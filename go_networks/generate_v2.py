"""
Generate GO Networks from a list of GO terms and the Sif dump.
"""
import json
import logging
import pickle
from collections import defaultdict
from datetime import datetime
from itertools import combinations
from textwrap import dedent
from typing import Dict, Iterator, Optional, Set, Tuple, Union, List

import ndex2.client
import pandas as pd
import pystow
from neo4j.exceptions import CypherSyntaxError
from indra.ontology.bio import bio_ontology
from indra_cogex.client.neo4j_client import Neo4jClient
from indra_cogex.representation import Node
from indra_db.client.principal.curation import get_curations
from ndex2 import NiceCXNetwork, create_nice_cx_from_server
from tqdm import tqdm

from go_networks.util import (
    DIRECTED_TYPES,
    get_ndex_web_client,
    get_networks_in_set,
    NDEX_ARGS,
)
from go_networks.network_assembly import GoNetworkAssembler, get_cx_layout


# Derived types
Term = Tuple[str, str]
PropAgg = Tuple[Dict[str, int], Dict[str, int], Dict[str, int]]
PropDict = Dict[Tuple[str, str], PropAgg]
Go2Genes = Dict[str, Set[str]]
EvCountDict = Dict[str, Dict[str, int]]
NameEntityMap = Dict[str, Tuple[str, str]]

# Constants
cache = pystow.module("go_networks")
GO_MAPPINGS = cache.join(name="go_mappings.pkl")
COGEX_SIF = cache.join(name="cogex_sif.pkl")
PROPS_FILE = cache.join(name="props.pkl")
GO_NETWORKS = cache.join(name="networks.pkl")
NCX_CACHE = cache.module("ncx_cache")
DEFAULT_NDEX_SERVER = "http://ndexbio.org"
TEST_GO_ID = None

MIN_GENE_COUNT = 5
MAX_GENE_COUNT = 200
DB_SOURCES = ["biogrid", "hprd", "signor", "phosphoelm", "signor", "biopax"]

logger = logging.getLogger(__name__)


def get_curation_set() -> Set[int]:
    def _is_incorrect(stmt_hash, evid_hash):
        # Evidence is incorrect if it was only curated as incorrect
        if (
            evid_hash in correct_stmt_evid[stmt_hash]["incorrect"]
            and evid_hash not in correct_stmt_evid[stmt_hash]["correct"]
        ):
            return True
        return False

    correct_tags = {"correct", "hypothesis", "act_vs_amt"}

    # Download curations
    logger.info("Getting curations from the database")
    try:
        curations = get_curations()
    except Exception as e:
        logger.error(f"Could not get curations: {e}")
        return set()

    # If a statement only has incorrect curations, it is overall incorrect
    # If a statement has any correct curations, it is overall correct
    correct = {c["pa_hash"] for c in curations if c["tag"] in correct_tags}
    incorrect = {c["pa_hash"] for c in curations if c["pa_hash"] not in correct}

    correct_stmt_evid = {}
    for c in curations:
        pa_hash = c["pa_hash"]
        if pa_hash in correct:
            if pa_hash not in correct_stmt_evid:
                correct_stmt_evid[pa_hash] = defaultdict(set)
            if c["tag"] in correct_tags:
                correct_stmt_evid[pa_hash]["correct"].add(c["source_hash"])
            else:
                correct_stmt_evid[pa_hash]["incorrect"].add(c["source_hash"])

    hashes_out = set()
    for c in curations:
        pa_hash = c["pa_hash"]
        src_hash = c["source_hash"]
        if pa_hash not in incorrect:
            # Check evidences in detail
            if pa_hash in correct_stmt_evid:
                if _is_incorrect(pa_hash, src_hash):
                    hashes_out.add(pa_hash)
                else:
                    pass
            else:
                pass
        else:
            hashes_out.add(pa_hash)

    logger.info(f"Found {len(hashes_out)} hashes to filter out")

    return hashes_out


def get_sif_from_cogex(limit: Optional[int] = None, with_apoc: bool = True) -> pd.DataFrame:
    """Get the SIF from the database

    Conditions:
        - No self loops, i.e. A != B
        - Only HGNC nodes
        - Skip relations where evidence_count == 1 AND source is a reader
        - Skip relations where stmt_type == 'Complex' AND sparser is the only
          source, for any evidence count

    Parameters
    ----------
    limit :
        Limit the number of edges to return. Useful for testing or debugging.
    with_apoc :
        If True, use the apoc library to perform some of the filtering directly
        in the query. If False, perform that filtering in Python.

    Returns
    -------
    :
        A pandas DataFrame with the SIF data
    """
    query_fmt = dedent("""\
    MATCH (gene1:BioEntity)-[r:indra_rel]-(gene2:BioEntity)
    {apoc_with_clause}
    WHERE
        // Filter out self loops
        gene1 <> gene2 AND
        // Filter out non-HGNC nodes
        gene1.id CONTAINS 'hgnc' AND
        gene2.id CONTAINS 'hgnc'
        {apoc_filter}
    RETURN gene1, gene2, r.belief, r.evidence_count, r.source_counts, r.stmt_hash, r.stmt_type
    """)
    if limit is not None and isinstance(limit, int):
        query_fmt += f"LIMIT {limit}"

    if with_apoc:
        apoc_with_clause = (
            "WITH gene1, gene2, r, "
            "apoc.convert.fromJsonMap(r.source_counts) AS source_counts"
        )
        apoc_filter = dedent(f"""\
            // Filter out Complex statements with only sparser as source
            AND NOT (r.stmt_type = 'Complex' AND keys(source_counts) = ['sparser'])
            // Filter out relations with only one evidence and from a reader
            AND NOT (
                r.evidence_count = 1 AND
                NOT apoc.coll.intersection(
                    keys(source_counts),
                    {DB_SOURCES}
                )
            )""")
    else:
        apoc_filter = ""
        apoc_with_clause = ""

    query = query_fmt.format(
        apoc_with_clause=apoc_with_clause, apoc_filter=apoc_filter
    )

    n4j_client = Neo4jClient()
    logger.info("Getting SIF interaction data from database")
    if with_apoc:
        try:
            results = n4j_client.query_tx(query)
        except CypherSyntaxError:
            logger.warning("Failed to get SIF from database with apoc query, trying without")
            with_apoc = False
            query = query_fmt.format(apoc_with_clause="", apoc_filter="")

    if not with_apoc:
        results = n4j_client.query_tx(query)

    res_tuples = []
    logger.info("Generating SIF from database results")
    for r in tqdm(results):
        gene1 = n4j_client.neo4j_to_node(r[0])
        gene2 = n4j_client.neo4j_to_node(r[1])
        if not with_apoc:
            # Filter out Complex statements with only sparser as source
            source_counts = json.loads(r[4])
            if r[6] == "Complex" and {"sparser"} == set(source_counts):
                continue
            # Filter out relations with only one evidence and from a reader
            # (i.e. not in DB_SOURCES)
            if int(r[3]) == 1 and not set(source_counts).intersection(DB_SOURCES):
                continue
        res_tuples.append(
            (
                gene1.db_ns,  # agA_ns
                gene1.db_id,  # agA_id
                gene1.data["name"],  # agA_name
                gene2.db_ns,  # agB_ns
                gene2.db_id,  # agB_id
                gene2.data["name"],  # agB_name
                r[2],  # belief
                r[3],  # evidence_count
                r[4],  # source_counts
                r[5],  # stmt_hash
                r[6],  # stmt_type
            )
        )

    logger.info("Converting to DataFrame")
    df = pd.DataFrame(
        res_tuples,
        columns=[
            "agA_ns",
            "agA_id",
            "agA_name",
            "agB_ns",
            "agB_id",
            "agB_name",
            "belief",
            "evidence_count",
            "source_counts",
            "stmt_hash",
            "stmt_type",
        ],
    )
    logger.info("Setting data types")
    df["belief"] = df["belief"].astype(pd.Float64Dtype())
    df["evidence_count"] = df["evidence_count"].astype(pd.Int64Dtype())
    df["stmt_hash"] = df["stmt_hash"].astype(pd.Int64Dtype())

    # Filter out curated incorrect statements
    wrong_hashes = get_curation_set()
    if wrong_hashes:
        logger.info("Filtering out statements curated as incorrect")
        df = df[~df.stmt_hash.isin(wrong_hashes)]
    return df


def get_sif(regenerate: bool = False) -> pd.DataFrame:
    if not regenerate and COGEX_SIF.exists():
        logger.info("Loading SIF from cache")
        with COGEX_SIF.open("rb") as bf:
            cogex_sif = pickle.load(file=bf)
        assert isinstance(cogex_sif, pd.DataFrame), "Cached SIF is not a dataframe"

    else:
        cogex_sif = get_sif_from_cogex()
    return cogex_sif


def generate_props(
    regenerate: bool = False,
) -> PropDict:
    """Generate properties per pair from the Sif dump

    For each pair of genes (A,B) (excluding self loops), generate the
    following properties:
        - "SOURCE => TARGET": aggregate number of evidences by statement
          type for A->B statements
        - "TARGET => SOURCE": aggregate number of evidences by statement
          type for B->A statements
        - "SOURCE - TARGET": aggregate number of evidences by statement type
          for A-B undirected statements

    Parameters
    ----------
    regenerate :
        Whether to regenerate the props. If True, the cache will be overwritten.

    Returns
    -------
    :
        The properties dictionary
    """
    # If the props file exists, load it, unless we want to regenerate
    if not regenerate and PROPS_FILE.exists():
        logger.info(f"Loading property lookup from {PROPS_FILE}")
        with PROPS_FILE.open(mode="rb") as fr:
            props_by_pair = pickle.load(fr)
    else:
        logger.info("Generating property lookup")

        # Load SIF
        sif_df = get_sif(regenerate)

        hashes_by_pair = defaultdict(set)
        props_by_hash = {}

        def get_direction(row, pair):
            directed = row.stmt_type in DIRECTED_TYPES
            direction = row.agA_name == pair[0]
            if directed:
                if direction:
                    return "forward"
                else:
                    return "reverse"
            else:
                return "undirected"

        # For each pair of genes, generate the properties
        logger.info("Generating properties by hash")
        for _, row in tqdm(sif_df.iterrows(), total=sif_df.shape[0]):
            pair = tuple(sorted([row.agA_name, row.agB_name]))
            hashes_by_pair[pair].add(row.stmt_hash)
            if row.stmt_hash not in props_by_hash:
                props_by_hash[row.stmt_hash] = {
                    "ev_count": row.evidence_count,
                    "stmt_type": row.stmt_type,
                    "direction": get_direction(row, pair),
                }
        hashes_by_pair = dict(hashes_by_pair)

        def aggregate_props(props) -> PropAgg:
            ev_forward = defaultdict(int)
            ev_reverse = defaultdict(int)
            ev_undirected = defaultdict(int)
            for prop in props:
                if prop["direction"] == "forward":
                    ev_forward[prop["stmt_type"]] += prop["ev_count"]
                elif prop["direction"] == "reverse":
                    ev_reverse[prop["stmt_type"]] += prop["ev_count"]
                else:
                    ev_undirected[prop["stmt_type"]] += prop["ev_count"]
            return dict(ev_forward), dict(ev_reverse), dict(ev_undirected)

        props_by_pair = {}
        logger.info("Aggregating properties by pair")
        for pair, hashes in tqdm(hashes_by_pair.items()):
            props_by_pair[pair] = aggregate_props([props_by_hash[h] for h in hashes])

        # Write to file if provided
        with PROPS_FILE.open(mode="wb") as fo:
            logger.info(f"Saving property lookup to {PROPS_FILE}")
            pickle.dump(obj=props_by_pair, file=fo)

    return props_by_pair


def go_term_gene_query() -> Iterator[Tuple[str, str]]:
    """An Iterator of go-term gene pairs

    Returns
    -------
    :
        An iterable of pairs of go term node with associated gene node
    """
    query = (
        "MATCH (gene:BioEntity)-[:associated_with]->(term:BioEntity) "
        "RETURN term.id, gene.name"
    )
    client = Neo4jClient()
    return client.query_tx(query)


def genes_by_go_id(regenerate: bool = False) -> Dict[str, Set[str]]:
    """Map go ids to gene symbols

    Parameters
    ----------
    regenerate :
        Whether to regenerate the cache. If True, the cache will be overwritten.

    Returns
    -------
    :
        A dictionary mapping go ids to sets of gene symbols
    """
    # Get all go IDs from the database, unless they're cached
    if GO_MAPPINGS.exists() and not regenerate:
        logger.info("Loading GO mapping from cache")
        with GO_MAPPINGS.open(mode="rb") as fr:
            return pickle.load(fr)

    logger.info("Loading GO mapping from database")

    # Set initial mapping
    genes_by_go = defaultdict(set)
    for go_curie, gene_symbol in tqdm(go_term_gene_query(),
                                      desc="Loading GO-gene associations"):
        # Use upper case for GO IDs so we get GO:0001234 instead of go:0001234 so that
        # it corresponds to the bio_ontology
        genes_by_go[go_curie.upper()].add(gene_symbol)

    # Load bio ontology
    logger.info("Adding genes of child terms to the parent terms")
    bio_ontology.initialize()

    # For each term, add the genes associated with its children as well
    for go_id in tqdm(set(genes_by_go.keys()), desc="Adding genes of child terms"):
        # children come out as ('GO', 'GO:0001234'), use only GO:0001234
        for db_ns, go_child in bio_ontology.get_children("GO", go_id.upper()):
            genes_by_go[go_id] |= genes_by_go[go_child]

    # Reset defaultdict to dict
    genes_by_go = dict(genes_by_go)

    # Save to cache
    with GO_MAPPINGS.open(mode="wb") as fw:
        logger.info("Caching GO mappings")
        pickle.dump(obj=genes_by_go, file=fw)

    return genes_by_go


def build_networks(
    go2genes_map: Go2Genes,
    pair_props: PropDict,
) -> Dict[str, Dict[str, Union[NiceCXNetwork, float]]]:
    """Build networks per go-id associated genes

    Parameters
    ----------
    go2genes_map :
        A dict mapping GO ID to a list of genes
    pair_props :
        Lookup for edges

    Returns
    -------
    :
        Dict of assembled networks by go id
    """
    networks = {}
    skipped = 0
    # Only pass the relevant parts of the pair_props dict
    for go_id, gene_set in tqdm(go2genes_map.items(), total=len(go2genes_map)):
        if TEST_GO_ID and go_id != TEST_GO_ID:
            continue

        def _key(g1, g2):
            return tuple(sorted([g1, g2]))

        # Get relevant pairs from pair_properties
        prop_dict = {
            _key(g1, g2): pair_props[_key(g1, g2)]
            for g1, g2 in combinations(gene_set, 2)
            if _key(g1, g2) in pair_props
        }

        if not prop_dict:
            tqdm.write(f"No statements for ID {go_id}")
            skipped += 1
            continue

        gna = GoNetworkAssembler(
            identifier=go_id,
            entity_list=list(gene_set),
            pair_properties=prop_dict,
        )
        gna.assemble()
        networks[go_id] = {
            "network": gna.network,
            "max_score": max(gna.rel_scores),
            "min_score": min(gna.rel_scores),
        }

    logger.info(f"Skipped {skipped} networks without statements")
    return networks


def filter_self_loops(df):
    """Remove self-loops from a dataframe

    Parameters
    ----------
    df :
        The dataframe to filter

    Returns
    -------
    :
        The filtered dataframe
    """
    return df[df.agA_name != df.agB_name]


def filter_go_ids(go2genes_map) -> Dict[str, Set[str]]:
    return {
        go_id: genes
        for go_id, genes in tqdm(go2genes_map.items(), desc="Filtering GO IDs")
        if MIN_GENE_COUNT <= len(genes) <= MAX_GENE_COUNT
    }


def generate(
    regenerate: bool = False,
) -> Dict[str, Dict[str, Union[NiceCXNetwork, float]]]:
    """Generate new GO networks from INDRA statements

    Parameters
    ----------
    regenerate :
        If True, regenerate the props from scratch and overwrite the existing
        props file, if it exists.
    """
    # Make genes by GO ID dict
    go2genes_map = genes_by_go_id(regenerate=regenerate)

    # Filter GO IDs
    go2genes_map = filter_go_ids(go2genes_map)

    # Generate properties
    sif_props = generate_props(regenerate=regenerate)

    # Iterate by GO ID and for each list of genes, build a network
    return build_networks(go2genes_map, sif_props)


def download_ncx_from_uuids(ncx_uuids):
    # Download the set of NiceCXnetworks given by the UUIDS
    res = {}
    for ncx_uuid in tqdm(ncx_uuids, desc="Downloading NCX"):
        ncx = create_nice_cx_from_server(uuid=ncx_uuid, **NDEX_ARGS)
        res[ncx_uuid] = ncx

    return res


def get_ncx_cache_from_set(ncx_set_uuid: str, refresh=False):
    """Get the set of ncx networks from a network set and cache them"""
    ncx_cache = {}
    ndex_web_client = get_ndex_web_client()
    ncx_uuids = get_networks_in_set(ncx_set_uuid, client=ndex_web_client)
    for ncx_uuid in tqdm(ncx_uuids,
                         desc=f"Downloading NCX from set {ncx_set_uuid}"):
        ncx_file = NCX_CACHE.join(name=f"{ncx_uuid}.pkl")
        if not refresh and ncx_file.is_file():
            with ncx_file.open("rb") as f:
                ncx = pickle.load(f)
        else:
            ncx = create_nice_cx_from_server(uuid=ncx_uuid, **NDEX_ARGS)
            with ncx_file.open("wb") as f:
                pickle.dump(obj=ncx, file=f)

        # network_info = ndex_web_client.get_network_summary(ncx_uuid)

        # Get the go id
        # go_id = None
        # for prop_dict in network_info["properties"]:
        #     if prop_dict["predicateString"] == "GO ID":
        #         go_id = prop_dict["value"]
        #         break
        #     else:
        #         continue
        # if go_id is None:
        #     logger.warning(f"No GO ID found for network {ncx_uuid}")
        #     continue

        # ncx_cache[go_id] = {"uuid": ncx_uuid, "ncx": ncx}
        ncx_cache[ncx_uuid] = ncx

    return ncx_cache


def get_go_uuid_mapping(set_uuid: str) -> Dict[str, str]:
    """Get mapping of GO IDs to UUID from NDEx

    Parameters
    ----------
    set_uuid :
        The UUID of the GO set

    Returns
    -------
    :
        A dict mapping GO IDs to UUIDs of networks
    """
    logger.info(f"Getting GO IDs for {set_uuid}")

    # Loop the networks in the set and get the go id for each network
    ndex_web_client = get_ndex_web_client()
    uuid_set = get_networks_in_set(network_set_id=set_uuid, client=ndex_web_client)
    go_uuid_mapping = {}

    logger.info("Getting GO ID-uuid mapping")
    for cx_uuid in tqdm(uuid_set):
        # Get the info for the network
        network_info = ndex_web_client.get_network_summary(cx_uuid)

        # Get the go id
        go_id = None
        for prop_dict in network_info["properties"]:
            if prop_dict["predicateString"] == "GO ID":
                go_id = prop_dict["value"]
                go_uuid_mapping[go_id] = cx_uuid
                break
            else:
                continue
        if go_id is None:
            logger.warning(f"No GO ID found for network {cx_uuid}")
    return go_uuid_mapping


def update_coordinates_for_network_set(set_uuid: str):
    """Update the coordinates for a network set

    Returns the set of UUIDs that failed a second time
    """
    ndex_web_client = get_ndex_web_client()

    # Get the network UUIDs for the given network set
    uuid_set = get_networks_in_set(network_set_id=set_uuid, client=ndex_web_client)

    failed = set()
    for network_uuid in tqdm(
            uuid_set, total=len(uuid_set), desc="Updating coordinates"
    ):
        try:
            update_coordinates_for_network(ncx_uuid=network_uuid,
                                           ndex_client=ndex_web_client)
        except Exception:
            failed.add(network_uuid)

    second_fail = set()
    if failed:
        logger.info(f"{len(failed)} updates failed, retrying after 10 s sleep")
        from time import sleep
        sleep(10)
        for network_uuid in tqdm(
                failed, total=len(failed), desc="Retrying updating coordinates"
        ):
            try:
                update_coordinates_for_network(ncx_uuid=network_uuid,
                                               ndex_client=ndex_web_client)
            except Exception as err:
                tqdm.write(
                    f"Failed a second time for uuid {network_uuid}: {err}"
                )
                second_fail.add(network_uuid)
    if second_fail:
        logger.warning(f"Returning {len(second_fail)} network UUIDs that "
                       f"failed to update a second time.")
    return second_fail


def update_coordinates_for_network(
    ncx_uuid: str,
    ndex_client: ndex2.client.Ndex2
):
    """Given an ncx uuid update the coordinates of the nodes in the graph"""
    # Get the NiceCXNetwork
    ncx: NiceCXNetwork = create_nice_cx_from_server(uuid=ncx_uuid, **NDEX_ARGS)

    # Get the coordinates
    node_layout_by_name = get_cx_layout(network=ncx)

    # Update the coordinates in the network
    node_name_id = {(nd['n'], _id) for _id, nd in ncx.nodes.items()}
    layout_aspect = []
    for hgnc_symb, node_id in node_name_id:
        x, y = node_layout_by_name[hgnc_symb]
        layout_aspect.append({"node": node_id, "x": x, "y": y})
    ncx.set_opaque_aspect("cartesianLayout", layout_aspect)

    # Update the network
    ndex_client.update_cx_network(cx_stream=ncx.to_cx_stream(),
                                  network_id=ncx_uuid)


def add_to_new_set(add_from_set_uuid: str,
                   new_set_uuid: str,
                   ndex_client: Optional[ndex2.Ndex2] = None):
    """Add networks from one set to another"""
    # Get client if not provided
    t = tqdm(desc="Adding networks to new set", total=3)
    ndex_client = get_ndex_web_client() if ndex_client is None else ndex_client
    t.update()

    # Get the networks in the input set
    network_uuids = get_networks_in_set(network_set_id=add_from_set_uuid,
                                        client=ndex_client)
    t.update()

    # Add the networks to the new set
    # NOTE: the docstring typing in 'add_networks_to_networkset' says it
    # returns None but a quick look at the code shows that it returns
    # res.text or res.json() from a 'requests' response
    res = ndex_client.add_networks_to_networkset(set_id=new_set_uuid,
                                                 networks=network_uuids)
    t.update()
    t.close()
    return res


def get_networks_and_date_for_user(
        user: str, ndex_client: Optional[ndex2.Ndex2] = None
) -> List[Dict[str, Union[datetime, str]]]:
    """Get network UUIDs together with their creation and modification date"""
    def _get_dt(raw_ts: int) -> datetime:
        str_ts = str(raw_ts)
        return datetime.fromtimestamp(float(f"{str_ts[:10]}.{str_ts[10:]}"))

    def _get_list(sums):
        out = []
        for s in sums:
            out.append({"uuid": s["externalId"],
                        "modificationTime": _get_dt(s["modificationTime"]),
                        "creationTime": _get_dt(s["creationTime"])})
        return out

    # Get client if not provided
    if ndex_client is None:
        ndex_client = get_ndex_web_client()

    res = []
    summaries = ndex_client.get_user_network_summaries(username=user)
    res += _get_list(summaries)

    # If the returned list has 1000 entries, it was likely cut off and we
    # need to make a new request with a new offset
    offset = 1000
    while len(summaries) == 1000:
        summaries = ndex_client.get_user_network_summaries(username=user,
                                                           offset=offset)
        res += _get_list(summaries)
        offset += 1000

    # Check for double counting
    uuids = set()
    for d in res:
        uuids.add(d["uuid"])

    assert len(uuids) == len(res), "Double counted UUIDs"
    return res


def get_networks_before_date_for_user(
        user: str, cutoff_date: datetime, ndex_client: Optional[ndex2.Ndex2] =
        None
):
    """Get all networks created and modified before the provided date"""
    # Check that the cutoff date is valid
    utcnow = datetime.utcnow()
    if cutoff_date > utcnow:
        raise ValueError(
            f"Cutoff date is in the future (year={cutoff_date.year}, month="
            f"{cutoff_date.month}, day={cutoff_date.day}), please correct "
            f"the cutoff date to be in the past"
        )

    # Get client if not provided
    if ndex_client is None:
        ndex_client = get_ndex_web_client()

    # Get the networks UUIDs and their datetimes
    network_dt_list = get_networks_and_date_for_user(user=user,
                                                     ndex_client=ndex_client)

    res = []
    for dt_dict in network_dt_list:
        if dt_dict["modificationTime"] < cutoff_date and \
                dt_dict["creationTime"] < cutoff_date:
            res.append(dt_dict)

    return res


def delete_networks(network_uuids: Set[str], ndex_client: ndex2.Ndex2):
    """Delete the networks provided by the UUIDs, return those that failed"""
    failed = set()
    for uuid in tqdm(network_uuids, desc="Deleting networks"):
        try:
            ndex_client.delete_network(uuid)
        except Exception:
            failed.add(uuid)

    if failed:
        logger.warning(f"{len(failed)} deletions failed")

    return failed


def format_and_update_network(
    ncx: NiceCXNetwork,
    network_set_id: str,
    style_ncx: NiceCXNetwork,
    ndex_client: ndex2.client.Ndex2,
    cx_uuid: Optional[str] = None,
) -> Tuple[str, Dict[str, bool]]:
    """Take a NiceCXNetwork and upload it to NDEx."""
    # Fixme: NiceCXNetwork.apply_style_from_network() does not exist in
    #  ndex2==2.0.1 but ==2.0.1 is the requirement for INDRA
    #  This method was added in 3.1.0:
    #  https://github.com/ndexbio/ndex2-client/issues/43
    ncx.apply_style_from_network(style_ncx)
    failed_public = False
    failed_update = False
    # If we have a UUID, update the network
    if cx_uuid:
        try:
            ndex_client.update_cx_network(
                cx_stream=ncx.to_cx_stream(), network_id=cx_uuid
            )
        except Exception as e:
            logger.warning(f"Failed to update network {cx_uuid}: {e}")
            failed_update = True
        finally:
            network_id = cx_uuid
    # If there is no UUID, create a new network
    else:
        network_url = ncx.upload_to(client=ndex_client)
        network_id = network_url.split("/")[-1]
        try:
            ndex_client.make_network_public(network_id)
        except Exception as e:
            logger.warning(f"Failed to make network {network_id} public: {e}")
            failed_public = True

        try:
            ndex_client.add_networks_to_networkset(network_set_id, [network_id])
        except Exception as e:
            logger.warning(
                f"Failed to add network {network_id} to network set "
                f"{network_set_id}: {e}"
            )

    return network_id, {"public": failed_public, "update": failed_update}


def _update_style_network(style_ncx: NiceCXNetwork, min_score: float, max_score: float):
    if NiceCXNetwork.CY_VISUAL_PROPERTIES in style_ncx.opaqueAspects:
        vis_prop = NiceCXNetwork.CY_VISUAL_PROPERTIES
        for style_dict in style_ncx.opaqueAspects[vis_prop]:
            if style_dict["properties_of"] == "edges:default":
                style_str = style_dict["mappings"]["EDGE_WIDTH"]["definition"]
                styles = style_str.split(",")
                min_ix = -1
                max_ix = -1
                for st in styles:
                    if "OV=0=" in st:
                        min_ix = styles.index(st)
                    if "OV=1=" in st:
                        max_ix = styles.index(st)
                    if min_ix != -1 and max_ix != -1:
                        break
                styles[min_ix] = f"OV=0={min_score}"
                styles[max_ix] = f"OV=1={max_score}"
                style_dict["mappings"]["EDGE_WIDTH"]["definition"] = ",".join(styles)
                break

        assert f"OV=0={min_score}" in style_dict["mappings"]["EDGE_WIDTH"]["definition"]
        assert f"OV=1={max_score}" in style_dict["mappings"]["EDGE_WIDTH"]["definition"]


def main(
    network_set: str,
    style_network: str,
    regenerate: bool,
    test_go_term: Optional[str] = None,
    ndex_server_style: str = DEFAULT_NDEX_SERVER,
):
    global TEST_GO_ID
    logger.info(f"Using network set id {network_set}")
    logger.info(f"Using network style {style_network}")
    if test_go_term:
        logger.info(f"Testing GO term {test_go_term}")
        TEST_GO_ID = test_go_term

    logger.info(f"Using ndex server {ndex_server_style} for style")

    networks = generate(regenerate=regenerate)

    # Only cache the networks if we're not testing
    if not test_go_term:
        with GO_NETWORKS.open("wb") as f:
            logger.info(f"Writing networks to {GO_NETWORKS}")
            pickle.dump(networks, f)
    else:
        logger.info("TEST: Not caching networks")

    style_ncx = create_nice_cx_from_server(server=ndex_server_style, uuid=style_network)

    ndex_web_client = get_ndex_web_client()

    go_uuid_mapping = get_go_uuid_mapping(network_set)

    failed_to_set_public = []
    failed_to_update = []
    logger.info(f"Uploading {len(networks)} networks to NDEx")
    for go_id, network_dict in tqdm(sorted(networks.items(), key=lambda x: x[0])):
        network = network_dict["network"]
        min_score = network_dict["min_score"]
        max_score = network_dict["max_score"]

        # Update style network
        _update_style_network(style_ncx, min_score=min_score, max_score=max_score)

        # Get uuid for GO term
        go_uuid = go_uuid_mapping.get(go_id)

        # Update/upload network
        network_id, failed = format_and_update_network(
            ncx=network,
            network_set_id=network_set,
            style_ncx=style_ncx,
            ndex_client=ndex_web_client,
            cx_uuid=go_uuid,
        )
        if failed["public"]:
            failed_to_set_public.append(network_id)
        if failed["update"]:
            failed_to_update.append(network_id)

        if TEST_GO_ID:
            print(f"Testing network uuid {network_id}")

    if failed_to_set_public:
        logger.warning(f"Failed to set {len(failed_to_set_public)} networks public")
        logger.info("Retrying to set networks public...")

        # Retry to set public
        for failed_uuid in failed_to_set_public:
            try:
                ndex_web_client.make_network_public(network_id=failed_uuid)
            except Exception as e:
                logger.warning(f"Failed to set network {failed_uuid} public again: {e}")

    if failed_to_update:
        logger.warning(f"Failed to update {len(failed_to_update)} networks")
