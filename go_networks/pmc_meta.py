"""Get PMC XML from the DB and extract metadata

# DATA COLUMNS USEFUL FOR BIOFACTOID

* journal
* article
* email
* corresponding_author (as Boolean)
* year
* indra_statement_count

# XPATH STATEMENTS

## Authors

.//contrib[@contrib-type="author"]

## Corresponding Author

.//contrib[@contrib-type="author" and @corresp="yes"]

## Author Sub-Properties (Per Author)

For what I have done in the past, I have just taken the first entries

### Last Names

.//surname

### First Name

.//given-names

### Email

.//email
"""
import codecs
import logging

from lxml import etree

from indra_db_lite.construction import query_to_csv
from gzip import decompress

logger = logging.getLogger(__name__)


def pmc_xml_to_csv(out_path):
    """Get PMC XML from the DB and extract metadata

    Parameters
    ----------
    out_path : str
        Path to the output CSV file.
    """
    query = """
    SELECT
        text_ref_id, encode(content, 'hex')
    FROM 
        text_content
    WHERE
        format = 'xml'
        AND
        text_ref_id IN (
            SELECT
                id
            FROM
                text_ref
            WHERE
                pmcid IS NOT NULL
        )
    """
    query_to_csv(query, out_path)


def text_ref_id_pmc_id_dump(out_path):
    """Get text_ref_id and pmcid from the DB

    Parameters
    ----------
    out_path : str
        Path to the output CSV file.
    """
    query = """
    SELECT
        id, pmcid
    FROM 
        text_ref
    WHERE
        pmcid IS NOT NULL
    """
    query_to_csv(query, out_path)


def hex_bin_to_str(raw_hex_bin: str) -> str:
    """Convert hex-encoded raw string to string

    Parameters
    ----------
    raw_hex_bin :
        Hex-encoded bytes as a plain string.

    Returns
    -------
    str
        String
    """
    decode_hex = codecs.getdecoder("hex_codec")
    # It's a plain text string containing hex-encoded bytes, so it's first
    # two characters are escaping the hex-encoding: '\\x1f8b0808......'
    hex_str = decode_hex(raw_hex_bin[2:])[0]
    return decompress(hex_str).decode()


def extract_info_from_pmc_xml(xml_str: str) -> dict:
    """Extract metadata from PMC XML

    * journal
    * article
    * email
    * corresponding_author (as Boolean)
    * year

    Parameters
    ----------
    xml_str :
        PMC XML as a string

    Returns
    -------
    :
    """
    tree = etree.fromstring(xml_str)

    # Authors
    def _get_corresponding_author(root):
        authors = root.xpath(".//contrib[@contrib-type='author']")
        # Find the corresponding author and return the email
        for author in authors:
            if author.attrib.get('corresp', 'no') == 'yes':
                return author.xpath(".//email")[0].text

    # Corresponding author
    email = _get_corresponding_author(tree)

    # Journal name
    journal = tree.xpath(".//journal-title")[0].text

    # Article title
    # FixMe: sometimes the title and the abstract are in the article-title tag
    article_title = tree.xpath(".//article-title")[0].text

    # Year
    year = tree.xpath(".//pub-date[@pub-type='epub']")[0].xpath(".//year")[0].text

    return {
        'journal': journal,
        'article': article_title,
        'email': email,
        'corresponding_author': email is not None,
        'year': year,
    }
