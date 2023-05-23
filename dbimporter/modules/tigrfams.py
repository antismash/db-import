# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Database importing logic for TIGRFam domains """

from antismash.detection import tigrfam
from antismash.detection.tigrfam.tigr_domain import TIGRDomain

from dbimporter.common.record_data import RecordData


def import_results(data: RecordData) -> None:
    """ Import all TIGRFam results for a record """
    results = data.module_results.get(tigrfam.__name__)
    if not results:
        return
    for domain in data.record.get_antismash_domains_by_tool(tigrfam.__name__.split(".")[-1]):
        handle_tigrfam_domain(data, domain)


def handle_tigrfam_domain(data: RecordData, feature: TIGRDomain) -> int:
    """Handle TIGRfam domain features."""
    params = {}

    params['location'] = str(feature.location)

    params['score'] = feature.score
    params['evalue'] = feature.evalue
    params['translation'] = feature.translation
    params['locus_tag'] = feature.locus_tag
    params['detection'] = feature.detection
    params['database'] = feature.database
    params['tigrfam_id'] = get_tigrfam_id(data.cursor, feature.identifier)
    params['location'] = str(feature.location)
    params['cds_id'] = data.feature_mapping[data.record.get_cds_by_name(feature.locus_tag)]

    return data.insert("""
INSERT INTO antismash.tigrfam_domains (
    database,
    detection,
    score,
    evalue,
    translation,
    tigrfam_id,
    location,
    cds_id
) VALUES (
    %(database)s,
    %(detection)s,
    %(score)s,
    %(evalue)s,
    %(translation)s,
    %(tigrfam_id)s,
    %(location)s,
    %(cds_id)s
)
RETURNING tigrfam_domain_id""", params)


def get_tigrfam_id(cur, identifier: str) -> id:
    """ Get a TIGRFam ID.
        Really just ensures that the identifier is present, as if it is present
        the return value and the identifier argument will be identical.
    """
    assert identifier
    cur.execute("SELECT tigrfam_id FROM antismash.tigrfams WHERE tigrfam_id = %s", (identifier,))
    ret = cur.fetchone()
    if ret is None:
        raise ValueError("Invalid tigrfam_id {!r}".format(identifier))

    return ret[0]
