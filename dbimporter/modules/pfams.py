# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Database importing logic for Pfam domains """

from antismash.common.secmet.features import PFAMDomain

from dbimporter.common.record_data import RecordData


def import_results(data: RecordData) -> None:
    """ Import all Pfam results for a record """
    for pfam in data.record.get_pfam_domains():
        handle_pfamdomain(data, pfam)


def handle_pfamdomain(data: RecordData, feature: PFAMDomain) -> None:
    """Handle PFAM_domain features."""
    params = {}

    params['location'] = str(feature.location)

    params['score'] = feature.score
    params['evalue'] = feature.evalue
    params['translation'] = feature.translation
    params['locus_tag'] = feature.locus_tag
    params['detection'] = feature.detection
    params['database'] = feature.database
    params['pfam_id'] = get_pfam_id(data.cursor, feature.identifier)
    params['location'] = str(feature.location)
    params['cds_id'] = data.feature_mapping[data.record.get_cds_by_name(feature.locus_tag)]

    pfam_id = data.insert("""
INSERT INTO antismash.pfam_domains (
    database,
    detection,
    score,
    evalue,
    translation,
    pfam_id,
    location,
    cds_id
) VALUES (
    %(database)s,
    %(detection)s,
    %(score)s,
    %(evalue)s,
    %(translation)s,
    %(pfam_id)s,
    %(location)s,
    %(cds_id)s
)
RETURNING pfam_domain_id""", params)

    if not feature.gene_ontologies:
        return

    statement = """
INSERT INTO antismash.pfam_go_entries (pfam_domain_id, go_id)
VALUES (%s, %s)"""
    for go_identifier in feature.gene_ontologies.ids:
        go_id = get_go_id(data.cursor, go_identifier)
        data.insert(statement, (pfam_id, go_id))


def get_pfam_id(cur, identifier: str) -> id:
    """Get the pfam_id for a domain by db_xref string."""
    assert identifier
    cur.execute("SELECT pfam_id FROM antismash.pfams WHERE pfam_id = %s", (identifier,))
    ret = cur.fetchone()
    if ret is None:
        raise ValueError("Invalid pfam_id {!r}".format(identifier))

    return ret[0]


def get_go_id(cursor, go_identifier: str) -> int:
    """Get a gene ontology entry ID from the GO identifier"""
    cursor.execute("SELECT go_id FROM antismash.gene_ontologies WHERE identifier = %s", (go_identifier,))
    ret = cursor.fetchone()
    if ret is None:
        raise ValueError("unknown GO identifier %s" % go_identifier)
    return ret[0]
