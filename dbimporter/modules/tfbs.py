# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Database importing logic for the TFBSFinder module from antiSMASH"""

from antismash.modules import tfbs_finder

from dbimporter.common import RecordData

_CONFIDENCE_IDS = {}
_REGULATOR_IDS = {}


def get_regulator_id(cursor, name):
    """Get the TFBS regulator_id given the name of the regulator."""
    if name not in _REGULATOR_IDS:
        cursor.execute("SELECT regulator_id FROM antismash.regulators WHERE name = %s", (name,))
        ret = cursor.fetchone()
        if ret is None:
            raise ValueError("could not find matching regulator entry for: %s" % name)
        _REGULATOR_IDS[name] = ret[0]
    return _REGULATOR_IDS[name]


def get_confidence_id(cursor, name):
    """Get the TFBS confidence level id given the confidence name."""
    if name not in _CONFIDENCE_IDS:
        cursor.execute("SELECT confidence_id FROM antismash.regulator_confidence WHERE name = %s", (name,))
        ret = cursor.fetchone()
        if ret is None:
            raise ValueError("could not find matching regulator confidence entry for: %s" % name)
        _CONFIDENCE_IDS[name] = ret[0]
    return _CONFIDENCE_IDS[name]


def import_results(data: RecordData) -> None:
    results = data.module_results[tfbs_finder.__name__]
    assert isinstance(results, tfbs_finder.TFBSFinderResults)
    for region in data.record.get_regions():
        hits = results.get_hits_by_region(region.get_region_number())
        region_id = data.feature_mapping[region]
        for hit in hits:
            params = {
                "regulator_id": get_regulator_id(data.cursor, hit.name),
                "region_id": region_id,
                "score": hit.score,
                "start": hit.start,
                "confidence_id": get_confidence_id(data.cursor, str(hit.confidence).lower()),
            }
            data.insert("""
INSERT INTO antismash.binding_sites (regulator_id, region_id, score, start, confidence_id)
VALUES (%(regulator_id)s, %(region_id)s, %(score)s, %(start)s, %(confidence_id)s)""", params)
