# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Database importing logic for the ClusterCompare module from antiSMASH"""

from antismash.modules import cluster_compare

from dbimporter.common import RecordData

SCORE_THRESHOLD = 0.2
_INSERT_STMT = """
INSERT INTO antismash.cluster_compare_hits (region_id, protocluster_id, reference_accession, description, score, identity_metric, order_metric, components_metric)
VALUES (%(region_id)s, %(protocluster_id)s, %(reference_accession)s, %(description)s, %(score)s, %(identity)s, %(order)s, %(components)s)
"""


def insert_hit(data: RecordData, scorer: cluster_compare.data_structures.ReferenceScorer,
               region_id: int = None, protocluster_id: int = None) -> None:
    assert isinstance(scorer, cluster_compare.data_structures.ReferenceScorer)
    params = {
        "region_id": region_id,
        "protocluster_id": protocluster_id,
        "reference_accession": scorer.reference.accession,
        "description": scorer.reference.description,
        "score": scorer.final_score,
        "identity": scorer.identity,
        "order": scorer.order,
        "components": scorer.component,
    }
    data.insert(_INSERT_STMT, params)


def import_results(data: RecordData) -> None:
    results = data.module_results.get(cluster_compare.__name__)
    if not results:
        return
    assert isinstance(results, cluster_compare.ClusterCompareResults)
    mibig_hits = results.by_database["MIBiG"]
    for region in data.record.get_regions():
        region_id = data.feature_mapping[region]
        all_mode_hits = mibig_hits.by_region[region.get_region_number()]
        # start with the query region to reference region mode
        hits = all_mode_hits["RegionToRegion_RiQ"]
        assert isinstance(hits, cluster_compare.results.VariantResults)
        for scorer in hits.details.details:
            if scorer.final_score < SCORE_THRESHOLD:
                continue
            insert_hit(data, scorer, region_id=region_id)
        # then the query protocluster to reference region mode
        hits = all_mode_hits["ProtoToRegion_RiQ"]
        for protocluster in region.get_unique_protoclusters():
            protocluster_id = data.feature_mapping[protocluster]
            proto_score = hits.details.details.get(protocluster.get_protocluster_number(), {})
            for scorer in proto_score.values():
                if scorer.final_score < SCORE_THRESHOLD:
                    continue
                insert_hit(data, scorer, protocluster_id=protocluster_id)
