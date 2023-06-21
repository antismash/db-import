# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Database importing logic for the ClusterBlast variants """

from antismash.common.secmet import Region
from antismash.modules import clusterblast

from dbimporter.common import RecordData


def store_clusterblast(data: RecordData, results, algorithm):
    """Store XClusterBlast results in the database."""
    data.cursor.execute("SELECT algorithm_id FROM antismash.clusterblast_algorithms WHERE name = %s", (algorithm, ))
    ret = data.cursor.fetchone()
    if ret is None:
        raise ValueError('Did not find algorithm_id for {!r}!'.format(algorithm))
    algorithm_id = ret[0]

    region = data.current_region
    region_string = f"c{int(region.location.start)}-{int(region.location.end)}"

    assert isinstance(results, clusterblast.results.RegionResult), type(results)
    # limit to the number of drawn hits, normally set with DEFAULT_AS_OPTIONS.cb_nclusters
    for i, hit in enumerate(results.ranking[:len(results.svg_builder.hits)]):
        ref_cluster, _ = hit
        # skip hits that refer back to this specific region, but allow other hits to this record
        if ref_cluster.accession == data.record.id and ref_cluster.cluster_label == region_string:
            continue
        params = {
            "algorithm_id": algorithm_id,
            "region_id": data.current_region_id,
            "rank": i + 1,
            "acc": ref_cluster.accession + (("_" + ref_cluster.cluster_label) if algorithm != "knownclusterblast" else ""),
            "description": ref_cluster.description,
            "similarity": results.svg_builder.hits[i].similarity,
        }
        data.insert("""
INSERT INTO antismash.clusterblast_hits (rank, region_id, acc, description, similarity, algorithm_id)
VALUES (%(rank)s, %(region_id)s, %(acc)s, %(description)s, %(similarity)s, %(algorithm_id)s)
RETURNING clusterblast_hit_id
    """, params)


def import_region_results(data: RecordData, region: Region, deferred: bool = False) -> None:
    clusterblast_results = data.module_results.get(clusterblast.__name__)

    if not clusterblast_results:
        return

    assert data.current_region == region

    region_index = region.get_region_number() - 1  # 1-indexed to 0-indexed for a list
    if deferred:
        if clusterblast_results.general:
            store_clusterblast(data, clusterblast_results.general.region_results[region_index], 'clusterblast')
    else:
        if clusterblast_results.knowncluster:
            store_clusterblast(data, clusterblast_results.knowncluster.region_results[region_index], 'knownclusterblast')
        if clusterblast_results.subcluster:
            store_clusterblast(data, clusterblast_results.subcluster.region_results[region_index], 'subclusterblast')
