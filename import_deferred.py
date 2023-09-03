#!/usr/bin/env python3
""" Import the sections of antiSMASH results for a record that contain crosslinks
    to other records within the database, requiring them to have been already
    committed.
"""
from argparse import ArgumentParser
import json
import os
import time
import traceback

# pylint: disable=line-too-long,missing-docstring

import antismash
from antismash.common.secmet import Record
import psycopg2
import psycopg2.extensions

from dbimporter.common.data import read_json
from dbimporter.common.record_data import RecordData
from dbimporter.common import (
    comparippson,
    getters,
    preparation,
)
from dbimporter.modules import (
    clusterblast,
)

psycopg2.extensions.register_type(psycopg2.extensions.UNICODE)
psycopg2.extensions.register_type(psycopg2.extensions.UNICODEARRAY)

DB_CONNECTION = "host='localhost' port=5432 user='postgres' password='secret' dbname='antismash'"
REPORTED_TYPES = set()

DEFAULT_AS_OPTIONS = antismash.config.build_config(["--minimal"], modules=antismash.main.get_all_modules())
DEFAULT_AS_OPTIONS.all_enabled_modules = []


class ExistingRecordError(ValueError):
    pass


class MissingAssemblyIdError(ValueError):
    pass


def main(filename, db_connection):
    """Run the import."""
    connection = psycopg2.connect(db_connection)
    connection.autocommit = False

    raw_data, results = read_json(filename)
    with connection.cursor() as cursor:
        try:
            assembly_id = getters.get_assembly_id(results.records[0])
            if not assembly_id:
                short_name, _ = os.path.splitext(os.path.basename(filename))
                id_parts = short_name.split("_")
                if id_parts[0] not in ("GCF", "GCA"):
                    raise MissingAssemblyIdError("assembly ID does begin with 'GCF'/'GCA'")
                assembly_id = "_".join(id_parts[:2])

            print("assembly_id:", assembly_id, end="\t")
            if assembly_id:
                input_basename = os.path.basename(filename)
                cursor.execute("SELECT (assembly_id) FROM antismash.filenames WHERE base_filename = %s AND assembly_id = %s", (input_basename, assembly_id))
                if cursor.fetchone() is None:
                    raise ValueError("No existing import exists to add deferred annotations to")
            record_no = 0
            for rec, module_results in zip(results.records, results.results):
                raw_record = raw_data["records"][record_no]
                record_no += 1
                preparation.prepare_record(rec, raw_record["areas"], module_results)
                load_deferred_sections(rec, module_results, cursor, assembly_id, record_no)
            connection.commit()
            print(assembly_id, "deferred changes committed", end="\t")
        except ExistingRecordError:
            connection.rollback()
        except Exception:
            connection.rollback()
            raise
    connection.close()


def load_deferred_sections(rec: Record, module_results, cur, assembly_id, record_no: int):
    """Import deferred portions of records."""
    if not rec.get_regions():
        return
    seq_id = rec.annotations['accessions'][0]

    data = RecordData(cur, rec, seq_id, assembly_id, module_results, record_no)

    for region in sorted(rec.get_regions()):
        handle_region(data, seq_id, region)

    comparippson.import_results(data)


def handle_region(data: RecordData, sequence_id, region):
    """Import deferred portions of regions."""
    assert region
    data.current_region = region

    clusterblast.import_region_results(data, region, deferred=True)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--db', default=DB_CONNECTION, help="DB connection string to use (default: %(default)s)")
    parser.add_argument('--from-filelist', action="store_true", default=False, help="Passed filename is a list of filenames")
    parser.add_argument('filenames', nargs="*")
    args = parser.parse_args()
    total_duration = 0
    total_imports = 0
    successful_imports = 0

    filenames: list[str] = args.filenames
    if args.from_filelist:
        filenames = []
        for filename in args.filenames:
            with open(filename, 'r', encoding="utf-8") as handle:
                content = handle.read()
                filenames.extend(content.splitlines())

    for filename in filenames:
        start_time = time.time()
        try:
            main(filename, args.db)
            successful_imports += 1
        except MissingAssemblyIdError as err:
            print("failed to import", filename, ":", err)
        except Exception as err:
            print("failed to import", filename, ":", err)
            traceback.print_exc()
        finally:
            end_time = time.time()
            import_duration = end_time - start_time
            print("took", round(import_duration, 2), "seconds", end="\t")
            total_imports += 1
            total_duration += import_duration
            print("average:", round(total_duration/total_imports, 2), f"for {total_imports} total imports ({successful_imports} successful)")
