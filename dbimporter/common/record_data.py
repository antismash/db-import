# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from antismash.common.secmet import Feature, Record, Region
from antismash.common.module_results import ModuleResults

class RecordData:
    def __init__(self, cursor, record: Record, record_id: int, assembly_id: int,
                 module_results: dict[str, ModuleResults], record_no: int):
        self.cursor = cursor
        self.record = record
        self.record_id = record_id
        assert record_id
        self.assembly_id = assembly_id
        self.module_results = module_results
        self.record_no = record_no

        self._current_region = None
        self._current_region_id = None
        self.feature_mapping: dict[Feature, int] = {}

    @property
    def current_region(self):
        assert self._current_region
        return self._current_region

    @current_region.setter
    def current_region(self, region):
        assert isinstance(region, Region)
        self._current_region = region
        self._current_region_id = self.feature_mapping[region]

    @property
    def current_region_id(self):
        assert self._current_region_id
        return self._current_region_id

    def insert(self, statement, values):
        try:
            self.cursor.execute(statement, values)
        except:
            print("failed insertion:")
            print("  statement:")
            print("   ", "\n    ".join(statement.strip().splitlines()))
            print("  values:", values)
            raise
        if "RETURNING" in statement:
            return self.cursor.fetchone()[0]
        return None
