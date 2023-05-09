# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A package with shared database and record data fetching functions """

from typing import Optional

from antismash.common.secmet import Record


def get_assembly_id(rec: Record) -> Optional[str]:
    """Extract the NCBI assembly ID from a record."""
    for ref in rec.dbxrefs:
        if not ref.startswith('Assembly:'):
            continue
        return ref[9:]

    return None
