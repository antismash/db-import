# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A package with shared file data loading functions """

import bz2
from io import StringIO
import json
from typing import Any

from antismash.common.serialiser import AntismashResults

RawJson = dict[str, Any]

def read_json(filename: str) -> tuple[RawJson, AntismashResults]:
    """Read the antiSMASH json file and return both the raw dict and the AntismashResults object"""
    if filename.endswith(".bz2"):
        with bz2.open(filename, mode='rt', encoding="utf-8") as hdl:
            handle = StringIO(hdl.read())
    else:
        with open(filename, encoding="utf-8") as hdl:
            handle = StringIO(hdl.read())

    raw_data: dict[str, Any] = json.load(handle)
    handle.seek(0)
    results = AntismashResults.from_file(handle)

    return raw_data, results