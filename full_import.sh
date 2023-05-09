#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

touch imported.txt

readonly SEQDIR=$1

for infile in $(find ${SEQDIR} -name "*.json"); do
    grep -q ${infile} imported.txt && echo "Skipping ${infile}" && continue
    echo ${infile}
    ./import_json.py ${infile}
    echo ${infile} >> imported.txt
done
