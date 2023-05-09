#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

touch imported.txt
touch deferred_imported.txt

readonly SEQDIR=$1

for infile in $(find ${SEQDIR} -name "*.json"); do
    grep -q ${infile} imported.txt && echo "Skipping ${infile}" && continue
    echo ${infile}
    ./import_json.py ${infile}
    echo ${infile} >> imported.txt
done

for infile in $(find ${SEQDIR} -name "*.json"); do
    grep -q ${infile} deferred_imported.txt && echo "Skipping ${infile}" && continue
    echo ${infile}
    ./import_deferred.py ${infile}
    echo ${infile} >> deferred_imported.txt
done
