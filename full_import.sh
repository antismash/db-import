#!/bin/bash

set -o errexit
set -o pipefail
set -o nounset

touch imported.txt
touch deferred_imported.txt

ERROR_FILE=import_errors.txt

rm -f $ERROR_FILE

IMPORTDIR=`dirname "$0"`
readonly SEQDIR=$1
readonly TAXONOMY=$2

echo "Importing base results"
for infile in $(find ${SEQDIR} -name "*.json"); do
    grep -q ${infile} imported.txt && echo "Skipping ${infile}" && continue
    echo "importing ${infile}"
    if $IMPORTDIR/import_json.py --taxonomy ${TAXONOMY} ${infile}; then
        echo ${infile} >> imported.txt
    else
        echo ${infile} >> $ERROR_FILE
        false  # marks the loop as a failure, but doesn't _exit_ the loop
    fi
done || { echo "Skipping deferred imports due to at least one base result import failure (see '$ERROR_FILE')"; exit 1;}


echo "Importing deferred results"
for infile in $(find ${SEQDIR} -name "*.json"); do
    grep -q ${infile} deferred_imported.txt && echo "Skipping ${infile}" && continue
    echo "importing deferred portions of ${infile}"
    if $IMPORTDIR/import_deferred.py ${infile}; then
        echo ${infile} >> deferred_imported.txt
    else
        echo "deferred" ${infile} >> import_errors.txt
        false
    fi
done || { echo "Not all deferred imports successful (see '$ERROR_FILE')"; exit 1;}
