#!/bin/bash -lx

RUN_NAME="DEM06"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

rm -rf POST_* ${RUN_NAME}* &> /dev/null
time -p ${MFIX} inputs DES_ONEWAY_COUPLED=.F.

if ! [ -z "${MFIX_BENCHMARKS_HOME}" ] && ! [ -z "${FCOMPARE}" ]; then
  ${FCOMPARE} ${MFIX_BENCHMARKS_HOME}/DEM06-y_plt00350 DEM0600350/
fi

post_dats=POST*.dat
for result in ${post_dats}; do
    numdiff -a 0.0 AUTOTEST/${result} ${result}
done
