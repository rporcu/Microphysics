#!/bin/bash -lx

RUN_NAME="DEM06"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

if [ -n "$2" ]; then
    FCOMPARE=$2/plt_compare_diff_grids
fi

MFIX_BENCHMARKS_HOME=${MFIX_BENCHMARKS_HOME:-}
FCOMPARE=${FCOMPARE:-}

rm -rf POST_* ${RUN_NAME}* &> /dev/null
time -p ${MFIX} inputs DES_ONEWAY_COUPLED=.F.

if ! [ -z "${MFIX_BENCHMARKS_HOME}" ] && ! [ -z "${FCOMPARE}" ]; then
  ${FCOMPARE} --infile1 "${MFIX_BENCHMARKS_HOME}/DEM06-z_plt00350" --infile2 DEM0600350/
fi

post_dats=POST*.dat
for result in ${post_dats}; do
    diff "AUTOTEST/${result}" "${result}"
done
