#!/bin/bash -exl

# set case directory
RUN_NAME="FLD01"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

if [ -n "$2" ]; then
    FCOMPARE=$2/plt_compare_diff_grids
    FEXTRACT=$2/fextract
fi
if [ -z "${FEXTRACT}" ]; then
    echo "FEXTRACT is not set. Aborting test."
    exit 1
fi

rm -rf POST_* ${RUN_NAME}* &> /dev/null
time -p ${MFIX} inputs

${FEXTRACT} -p FLD0100000/ -d 1 -v w_g -s POST_WG.dat
${FEXTRACT} -p FLD0100000/ -d 3 -v p_g -s POST_PG.dat

post_dats=POST*.dat
for result in ${post_dats}; do
  numdiff -a 0.0 AUTOTEST/${result} ${result}
done

if ! [ -z "${MFIX_BENCHMARKS_HOME}" ] && ! [ -z "${FCOMPARE}" ]; then
    ${FCOMPARE} --infile1 ${MFIX_BENCHMARKS_HOME}/FLD01-z_FLD01-z_plt00000 --infile2 FLD0100000/
fi
