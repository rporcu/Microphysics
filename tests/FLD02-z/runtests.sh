#!/bin/bash -lex

RUN_NAME="FLD02"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

rm -rf POST_* ${RUN_NAME}* &> /dev/null
time -p ${MFIX} inputs

if ! [ -z "${MFIX_BENCHMARKS_HOME}" ] && ! [ -z "${FCOMPARE}" ]; then
    ${FCOMPARE} ${MFIX_BENCHMARKS_HOME}/FLD02-z_FLD02-z_plt00000 FLD0200000/
fi

if ! [ -z "${FEXTRACT}" ]; then
    ${FEXTRACT} -p FLD0200000/ -d 1 -v w_g && mv FLD0200000.slice POST_WG.dat
    ${FEXTRACT} -p FLD0200000/ -d 3 -v p_g && mv FLD0200000.slice POST_PG.dat

    post_dats=POST*.dat
    for result in ${post_dats}; do
        numdiff -a 0.0 AUTOTEST/${result} ${result}
    done
fi
