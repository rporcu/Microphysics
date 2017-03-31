#!/bin/bash -lex

RUN_NAME="FLD02"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

rm -rf POST_* ${RUN_NAME}* &> /dev/null
time -p ${MFIX} inputs

if ! [ -z "${MFIX_BENCHMARKS_HOME}" ] && ! [ -z "${FCOMPARE}" ]; then
    ${FCOMPARE} ${MFIX_BENCHMARKS_HOME}/FLD02-x_FLD02-x_plt00000 FLD0200000/
fi

if ! [ -z "${FEXTRACT}" ]; then
    ${FEXTRACT} -p FLD0100000/ -d 2 -v u_g && mv FLD0100000.slice POST_UG.dat
    ${FEXTRACT} -p FLD0100000/ -d 1 -v p_g && mv FLD0100000.slice POST_PG.dat

    post_dats=POST*.dat
    for result in ${post_dats}; do
        numdiff -a 0.0 AUTOTEST/${result} ${result}
    done
fi
