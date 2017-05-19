#!/bin/bash -exl

# set case directory
RUN_NAME="FLD03"

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
time -p mpirun -np 4 ${MFIX} inputs

${FEXTRACT} -p FLD0300000/ -d 1 -v w_g && mv FLD0300000.slice POST_WG_MPI.dat
${FEXTRACT} -p FLD0300000/ -d 3 -v p_g && mv FLD0300000.slice POST_PG_MPI.dat

post_dats=POST*.dat
for result in ${post_dats}; do
    numdiff -a 0.0 -r 0.0 AUTOTEST/${result} ${result}
done

if ! [ -z "${MFIX_BENCHMARKS_HOME}" ] && ! [ -z "${FCOMPARE}" ]; then
    ${FCOMPARE} --infile1 ${MFIX_BENCHMARKS_HOME}/FLD03-z_FLD03-z_plt00000 --infile2 FLD0300000/
fi
