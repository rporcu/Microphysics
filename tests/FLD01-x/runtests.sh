#!/bin/bash -exl

# set case directory
RUN_NAME="FLD01"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

rm -f POST_* &> /dev/null
rm -rf ${RUN_NAME}* &> /dev/null
time -p ${MFIX} inputs

${FEXTRACT} -p FLD0100000/ -d 2 -v u_g && mv FLD0100000.slice POST_UG.dat
${FEXTRACT} -p FLD0100000/ -d 1 -v p_g && mv FLD0100000.slice POST_PG.dat

post_dats=POST*.dat
for result in ${post_dats}; do
    numdiff -a 0.0 AUTOTEST/${result} ${result}
done
