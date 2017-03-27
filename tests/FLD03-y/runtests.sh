#!/bin/bash -exl

# set case directory
RUN_NAME="FLD03"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

rm -f POST_* &> /dev/null
rm -rf ${RUN_NAME}* &> /dev/null
time -p ${MFIX} inputs

${FEXTRACT} -p FLD0300000/ -d 3 -v v_g && mv FLD0300000.slice POST_VG.dat
${FEXTRACT} -p FLD0300000/ -d 2 -v p_g && mv FLD0300000.slice POST_PG.dat
