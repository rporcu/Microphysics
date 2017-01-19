#!/bin/bash -exl

# set case directory
RUN_NAME="FLD03"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

rm -f POST_* &> /dev/null

for GMAX in 8; do
   rm -f ${RUN_NAME}* &> /dev/null
   time -p ${MFIX} inputs IMAX=${GMAX} JMAX=${GMAX} 
done
