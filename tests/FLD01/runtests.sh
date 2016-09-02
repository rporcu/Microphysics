#!/bin/bash -exl

# set case directory
RUN_NAME="FLD01"

for GMAX in 16 32 64 128; do
   rm -f ${RUN_NAME}* &> /dev/null
   time -p ./mfix IMAX=${GMAX} JMAX=${GMAX}
done
