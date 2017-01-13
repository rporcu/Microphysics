#!/bin/bash -exl

# set case directory
RUN_NAME="FLD01"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

rm -f POST_* &> /dev/null

for GMAX in 8 16 32; do
   rm -f ${RUN_NAME}* &> /dev/null
   time -p ${MFIX} inputs IMAX=${GMAX} JMAX=${GMAX} 
done

post_dats=AUTOTEST/POST*.dat

for test_post_file in ${post_dats}; do
    numdiff -a 0.000001 -r 0.05 ${test_post_file} \
      $(basename ${test_post_file})
done
