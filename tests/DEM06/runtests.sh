#!/bin/bash -lx

RUN_NAME="DEM06"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

rm -f POST_* &> /dev/null

rm -f ${RUN_NAME}* &> /dev/null
time -p ${MFIX} DES_ONEWAY_COUPLED=.F. \
    DES_INTERP_ON=.F. DES_INTERP_MEAN_FIELDS=.F.

rm -f ${RUN_NAME}* &> /dev/null
time -p ${MFIX} DES_ONEWAY_COUPLED=.F. \
    DES_INTERP_ON=.T. DES_INTERP_MEAN_FIELDS=.T. \
    DES_INTERP_SCHEME=\'GARG_2012\'

post_dats=AUTOTEST/POST*.dat

for test_post_file in ${post_dats}; do
    numdiff -a 0.000001 -r 0.05 ${test_post_file} $(basename ${test_post_file})
done
