#!/bin/bash -lx

RUN_NAME="DEM06"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

rm -rf POST_* ${RUN_NAME}* &> /dev/null
time -p ${MFIX} inputs DES_ONEWAY_COUPLED=.F.
