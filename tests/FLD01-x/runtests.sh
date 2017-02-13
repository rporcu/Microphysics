#!/bin/bash -exl

# set case directory
RUN_NAME="FLD01"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

rm -f fort* POST_* &> /dev/null

rm -f ${RUN_NAME}* &> /dev/null
time -p ${MFIX} inputs

