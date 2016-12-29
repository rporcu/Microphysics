#!/bin/bash -lex

RUN_NAME="DEM03"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

rm -f POST_* &> /dev/null

DES_IM=ADAMS_BASHFORTH
for DES_ETA in 1.0; do
  rm -f ${RUN_NAME}* &> /dev/null
  time -p ${MFIX} DES_INTG_METHOD=\"${DES_IM}\" \
    DES_EN_INPUT\(1\)=${DES_ETA} \
    DES_EN_INPUT\(2\)=${DES_ETA} \
    DES_EN_INPUT\(3\)=${DES_ETA} \
    DES_EN_WALL_INPUT\(1\)=${DES_ETA} \
    DES_EN_WALL_INPUT\(2\)=${DES_ETA}
done

