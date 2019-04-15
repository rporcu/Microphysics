#!/bin/bash -lex

RUN_NAME="DEM05"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

rm -rf const_plt* POST_* &> /dev/null

rm -f ${RUN_NAME}* &> /dev/null
time -p ${MFIX} inputs \
  KN=1.72d7 KT_FAC="@(1.48/1.72)" KN_W=1.72d7 KT_W_FAC="@(1.48/1.72)" \
  "DES_EN_INPUT(1:3)=3*1.0 DES_EN_WALL_INPUT(1:2)=2*1.0"

post_dats=AUTOTEST/POST*.dat

for test_post_file in ${post_dats}; do
    diff ${test_post_file} $(basename ${test_post_file})
done
