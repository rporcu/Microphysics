#!/bin/bash -lex

RUN_NAME="DEM01"

MFIX=./mfix
if [ -n "$1" ]; then
   MFIX=$1
fi

rm -f POST_* &> /dev/null

DES_KN=10000
for DES_ETA in 0.9 0.8 0.7 0.6; do
  rm -f ${RUN_NAME}* &> /dev/null
  time -p ${MFIX} inputs \
    DES_EN_INPUT=${DES_ETA} DES_EN_WALL_INPUT=${DES_ETA} \
    KN=${DES_KN} KN_W=${DES_KN}
done

for DES_KN in 25000 50000 100000; do
  for DES_ETA in 1.0 0.9 0.8 0.7 0.6; do
    rm -f ${RUN_NAME}* &> /dev/null
    time -p ${MFIX} inputs \
      DES_EN_INPUT=${DES_ETA} DES_EN_WALL_INPUT=${DES_ETA} \
      KN=${DES_KN} KN_W=${DES_KN}
  done
done

post_dats=AUTOTEST/POST*.dat

for test_post_file in ${post_dats}; do
    numdiff -a 0.000001 -r 0.05 ${test_post_file} $(basename ${test_post_file})
done
