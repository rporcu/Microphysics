#!/bin/bash -lex

RUN_NAME="DEM04"

DES_IM=ADAMS_BASHFORTH
for DES_MEW in 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0; do
  rm -f ${RUN_NAME}* &> /dev/null
  time -p ./mfix DES_INTG_METHOD=\"${DES_IM}\" \
    MEW=${DES_MEW} MEW_W=${DES_MEW}
done

post_dats=AUTOTEST/POST*.dat

for test_post_file in ${post_dats}; do
    numdiff -a 0.000001 -r 0.05 ${test_post_file} $(basename ${test_post_file}) || echo "Post results differ"
done
