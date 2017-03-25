#!/bin/bash -lex

RUN_NAME="FLD02"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

rm -f POST_* &> /dev/null

for DELP_Y in -3.0 -2.0 -1.0 0.0 1.0 2.0 3.0; do
  rm -rf ${RUN_NAME}* &> /dev/null
  time -p ${MFIX} inputs DELP_Y=${DELP_Y}
  ${FEXTRACT} -p FLD0200000/ -d 3 -v v_g
  cat FLD0200000.slice >> POST_VG.dat
  echo -e "\n\n\n" >> POST_VG.dat
done

post_dats=POST*.dat
for result in ${post_dats}; do
    numdiff -a 0.0 AUTOTEST/${result} ${result}
done
