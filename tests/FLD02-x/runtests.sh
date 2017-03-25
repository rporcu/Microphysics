#!/bin/bash -lex

RUN_NAME="FLD02"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

rm -f POST_* &> /dev/null

for DELP_X in -3.0 -2.0 -1.0 0.0 1.0 2.0 3.0; do
  rm -rf ${RUN_NAME}* &> /dev/null
  time -p ${MFIX} inputs DELP_X=${DELP_X}
  ${FEXTRACT} -p FLD0200000/ -d 2 -v u_g
  cat FLD0200000.slice >> POST_UG.dat
  echo -e "\n\n\n" >> POST_UG.dat
done

post_dats=POST*.dat
for result in ${post_dats}; do
    numdiff -a 0.0 AUTOTEST/${result} ${result}
done
