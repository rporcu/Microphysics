#!/bin/bash -lex

RUN_NAME="DEM05"

write_data() {
    echo "   " >> $1
    echo "   Particle $2" >> $1
    echo "   " >> $1

    ${FJOIN_PAR} -f DEM05_par --end 25 --id $2 --var $3 --var $4 --var $5 --format 6 --dt 0.0001 >> $1

}

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

if [ -n "$2" ]; then
    FJOIN_PAR=$2/fjoin_par
fi

rm -f ${RUN_NAME}* const_plt* POST_* &> /dev/null

time -p ${MFIX} inputs \
  KN=1.72d7 KT_FAC="@(1.48/1.72)" KN_W=1.72d7 KT_W_FAC="@(1.48/1.72)" \
  "DES_EN_INPUT(1:3)=3*1.0 DES_EN_WALL_INPUT(1:2)=2*1.0"


for pID in {1..62}; do
    write_data POST_VEL.dat ${pID} 11 10 12
done

post_dats=AUTOTEST/POST*.dat

for test_post_file in ${post_dats}; do
    diff ${test_post_file} $(basename ${test_post_file})
done
