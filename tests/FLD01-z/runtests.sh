#!/bin/bash -exl

# set case directory
RUN_NAME="FLD01"

MFIX=./mfix
if [ -n "$1" ]; then
    MFIX=$1
fi

rm -f POST_* &> /dev/null

rm -rf ${RUN_NAME}* &> /dev/null
time -p ${MFIX} inputs

rm -rf ${RUN_NAME}* &> /dev/null
time -p ${MFIX} inputs
~/packages/amrex/amrex-source/Tools/Postprocessing/F_Src/fextract.Linux.gfortran.exe -p FLD0100000/ -d 2 -v u_g && mv FLD0100000.slice POST_UG.dat
