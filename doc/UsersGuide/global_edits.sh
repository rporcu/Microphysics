#!/bin/bash

# This script with transverse all the tex files in the diretory and replace 
# the text string A with the text string B with sed. The command for doing so 
# with sed is:
#               sed -i 's/A/B/g' <filename> 
#


for d in $(find . -maxdepth 2 -name "*.tex" )
do
   echo $d
   #sed -i 's|\\AMReX|\\amrex|g' $d
   #sed -i 's|\\nyx|\\amrex|g' $d
   #sed -i 's/[[:space:]]amrex/ AMReX/g' $d
   #sed -i 's/[[:space:]]Nyx/ MFIX-Exa/g' $d
   sed -i 's|/Nyx|/mfix|g' $d
   
done 
