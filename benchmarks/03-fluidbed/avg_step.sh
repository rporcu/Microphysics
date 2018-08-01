#/bin/bash -l

FILE=$1

grep "Time per step" $FILE | sed 's/Time per step//g' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }'
