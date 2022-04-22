#!/bin/sh
#set -o noclobber
FILE=status.txt

#delete the file if it exists
if [ -f "$FILE" ]; then
    rm status.txt
fi

#iterate and run the optimization 1 step at a time
for i in {1..10}
do
    echo "Global iteration $i" >> status.txt
    mpiexec_mpt -n 192 python optimization.py 2>&1 > output.txt
done