#!/bin/bash
export OMP_NUM_THREADS=8


# STDOUT_FNAME=runlog.txt 
STDOUT_FNAME=runlog_$(date "+%Y.%m.%d-%H.%M.%S").txt 
mkdir -p log
echo 'Number of threads used = '$OMP_NUM_THREADS > 'log/'$STDOUT_FNAME
./a.out >>'log/'$STDOUT_FNAME 2>&1 &
