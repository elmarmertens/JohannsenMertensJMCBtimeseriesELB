#!/bin/bash
# batchjob exec from to

# make cleanall
export OMP_NUM_THREADS=16
FCmode=prod

# launch
for datalabel in spectreTB3MSGS020510OutputGapHeadline2018Q4 spectreTB3MSGS020510Headline2018Q4
do
    echo $datalabel
    make -f irf.makefile run FCmode=$FCmode datalabel=$datalabel Nparticles=10000
done
