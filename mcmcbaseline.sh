#!/bin/bash
# batchjob exec from to

# make cleanall # comment in to clean reset directory
export OMP_NUM_THREADS=4
FCmode=prod
T=236
p=2

Nsim=10000
burnin=10000


THIS=spectre2018
make FCmode=$FCmode THIS=$THIS

datalabel=spectreTB3MSGS020510OutputGapHeadline2018Q4 # set to spectreTB3MSGS020510Headline2018Q4 to use unemployment gap instead of output gap
make run FCmode=$FCmode THIS=$THIS datalabel=$datalabel ZLBrejection=1 Nsim=$Nsim burnin=$burnin p=$p T=$T 
 
