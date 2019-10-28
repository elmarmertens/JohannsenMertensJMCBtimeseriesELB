#!/bin/bash
# batchjob exec from to

# make cleanall # comment in to clean reset directory
export OMP_NUM_THREADS=4
FCmode=prod
T=236
p=2

Nsim=10000
burnin=10000


for THIS in spectre2018 spectre2018constvar spectre2018rbarSV
do
    make FCmode=$FCmode THIS=$THIS

    # launch
    for datalabel in spectreTB3MSGS020510OutputGapHeadline2018Q4 spectreTB3MSGS020510Headline2018Q4
    do
	make run FCmode=$FCmode THIS=$THIS datalabel=$datalabel ZLBrejection=1 Nsim=$Nsim burnin=$burnin p=$p T=$T > mcmc.$THIS.$datalabel.log 
    done

done

# calling nomas version (w/o long-term bonds)
THIS=nomas2018
Nyield=0
gapOrderKey=123
make FCmode=$FCmode THIS=$THIS

    # launch
    for datalabel in spectreTB3MSGS020510OutputGapHeadline2018Q4 spectreTB3MSGS020510Headline2018Q4
    do
	make run FCmode=$FCmode THIS=$THIS datalabel=$datalabel ZLBrejection=1 Nsim=$Nsim burnin=$burnin p=$p T=$T Nyield=$Nyield gapOrderKey=$gapOrderKey > mcmc.$THIS.$datalabel.log 
    done
