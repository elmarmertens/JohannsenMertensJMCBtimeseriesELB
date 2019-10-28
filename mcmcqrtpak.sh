#!/bin/bash
# batchjob exec from to

export OMP_NUM_THREADS=4
FCmode=prod
p=2

Nsim=10000
burnin=10000


for THIS in spectre2018 spectre2018constvar spectre2018rbarSV
do
    make FCmode=$FCmode THIS=$THIS

    # launch
    for datalabel in spectreTB3MSGS020510OutputGapHeadline2018Q4 spectreTB3MSGS020510Headline2018Q4
    do
        for T in {100..236}
        do
	    make run FCmode=$FCmode THIS=$THIS datalabel=$datalabel ZLBrejection=1 Nsim=$Nsim burnin=$burnin p=$p T=$T > mcmc.$THIS.$datalabel.$T.log
        done     
    done

done

# calling nomas version (w/o long-term bonds)
THIS=nomas2018
Nyield=0
gapOrderKey=123
make FCmode=$FCmode THIS=$THIS

for datalabel in spectreTB3MSGS020510OutputGapHeadline2018Q4 spectreTB3MSGS020510Headline2018Q4
do
    for T in {100..236}
    do
	make run FCmode=$FCmode THIS=$THIS datalabel=$datalabel ZLBrejection=1 Nsim=$Nsim burnin=$burnin p=$p T=$T Nyield=$Nyield gapOrderKey=$gapOrderKey > mcmc.$THIS.$datalabel.$T.log
    done
done

# using data for camparison with Wu-Xia
THIS=spectre2018
for datalabel in spectreTB3MSSVEN020510EOQOutputGapHeadline2018Q4 spectreTB3MSSVEN020510EOQHeadline2018Q4
do
    for T in {60..116}
    do
	make run FCmode=$FCmode THIS=$THIS datalabel=$datalabel ZLBrejection=1 Nsim=$Nsim burnin=$burnin p=$p T=$T T0=120 > mcmc.$THIS.$datalabel.$T.$T0log
    done     
done
