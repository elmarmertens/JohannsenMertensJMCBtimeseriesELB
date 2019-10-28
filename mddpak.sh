#!/bin/bash
# batchjob exec from to

export OMP_NUM_THREADS=4
FCmode=prod

make THIS=spectre2018mdd
make THIS=spectre2018constvarmdd
make THIS=spectre2018rbarSVmdd

for datalabel in spectreTB3MSGS020510OutputGapHeadline2018Q4
do
    for i in {1..10}
    do
	./spectre2018mdd 1 0 1 0 0 236 0 1234 $datalabel 3 2 10000 1000 100000 4 $i > mdd_$i.$datalabel.log
	./spectre2018constvarmdd 1 0 1 0 0 236 0 1234 $datalabel 3 2 10000 1000 100000 4 $i > constvarmdd_$i.$datalabel.log
	./spectre2018rbarSVmdd 1 0 1 0 0 236 0 1234 $datalabel 3 2 10000 1000 100000 4 $i > rbarSVmdd_$i.$datalabel.log
    done
done
