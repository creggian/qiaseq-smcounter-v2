#!/bin/bash

python sm_counter_v2.py --runPath /home/qiauser/test_v2/ --bamName /srv/qgen/test_smcounter-v2/NB956-240-3-10_S1.highconfidence.bam \
       --bedName /srv/qgen/test_smcounter-v2/high.confidence.variants.bed --prefix NB956-240-3-10_S1.test --nCPU 2 --minBQ 25 \
       --minMQ 50 --hpLen 8 --mismatchThr 6.0 --primerDist 2 --mtThreshold 0.8 --rpb 7.6 --primerSide 1 --minAltUMI 3 --maxAltAllele 2 \
       --refGenome /srv/qgen/data/genome/ucsc.hg19.fa --repBed /srv/qgen/data/annotation/simpleRepeat.full.bed \
       --srBed /srv/qgen/data/annotation/SR_LC_SL.full.bed

python compare_outlong.py /home/qiauser/test_v2/intermediate/NB956-240-3-10_S1.test.VariantList.long.txt \
       /srv/qgen/test_smcounter-v2/NB956-240-3-10_S1.highconfidence.VariantList.long.txt True True
