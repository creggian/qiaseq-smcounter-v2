# qiaseq-smcounter-v2
```
usage: smCounter.v2.4.py [-h] [--runPath RUNPATH] [--bedName BEDNAME]
                         [--bamName BAMNAME] [--prefix PREFIX] [--nCPU NCPU]
                         [--minBQ MINBQ] [--minMQ MINMQ] [--hpLen HPLEN]
                         [--mismatchThr MISMATCHTHR] [--primerDist PRIMERDIST]
                         [--mtThreshold MTTHRESHOLD] [--rpb RPB] [--isRna]
                         [--primerSide PRIMERSIDE] [--minAltUMI MINALTUMI]
                         [--maxAltAllele MAXALTALLELE] [--refGenome REFGENOME]
                         [--repBed REPBED] [--srBed SRBED]

Variant calling using molecular barcodes

optional arguments:
  -h, --help            show this help message and exit
  --runPath RUNPATH     path to working directory
  --bedName BEDNAME     BED file
  --bamName BAMNAME     BAM file
  --prefix PREFIX       file name prefix
  --nCPU NCPU           number of CPU to use in parallel
  --minBQ MINBQ         minimum base quality allowed for analysis
  --minMQ MINMQ         minimum mapping quality allowed for analysis
  --hpLen HPLEN         Minimum length for homopolymers
  --mismatchThr MISMATCHTHR
                        average number of mismatches per 100 bases allowed
  --primerDist PRIMERDIST
                        filter variants that are within X bases to primer
  --mtThreshold MTTHRESHOLD
                        threshold on read proportion to determine MT level
                        consensus
  --rpb RPB             mean read pairs per UMI; default at 0 and will be
                        calculated
  --isRna               RNAseq varinat calling only; default is DNAseq
  --primerSide PRIMERSIDE
                        read end that includes the primer; default is 1
  --minAltUMI MINALTUMI
                        minimum requirement of ALT UMIs; default is 3
  --maxAltAllele MAXALTALLELE
                        maximum ALT alleles that meet minAltUMI to be
                        reported; default is 2 (tri-allelic variants)
  --refGenome REFGENOME
                        Path to the reference fasta file
  --repBed REPBED       Path to the simpleRepeat bed file
  --srBed SRBED         Path to the full repeat bed file
```

[![Build Status](https://travis-ci.org/qiaseq/qiaseq-smcounter-v2.svg?branch=master)](https://travis-ci.org/qiaseq/qiaseq-smcounter-v2)
