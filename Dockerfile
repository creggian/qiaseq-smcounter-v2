####################################################################################
########        This is a Dockerfile to describe QIAGEN's read processing  #########
########        runtime framework for spe-dna panels                       #########
####################################################################################

# Using a biocontainer base image
# Please see below for further details : 
# https://github.com/BioContainers/containers/blob/master/biocontainers/Dockerfile
FROM biocontainers/biocontainers:latest

MAINTAINER Raghavendra Padmanabhan <raghavendra.padmanabhan@qiagen.com>

################ Create appropriate directory structure for code to run ################
USER root
RUN mkdir -p /srv/qgen/code && \
    mkdir -p /srv/qgen/bin/downloads && \
    mkdir -p /srv/qgen/data/genome && \
    mkdir -p /srv/qgen/data/annotation && \
    mkdir -p /srv/qgen/example/


################ Update package repository and install dependencies using apt-get ################
RUN apt-get -y update && \
    apt-get -y install r-base

################ Install various version specific 3rd party tools ################
RUN conda install bedtools=2.25.0

################ Install python modules ################
## Install some modules with conda
RUN conda install scipy MySQL-python openpyxl pysam=0.9.0
RUN pip install statistics
    
################ R packages ################
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('plyr')"


################ Add latest samtools version for sort by Tag feature ################
RUN wget https://github.com/samtools/samtools/releases/download/1.5/samtools-1.5.tar.bz2 -O /srv/qgen/bin/downloads/samtools-1.5.tar.bz2 && \
    cd /srv/qgen/bin/downloads/ && \
    tar -xvf samtools-1.5.tar.bz2 && \
    cd samtools-1.5  && \
    mkdir -p /srv/qgen/bin/samtools-1.5 && \
    ./configure --prefix /srv/qgen/bin/samtools-1.5 && \
    make && \
    make install 

################ Add data directory ################
## Download genome files
RUN wget https://storage.googleapis.com/qiaseq-dna/data/genome/ucsc.hg19.dict \
         https://storage.googleapis.com/qiaseq-dna/data/genome/ucsc.hg19.fa.gz -P /srv/qgen/data/genome/
	 
RUN cd /srv/qgen/data/genome && \
    gunzip ucsc.hg19.fa.gz  && \
    ## Index the fasta using samtools
    /srv/qgen/bin/samtools-1.5/bin/samtools faidx /srv/qgen/data/genome/ucsc.hg19.fa
    
## Download Annotation files
RUN wget https://storage.googleapis.com/qiaseq-dna/data/annotation/bkg.error.v2.RData \
	 https://storage.googleapis.com/qiaseq-dna/data/annotation/SR_LC_SL.full.bed \
	 https://storage.googleapis.com/qiaseq-dna/data/annotation/simpleRepeat.full.bed \
	  -P /srv/qgen/data/annotation/

## Add test files for smCounterv2
RUN mkdir -p /srv/qgen/test_smcounter-v2/
RUN wget https://storage.googleapis.com/qiaseq-dna/test_files/high.confidence.variants.bed \
         https://storage.googleapis.com/qiaseq-dna/test_files/NB956-240-3-10_S1.highconfidence.bam \
	 https://storage.googleapis.com/qiaseq-dna/test_files/NB956-240-3-10_S1.highconfidence.VariantList.long.txt \
	 https://storage.googleapis.com/qiaseq-dna/test_files/NB956-240-3-10_S1.highconfidence.bam.bai \
	 -P /srv/qgen/test_smcounter-v2/

################ Update Environment Variables ################
ENV PYTHONPATH $PYTHONPATH:/opt/conda/lib/python2.7/site-packages/:/srv/qgen/code/qiaseq-dna/

