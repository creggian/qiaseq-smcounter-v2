#!/usr/bin/Rscript
# vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
# integrate prior info of bkg error distribution (based on duplex data) and data specific error distribution
# modified the prior distribution using more stringint consensus criteria. 
# different prior distributions for including/excluding 1 read pair UMIs
# Chang Xu, 16MAY2017

rm(list=ls())
library(plyr)
#bkgErrorDistSimulation <- '/qgen/home/xuc/frequentlyUsedFiles/bkg.error.v2.RData'
#bkgErrorDistSimulation <- '/mnt/webserver/datadisk/varcall/frequentlyUsedFiles/bkg.error.v2.RData'
bkgErrorDistSimulation <- '/home/xuc/frequentlyUsedFiles/bkg.error.v2.RData'

args <- commandArgs(TRUE)
wd <- args[1]
outlong <- args[2]
bkgfile <- args[3]
seed <- as.numeric(args[4])
nsim <- as.numeric(args[5])
outfile <- args[6]
rpb <- as.numeric(args[7])
minAltUMI <- as.numeric(args[8])
min.mtDepth <- 1000

# set working directory
setwd(wd)
set.seed(seed)

# function to calculate standard deviation
beta.sd <- function(a,b) sqrt(a*b) / ((a+b) * sqrt(a+b+1))
# function to estimate a
calc.a <- function(mu, sigma) mu * (mu*(1-mu) / sigma^2 - 1)
# function to estimate b
calc.b <- function(mu, sigma) (1-mu) * (mu*(1-mu) / sigma^2 - 1)
# function to compute p values
calc.pval <- function(TYPE, REF, ALT, sForUMT, sRevUMT, sForVMT, sRevVMT, p.high.final, p.low.final){
  totalN <- sForUMT + sRevUMT
  totalX <- sForVMT + sRevVMT
  if(totalX >= minAltUMI){    
    if(TYPE=='INDEL' | (REF=='A' & ALT=='G') | (REF=='G' & ALT=='A') | (REF=='C' & ALT=='T') | (REF=='T' & ALT=='C')){
      pr <- p.high.final
    } else{
      pr <- p.low.final
    }
    tmp <- pbinom(q=totalX-1, size=totalN, prob=pr, lower.tail=F)
    pval <- ifelse(totalN==0, 1, mean(tmp, na.rm=T))
  } else{
    pval <- 1.0 
  }
  return(pval)
}

# define constants
cols <- c('chrom', 'pos', 'ref', 'AG', 'GA', 'CT', 'TC', 'AC', 'AT', 'CA', 'CG', 'GC', 'GT', 'TA', 'TG', 'neg.strand', 'pos.strand', 'all.smt')
out <- NULL

# read in smCounter output 
dat <- read.delim(outlong, header=T, stringsAsFactors=F)

# read in prior information
load(bkgErrorDistSimulation)
if(rpb >= 3.0){
  top4 <- bkg.error$top4.exclude.1rpUMI
} else{
  top4 <- bkg.error$top4.include.1rpUMI
}

a.ga.orig <- top4$shape1[2]
b.ga.orig <- top4$shape2[2]
sigma.high <- beta.sd(a.ga.orig, b.ga.orig)

a.ct.orig <- top4$shape1[3]
b.ct.orig <- top4$shape2[3]

# proportion of zeros
p0.high <- top4$p0[2]
p0.low <- top4$p0[3]
n0.high <- floor(nsim * p0.high)
n0.low <- floor(nsim * p0.low)

# read in data-specific background error file
bkg <- read.delim(bkgfile, header=T, stringsAsFactors=F, sep='\t')
colnames(bkg) <- cols

################### bkg errors from the readset ##################
# A/G error rate from data 
tmp <- bkg[(bkg$ref=='A' & bkg$neg.strand > min.mtDepth) | (bkg$ref=='T' & bkg$pos.strand > min.mtDepth), ]
tmp$all <- ifelse(tmp$ref=='A', tmp$neg.strand, tmp$pos.strand)
d.ag <- tmp$AG / tmp$all
tmp <- tmp[d.ag < 0.01,]
mean.ag <- sum(tmp$AG) / sum(tmp$all)
n.ag <- nrow(tmp)

# G/A error rate from data
tmp <- bkg[(bkg$ref=='G' & bkg$neg.strand > min.mtDepth) | (bkg$ref=='C' & bkg$pos.strand > min.mtDepth), ]
tmp$all <- ifelse(tmp$ref=='G', tmp$neg.strand, tmp$pos.strand)
d.ga <- tmp$GA / tmp$all
tmp <- tmp[d.ga < 0.01,]
mean.ga <- sum(tmp$GA) / sum(tmp$all)
n.ga <- nrow(tmp)

# highest error rate
mu.high <- max(mean.ag, mean.ga)
n.high <- min(n.ag, n.ga)

if(is.na(mu.high) | is.na(n.high) | n.high < 100) {
   p.high <- rbeta(n=nsim, shape1=a.ga.orig, shape2=b.ga.orig)
} else{
   a.high <- calc.a(mu.high, sigma.high)
   b.high <- calc.b(mu.high, sigma.high)
   p.high <- rbeta(n=nsim, shape1=a.high, shape2=b.high)
}

p.low <- c(rbeta(n=nsim-n0.low, shape1=a.ct.orig, shape2=b.ct.orig), rep(0, n0.low))

# compute p-values
dat$sForUMT <- as.numeric(dat$sForUMT)
dat$sRevUMT <- as.numeric(dat$sRevUMT)
dat$sForVMT <- as.numeric(dat$sForVMT)
dat$sRevVMT <- as.numeric(dat$sRevVMT)

tmp <- subset(dat, select=c(TYPE, REF, ALT, sForUMT, sForVMT, sRevUMT, sRevVMT))
pval <- mdply(tmp, calc.pval, p.high.final=p.high, p.low.final=p.low)
# set mininum at 1e-200 to avoid log(0)
raw.pval <- pmax(1e-200, pval$V1)
# take -log10
dat$logpval <- round(-log10(raw.pval), 2)

# save to disk
dat <- subset(dat, select=c(CHROM, POS, REF, ALT, TYPE, sUMT, sForUMT, sRevUMT, sVMT, sForVMT, sRevVMT, sVMF, sForVMF, sRevVMF, VDP, VAF, RefForPrimer, RefRevPrimer, primerOR, pLowQ, hqUmiEff, allUmiEff, refMeanRpb, altMeanRpb, rpbEffectSize, repType, hpInfo, simpleRepeatInfo, tandemRepeatInfo, DP, FR, MT, UFR, sUMT_A, sUMT_T, sUMT_G, sUMT_C, logpval, FILTER))
write.table(dat, outfile, sep='\t', row.names=F, col.names=T, quote=F)


