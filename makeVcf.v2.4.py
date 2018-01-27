#!/usr/bin/python
# vim: tabstop=9 expandtab shiftwidth=3 softtabstop=3
# convert long text output to VCF format (v2.4). Merges multi-allelic variants to a single row
# currently, only allow multi-allelic SNPs
# Chang Xu. 20OCT2017

import os
import sys
from operator import itemgetter

cutoff = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]
rangeCutoff = range(len(cutoff))
minCutoff = cutoff[0]
FORMAT = 'GT:DP:VDP:VAF:UMT:VMT:VMF:LogP'
ID = '.'
headerAll = ['CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'sUMT', 'sForUMT', 'sRevUMT', 'sVMT', 'sForVMT', 'sRevVMT', 'sVMF', 'sForVMF', 'sRevVMF', 'VDP', 'VAF', 'RefForPrimer', 'RefRevPrimer', 'primerOR', 'pLowQ', 'hqUmiEff', 'allUmiEff', 'refMeanRpb', 'altMeanRpb', 'rpbEffectSize', 'repType', 'hpInfo', 'simpleRepeatInfo', 'tandemRepeatInfo', 'DP', 'FR', 'MT', 'UFR', 'sUMT_A', 'sUMT_T', 'sUMT_G', 'sUMT_C', 'logpval', 'FILTER']

#--------------------------------------------------------------------------------------
# function to handle normal variants
#--------------------------------------------------------------------------------------
def biAllelicVar(alleles, RepRegion, outvcf):
   chrom, pos, ref, alt, typ, dp, vdp, vaf, umt, vmt, vmf, qual, fqual, fltr = alleles[0]
   for i in rangeCutoff:
      if fqual >= cutoff[i]:
         INFO = ';'.join(['SAMPLE=' + sampleName, 'TYPE=' + typ, 'RepRegion=' + RepRegion]) 
         SAMPLE = ':'.join(['0/1', dp, vdp, vaf, umt, vmt, vmf, qual])
         vcfLine = '\t'.join([chrom, pos, ID, ref, alt, qual, fltr, INFO, FORMAT, SAMPLE]) + '\n'
         outvcf[i].write(vcfLine)

#--------------------------------------------------------------------------------------
# function to handle multi-allelic variants
#--------------------------------------------------------------------------------------
def multiAllelicVar(alleles, RepRegion, outvcf):
   for i in rangeCutoff:
      tmpAlleles = [x for x in alleles if x[-2] >= cutoff[i] and x[-1] == 'PASS']
      lenTmpAlleles = len(tmpAlleles)
      if lenTmpAlleles == 0:
         pass
      elif lenTmpAlleles == 1:
         chrom, pos, ref, alt, typ, dp, vdp, vaf, umt, vmt, vmf, qual, fqual, fltr = tmpAlleles[0]
         INFO = ';'.join(['SAMPLE=' + sampleName, 'TYPE=' + typ, 'RepRegion=' + RepRegion]) 
         SAMPLE = ':'.join(['0/1', dp, vdp, vaf, umt, vmt, vmf, qual])
         vcfLine = '\t'.join([chrom, pos, ID, ref, alt, qual, fltr, INFO, FORMAT, SAMPLE]) + '\n'
         outvcf[i].write(vcfLine)
      else:
         VDPs, VAFs, VMTs, VMFs, QUALs, fQUALs, TYPEs, REFs, ALTs= [], [], [], [], [], [], [], [], []
         for allele in tmpAlleles:
            chrom, pos, ref, alt, typ, dp, vdp, vaf, umt, vmt, vmf, qual, fqual, fltr = allele
            VDPs.append(vdp)
            VAFs.append(vaf)
            VMTs.append(vmt)
            VMFs.append(vmf)
            QUALs.append(qual)
            TYPEs.append(typ)
            REFs.append(ref)
            ALTs.append(alt)
            fQUALs.append(fqual)

         # align multiple alleles to the same REF if necessary
         if all(x==REFs[0] for x in REFs):
            finalRef = REFs[0]
            finalAlt = ','.join(ALTs)
         else:
            # Assumption: the first bases are the same 
            finalRef = max(REFs, key=len)
            for j in range(len(ALTs)):
               ALTs[j] = ALTs[j] if REFs[j] == finalRef else ALTs[j] + finalRef[len(REFs[j]):]
            finalAlt = ','.join(ALTs)

         newQual = str(min(fQUALs))
         allTypes = ','.join(TYPEs)
         allVDPs = ','.join(VDPs)
         allVAFs = ','.join(VAFs)
         allVMTs = ','.join(VMTs)
         allVMFs = ','.join(VMFs)
         allQUALs = ','.join(QUALs)

         INFO = ';'.join(['SAMPLE=' + sampleName, 'TYPE=' + allTypes, 'RepRegion=' + RepRegion]) 
         SAMPLE = ':'.join(['0/1', dp, allVDPs, allVAFs, umt, allVMTs, allVMFs, allQUALs])
         vcfLine = '\t'.join([chrom, pos, ID, finalRef, finalAlt, newQual, 'PASS', INFO, FORMAT, SAMPLE]) + '\n'
         outvcf[i].write(vcfLine)

#--------------------------------------------------------------------------------------
# main function
#--------------------------------------------------------------------------------------
def main(runPath, outlong, sampleName):
   # change working directory to runDir and make output directories 
   os.chdir(runPath)
   if not os.path.exists('vcf'):
      os.makedirs('vcf')

   # create a vcf file at each cutoff
   base = os.path.basename(outlong)
   header_vcf = '##fileformat=VCFv4.2' + '\n' + \
         '##FILTER=<ID=LM,Description="Low coverage (fewer than 5 barcodes)">' + '\n' + \
         '##FILTER=<ID=RepT,Description="Variant in tandem repeat (TFR) regions">' + '\n' + \
         '##FILTER=<ID=RepS,Description="Variant in simple repeats (RepeatMasker) regions">' + '\n' + \
         '##FILTER=<ID=HP,Description="Inside or flanked by homopolymer regions">' + '\n' + \
         '##FILTER=<ID=LowC,Description="Variant in Low complexity regions, as defined in RepeatMasker">' + '\n' + \
         '##FILTER=<ID=SL,Description="Variant in micro-satelite regions, as defined in RepeatMasker">' + '\n' + \
         '##FILTER=<ID=SB,Description="Strand Bias">' + '\n' + \
         '##FILTER=<ID=DP,Description="Too many discordant pairs">' + '\n' + \
         '##FILTER=<ID=MM,Description="Too many mismatches in a read. Default threshold is 6.5 per 100 bases">' + '\n' + \
         '##FILTER=<ID=LowQ,Description="Low base quality">' + '\n' + \
         '##FILTER=<ID=RBCP,Description="Variant are clustered at the end of barcode-side reads">' + '\n' + \
         '##FILTER=<ID=RPCP,Description="Variant are clustered at the end of primer-side reads">' + '\n' + \
         '##FILTER=<ID=PB,Description="Primer bias filter. odds ratio > 10 or < 0.1">' + '\n' + \
         '##FILTER=<ID=PrimerCP,Description="variant is clustered within 2 bases from primer sequence due to possible primer dimers">' + '\n' + \
         '##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample name">' + '\n' + \
         '##INFO=<ID=TYPE,Number=.,Type=String,Description="Variant type: SNP/INDEL/COMPLEX">' + '\n' + \
         '##INFO=<ID=RepRegion,Number=.,Type=String,Description="Repetitive region">' + '\n' + \
         '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' + '\n' + \
         '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total read depth">' + '\n' + \
         '##FORMAT=<ID=VDP,Number=.,Type=Integer,Description="Variant read depth">' + '\n' + \
         '##FORMAT=<ID=VAF,Number=.,Type=Float,Description="Variant read allele frequency">' + '\n' + \
         '##FORMAT=<ID=UMT,Number=1,Type=Integer,Description="Total used UMI depth">' + '\n' + \
         '##FORMAT=<ID=VMT,Number=.,Type=Integer,Description="Variant UMI depth">' + '\n' + \
         '##FORMAT=<ID=VMF,Number=.,Type=Float,Description="Variant UMI allele frequency">' + '\n' + \
         '##FORMAT=<ID=LogP,Number=.,Type=Float,Description="Log P value. For bi- and multi-allelic variants, QUAL is the minimum of LogP">' + '\n' + \
         '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + [sampleName]) + '\n'

   outvcf, alleles = [], []
   lastCHROM, lastPOS = '', ''

   for i in rangeCutoff:
      vcfName = 'vcf/' + sampleName + '.' + str(cutoff[i]) + '.vcf'
      outvcf.append(open(vcfName, 'w'))
      outvcf[i].write(header_vcf)

   cnt = 1
   with open(outlong, 'r') as f:
      next(f)
      for line in f:
         cnt += 1
         CHROM, POS, REF, ALT, TYPE, sUMT, sForUMT, sRevUMT, sVMT, sForVMT, sRevVMT, sVMF, sForVMF, sRevVMF, VDP, VAF, RefForPrimer, RefRevPrimer, primerOR, pLowQ, hqUmiEff, allUmiEff, refMeanRpb, altMeanRpb, rpbEffectSize, repType, hpInfo, simpleRepeatInfo, tandemRepeatInfo, DP, FR, MT, UFR, sUMT_A, sUMT_T, sUMT_G, sUMT_C, logpval, FILTER = line.strip().split('\t')

         if ALT == 'DEL': 
            continue

         QUAL = logpval if logpval != 'NA' else '0.00'

         try:
            fQUAL = float(QUAL)
         except ValueError:
            fQUAL = 0.00

         if fQUAL < minCutoff:
            lastCHROM = '.'
            continue
         try:
            VAF = str(float(VAF)/100)
         except ValueError:
            VAF = '-1'
         try:
            sVMF = str(float(sVMF)/100)
         except ValueError:
            sVMF = '-1'

         # rep types are separeted by ";" in the long output. Replace to "," to comply with VCF format
         RepRegion = repType.replace(';', ',')

         currentAllele = (CHROM, POS, REF, ALT, TYPE, DP, VDP, VAF, sUMT, sVMT, sVMF, QUAL, fQUAL, FILTER)
         lenAlleles = len(alleles)

         # if current chrom and position equal to last line, append it for potential multi-allelic output
         if lenAlleles == 0 or (CHROM == lastCHROM and POS == lastPOS):
            alleles.append(currentAllele)

         # for new chrom or position, if last variant is not multi-allelic, write to vcf directly
         elif lenAlleles == 1:
            biAllelicVar(alleles, RepRegion, outvcf)
            alleles = [currentAllele]

         # if last variant is possible multi-allelic, combine and write as one 
         else: 
            multiAllelicVar(alleles, RepRegion, outvcf)
            alleles = [currentAllele]

         lastCHROM, lastPOS = CHROM, POS

   # take care of the last line
   lenAlleles = len(alleles)
   if lenAlleles == 1:
      biAllelicVar(alleles, RepRegion, outvcf)
   elif lenAlleles >= 2:
      multiAllelicVar(alleles, RepRegion, outvcf)
   else:
      pass

   # close all output vcf
   for vcf in outvcf:
      vcf.close()

#----------------------------------------------------------------------------------------------
#pythonism to run from the command line
#----------------------------------------------------------------------------------------------
if __name__ == "__main__":
   runPath = sys.argv[1]
   outlong = sys.argv[2]
   sampleName = sys.argv[3]
   main(runPath, outlong, sampleName)





