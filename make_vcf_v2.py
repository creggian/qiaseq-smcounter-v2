#!/usr/bin/python
# vim: tabstop=9 expandtab shiftwidth=3 softtabstop=3
# convert long text output to VCF format (v2.4). Merges multi-allelic variants to a single row
# currently, only allow multi-allelic SNPs
# Chang Xu. 20OCT2017

import os
import sys
from operator import itemgetter


def assign_gt(alt,chrom,vmf):
   ''' Function for faking the Genotype i.e. GT field
   for downstream tools
   :param alts (str) alternative allele(s)
   :param chrom (str) chromosome the variant is on
   :param vmf (str) variant minor allele frequency (comma seperated for multi-allelic sites)
   '''
   alts = alt.split(",")
   if len(alts) >= 2: ## Treat all multiallelic sites as heterozygotes for the first 2 variant alleles
      genotype = '1/2'
   elif chrom == "chrY" or chrom == "chrM":
      genotype = '1'
   elif float(vmf) > 0.95 : ## Treat as Heterozygous
      genotype = '1/1'
   else:
      genotype = '0/1'

   return genotype

def assign_ad(uumi,vumi):
   ''' Function for faking the Allele Depth i.e. AD field
   for downstream tools
   :param uumi (str) total umis at the variant site
   :param vumi (str) umis corresponding to the non-reference allele(s) at the variant site (comma seperated for multi-allelic sites)
   '''    
   vumis = vumi.split(',')   
   refumi = int(uumi)
   for umi in vumis:
      refumi = refumi - int(umi)
   refumi = str(refumi)
   ad = refumi + ',' + ','.join(vumis)
   return ad 
      
#--------------------------------------------------------------------------------------
# function to handle normal variants
#--------------------------------------------------------------------------------------
def biAllelicVar(alleles, RepRegion, outVcf, outVariants):
   ID = '.'
   chrom, pos, ref, alt, typ, dp, vdp, vaf, umt, vmt, vmf, qual, fqual, fltr = alleles[0]

   INFO = ';'.join(
      ['TYPE=' +typ,'RepRegion=' + RepRegion,'DP='+dp,'UMT='+umt,'VMT='+vmt,
      'VMF='+vmf]
      ) 
   FORMAT = 'GT:AD:VF'
   gt = assign_gt(alt,chrom,vmf)      
   ad = assign_ad(umt,vmt)         
   SAMPLE = ':'.join([gt,ad,vmf])
   vcfLine = '\t'.join([chrom, pos, ID, ref, alt, qual, fltr, INFO, FORMAT, SAMPLE]) + '\n'
   outVcf.write(vcfLine)
   cutVarLine = '\t'.join([chrom, pos, ref, alt, typ, dp, vdp, vaf, umt, vmt, vmf, qual, fltr]) + '\n'
   outVariants.write(cutVarLine)

#--------------------------------------------------------------------------------------
# function to handle multi-allelic variants
#--------------------------------------------------------------------------------------
def multiAllelicVar(alleles, RepRegion, outVcf, outVariants):
  ID = '.'
  tmpAlleles = [x for x in alleles if x[-1] == 'PASS']
  lenTmpAlleles = len(tmpAlleles)
  if lenTmpAlleles == 0:
     pass
  elif lenTmpAlleles == 1:
     chrom, pos, ref, alt, typ, dp, vdp, vaf, umt, vmt, vmf, qual, fqual, fltr = tmpAlleles[0]
     INFO = ';'.join(
            ['TYPE=' +typ,'RepRegion=' + RepRegion,'DP='+dp,'UMT='+umt,'VMT='+vmt,
             'VMF='+vmf]
         )
     FORMAT = 'GT:AD:VF'
     gt = assign_gt(alt,chrom,vmf)
     ad = assign_ad(umt,vmt)
     SAMPLE = ':'.join([gt,ad,vmf])     
     vcfLine = '\t'.join([chrom, pos, ID, ref, alt, qual, fltr, INFO, FORMAT, SAMPLE]) + '\n'
     outVcf.write(vcfLine)
     cutVarLine = '\t'.join([chrom, pos, ref, alt, typ, dp, vdp, vaf, umt, vmt, vmf, qual, fltr]) + '\n'
     outVariants.write(cutVarLine)
  else:
     VDPs, VAFs, VMTs, UMTs, VMFs, QUALs, fQUALs, TYPEs, REFs, ALTs, DPs = [], [], [], [], [], [], [], [], [], [], []
     for allele in tmpAlleles:
        chrom, pos, ref, alt, typ, dp, vdp, vaf, umt, vmt, vmf, qual, fqual, fltr = allele
        VDPs.append(vdp)
        VAFs.append(vaf)
        VMTs.append(vmt)
        UMTs.append(umt)
        VMFs.append(vmf)
        QUALs.append(qual)
        TYPEs.append(typ)
        REFs.append(ref)
        ALTs.append(alt)
        fQUALs.append(fqual)
        DPs.append(dp)

     # debug check
     assert len(set(UMTs)) == 1, "The number of used UMIs at a site should be the same across all alleles"

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
     allDPs = ','.join(DPs)

     INFO = ';'.join(
            ['TYPE=' +allTypes,'RepRegion=' + RepRegion,'DP='+allDPs,'UMT='+umt,'VMT='+allVMTs,
             'VMF='+allVMFs]
         )          
     FORMAT = 'GT:AD:VF' 
     gt = assign_gt(finalAlt,chrom,allVMFs)
     ad = assign_ad(umt,allVMTs)   
     SAMPLE = ':'.join([gt,ad,allVMFs])
     vcfLine = '\t'.join([chrom, pos, ID, finalRef, finalAlt, newQual, 'PASS', INFO, FORMAT, SAMPLE]) + '\n'     
     outVcf.write(vcfLine)
     cutVarLine = '\t'.join([chrom, pos, finalRef, finalAlt, allTypes, allDPs, allVDPs, allVAFs, umt, allVMTs, allVMFs, newQual,'PASS']) + '\n'
     outVariants.write(cutVarLine)

#--------------------------------------------------------------------------------------
# main function
#--------------------------------------------------------------------------------------
def main(runPath, outlong, sampleName):
   # change working directory to runDir
   os.chdir(runPath)
   outAll = open(sampleName + '.smCounter.all.txt', 'w')
   outVariants = open(sampleName + '.smCounter.cut.txt','w')
   outVcf = open(sampleName + '.smCounter.cut.vcf','w')
   outLowPi = open(sampleName + '.smCounter.lowQ.txt','w')
   
   cutoff = 6
   minCutoff = {'INDEL': 2,'SNP':2} ## Cutoff for the low-PI file

   
   ID = '.'
   headerAll = ['CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'sUMT', 'sForUMT', 'sRevUMT', 'sVMT', 'sForVMT', 'sRevVMT', 'sVMF', 'sForVMF', 'sRevVMF', 'VDP', 'VAF', 'RefForPrimer', 'RefRevPrimer', 'primerOR', 'pLowQ', 'hqUmiEff', 'allUmiEff', 'refMeanRpb', 'altMeanRpb', 'rpbEffectSize', 'repType', 'hpInfo', 'simpleRepeatInfo', 'tandemRepeatInfo', 'DP', 'FR', 'MT', 'UFR', 'sUMT_A', 'sUMT_T', 'sUMT_G', 'sUMT_C', 'logpval', 'FILTER']
   headerVariants = ['CHROM','POS','REF','ALT','TYPE','DP','VDP','VAF','sUMT','sVMT','sVMF','QUAL','FILTER']
   headerLowPi = [sampleName] + headerVariants   
   headerVcf = '##fileformat=VCFv4.2' + '\n' + \
         '##reference=GRCh37' + '\n' + \
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
         '##INFO=<ID=TYPE,Number=.,Type=String,Description="Variant type: SNP/INDEL/COMPLEX">' + '\n' + \
         '##INFO=<ID=RepRegion,Number=.,Type=String,Description="Repetitive region">' + '\n' + \
         '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">' + '\n' + \
         '##INFO=<ID=UMT,Number=1,Type=Integer,Description="Total used UMI depth">' + '\n' + \
         '##INFO=<ID=VMT,Number=.,Type=Integer,Description="Variant UMI depth">' + '\n' + \
         '##INFO=<ID=VMF,Number=.,Type=Float,Description="Variant UMI allele frequency">' + '\n' + \
         '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' + '\n' + \
         '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Filtered allelic MT depths for the ref and alt alleles">' + '\n' + \
         '##FORMAT=<ID=VF,Number=.,Type=Float,Description="Variant UMI allele frequency, same as VMF">' + '\n' + \
         '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + [sampleName]) + '\n'

   alleles = []
   lastCHROM, lastPOS = '', ''

   outAll.write('\t'.join(headerAll)+'\n')
   outVariants.write('\t'.join(headerVariants)+'\n')
   outLowPi.write('\t'.join(headerLowPi)+'\n')
   outVcf.write(headerVcf)

   cnt = 1
   with open(outlong, 'r') as f:
      next(f)
      for line in f:
         outAll.write(line)
         cnt += 1
         CHROM, POS, REF, ALT, TYPE, sUMT, sForUMT, sRevUMT, sVMT, sForVMT, sRevVMT, sVMF, sForVMF, sRevVMF, VDP, VAF, RefForPrimer, RefRevPrimer, primerOR, pLowQ, hqUmiEff, allUmiEff, refMeanRpb, altMeanRpb, rpbEffectSize, repType, hpInfo, simpleRepeatInfo, tandemRepeatInfo, DP, FR, MT, UFR, sUMT_A, sUMT_T, sUMT_G, sUMT_C, logpval, FILTER = line.strip().split('\t')
         
         if TYPE == '0':
            continue
         
         if ALT == 'DEL': 
            continue

         QUAL = logpval if logpval != 'NA' else '0.00'

         try:
            fQUAL = float(QUAL)
         except ValueError:
            fQUAL = 0.00

         if fQUAL < minCutoff[TYPE.upper()]:
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
         tempVar = (CHROM, POS, REF, ALT, TYPE, DP, VDP, VAF, sUMT, sVMT, sVMF, QUAL, FILTER)
         lenAlleles = len(alleles)

         if fQUAL < cutoff: ## Write to low-PI file
            outLowPi.write(sampleName+'\t'+'\t'.join(tempVar)+'\n')
            continue
            
         # if current chrom and position equal to last line, append it for potential multi-allelic output
         if lenAlleles == 0 or (CHROM == lastCHROM and POS == lastPOS):
            alleles.append(currentAllele)

         # for new chrom or position, if last variant is not multi-allelic, write to vcf directly
         elif lenAlleles == 1:
            biAllelicVar(alleles, RepRegion, outVcf, outVariants)
            alleles = [currentAllele]

         # if last variant is possible multi-allelic, combine and write as one 
         else:            
            multiAllelicVar(alleles, RepRegion, outVcf, outVariants)
            alleles = [currentAllele]

         lastCHROM, lastPOS = CHROM, POS

   # take care of the last line
   lenAlleles = len(alleles)
   if lenAlleles == 1:
      biAllelicVar(alleles, RepRegion, outVcf, outVariants)
   elif lenAlleles >= 2:
      multiAllelicVar(alleles, RepRegion, outVcf, outVariants)
   else:
      pass

   # close all output file handles
   outAll.close()
   outVariants.close()
   outVcf.close()
   outLowPi.close()

#----------------------------------------------------------------------------------------------
#pythonism to run from the command line
#----------------------------------------------------------------------------------------------
if __name__ == "__main__":
   runPath = sys.argv[1]
   outlong = sys.argv[2]
   sampleName = sys.argv[3]
   main(runPath, outlong, sampleName)





