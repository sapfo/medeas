#!/usr/bin/env python
from __future__ import division

'''
genotypes:
3: missing
0: homozygote for ref allele1: ref1, ref1
1: heterozygote: ref1, ref2
2: homozygote for ref allele2: ref2, ref2
'''
import sys
import numpy as np
from plinkio import plinkfile
import argparse

#plink_file = plinkfile.open( "/home/sapfo/Projects/australian/mds/masked_yorgos/yorgos_autosomes_pruned_mind08_geno01" )
parser = argparse.ArgumentParser()

parser.add_argument("--pathin", help = "input path where the ped file is", required = False, type = str, default = "/home/sapfo/Projects/australian/mds/masked_yorgos") #or testdata
parser.add_argument("--pathout", help = "output path where the asd file should go", required = False, type = str, default = "asdfiles") 
parser.add_argument("--filenamein", help = "basename of the ped file", required = True, type = str) 
parser.add_argument("--nsamples",help = "nsamples", required = False, type = int, default = None)

args = parser.parse_args()                                                       

pathin = args.pathin
pathout = args.pathout
#pathin = "/home/sapfo/Projects/australian/mds/masked_yorgos"
#pathout = "asdfiles"

filenamein = args.filenamein #"yorgos_autosomes_pruned_mind04_geno01"
 
print "reading in %s/%s"%(pathin,filenamein)

plink_file = plinkfile.open( "%s/%s"%(pathin,filenamein))
fam_file = pathin+'/'+filenamein+'.fam'

asd_file = pathout+'/'+filenamein+".asd"
label_file = pathout+'/'+filenamein+'.asd.labels'

if not plink_file.one_locus_per_row( ):
     print( "This script requires that snps are rows and samples columns." )
     exit( 1 )

sample_list = plink_file.get_samples( )
locus_list = plink_file.get_loci( )

#for locus, row in zip( locus_list, plink_file ):
#  for sample, genotype in zip( sample_list, row ):
 #       print( "Individual {0} has genotype {1} for snp {2}.".format( sample.iid, genotype, locus.name ) )

verbose = False
cutoff_snps = 1e3
cutoff_samples = 1000

nsnp = 0


genoT = np.array([list(row) for row in plink_file]) # row is a specific snp
geno = genoT.T

print "shape of data ", np.shape(geno)


#a1 = ar[0][np.where((ar[0]!=3) & (ar[1]!=3))]
#a1 = ar[0][np.where((ar[0]!=3) & (ar[1]!=3))]

nsamples = args.nsamples
if nsamples == None:
     nsamples = np.shape(geno)[0]

nsnps = np.shape(geno)[1]

print "nsites: ", nsnps
print "nsamples: ", nsamples

ASD = np.zeros((nsamples,nsamples))
print "computing asd file ..."
for i1 in range(nsamples):
     for i2 in range(i1,nsamples):
          if verbose: print i1,i2

          a1 = geno[i1][np.where((geno[i1]!=3) & (geno[i2]!=3))]
          a2 = geno[i2][np.where((geno[i1]!=3) & (geno[i2]!=3))]
          
          ASD[i1,i2] = np.sum(np.abs(a1-a2))/len(a1)/2
          ASD[i2,i1] = ASD[i1,i2]

 
print "ASD file: ", ASD
print "asdfile is being saved under: ", asd_file

np.savetxt(asd_file, ASD)

# save labels
labels = [xx.split()[0]+' '+xx.split()[1]+'\n' for xx in open(fam_file).readlines()[:nsamples]]
#labels += [xx.split()[1]+'\n' for xx in open(fam_file).readlines()[:nsamples]]+['\n']

print "labelfile is being saved under: ", label_file

labelout = open(label_file,"w")
labelout.write("".join(labels))
labelout.close()

if len(labels)>10:
     print "labels: ",np.array(labels)[:10]
else: 
     print "labels: ",np.array(labels)
  

np.savetxt(asd_file, ASD)
