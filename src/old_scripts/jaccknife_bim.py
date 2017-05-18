#!/usr/bin/env python

'''
two resampling options:
- jackknife, where you just leave one out
- bootstrap where you resample

obviously....it would be faster to compute the statistics per block and play with that.

but being lazy....the easier is to remove a block at a time, and then use plink to filter out those sites.

so in short. take the bim file. 
define the intervals. 
create nblock filter files.

use plink to create all bed files. 

create all pairwise Ds.

e.g. with plink:

plink1.9 --bfile yorgos_autosomes_pruned_mind08_geno01 --exclude exclude_blockXX.bim --make-bed --out yorgos_autosomes_pruned_mind08_geno01_wo_blockXX
'''
from __future__ import division

import argparse,sys
import os

parser = argparse.ArgumentParser()

parser.add_argument("--bimfile", help = "bim file from which to create subsets", type = str, required = True)

parser.add_argument("--nblocks",help = "number of blocks for the jackknife", type = int, required = True)

parser.add_argument("--verbose",help = "increase output verbosity", action = "store_true", required = False)    

parser.add_argument("--pathout", help = "path where the blocks should be added", type = str, required = False, default = None)

parser.add_argument("--basebim", help = "basename for the output bim", type = str, required = False, default = None)

parser.add_argument("--runplink",help = "runplink at the end as well", action = "store_true", required = False) 
   
parser.add_argument("--asdpath", help = "path for the asdfile", type = str, required = False, default = "/home/sapfo/Dropbox/mds_genealogy/applications/asdfiles")


args = parser.parse_args()

bimfile = open(args.bimfile).readlines()

pathout = args.pathout
 
if pathout == None:
    pathout = os.path.dirname(args.bimfile)

basebim = args.basebim

if basebim == None:
    basebim = ''.join(os.path.basename(args.bimfile).split('.')[:-1])

nblocks = args.nblocks
nlines = len(bimfile)

print "nblocks: ",nblocks
print "nlines: ",nlines
print "path out: ",pathout
print "basename: ",basebim


sob = int(nlines/nblocks)

print "size of blocks: ",sob

inter_starts = range(0,nlines,sob)

print "size of last block: ", nlines-inter_starts[-1]

if (nlines-inter_starts[-1])/sob <0.9:
    print "the last block in proportion is: ",(nlines-inter_starts[-1])/sob
    inter_starts = inter_starts[:-1]
    print "so the last block will be a bit larger: ",(nlines-inter_starts[-1])

print "actual number of blocks ",len(inter_starts)

inter = inter_starts+[nlines]

print "example file of bimfile: ",[bimfile[0]]

for start_ii,start in enumerate(inter_starts):
    
    STR = ''.join(bimfile[inter[start_ii]:inter[start_ii+1]])
    
    filename = pathout+'/'+"temporary_"+basebim+'_nblocks%s_excl_block%s.bim'%(nblocks,start_ii)
    print "writing file under: ",filename

    out = open(filename,'w')
    out.write(STR)
    out.close()

# runplink
if args.runplink:
    asdcommands = ""
    for start_ii,start in enumerate(inter_starts):

        os.chdir(pathout)

        bimexclude = "temporary_"+basebim+'_nblocks%s_excl_block%s.bim'%(nblocks,start_ii)

        dico = {"basebim":basebim,"bimexclude":bimexclude,"iteration":start_ii,"nblocks":nblocks,"asdpath":args.asdpath,"pathout":pathout}

        plinkout = "%(basebim)s_nblocks%(nblocks)s_wo_block%(iteration)s"%dico

        dico["plinkout"]=plinkout

        plinkcommand = "plink1.9 --bfile %(basebim)s --exclude %(bimexclude)s  --make-bed --out %(plinkout)s"%dico

        print "run plink in this directory :", pathout

        print "plinkcommand: ",plinkcommand
        
        os.system(plinkcommand)


        asdcommands += "./plink2asd.py --pathin %(pathout)s --pathout %(asdpath)s --filenamein %(plinkout)s; rm %(pathout)s/%(plinkout)s.bed %(pathout)s/%(plinkout)s.fam;\n"%dico
    
        
    print asdcommands
    shfile = "%(pathout)s/run_plink2asd_nblocks%(nblocks)s.sh"%dico
    print "written to %s file"%shfile
    

    sh = open(shfile,'w')
    sh.write(asdcommands)
    sh.close()
    
    

#    run_plink2asd_nblocks3.sh

        #try:
        #    os.system("")

    
