#!/usr/bin/env python

"""the goal here is to parse a ped file (NOT PARSING A PED FILE, JUST ASD FILE), choose the populations, stick to data with non missing sites

and 
1. compute MDS.
2. compute PCA.

relax the missingness assumption start over.

in each case compute 
eigenvalues and eigenvectors and from that
T, D and D_prime and ...D_primeprime.

considering either 1, 2, 3, 4 populations.

Suggested. 
Africans
East Asians
Europeans

and or
Australians (WCD)
Papuans (HGDP? or and others)


a. load matrix
b. cmdscale
c. plot dim1, dim2; dim1, dim3; dim2, dim3
d. plot eigenvalues
e. plot eigenvectors 1, 2, 3
f. assuming a certain number of population compute T, D, Dprime etc (need to solve lambda2DT_pop2 etc)."""

import matplotlib
matplotlib.use('Agg')

import numpy as np
import python_cmdscale
import exp
import pylab as py
import sys
import argparse
import matplotlib.cm as cm

                                               
parser = argparse.ArgumentParser()
parser.add_argument("--verbose", help = "increase output verbosity",action="store_true")
parser.add_argument("--asdfile", help = "input asd file", required = True, type = str) 
parser.add_argument("--basename", help = "output basename",required = False, type = str, default = False) 

parser.add_argument("--labelfile", help = "input label file",required = False, type = str, default = None) 

parser.add_argument('--filter', nargs='+', help='<Required> Set flag', required=False, default = 'all')

parser.add_argument('--out', help='pdf out', required=False, default = None)

args = parser.parse_args()                                                       
                                                                                
verbose = args.verbose
basename = args.basename
if not basename:    
    basename = args.asdfile[:-4]
    
    
asdfile = args.asdfile
D = np.loadtxt(asdfile)

if args.labelfile == None:
    args.labelfile = asdfile+".labels"

labels = np.loadtxt(args.labelfile,dtype=str)

#if verbose:
 
print "2 first labels: ",labels[:2]
print "shapes of labels: ",np.shape(labels)
print "labels type: ",type(labels)

if len(labels[0])==2:   
    if len(np.shape(labels))>1:
        pops = labels[:,0]
    else:
        pops = labels
    
else: 
    print "expecting population labels in one row and sample labels in the next row"
    pops = labels[0,:]

upops = list(set(pops))
upops.sort()
print "unique populations in asdfile before filtering: ",' '.join(upops)

if args.filter != "all": 
    pops_toinclude = args.filter
    if verbose: print "populations to be included: ",pops_toinclude
    indices_toinclude = []

    for pop_toinclude in pops_toinclude:
        indices_toinclude += [i for i, x in enumerate(pops) if x == pop_toinclude]

    if verbose: print "indices_toinclude: ",indices_toinclude

    indices_toremove = np.delete(np.arange(len(labels)),indices_toinclude,axis=0)
    
    if verbose: print "indices_toremove: ",indices_toremove

    pops = np.delete(pops,indices_toremove,axis=0)
    D = np.delete(D,indices_toremove,axis=0)
    D = np.delete(D,indices_toremove,axis=1)

print "unique pops after filtering: ",set(pops)
print "shape of D: ",np.shape(D)

nsamples = np.shape(D)[0]

#if pops_toinclude == None:
# the population labels can be either 1. over two lines: first line population, second line samples. or 2. two columns: first column population, second column sample. 

unique_pops = list(set(pops))
unique_pops.sort()
nb_pops = len(unique_pops)

if nb_pops == 2:
    n1 = list(pops).count(unique_pops[0])
    n2 = list(pops).count(unique_pops[1])
    print "n_%s = %s"%(unique_pops[0],n1)
    print "n_%s = %s"%(unique_pops[1],n2)

#set_cmap('Set2')
colors = cm.jet(np.linspace(0,1,nb_pops))

#colors = len(pops)*['peru','dodgerblue','orange','dargrey','olive','purple','orange','yellow','aqua','tomato','cyan','steelblue']

shapes = len(pops)*['o','D','v']

dico_colors = dict(zip(unique_pops,colors[:nb_pops]))
dico_shapes = dict(zip(unique_pops,shapes[:nb_pops]))

print "dico_colors: ",dico_colors
print "dico_shapes: ",dico_shapes

pops_colors = [dico_colors[pp] for pp in pops]
pops_shapes = [dico_shapes[pp] for pp in pops]

print "unique population labels: ",unique_pops
print "number of unique population labels: ",len(unique_pops)


print "computing cmdscale..."
evals, evecs, Y = python_cmdscale.cmdscale(D)

n = len(evals)

print "evals: ",evals
print "number of evals (or number of individuals): ",n

if nb_pops == 2:    
    params = exp.T_D_two_pops(eigenvalues=evals,n1=n1,n2=n2,gen_time=29,Ne=10000,diploid=2,verbose=True)

    dico_parameters = dict(zip(['nb_pops','T','D','Drescaled'],[nb_pops]+list(params)))

    print "assuming %(nb_pops)s populations, we find that T = %(T)s , D = %(D)s, D(years) = %(Drescaled)s"%dico_parameters


print "plotting stuff..."
fig = py.figure(figsize=(25,12))
#title = "Populations: %s"%(' '.join([unique_pops[ii]+' ('+colors[ii]+')' for ii in range(nb_pops)]))

title = "Populations: %s"%(' '.join([unique_pops[ii] for ii in range(nb_pops)]))

print title
py.suptitle(title) 
py.subplot(2,4,1)
py.title("classic mds (dim1, dim2)")

legend_labels = []
legend_points = []
for ii,pop in enumerate(pops): 
    # collect the first appearance of each pops
    print "pop: ",pop
    if pop in legend_labels:
        py.scatter(Y[ii,0],Y[ii,1],marker=pops_shapes[ii],color=pops_colors[ii])
    else:
        legend_points.append(py.scatter(Y[ii,0],Y[ii,1],marker=pops_shapes[ii],color=pops_colors[ii],label=pop))
        legend_labels.append(pop)

print legend_labels

py.xlabel("dim 1 (%.2f %%)"%(1.*evals[0]/np.average(evals[:-1])))
py.ylabel("dim 2 (%.2f %%)"%(1.*evals[1]/np.average(evals[:-1])))
#py.xlim(-0.15,0.2)
#py.ylim(-0.1,0.15)

py.subplot(2,4,2)
py.title("classic mds (dim1, dim2)")

for ii in range(nsamples): py.scatter(Y[ii,0],Y[ii,2],marker=pops_shapes[ii],color=pops_colors[ii])

py.xlabel("dim 1 (%.2f %%)"%(1.*evals[0]/np.average(evals[:-1])))
py.ylabel("dim 3 (%.2f %%)"%(1.*evals[2]/np.average(evals[:-1])))

py.subplot(2,4,3)
py.title("classic mds (dim2, dim3)")
for ii in range(nsamples): py.scatter(Y[ii,1],Y[ii,2],marker=pops_shapes[ii],color=pops_colors[ii])

py.xlabel("dim 2 (%.2f %%)"%(1.*evals[2]/np.average(evals[:-1])))
py.ylabel("dim 3 (%.2f %%)"%(1.*evals[2]/np.average(evals[:-1])))

py.subplot(2,4,4)
py.title("eigenvalues")
py.plot(range(1,n),evals[:-1],'o',color='black') 
py.xlim(-1,len(evals)+1)

py.subplot(2,4,5)              
py.title("eigenvector 1")
for ii in range(n):py.scatter(ii+1,evecs.T[0][ii],marker=pops_shapes[ii],color=pops_colors[ii])
py.xlabel("individuals (numbered)")

py.subplot(2,4,6)
py.title("eigenvector 2") 
for ii in range(n):py.scatter(ii+1,evecs.T[1][ii],marker=pops_shapes[ii],color=pops_colors[ii])
py.xlabel("individuals (numbered)")

py.subplot(2,4,7) 
py.title("eigenvector 3")
for ii in range(n):py.scatter(ii+1,evecs.T[2][ii],marker=pops_shapes[ii],color=pops_colors[ii])
py.xlabel("individuals (numbered)")

fig.legend(legend_points,legend_labels,'lower right')

if args.out == None:
    pdfout = basename+'.pdf'

else:
    if args.out.rstrip()[-3:]=='pdf':        
        pdfout = args.out
    else:
        pdfout = args.out+'.pdf'

print "saving file under %s"%pdfout
py.savefig(pdfout)
#py.show()








