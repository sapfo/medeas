#!/usr/bin/env python

"""
Created Wed Oct  7 11:49:12 CEST 2015

@author: sapfo
"""
from __future__ import division

import numpy as np
import matplotlib.pyplot as pl
import os

def parse_ms(n=10,nreps=10,nsites=1,out = ['ms 7 1 -s 1 -I 2 5 2 -ej 0.100000 1 2 -L \n', '35929 49668 37267\n', '\n', '//\n', 'time:\t0.339923\t1.348038\n', 'segsites: 1\n', 'positions: 0.0475 \n', '1\n', '1\n', '1\n', '0\n', '1\n', '1\n', '1\n'],verbose=0):
    tree_lengths = []
    data = np.zeros((n,nreps))
    dd = []
    site_counter = -1
    for ll in out[3:]: # starts at first locus
        ll=ll.rstrip()
        #print 'llori: ',ll
        if ll=='//': # new locus
            #print "dd ",dd
            if dd!=[]:data[:,site_counter]=dd 
            #print ll
            dd = []
            site_counter += 1
        elif ll[0:4]=="time": # the timeline gives the probability of the tree folloyed by the tree length
            tt=ll.split()[-1]            
            tree_lengths.append(float(tt))
        elif ':' not in ll: # ll is the data here: 0 or 1 (considering single locus)
            if ll!='':                
                #print 'll: ',ll
                dd.append(ll[0])                                 
    
    if verbose: 
        print "tree_lengths: ",tree_lengths
    return data,tree_lengths

    
def ms_one_pops(n=5, nreps=20,verbose=1):
    '''
    ms nsam nreps -s segsites -I npop n1 n2 ... [4N0m] -ej t i j
    -I  assume island model with symmetric migration rate, 4N0m, and sample configuration n1, n2 . . . .
    -s j make samples with fixed number of segregating sites, j.
    -ej move all lineages in subpopulation i to subpopulation j at time t.    
    '''

    params = dict(zip(['n','nreps'],[n,nreps]))

    mscommand = '''ms %(n)i %(nreps)i -s 1 -L'''%params
    if verbose: print "mscommand:\n",mscommand
    
    out = os.popen(mscommand).readlines()
        
    data,tree_lengths = parse_ms(n=n,nreps=nreps,nsites=1,out=out,verbose=verbose)
    return params,data,tree_lengths

def ms_two_pops(n1=5, n2=2, D=0.1, nreps=20,verbose=1):
    '''
    ms nsam nreps -s segsites -I npop n1 n2 ... [4N0m] -ej t i j
    -I  assume island model with symmetric migration rate, 4N0m, and sample configuration n1, n2 . . . .
    -s j make samples with fixed number of segregating sites, j.
    -ej move all lineages in subpopulation i to subpopulation j at time t.    
    '''

    n = n1+n2

    params = dict(zip(['n','n1','n2','D','nreps'],[n,n1,n2,D,nreps]))

    mscommand = '''ms %(n)i %(nreps)i -s 1 -I 2 %(n1)i %(n2)i -ej %(D)f 1 2 -L'''%params
    if verbose: print "mscommand:\n",mscommand
    
    out = os.popen(mscommand).readlines()
        
    data,tree_lengths = parse_ms(n=n,nreps=nreps,nsites=1,out=out)
    return params,data,tree_lengths


def ms_three_pops(n1=5, n2=2, n3=10, D=0.2, D1=0.1, nreps=20,verbose=1):
    '''
    ms nsam nreps -s segsites -I npop n1 n2 ... [4N0m] -ej t i j
    -I  assume island model with symmetric migration rate, 4N0m, and sample configuration n1, n2 . . . .
    -s j make samples with fixed number of segregating sites, j.
    -ej move all lineages in subpopulation i to subpopulation j at time t.    
    assuming here D1<D
    '''

    n = n1+n2+n3

    params = dict(zip(['n','n1','n2','n3','D1','D','nreps'],[n,n1,n2,n3,D1,D,nreps]))

    mscommand = '''ms %(n)i %(nreps)i -s 1 -I 3 %(n1)i %(n2)i %(n3)i -ej %(D1)f 1 2 -ej %(D)f 2 3 -L'''%params
    if verbose: print "mscommand:\n",mscommand
    
    out = os.popen(mscommand).readlines()
        
    data,tree_lengths = parse_ms(n=n,nreps=nreps,nsites=1,out=out)
    return params,data,tree_lengths


def distance_matrix(data=np.array([[0,1],[1,0],[1,1]]),verbose=0, n1=None, n2=None, rescale=False):
    data = np.array(data)    
    n = len(data)
    nssites = len(data[0])
    if verbose: print "expecting %s samples"%n
    Delta = np.zeros((n,n))
    for index_i,indiv_i in enumerate(data):
        for index_j,indiv_j in enumerate(data):
            
            if index_i > index_j:
                Delta[index_i,index_j]=Delta[index_j,index_i]
            else:
                #print indiv_i; print indiv_j; print index_i,index_j,print np.sum(indiv_i*indiv_j+((1-indiv_i)*(1-indiv_j)))/nssites
                Delta[index_i,index_j]=1-(1.*np.sum(indiv_i*indiv_j+((1-indiv_i)*(1-indiv_j)))/nssites) #1 - nmatches
                if rescale=='diag' or rescale=='all':
                    if verbose: print "i am being rescaled"

                    if (index_i<n1) and (index_j<n1):
                        Delta[index_i,index_j] *= np.sqrt(n1/(n1+n2))
                    elif (index_i>=n1) and (index_j>=n1):
                        Delta[index_i, index_j] *= np.sqrt(n2/(n1+n2))

                    if rescale =='all':
                        if (index_i<n1) and (index_j>=n1):
                            print "I AM BEING RESCALED!!!"
                            Delta[index_i,index_j] *= np.sqrt(1/2.-1./n)
                        
                
    if verbose: print Delta        
    return Delta

#def ms_three_pops(n1=2, n2=3, n3=3, D=0.1, D1=0.2, T=0.3):

#params, data, tree_lengths = ms_two_pops(n1=80, n2=20, D=0.1, nreps=100,verbose=1)
'''
params, data, tree_lengths = ms_three_pops(n1=2, n2=3, n3=2, D=1, D1=4, nreps=10,verbose=1)
Delta = distance_matrix(data=data)

#print params
#print data
#print tree_lengths
print np.average(tree_lengths)
print Delta
'''
