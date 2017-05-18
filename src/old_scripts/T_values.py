#!/usr/bin/env python
"""
Created Wed Oct  7 15:04:36 CEST 2015

@author: sapfo
"""
import simul_ms
import python_cmdscale
import exp

import numpy as np
import pylab as py

'''
We want to pick n1, n2, D, T?
Simulate data
Compute the distance matrix
MDS the distance matrix
Get coordinates
Get eigenvalues, eigenvectors
Plot comparing with the other eigenvalues
'''
#################### FIXED #############
n = 100

n1 = 21
n2 = 20
#n3 = 20

D = 10
#D1 = 0.2 #(D1<D)

nreps = 10000
## simulate data
rescaling = 2.0 

#### various popgen numrical values

def TMRCA(n):
    if n!=0:
        return 2*(1-1./n)
    else:
        return 0
        
def Ttotal(n):
    tt = 0
    for ii in range(1,n-1):
        tt+=(1./ii)
    return 2*tt

########### 2 populations ##############
print "########### 2 populations ##############"

#ms simul
params_2pops,data_2pops,tree_lengths_2pops =  simul_ms.ms_two_pops(n1=n1, n2=n2, D=1./rescaling*D,nreps=nreps,verbose=0)
avg_tree_length_2pops = rescaling*np.average(tree_lengths_2pops)

print "tree_lengths: ",rescaling*np.array(tree_lengths_2pops)
print "average tree length (ms): ",rescaling*np.average(tree_lengths_2pops)
print "TMRCA pop1: ",TMRCA(n1)
print "TMRCA pop2: ",TMRCA(n2)
print "Ttotal pop1: ",Ttotal(n1)
print "Ttotal pop2: ",Ttotal(n2)
print "Total pop1 + pop2:",Ttotal(n1)+Ttotal(n2)


