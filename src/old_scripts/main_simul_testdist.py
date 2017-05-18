#!/usr/bin/env python
"""
Created Wed Oct  7 15:04:36 CEST 2015

@author: sapfo
"""
import matplotlib 
matplotlib.use("Agg")

import simul_ms
import python_cmdscale
import exp
import sys

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

n1 = 5
n2 = 5
n3 = 5

D = 0.4
D1 = 0.1 #(D1<D)

nreps = 10000
## simulate data
rescaling = 2.0 

########### 2 populations ##############
print "########### 2 populations ##############"

#ms simul
params_2pops,data_2pops,tree_lengths_2pops =  simul_ms.ms_two_pops(n1=n1, n2=n2, D=1./rescaling*D,nreps=nreps,verbose=0)
avg_tree_length_2pops = rescaling*np.average(tree_lengths_2pops)

Delta_2pops = simul_ms.distance_matrix(data=data_2pops,rescale="all", n1=n1,n2=n2,verbose=1)


#cmdscale
evals_2pops, evecs_2pops, Y_2pops = python_cmdscale.cmdscale(Delta_2pops)
exp.T_D_two_pops(eigenvalues = evals_2pops,n1=n1,n2=n2,diploid=2)

# analytical
params_exp_2pops,evals_exp_2pops, evec_exp_2pops = exp.two_pops(n1=n1, n2=n2, D=D, T=avg_tree_length_2pops)


print "params_2pops (ms): ",params_2pops
print "params_exp_2pops: ",params_exp_2pops
print "average tree length (ms): ",rescaling*np.average(tree_lengths_2pops)
#print "expected tree length (coal): ",exp_tree_length

print "expected lambda1 (analytical): ",evals_exp_2pops[0]
print "observed lambda1 (cmdscale): ",evals_2pops[0]

print "expected lambda2 (analytical): ",evals_exp_2pops[1]
print "observed lambda2 (cmdscale): ",evals_2pops[1]
print "average observed lambda2...n-1 (cmdscale): ",np.average(evals_2pops[1:-1])

print evals_exp_2pops[:10]
print evals_2pops[:10]

#print "observed lambda1 (mds): ",evals[0]
#print "observed average lambdas (mds): ",np.average(evals[:-1])

## observed and expected distance between the two populations

dist_obs = np.abs(np.average(Y_2pops[:n1,0])-np.average(Y_2pops[n1:,0]))
dist_exp = exp.dist_two_populations(T=avg_tree_length_2pops, D=D, n1=n1, n2=n2)
dist_exp_resc_diag = exp.dist_two_populations_resc_diag(T=avg_tree_length_2pops, D=D, n1=n1, n2=n2)
dist_exp_resc_all = exp.dist_two_populations_resc_all(T=avg_tree_length_2pops, D=D, n1=n1, n2=n2)

print "dist_obs: ",dist_obs
print "dist_exp: ",dist_exp
print "dist_exp_resc_diag: ",dist_exp_resc_diag
print "dist_exp_resc_all: ",dist_exp_resc_all
sys.exit()

### plotting two population ###
py.figure()
py.plot(Y_2pops[:n1,0],Y_2pops[:n1,1],'x',color='orange')
py.plot(Y_2pops[n1:,0],Y_2pops[n1:,1],'o',color='blue')

py.title("simulations 2 pops n1 = %s, n2 = %s, D = %s, nreps = %s "%(n1,n2,D,nreps))
py.xlabel("dim 1")
py.ylabel("dim 2")

#py.xlabel("dim 1 (%.2f %%)"%(1.*evals[0]/np.average(evals[:-1])))
#py.ylabel("dim 2 (%.2f %%)"%(1.*evals[1]/np.average(evals[:-1])))
#py.show()

########### 3 populations ##############
print "########### 3 populations ##############"
nreps = 100
#ms simul
params_3pops,data_3pops,tree_lengths_3pops =  simul_ms.ms_three_pops(n1=n1, n2=n2, n3=n3, D=1./rescaling*D, D1 = 1./rescaling*D1,nreps=nreps,verbose=0)
avg_tree_length_3pops = rescaling*np.average(tree_lengths_3pops)
Delta_3pops = simul_ms.distance_matrix(data=data_3pops,verbose=0)

#cmdscale
evals_3pops, evecs_3pops, Y_3pops = python_cmdscale.cmdscale(Delta_3pops)

try:
    Texp,Dexp,D1exp,Drescaledexp,D1rescaledexp = exp.T_D_D1_three_pops(eigenvalues = evals_3pops,n1=n1,n2=n2,n3=n3,diploid=2)
except:
    Texp,Dexp,D1exp,Drescaledexp,D1rescaledexp= 1,1,1,1,1

print "average tree length (ms): ",rescaling*np.average(tree_lengths_3pops)
print "params_3pops (ms): ",params_3pops

# analytical
params_exp_3pops,evals_exp_3pops, evec_exp_3pops = exp.three_pops(n1=n1, n2=n2, n3=n3, D=D, D1=D1, T=avg_tree_length_3pops)

print "params_3pops (ms): ",params_3pops
print "params_exp_3pops: ",params_exp_3pops
print "average tree length (ms): ",rescaling*np.average(tree_lengths_3pops)
#print "expected tree length (coal): ",exp_tree_length

print "expected lambda1 (analytical): ",evals_exp_3pops[0]
print "observed lambda1 (cmdscale): ",evals_3pops[0]
print ""
print "expected lambda2 (analytical): ",evals_exp_3pops[1]
print "observed lambda2 (cmdscale): ",evals_3pops[1]
print ""
print "expected lambda3 (analytical): ",evals_exp_3pops[2]
print "observed lambda3 (cmdscale): ",evals_3pops[2]

print "average observed lambda3...n-1 (cmdscale): ",np.average(evals_3pops[2:-1])

print evals_exp_3pops[:10]
print evals_3pops[:10]

sys.exit()
### plotting three population ###
py.figure()
py.plot(Y_3pops[:,0][:n1],Y_3pops[:,1][:n1],'D',color='orange')
py.plot(Y_3pops[:,0][n1:n1+n2],Y_3pops[:,1][n1:n1+n2],'o',color='blue')
py.plot(Y_3pops[:,0][n1+n2:],Y_3pops[:,1][n1+n2:],'v',color='green')

py.title("simulations 3 pops n1 = %(n1)s, n2 = %(n2)s, n3 =  %(n3)s, D = %(D)s, D1 = %(D1)s, nreps = %(nreps)s "%params_3pops)
py.xlabel("dim 1")
py.ylabel("dim 2")
py.show()

########### 4 populations and above ##############


