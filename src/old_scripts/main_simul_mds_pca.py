#!/usr/bin/env python
"""
Created Wed Oct  7 15:04:36 CEST 2015

@author: sapfo
"""
import matplotlib 
#matplotlib.use('Agg')

import simul_ms
import python_cmdscale
import python_pca
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
n = 10

n1 = 5
n2 = 5
n3 = 5

D = 0.4
D1 = 0.1 #(D1<D)

nreps = 100
## simulate data
rescaling = 2.0 

verbose = True
########### 1 population ##############
print "########### 1 population ##############"
## expected tree length for one population
exp_tree_length = 0
for i in range(2,n+1):
    exp_tree_length += 2./(i-1)

nsnps = [100,1000]

T_mds = {}
T_pca = {}
for nsnp in nsnps:

    T_mds[nsnp] = []
    T_pca[nsnp] = []

    for iteration in range(nreps):
        params,data,tree_lengths =  simul_ms.ms_one_pops(n=n,nreps=nsnp,verbose=0)

        Delta = simul_ms.distance_matrix(data=data,verbose=0)

        #print "data: ",data
        #print "size of data: ",data.shape
        #print "Delta: ",Delta

        evals_mds, evecs_mds, Y_mds = python_cmdscale.cmdscale(Delta)

        evals_pca, evecs_pca, Y_pca = python_pca.PCA(data.T)

        #print "params: ",params
        if verbose: print "average tree length (computed with ms): ",rescaling*np.average(tree_lengths)

        if verbose: print "expected tree length (analytical coal): ",exp_tree_length

        # mds expected total tree length, bias, rmse
        t_mds = (2./(np.average(evals_mds[:-1])))**(1/2.)
        T_mds[nsnp].append(t_mds)

        if verbose: print "expected T (mds) from eigenvalues: ",T_mds

        # pca expected tree length, bias, rmse
        t_pca = 1./np.average(evals_pca[:-1])
        T_pca[nsnp].append(t_pca)

        if verbose: print "expected T (pca) from eigenvalues: ",T_pca


    #print "expected lambda1 (mds) for (Ivan analytical): ",2./((exp_tree_length)**2)
    #print "expected lambda1 (pca) for (Ivan analytical): ",1./((exp_tree_length))


    #print "observed lambda1 (mds procedure): ",evals_mds[0]
    #print "observed lambda1 (pca procedure): ",evals_pca[0]
    #print "observed average lambdas (mds): ",np.average(evals_mds[:-1])
    #print "observed average lambdas (pca): ",np.average(evals_pca[:-1])

    #print "evals (first 10): ",evals_mds[:10]

Bias_pca = []
Bias_mds = []
Var_pca = []
Var_mds = []

for nsnp in nsnps:
    Bias_mds.append(np.average(T_mds[nsnp])-exp_tree_length)
    Bias_pca.append(np.average(T_pca[nsnp])-exp_tree_length)

    Var_mds.append(np.var(T_mds[nsnp]))
    Var_pca.append(np.var(T_pca[nsnp]))

print Bias_mds
print Bias_pca

print Var_mds
print Var_pca



fig = py.figure()
py.suptitle("1 population: estimate of T, mds vs pca")
ax1 = fig.add_subplot(3,1,1)

py.title("T")
py.plot(nsnps,[T_mds[ss][0] for ss in nsnps],'o',color='grey',label='mds'\
)
py.plot(1.1*np.array((nsnps)),[T_pca[ss][0] for ss in nsnps],'^',color='black', label = 'pca')

py.plot(nsnps,[T_mds[ss] for ss in nsnps],'o',color='grey')
py.plot(1.1*np.array((nsnps)),[T_pca[ss] for ss in nsnps],'^',color='black')
py.hlines(exp_tree_length,min(nsnps),max(nsnps)*1.1,'r')

ax1.set_xscale('log')
py.legend(loc='best')

ax2=fig.add_subplot(3,1,2)
py.title("Bias")
py.plot(nsnps,Bias_mds,'o-',color = "grey",label='mds')
py.plot(1.1*np.array(nsnps), Bias_pca,'^-',color = "black",label='pca')

ax2.set_xscale('log')

ax3 = fig.add_subplot(3,1,3)
py.title("Variance")
py.plot(nsnps,Var_mds,'o-', color = 'grey',label='mds')
py.plot(1.1*np.array(nsnps),Var_pca,'^-', color = 'black',label='pca')

ax3.set_xscale('log')

py.savefig("1pop_mds_pca_T.pdf")
py.show()

sys.exit()


### plotting one population ###
py.plot(Y[:,0],(Y[:,1]),'o',color='blue')
py.title("simulations 1 population n = %s, nreps = %s "%(n,nreps))
py.xlabel("dim 1 (%.2f %%)"%(1.*evals[0]/np.average(evals[:-1])))
py.ylabel("dim 2 (%.2f %%)"%(1.*evals[1]/np.average(evals[:-1])))
########### 2 populations ##############
print "########### 2 populations ##############"

#ms simul
params_2pops,data_2pops,tree_lengths_2pops =  simul_ms.ms_two_pops(n1=n1, n2=n2, D=1./rescaling*D,nreps=nreps,verbose=0)
avg_tree_length_2pops = rescaling*np.average(tree_lengths_2pops)
Delta_2pops = simul_ms.distance_matrix(data=data_2pops,verbose=0)


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


### plotting two population ###
py.figure()
py.plot(Y_2pops[:,0][:n1],Y_2pops[:,1][:n1],'x',color='orange')
py.plot(Y_2pops[:,0][n1:],Y_2pops[:,1][n1:],'o',color='blue')

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


