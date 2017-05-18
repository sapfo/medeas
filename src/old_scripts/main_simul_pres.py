#!/usr/bin/env python
"""
Created Wed Oct  7 15:04:36 CEST 2015

@author: sapfo
"""
import matplotlib
matplotlib.use('Agg')

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

nreps = 2
## simulate data
rescaling = 2.0 

#figure property
#matplotlib.rc('xtick', labelsize=12) 
#matplotlib.rc('ytick', labelsize=12)

'''
########### 1 population ##############
print "########### 1 population ##############"

n = 3
## expected tree length for one population
exp_tree_length = 0
for i in range(2,n+1):
    exp_tree_length += 2./(i-1)

npoints = 30

lambdas = []
Xs = np.logspace(1,4,npoints)
#Xs = [1000]
for nrep in Xs:
   
    params,data,tree_lengths =  simul_ms.ms_one_pops(n=n,nreps=nrep,verbose=0)
    Delta = simul_ms.distance_matrix(data=data,verbose=0)
    print "Delta: ",Delta
    evals, evecs, Y = python_cmdscale.cmdscale(Delta)       

    print "params: ",params
    print "average tree length (ms): ",rescaling*np.average(tree_lengths)
    print "expected tree length (coal): ",exp_tree_length
    print "expected lambda1 (analytical): ",2./((rescaling*np.average(tree_lengths))**2)
    print "observed lambda1 (mds): ",evals[0]
    print "observed average lambdas (mds): ",np.average(evals[:-1])
    print "evals (first 10): ",evals[:10]
        
    lambdas.append(np.average(evals[:-1]))
    
if npoints==1:                        
   ### plotting one population ###
    py.figure(1)
    py.subplot(1,1,1)
    py.plot(Y[:,0],(Y[:,1]),'o',color='black')
    py.title("simulations 1 population n = %s, nreps = %s "%(n,nrep))
    py.xlabel("dim 1 (%.2f %%)"%(1.*evals[0]/np.average(evals[:-1])))
    py.ylabel("dim 2 (%.2f %%)"%(1.*evals[1]/np.average(evals[:-1])))
    ''py.subplot(1,3,2)
    py.plot(Y[:,0],(Y[:,2]),'o',color='black')
    py.title("simulations 1 population n = %s, nreps = %s "%(n,nrep))        
    py.xlabel("dim 1 (%.2f %%)"%(1.*evals[0]/np.average(evals[:-1])))
    py.ylabel("dim 3 (%.2f %%)"%(1.*evals[2]/np.average(evals[:-1])))
    py.subplot(1,3,3)
    py.plot(Y[:,1],(Y[:,2]),'o',color='black')
    py.title("simulations 1 population n = %s, nreps = %s "%(n,nrep))
    py.xlabel("dim 2 (%.2f %%)"%(1.*evals[1]/np.average(evals[:-1])))
    py.ylabel("dim 3 (%.2f %%)"%(1.*evals[2]/np.average(evals[:-1])))'

    py.savefig("1pop_projections_%s.pdf"%np.random.randint(1,100,1)[0])

        
if npoints>5:
    fig = py.figure(2)        
    ax = fig.add_subplot(1,1,1)

    py.plot(Xs,lambdas,'o-',color='black')
    py.hlines(2./(exp_tree_length**2),min(Xs),max(Xs),'r')
    ax.set_xscale('log')
    py.savefig("1pop_lambdas_%s.pdf"%np.random.randint(1,100,1)[0])

print lambdas
        
#py.show()
'''

'''
########### 2 populations ##############
print "########### 2 populations ##############"
npoints = 20
Xs = np.logspace(1,4,npoints)
Ds = []; Ts = []
D = 0.2

for nreps in Xs:
    #ms simul
    params_2pops,data_2pops,tree_lengths_2pops =  simul_ms.ms_two_pops(n1=n1, n2=n2, D=1./rescaling*D,nreps=nreps,verbose=0)
    avg_tree_length_2pops = rescaling*np.average(tree_lengths_2pops)
    Delta_2pops = simul_ms.distance_matrix(data=data_2pops,verbose=0)


    #cmdscale
    evals_2pops, evecs_2pops, Y_2pops = python_cmdscale.cmdscale(Delta_2pops)
    Texp,Dexp,Dexprescaled = exp.T_D_two_pops(eigenvalues = evals_2pops,n1=n1,n2=n2,diploid=2)
    Ts.append(Texp)
    Ds.append(Dexp)
    
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
    
fig = py.figure(1,figsize=(5,2.5))        
ax = fig.add_subplot(1,1,1)
py.plot(Xs,Ds,'o-',color='black')
py.hlines(D,min(Xs),max(Xs),'r')
ax.set_xscale('log')
 
ax = fig.add_subplot(2,1,2)
py.plot(Xs,Ts,'o-',color='black')
py.hlines(avg_tree_length_2pops,min(Xs),max(Xs),'r')
ax.set_xscale('log')
  
py.savefig("2pop_Ds_Ts_%s.pdf"%(np.random.randint(1,100,1)[0]))
'''
    
#print "observed lambda1 (mds): ",evals[0]
#print "observed average lambdas (mds): ",np.average(evals[:-1])

### plotting two population ###
if 0:
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

npoints = 20
n1 = 20; n2 = 20; n3 = 20;

Xs = np.logspace(1,4,npoints)
Ds = []; Ts = []
D = 0.4
D1 = 0.2

nreps = 100

Ds = []
D1s =  []
Ts = []
#ms simul
for nreps in Xs:
    params_3pops,data_3pops,tree_lengths_3pops =  simul_ms.ms_three_pops(n1=n1, n2=n2, n3=n3, D=1./rescaling*D, D1 = 1./rescaling*D1,nreps=nreps,verbose=0)
    avg_tree_length_3pops = rescaling*np.average(tree_lengths_3pops)
    Delta_3pops = simul_ms.distance_matrix(data=data_3pops,verbose=0)

    #cmdscale
    evals_3pops, evecs_3pops, Y_3pops = python_cmdscale.cmdscale(Delta_3pops)

    #try
    Texp,Dexp,D1exp,Drescaledexp,D1rescaledexp = exp.T_D_D1_three_pops(eigenvalues = evals_3pops,n1=n1,n2=n2,n3=n3,diploid=2)
    
    #except:
    #    Texp,Dexp,D1exp,Drescaledexp,D1rescaledexp= 1,1,1,1,1

    Ds.append(Dexp)
    D1s.append(D1exp)
    '''
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
    print evals_3pops[:10]'''

fig = py.figure(1,figsize=(5,3))        
ax = fig.add_subplot(2,1,1)
py.plot(Xs,Ds,'o-',color='black')
py.hlines(D,min(Xs),max(Xs),'r')
ax.set_xscale('log')
 
ax = fig.add_subplot(2,1,2)
py.plot(Xs,D1s,'o-',color='black')
py.hlines(D1,min(Xs),max(Xs),'r')
ax.set_xscale('log')
  
py.savefig("3pop_Ds_D1s_%s.pdf"%(np.random.randint(1,100,1)[0]))


sys.exit()

'''### plotting three population ###
py.figure()
py.plot(Y_3pops[:,0][:n1],Y_3pops[:,1][:n1],'D',color='orange')
py.plot(Y_3pops[:,0][n1:n1+n2],Y_3pops[:,1][n1:n1+n2],'o',color='blue')
py.plot(Y_3pops[:,0][n1+n2:],Y_3pops[:,1][n1+n2:],'v',color='green')

py.title("simulations 3 pops n1 = %(n1)s, n2 = %(n2)s, n3 =  %(n3)s, D = %(D)s, D1 = %(D1)s, nreps = %(nreps)s "%params_3pops)
py.xlabel("dim 1")
py.ylabel("dim 2")
py.show()'''

########### 4 populations and above ##############


