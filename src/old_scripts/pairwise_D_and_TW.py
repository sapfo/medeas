# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 13:26:45 2016

@author: ivan
"""

import numpy as np
import matplotlib.pyplot as plt

from collections import Counter, OrderedDict

def lamb_theor(D, T, n, n1, n2):
    return 2*(2*n1*n2*D*(D+2) + n)/(n*T**2)

def D_theor(lamb, T, n, n1, n2):
    return (np.sqrt(1+n*(T**2*lamb/2 - 1)/(2*n1*n2)) - 1)

def D_theor_P(lamb, T, n, n1, n2):
    return n*(T*lamb - 1)/(2*n1*n2)

def mu(m,n):
    return (np.sqrt(n-1)+np.sqrt(m))**2 / n

def sigma(m,n):
    return (np.sqrt(n-1)+np.sqrt(m))*(1/np.sqrt(n-1)+1/np.sqrt(m))**(1/3)/n


def mark(ls,m):
    return (m+1)*np.sum(ls)**2/((m-1)*np.sum(ls**2)-np.sum(ls)**2)

def parse_file(fname):

    dists = np.loadtxt(fname+'.asd')
    
    with open(fname+'.asd.labels') as f:
        pops_lst = [line.strip()[:3] for line in f.readlines()]
    
#    for i,pop in enumerate(pops_lst):
#        if pop == 'Eng' or pop == 'Bri':
#            pops_lst[i] = 'Eur'
    
    pops_lst = np.array(pops_lst)
#    no_ind = np.where(pops_lst != 'Ind')[0]
    no_ind = np.where(np.logical_and(pops_lst != 'Ind', pops_lst != 'Eng'))[0]
    dists = dists[no_ind][:,no_ind]
    
    pops_lst = pops_lst[no_ind]
    
    pops = OrderedDict(sorted(Counter(pops_lst).items()))
    
    return dists, pops, pops_lst

def pair_d(pop1, pop2, pops_lst, dists):
    ind = np.where(np.logical_or(pops_lst == pop1, pops_lst == pop2))[0]
    return dists[ind][:,ind]
    
def find_eig(dst):

    n = len(dst)

    a = -dst**2/2
    at0 = np.sum(a, 0)/n
    att = np.sum(a)/n**2
 
    one = np.ones((n,))
    b = a - np.outer(at0, one) - np.outer(one, at0) + att
 
    lambdas, vecs = np.linalg.eig(b)
 
    arr = np.hstack((lambdas.reshape((n, 1)), vecs.T))
    arr = sorted(arr, key=lambda x: x[0], reverse=True)
    arr = arr[:3]
     
    for i, v in enumerate(arr):
        arr[i] = np.sqrt(v[0])*v[1:]
    
    return np.array(sorted(lambdas)[::-1]), arr[0], arr[1], arr[2]


def find_eig_P(dst):

    n = len(dst)

    a = -dst/2
    at0 = np.sum(a, 0)/n
    att = np.sum(a)/n**2
 
    one = np.ones((n,))
    b = a - np.outer(at0, one) - np.outer(one, at0) + att
 
    lambdas, vecs = np.linalg.eig(b)
 
    arr = np.hstack((lambdas.reshape((n, 1)), vecs.T))
    arr = sorted(arr, key=lambda x: x[0], reverse=True)
    arr = arr[:3]
     
    for i, v in enumerate(arr):
        arr[i] = np.sqrt(v[0])*v[1:]
    
    return np.array(sorted(lambdas)[::-1]), arr[0], arr[1], arr[2]


def TW_test_one(lambs, verbose=False):
    lambs = lambs[:-1]
    
    m = len(lambs)
    n = mark(lambs, m)
    MU = mu(m,n)
    S = sigma(m,n)
    l = (m-1)*lambs[0]/np.sum(lambs)
    if verbose: print('TW statistics ', (l-MU)/S)
    
    return (l-MU)/S

def find_D_T_TW(lambdas, lambdas_P, pop1, pop2, pops, verbose=False):
    
    n1 = pops[pop1]
    n2 = pops[pop2]
    n = n1+n2    
    
    T = np.sqrt(2/np.average(lambdas[1:-1]))
    D = D_theor(lambdas[0], T, n, n1, n2)
    
    T_P = 1/np.average(lambdas_P[1:-1])
    D_P = D_theor_P(lambdas_P[0], T_P, n, n1, n2)
    
    if verbose:
        print('T for pair {0} and {1} is {2}'.format(pop1, pop2, T))
        print('D for pair {0} and {1} is {2}'.format(pop1, pop2, D))
        print('T_PCA for pair {0} and {1} is {2}'.format(pop1, pop2, T_P))
        print('D_PCA for pair {0} and {1} is {2}'.format(pop1, pop2, D_P))
    if verbose: print('First eigenvalue ', end='')
    TW_test_one(lambdas, verbose)
    if verbose: print('Second eigenvalue ', end='')
    TW_test_one(lambdas[1:], verbose)
    if verbose: print('Third eigenvalue ', end='')
    TW_test_one(lambdas[2:], verbose)
    
    return D, T
    
def plot_lambdas(lambdas, T, pop1, pop2):

    plt.figure()
    plt.plot(lambdas, 'bo')
    plt.plot([2,len(lambdas)], [2/T**2]*2, 'r-', label='average small eigenvalues')
    plt.xlabel('number of eigenvalue')
    plt.ylabel('eigenvalue')
    plt.title(pop1+' and '+ pop2)
    plt.legend()
    plt.savefig('lambdas_{}{}.pdf'.format(pop1, pop2))

def plot_scatter(pop1, pop2, X, Y, pops):
    
    n1 = pops[pop1]
    
    plt.figure()
    plt.plot(X[:n1], Y[:n1], 'bo', label=pop1)
    plt.plot(X[n1:], Y[n1:], 'ro', label=pop2)
    plt.legend()
    print('saving...')
    plt.savefig('scatter_{}{}.pdf'.format(pop1, pop2))

def proc_file(fname, show=False):
    dists, pops, pops_lst = parse_file(fname)
    names = list(pops.keys())
    
    Ds = {}
    Ts = {}

    for i, p1 in enumerate(names):
        for p2 in names[i+1:]:
            matr = pair_d(p1, p2, pops_lst, dists)
            ls, X, Y, Z = find_eig(matr)
            ls_P, X_P, Y_P, Z_P = find_eig_P(matr)
            D, T = find_D_T_TW(ls, ls_P, p1, p2, pops, verbose=True)
            if show:
                plot_lambdas(ls, T, p1, p2)
                plot_scatter(p1, p2, X, Y, pops)

            Ds[p1+p2] = D
            Ts[p1+p2] = T
            
    
    return Ds, Ts

mind = 8

for nblocks in [1]: #,2,3,5,10,20,50]: #range(2,9,2):

    fbase = '../applications/asdfiles/yorgos_autosomes_pruned_mind0{}_indep5052_noPILCAIWPA_nblocks{}_wo_block'.format(mind, nblocks)
    #fname_m = '../applications/asdfiles/yorgos_autosomes_pruned_mind0{}_geno01'.format(mind)
    
    Ds, Ts = proc_file(fbase+str(0), show=True)
    
    Dsb = {k:[v] for k,v in Ds.items()}
    Tsb = {k:[v] for k,v in Ts.items()}
    
    for i in range(1, nblocks):
        fname = fbase + str(i)
        nD, nT = proc_file(fname, show=True)
        for pair in Ds:
            Dsb[pair].append(nD[pair])
            Tsb[pair].append(nT[pair])
    
    
    print('============= mind0{}nblocks{} =============\n'.format(mind, nblocks))

    outfile = '../applications/results_for_mind0{}_indep5052_nblocks{}.txt'.format(mind, nblocks) 
    
    with open(outfile, 'w') as outf:
    
        for pair in Ds:
            #print('The D for {} is {:5f}'.format(pair, Ds[pair]))
            #print('The T for {} is {:5f}'.format(pair, Ts[pair]))
            print('The average D for {} is {:.3f}, error {:.2f} %'.format(
                   pair, np.mean(Dsb[pair]), 100*np.sqrt((nblocks-1)*np.var(Dsb[pair]))/np.mean(Dsb[pair])),
                                 file=outf)
            print('The average T for {} is {:.3f}, error {:.2f} %'.format(
                   pair, np.mean(Tsb[pair]), 100*np.sqrt((nblocks-1)*np.var(Tsb[pair]))/np.mean(Tsb[pair])),
                                 file=outf)
        
            print('', file=outf)

#plt.show()
