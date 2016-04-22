# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 16:05:19 2016

@author: ivan
"""

import numpy as np
import matplotlib.pyplot as plt

L = 100000
n = 100
p = 0.05

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
    
def sample(L):
    return np.random.randint(0,2,L)
    
base1 = sample(L)
base2 = sample(L)

def make_pop(base, p, size):
    result = []
    for _i in range(size):
        ind = np.copy(base)
        for i, elem in enumerate(ind):
            if np.random.random() < p:
                ind[i] = 1 - elem
        
        result.append(ind)
        
    return result

def process():
    base1 = sample(L)
    base2 = sample(L)
       
    inds = make_pop(base1, p, n//2) + make_pop(base2, p, n//2)
    
    dists = np.zeros((n,n))
    dists2 = np.zeros((n,n))
    
    for i in range(n):
        for j in range(n):
            dists[i,j] = np.sum(np.abs(inds[i]-inds[j]))/L
            dists2[i,j] = np.sqrt(np.sum((inds[i]-inds[j])**2)/L)
            
    #plt.pcolor(dists)
    #plt.colorbar()
    
    ls, x, Y, Z = find_eig(dists)
    ls, X, Y, Z = find_eig(dists2)
    
    arr1 = x[:n//2]
    arr2 = x[n//2:]
    m1 = np.mean(arr1)
    m2 = np.mean(arr2)
    v1 = np.std(arr1)
    v2 = np.std(arr2)
    
    arr1 = X[:n//2]
    arr2 = X[n//2:]
    M1 = np.mean(arr1)
    M2 = np.mean(arr2)
    V1 = np.std(arr1)
    V2 = np.std(arr2)
    
    #xx = np.linspace(min(X), max(X), 300)
    #yy1 = np.exp(-(xx-m1)**2/(2*v1**2))/(v1*np.sqrt(2*np.pi))
    #yy2 = np.exp(-(xx-m2)**2/(2*v2**2))/(v2*np.sqrt(2*np.pi))
    
    #plt.figure()
    #plt.plot(X, [0.01]*n, 'bo')
    #plt.plot(xx, yy1+yy2, 'r-')
    
    return (m1 - m2)/v1, (M1 - M2)/V2

for pr in range(10):
    print(process())

#plt.show()
