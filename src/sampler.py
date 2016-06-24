# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 16:05:19 2016

@author: ivan
"""

import numpy as np
import matplotlib.pyplot as plt

L = 10000
n = 100
p = 0.05

numtests = 30

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

def process(illustrate=False):
    base1 = sample(L)
    base2 = sample(L)
       
    inds = make_pop(base1, p, n//2) + make_pop(base2, p, n//2)
    
    dists = np.zeros((n,n))
    dists2 = np.zeros((n,n))
    
    for i in range(n):
        for j in range(n):
            dists[i,j] = np.sum(np.abs(inds[i]-inds[j]))/L
            dists2[i,j] = np.sqrt(np.sum((inds[i]-inds[j])**2)/L)
    
    if illustrate:    
        plt.pcolor(dists)
        plt.colorbar()
    
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
    
    if illustrate:
        xx = np.linspace(min(x)*1.3, max(x)*1.3, 2500)
        yy1 = np.exp(-(xx-m1)**2/(2*v1**2))/(v1*np.sqrt(2*np.pi))
        yy2 = np.exp(-(xx-m2)**2/(2*v2**2))/(v2*np.sqrt(2*np.pi))
        
        plt.figure()
        plt.plot(x, [1]*n, 'bo')
        plt.plot(xx, yy1+yy2, 'r-')
    
    s2 = np.std(dists2[:n//2,n//2:])
    d2 = np.mean(dists2[:n//2,n//2:])
    s = np.std(dists[:n//2,n//2:])
    d = np.mean(dists[:n//2,n//2:])

    
    return v1/np.abs(m1 - m2), V1/np.abs(M1 - M2), s/d, s2/d2

rat_mat = []
rat_mds = []
for pr in range(numtests):
    s_m, s_m2, s_d, s_d2 = process(pr==numtests-1)
    print('---------------------------------------------------')
    print('sigma/mean for matrix elements (sqrt)   :', s_d2)
    print('sigma/mean for matrix elements (linear) :', s_d)
    print('sigma/distance for mds plot (sqrt)      :', s_m2)
    print('sigma/distance for mds plot (linear)    :', s_m)
    print('ratio of first two                      :', s_d2/s_d)
    print('ratio of last two                       :', s_m2/s_m)
    rat_mat.append(s_d2/s_d)
    rat_mds.append(s_m2/s_m)

print('---------------------------------------------------')
print('Grand mean for sqrt/linaer in mds       :', np.mean(rat_mds))
print('Grand sigma for sqrt/linaer in mds      :', np.std(rat_mds))
print('Grand mean for sqrt/linaer in matrix    :', np.mean(rat_mat))
print('Grand sigma for sqrt/linaer in matrix   :', np.std(rat_mat))


plt.show()
