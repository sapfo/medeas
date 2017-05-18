# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 21:25:12 2015

@author: ivan
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Ellipse

from collections import Counter, OrderedDict

c = 'k 0.4 r c m g b y 0.8'.split()
#c = 'g g r c m g b y 0.8'.split()

# ========== parsing file =================

with open('aboriginals.ibs') as f:
    data = f.readlines()

data = [line.strip().split() for line in data]

good = ['WCD', 'PIL', 'NGA', 'ENY']
bad = ['WON']
good_pos = []
for i, pop in enumerate(np.array(data).T[0]):
    if pop not in bad:
    #if pop in good:
        good_pos.append(i)

good_pos_T = [0,1] + list(np.array(good_pos)+2)

good_pos = np.array(good_pos)
good_pos_T = np.array(good_pos_T)

#data = np.array(data)[good_pos].T[good_pos_T].T

pops = [line[0] for line in data]
pops = OrderedDict(sorted(Counter(pops).items()))

names = list(pops.keys())
sizes = list(pops.values())
starts = [0] + list(np.cumsum(sizes))

dists = [list(map(float, line[2:])) for line in data]
dists = np.array(dists)

n = len(dists)

# ========== finding eigensystem ==========

cc = 0
cc = 0.19
dists = dists*(1-cc) + cc*np.ones((n,n))

a = -dists**2/2

at0 = np.sum(a, 0)/n
att = np.sum(a)/n**2
 
one = np.ones((n,))
b = a - np.outer(at0, one) - np.outer(one, at0) + att
 
lambdas, vecs = np.linalg.eig(b)
 
arr = np.hstack((lambdas.reshape((n, 1)), vecs.T))
arr = sorted(arr, key=lambda x: x[0], reverse=True)
arr = arr[:4] # take four largest eigenvalues and their vectors
 
for i, v in enumerate(arr):
    arr[i] = np.sqrt(v[0])*v[1:]

T = np.sqrt(2/np.average(lambdas[2:-1]))
print('T = ', T)

# ========== testing clusters ============

ind_c = []

my_c = 'b m k'.split()

ind_c.append(np.where(arr[0]>0.05))
ind_c.append(np.where(arr[0]<0.05))
#ind_c.append(np.where((arr[0]<0.05) & (arr[1]>-0.05)))

arr0 = {}
arr1 = {}
av0 = {}
av1 = {}
cv = {}
st0 = {}
st1 = {}
al = {}

for p in range(2):
    arr0[p]= arr[0][ind_c[p]]
    arr1[p] = arr[1][ind_c[p]]
    av0[p] = np.average(arr0[p])
    av1[p] = np.average(arr1[p])
    cv[p] = np.cov(arr0[p], arr1[p])
    ls, vs = np.linalg.eig(cv[p])
    st0[p], st1[p] = np.sqrt([max(ls), min(ls)])
    maxv = vs[np.where(ls == max(ls))][0]
    al[p] = - np.arctan2(maxv[1], maxv[0])*180/np.pi

fig = plt.figure()
ax1 = fig.add_subplot(111, aspect='equal')

for i,s in enumerate(sizes):
    st = starts[i]
    ax1.plot(arr[0][st:st+s], arr[1][st:st+s],marker='o', color=c[i],
                  linewidth=0, label = names[i])

sig2 = 2*np.sqrt(5.991)
sig3 = 2*np.sqrt(9.210)

ells = [Ellipse(xy=[av0[p], av1[p]], width=sig2*st0[p], 
                height=sig2*st1[p], angle=al[p])
        for p in range(2)]

ells2 = [Ellipse(xy=[av0[p], av1[p]], width=sig3*st0[p], 
                 height=sig3*st1[p], angle = al[p])
        for p in range(2)]

for i,e in enumerate(ells+ells2):
    ax1.add_artist(e)
    e.set_clip_box(ax1.bbox)
    e.set_alpha(0.3)
    e.set_facecolor(my_c[i % 3])

plt.xlabel('dimension 1')
plt.ylabel('dimension 2')

# ========== significance test ===========

def mu(m,n):
    return (np.sqrt(n-1)+np.sqrt(m))**2 / n

def sigma(m,n):
    return (np.sqrt(n-1)+np.sqrt(m))*(1/np.sqrt(n-1)+1/np.sqrt(m))**(1/3)/n


def mark(ls,m):
    return (m+1)*np.sum(ls)**2/((m-1)*np.sum(ls**2)-np.sum(ls)**2)


lambdas = np.array(sorted(lambdas)[::-1])
lambdas = lambdas[:-1]

m = len(lambdas)
n = mark(lambdas, m)
MU = mu(m,n)
S = sigma(m,n)
l = (m-1)*lambdas[0]/np.sum(lambdas)
print((l-MU)/S)

lambdas1 = lambdas[1:]
m = len(lambdas1)
n = mark(lambdas1, m)
MU = mu(m,n)
S = sigma(m,n)
l = (m-1)*lambdas1[0]/np.sum(lambdas1)
print((l-MU)/S)

plt.show()
