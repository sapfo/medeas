import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from collections import Counter, OrderedDict

c = 'k 0.4 r c m g b y 0.8'.split()

# ========== parsing file =================

with open('/home/ivan/aboriginals.ibs') as f:
    data = f.readlines()

data = [line.strip().split() for line in data]

pops = [line[0] for line in data]
pops = OrderedDict(sorted(Counter(pops).items()))

names = list(pops.keys())
sizes = list(pops.values())
starts = [0] + list(np.cumsum(sizes))

dists = [list(map(float, line[2:])) for line in data]
dists = np.array(dists)
n = len(dists)

# ========== finding eigensystem ==========

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

# ========== 2D plots ====================
dpi = 80
fig = plt.figure(figsize=(800/dpi, 600/dpi), dpi=dpi)
plt.plot(list(reversed(sorted(lambdas))), 'bo')
plt.xlabel('number of eigenvalue')
plt.ylabel('eigenvalue')

dpi = 60
f, ax = plt.subplots(4, 4, sharex='col', sharey='row', 
                     figsize=(1400/dpi, 1000/dpi), dpi=dpi)
for p in range(4):
    for q in range(4):
        for i,s in enumerate(sizes):
            st = starts[i]
            ax[q][p].plot(arr[p][st:st+s], arr[q][st:st+s],marker='o', color=c[i],
                     linewidth=0, label = names[i])
            if q == 3: ax[q][p].set_xlabel('dimension ' + str(p+1))
            if p == 0: ax[q][p].set_ylabel('dimension ' + str(q+1))

f.subplots_adjust(hspace=0.1, wspace=0.1)
f.tight_layout()
ax[0][3].legend()

# ========== 3D plot =====================
dpi = 80
fig = plt.figure(figsize=(800/dpi, 600/dpi), dpi=dpi)
ax = fig.gca(projection='3d')
for i,s in enumerate(sizes):
    st = starts[i]
    ax.scatter(arr[0][st:st+s], arr[1][st:st+s], arr[2][st:st+s],
               marker='o', color=c[i],
             label = names[i], alpha = 0.9)
plt.xlabel('dimension 1')
plt.ylabel('dimension 2')
ax.set_zlabel('dimension 3')
plt.legend(loc=2)

plt.show()