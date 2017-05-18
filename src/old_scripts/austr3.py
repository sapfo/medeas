import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from collections import Counter, OrderedDict

c = 'k 0.4 r c m g b y 0.8'.split()
#c = 'g g r c m g b y 0.8'.split()

# ========== parsing file =================

with open('/home/ivan/aboriginals.ibs') as f:
    data = f.readlines()

data = [line.strip().split() for line in data]

good = ['WCD', 'PIL', 'NGA', 'ENY']
good_pos = []
for i, pop in enumerate(np.array(data).T[0]):
    if pop in good:
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

T = np.sqrt(2/np.average(lambdas[4:]))

# ========== 2D plots ====================
plt.pcolor(dists)
plt.colorbar()
plt.xlim(0,45)
plt.ylim(0,45)
plt.xlabel('individuum numner')
plt.ylabel('individuum numner')

dpi = 80
fig = plt.figure(figsize=(800/dpi, 600/dpi), dpi=dpi)
plt.plot(list(reversed(sorted(lambdas))), 'bo')
plt.plot([4,45], [2/T**2]*2, 'r-', label='average small eigenvalues')
plt.xlabel('number of eigenvalue')
plt.ylabel('eigenvalue')
plt.legend()

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

arr1 = arr[0][np.where(arr[0]>0.04)]
arr2 = arr[0][np.where(arr[0]<0.04)]

m1 = np.mean(arr1)
m2 = np.mean(arr2)

v1 = np.std(arr1)
v2 = np.std(arr2)

xx = np.linspace(min(arr[0]), max(arr[0]), 100)
yy1 = np.exp(-(xx-m1)**2/(2*v1**2))/(v1*np.sqrt(2*np.pi))
yy2 = np.exp(-(xx-m2)**2/(2*v2**2))/(v2*np.sqrt(2*np.pi))

plt.figure()
plt.plot(arr[0], [0.01]*len(arr[0]), 'bo')
plt.plot(xx, yy1+yy2, 'r-')

# ========== significance test ===========

def mu(m,n):
    return (np.sqrt(n-1)+np.sqrt(m))**2 / n

def sigma(m,n):
    return (np.sqrt(n-1)+np.sqrt(m))*(1/np.sqrt(n-1)+1/np.sqrt(m))**(1/3)/n


def mark(ls,m):
    return (m+1)*np.sum(ls)**2/((m-1)*np.sum(ls**2)-np.sum(ls)**2)


lambdas = np.array(sorted(lambdas)[::-1])
#lambdas = lambdas[:-2]

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
