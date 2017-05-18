# -*- coding: utf-8 -*-
"""
Created on Mon May  4 12:22:31 2015

@author: ivan
"""

import numpy as np
import matplotlib.pyplot as pl
from textwrap import dedent
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.legend_handler import HandlerLine2D
import random

T = 2
RANDOM = 3

TEST = '4pops-grad'
THEORY = True

# --------------- setup ---------------

PRESETS = {

'2pops'     : ([20, 40],             # population sizes
               ['b', 'g'],           # colors for plotting
               '''
               d[0, 1] = 8
               '''),                 # distances

'3pops'     : ([50, 30, 20],
               ['b', 'g', 'r'],
               '''
               d[0, 1] = 3
               d[0, 2] = d[1, 2] = 8
               '''),

'4pops-pair' : ([110, 20, 30, 40],
                ['b', 'g', 'r', 'c'],
                '''
                d[0, 1] = 3
                d[2, 3] = 5
                d[0, 2] = d[0, 3] = d[1, 2] = d[1, 3] = 8
                '''),

'4pops-grad' : ([110, 20, 30, 40],
                ['b', 'g', 'r', 'c'],
                '''
                d[0, 1] = 3
                d[0, 2] = d[1, 2] = 5
                d[0, 3] = d[1, 3] = d[2, 3] = 8
                '''),

}

ns, cs, ds = PRESETS[TEST]
N = len(ns)                        # number of pops
n = sum(ns)                        # total number of individs

d = np.zeros((N, N))
exec(dedent(ds), globals())

for p in range(N):
    for q in range(p, N):
        d[q, p] = d[p, q]            # distances must be symmetric

# --------------- theoretical results ---------------

def vecs2(n1, n2, D, n=n):
    '''
    Calculate the largest eigenvector for two populations
    '''
    lambda_1 = (2/n/T**2)*(2*n1*n2*D*(D+2) + n)

    x = np.sqrt(n2/n/n1)
    y = np.sqrt(n1/n/n2)

    return [-np.sqrt(lambda_1)*x, np.sqrt(lambda_1)*y]

def vecs3(n1, n2, n3, D, D1, n=n):
    '''
    Calculate two largest eigenvectors for three populations
    '''
    dt1 = D1*(D1+2)
    dt = D*(D+2)

    R = np.sqrt(dt**2    * (n1+n2)**2*n3**2
               +dt1**2   * n1*n2*(n1+n3)*(n2+n3)
               -2*dt*dt1 * n1*n2*n3*(n+n3))

    lambdat_p = -(dt1*n1*n2 + dt*(n1+n2)*n3 + R)/n
    lambdat_m = -(dt1*n1*n2 + dt*(n1+n2)*n3 - R)/n
    lambda_p =   (2/T**2)*(1 - lambdat_p)
    lambda_m =   (2/T**2)*(1 - lambdat_m)

    xp = -(-dt*(n1+n2)*n3 + dt1*(n2+n3)*n1 + R)/(dt1*(n1-n2)*n1)
    yp =  (-dt*(n1+n2)*n3 + dt1*(n1+n3)*n2 + R)/(dt1*(n1-n2)*n2)
    xm =  ( dt*(n1+n2)*n3 - dt1*(n2+n3)*n1 + R)/(dt1*(n1-n2)*n1)
    ym = -( dt*(n1+n2)*n3 - dt1*(n1+n3)*n2 + R)/(dt1*(n1-n2)*n2)

    norm_p = np.sqrt(lambda_p)/np.sqrt(n1*xp**2 + n2*yp**2 + n3)
    X = [xp*norm_p, yp*norm_p, norm_p]
    norm_m = np.sqrt(lambda_m)/np.sqrt(n1*xm**2 + n2*ym**2 + n3)
    Y = [xm*norm_m, ym*norm_m, norm_m]

    return X, Y

def vecsn(ns=ns, d=d):
    '''
    Calculate three largest eigenvectors for the general case
    '''
    n = sum(ns)
    N = len(ns)

    a = (d+1)**2

    one = np.ones((N,))
    an = np.outer(ns, one)

    at = (np.sum(an*a, 0)-1)/n
    att = (np.sum(an*an.T*a) -n)/n**2
    b = a - np.outer(at, one) - np.outer(one, at) + att

    lambdas, vecs = np.linalg.eig(b*np.outer(one, ns))
    lambdas = (2/T**2)*(1-lambdas)

    arr = np.hstack((lambdas.reshape((N, 1)), vecs.T))
    arr = sorted(arr, key=lambda x: x[0], reverse=True)
    arr = arr[:3]

    for i, v in enumerate(arr):
        arr[i] = np.sqrt(v[0])*v[1:]/np.sqrt(sum(v[1:]**2*ns))

    return tuple(arr)

# --------------- dispatching ---------------

if N == 2:
    X = vecs2(ns[0], ns[1], d[0, 1])
    Y = [0, 0]
    Z = [0, 0]

elif N == 3:
    X, Y = vecs3(ns[0], ns[1], ns[2], d[1, 2], d[0, 1])
    Z = [0, 0, 0]

else:
    X, Y, Z = vecsn()

# --------------- making matrix --------------- < SIMULATION STARTS HERE

delta = np.zeros((n, n))
for p in range(N):
    for q in range(N):
        for i in range(ns[p]):
            for j in range(ns[q]):
                ii = i + sum(ns[:p], 0)
                jj = j + sum(ns[:q], 0)
                if p == q:
                    if i != j:
                        delta[ii, jj] = 1
                else:
                    delta[ii, jj] = d[p, q] + 1

if RANDOM:
    for i in range(n):
        for j in range(i+1, n):
            dst = RANDOM*random.random()
            delta[i, j] += dst # add small random distance between idivids
            delta[j, i] += dst

delta = delta*2/T

# --------------- shuffling individs ---------------

perm = list(range(n))
random.shuffle(perm) # make a random permutation

iperm = [perm.index(i) for i in range(n)] # and its inverse

old_delta = delta.copy()
delta = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        delta[i, j] = old_delta[perm[i], perm[j]]

# --------------- finding the largest eigenvectors ---------------

a = -delta**2/2

at0 = np.sum(a, 0)/n
att = np.sum(a)/n**2

one = np.ones((n,))
b = a - np.outer(at0, one) - np.outer(one, at0) + att

lambdas, vecs = np.linalg.eig(b)

arr = np.hstack((lambdas.reshape((n, 1)), vecs.T))
arr = sorted(arr, key=lambda x: x[0], reverse=True)
arr = arr[:3] # take three largest eigenvalues and their vectors

for i, v in enumerate(arr):
    arr[i] = np.sqrt(v[0])*v[1:]

arr = np.array(arr).T
arr = [arr[iperm[i]] for i in range(n)] # restore order for coloring
arr = np.array(arr).T

# --------------- plotting --------------- < SIMULATION ENDS HERE

pl.pcolor(delta, cmap='gray')
pl.colorbar()
pl.title(r'Dissimilarity matrix $\delta_{rs}$', size=14)
pl.xlabel(r'$r$', size=14)
pl.ylabel(r'$s$', size=14)
pl.savefig('dissim-'+TEST+'.pdf')

csign = lambda arr: np.copysign(1, arr)

fig = pl.figure()
ax = fig.gca(projection='3d')
for p in range(N):
    s = slice(sum(ns[:p], 0), sum(ns[:p], 0) + ns[p])
    Xs = arr[0][s]
    Ys = arr[1][s]
    Zs = arr[2][s]
    Xs = Xs * csign(X[p]) * csign(sum(Xs)) # syncronize projection side
    Ys = Ys * csign(Y[p]) * csign(sum(Ys)) # for theory and simulation
    Zs = Zs * csign(Z[p]) * csign(sum(Zs))
    ax.scatter(Xs, Ys, Zs, marker='o', color=cs[p], alpha=0.7)
    if THEORY:
        ax.scatter(X[p], Y[p], Z[p], marker='D', color='k')
        # theoretical data always shown by black diamonds

box = [np.min([X, Y, Z])-1, np.max([X, Y, Z])+1]
ax.set_xlim(box)
ax.set_ylim(box)
ax.set_zlim(box)

pl.figure()
lines = list(range(N+THEORY))
for p in range(N):
    s = slice(sum(ns[:p], 0), sum(ns[:p], 0) + ns[p])
    Xs = arr[0][s]
    Ys = arr[1][s]
    Xs = Xs * csign(X[p]) * csign(sum(Xs))
    Ys = Ys * csign(Y[p]) * csign(sum(Ys))
    lines[p], = pl.plot(Xs, Ys, 'o', color=cs[p], linewidth=0, alpha=0.7,
                     label=r'population #{:d}'.format(p+1))
    if THEORY:
        lines[N], = pl.plot(X[p], Y[p], 'D', color='k', markersize=8,
                     label=r'theory' if p == N-1 else None)
# a nicer legend
hmap = {line: HandlerLine2D(numpoints=1) for line in lines}
pl.legend(handler_map=hmap, loc='lower left')

pl.xlim([min(X)-RANDOM, max(X)+RANDOM])
pl.ylim([min(Y)-RANDOM, max(Y)+RANDOM])
pl.savefig('2D-'+TEST+'.pdf')

# --------------- selscted set of purely theoretical plots ---------------

if THEORY:
    if N == 2:
        pl.figure()
        D = np.linspace(0.1, 8, 100)
        for n12 in np.linspace(1/n, 1/2, 5):
            pl.plot(D, vecs2(n*n12, n*(1-n12), D)[1]
                     - vecs2(n*n12, n*(1-n12), D)[0],
                     label=r'$n_1/n = $ {:.4f}'.format(n12))
        pl.xlabel(r'real distance, $D$', size=14)
        pl.ylabel(r'inferred distance, $\delta_{12}$', size=14)
        pl.legend()
        pl.savefig('plot1-'+TEST+'.pdf')

        pl.figure()
        D = 8
        n12 = np.linspace(1/n, 1-1/n, 100)
        for n in range(1, 15, 3):
            pl.plot(n12, vecs2(n*n12, n*(1-n12), D, n)[1]
                       - vecs2(n*n12, n*(1-n12), D, n)[0],
                     label=r'$n = $ {:.1f}'.format(n))
        pl.xlabel(r'relative size of first population, $n_1/n$', size=14)
        pl.ylabel(r'inferred distance, $\delta_{12}$', size=14)
        pl.legend()
        pl.savefig('plot2-'+TEST+'.pdf')

        pl.figure()
        n12 = 0.3
        n = np.linspace(1, 20, 100)
        for D in range(1, 8, 2):
            pl.plot(n, vecs2(n*n12, n*(1-n12), D, n)[1]
                     - vecs2(n*n12, n*(1-n12), D, n)[0],
                     label=r'$D = $ {:.1f}'.format(D))
        pl.xlabel(r'total number of individs, $n$', size=14)
        pl.ylabel(r'inferred distance, $\delta_{12}$', size=14)
        pl.legend()
        pl.savefig('plot3-'+TEST+'.pdf')

    elif N == 3:
        n1n = np.linspace(1/n+0.01, 1-1/n-0.015, 100)
        n12 = np.linspace(1/n+0.015, 1-1/n-0.01, 100)
        n1n, n12 = np.meshgrid(n1n, n12)
        n1 = n1n*n
        n2 = n12*(1-n1n)*n
        n3 = n - n1 - n2

        vx, vy = vecs3(n1, n2, n3, d[1, 2], d[0, 1])
        dist01 = np.sqrt((vx[0] - vx[1])**2 + (vy[0] - vy[1])**2)

        fig = pl.figure()
        ax = fig.gca(projection='3d')
        ax.plot_surface(n1, n2, dist01, cmap='coolwarm',
                        rstride=1, cstride=1, linewidth=0)
        ax.set_xlabel(r'population 1 size, $n_{1}$')
        ax.set_ylabel(r'population 2 size, $n_{2}$')
        ax.set_zlabel(r'inferred distance $\delta_{12}$')

        n1n = np.linspace(1/n+0.01, 1/2-1/n-0.015, 300)
        n2n = np.linspace(1/n+0.015, 1/2-1/n-0.01, 30)
        n1n, n2n = np.meshgrid(n1n, n2n)

        vx, vy = vecs3(n1n*n, n2n*n, (1-n1n-n2n)*n, d[1, 2], d[0, 1])
        dist02 = np.sqrt((vx[0] - vx[2])**2 + (vy[0] - vy[2])**2)

        fig = pl.figure()
        ax = fig.gca(projection='3d')
        ax.plot_surface(n1n, n2n, dist02, cmap='coolwarm',
                        rstride=1, cstride=1, linewidth=0)
        ax.set_xlabel(r'population 1 ratio, $n_{1}/n$')
        ax.set_ylabel(r'population 2 ratio, $n_{2}/n$')
        ax.set_zlabel(r'inferred distance $\delta_{13}$')

    else:
        n1ns = np.linspace(1/n+0.01, 1/4-1/n-0.015, 50)
        n2ns = np.linspace(1/n+0.015, 1/4-1/n-0.01, 50)
        n1np, n2np = np.meshgrid(n1ns, n2ns)
        #sorry for this, it is not easy to bradcast through genric function
        dist01 = np.zeros((len(n1ns), len(n2ns)))
        dist02 = np.zeros((len(n1ns), len(n2ns)))
        for i, n1n in enumerate(n1ns):
            for j, n2n in enumerate(n2ns):
                vx, vy, vz = vecsn([n1n*n, n2n*n, 50, (1-n1n-n2n)*n-50])
                dist01[i,j] = np.sqrt((vx[0] - vx[1])**2 
                                    + (vy[0] - vy[1])**2)
                dist02[i,j] = np.sqrt((vx[0] - vx[2])**2 
                                    + (vy[0] - vy[2])**2)

        fig = pl.figure()
        ax = fig.gca(projection='3d')
        ax.plot_surface(n1np, n2np, dist01, cmap='coolwarm',
                        rstride=1, cstride=1, linewidth=0)
        ax.set_xlabel(r'population 1 ratio, $n_{1}/n$')
        ax.set_ylabel(r'population 2 ratio, $n_{2}/n$')
        ax.set_zlabel(r'inferred 2D distance $\delta_{12}$')

        fig = pl.figure()
        ax = fig.gca(projection='3d')
        ax.plot_surface(n1np, n2np, dist02, cmap='coolwarm',
                        rstride=1, cstride=1, linewidth=0)
        ax.set_xlabel(r'population 1 ratio, $n_{1}/n$')
        ax.set_ylabel(r'population 2 ratio, $n_{2}/n$')
        ax.set_zlabel(r'inferred 2D distance $\delta_{13}$')        

pl.show()
