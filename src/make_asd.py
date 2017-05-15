# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 15:29:37 2016

@author: ivan
"""

import sys
import numpy as np
from multiprocessing import Process, Queue, cpu_count

TESTING = True

if TESTING:
    import matplotlib.pyplot as plt

N = 294


def dist_and_norm(a, b, dist_func):
    filt = np.logical_and(a, b)
    dst = filt * dist_func(a-b)
    return np.sum(dst), np.sum(filt)


def compute(i, data, dist_func):
    dists = []
    norms = []
    for j in range(i+1, N):
        dist, norm = dist_and_norm(data[i], data[j], dist_func)
        dists.append(dist)
        norms.append(norm)
    return dists, norms


def work(tasks, results, data, dist_func):
    while True:
        i = tasks.get()
        if i < 0:
            return
        results.put((i, compute(i, data, dist_func)))


def process(data, dist_func):
    dists = np.zeros((N, N))
    norms = np.zeros((N, N))

    for i in range(N):
        tasks.put(i)
    for _ in range(NPROC):
        tasks.put(-1)

    procs = [Process(target=work, args=(tasks, results, data, dist_func))
             for _ in range(NPROC)]
    for proc in procs:
        proc.start()

    rest = N
    while rest:
        i, (dist, norm) = results.get()
        rest -= 1
        dists[i, i+1:] = dist
        norms[i, i+1:] = norm

    for proc in procs:
        proc.join()

    return dists, norms


def dens(x, a, b):
    print('arg1', (2*x - a - b)/(b - a))
    arc1 = np.arcsin((2*x - a - b)/(b - a))
    print('!!!!!!!res1', arc1)
    print('arg2', ((a+b)*x-2*a*b)/x/(b - a))
    arc2 = np.arcsin(((a+b)*x-2*a*b)/x/(b - a))
    print('!!!!!!!res2', arc2)
    res = (np.sqrt((b-x)*(x-a))
            + (a+b)*arc1/2
            - np.sqrt(a*b)*arc2
            + np.pi*((a+b)/2 - np.sqrt(a*b))/2)
    print('RESULT', res)
    return res

def dens_fit(x, T, L):
    a = 2/T**2 * (1 - np.sqrt(N/L))**2
    b = 2/T**2 * (1 + np.sqrt(N/L))**2
    print("CALLED:", T, L, "with", a, b, x)
    if x <= a+0.001:
        return 0.0
    if x >= b-0.001:
        return N*1.0
    return N * dens(x, a, b) / dens(b-0.001, a, b)


def find_TW(lambdas):
    m = len(lambdas)
    n = 1818
    mu = (np.sqrt(n - 1) + np.sqrt(m))**2 / n
    sigma = (1/np.sqrt(n-1) + 1/np.sqrt(m))**(1/3) * (np.sqrt(n - 1) + np.sqrt(m)) / n
    l = m*lambdas[0]/np.sum(lambdas)
    print(l, mu, sigma)
    s = (l - mu) / sigma
    return s

def find_TW_1(lambdas):
    m = len(lambdas)
    n = 2015
    mu = np.sqrt(2*m) / np.sqrt(n)
    sigma = 1/np.sqrt(2*m**(1/3)) / np.sqrt(n)
    l = m*lambdas[0]/np.sum(lambdas)
    print(l, mu, sigma)
    s = (l - mu) / sigma
    return s


def find_TW_2(lambdas):
    m = len(lambdas)
    n = 2015
    mu = (np.sqrt(n - 1) + np.sqrt(N))**2 / n
    sigma = (1/np.sqrt(n-1) + 1/np.sqrt(N))**(1/3) * (np.sqrt(n - 1) + np.sqrt(N)) / n
    l = m*lambdas[0]/np.sum(lambdas)
    print(l, mu, sigma)
    s = (l - mu) / sigma
    return s


if __name__ == '__main__':

    # On Windows, processes execute the whole file before forking
    # therefore we protect this code with if __name__ == '__main__'
    # Need to think how to avoid copying ``data`` on forking.
    # Maybe process input file in chunks?
    # On POSIX everything is already fine because of "copy-on-write"

    def test_func(dst):
        if dst > 1:
            return dst**2
        else:
            return 0

    dist_func = np.vectorize(test_func)
    dist_func = np.abs
    dist_func = np.square
    dist_func = lambda x: np.abs(x)**2

    # ---------- constants

    N = 294
    sites = 24112  # TODO: read this from file
    cut = N*4
    name = sys.argv[1]  # C:\Users\levkivskyi\PycharmProjects\medeas\test\...

    NPROC = cpu_count()

    # this should be large to avoid overhead of spawning new processes
    # or we need to reuse them somehow
    MAXSIZE = 100*2**20  # 100 MB

    # ---------- global data

    tot_dists = np.zeros((N, N))
    tot_norms = np.zeros((N, N))
    delta = np.zeros((N, N))
    tasks = Queue()
    results = Queue()

    with open(name) as f:
        while True:
            data_lines = f.readlines(MAXSIZE)
            if not data_lines:
                break
            data = np.array([np.fromstring(line[-cut:-1], sep=' ',
                                           dtype='int8')
                             for line in data_lines])
            data = data[:, ::2] + data[:, 1::2]
            data = data.T.copy()
            dists, norms = process(data, dist_func)
            tot_dists += dists
            tot_norms += norms

    for i in range(N):
        delta[i, i] = 0

    for i in range(N):
        for j in range(i+1, N):
            delta[i, j] = delta[j, i] = tot_dists[i, j]/tot_norms[i, j]
    delta = delta**.5

    print('woohoo!')
    if TESTING:
        plt.pcolor(delta)
        plt.colorbar()
        plt.xlabel('i individual')
        plt.ylabel('j individual')
        plt.title('Distance matrix')

        a = -delta**2/2

        at0 = np.sum(a, 0)/N
        att = np.sum(a)/N**2

        one = np.ones((N,))
        b = a - np.outer(at0, one) - np.outer(one, at0) + att

        lambdas, vecs = np.linalg.eig(b)


        # finding L and T from fit
        plt.figure()
        lambdas_s = np.array(sorted(lambdas))
        plt.plot(lambdas_s, range(len(lambdas_s)))
        lambdas_se = np.linspace(lambdas.min(), lambdas.max(), 5000)
        L = 1818  # Just an attempt
        t = 2.85
        a = 2/t**2 * (1 - np.sqrt(N/L))**2
        b = 2/t**2 * (1 + np.sqrt(N/L))**2
        pref = N/(2*np.pi*N/L)
        l_dens = [4*pref*dens(b, a, b) if l > b else
                  4*pref*dens(l, a, b) if l > a else 0 for l in lambdas_se]
        plt.plot(lambdas_se, l_dens)

        from scipy.optimize import curve_fit

        popt, pcov = curve_fit(np.vectorize(dens_fit, otypes=[np.float]),
                               lambdas_s, range(len(lambdas_s)),
                               p0 = (1, N+1), bounds=([0.1, N], [10, 100*N]))
        print('FIT: ', popt, pcov)
        #t, L = popt
        plt.figure()

        plt.plot(lambdas_s, range(len(lambdas_s)))
        l_dens_fit = [dens_fit(l, *popt) for l in lambdas_se]
        plt.plot(lambdas_se, l_dens_fit)


        arr = np.hstack((lambdas.reshape((N, 1)), vecs.T)).copy()
        arr = sorted(arr, key=lambda x: x[0], reverse=True)
        for i, v in enumerate(arr.copy()):
            arr[i] = np.sqrt(v[0])*v[1:]
        arr = arr[:14]
        arr = np.array(arr)
        arr = arr.T
        lambdas = list(sorted(lambdas, reverse=True))
        lambdas = np.array(lambdas)
        plt.figure()
        plt.plot(lambdas, 'b.')
        plt.figure()
        plt.hist(lambdas[5:], 50)
        #plt.figure()
        #plt.hist(np.diff(lambdas[8:])[1:], 50)

        s1 = np.sum(lambdas)
        s2 = np.sum(lambdas**2)

        print((N+1)*s1**2/((N-1)*s2 - s1**2))
        print(s1)
        print(s2)

        print('-'*20)
        s = 5
        i = 1
        while s > 3.2724:  # for p-value 0.001
            s = find_TW(lambdas)
            print(f'eigenvalue {i}: TW = {s}')
            i += 1
            lambdas = lambdas[1:]
        print('-'*20)
        print('shape', arr.shape)
        print('-'*20)
        lambdas1 = lambdas[5:]

        s1 = np.sum(lambdas1)
        s2 = np.sum(lambdas1**2)
        print((N+1)*s1**2/((N-1)*s2 - s1**2))


        flat = delta.reshape(N**2,)
        plt.figure()
        plt.hist(flat, 1250)
        npop =8
        from sklearn.cluster import AgglomerativeClustering as AC
        clusterer = AC(n_clusters=npop)
        labs = clusterer.fit_predict(arr)
        labels = [hex(l)[-1].upper() for l in labs]

        with open('../test/bla.tfam') as f:
            new_data = f.readlines()
        labels_0 = [l.split()[0].strip('"') for l in new_data]

        for p, q in [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)]:
            fig, ax = plt.subplots()
            ax.scatter(arr.T[p], arr.T[q])
            for i, txt in enumerate(labels):
                ax.annotate(txt, (arr.T[p, i], arr.T[q, i]))

        from collections import defaultdict
        dd = defaultdict(list)
        dd = {}
        ds = np.zeros((npop, npop))
        stds = np.zeros((npop, npop))

        ulabs = set(labs)
        small = lambdas[:-1]
        T = np.sqrt(2/np.average(small))

        for i in ulabs:
            for j in ulabs:
                i_s = np.where(labs == i)
                j_s = np.where(labs == j)
                block = delta[i_s].T[j_s]*T/2
                block = block.reshape(block.shape[0]*block.shape[1],)
                if i == j:
                    block = [el for el in block if el > 1e-5]
                ds[i, j] = np.average(block)
                stds[i, j] = np.std(block)
                dd[i, j] = block

        plt.figure()
        plt.pcolor(ds)
        plt.colorbar()
        plt.figure()
        plt.pcolor(stds)
        plt.colorbar()
        #from pprint import pprint
        #pprint(dd)
        #pprint(ds)
        plt.show()


def inv_filter(p1, p2, p3): # This should be moved to other place actually
    M = np.eye(3)
    for i in range(3):
        M[i, 0] += 2*p1*p2
        M[i, 1] += 2*p1*p3
        M[i, 2] += 2*p2*p3
    M[0, 0] -= p1+p2
    M[1, 1] -= p1+p3
    M[2, 2] -= p2+p3
    M[0, 1] -= p3
    M[0, 2] -= p3
    M[1, 0] -= p2
    M[1, 2] -= p2
    M[2, 0] -= p1
    M[2, 1] -= p1
    return np.matrix(M)**-1