# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 15:29:37 2016

@author: ivan
"""

TESTING = True

import sys
import numpy as np
from multiprocessing import Process, Queue, cpu_count

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

if __name__ == '__main__':

    # On Windows, processes execute the wholw file before forking
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
    #dist_func = np.square

    # ---------- constants

    N = 294
    sites = 24112 # TODO: read this from file
    cut = N*4
    name = sys.argv[1] #'C:\\Users\\levkivskyi\\PycharmProjects\\medeas\\test\\bla.tped'

    NPROC = cpu_count()

    # this should be large to avoid overhead of spawning new processes
    # or we need to reuse them somehow
    MAXSIZE = 100*2**20 # 100 MB

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
            data = np.array([np.fromstring(line[-cut:-1], sep=' ', dtype='int8')
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
    delta = np.sqrt(delta)

    print('woohoo!')
    if TESTING:
        plt.pcolor(delta)
        plt.colorbar()
        plt.xlabel('i individual')
        plt.ylabel('j individual')
        plt.title('Distance matrix')
        plt.show()
