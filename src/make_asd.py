# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 15:29:37 2016

@author: ivan
"""

import numpy as np
from multiprocessing import Process, Queue, cpu_count

# ---------- constants

N = 294
sites = 24112
cut = N*4
name = '../test/bla.tped'

NPROC = cpu_count()

# this should be large to avoid overhead of spawning new processes
# or we need to reuse them somehow
MAXSIZE = 100*2**20 # 100 MB

# ---------- distance functions

def test_func(dst):
    if dst > 1:
        return dst**2
    else:
        return 0

dist_func = np.vectorize(test_func)
dist_func = np.abs
dist_func = np.square

# ---------- global data

tot_dists = np.zeros((N, N))
tot_norms = np.zeros((N, N))
tasks = Queue()
results = Queue()

def dist_and_norm(a, b, dist_func):
    filt = np.logical_and(a, b)
    dst = filt * dist_func(a-b)
    return np.sum(dst), np.sum(filt)

def compute(i):
    dists = []
    norms = []
    for j in range(i+1, N):
        dist, norm = dist_and_norm(data[i], data[j], dist_func)
        dists.append(dist)
        norms.append(norm)
    return dists, norms

def work(tasks, results):
    while True:
        i = tasks.get()
        if i < 0:
            return
        results.put((i, compute(i)))

def process():
    dists = np.zeros((N, N))
    norms = np.zeros((N, N))

    for i in range(N):
        tasks.put(i)
    for _ in range(NPROC):
        tasks.put(-1)

    procs = [Process(target=work, args=(tasks, results))
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

# ---------- main part

with open(name) as f:
    while True:
        data_lines = f.readlines(MAXSIZE)
        if not data_lines:
            break
        data = np.array([np.fromstring(line[-cut:-1], sep=' ', dtype='int8')
                         for line in data_lines])
        data = data[:, ::2] + data[:, 1::2]
        data = data.T.copy()
        dists, norms = process() # TODO: avoid being beaten :)
        tot_dists += dists
        tot_norms += norms

delta = (tot_dists + 1e-16)/(tot_norms + 1e-8)  # TODO: do this cleanly
for i in range(N):
    delta[i, i] = 0

delta = delta + delta.T
delta = np.sqrt(delta)
