# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 15:29:37 2016

@author: ivan
"""

import numpy as np
from multiprocessing import Process, Queue, current_process, cpu_count
from queue import Empty

from time import time, sleep
start = time()

N = 294
sites = 24112
cut = N*4 - 1
name = '../test/bla.tped'

NPROC = cpu_count()
dist_func = np.abs

def test_func(dst):
    if dst > 1:
        return dst**2
    else:
        return 0

#dist_func = np.vectorize(test_func)

def dist(a, b, dist_func):
    filt = np.logical_and(a, b)
    dst = filt * dist_func(a-b)
    return np.sum(dst)/np.sum(filt)

with open(name) as f:
    data = f.read().splitlines()

data = np.array([np.fromstring(l[-cut:], sep=' ', dtype='int8') for l in data])

data = data[:,::2] + data[:,1::2]

delta = np.zeros((N, N))

tasks = Queue()
for i in range(N):
    tasks.put(i)

for _ in range(NPROC):
    tasks.put(-1)

results = Queue()

def compute(i):
    print(current_process().name + ': ' + str(i))
    res = []
    for j in range(i+1, N):
        res.append(dist(data[:,i], data[:,j], dist_func))
    return res

def work(t, r):
    while True:
        i = t.get()
        if i < 0:
            print(current_process().name + ' exiting on empty')
            return
        r.put([i, compute(i)])

procs = [Process(target=work, args=(tasks, results)) for _ in range(NPROC)]

for proc in procs:
    proc.start()

rest = N
while rest:
    i, res = results.get()
    rest -= 1
    delta[i, i+1:] = res

for proc in procs:
    proc.join()

delta = delta + delta.T
print(time()-start)
