# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 15:29:37 2016

@author: ivan
"""

import numpy as np
import time

start = time.time()

N = 294
sites = 24112
cut = N*4 - 1
name = '../test/bla.tped'

def fun(a, b):
    filt = np.logical_and(a, b)
    dst = filt * np.abs(a-b)
    norm = np.count_nonzero(filt)
    return np.sum(dst)/norm

with open(name) as f:
    data = f.read().splitlines()

data = np.array([np.fromstring(l[-cut:], sep=' ', dtype='int8') for l in data])

finish = time.time()
print(finish-start)

delta = np.zeros((2*N, 2*N))

for i in range(2*N):
    print(i)
    for j in range(i+1, 2*N):
        delta[i,j] = fun(data[:,i], data[:,j])

delta = delta + delta.T

finish = time.time()
print(finish-start)
