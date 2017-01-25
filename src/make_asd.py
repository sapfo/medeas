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
    dst = np.sign(a)**2 * np.sign(b)**2 * np.abs(a-b)
    norm = np.sign(a)**2 * np.sign(b)**2
    return np.sum(dst)/np.sum(norm)

with open(name) as f:
    data = f.read().splitlines()

data = np.array([np.fromstring(l[-cut:], sep=' ', dtype='int8') for l in data])
data = data.T

finish = time.time()
print(finish-start)

delta = np.zeros((N, N))

for i in range(N):
    print(i)
    for j in range(i+1, N):
        delta[i,j] = fun(data[i], data[j])

delta = delta + delta.T

finish = time.time()
print(finish-start)