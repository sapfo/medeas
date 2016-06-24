# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 15:29:37 2016

@author: ivan
"""

import numpy as np
import time

start = time.time()

name = 'test.what'

def fun(a, b):
    return (a - b)**2

lines = 300
symbs = 300000

with open(name) as f:
    data = np.fromstring(f.read(), sep=' ', dtype='int8').reshape((lines, symbs))

finish = time.time()
print(finish-start)

delta = np.zeros((lines, lines))

for i in range(lines):
    print(i)
    for j in range(i+1, lines):
        delta[i,j] = np.sum(fun(data[i], data[j]))

delta = delta + delta.T

finish = time.time()
print(finish-start)