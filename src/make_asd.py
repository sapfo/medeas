# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 15:29:37 2016

@author: ivan
"""

import numpy as np
import time

start = time.time()

name = 'test.what'

lines = 300
symbs = 300000

f = open(name, 'r')

delta = np.zeros((lines, lines))
for i in range(lines):
    f.seek(i*2*symbs)
    a = np.fromstring(f.readline(), sep=' ', dtype='int8')
    for j in range(i+1, lines):
        f.seek(j*2*symbs)
        b = np.fromstring(f.readline(), sep=' ', dtype='int8')
        delta[i,j] = np.sqrt(np.sum((a-b)**2)/lines)

    print(i)

f.close()
delta = delta + delta.T

finish = time.time()

print(finish-start)