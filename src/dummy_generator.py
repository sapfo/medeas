# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 14:52:13 2016

@author: ivan
"""

import random

name = 'test.what'
lines = 900
symbs = 600000

with open(name, 'w') as f:
    for i in range(lines):
        f.write(' '.join(map(str, (random.choice([0, 1, 2]) for j in range(symbs)))))
        f.write('\n')
        print(i/lines*100)
