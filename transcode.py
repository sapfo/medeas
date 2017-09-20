#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 12:11:31 2017

@author: ivan
"""

import sys

n = 100

with open(sys.argv[1]) as f:
    data = f.readlines()

data = [d.strip() for d in data]
data = [d for d in data[1:] if d]

tot_blocks = len(data)//n
blocks = [None] * tot_blocks
print(f'Total blocks: {tot_blocks}')

for i in range(tot_blocks):
    blocks[i] = data[i*n:(i+1)*n]

with open(sys.argv[2], 'w') as f:
    for block in blocks:
        for snp in range(len(block[0])):
            for ind in range(n):
                f.write(block[ind][snp])
                if ind < n - 1:
                    f.write(' ')
            f.write('\n')
