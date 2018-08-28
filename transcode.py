#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 12:11:31 2017

@author: ivan
"""

# According to scrm, and we follow the convention here, a locus is
# a physical region which  might contain various snp.

import sys

nb_individual = int(sys.argv[4])
print(f'n = {nb_individual}')

with open(sys.argv[1]) as f:
    data = f.readlines()

data = [d.strip() for d in data]
seed = data[0]
with open(sys.argv[3], 'w') as f:
    f.write(seed)
data = [d for d in data[1:] if d]

nb_loci = len(data) // nb_individual
loci = [None] * nb_loci
print(f'number of loci: {nb_loci}')

for i in range(nb_loci):
    loci[i] = data[i * nb_individual:(i + 1) * nb_individual]

with open(sys.argv[2], 'w') as f:
    for locus in loci:
        for snp in range(len(locus[0])):
            for individual in range(nb_individual):
                s = locus[individual][snp]
                f.write(str(int(s) + 1))
                if individual < nb_individual - 1:
                    f.write(' ')
            f.write('\n')
