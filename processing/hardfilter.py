#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 14:56:35 2017

@author: ivan
"""

import numpy as np

mu = {}
sigma = {}
for group in range(1, 5):
    print(f'Initial loading, group {group}')
    freqs = np.loadtxt(f'abo.all.freqs.{group}.txt', dtype='int8')
    print(freqs.shape)
    mu[group] = np.mean(freqs)
    sigma[group] = np.std(freqs)

def process_chromosome(num: int) -> None:
    print(f'Started processing chromosome #{num}')
    
    ancestry = np.genfromtxt(f'eur.chi.pap.wcd.abo.chr{num}.g10.txt.0.Viterbi.txt',
                             dtype='int8')
    print(f'Loaded ancestry for chromosome #{num}')
    snps = np.loadtxt(f'abo.all.chr{num}.stped', dtype='int8')
    cond = np.ones((snps.shape[0],))

    group: int
    for group in range(1, 5):
        is_ancestry = 1 - np.sign(ancestry - group) ** 2
        print(f'Processed selectors for group {group}')
        freqs = np.sum(is_ancestry, axis=1)
        print(f'Processed filters for group {group}')
        cond = np.logical_and(cond, freqs < mu[group] + 2*sigma[group])
        cond = np.logical_and(cond, freqs > mu[group] - 2*sigma[group])
    cond = np.logical_and(cond, freqs > 58)
        
    snps = snps[np.where(cond)]
    np.savetxt(f'abo.all.chr{num}.hardfilter.stped', snps, fmt='%1d')

for num in range(1, 23):
    process_chromosome(num)