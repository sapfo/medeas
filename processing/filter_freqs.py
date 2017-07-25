#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 12:58:54 2017

@author: ivan
"""

import numpy as np

def load_freqs(freq_pattern):
    mu = {}
    sigma = {}
    for group in range(1, 5):
        print(f'Initial loading, group {group}')
        freqs = np.loadtxt(freq_pattern.format(group), dtype='int8')
        mu[group] = np.mean(freqs)
        sigma[group] = np.std(freqs)
    return mu, sigma

def soft_filter(ancestry_file: str, snp_file: str, outfile: str, ancestry_outfile: str,
                mu, sigma, width: float) -> None:

    print(f'Started processing', ancestry_file)
    
    ancestry = np.genfromtxt(ancestry_file, dtype='int8')
    print(f'Loaded ancestry from', ancestry_file)
    snps = np.loadtxt(snp_file, dtype='int8')
    cond = np.ones((snps.shape[0],))

    group: int
    for group in range(1, 5):
        is_ancestry = 1 - np.sign(ancestry - group) ** 2
        print(f'Processed selectors for group {group}')
        freqs = np.sum(is_ancestry, axis=1)
        print(f'Processed filters for group {group}')
        cond = np.logical_and(cond, freqs < mu[group] + width*sigma[group])
        cond = np.logical_and(cond, freqs > mu[group] - width*sigma[group])
        
    snps = snps[np.where(cond)]
    np.savetxt(outfile, snps, fmt='%1d')
    ancestry = ancestry[np.where(cond)]
    np.savetxt(ancestry_outfile, ancestry, fmt='%1d')
