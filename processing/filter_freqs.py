#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 12:58:54 2017

@author: ivan
"""

import numpy as np
import pickle
from typing import Tuple, Dict


def load_freqs(freq_pattern: str) -> Tuple[Dict[int, float], Dict[int, float]]:
    """Calculate averages and standard deviations from the per site
    ancestry frequencies from files giben by 'freq_pattern'.
    """
    mu = {}
    sigma = {}
    for group in range(1, 5):
        print(f'Initial loading, group {group}')
        with open(freq_pattern.format(group), 'rb') as f:
            freqs = pickle.load(f)
        mu[group] = np.mean(freqs)
        sigma[group] = np.std(freqs)
    return mu, sigma


def soft_filter(ancestry_file: str, snp_file: str, outfile: str,
                ancestry_outfile: str,
                mu: Dict[int, float], sigma: Dict[int, float],
                width: float) -> None:
    """Apply soft filter by per site ancestry data.

    This removes the sites where ancestries are anomalous, i.e. are outside
    of the interval [mu - width * sigma, mu - width * sigma].
    """
    print(f'Started processing', ancestry_file)
    with open(ancestry_file, 'rb') as f:
        ancestry = pickle.load(f)
    print(f'Loaded ancestry from', ancestry_file)
    with open(snp_file, 'rb') as f:
        snps = pickle.load(f)
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
    with open(outfile, 'wb') as f:
        pickle.dump(snps, f)
    ancestry = ancestry[np.where(cond)]
    with open(ancestry_outfile, 'wb') as f:
        pickle.dump(ancestry, f)
