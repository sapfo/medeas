#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 14:56:35 2017

@author: ivan
"""

import numpy as np


def hard_filter(ancestry_file: str, new_ancestry_file: str,
                infile: str, outfile: str,
                group: int, ratio: float) -> None:
    print(f'Started processing', ancestry_file)
    
    ancestry = np.genfromtxt(ancestry_file, dtype='int8')
    snps = np.loadtxt(infile, dtype='int8')
    if ancestry.shape[0]:  # guard for empty files
        is_ancestry = 1 - np.sign(ancestry - group) ** 2
        freqs = np.sum(is_ancestry, axis=1)
        print(f'Loaded ancestry from', ancestry_file)
        max_freq = is_ancestry.shape[1]
        cond = freqs > max_freq * ratio

        rows = np.where(cond)
        snps = snps[rows]
        ancestry = ancestry[rows]
    np.savetxt(outfile, snps, fmt='%1d')
    np.savetxt(new_ancestry_file, ancestry, fmt='%1d')
