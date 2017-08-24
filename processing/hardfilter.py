#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 14:56:35 2017

@author: ivan
"""

import numpy as np
from typing import List
import pickle


def hard_filter(ancestry_file: str, new_ancestry_file: str,
                infile: str, outfile: str,
                groups: List[int], ratio: float) -> None:
    print(f'Started processing', ancestry_file)
    
    with open(ancestry_file, 'rb') as f:
        ancestry = pickle.load(f)
    with open(infile, 'rb') as f:
        snps = pickle.load(f)
    if ancestry.shape[0]:  # guard for empty files
        is_ancestry = np.zeros(ancestry.shape)
        for group in groups:
            is_ancestry = np.logical_or(is_ancestry,
                                      1 - np.sign(ancestry - group) ** 2)
        freqs = np.sum(is_ancestry, axis=1)
        print(f'Loaded ancestry from', ancestry_file)
        max_freq = is_ancestry.shape[1]
        cond = freqs > max_freq * ratio

        rows = np.where(cond)
        snps = snps[rows]
        ancestry = ancestry[rows]
    with open(outfile, 'wb') as f:
        pickle.dump(snps, f)
    with open(new_ancestry_file, 'wb') as f:
        pickle.dump(ancestry, f)
