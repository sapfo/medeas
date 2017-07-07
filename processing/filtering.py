#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 16:27:49 2017

@author: ivan
"""

import numpy as np
import sys

num = int(sys.argv[1])

def process_chromosome(num: int) -> None:
    print(f'Started processing chromosome #{num}')
    
    full_data = np.genfromtxt(f'eur.chi.pap.wcd.abo.chr{num}.txt',
                              dtype='int8', delimiter=1)
    print(f'Loaded SNPs on chromosome #{num}')
    
    ancestry = np.genfromtxt(f'eur.chi.pap.wcd.abo.chr{num}.g10.txt.0.Viterbi.txt',
                             dtype='int8')
    print(f'Loaded ancestry for chromosome #{num}')
    
    with open('haplotype_labels.txt') as f:
        labels = f.readlines()
    labels = np.array([l.split()[0] for l in labels])
    abo_columns = np.where(labels == 'ABO')[0]
    wcd_columns = np.where(labels == 'WCD')[0]
    abo_data = full_data[:, abo_columns]
    wcd_data = full_data[:, wcd_columns]
    print('Selected columns from data')
    
    is_wcd_ancestry = 1 - np.sign(ancestry - 4) ** 2 # 4 is fot WCD ancestry
    abo_filtered = (abo_data + 1)[:-1] * is_wcd_ancestry  # -1 ONLY FOR chr=21
    print('Filtered data')
    
    out_data = np.concatenate((abo_filtered, wcd_data[:-1] + 1), axis=1) # -1 ONLY FOR chr=21
    np.savetxt(f'abo.all.chr{num}.stped', out_data, fmt='%1d')
    print(f'Saved data for chromosome #{num}')

for num in [21]:
    process_chromosome(num)