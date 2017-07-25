#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 16:27:49 2017

@author: ivan
"""

import numpy as np
from typing import List

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

def set_missing(ancestry_file: str, infile: str, outfile: str,
                all_labels: List[str], groups: List[int]) -> None:
    print("Processing file", ancestry_file)
    all_labels = np.array(all_labels)
    columns = np.where(all_labels == 'ABO')
    full_data = np.genfromtxt(infile, dtype='int8')
    ancestry = np.genfromtxt(ancestry_file, dtype='int8')
    if ancestry.shape[0]: # guard for empty files
        for group in groups:
            is_not_ancestry = np.sign(ancestry - group) ** 2
            full_data.T[columns] = full_data.T[columns] * is_not_ancestry.T
    np.savetxt(outfile, full_data, fmt='%1d')

def filter_sparse(infile: str, outfile: str, ratio: float,
                  labels: List[str]) -> 'np.ndarray[str]':
    labels = np.array(labels)
    data = np.genfromtxt(infile, dtype='int8')
    total = data.shape[0]
    non_missing = np.sum(np.sign(data), axis=0)
    columns = np.where(non_missing > ratio * total)
    data = data.T[columns].T
    labels = labels[columns]
    np.savetxt(outfile, data, fmt='%1d')
    return labels

def filter_manual(infile: str, outfile: str, pops: List[str],
                  labels: List[str]) -> 'np.ndarray[str]':
    short_labs = np.array([l.split()[0] for l in labels])
    print(short_labs)
    labels = np.array(labels)
    data = np.genfromtxt(infile, dtype='int8')

    mask = np.ones((data.shape[1],))
    for pop in pops:
        mask[np.where(short_labs == pop)] = 0
    columns = np.where(mask)

    data = data.T[columns].T
    labels = labels[columns]
    print(labels)
    np.savetxt(outfile, data, fmt='%1d')
    return labels