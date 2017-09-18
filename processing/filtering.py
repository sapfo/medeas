#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 16:27:49 2017

@author: ivan
"""

import numpy as np
from typing import List
import pickle

chromosomes = range(1, 23)


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
    with open(infile, 'rb') as f:
        full_data = pickle.load(f)
    with open(ancestry_file, 'rb')as f:
        ancestry = pickle.load(f)
    if ancestry.shape[0]: # guard for empty files
        for group in groups:
            is_not_ancestry = np.sign(ancestry - group) ** 2
            full_data.T[columns] = full_data.T[columns] * is_not_ancestry.T
    with open(outfile, 'wb') as f:
        pickle.dump(full_data, f)

def filter_sparse(infile_pattern: str, outfile_pattern: str, ratio: float,
                  labels: List[str]) -> 'np.ndarray[str]':
    labels = np.array(labels)
    total = 0
    non_missing = None
    for n in chromosomes:
        print("First read: chromosome", n)
        with open(infile_pattern.format(n), 'rb') as f:
            data = pickle.load(f)
        total += data.shape[0]
        non_missing = (non_missing if non_missing is not None else
                       np.zeros((data.shape[1],))) + np.sum(np.sign(data), axis=0)
    columns = np.where(non_missing > ratio * total)
    for n in chromosomes:
        with open(infile_pattern.format(n), 'rb') as f:
            data = pickle.load(f)
        data = data.T[columns].T
        print("Writing file: chromosome", n)
        with open(outfile_pattern.format(n), 'wb') as f:
            pickle.dump(data, f)
    labels = labels[columns]
    return labels

def filter_manual(infile: str, outfile: str, pops: List[str],
                  labels: List[str]) -> 'np.ndarray[str]':
    short_labs = np.array([l.split()[0] for l in labels])
    print('Manual drop: ', pops, infile)
    labels = np.array(labels)
    with open(infile, 'rb') as f:
        data = pickle.load(f)

    mask = np.ones((data.shape[1],))
    for pop in pops:
        mask[np.where(short_labs == pop)] = 0
    columns = np.where(mask)

    data = data.T[columns].T
    labels = labels[columns]
    with open(outfile, 'wb') as f:
        pickle.dump(data, f)
    return labels