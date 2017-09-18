#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 12:20:25 2017

@author: ivan
"""

# TODO: use argparse and/or config file

from typing import List, Iterable

import sys
import numpy as np
import matplotlib.pyplot as plt

from src.make_asd import asd_main
from src.mds import calc_mds
from src.lambda_analyze import find_T_and_L, find_K
from src.clustering import perform_clustering, find_tree, find_distances, validate_dists
from processing.recode import recode_wide
from processing.freqs import calculate_freqs
from processing.filter_freqs import load_freqs, soft_filter
from processing.hardfilter import hard_filter
from processing.filtering import set_missing, filter_sparse, filter_manual
from concurrent.futures import ProcessPoolExecutor as PPE
from multiprocessing import cpu_count

ancestry_pattern = sys.argv[1]
snps_pattern = sys.argv[2] if len(sys.argv) > 2 else None

labels_file = sys.argv[3]
chromosomes = range(1, 23)

snps_pattern_stped = snps_pattern + '.stped'
group_pattern = 'group.{}'
freqs_pattern = group_pattern + '.freqs'
filtered_pattern = snps_pattern_stped + '.filtered'
filtered_ancestry_pattern = ancestry_pattern + '.filtered'
hard_filtered_pattern = filtered_pattern + '.hard'
hard_filtered_ancestry_pattern = filtered_ancestry_pattern + '.hard'

# TODO: Use temporary folder (or read these from argv[]?)
missingnes_pattern = 'processed.chr.{}.stped'
big_file_name = '../raw_data_copy/big_file.stped'

final_big_file = big_file_name + '.final'
asd_pattern = 'pp.{}.asd.data'
vec_pattern = 'pp.{}.vecs.data'

def do(task, count=cpu_count()):
    executor = PPE(count)
    return list(executor.map(task, chromosomes))

def read_labs(file: str) -> List[str]:
    with open(file) as f:
        return f.readlines()

def save_labs(labs: Iterable[str], file: str) -> None:
    with open(file, 'w') as f:
        return f.writelines(labs)

# TODO: add options to read and write text files

if '-recode' in sys.argv:
    def recode(n):
        recode_wide(snps_pattern.format(n), snps_pattern_stped.format(n),
                    ancestry_pattern.format(n),
                    ancestry_pattern.format(n) + '.recoded')
    do(recode, 2)

if '-freqs' in sys.argv:
    calculate_freqs(ancestry_pattern + '.recoded', group_pattern)

if '-manual' in sys.argv:
    for n in chromosomes:
        new_labs = filter_manual(snps_pattern_stped.format(n),
                                 snps_pattern_stped.format(n) + '.selected',
                                 ['ABO', 'WCD'], read_labs(labels_file))
    save_labs(new_labs, labels_file + '.selected')

if '-filter' in sys.argv:
    mu, sigma = load_freqs(freqs_pattern)
    def filterf(n):
        soft_filter(ancestry_pattern.format(n) + '.recoded',
                    snps_pattern_stped.format(n) + '.selected',
                    filtered_pattern.format(n),
                    filtered_ancestry_pattern.format(n), mu, sigma, 5),
    do(filterf)

if '-hard' in sys.argv:
    def hardfilt(n):
        hard_filter(filtered_ancestry_pattern.format(n),
                    hard_filtered_ancestry_pattern.format(n),
                    filtered_pattern.format(n),
                    hard_filtered_pattern.format(n),
                    [4, 3, 2, 1], 0.1) # PAP or WCD
    do(hardfilt)

labs = np.array([l.split()[0] for l in read_labs(labels_file)])

if '-miss' in sys.argv:
    def setmiss(n):
        set_missing(hard_filtered_ancestry_pattern.format(n),
                    hard_filtered_pattern.format(n),
                    missingnes_pattern.format(n),
                    labs, [])  # EUR and CHI
    do(setmiss)

if '-sparse' in sys.argv:
        new_labs = filter_sparse(missingnes_pattern, missingnes_pattern + '.filtered',
                                 0.3, read_labs(labels_file + '.selected'))
        save_labs(new_labs, labels_file + '.filtered')

if '-asd' in sys.argv:
    asd_main(1, missingnes_pattern + '.filtered', asd_pattern.format(1))
    asd_main(2, missingnes_pattern + '.filtered', asd_pattern.format(2))

if '-analyze' in sys.argv:
    calc_mds(asd_pattern.format(1), vec_pattern.format(1))
    calc_mds(asd_pattern.format(2), vec_pattern.format(2))
    T, L = find_T_and_L(vec_pattern.format(2))
    K = find_K(vec_pattern.format(2), L)
    print('Number of clusters found:', K)
    labels, short_array, lambdas, res_labels = perform_clustering(K,
                                                      vec_pattern.format(1),
                                                      labels_file + '.filtered')
    outgroups = ['PAP']
    tree, ns, blocks = find_tree(K, asd_pattern.format(1), labels, short_array,
                                 outgroups, res_labels)
    dists = find_distances(K, T, tree, ns, lambdas, blocks)
    print('Found distances:', dists)

# TODO: Implement statistical bootstrap
# TODO: Refactor main into four parts: actual main, prepare.py, single_pass.py, bootstrap.py
