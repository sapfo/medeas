#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 12:20:25 2017

@author: ivan

This file is the 'medeas' main entry point.

Here is the summary of how things are organized.
On large scale, there are two folders: 'procesing' and 'src'.
The first one contains files for preprocessing of the data before
actuall analyzis, this includes various filtering steps, format
transforms, generation of auxilliary files, setting missingness etc.
(see function docstrings for details). The second folder contains modules
for actual analyzis using the eigenvalue statistics.

Here is the summary of steps in normal data pre-processing:
1. Convert the input files with SNP data and ancestry data into internal
   binary format (currently just stdlib pickle).
2. Calculate ancestry frequencies per site.
3. Manually remove some populations/individuals (if desired).
3. Filter sites depending on ancetry frequencies, this removes the sites
   where ancestry is anomaluos (i.e. to far from average).
4. Filter sites by keeping only those where amount of given ancestry is
   larger than given level.
5. Set SNP data with given ancestry as missing.
6. Filter individuals by keeping only those with given amount of non-missing
   sites.

Here is the summary of actual analysis workflow:
1. Calculate distance matrixes for p=1 and p=2 (together with several
   bootstrapped -- resampled with replacement -- matrixes).
2. Calculate the eigensystem for two matrixes.
3. Find T (total tree length) and L (effective number of markers) by fitting
   the bulk distribution of eigenvalues for p=2.
4. Find the number of populations (clusters) K using the Tracy-Widom
   statistics.
5. Perform agglomerative clustering using K and p=1 distances.
6. Find the matrix of distances between centers of mass of every cluster.
7. Use neighbour join and an outgroup to find the tree topology from
   the above matrix.
8. Use the tree to construct the constraints and equations for split times D.
9. Find the split times, repeat steps 5-9 to bootstrap
   the confidence interval.
"""

# TODO: use argparse and/or config file

SIMULATION = False

import options
options.TESTING = False  # Shows some debugging info and many plots
options.FST = True  # Also calculates F_ST
options.BOOTRUNS = BOOTRUNS = 10  # How many bootstrap runs we need.

from typing import List, Iterable, Callable

import sys
import numpy as np
import matplotlib.pyplot as plt
from skbio.tree import TreeNode

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


def do(task: Callable[..., None], count: int = cpu_count()):
    """A simple helper to paralelize given task across chromosomes."""
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
                                 ['ABO', 'WCD', 'PAP'], read_labs(labels_file))
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
    if SIMULATION:
        asd_main(1, '/Users/ivan/scrm/tmpres.txt', asd_pattern.format(1), txt_format=True)
        asd_main(2, '/Users/ivan/scrm/tmpres.txt', asd_pattern.format(2), txt_format=True)
    else:
        asd_main(1, missingnes_pattern + '.filtered', asd_pattern.format(1))
        asd_main(2, missingnes_pattern + '.filtered', asd_pattern.format(2))

if '-analyze' in sys.argv:
    calc_mds(asd_pattern.format(1), vec_pattern.format(1))
    for boot in range(BOOTRUNS):
        suffix = f'.boot.{boot}'
        calc_mds(asd_pattern.format(1) + suffix, vec_pattern.format(1) + suffix)
    calc_mds(asd_pattern.format(2), vec_pattern.format(2))
    T, L = find_T_and_L(vec_pattern.format(2))
    K = find_K(vec_pattern.format(2), L, T)
    print('Number of clusters found:', K)
    K_over = 2
    if K_over:
        print(f'OVERRIDING WITH: K = {K_over}')
        K = K_over

    res = []
    if SIMULATION:
        new_labels_file = '/Users/ivan/scrm/fake_labs.txt'
    else:
        new_labels_file =  labels_file + '.filtered'

    def run_once(boot: int) -> None:
        suffix = '' if boot == -1 else f'.boot.{boot}'

        labels, short_array, lambdas, res_labels = perform_clustering(K,
                                                          vec_pattern.format(1) + suffix,
                                                          new_labels_file)
        outgroups = ['PAP']
        tree, ns, blocks = find_tree(K, asd_pattern.format(1) + suffix, labels, short_array,
                                     outgroups, res_labels)

        for _ in range(min(10 + 2**K, 100)):
            dists, constraints = find_distances(K, T, tree, ns, lambdas, blocks)
            print('Found distances:', dists.x)
            if validate_dists(dists.x, constraints):
                print('OK')
                res.append(dists.x)
            else:
                print('Invalid')

    for boot in range(-1, BOOTRUNS):
        run_once(boot)
    if res:
        Dst = np.mean([d[0] for d in res])
        delta_Dst = np.std([d[0] for d in res])
    else:
        Dst = delta_Dst = 0
    with open('res10.txt', 'a') as f:
        f.write(f' {Dst} {delta_Dst} {K} {T} {L}\n')

# TODO: Refactor main into four parts: actual main, prepare.py, single_pass.py, bootstrap.py
