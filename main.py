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


Most of the step above are performed if an option is present in the input file
Pre-processing step:
1: -recode
2: -freqs
3: -manual
4: -filter
5 and 6 are related to -hard -miss and -sparse... not so clear what exactly

Actual analysis
1: -asd
2  -> 9: -analyse

The tag SIMULATION should be changed manually if the data comes from simulation.
In this case, the genotypes result are directly read from the file _tmp_res.txt without
the pre-processing step .

Argument used:
    #argv[1] : ancestry pattern (painting)
    #argv[2] : SNP data
    #argv[3] : labels file
    #argv[4:-4]: flags
    #argv[-4]: true nssites
    #argv[-3]: true split time
    #argv[-2]: folder where we want to save the data
    #argv[-1]: filename where to store the summary of the simulation

the pattern are the file name and should be passed to it with {}
that will be replaced with the number of the chromosome

"""

# TODO: use argparse and/or config file



import options

BOOTRUNS = options.BOOTRUNS
SIMULATION = options.SIMULATION

from typing import List, Iterable, Callable

import sys
import os
import numpy as np
import matplotlib
matplotlib.use(options.MATHPLOTLIB_BACKEND)


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

folder = sys.argv[-2]
snps_pattern_stped = snps_pattern + '.stped'
group_pattern = os.path.join(folder, 'group.{}')
freqs_pattern = group_pattern + '.freqs'
filtered_pattern = snps_pattern_stped + '.filtered'
filtered_ancestry_pattern = ancestry_pattern + '.filtered'
hard_filtered_pattern = filtered_pattern + '.hard'
hard_filtered_ancestry_pattern = filtered_ancestry_pattern + '.hard'

missingnes_pattern = os.path.join(folder,'processed','chr.{}.stped')
directory = os.path.join(folder, 'processed')
if not os.path.exists(directory):
    os.makedirs(directory)
asd_pattern = os.path.join(folder, 'pp.{}.asd.data')
vec_pattern = os.path.join(folder, 'pp.{}.vecs.data')


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
    do(recode, 6)

if '-freqs' in sys.argv:
    calculate_freqs(ancestry_pattern + '.recoded', group_pattern)

if '-manual' in sys.argv:
    for n in chromosomes:
        new_labs = filter_manual(snps_pattern_stped.format(n),
                                 snps_pattern_stped.format(n) + '.selected',
                                 ['BRI', 'CHI'], read_labs(labels_file))
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

if not SIMULATION:
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
        scrm_file_clean_result = os.path.join(folder, '_tmp_res.txt')
        asd_main(1, scrm_file_clean_result, asd_pattern.format(1), txt_format=True)
        asd_main(2, scrm_file_clean_result, asd_pattern.format(2), txt_format=True)
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
    print(f'Extrapolated value for the total tree length T: {T}')
    print(f'Extrapolated value for number of loci L:, {L}')
    T = 6.82
    K = find_K(vec_pattern.format(2), L, T)
    K_inf = K
    print('Number of clusters found:', K)
    K_over = options.K_OVERRIDE
    if K_over:
        print(f'OVERRIDING WITH: K = {K_over}')
        K = K_over

    res = []
    if SIMULATION:
        new_labels_file = os.path.join(folder, 'fake_labs.txt')
    else:
        new_labels_file =  labels_file + '.filtered'

    def run_once(boot: int) -> None:
        sL = sys.argv[-4]
        sD = sys.argv[-3]
        suffix = '' if boot == -1 else f'.boot.{boot}'

        if boot == -1 and SIMULATION:
            labels, short_array, lambdas, res_labels = perform_clustering(K,
                                                              vec_pattern.format(2) + suffix,
                                                              new_labels_file)
            np.savetxt(os.path.join(folder, f'XY_p2_L_{sL}_D_{sD}.txt'),
                       short_array[:, :2])
            np.savetxt(os.path.join(folder, f'lambdas_L_{sL}_D_{sD}.txt'), sorted(lambdas, reverse=True))
        labels, short_array, lambdas, res_labels = perform_clustering(K,
                                                          vec_pattern.format(1) + suffix,
                                                          new_labels_file)
        if boot == -1:
            np.savetxt(os.path.join(folder, f'labels_L_{sL}_D_{sD}.txt'), labels)
            np.savetxt(os.path.join(folder, f'XY_p1_L_{sL}_D_{sD}.txt'),
                       short_array[:, :2])

        outgroups = ['PAP']
        tree, ns, blocks = find_tree(K, asd_pattern.format(1) + suffix, labels, short_array,
                                     outgroups, res_labels)


        for _ in range(1):
            dists, constraints = find_distances(K, T, tree, ns, lambdas, blocks)
            print('Found distances:', dists.x)
            if validate_dists(dists.x, constraints):
                print('OK')
                res.append(dists.x)
            else:
                print('Invalid')
        with open('all_distance.txt','a') as f:
            f.write(str(dists.x[0])+'\n')

    with open('all_distance.txt','w') as f:
        f.write('')
    for boot in range(-1, BOOTRUNS):
        run_once(boot)

    if res:
        Dst = np.mean([d for d in res])
        delta_Dst = np.std([d for d in res])
    else:
        Dst = delta_Dst = 0
    with open(os.path.join(folder, sys.argv[-1]), 'a') as f:
    #    f.write('Dst delta_Dst K_inf T L\n')
        f.write(f' {Dst} {delta_Dst} {K_inf} {T} {L}\n')
    with open(os.path.join(folder, "all_extrapolated_distances.txt"), 'a') as f:
        np.savetxt(f,res)

# TODO: Refactor main into four parts: actual main, prepare.py, single_pass.py, bootstrap.py
