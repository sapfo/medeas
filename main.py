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



the pattern are the file name and should be passed to it with {}
that will be replaced with the number of the chromosome

"""



import options

BOOTRUNS = options.BOOTRUNS
VERBOSE = options.VERBOSE

import os
import numpy as np
import matplotlib
import sys
matplotlib.use(options.MATHPLOTLIB_BACKEND)
from src.simulation import simulation_info
from src.make_asd import asd_main
from src.mds import calc_mds
from src.lambda_analyze import find_T_and_L, find_K
from single_pass import run_once
from prepare import preprocess_data

Simulation = simulation_info()


if not Simulation.simulation and not Simulation.skip_preprocessing:
    preprocess_data(Simulation.snps_pattern, Simulation.ancestry_pattern, Simulation.output_folder, Simulation.output_file, Simulation.chromosomes, Simulation.labels_file)



if not Simulation.skip_calculate_matrix:
    if Simulation.simulation:
        scrm_file_clean_result = Simulation.snps_pattern
        asd_main(1, Simulation.snps_pattern, Simulation.asd_pattern.format(1), Simulation.chromosomes, Simulation.bootsize, Simulation.labels_file, txt_format=True)
        asd_main(2, Simulation.snps_pattern, Simulation.asd_pattern.format(2), Simulation.chromosomes, Simulation.bootsize, Simulation.labels_file, txt_format=True)
    else:
        asd_main(1, Simulation.output_file, Simulation.asd_pattern.format(1), Simulation.chromosomes, Simulation.bootsize, Simulation.labels_file)
        asd_main(2, Simulation.output_file, Simulation.asd_pattern.format(2), Simulation.chromosomes, Simulation.bootsize, Simulation.labels_file)

    calc_mds(Simulation.asd_pattern.format(1), Simulation.vec_pattern.format(1))
    calc_mds(Simulation.asd_pattern.format(2), Simulation.vec_pattern.format(2))
    for boot in range(BOOTRUNS):
        suffix = f'.boot.{boot}'
        calc_mds(Simulation.asd_pattern.format(1) + suffix, Simulation.vec_pattern.format(1) + suffix)
        calc_mds(Simulation.asd_pattern.format(2) + suffix, Simulation.vec_pattern.format(2) + suffix)


if not Simulation.skip_analysis:
    T, L = find_T_and_L(Simulation.vec_pattern.format(2))
    K = find_K(Simulation.vec_pattern.format(2), L, T)
    K_over = Simulation.K
    if K_over:
        print(f'OVERRIDING  K = {K} WITH: K = {K_over}')
        K = K_over

    all_res = []

    for boot in range(-1, BOOTRUNS):
        boot_res = run_once(boot, Simulation.outgroups, K, T, Simulation.asd_pattern, Simulation.vec_pattern, Simulation.labels_file)
        for res in boot_res:
            all_res.append(res)


    with open(os.path.join(Simulation.output_folder, "all_extrapolated_distances.txt"), 'w') as f:
        np.savetxt(f, all_res)

