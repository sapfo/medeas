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

import argparse
parser = argparse.ArgumentParser()


parser.add_argument("-sf","--snps_file", help="Prepare the data wihtout performing the actual analysis on them",
                    required=True)
parser.add_argument("-lf","--labels_file", help="File containing the labels",
                    required=True)
parser.add_argument("-of","--output_folder", help="Folder where results and temporal data should be store")
parser.add_argument("-K", help="Number of population. If K=0 (default), the soft detect automatically the number of population",
                    type = int, default=0)
parser.add_argument("--n_chromosome", help="Number of chromosome If they are stored in different file",
                    type=int,default=1)
parser.add_argument("--outgroup", help="Who is the outgroup in your data", nargs='+')


parser.add_argument("-af","--ancestry_file", help="File containing the ancestry of each locus")



parser.add_argument("-bws","--boot_window_size",
                    help="How many markers do we have in each bootstraping windows",
                    type = int, default=100)

parser.add_argument("--simulation", help="Does the data come from a simulation",
                    action="store_true")
parser.add_argument("--skip_calculate_matrix", help="Skip the computation of the distance matrices and the related MDS matrix",
                    action="store_true")
parser.add_argument("--skip_preprocessing", help="Directly proceed to analysis without preparing the data",
                    action="store_true")
parser.add_argument("--skip_analysis", help="Prepare the data wihtout performing the actual analysis on them",
                    action="store_true")

args = parser.parse_args()
ancestry_pattern =args.ancestry_file
snps_pattern = args.snps_file
output_folder = args.output_folder
labels_file = args.labels_file
chromosomes = range(1, args.n_chromosome+1)
bootsize = args.boot_window_size

import options

BOOTRUNS = options.BOOTRUNS
VERBOSE = options.VERBOSE




import os
import numpy as np
import matplotlib
matplotlib.use(options.MATHPLOTLIB_BACKEND)


from src.make_asd import asd_main
from src.mds import calc_mds
from src.lambda_analyze import find_T_and_L, find_K
from prepare import preprocess_data
from single_pass import run_once


output_file = os.path.join(output_folder, 'processed', 'chr.{}.stped')

if not args.simulation and not args.skip_preprocessing:
    preprocess_data(snps_pattern,ancestry_pattern, output_folder,output_file,chromosomes,labels_file)

asd_folder = "asd_matrices"
mds_folder = "MDS_eigensystem"
asd_full_path = os.path.join(output_folder,asd_folder)
mds_full_path = os.path.join(output_folder,mds_folder)
all_path = [asd_full_path,mds_full_path]
for path in all_path:
    if not os.path.exists(path):
        os.makedirs(path)
asd_pattern = os.path.join(asd_full_path, 'p{}.asd.data')
vec_pattern = os.path.join(mds_full_path, 'p{}.vecs.data')

if not args.skip_calculate_matrix:
    if args.simulation:
        scrm_file_clean_result = snps_pattern
        asd_main(1, scrm_file_clean_result, asd_pattern.format(1), chromosomes,bootsize, txt_format=True)
        asd_main(2, scrm_file_clean_result, asd_pattern.format(2), chromosomes,bootsize, txt_format=True)
    else:
        asd_main(1, output_file, asd_pattern.format(1),chromosomes,bootsize)
        asd_main(2, output_file, asd_pattern.format(2),chromosomes,bootsize)

    calc_mds(asd_pattern.format(1), vec_pattern.format(1))
    calc_mds(asd_pattern.format(2), vec_pattern.format(2))
    for boot in range(BOOTRUNS):
        suffix = f'.boot.{boot}'
        calc_mds(asd_pattern.format(1) + suffix, vec_pattern.format(1) + suffix)
        calc_mds(asd_pattern.format(2) + suffix, vec_pattern.format(2) + suffix)


if not args.skip_analysis:
    T, L, shift = find_T_and_L(vec_pattern.format(2))
    K = find_K(vec_pattern.format(2), L, T,shift)
    K_over = args.K
    if K_over:
        print(f'OVERRIDING  K = {K} WITH: K = {K_over}')
        K = K_over

    all_res = []

    outgroups = args.outgroup

    for boot in range(-1, BOOTRUNS):
        boot_res = run_once(boot,outgroups,K,T,asd_pattern, vec_pattern,labels_file)
        for res in boot_res:
            all_res.append(res)


    with open(os.path.join(output_folder, "all_extrapolated_distances.txt"), 'w') as f:
        np.savetxt(f, all_res)

