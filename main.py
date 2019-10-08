#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 12:20:25 2017

@author: ivan


"""

import matplotlib
matplotlib.use('Agg')
import numpy as np

from src.simulation import SimulationInfo
import pickle
from src.make_asd import compute_asd_matrix
from src.mds import calc_mds

from single_pass import run_once
from src.clustering import find_tree, get_mds_coordinate, set_tree_from_input, build_population_dimension
import sys
from src.extrapolate_split_time import extrapolate_split_time

simulation = SimulationInfo()
if not simulation.skip_calculate_matrix:
    try:
        compute_asd_matrix(simulation)
    except:
        sys.exit("Error: A problem occurs when computing the distance matrix. Please check that your genotype matrix is in the right format.")
    simulation.export_sfs()

with open(simulation.asd_pattern.format(1), 'rb') as f:
     delta = pickle.load(f)


simulation.plot_distance_matrix(delta)


if simulation.output_level > 1:
    print(f"number of individual in the distance matrix: {len(delta)}")
calc_mds(simulation.asd_pattern.format(1), simulation.vec_pattern.format(1))

calc_mds(simulation.asd_pattern.format(2), simulation.vec_pattern.format(2))
simulation.plot_eigenvalues()

for boot in range(simulation.bootstrap_number):
    suffix = f'.boot.{boot}'
    calc_mds(simulation.asd_pattern.format(1) + suffix, simulation.vec_pattern.format(1) + suffix)
    calc_mds(simulation.asd_pattern.format(2) + suffix, simulation.vec_pattern.format(2) + suffix)

coordinates_mds = get_mds_coordinate(simulation, 1)
simulation.plot_mds(coordinates_mds,"MDS_")
ns = build_population_dimension(simulation.K,simulation.numerical_labels)
simulation.ns = ns

coordinates_pca = get_mds_coordinate(simulation, 2)
simulation.plot_mds(coordinates_pca, "PCA_")
if simulation.K > 1:
    if simulation.topology == None:
        tree = find_tree(simulation.K, simulation.numerical_labels, coordinates_mds)
    else:
        tree = set_tree_from_input(simulation.asd_pattern.format(1), simulation)
else:
    exit("Error: Unable to perform tree reconstruction. You should have more than one population")
simulation.set_tree(tree)
simulation.save_tree()


simulation.all_between_pop_coalescence_time = []
simulation.all_within_pop_coalescence_time = []
simulation.all_T = []
if simulation.pops_contain_at_least_2_individual():
    for boot in range(-1, simulation.bootstrap_number):
        (distances, effective_size, T) = run_once(boot, simulation)
        if distances is not None:
            simulation.all_between_pop_coalescence_time.append(distances)
            simulation.all_within_pop_coalescence_time.append(effective_size)
            simulation.all_T.append(T)
else:
    exit("Error: Unable to perform the time inference.Each population should have more than one individual")

simulation.all_split_time = []
simulation.all_effective_size = []

for within_pop_coalescence_time, between_pop_coalescence_time in zip(simulation.all_within_pop_coalescence_time,
                                                                     simulation.all_between_pop_coalescence_time):
    effective_sizes, split_times = extrapolate_split_time(simulation.tree,
                                                          simulation.split_index_matrix,
                                                          within_pop_coalescence_time,
                                                          between_pop_coalescence_time)
    simulation.all_effective_size.append(effective_sizes)
    simulation.all_split_time.append(split_times)


simulation.generate_final_output()

