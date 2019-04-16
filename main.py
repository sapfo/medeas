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
from src.lambda_analyze import find_T_and_t_within
from single_pass import run_once
from src.clustering import perform_clustering, find_tree, get_mds_coordinate, set_tree_from_input, build_population_dimension

simulation = SimulationInfo()

if not simulation.skip_calculate_matrix:
    compute_asd_matrix(1, simulation, txt_format=simulation.simulation)
    compute_asd_matrix(2, simulation, txt_format=simulation.simulation)

with open(simulation.asd_pattern.format(1), 'rb') as f:
     delta = pickle.load(f)
simulation.plot_distance_matrix(delta)

calc_mds(simulation.asd_pattern.format(1), simulation.vec_pattern.format(1))

calc_mds(simulation.asd_pattern.format(2), simulation.vec_pattern.format(2))
simulation.plot_eigenvalues()

for boot in range(simulation.bootstrap_number):
    suffix = f'.boot.{boot}'
    calc_mds(simulation.asd_pattern.format(1) + suffix, simulation.vec_pattern.format(1) + suffix)
    calc_mds(simulation.asd_pattern.format(2) + suffix, simulation.vec_pattern.format(2) + suffix)

T,ts = find_T_and_t_within(simulation.asd_pattern.format(2), simulation.labels)


coordinates_mds = get_mds_coordinate(simulation, 1)
simulation.plot_mds(coordinates_mds,"MDS_")
ns = build_population_dimension(simulation.K,simulation.numerical_labels)
simulation.ns = ns

coordinates_pca = get_mds_coordinate(simulation, 2)
simulation.plot_mds(coordinates_pca, "PCA_")
if simulation.topology == None:
    tree = find_tree(simulation.K, simulation.asd_pattern.format(1), simulation.numerical_labels, coordinates_mds, simulation)
else:
    tree = set_tree_from_input(simulation.asd_pattern.format(1), simulation)

simulation.tree = tree
simulation.plot_tree()


simulation.all_distance = []
simulation.all_effective_size = []
for boot in range(-1, simulation.bootstrap_number):
    (distances, effective_size) = run_once(boot, simulation)
    if distances is not None:
        simulation.all_distance.append(distances)
        simulation.all_effective_size.append(effective_size)


simulation.generate_final_output()

