#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 12:20:25 2017

@author: ivan


"""

import matplotlib
matplotlib.use('Agg')

from src.simulation import SimulationInfo
from src.make_asd import compute_asd_matrix
from src.mds import calc_mds
from src.lambda_analyze import find_T_and_L, find_K
from single_pass import run_once
from src.clustering import perform_clustering, find_tree, get_coordinate
simulation = SimulationInfo()

if not simulation.skip_calculate_matrix:
    compute_asd_matrix(1, simulation, txt_format=simulation.simulation)
    compute_asd_matrix(2, simulation, txt_format=simulation.simulation)

    calc_mds(simulation.asd_pattern.format(1), simulation.vec_pattern.format(1))
    calc_mds(simulation.asd_pattern.format(2), simulation.vec_pattern.format(2))

    for boot in range(simulation.bootstrap_number):
        suffix = f'.boot.{boot}'
        calc_mds(simulation.asd_pattern.format(1) + suffix, simulation.vec_pattern.format(1) + suffix)
        calc_mds(simulation.asd_pattern.format(2) + suffix, simulation.vec_pattern.format(2) + suffix)



T, L = find_T_and_L(simulation, simulation.vec_pattern.format(2))
K = find_K(simulation.vec_pattern.format(2), L, T, simulation)

if simulation.K:
    print(f'OVERRIDING  K = {K} WITH: K = {simulation.K}')
    K = simulation.K
else:
    simulation.K = K

coordinates_mds = get_coordinate(simulation, 1)
inferred_labels = perform_clustering(simulation.K, coordinates_mds, simulation)


simulation.plot_mds(coordinates_mds,inferred_labels,"MDS_")
coordinates_pca = get_coordinate(simulation, 2)
simulation.plot_mds(coordinates_pca,inferred_labels,"PCA_")

tree, ns, blocks = find_tree(K, simulation.asd_pattern.format(1), inferred_labels, coordinates_mds, simulation)
simulation.tree = tree
simulation.ns = ns
simulation.blocks = blocks

simulation.plot_tree()


simulation.all_res = []

for boot in range(-1, simulation.bootstrap_number):
    boot_res = run_once(boot, simulation.K, simulation)
    for res in boot_res:
        simulation.all_res.append(res)

simulation.generate_final_output()

