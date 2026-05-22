#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 12:20:25 2017

@author: ivan


"""

import matplotlib
matplotlib.use('Agg')
import numpy as np
import os

from src.simulation import SimulationInfo
import pickle
from src.make_asd import compute_asd_matrix
from src.make_asd import compute_asd_matrix_plink
from src.mds import calc_mds

from single_pass import run_once
from src.clustering import find_tree, get_mds_coordinate, set_tree_from_input, build_population_dimension
import sys
from src.extrapolate_split_time import extrapolate_split_time

def main():
    simulation = SimulationInfo()

    required_asd = [simulation.asd_pattern.format(1), simulation.asd_pattern.format(2)]

    if simulation.skip_calculate_matrix:
        missing = [p for p in required_asd if not os.path.isfile(p)]
        if missing:
            msg = [
                "Error: --skip-calculate-matrix was requested, but required ASD matrices are missing:",
            ]
            msg.extend([f"  - {p}" for p in missing])
            msg.append(
                "Run once without --skip_calculate_matrix to generate these files, or point -of to an existing output folder."
            )
            sys.exit("\n".join(msg))

    else:
        try:
            if simulation.use_plink:
                compute_asd_matrix_plink(simulation)
            else:
                compute_asd_matrix(simulation)
        except Exception as e:
            import traceback; traceback.print_exc()
            sys.exit("Error: A problem occurs when computing the distance matrix. Please check that your genotype matrix is in the right format.")

        if not simulation.use_plink:
            if simulation.detailed_output:
                simulation.export_sfs()
            else:
                simulation.export_sfs_simple()

    # Loading delta the distance matrix for p = 1
    if simulation.use_plink:
        delta = np.loadtxt(simulation.asd_pattern.format(1) + ".mdist.gz")
    else:
        with open(simulation.asd_pattern.format(1), 'rb') as f:
            delta = pickle.load(f)
         

    # plot distance histogram
    if simulation.detailed_output:
        simulation.plot_distance_matrix(delta)

    if simulation.output_level > 1:
        print(f"number of individual in the distance matrix: {len(delta)}")

    # Compute mds for p = 1 and for p = 2 with all data
    if simulation.use_plink:
        calc_mds(simulation.asd_pattern.format(1) + ".mdist.gz", simulation.vec_pattern.format(1), "MDS")
        calc_mds(simulation.asd_pattern.format(2) + ".mdist.gz", simulation.vec_pattern.format(2), "PCA")
    else:
        calc_mds(simulation.asd_pattern.format(1), simulation.vec_pattern.format(1), "MDS")
        calc_mds(simulation.asd_pattern.format(2), simulation.vec_pattern.format(2), "PCA")

    if simulation.detailed_output:
        simulation.plot_eigenvalues()
    else:
        simulation.plot_eigenvalues_simple()

    if not simulation.no_mds:
        coordinates_mds = get_mds_coordinate(simulation, 1)
        simulation.plot_mds(coordinates_mds, "MDS_", "mds_plot")
        with open(os.path.join(simulation.output_folder, "MDS_coordinate.txt"), 'w') as f:
            np.savetxt(f, coordinates_mds)

    if not simulation.no_pca:
        coordinates_pca = get_mds_coordinate(simulation, 2)
        simulation.plot_mds(coordinates_pca, "PCA_", "pca_plot")
        with open(os.path.join(simulation.output_folder, "PCA_coordinate.txt"), 'w') as f:
            np.savetxt(f, coordinates_pca)

    ## stop here if no split estimation is requested
    if simulation.no_split:
        simulation.generate_mds_pca_only_output()
        return


    ## split time estimation
    # Compute mds for p = 1 and for p = 2 for all bootstrap replicate
    for boot in range(simulation.bootstrap_number):
        suffix = f'.boot.{boot}'
        calc_mds(simulation.asd_pattern.format(1) + suffix, simulation.vec_pattern.format(1) +  suffix, f"MDS bootstrap {boot}")
        calc_mds(simulation.asd_pattern.format(2) + suffix, simulation.vec_pattern.format(2) + suffix, f"PCA bootstrap {boot}")

    # ns is the vector of population sample sizes
    ns = build_population_dimension(simulation.K, simulation.numerical_labels)
    simulation.ns = ns

    if simulation.K > 1:
        if simulation.topology is None:
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

    for within_pop_coalescence_time, between_pop_coalescence_time in zip(
        simulation.all_within_pop_coalescence_time,
        simulation.all_between_pop_coalescence_time,
    ):
        effective_sizes, split_times = extrapolate_split_time(
            simulation.tree,
            simulation.split_index_matrix,
            within_pop_coalescence_time,
            between_pop_coalescence_time,
        )
        simulation.all_effective_size.append(effective_sizes)
        simulation.all_split_time.append(split_times)

    simulation.generate_final_output()


if __name__ == "__main__":
    main()
