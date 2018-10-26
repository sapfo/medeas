#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 12:20:25 2017

@author: ivan


"""

import matplotlib
matplotlib.use('Agg')

from src.simulation import SimulationInfo
from src.make_asd import asd_main
from src.mds import calc_mds
from src.lambda_analyze import find_T_and_L, find_K
from single_pass import run_once
from prepare import preprocess_data



simulation = SimulationInfo()


if not simulation.simulation and not simulation.skip_preprocessing:
    preprocess_data(simulation.snps_pattern, simulation.ancestry_pattern, simulation.output_folder, simulation.output_file, simulation.chromosomes, simulation.labels_file)



if not simulation.skip_calculate_matrix:
    if simulation.simulation:
        scrm_file_clean_result = simulation.snps_pattern
        asd_main(1, simulation, txt_format=True)
        asd_main(2, simulation, txt_format=True)
    else:
        asd_main(1, simulation)
        asd_main(2, simulation)

    calc_mds(simulation.asd_pattern.format(1), simulation.vec_pattern.format(1))
    calc_mds(simulation.asd_pattern.format(2), simulation.vec_pattern.format(2))
    for boot in range(simulation.bootstrap_number):
        suffix = f'.boot.{boot}'
        calc_mds(simulation.asd_pattern.format(1) + suffix, simulation.vec_pattern.format(1) + suffix)
        calc_mds(simulation.asd_pattern.format(2) + suffix, simulation.vec_pattern.format(2) + suffix)


if not simulation.skip_analysis:
    T, L = find_T_and_L(simulation, simulation.vec_pattern.format(2))
    K = find_K(simulation.vec_pattern.format(2), L, T, simulation)
    K_over = simulation.K
    if K_over:
        print(f'OVERRIDING  K = {K} WITH: K = {K_over}')
        K = K_over

    simulation.all_res = []

    for boot in range(-1, simulation.bootstrap_number):
        boot_res = run_once(boot, K, T, simulation)
        for res in boot_res:
            simulation.all_res.append(res)

    simulation.generate_output()

