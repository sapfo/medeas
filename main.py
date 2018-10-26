#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 12:20:25 2017

@author: ivan


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
        asd_main(1, Simulation,  txt_format=True)
        asd_main(2, Simulation,  txt_format=True)
    else:
        asd_main(1, Simulation)
        asd_main(2, Simulation)

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

    Simulation.all_res = []

    for boot in range(-1, BOOTRUNS):
        boot_res = run_once(boot, K, T, Simulation)
        for res in boot_res:
            Simulation.all_res.append(res)

    Simulation.generate_output()

