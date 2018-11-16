# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 15:29:37 2018

@author: Frederic
"""
import pickle
from src.clustering import find_distances, validate_dists, build_distance_subblock
from src.lambda_analyze import find_T_and_L

def run_once(boot: int, simulation) -> None:
    suffix = f'.boot.{boot}' if boot > -1 else ''
    with open(simulation.vec_pattern.format(1) + suffix, 'rb') as f:
        lambdas, vec = pickle.load(f)
    with open(simulation.asd_pattern.format(1) + suffix, 'rb') as f:
        delta = pickle.load(f)
    T, L = find_T_and_L(simulation, simulation.vec_pattern.format(2) + suffix)
    distance_subblocks = build_distance_subblock(simulation.K, simulation.used_labels, delta)
    res = []
    for _ in range(min(10 + 2 ** simulation.K, 100)):
        dists, constraints = find_distances(simulation.K, T, simulation.tree, simulation.ns, lambdas, distance_subblocks, simulation)
        if simulation.output_level >= 1:
            print('Found distances:', dists.x)
        if validate_dists(dists.x, constraints):
            if simulation.output_level >= 1:
                print('OK')
            res.append(dists.x)
        else:
            if simulation.output_level >= 1:
                print('Invalid')
    return(res)
