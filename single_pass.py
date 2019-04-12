# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 15:29:37 2018

@author: Frederic
"""
import pickle
from src.clustering import find_distances, validate_dists, build_distance_subblock
from src.lambda_analyze import find_T_and_t_within

def run_once(boot: int, simulation) -> None:
    suffix = f'.boot.{boot}' if boot > -1 else ''
    with open(simulation.vec_pattern.format(1) + suffix, 'rb') as f:
        lambdas, vec = pickle.load(f)
    with open(simulation.asd_pattern.format(1) + suffix, 'rb') as f:
        delta = pickle.load(f)
    T, t_within = find_T_and_t_within(simulation.asd_pattern.format(1) + suffix, simulation.labels)
    distance_subblocks = build_distance_subblock(simulation.K, simulation.numerical_labels, delta)
    distance_validity = False
    nbLoop = 0
    maxNbLoop = min(10 + 2 ** simulation.K, 100)
    while not distance_validity:
        if nbLoop == maxNbLoop: break
        nbLoop = nbLoop + 1
        dists, constraints = find_distances(simulation.K, T,t_within, simulation.tree, simulation.ns, lambdas,
                                            distance_subblocks, simulation.output_level)
        distance_validity = validate_dists(dists.x, constraints)
        if simulation.output_level >= 1:
            print('Found distances:', dists.x)
            if distance_validity:
                print('Valide distance')
            else:
                print('Invalid distance')

    else:
        return(dists.x)

    print('Unable to find valid distances for this bootstrap sample')
    return(None)
