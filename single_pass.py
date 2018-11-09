# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 15:29:37 2018

@author: Frederic
"""

from src.clustering import find_distances, validate_dists

def run_once(boot: int, K: int, T: float, simulation) -> None:
    suffix = f'.boot.{boot}' if boot > -1 else ''


    res = []
    for _ in range(min(10 + 2 ** K, 100)):
        dists, constraints = find_distances(K, T, simulation.tree, simulation.ns, simulation.lambdas, simulation.blocks, simulation)
        if simulation.output_level >= 1:
            print('Found distances:', dists.x)
        if validate_dists(dists.x, constraints):
            if simulation.output_level >= 1:
                print('OK')
            res.append(dists.x)
        else:
            if simulation.output_level >= 1:
                print('Invalid')
    simulation.generate_ouput_single_bootstrap(boot)
    return(res)
