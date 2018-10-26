# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 15:29:37 2018

@author: Frederic
"""
import options
from src.clustering import perform_clustering, find_tree, find_distances, validate_dists

VERBOSE = options.VERBOSE

def run_once(boot: int, K: int, T: float, simulation) -> None:
    outgroups = simulation.outgroups
    suffix = f'.boot.{boot}' if boot > -1 else ''
    labels, short_array, lambdas, res_labels = perform_clustering(K, simulation)

    tree, ns, blocks = find_tree(K, simulation.asd_pattern.format(1) + suffix, labels, short_array,
                                 outgroups, res_labels)

    res = []
    for _ in range(min(10 + 2 ** K, 100)):
        dists, constraints = find_distances(K, T, tree, ns, lambdas, blocks)
        if VERBOSE >= 1:
            print('Found distances:', dists.x)
        if validate_dists(dists.x, constraints):
            if VERBOSE >= 1:
                print('OK')
            res.append(dists.x)
        else:
            if VERBOSE >= 1:
                print('Invalid')
    return(res)
