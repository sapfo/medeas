# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 15:29:37 2016

@author: ivan
"""

import sys
import numpy as np
from multiprocessing import Process, Queue, cpu_count
from typing import Tuple, Callable, List, IO
import pickle
from random import randint
import matplotlib.pyplot as plt


NPROC = cpu_count()


def dist_and_norm(a: 'np.ndarray[int]', b: 'np.ndarray[int]',
                  dist_func: Callable[[np.ndarray], np.ndarray]
                  ) -> Tuple[float, int]:
    """Calculate un-normalized distance and norm between vectors 'a' and 'b'.
    Norm is number of sites 'i' where both 'a[i]' and 'b[i]' are non-zero.
    Other sites does not contribute to distance.
    """
    filt = np.logical_and(a, b)
    dst = filt * dist_func(a - b)
    return np.sum(dst), np.sum(filt)


def compute(i: int, data: 'np.ndarray[int]',
            dist_func: Callable[[np.ndarray], np.ndarray],
            N: int) -> Tuple[List[float], List[int]]:
    """Compute all distances and norms for 'i'th row in 'data'."""
    dists: List[float] = []
    norms: List[int] = []
    #if TESTING:
       # print(f'Processing row #{i}')
    for j in range(i + 1, N):
        dist, norm = dist_and_norm(data[i], data[j], dist_func)
        dists.append(dist)
        norms.append(norm)
    return dists, norms


def work(tasks: 'Queue[int]',
         results: 'Queue[Tuple[int, Tuple[List[float], List[int]]]]',
         data: 'np.ndarray[int]', dist_func: Callable[[np.ndarray], np.ndarray],
         N: int) -> None:
    """Compute distaces and norms for rows from 'tasks'."""
    while True:
        i = tasks.get()
        if i < 0:
            return
        results.put((i, compute(i, data, dist_func, N)))


def process(data: 'np.ndarray[int]',
            dist_func: Callable[[np.ndarray], np.ndarray],
            tasks, results, N) -> Tuple['np.ndarray[float]', 'np.ndarray[int]']:
    """Calculate matrices of un-normalized distances and norms for 'data'
    using given distance function.
    """
    dists = np.zeros((N, N))
    norms = np.zeros((N, N))

    for i in range(N):
        tasks.put(i)
    for _ in range(NPROC):
        tasks.put(-1)

    procs = [Process(target=work, args=(tasks, results, data, dist_func, N))
             for _ in range(NPROC)]
    for proc in procs:
        proc.start()

    rest = N
    while rest:
        i, (dist, norm) = results.get()
        rest -= 1
        dists[i, i + 1:] = dist
        norms[i, i + 1:] = norm

    for proc in procs:
        proc.join()

    return dists, norms


def  compute_asd_matrix(simulation) -> None:
    """Calculate the distance matrix for two Minkowsky parapeter pp.
    Take care that here we are assuming haploid individuals!
    'name' is the input file name (or format string) with SNP data.
    'out_name' is the output binary file (pickle).
    """

    # On Windows, processes execute the whole file before forking
    # therefore we protect this code with if __name__ == '__main__'
    # Need to think how to avoid copying ``data`` on forking.
    # Maybe process input file in chunks?
    # On POSIX everything is already fine because of "copy-on-write"


    name =  simulation.snps_pattern
    bootsize = simulation.bootsize
    label = simulation.labels_file
    dist_func = lambda x: np.abs(x)

    # ---------- constants

    # this should be large to avoid overhead of spawning new processes
    # or we need to reuse them somehow
    MAXSIZE = 5000 * 2 ** 20  # 5 Gb


    # ---------- global data

    with open(name) as f:
        N = len(f.readline()) // 2
        print(f'nb individual = {N}')

    tot_dists = np.zeros((N, N))
    tot_norms = np.zeros((N, N))
    delta = np.zeros((N, N))
    tasks = Queue()
    results = Queue()

    f: IO
    chunk_data: List[Tuple['np.ndarray[float]', 'np.ndarray[float]']] = []
    remainder = None  # np.zeros((1, N))

    def process_chunks(data) -> None:
        nonlocal remainder, tot_dists, tot_norms
        start_i = 0
        end_i = bootsize - len(remainder) if remainder is not None else bootsize
        while end_i <= len(data):
            if simulation.output_level >= 1:
                print(f'Processing site {start_i}')
            chunk = data[start_i:end_i]
            if remainder is not None and start_i == 0:
                chunk = np.vstack((remainder, chunk))
            datac = chunk.T.copy()
            dists, norms = process(datac, dist_func, tasks, results, N)
            tot_dists += dists
            tot_norms += norms
            chunk_data.append((dists, norms))
            start_i += bootsize
            end_i = start_i + bootsize
        remainder = data[start_i:]

    pre_sfs = np.empty((0))
    with open(name) as f:
        while True:
            data_lines = f.readlines(MAXSIZE)
            print('Chunk loading started')
            if not data_lines:
                break
            data = np.array([np.fromstring(line, sep=' ',  # line[-cut:-1]
                                           dtype='int8')
                             for line in data_lines])
            # data = data[:, ::2] + data[:, 1::2]
            print('Chunk loaded')

            nb_mut = np.sum(data==1,axis=1)
            nb_missing = np.sum(data==0,axis=1)
            nb_mut[np.where((N - nb_missing)==0)] = 0 ## removing in a stupid way site with no data at all
            nb_missing[np.where((N - nb_missing)==0)] = 1 ## removing in a stupid way site with no data at all
            freq = nb_mut/(N - nb_missing)
            nb_other_mut = np.random.binomial(nb_missing,freq,len(nb_mut))
            nb_mut = nb_mut + nb_other_mut
            nb_mut[nb_mut > N/2] = N - nb_mut[nb_mut > N/2] # Folding the SFS, since 0 and 1 are likely to be not well defined
            pre_sfs = np.append(pre_sfs, nb_mut)
            process_chunks(data)
    sfs = np.unique(pre_sfs, return_counts = True)
    simulation.sfs = np.array(sfs)



    for i in range(N):
        delta[i, i] = 0

    for i in range(N):
        for j in range(i + 1, N):
            delta[i, j] = delta[j, i] = tot_dists[i, j] / tot_norms[i, j]
    pp = 1
    delta_1 = np.copy(delta)
    delta_1 = delta_1 ** (1 / pp)
    out_name =  simulation.asd_pattern.format(pp)
    with open(out_name, 'wb') as f:
        pickle.dump(delta_1, f)
    pp = 2
    delta_2 = np.copy(delta)
    delta_2 = delta_2 ** (1 / pp)
    out_name =  simulation.asd_pattern.format(pp)
    with open(out_name, 'wb') as f:
        pickle.dump(delta_2, f)


    # bootstraping ---------------

    clen = len(chunk_data)
    print(f'NUMBER OF BLOCKS: {clen}')
    for boot in range(simulation.bootstrap_number):
        chunk_res = chunk_data.copy()
        for i in range(clen):
            chunk_res[i] = chunk_data[randint(0, clen - 1)]
        delta = np.zeros((N, N))
        tot_dists = np.sum(np.array([c[0] for c in chunk_res]), axis=0)
        tot_norms = np.sum(np.array([c[1] for c in chunk_res]), axis=0)
        for i in range(N):
            for j in range(i + 1, N):
                delta[i, j] = delta[j, i] = tot_dists[i, j] / tot_norms[i, j]

        pp = 1
        delta_1 = np.copy(delta)
        delta_1 = delta_1 ** (1 / pp)
        out_name = simulation.asd_pattern.format(pp)
        with open(out_name + f'.boot.{boot}', 'wb') as f:
            pickle.dump(delta_1, f)

        pp = 2
        delta_2 = np.copy(delta)
        delta_2 = delta_2 ** (1 / pp)
        out_name = simulation.asd_pattern.format(pp)
        with open(out_name + f'.boot.{boot}', 'wb') as f:
            pickle.dump(delta_2, f)

    print('Distance matrix computed')
    simulation.plot_distance_matrix(delta)




if __name__ == '__main__':
    pp = int(sys.argv[2])
    name = sys.argv[1]  # C:\Users\levkivskyi\PycharmProjects\medeas\test\...
    out_name = f'out.pp.{pp}.asd'

    compute_asd_matrix(pp, name, out_name)
