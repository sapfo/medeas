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

from options import TESTING
from options import BOOTRUNS
from options import VERBOSE

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


def asd_main(pp: int, name: str, out_name: str, chromosomes: range,
             bootsize: int,
             txt_format: bool = False
             ) -> None:
    """Calculate the distance matrix with Minkowski parameter 'pp'.

    'name' is the input file name (or format string) with SNP data.
    'out_name' is the output binary file (pickle). If 'txt_format' is
    True, the SNP data will be read from a single text file (useful for
    tests with scrm), otherwise the data will be read from (binary)
    chromosome files.
    """

    # On Windows, processes execute the whole file before forking
    # therefore we protect this code with if __name__ == '__main__'
    # Need to think how to avoid copying ``data`` on forking.
    # Maybe process input file in chunks?
    # On POSIX everything is already fine because of "copy-on-write"
    def test_func(dst: int) -> int:
        if dst > 1:
            return dst ** 2
        else:
            return 0

    dist_func = np.vectorize(test_func)
    dist_func = np.abs
    dist_func = np.square
    dist_func = lambda x: np.abs(x) ** pp

    # ---------- constants

    # this should be large to avoid overhead of spawning new processes
    # or we need to reuse them somehow
    MAXSIZE = 200 * 2 ** 20  # 200 MB


    # ---------- global data

    if txt_format:
        with open(name) as f:
            N = len(f.readline()) // 2
            print(f'N = {N}')
    else:
        with open(name.format(1), 'rb') as f:
            N = pickle.load(f).shape[1]
            print(f'N = {N}')

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
            if VERBOSE >= 1:
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

    if txt_format:
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
                process_chunks(data)

    else:
        for n in chromosomes:
            with open(name.format(n), 'rb') as f:
                data = pickle.load(f)
            print(f'Loaded data for chromosome: {n}')
            process_chunks(data)

    for i in range(N):
        delta[i, i] = 0

    for i in range(N):
        for j in range(i + 1, N):
            delta[i, j] = delta[j, i] = tot_dists[i, j] / tot_norms[i, j]
    delta = delta ** (1 / pp)
    if TESTING:
        plt.pcolor(delta)
        plt.show()
    with open(out_name, 'wb') as f:
        pickle.dump(delta, f)

    # bootstraping ---------------

    clen = len(chunk_data)
    print(f'NUMBER OF BLOCKS: {clen}')
    for boot in range(BOOTRUNS):
        chunk_res = chunk_data.copy()
        for i in range(clen):
            chunk_res[i] = chunk_data[randint(0, clen - 1)]
        delta = np.zeros((N, N))
        tot_dists = np.sum(np.array([c[0] for c in chunk_res]), axis=0)
        tot_norms = np.sum(np.array([c[1] for c in chunk_res]), axis=0)
        for i in range(N):
            for j in range(i + 1, N):
                delta[i, j] = delta[j, i] = tot_dists[i, j] / tot_norms[i, j]
        delta = delta ** (1 / pp)
        if TESTING:
            plt.pcolor(delta)
            plt.show()
        with open(out_name + f'.boot.{boot}', 'wb') as f:
            pickle.dump(delta, f)

    print('Distance matrix computed')
    if TESTING:
        plt.pcolor(delta)
        plt.colorbar()
        plt.xlabel('i individual')
        plt.ylabel('j individual')
        plt.title('Distance matrix')
        flat = delta.reshape(N ** 2, )
        plt.figure()
        plt.hist(flat, 50)
        plt.title('Distances')
        plt.figure()
        norm_counts = tot_norms.reshape((N ** 2,))
        plt.hist(norm_counts, 50)
        plt.title('SNP counts between pairs')
        plt.figure()
        plt.pcolor(tot_dists)
        plt.colorbar()
        plt.title('Unnormalized distances')
        plt.figure()
        plt.pcolor(tot_norms)
        plt.colorbar()
        plt.title('Norms')
        plt.show()


def inv_filter(p1, p2, p3):  # This should be moved to other place actually
    M = np.eye(3)
    for i in range(3):
        M[i, 0] += 2 * p1 * p2
        M[i, 1] += 2 * p1 * p3
        M[i, 2] += 2 * p2 * p3
    M[0, 0] -= p1 + p2
    M[1, 1] -= p1 + p3
    M[2, 2] -= p2 + p3
    M[0, 1] -= p3
    M[0, 2] -= p3
    M[1, 0] -= p2
    M[1, 2] -= p2
    M[2, 0] -= p1
    M[2, 1] -= p1
    return np.matrix(M) ** -1


if __name__ == '__main__':
    pp = int(sys.argv[2])
    name = sys.argv[1]  # C:\Users\levkivskyi\PycharmProjects\medeas\test\...
    out_name = f'out.pp.{pp}.asd'

    asd_main(pp, name, out_name)
