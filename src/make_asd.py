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

TESTING = True

if TESTING:
    import matplotlib.pyplot as plt

N = 294


def dist_and_norm(a: 'np.ndarray[int]', b: 'np.ndarray[int]',
                  dist_func: Callable[[np.ndarray], np.ndarray]
                  ) -> Tuple[float, int]:
    """Calculate un-normalized distance and norm between vectors 'a' and 'b'.
    Norm is number of sites 'i' where both 'a[i]' and 'b[i]' are non-zero.
    Other sites does not contribute to distance.
    """
    filt = np.logical_and(a, b)
    dst = filt * dist_func(a-b)
    return np.sum(dst), np.sum(filt)


def compute(i: int, data: 'np.ndarray[int]',
            dist_func: Callable[[np.ndarray], np.ndarray]
            ) -> Tuple[List[float], List[int]]:
    """Compute all distances and norms for 'i'th row in 'data'."""
    dists: List[float] = []
    norms: List[int] = []
    for j in range(i+1, N):
        dist, norm = dist_and_norm(data[i], data[j], dist_func)
        dists.append(dist)
        norms.append(norm)
    return dists, norms


def work(tasks: 'Queue[int]',
         results: 'Queue[Tuple[int, Tuple[List[float], List[int]]]]',
         data: 'np.ndarray[int]', dist_func: Callable[[np.ndarray], np.ndarray]
         ) -> None:
    """Compute distaces and norms for rows from 'tasks'."""
    while True:
        i = tasks.get()
        if i < 0:
            return
        results.put((i, compute(i, data, dist_func)))


def process(data: 'np.ndarray[int]',
            dist_func: Callable[[np.ndarray], np.ndarray]
            ) -> Tuple['np.ndarray[float]', 'np.ndarray[int]']:
    """Calculate matrices of un-normalized distances and norms for 'data'
    using given distance function.
    """
    dists = np.zeros((N, N))
    norms = np.zeros((N, N))

    for i in range(N):
        tasks.put(i)
    for _ in range(NPROC):
        tasks.put(-1)

    procs = [Process(target=work, args=(tasks, results, data, dist_func))
             for _ in range(NPROC)]
    for proc in procs:
        proc.start()

    rest = N
    while rest:
        i, (dist, norm) = results.get()
        rest -= 1
        dists[i, i+1:] = dist
        norms[i, i+1:] = norm

    for proc in procs:
        proc.join()

    return dists, norms


if __name__ == '__main__':

    # On Windows, processes execute the whole file before forking
    # therefore we protect this code with if __name__ == '__main__'
    # Need to think how to avoid copying ``data`` on forking.
    # Maybe process input file in chunks?
    # On POSIX everything is already fine because of "copy-on-write"

    def test_func(dst: int) -> int:
        if dst > 1:
            return dst**2
        else:
            return 0

    dist_func = np.vectorize(test_func)
    dist_func = np.abs
    dist_func = np.square
    dist_func = lambda x: np.abs(x)**2

    # ---------- constants

    N = 294
    sites = 24112  # TODO: read this from file
    cut = N*4
    name = sys.argv[1]  # C:\Users\levkivskyi\PycharmProjects\medeas\test\...

    NPROC = cpu_count()

    # this should be large to avoid overhead of spawning new processes
    # or we need to reuse them somehow
    MAXSIZE = 100*2**20  # 100 MB

    # ---------- global data

    tot_dists = np.zeros((N, N))
    tot_norms = np.zeros((N, N))
    delta = np.zeros((N, N))
    tasks = Queue()
    results = Queue()

    f: IO[str]
    with open(name) as f:
        while True:
            data_lines = f.readlines(MAXSIZE)
            if not data_lines:
                break
            data = np.array([np.fromstring(line[-cut:-1], sep=' ',
                                           dtype='int8')
                             for line in data_lines])
            data = data[:, ::2] + data[:, 1::2]
            data = data.T.copy()
            dists, norms = process(data, dist_func)
            tot_dists += dists
            tot_norms += norms

    for i in range(N):
        delta[i, i] = 0

    for i in range(N):
        for j in range(i+1, N):
            delta[i, j] = delta[j, i] = tot_dists[i, j]/tot_norms[i, j]
    delta = delta**.5
    with open('temp_asd.asd', 'wb') as f:
        pickle.dump(delta, f)

    print('Distance matrix computed')
    if TESTING:
        plt.pcolor(delta)
        plt.colorbar()
        plt.xlabel('i individual')
        plt.ylabel('j individual')
        plt.title('Distance matrix')
        flat = delta.reshape(N**2,)
        plt.figure()
        plt.hist(flat, 1250)
        plt.show()


def inv_filter(p1, p2, p3): # This should be moved to other place actually
    M = np.eye(3)
    for i in range(3):
        M[i, 0] += 2*p1*p2
        M[i, 1] += 2*p1*p3
        M[i, 2] += 2*p2*p3
    M[0, 0] -= p1+p2
    M[1, 1] -= p1+p3
    M[2, 2] -= p2+p3
    M[0, 1] -= p3
    M[0, 2] -= p3
    M[1, 0] -= p2
    M[1, 2] -= p2
    M[2, 0] -= p1
    M[2, 1] -= p1
    return np.matrix(M)**-1