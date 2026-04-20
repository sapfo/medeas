# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 15:29:37 2016

@author: ivan
"""

import sys
import numpy as np
from multiprocessing import Process, Queue
from typing import Tuple, Callable, List, IO
import pickle
from random import randint
import matplotlib.pyplot as plt


def _iter_haploid_chunks_from_bed(bfile_prefix: str, chunk_snps: int = 2000, max_snps: int = None):
    """Yield haploid-coded genotype chunks from a PLINK BED prefix.

    Output chunk shape is (n_snps_chunk, 2*n_samples), values in {0,1,2}.
    """
    fam_file = bfile_prefix + ".fam"
    bim_file = bfile_prefix + ".bim"
    bed_file = bfile_prefix + ".bed"

    n_samples = sum(1 for line in open(fam_file) if line.strip())
    n_snps = sum(1 for line in open(bim_file) if line.strip())
    if max_snps is not None:
        n_snps = min(n_snps, max_snps)
    n_haploid = 2 * n_samples
    n_bytes_per_snp = (n_samples + 3) // 4

    rng = np.random.default_rng()
    with open(bed_file, "rb") as bed:
        magic = bed.read(3)
        if magic != b"\x6c\x1b\x01":
            raise ValueError(
                f"Not a valid PLINK BED file (unexpected magic bytes): {bed_file}"
            )

        done = 0
        while done < n_snps:
            current = min(chunk_snps, n_snps - done)
            raw = np.frombuffer(bed.read(current * n_bytes_per_snp), dtype=np.uint8)
            if raw.size != current * n_bytes_per_snp:
                raise ValueError("Unexpected end of BED file while reading SNP data")
            raw = raw.reshape(current, n_bytes_per_snp)

            bits = np.unpackbits(raw, axis=1, bitorder="little")
            bits = bits[:, : 2 * n_samples].reshape(current, n_samples, 2)
            codes = bits[:, :, 0] + 2 * bits[:, :, 1]

            hap1 = np.zeros((current, n_samples), dtype=np.int8)
            hap2 = np.zeros((current, n_samples), dtype=np.int8)

            mask_hom_ref = (codes == 0)
            hap1[mask_hom_ref] = 1
            hap2[mask_hom_ref] = 1

            mask_hom_alt = (codes == 3)
            hap1[mask_hom_alt] = 2
            hap2[mask_hom_alt] = 2

            mask_het = (codes == 2)
            flip = rng.random((current, n_samples)) > 0.5
            hap1[mask_het & ~flip] = 1
            hap2[mask_het & ~flip] = 2
            hap1[mask_het & flip] = 2
            hap2[mask_het & flip] = 1

            data = np.empty((current, n_haploid), dtype=np.int8)
            data[:, 0::2] = hap1
            data[:, 1::2] = hap2
            done += current
            yield data, done, n_snps





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
    """Compute distances and norms for rows from 'tasks'."""
    while True:
        i = tasks.get()
        if i < 0:
            return
        results.put((i, compute(i, data, dist_func, N)))


def process(data: 'np.ndarray[int]',
            dist_func: Callable[[np.ndarray], np.ndarray],
            tasks, results, N, NPROC) -> Tuple['np.ndarray[float]', 'np.ndarray[int]']:
    """Calculate matrices of un-normalized distances and norms for 'data'
    using given distance function.
    """
    # Fast path for the medeas haploid encoding {0,1,2}:
    # - norm(i,j): number of sites where both are non-missing
    # - dist(i,j): mismatching non-missing sites (|1-2|=1)
    data_i16 = data.astype(np.int16, copy=False)
    if np.all((data_i16 == 0) | (data_i16 == 1) | (data_i16 == 2)):
        non_missing = (data_i16 != 0).astype(np.int32)
        norms_full = non_missing @ non_missing.T

        is1 = (data_i16 == 1).astype(np.int32)
        is2 = (data_i16 == 2).astype(np.int32)
        matches_full = (is1 @ is1.T) + (is2 @ is2.T)
        dists_full = norms_full - matches_full

        dists = np.triu(dists_full, k=1).astype(float, copy=False)
        norms = np.triu(norms_full, k=1).astype(float, copy=False)
        return dists, norms

    # Fallback for other encodings.
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
    """Calculate the distance matrix for two Minkowsky parameter pp.
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
    bfile_prefix = getattr(simulation, "bfile_prefix", None)
    bootsize = simulation.bootsize
    no_split = simulation.no_split
    max_snps = getattr(simulation, "max_snps", None)
    dist_func = lambda x: np.abs(x)

    # For no-split runs we do not need bootstrap chunks, so use larger
    # windows to reduce Python overhead and improve BLAS throughput.
#    if no_split:
#        bootsize = max(bootsize, 5000)

    # ---------- constants

    # this should be large to avoid overhead of spawning new processes
    # or we need to reuse them somehow
    MAXSIZE = 5000 * 2 ** 20  # 5 Gb


    # ---------- global data

    print('Distance matrix computing started')

    if bfile_prefix is not None:
        fam_file = bfile_prefix + ".fam"
        n_samples = sum(1 for line in open(fam_file) if line.strip())
        N = 2 * n_samples
    else:
        with open(name) as f:
            N = len(f.readline()) // 2
            # print(f'nb individual = {N}')

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
        if remainder is not None and len(remainder) > 0:
            data = np.vstack((remainder, data))
        remainder = None

        start_i = 0
        while start_i + bootsize <= len(data):
            if simulation.output_level >= 1:
                print(f'Processing site {start_i}')
            chunk = data[start_i:start_i + bootsize]
            datac = chunk.T.copy()

            dists, norms = process(datac, dist_func, tasks, results, N, simulation.NCORE)
            tot_dists += dists
            tot_norms += norms
            if not no_split:
                chunk_data.append((dists, norms))
            start_i += bootsize

        print(f'Chunk processing done, {len(data) - start_i} sites remain for next chunk')

        remainder = data[start_i:]

    # Frequency histogram for folded SFS to avoid repeated np.append reallocations.
    sfs_counts = np.zeros((N + 1,), dtype=np.int64)
    if bfile_prefix is not None:
        for data, done, total in _iter_haploid_chunks_from_bed(bfile_prefix, max_snps=max_snps):
            print(f'   Chunk loaded ({data.shape[0]} sites, {data.shape[1]} individuals) [{done}/{total}]')

            nb_mut = np.sum(data == 1, axis=1)
            nb_missing = np.sum(data == 0, axis=1)
            nb_mut[np.where((N - nb_missing) == 0)] = 0
            nb_missing[np.where((N - nb_missing) == 0)] = 1
            freq = nb_mut / (N - nb_missing)
            nb_other_mut = np.random.binomial(nb_missing, freq, len(nb_mut))
            nb_mut = nb_mut + nb_other_mut
            nb_mut[nb_mut > N / 2] = N - nb_mut[nb_mut > N / 2]
            vals, cnts = np.unique(nb_mut.astype(int), return_counts=True)
            sfs_counts[vals] += cnts
            process_chunks(data)
    else:
        processed_snps = 0
        with open(name) as f:
            while True:
                if max_snps is not None and processed_snps >= max_snps:
                    break
                data_lines = f.readlines(MAXSIZE)
                if not data_lines:
                    break
                if max_snps is not None:
                    remaining = max_snps - processed_snps
                    if remaining <= 0:
                        break
                    if len(data_lines) > remaining:
                        data_lines = data_lines[:remaining]
                    processed_snps += len(data_lines)
                print('   Chunk loading started')
                data = np.array([np.fromstring(line, sep=' ',  # line[-cut:-1]
                                               dtype='int8')
                                 for line in data_lines])
                # data = data[:, ::2] + data[:, 1::2]
                print(f'   Chunk loaded ({len(data_lines)} sites, {data.shape[1]} individuals)')

                nb_mut = np.sum(data==1,axis=1)
                nb_missing = np.sum(data==0,axis=1)
                nb_mut[np.where((N - nb_missing)==0)] = 0 ## removing in a stupid way site with no data at all
                nb_missing[np.where((N - nb_missing)==0)] = 1 ## removing in a stupid way site with no data at all
                freq = nb_mut/(N - nb_missing)
                nb_other_mut = np.random.binomial(nb_missing,freq,len(nb_mut))
                nb_mut = nb_mut + nb_other_mut
                nb_mut[nb_mut > N/2] = N - nb_mut[nb_mut > N/2] # Folding the SFS, since 0 and 1 are likely to be not well defined
                vals, cnts = np.unique(nb_mut.astype(int), return_counts=True)
                sfs_counts[vals] += cnts
                process_chunks(data)

    # Process trailing sites that did not fill a full window.
    if remainder is not None and len(remainder) > 0:
        datac = remainder.T.copy()
        dists, norms = process(datac, dist_func, tasks, results, N, simulation.NCORE)
        tot_dists += dists
        tot_norms += norms
        if not no_split:
            chunk_data.append((dists, norms))

    freq_values = np.nonzero(sfs_counts)[0]
    simulation.sfs = np.array((freq_values, sfs_counts[freq_values]))

    upper = np.triu_indices(N, k=1)
    valid = tot_norms[upper] > 0
    delta_vals = np.zeros_like(tot_dists[upper], dtype=float)
    delta_vals[valid] = tot_dists[upper][valid] / tot_norms[upper][valid]
    delta[upper] = delta_vals
    delta[(upper[1], upper[0])] = delta_vals
    np.fill_diagonal(delta, 0)

    zero_norm_pairs = int(np.size(valid) - np.count_nonzero(valid))
    if zero_norm_pairs > 0 and simulation.output_level >= 1:
        print(f'Warning: {zero_norm_pairs} pairs had zero overlap; distances set to 0.')
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


    # bootstrapping ---------------

    if not no_split:
        clen = len(chunk_data)
        print(f'NUMBER OF BLOCKS: {clen}')
        for boot in range(simulation.bootstrap_number):
            chunk_res = chunk_data.copy()
            for i in range(clen):
                chunk_res[i] = chunk_data[randint(0, clen - 1)]
            delta = np.zeros((N, N))
            tot_dists = np.sum(np.array([c[0] for c in chunk_res]), axis=0)
            tot_norms = np.sum(np.array([c[1] for c in chunk_res]), axis=0)
            upper = np.triu_indices(N, k=1)
            valid = tot_norms[upper] > 0
            delta_vals = np.zeros_like(tot_dists[upper], dtype=float)
            delta_vals[valid] = tot_dists[upper][valid] / tot_norms[upper][valid]
            delta[upper] = delta_vals
            delta[(upper[1], upper[0])] = delta_vals
            np.fill_diagonal(delta, 0)

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

