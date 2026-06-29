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
import subprocess
import os
import pandas as pd


def _iter_chunks_from_medeas(name: str):
    """Yield genotype chunks from a MEDEAS text file.

    Each line is one SNP; values are space-separated integers in {0, 1, 2}.
    Yields (data, done, total) where total is -1 (unknown ahead of time).
    """
    MAXSIZE = 5000 * 2 ** 20  # 5 GB read buffer
    processed = 0
    with open(name) as f:
        while True:
            lines = f.readlines(MAXSIZE)
            if not lines:
                break
            processed += len(lines)
            data = np.array(
                [np.fromstring(line, sep=' ', dtype='int8') for line in lines]
            )
            yield data, processed, -1


def _iter_chunks_from_plink(bfile_prefix: str, chunk_snps: int = 2000):
    """Yield pseudo-haploid genotype chunks from a PLINK BED prefix.

    Output chunk shape is (n_snps_chunk, n_samples), values in {0,1,2}.
    Each diploid sample contributes a single pseudo-haploid genotype.
    """
    fam_file = bfile_prefix + ".fam"
    bim_file = bfile_prefix + ".bim"
    bed_file = bfile_prefix + ".bed"

    ## medeas format:
    ## 0: missing data.
    ## 1: reference allele.
    ## 2: alternative allele.

    ## plink format (*.bed):
    ## Bits 00: homozygous A1
    ## Bits 01: missing
    ## Bits 10: heterozygous
    ## Bits 11: homozygous A2

    ## conversion
    ## code 0 -> 1, 
    ## code 1 -> 0 (missing).
    ## code 2 -> random 1 or 2, 
    ## code 3 -> 2, 

    n_samples = sum(1 for line in open(fam_file) if line.strip())
    n_snps = sum(1 for line in open(bim_file) if line.strip())
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

            data = np.zeros((current, n_samples), dtype=np.int8)

            ## homozygous reference allele
            data[codes == 0] = 1

            ## homozygous alternative allele
            data[codes == 3] = 2

            ## heterozygote: randomly assign to reference or alternative allele
            mask_het = (codes == 2)
            data[mask_het] = rng.integers(1, 3, size=np.count_nonzero(mask_het), dtype=np.int8)

            done += current
            yield data, done, n_snps


def _iter_chunks_from_vcf(variant_file: str, mode: str = "vcf_random", chunk_snps: int = 2000):
    """Yield pseudo-haploid genotype chunks from VCF/BCF (gzipped or plain).

    Modes:
    - vcf_random: always pseudo-haploidize by random allele selection
    - vcf_phased1: if phased use haplotype 1, else random
    - vcf_phased2: if phased use haplotype 2, else random
    """
    try:
        import pysam
    except ImportError as exc:
        raise ImportError(
            "pysam is required for --vcf/--vcf_phased1/--vcf_phased2 input. "
            "Install it with 'pip install pysam'."
        ) from exc

    if mode not in {"vcf_random", "vcf_phased1", "vcf_phased2"}:
        raise ValueError(f"Unsupported variant decoding mode: {mode}")

    rng = np.random.default_rng()
    vf = pysam.VariantFile(variant_file)
    n_samples = len(vf.header.samples)
    if n_samples == 0:
        vf.close()
        raise ValueError(f"Variant file has no sample columns: {variant_file}")

    total_snps = -1  # record count unknown without scanning the whole file

    chunk = np.empty((chunk_snps, n_samples), dtype=np.int8)
    row_i = 0
    done = 0

    if mode == "vcf_random":
        def _pick_heterozygote(gt, sample, idx, rp):
            if rp[0] is None:
                rp[0] = rng.integers(0, 2, size=n_samples, dtype=np.int8)
            return gt[int(rp[0][idx])]
    elif mode == "vcf_phased1":
        def _pick_heterozygote(gt, sample, idx, rp):
            if sample.phased:
                return gt[0] ## haploid 1
            if rp[0] is None:
                rp[0] = rng.integers(0, 2, size=n_samples, dtype=np.int8)
            return gt[int(rp[0][idx])]
    elif mode == "vcf_phased2":
        def _pick_heterozygote(gt, sample, idx, rp):
            if sample.phased:
                return gt[1] ## haploid 2
            if rp[0] is None:
                rp[0] = rng.integers(0, 2, size=n_samples, dtype=np.int8)
            return gt[int(rp[0][idx])]

    for rec in vf:
        row = chunk[row_i]
        rp = [None]  # lazily filled once per record for random/fallback picks
        for idx, sample in enumerate(rec.samples.values()):
            gt = sample.get("GT")
            if gt is None or len(gt) == 0:
                allele = None
            elif len(gt) == 1:
                allele = gt[0]
            else:
                allele = _pick_heterozygote(gt, sample, idx, rp)

            if allele is None:
                row[idx] = 0
            elif allele == 0:
                row[idx] = 1
            else:
                row[idx] = 2

        row_i += 1
        done += 1

        if row_i == chunk_snps:
            emit = chunk
            chunk = np.empty((chunk_snps, n_samples), dtype=np.int8)
            yield emit, done, total_snps
            row_i = 0

    vf.close()

    if row_i > 0:
        yield chunk[:row_i], done, total_snps



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


def ibs_distance_plink(bfile: str, dist_out: str, plink_path: str,
                        plink_params: str, threads: int) -> str:
    """Compute pairwise IBS distance file using PLINK --distance."""
    prefix = dist_out
    cmd = [
        plink_path,
        "--bfile", bfile,
        "--allow-no-sex",
        "--distance", "gz", "1-ibs", "square", "flat-missing",
        "--threads", str(threads),
        "--out", prefix,
    ]
    if plink_params:
        cmd += plink_params.split()

    print(f"Running: {' '.join(cmd)}", file=sys.stderr)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.stderr.write(result.stderr)
        raise RuntimeError(f"PLINK distance computation failed (exit {result.returncode})")

    dist_file = f"{prefix}.mdist.gz"  # gz + square → gzipped text square matrix
    if os.path.isfile(dist_file):
        return dist_file

    raise RuntimeError(
        "PLINK distance output file not found after successful run. "
        f"Checked: {dist_file}"
    )

def _compute_asd_bootstraps(simulation, chunk_data: list, N: int) -> None:
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

def  compute_asd_matrix_plink(simulation, save_pickle = False) -> None:
    """Calculates the same distance matrix as compute_asd_matrix() with plink.
    Take care that here we are assuming haploid individuals!
    'snp_file' is the input file name (or format string) with SNP data.
    'out_name' is the output binary file (pickle).
    """
    print('Compute distance matrix with PLINK')
    bfile_prefix = simulation.snp_file[:-4]  # remove .bed suffix

    ## compute distance with plink
    pp = 1
    out_name =  simulation.asd_pattern.format(pp)
    delta1_path = ibs_distance_plink(bfile_prefix, out_name, simulation.plink_path, "", simulation.NCORE)

    ## remove intermediate PLINK output files
    for ext in ('.mdist.id', '.nosex', '.log'):
        p = out_name + ext
        if os.path.isfile(p):
            os.remove(p)

    ## read delta1 to compute delta2
    pp = 2
    delta1 = pd.read_csv(delta1_path, sep=r'\s+', header=None).to_numpy()
    delta2 = delta1 ** (1 / pp)
    np.savetxt(simulation.asd_pattern.format(pp) + ".mdist.gz", delta2)

    if save_pickle:
        with open(simulation.asd_pattern.format(1), 'wb') as f:
            pickle.dump(dist_matrix, f)
        with open(simulation.asd_pattern.format(2), 'wb') as f:
            pickle.dump(delta2, f)


def  compute_asd_matrix(simulation) -> None:
    """Calculate the distance matrix for two Minkowsky parameter pp.
    Take care that here we are assuming haploid individuals!
    'snp_file' is the input file name (or format string) with SNP data.
    'out_name' is the output binary file (pickle).
    """

    # On Windows, processes execute the whole file before forking
    # therefore we protect this code with if __name__ == '__main__'
    # Need to think how to avoid copying ``data`` on forking.
    # Maybe process input file in chunks?
    # On POSIX everything is already fine because of "copy-on-write"

    snp_file = simulation.snp_file
    snp_file_type = simulation.snp_file_type

    bootsize = simulation.bootsize
    no_split = simulation.no_split
    dist_func = lambda x: np.abs(x)

    # For no-split runs we do not need bootstrap chunks, so use larger
    # windows to reduce Python overhead and improve BLAS throughput.
#    if no_split:
#        bootsize = max(bootsize, 5000)

    # ---------- global data

    print('Distance matrix computing started')

    N = len(simulation.labels)

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
    if snp_file_type == "plink":
        bfile_prefix = snp_file[:-4]  # remove .bed suffix
        print(f'Loading data from PLINK files with prefix: {bfile_prefix}')
        chunk_iter = _iter_chunks_from_plink(bfile_prefix)
    elif snp_file_type.split("_")[0] == "vcf":
        print(f'Loading data from VCF file: {snp_file}')
        chunk_iter = _iter_chunks_from_vcf(snp_file, mode=snp_file_type)
    else:
        print(f'Loading data from MEDEAS text file: {snp_file}')
        chunk_iter = _iter_chunks_from_medeas(snp_file)

    for data, done, total in chunk_iter:
        total_str = str(total) if total != -1 else "?"
        print(f'   Chunk loaded ({data.shape[0]} sites, {data.shape[1]} individuals) [{done}/{total_str}]')

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
    print(f'Saving ASD matrix with pp={pp} to: {out_name}')
    with open(out_name, 'wb') as f:
        pickle.dump(delta_1, f)
    #np.savetxt(out_name + '.txt', delta_1, fmt='%.6f')
    pp = 2
    delta_2 = np.copy(delta)
    delta_2 = delta_2 ** (1 / pp)
    out_name =  simulation.asd_pattern.format(pp)
    with open(out_name, 'wb') as f:
        pickle.dump(delta_2, f)
    #np.savetxt(out_name + '.txt', delta_2, fmt='%.6f')

    if not no_split:
        _compute_asd_bootstraps(simulation, chunk_data, N)

    print('Distance matrix computed')
    # bootstrapping ---------------





