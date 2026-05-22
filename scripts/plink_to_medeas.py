#!/usr/bin/env python3

import argparse
import os
from collections import deque
from concurrent.futures import ThreadPoolExecutor
import numpy as np


def _convert_chunk_to_pseudohaploid(raw_chunk, n_samples, seed):
    """Convert one BED SNP chunk to one pseudo-haploid column per sample."""
    current = raw_chunk.shape[0]
    rng = np.random.default_rng(seed)

    bits = np.unpackbits(raw_chunk, axis=1, bitorder="little")
    bits = bits[:, : 2 * n_samples].reshape(current, n_samples, 2)
    codes = bits[:, :, 0] + 2 * bits[:, :, 1]

    data = np.zeros((current, n_samples), dtype=np.int8)

    ## homozygote
    data[codes == 0] = 1
    data[codes == 3] = 2

    ## heterozygote
    mask_het = (codes == 2)
    data[mask_het] = rng.integers(1, 3, size=np.count_nonzero(mask_het), dtype=np.int8)
    return data


def convert_plink_to_medeas(bfile_prefix, snp_file, labels_file, threads=1):
    """Convert PLINK binary files (.bed/.bim/.fam) to pseudo-haploid medeas format.

    Each diploid individual contributes one pseudo-haploid genotype column.
    Homozygous states are kept as-is, and heterozygous states are sampled
    uniformly as 1 or 2.

    Output format (snp_file):
      - One SNP per row, Individual per column, space-separated integers.
      - Values: 0=missing, 1=ref allele (A1), 2=alt allele (A2).
    """
    fam_file = bfile_prefix + ".fam"
    bim_file = bfile_prefix + ".bim"
    bed_file = bfile_prefix + ".bed"

    snp_dir = os.path.dirname(snp_file)
    if snp_dir:
        os.makedirs(snp_dir, exist_ok=True)

    for path in (fam_file, bim_file, bed_file):
        if not os.path.isfile(path):
            raise FileNotFoundError(f"PLINK file not found: {path}")

    family_ids = []
    with open(fam_file) as f:
        for line in f:
            parts = line.strip().split()
            if parts:
                family_ids.append(parts[0])
    n_samples = len(family_ids)

    with open(bim_file) as f:
        n_snps = sum(1 for line in f if line.strip())

    print(f"Converting PLINK data: {n_samples} samples, {n_snps} SNPs")

    n_bytes_per_snp = (n_samples + 3) // 4
    chunk_snps = 2000
    rng = np.random.default_rng()
    if threads < 0:
        raise ValueError("threads must be >= 0")
    effective_threads = os.cpu_count() if threads == 0 else threads
    if effective_threads is None or effective_threads < 1:
        effective_threads = 1

    print(f"Using {effective_threads} thread(s)")

    print(f"Writing: {snp_file}")
    with open(bed_file, "rb") as bed, open(snp_file, "w") as out:
        magic = bed.read(3)
        if magic != b"\x6c\x1b\x01":
            raise ValueError(
                f"Not a valid PLINK BED file (unexpected magic bytes): {bed_file}"
            )

        snp_done = 0
        if effective_threads == 1:
            while snp_done < n_snps:
                current = min(chunk_snps, n_snps - snp_done)
                raw = np.frombuffer(bed.read(current * n_bytes_per_snp), dtype=np.uint8)
                if raw.size != current * n_bytes_per_snp:
                    raise ValueError("Unexpected end of BED file while reading SNP data")
                raw = raw.reshape(current, n_bytes_per_snp)

                seed = int(rng.integers(0, np.iinfo(np.int64).max, dtype=np.int64))
                geno_pseudohaploid = _convert_chunk_to_pseudohaploid(raw, n_samples, seed)

                np.savetxt(out, geno_pseudohaploid, fmt="%d", delimiter=" ")
                snp_done += current
                if snp_done % 20000 == 0 or snp_done == n_snps:
                    print(f"Converted {snp_done}/{n_snps} SNPs")
        else:
            pending = deque()
            max_pending = max(2, effective_threads * 2)
            with ThreadPoolExecutor(max_workers=effective_threads) as executor:
                while snp_done < n_snps:
                    current = min(chunk_snps, n_snps - snp_done)
                    raw = np.frombuffer(bed.read(current * n_bytes_per_snp), dtype=np.uint8)
                    if raw.size != current * n_bytes_per_snp:
                        raise ValueError("Unexpected end of BED file while reading SNP data")
                    raw = raw.reshape(current, n_bytes_per_snp)

                    seed = int(rng.integers(0, np.iinfo(np.int64).max, dtype=np.int64))
                    future = executor.submit(_convert_chunk_to_pseudohaploid, raw, n_samples, seed)
                    pending.append((current, future))

                    if len(pending) >= max_pending:
                        written_current, written_future = pending.popleft()
                        np.savetxt(out, written_future.result(), fmt="%d", delimiter=" ")
                        snp_done += written_current
                        if snp_done % 20000 == 0 or snp_done == n_snps:
                            print(f"Converted {snp_done}/{n_snps} SNPs")

                while pending:
                    written_current, written_future = pending.popleft()
                    np.savetxt(out, written_future.result(), fmt="%d", delimiter=" ")
                    snp_done += written_current
                    if snp_done % 20000 == 0 or snp_done == n_snps:
                        print(f"Converted {snp_done}/{n_snps} SNPs")


    if labels_file is not None:
        print(f"Writing: {labels_file}")
        labels_dir = os.path.dirname(labels_file)
        if labels_dir:
            os.makedirs(labels_dir, exist_ok=True)

        with open(labels_file, "w") as f:
            for fid in family_ids:
                f.write(fid + "\n")

    print("PLINK conversion complete.")


def main():
    parser = argparse.ArgumentParser(
        description="Convert diploid PLINK .bed/.bim/.fam files to pseudo-haploid medeas SNP/labels files (N -> N)"
    )
    parser.add_argument("--bfile", required=True, metavar="PREFIX",
                        help="PLINK file prefix (without .bed/.bim/.fam)")
    parser.add_argument("--out", required=True, metavar="FILE",
                        help="Output SNP matrix file (medeas format)")
    parser.add_argument("--out_labels", default=None, metavar="FILE",
                        help="Output labels file; if given, population labels from the .fam file (first column) are written to this file")
    parser.add_argument("-t", "--threads", type=int, default=1, metavar="N",
                        help="Number of worker threads for chunk conversion (0 = all cores)")
    args = parser.parse_args()

    convert_plink_to_medeas(args.bfile, args.out, args.out_labels, threads=args.threads)


if __name__ == "__main__":
    main()
