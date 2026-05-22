#!/usr/bin/env python3
"""Convert medeas internal format directly to PLINK BED/BIM/FAM.

Medeas internal format:
- one SNP per row
- one haploid individual per column
- values: 0 = missing, 1 = reference allele, 2 = alternative allele

Each haploid column is converted to a homozygous diploid sample (N -> N).
Creates a temporary PLINK transposed dataset, then always emits only:
- <prefix>.bed
- <prefix>.bim
- <prefix>.fam
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor


ALLELE_MAP = {
    "0": "0",
    "1": "A",
    "2": "B",
}


def read_labels(labels_file):
    with open(labels_file) as f:
        labels = [line.strip() for line in f if line.strip()]
    if not labels:
        raise ValueError("Label file is empty")
    return labels


def write_tfam(labels, tfam_path, iids=None):
    with open(tfam_path, "w") as f:
        for sample_index, label in enumerate(labels, start=1):
            fid = label
            iid = iids[sample_index - 1] if iids is not None else f"{label}_{sample_index}"
            f.write(f"{fid} {iid} 0 0 0 -9\n")


def convert_row_to_tped_line(snp_index, line, expected_haploids, chrom):
    line = line.strip()
    if not line:
        return None
    genotypes = line.split()
    if len(genotypes) != expected_haploids:
        raise ValueError(
            f"Row {snp_index} has {len(genotypes)} columns but label file has {expected_haploids} haploid labels"
        )

    fields = []
    for geno in genotypes:
        try:
            allele = ALLELE_MAP[geno]
        except KeyError as exc:
            raise ValueError(
                f"Unexpected genotype code '{exc.args[0]}'; expected only 0, 1, or 2"
            ) from exc
        # homozygous diploid: repeat the same allele twice
        fields.extend([allele, allele])

    marker_id = f"snp{snp_index}"
    position = snp_index
    prefix_fields = [str(chrom), marker_id, "0", str(position)]
    return " ".join(prefix_fields + fields) + "\n"


def write_tped_multithread(snps_file, tped_path, expected_haploids, chrom, threads):
    with open(snps_file) as fin:
        lines = fin.readlines()

    worker_count = os.cpu_count() if threads == 0 else threads
    if worker_count is None or worker_count < 1:
        worker_count = 1

    with open(tped_path, "w") as fout:
        if worker_count == 1:
            for snp_index, line in enumerate(lines, start=1):
                tped_line = convert_row_to_tped_line(snp_index, line, expected_haploids, chrom)
                if tped_line is not None:
                    fout.write(tped_line)
            return

        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            futures = [
                executor.submit(convert_row_to_tped_line, idx, line, expected_haploids, chrom)
                for idx, line in enumerate(lines, start=1)
            ]
            for future in futures:
                tped_line = future.result()
                if tped_line is not None:
                    fout.write(tped_line)


def make_bed_from_tfile(tfile_prefix, out_prefix, plink_path, threads):
    if shutil.which(plink_path) is None:
        raise FileNotFoundError(f"PLINK executable not found: {plink_path}")

    cmd = [
        plink_path,
        "--tfile", tfile_prefix,
        "--make-bed",
        "--threads", str((os.cpu_count() if threads == 0 else threads) or 1),
        "--out", out_prefix,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.stderr.write(result.stderr)
        raise RuntimeError(f"PLINK --make-bed failed (exit {result.returncode})")


def convert_medeas_to_plink(snps_file, labels_file, out_prefix, chrom, plink_path, threads):
    out_dir = os.path.dirname(out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    if threads < 0:
        raise ValueError("threads must be >= 0 (0 means all cores)")

    if labels_file is not None:
        labels = read_labels(labels_file)
        iids = None
    else:
        with open(snps_file) as f:
            n = len(f.readline().split())
        labels = ["pop1"] * n
        iids = [f"ind_{i}" for i in range(1, n + 1)]

    expected_haploids = len(labels)

    with tempfile.TemporaryDirectory(prefix="medeas_to_plink_") as tmpdir:
        tmp_prefix = os.path.join(tmpdir, "tmp_plink")
        tfam_path = tmp_prefix + ".tfam"
        tped_path = tmp_prefix + ".tped"

        write_tfam(labels, tfam_path, iids=iids)
        write_tped_multithread(snps_file, tped_path, expected_haploids, chrom, threads)

        make_bed_from_tfile(tmp_prefix, out_prefix, plink_path, threads)


def main():
    parser = argparse.ArgumentParser(
        description="Convert haploid medeas files to homozygous diploid PLINK BED/BIM/FAM files (N -> N)"
    )
    parser.add_argument("--snps", required=True, metavar="FILE",
                        help="SNP file in MEDEAS format (one SNP per row, space-separated integers, one column per individual).")
    parser.add_argument("--labels", metavar="FILE",
                        help="Medeas label file (one population per line). "
                             "If omitted, all individuals are assigned to a single population 'pop1' "
                             "and named consecutively (ind_1, ind_2, ...).")
    parser.add_argument("--out", required=True, metavar="PREFIX",
                        help="Output PLINK prefix (writes .bed/.bim/.fam)")
    parser.add_argument("--chrom", default="1",
                        help="Chromosome code to use in TPED conversion (default: 1)")
    parser.add_argument("--plink-path", default="plink", metavar="FILE",
                        help="Path to the PLINK executable (default: plink)")
    parser.add_argument("-t", "--threads", type=int, default=1, metavar="N",
                        help="Worker threads for conversion (0 = all cores)")
    args = parser.parse_args()

    try:
        convert_medeas_to_plink(
            args.snps,
            args.labels,
            args.out,
            args.chrom,
            args.plink_path,
            args.threads,
        )
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
