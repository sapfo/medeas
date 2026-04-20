#!/usr/bin/env python3
"""Convert medeas internal format directly to PLINK BED/BIM/FAM.

Medeas internal format:
- one SNP per row
- one haploid individual per column
- values: 0 = missing, 1 = reference allele, 2 = alternative allele

This script pairs consecutive haploid columns into diploid samples, creates a
temporary PLINK transposed dataset, then always emits only:
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
    if len(labels) % 2 != 0:
        raise ValueError(
            "Medeas label file must contain an even number of haploid labels so they can be paired into diploid samples"
        )
    return labels


def write_tfam(labels, tfam_path):
    with open(tfam_path, "w") as f:
        for sample_index in range(0, len(labels), 2):
            label1 = labels[sample_index]
            label2 = labels[sample_index + 1]
            if label1 != label2:
                raise ValueError(
                    f"Expected consecutive haploid labels to match for diploid reconstruction, got '{label1}' and '{label2}' at lines {sample_index + 1}-{sample_index + 2}"
                )
            diploid_index = sample_index // 2 + 1
            fid = label1
            iid = f"{label1}_{diploid_index}"
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

    if len(genotypes) % 2 != 0:
        raise ValueError("Each SNP row must contain an even number of haploid genotype columns")

    fields = []
    for idx in range(0, len(genotypes), 2):
        try:
            a1 = ALLELE_MAP[genotypes[idx]]
            a2 = ALLELE_MAP[genotypes[idx + 1]]
        except KeyError as exc:
            raise ValueError(
                f"Unexpected genotype code '{exc.args[0]}'; expected only 0, 1, or 2"
            ) from exc
        fields.extend([a1, a2])

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

    labels = read_labels(labels_file)
    expected_haploids = len(labels)

    with tempfile.TemporaryDirectory(prefix="medeas_to_plink_") as tmpdir:
        tmp_prefix = os.path.join(tmpdir, "tmp_plink")
        tfam_path = tmp_prefix + ".tfam"
        tped_path = tmp_prefix + ".tped"

        write_tfam(labels, tfam_path)
        write_tped_multithread(snps_file, tped_path, expected_haploids, chrom, threads)

        make_bed_from_tfile(tmp_prefix, out_prefix, plink_path, threads)


def main():
    parser = argparse.ArgumentParser(
        description="Convert haploid medeas files to diploid PLINK BED/BIM/FAM files (2N -> N)"
    )
    parser.add_argument("-sf", "--snps_file", required=True,
                        help="Medeas genotype matrix (one SNP per row, haploid columns)")
    parser.add_argument("-lf", "--labels_file", required=True,
                        help="Medeas label file (one haploid label per line)")
    parser.add_argument("-o", "--out", required=True,
                        help="Output PLINK prefix (writes .bed/.bim/.fam)")
    parser.add_argument("--chrom", default="1",
                        help="Chromosome code to use in TPED conversion (default: 1)")
    parser.add_argument("--plink-path", default="plink",
                        help="Path to the PLINK executable (default: plink)")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Worker threads for conversion (0 = all cores)")
    args = parser.parse_args()

    try:
        convert_medeas_to_plink(
            args.snps_file,
            args.labels_file,
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
