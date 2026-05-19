#!/usr/bin/env python3
"""Convert MEDEAS internal format to a VCF file.

MEDEAS internal format:
- one SNP per row
- one haploid individual per column
- values: 0 = missing, 1 = reference allele, 2 = alternative allele

This script converts each haploid column into one homozygous diploid sample (N -> N):
- 0 -> ./.
- 1 -> 0/0
- 2 -> 1/1
"""

import argparse
import contextlib
import io
import os
import shutil
import subprocess
import sys


GENOTYPE_MAP = {
    "0": "./.",
    "1": "0/0",
    "2": "1/1",
}


def read_labels(labels_file):
    with open(labels_file) as f:
        labels = [line.strip() for line in f if line.strip()]
    if not labels:
        raise ValueError("Label file is empty")
    return labels


def make_unique_sample_ids(labels):
    counts = {}
    sample_ids = []
    had_duplicates = False

    for label in labels:
        counts[label] = counts.get(label, 0) + 1
        n = counts[label]
        if n == 1:
            sample_ids.append(label)
        else:
            had_duplicates = True
            sample_ids.append(f"{label}_{n}")

    if had_duplicates:
        print(
            "Warning: duplicate labels detected; sample IDs were made unique in VCF header.",
            file=sys.stderr,
        )

    return sample_ids


@contextlib.contextmanager
def open_text_output(path):
    """Open a text file for writing, using bgzip for .gz paths."""
    if path.endswith(".gz"):
        if shutil.which("bgzip") is None:
            raise RuntimeError(
                "bgzip not found on PATH. Install htslib (e.g. 'conda install -c bioconda htslib') "
                "to produce bgzip-compressed VCF files compatible with pysam/tabix."
            )
        with open(path, "wb") as raw_out:
            proc = subprocess.Popen(
                ["bgzip", "-c"],
                stdin=subprocess.PIPE,
                stdout=raw_out,
            )
            text_wrapper = io.TextIOWrapper(proc.stdin, encoding="utf-8")
            try:
                yield text_wrapper
            finally:
                text_wrapper.flush()
                proc.stdin.close()
                proc.wait()
                if proc.returncode != 0:
                    raise RuntimeError(f"bgzip exited with code {proc.returncode}")
    else:
        with open(path, "w") as f:
            yield f


def convert_row_to_vcf_line(snp_index, line, expected_haploids, chrom, ref, alt):
    stripped = line.strip()
    if not stripped:
        return None

    genotypes = stripped.split()
    if len(genotypes) != expected_haploids:
        raise ValueError(
            f"Row {snp_index} has {len(genotypes)} columns but label file has {expected_haploids} labels"
        )

    calls = []
    for geno in genotypes:
        try:
            calls.append(GENOTYPE_MAP[geno])
        except KeyError as exc:
            raise ValueError(
                f"Unexpected genotype code '{exc.args[0]}'; expected only 0, 1, or 2"
            ) from exc

    fixed_fields = [
        str(chrom),
        f"snp{snp_index}",
        ".",
        ref,
        alt,
        ".",
        "PASS",
        ".",
        "GT",
    ]
    return "\t".join(fixed_fields + calls) + "\n"


def _default_labels_out_path(out_file):
    if out_file.endswith(".vcf.gz"):
        return out_file[:-7] + ".labels.lab"
    if out_file.endswith(".vcf"):
        return out_file[:-4] + ".labels.lab"
    return out_file + ".labels.lab"


def convert_medeas_to_vcf(
    snps_file,
    labels_file,
    out_file,
    chrom="1",
    ref="A",
    alt="T",
    out_labels=None,
):
    if ref == alt:
        raise ValueError("REF and ALT alleles must be different")

    labels = read_labels(labels_file)
    expected_haploids = len(labels)
    sample_ids = make_unique_sample_ids(labels)

    out_dir = os.path.dirname(out_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    labels_out_path = out_labels if out_labels is not None else _default_labels_out_path(out_file)
    labels_out_dir = os.path.dirname(labels_out_path)
    if labels_out_dir:
        os.makedirs(labels_out_dir, exist_ok=True)

    with open(labels_out_path, "w") as flabel:
        flabel.write("\n".join(sample_ids) + "\n")

    with open_text_output(out_file) as fout:
        fout.write("##fileformat=VCFv4.2\n")
        fout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fout.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(sample_ids)
            + "\n"
        )

        with open(snps_file) as fin:
            for snp_index, line in enumerate(fin, start=1):
                vcf_line = convert_row_to_vcf_line(
                    snp_index=snp_index,
                    line=line,
                    expected_haploids=expected_haploids,
                    chrom=chrom,
                    ref=ref,
                    alt=alt,
                )
                if vcf_line is not None:
                    fout.write(vcf_line)

    return labels_out_path


def main():
    parser = argparse.ArgumentParser(
        description="Convert MEDEAS SNP/labels files to a homozygous diploid VCF file (N -> N)"
    )
    parser.add_argument(
        "--snps",
        required=True,
        help="SNP file in MEDEAS format (one SNP per row, space-separated integers)",
    )
    parser.add_argument(
        "--labels",
        required=True,
        help="Label file (one label per MEDEAS haploid column)",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output VCF path (.vcf or .vcf.gz)",
    )
    parser.add_argument(
        "--out-labels",
        default=None,
        help=(
            "Output labels path to use with medeas --labels. "
            "Default: <out without .vcf/.vcf.gz>.labels.lab"
        ),
    )
    parser.add_argument(
        "--chrom",
        default="1",
        help="Chromosome value to write in VCF CHROM column (default: 1)",
    )
    parser.add_argument(
        "--ref",
        default="A",
        help="Reference allele string for all sites (default: A)",
    )
    parser.add_argument(
        "--alt",
        default="T",
        help="Alternative allele string for all sites (default: T)",
    )

    args = parser.parse_args()

    try:
        labels_out_path = convert_medeas_to_vcf(
            snps_file=args.snps,
            labels_file=args.labels,
            out_file=args.out,
            chrom=args.chrom,
            ref=args.ref,
            alt=args.alt,
            out_labels=args.out_labels,
        )
        print(f"Wrote VCF: {args.out}")
        print(f"Wrote labels for medeas: {labels_out_path}")
        print(
            f"Run: medeas --vcf {args.out} --labels {labels_out_path} --out_dir <output_dir>"
        )
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
