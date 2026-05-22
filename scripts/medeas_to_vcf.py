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
import sys


GENOTYPE_MAP = {
    "0": "./.",
    "1": "0/0",
    "2": "1/1",
}



@contextlib.contextmanager
def open_text_output(path):
    """Open a text file for writing, using pysam.BGZFile for .gz paths."""
    if path.endswith(".gz"):
        import pysam
        with pysam.BGZFile(path, "wb") as raw:
            text_wrapper = io.TextIOWrapper(raw, encoding="utf-8")
            try:
                yield text_wrapper
            finally:
                text_wrapper.flush()
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
            f"Row {snp_index} has {len(genotypes)} columns but expected {expected_haploids}"
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
        str(snp_index),
        f"snp{snp_index}",
        ref,
        alt,
        ".",
        "PASS",
        ".",
        "GT",
    ]
    return "\t".join(fixed_fields + calls) + "\n"


def convert_medeas_to_vcf(
    snps_file,
    out_file,
    chrom="1",
    ref="A",
    alt="T",
):
    if ref == alt:
        raise ValueError("REF and ALT alleles must be different")

    with open(snps_file) as fin:
        for first_line in fin:
            if first_line.strip():
                expected_haploids = len(first_line.strip().split())
                break
        else:
            raise ValueError("SNP file is empty")

    sample_ids = [f"ind_{i}" for i in range(1, expected_haploids + 1)]

    out_dir = os.path.dirname(out_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

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


def main():
    parser = argparse.ArgumentParser(
        description="Convert MEDEAS SNP file to a homozygous diploid VCF file (N -> N)"
    )
    parser.add_argument(
        "--snps",
        required=True, metavar="FILE",
        help="SNP file in MEDEAS format (one SNP per row, space-separated integers)",
    )
    parser.add_argument(
        "--out",
        required=True, metavar="FILE",
        help="Output VCF file (.vcf or .vcf.gz)",
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
        convert_medeas_to_vcf(
            snps_file=args.snps,
            out_file=args.out,
            chrom=args.chrom,
            ref=args.ref,
            alt=args.alt,
        )
        print(f"Wrote VCF: {args.out}")
        print(f"Run: medeas --vcf {args.out} --out_dir <output_dir>")
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
