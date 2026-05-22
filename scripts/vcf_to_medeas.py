#!/usr/bin/env python3
"""Convert a VCF/BCF file to the medeas pseudo-haploid SNP matrix format.

Each diploid sample is reduced to one pseudo-haploid column:
  --vcf         always pseudo-haploidize by random allele selection
  --vcf_phased1 phased → haplotype 1, unphased → random
  --vcf_phased2 phased → haplotype 2, unphased → random

Output: one SNP per row, space-separated integers (0=missing, 1=ref, 2=alt).
"""

import argparse
import sys
import numpy as np


def _convert_vcf(vcf_file: str, mode: str, out_file: str) -> None:
    try:
        import pysam
    except ImportError as exc:
        raise ImportError(
            "pysam is required. Install it with 'pip install pysam'."
        ) from exc

    rng = np.random.default_rng()
    vf = pysam.VariantFile(vcf_file)
    n_samples = len(vf.header.samples)
    if n_samples == 0:
        vf.close()
        sys.exit(f"Error: no sample columns found in {vcf_file}")

    print(f"Converting {vcf_file} ({n_samples} samples, mode={mode}) -> {out_file}")

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

    n_written = 0
    with open(out_file, "w") as fout:
        for rec in vf:
            row = np.zeros(n_samples, dtype=np.int8)
            rp = [None]
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

            fout.write(" ".join(map(str, row)) + "\n")
            n_written += 1
            if n_written % 10000 == 0:
                print(f"  {n_written} SNPs written...")

    vf.close()
    print(f"Done: {n_written} SNPs written.")


def main():
    parser = argparse.ArgumentParser(
        description="Convert a VCF/BCF file to pseudo-haploid medeas SNP matrix."
    )

    input_group = parser.add_argument_group("Input source (required, choose exactly one)")
    input_source = input_group.add_mutually_exclusive_group(required=True)

    input_source.add_argument("--vcf", metavar="FILE",
                              help="VCF/BCF input. Random pseudo-haploidize.",
                              type=str, default=None)
    input_source.add_argument("--vcf-phased1", metavar="FILE",
                              help="VCF/BCF input. If phased, use haplotype 1; else random.",
                              type=str, default=None)
    input_source.add_argument("--vcf-phased2", metavar="FILE",
                              help="VCF/BCF input. If phased, use haplotype 2; else random.",
                              type=str, default=None)

    parser.add_argument("--out", required=True, metavar="FILE",
                        help="Output SNP matrix file (medeas format).")

    args = parser.parse_args()

    if args.vcf is not None:
        vcf_file = args.vcf
        mode = "vcf_random"
    elif args.vcf_phased1 is not None:
        vcf_file = args.vcf_phased1
        mode = "vcf_phased1"
    else:
        vcf_file = args.vcf_phased2
        mode = "vcf_phased2"

    _convert_vcf(vcf_file, mode, args.out)


if __name__ == "__main__":
    main()
