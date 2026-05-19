#!/usr/bin/env python3
"""Convert a VCF/BCF file to the medeas pseudo-haploid format.

Workflow:
  1. PLINK converts VCF -> BED (in a temporary directory).
  2. Each diploid individual is converted to one pseudo-haploid column;
     heterozygous sites are randomly sampled as 1 or 2.
  3. Population labels are written once per individual.

Population labels are taken from the FAM file produced by PLINK (= VCF sample
IDs) unless --labels is supplied, in which case that file is used verbatim.
"""

import argparse
import os
import subprocess
import sys
import tempfile

# Allow ``from plink_to_medeas import ...`` when running as a script.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plink_to_medeas import convert_plink_to_medeas  # noqa: E402


def convert_vcf_to_medeas(
    vcf_file,
    snp_file,
    labels_file,
    labels_input=None,
    plink_path="plink",
    allow_extra_chr=True,
    threads=1,
):
    """Convert vcf_file to the medeas pseudo-haploid format.

    if not os.path.isfile(vcf_file):
        raise FileNotFoundError(f"VCF file not found: {vcf_file}")

    if labels_input is not None and not os.path.isfile(labels_input):
        raise FileNotFoundError(f"Labels file not found: {labels_input}")

    with tempfile.TemporaryDirectory() as tmpdir:
        bfile = os.path.join(tmpdir, "tmp")

        # ── Step 1: VCF → BED ──────────────────────────────────────────────
        cmd = [plink_path, "--vcf", vcf_file, "--make-bed", "--out", bfile]
        if allow_extra_chr:
            cmd.append("--allow-extra-chr")
        cmd += ["--threads", str(max(1, threads if threads != 0 else (os.cpu_count() or 1)))]

        print("Running PLINK: " + " ".join(cmd))
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            sys.stderr.write(result.stderr)
            raise RuntimeError("PLINK VCF→BED conversion failed (see output above)")

        # ── Step 2: BED → medeas (SNP matrix + labels from FAM) ───────────
        convert_plink_to_medeas(bfile, snp_file, labels_file, threads=threads)

        # ── Step 3 (optional): override labels with user-supplied file ─────
        if labels_input is not None:
            with open(labels_input) as f:
                diploid_labels = [ln.strip() for ln in f if ln.strip()]

            fam_path = bfile + ".fam"
            with open(fam_path) as f:
                n_samples = sum(1 for ln in f if ln.strip())

            if len(diploid_labels) != n_samples:
                raise ValueError(
                    f"--labels file has {len(diploid_labels)} entries but VCF "
                    f"contains {n_samples} samples"
                )

            with open(labels_file, "w") as f:
                for lbl in diploid_labels:
                    f.write(lbl + "\n")

    print("VCF conversion complete.")


def main():
    parser = argparse.ArgumentParser(
        description=("Convert a VCF/BCF file to pseudo-haploid medeas SNP/labels files. ")
    )
    parser.add_argument("--vcf", required=True,
                        help="Input VCF or BCF file (plain or gzip-compressed).")
    parser.add_argument("--snp-out", required=True,
                        help="Output SNP matrix file (medeas format).")
    parser.add_argument("--labels", metavar="FILE", default=None,
                        help="Optional input labels file with one population label per diploid individual (same order as VCF samples).  When omitted, VCF sample IDs are used.")
    parser.add_argument("--labels-out", required=True, metavar="FILE",
                        help="Output labels file (one label per pseudo-haploid individual).")
    parser.add_argument("--plink-path", default="plink", metavar="PATH",
                        help="Path to the PLINK 1.9 executable (default: plink).")
    parser.add_argument("--no-allow-extra-chr", dest="allow_extra_chr", action="store_false", default=True,
                        help="Do not pass --allow-extra-chr to PLINK.")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Worker threads for SNP conversion (0 = all cores, default: 1).")

    args = parser.parse_args()

    convert_vcf_to_medeas(
        vcf_file=args.vcf,
        snp_file=args.snp_out,
        labels_file=args.labels_out,
        labels_input=args.labels,
        plink_path=args.plink_path,
        allow_extra_chr=args.allow_extra_chr,
        threads=args.threads,
    )


if __name__ == "__main__":
    main()
