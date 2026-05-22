#!/usr/bin/env python3
"""Convert a VCF/BCF file directly to PLINK BED/BIM/FAM.

Outputs only:
- <prefix>.bed
- <prefix>.bim
- <prefix>.fam
"""

import argparse
import os
import shutil
import subprocess
import sys


def apply_labels_to_fam(fam_path, labels_file):
    ## read labels
    with open(labels_file) as f:
        labels = [line.strip() for line in f if line.strip()]

    if not labels:
        raise ValueError("Labels file is empty")

    ## read FAM file
    with open(fam_path) as f:
        fam_rows = [line.strip().split() for line in f if line.strip()]

    if len(labels) != len(fam_rows):
        raise ValueError(
            f"Labels file has {len(labels)} entries but FAM has {len(fam_rows)} samples"
        )

    ## rewrite fam file
    with open(fam_path, "w") as f:
        for label, row in zip(labels, fam_rows):
            if len(row) < 6:
                raise ValueError("Unexpected FAM format: expected at least 6 columns")
            row[0] = label
            f.write(" ".join(row) + "\n")


def convert_vcf_to_plink(vcf_file, out_prefix, plink_path, allow_extra_chr, threads, labels_file=None):
    if not os.path.isfile(vcf_file):
        raise FileNotFoundError(f"VCF file not found: {vcf_file}")

    if labels_file is not None and not os.path.isfile(labels_file):
        raise FileNotFoundError(f"Labels file not found: {labels_file}")

    if shutil.which(plink_path) is None:
        raise FileNotFoundError(f"PLINK executable not found: {plink_path}")

    if threads < 0:
        raise ValueError("threads must be >= 0 (0 means all cores)")
    if threads == 0:
        threads = os.cpu_count()

    out_dir = os.path.dirname(out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    cmd = [
        plink_path,
        "--vcf", vcf_file,
        "--make-bed",
        "--threads", str(threads),
        "--out", out_prefix,
    ]
    if allow_extra_chr:
        cmd.append("--allow-extra-chr")

    print("Running PLINK: " + " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.stderr.write(result.stderr)
        raise RuntimeError(f"PLINK VCF->BED/BIM/FAM conversion failed (exit {result.returncode})")

    if labels_file is not None:
        fam_path = out_prefix + ".fam"
        print("Applying labels to FAM: " + labels_file)
        apply_labels_to_fam(fam_path, labels_file)

    print("VCF to PLINK conversion complete.")


def main():
    parser = argparse.ArgumentParser(
        description="Convert a VCF/BCF file to PLINK BED/BIM/FAM files (2N -> 2N)"
    )
    parser.add_argument("--vcf", required=True, metavar="FILE",
                        help="Input VCF or BCF file (plain or gzip-compressed)")
    parser.add_argument("--out", required=True, metavar="PREFIX",
                        help="Output PLINK prefix (writes .bed/.bim/.fam)")
    parser.add_argument("--labels", default=None, metavar="FILE",
                        help="Optional labels file with one label per line; used to set FID in output .fam")
    parser.add_argument("--plink-path", default="plink", metavar="EXEC",
                        help="Path to the PLINK executable (default: plink)")
    parser.add_argument("--no-allow-extra-chr", dest="allow_extra_chr", action="store_false", default=True,
                        help="Do not pass --allow-extra-chr to PLINK")
    parser.add_argument("-t", "--threads", type=int, default=1, metavar="N",
                        help="Worker threads for PLINK (0 = all cores)")
    args = parser.parse_args()

    try:
        convert_vcf_to_plink(
            vcf_file=args.vcf,
            out_prefix=args.out,
            plink_path=args.plink_path,
            allow_extra_chr=args.allow_extra_chr,
            threads=args.threads,
            labels_file=args.labels,
        )
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
