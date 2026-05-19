#!/usr/bin/env python3
"""Convert PLINK BED/BIM/FAM to VCF.

Outputs:
- <out>.vcf when --out ends with .vcf
- <out>.vcf.gz when --out ends with .vcf.gz

Internally PLINK writes an intermediate .vcf that is moved/compressed to the
requested destination.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile


def _read_labels(labels_file):
    with open(labels_file) as f:
        labels = [line.strip() for line in f if line.strip()]
    if not labels:
        raise ValueError("Labels file is empty")
    return labels


def _rewrite_vcf_header_samples(vcf_path, labels):
    tmp_path = vcf_path + ".tmp"
    replaced = False

    with open(vcf_path, "r") as fin, open(tmp_path, "w") as fout:
        for line in fin:
            if line.startswith("#CHROM\t"):
                parts = line.rstrip("\n").split("\t")
                sample_names = parts[9:]
                if len(labels) != len(sample_names):
                    raise ValueError(
                        f"Labels file has {len(labels)} entries but VCF has {len(sample_names)} samples"
                    )
                fout.write("\t".join(parts[:9] + labels) + "\n")
                replaced = True
            else:
                fout.write(line)

    if not replaced:
        os.remove(tmp_path)
        raise RuntimeError("VCF header line (#CHROM) not found")

    os.replace(tmp_path, vcf_path)


def _bgzip_copy(src_path, dst_path):
    if shutil.which("bgzip") is None:
        raise RuntimeError(
            "bgzip not found on PATH. Install htslib (e.g. 'conda install -c bioconda htslib') "
            "to produce bgzip-compressed VCF files compatible with pysam/tabix."
        )

    with open(src_path, "rb") as fin, open(dst_path, "wb") as fout:
        proc = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=fout)
        try:
            shutil.copyfileobj(fin, proc.stdin)
        finally:
            proc.stdin.close()
            proc.wait()
        if proc.returncode != 0:
            raise RuntimeError(f"bgzip compression failed (exit {proc.returncode})")


def convert_plink_to_vcf(
    bfile_prefix,
    out_vcf,
    labels_file=None,
    plink_path="plink",
    allow_extra_chr=True,
    threads=1,
    iid_only=True,
):
    bed = bfile_prefix + ".bed"
    bim = bfile_prefix + ".bim"
    fam = bfile_prefix + ".fam"

    for path in (bed, bim, fam):
        if not os.path.isfile(path):
            raise FileNotFoundError(f"PLINK file not found: {path}")

    if labels_file is not None and not os.path.isfile(labels_file):
        raise FileNotFoundError(f"Labels file not found: {labels_file}")

    if shutil.which(plink_path) is None:
        raise FileNotFoundError(f"PLINK executable not found: {plink_path}")

    if threads < 0:
        raise ValueError("threads must be >= 0 (0 means all cores)")
    if threads == 0:
        threads = os.cpu_count() or 1

    if not (out_vcf.endswith(".vcf") or out_vcf.endswith(".vcf.gz")):
        raise ValueError("--out must end with .vcf or .vcf.gz")

    out_dir = os.path.dirname(out_vcf)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="plink_to_vcf_") as tmpdir:
        tmp_prefix = os.path.join(tmpdir, "plink_export")
        tmp_vcf = tmp_prefix + ".vcf"

        cmd_vcf = [
            plink_path,
            "--bfile", bfile_prefix,
            "--recode", "vcf",
            "--threads", str(threads),
            "--out", tmp_prefix,
        ]
        if allow_extra_chr:
            cmd_vcf.append("--allow-extra-chr")

        # PLINK1.9 uses --recode vcf-iid for IID-only sample names.
        # Use it when requested and available by trying it first.
        if iid_only:
            cmd_vcf_iid = [
                plink_path,
                "--bfile", bfile_prefix,
                "--recode", "vcf-iid",
                "--threads", str(threads),
                "--out", tmp_prefix,
            ]
            if allow_extra_chr:
                cmd_vcf_iid.append("--allow-extra-chr")

            result = subprocess.run(cmd_vcf_iid, capture_output=True, text=True)
            if result.returncode != 0:
                # Fallback for PLINK variants lacking vcf-iid token.
                result = subprocess.run(cmd_vcf, capture_output=True, text=True)
        else:
            result = subprocess.run(cmd_vcf, capture_output=True, text=True)

        if result.returncode != 0:
            sys.stderr.write(result.stderr)
            raise RuntimeError(f"PLINK PLINK->VCF conversion failed (exit {result.returncode})")

        if not os.path.isfile(tmp_vcf):
            raise RuntimeError("Expected VCF output was not created by PLINK")

        if labels_file is not None:
            _rewrite_vcf_header_samples(tmp_vcf, _read_labels(labels_file))

        if out_vcf.endswith(".vcf"):
            shutil.move(tmp_vcf, out_vcf)
        else:
            _bgzip_copy(tmp_vcf, out_vcf)



def main():
    parser = argparse.ArgumentParser(
        description="Convert PLINK BED/BIM/FAM files to VCF (.vcf or .vcf.gz)"
    )
    parser.add_argument("--bfile", required=True,
                        help="PLINK file prefix (without .bed/.bim/.fam)")
    parser.add_argument("--out", required=True,
                        help="Output VCF path (.vcf or .vcf.gz)")
    parser.add_argument("--labels", default=None,
                        help="Optional labels file with one label per line; when provided, these labels replace the VCF sample names")
    parser.add_argument("--plink-path", default="plink",
                        help="Path to the PLINK executable (default: plink)")
    parser.add_argument("--no-allow-extra-chr", dest="allow_extra_chr", action="store_false", default=True,
                        help="Do not pass --allow-extra-chr to PLINK")
    parser.add_argument("--fid-iid", dest="iid_only", action="store_false", default=True,
                        help="Keep default PLINK sample naming (FID+IID) instead of IID-only")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Worker threads for PLINK (0 = all cores)")
    args = parser.parse_args()

    try:
        convert_plink_to_vcf(
            bfile_prefix=args.bfile,
            out_vcf=args.out,
            labels_file=args.labels,
            plink_path=args.plink_path,
            allow_extra_chr=args.allow_extra_chr,
            threads=args.threads,
            iid_only=args.iid_only,
        )
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
