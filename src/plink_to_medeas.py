import os
import numpy as np


def convert_plink_to_medeas(bfile_prefix, snp_file, labels_file):
    """Convert PLINK binary files (.bed/.bim/.fam) to the medeas haploid format.

    Mirrors the shell pipeline in launch.sh:
      1. PLINK --recode12 --transpose  -> alleles encoded as 1/2, missing as 0.
      2. Each diploid individual is split into two haploid pseudoindividuals
         (the two alleles), with heterozygous sites randomly phased.
      3. Population labels (FAM family-ID) are written once per haploid
         pseudoindividual (each label appears twice).

    Output format (snp_file):
      - One SNP per row, space-separated integers.
      - Column order: hap1_ind1 hap2_ind1 hap1_ind2 hap2_ind2 ...
      - Values: 0=missing, 1=ref allele (A1), 2=alt allele (A2).
    """
    fam_file = bfile_prefix + ".fam"
    bim_file = bfile_prefix + ".bim"
    bed_file = bfile_prefix + ".bed"

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

    print(f"Writing: {snp_file}")
    with open(bed_file, "rb") as bed, open(snp_file, "w") as out:
        magic = bed.read(3)
        if magic != b"\x6c\x1b\x01":
            raise ValueError(
                f"Not a valid PLINK BED file (unexpected magic bytes): {bed_file}"
            )

        snp_done = 0
        while snp_done < n_snps:
            current = min(chunk_snps, n_snps - snp_done)
            raw = np.frombuffer(bed.read(current * n_bytes_per_snp), dtype=np.uint8)
            if raw.size != current * n_bytes_per_snp:
                raise ValueError("Unexpected end of BED file while reading SNP data")
            raw = raw.reshape(current, n_bytes_per_snp)

            bits = np.unpackbits(raw, axis=1, bitorder="little")
            bits = bits[:, : 2 * n_samples].reshape(current, n_samples, 2)
            codes = bits[:, :, 0] + 2 * bits[:, :, 1]

            hap1 = np.zeros((current, n_samples), dtype=np.int8)
            hap2 = np.zeros((current, n_samples), dtype=np.int8)

            mask_hom_ref = codes == 0
            hap1[mask_hom_ref] = 1
            hap2[mask_hom_ref] = 1

            mask_hom_alt = codes == 3
            hap1[mask_hom_alt] = 2
            hap2[mask_hom_alt] = 2

            mask_het = codes == 2
            flip = rng.random((current, n_samples)) > 0.5
            hap1[mask_het & ~flip] = 1
            hap2[mask_het & ~flip] = 2
            hap1[mask_het & flip] = 2
            hap2[mask_het & flip] = 1

            geno_haploid = np.empty((current, 2 * n_samples), dtype=np.int8)
            geno_haploid[:, 0::2] = hap1
            geno_haploid[:, 1::2] = hap2

            np.savetxt(out, geno_haploid, fmt="%d", delimiter=" ")
            snp_done += current
            if snp_done % 20000 == 0 or snp_done == n_snps:
                print(f"Converted {snp_done}/{n_snps} SNPs")

    print(f"Writing: {labels_file}")
    with open(labels_file, "w") as f:
        for fid in family_ids:
            f.write(fid + "\n")
            f.write(fid + "\n")

    print("PLINK conversion complete.")
