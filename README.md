# MEDEAS

MEDEAS is a Python tool for inferring coalescence times within and between populations from genotype data. It computes a multi-dimensional scaling (MDS) projection and infers coalescence times from the eigenvalues of the MDS matrix.

The goal is achieved by computing a multi-dimensional scaling (MDS) projection in a K-dimensional space and then inferring coalescence times from the eigenvalues of the MDS. The population tree is inferred using the center of mass of the various clusters in the MDS projection.

Split-time estimation requires population labels. Without labels (or with `--no-split`), MEDEAS still produces MDS/PCA plots, and coordinate files.

MEDEAS works with haploid individuals. If diploid individuals are passed (using plink or vcf format) MEDEAS converts the diploid individuals to haploid individuals. MEDEAS requires two inputs to compute split times:
### Genotypes matrix format
MEDEAS accepts three genotype matrix formats:

- **MEDEAS format** (default): A plain-text file where each column (separated by a space) is one haploid individual and each row is one genomic position. No header or row names. Reference allele is encoded as 1, alternative allele as 2, and missing data as 0.
- **PLINK format** (`.bed/.bim/.fam`): Diploid unphased genotypes are automatically pseudo-haploidized by randomly selecting one allele per heterozygous site. If no label file is provided, population labels are taken from the first column (FID) of the `.fam` file.
- **VCF format** (`.vcf`, `.vcf.gz`, `.bcf`): Diploid genotypes are automatically pseudo-haploidized. With `--vcf`, alleles at heterozygous sites are chosen randomly. With `--vcf-phased1` or `--vcf-phased2`, the first or second haplotype is used for phased variants; unphased variants are randomized.

### Label file format

The label file should contain the population of each haploid individual, one by line, without any other information. The number of lines in this file should be the same as the number of columns in the genotype matrix file.

A plink file may also contain population information in the `.fam` file (FID; first column). If a plink file is used as input without a specific label file, then the labels from the `.fam` file are used.

---

## Dependencies

- Python ≥ 3.6
- numpy
- matplotlib
- scipy
- scikit-bio
- pysam *(required for VCF/BCF input)*
- plink 1.9 *(required for `--bfile` and `--use-plink`)*

---

## Installation

```bash
git clone git@github.com:sapfo/medeas.git
cd medeas
git checkout sam-dev
pip install .
```

---

## Usage

```bash
medeas --snps <genotypes_file> [--labels <labels_file>] --out-dir <output_folder>
```

### Input formats (choose exactly one)

| Flag | Description |
|------|-------------|
| `--snps FILE` | MEDEAS text format: one SNP per row, space-separated integers (0=missing, 1=ref, 2=alt), one column per haploid individual |
| `--bfile PREFIX` | PLINK binary prefix (`.bed/.bim/.fam`). Diploid individuals are pseudo-haploidized. |
| `--vcf FILE` | VCF/BCF (plain or gzipped). Each diploid sample is randomly pseudo-haploidized. |
| `--vcf-phased1 FILE` | VCF/BCF. Phased samples use haplotype 1; unphased are randomly pseudo-haploidized. |
| `--vcf-phased2 FILE` | VCF/BCF. Phased samples use haplotype 2; unphased are randomly pseudo-haploidized. |

### Key options

| Flag | Default | Description |
|------|---------|-------------|
| `--labels FILE` | — | Population label file: one label per individual, one per line, matching column order of the genotype file. Required for split-time estimation. |
| `--out-dir DIR` | '.' | Output folder (created if absent). |
| `--no-split` | off | Skip tree and split-time estimation; produce only distance, MDS, and PCA outputs. Implied when no labels are provided. |
| `-bws N` | 100 | SNPs per bootstrap window. Recommendation: ~1% of total SNPs. |
| `-bsn N` | 100 | Number of bootstrap replicates (0 = single run with all SNPs). |
| `--topology NEWICK` | — | Fix the tree topology (Newick, population indices in label order). |
| `--skip-calculate-matrix` | off | Reuse previously computed distance matrices (skips genotype file reading). |
| `--threads N` | 0 | Worker threads (0 = all available cores). |
| `--output-level N` | 0 | Verbosity: 0 minimal, 1 conventional, 2 maximum. |
| `--n-dims N` | all | Number of MDS/PCA dimensions to compute and plot. |
| `--plot-dims "p,q ..."` | auto | Pairs of dimensions to plot, 1-based (e.g. `"1,2 3,4"`). Default: consecutive pairs. |
| `--no-mds` | off | Skip MDS scatter plots. |
| `--no-pca` | off | Skip PCA scatter plots. |
| `--use-plink` | off | Use PLINK routines for MDS/PCA (requires `--bfile` and `--no-split`). |
| `--plink-path FILE` | 'plink' | Path to PLINK 1.9 executable. |
| `--detailed-output` | off | Generate additional diagnostic plots (distance matrix, eigenvalue histogram, SFS). |

A full list of parameters is available via `medeas --help`.

### Example

Running with split-time estimation:

```bash
medeas --snps test/snp.dat --labels test/pop_label.dat --out-dir results/
```

Running to obtain MDS/PCA only:

```bash
medeas --snps test/snp.dat --labels test/pop_label.dat --out-dir results/  --no-split
```

---

## Output files

| File / folder | Description |
|---------------|-------------|
| `mds_plot/` | MDS scatter plots |
| `pca_plot/` | PCA scatter plots |
| `MDS_coordinate.txt` | MDS coordinates (one row per individual) |
| `PCA_coordinate.txt` | PCA coordinates (one row per individual) |
| `eigenvalues.pdf` | Eigenvalue scree plot |
| `SFS.pdf` / `SFS.txt` | Folded site-frequency spectrum |
| `tree.txt` | Inferred population tree *(requires labels)* |
| `split_time.txt` | Split times per node, all bootstrap replicates *(requires labels)* |
| `effective_size.txt` | Effective sizes after split, all bootstrap replicates *(requires labels)* |
| `between_population_coalescence_time.txt` | Between-population coalescence times *(requires labels)* |
| `within_population_coalescence_time.txt` | Within-population coalescence times *(requires labels)* |
| `split_bootstraped_confidence_interval.txt` | 95% CI for split times *(requires labels)* |
| `population_bootstraped_confidence_interval.txt` | 95% CI for within-pop coalescence and effective sizes *(requires labels)* |
| `asd_matrices/` | Cached distance matrices (reused with `--skip-calculate-matrix`) |
| `MDS_eigensystem/` | Cached MDS eigensystems |

---

## Conversion scripts

MEDEAS accepts different input formats. The `scripts/` folder contains conversion scripts to convert between the different formats (e.g., useful for testing purposes). A full list of available parameters for each converter can be obtained by running the script with `--help`.

### `vcf_to_medeas.py` — VCF → MEDEAS format

Converts a VCF/BCF directly to the MEDEAS SNP matrix. Supports the same three pseudo-haploidization modes as the main tool.

```bash
vcf_to_medeas --vcf input.vcf.gz --out snp.dat
vcf_to_medeas --vcf-phased1 input.vcf.gz --out snp.dat
vcf_to_medeas --vcf-phased2 input.vcf.gz --out snp.dat
```

### `plink_to_medeas.py` — PLINK → MEDEAS format

Converts a PLINK to MEDEAS SNP format. Pseudo-haplodization is random. Labels may be overtaken from the .fam file (FID; first column).

```bash
plink_to_medeas --bfile prefix --out snp.dat --out-labels labels.dat
```

### `vcf_to_plink.py` — VCF → MEDEAS format

Converts a VCF/BCF to PLINK.

```bash
vcf_to_plink --vcf input.vcf.gz --out prefix
```

### `medeas_to_plink.py` — MEDEAS format → PLINK

Converts a MEDEAS format to a VCF/BCF.

```bash
medeas_to_plink --snps snp.dat --labels labels.dat --out prefix
```

### `plink_to_vcf.py` — PLINK → VCF

Converts a PLINK to a VCF/BCF.

```bash
plink_to_vcf --bfile prefix --out out.vcf.gz
```

### `medeas_to_vcf.py` — MEDEAS format → VCF

Converts a MEDEAS format to a VCF/BCF.

```bash
medeas_to_vcf --snps snp.dat --out out.vcf.gz
```


---

## Publication

Pending.
