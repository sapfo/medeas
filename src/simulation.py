import argparse
import copy
import os
import tempfile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.ticker import MaxNLocator
import subprocess
import datetime
import sys
from multiprocessing import cpu_count


from skbio.tree import TreeNode
import pickle
from src.clustering import get_mds_coordinate
from src.clustering import build_split_index_matrix
from src.make_asd import compute_asd_matrix


def _mds_plink(bfile: str, k: int, dist_path: str, plink_path: str,
               plink_params: str, threads: int, tmpdir: str):
    """Run PLINK MDS from a precomputed IBS distance matrix.

    Returns (coords, eigenvalues, samples)."""
    prefix = os.path.join(tmpdir, "plink_mds")

    mds_cmd = [
        plink_path,
        "--bfile", bfile,
        "--read-dists", dist_path,
        "--allow-no-sex",
        "--cluster",
        "--mds-plot", str(k), "eigvals",
        "--threads", str(threads),
        "--out", prefix,
    ]
    if plink_params:
        mds_cmd += plink_params.split()

    print(f"Running: {' '.join(mds_cmd)}", file=sys.stderr)
    r = subprocess.run(mds_cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write(r.stderr)
        raise RuntimeError(f"PLINK MDS step failed (exit {r.returncode})")

    mds_df = pd.read_csv(f"{prefix}.mds", sep=r"\s+")
    samples = mds_df[["FID", "IID"]].copy()
    coord_cols = [c for c in mds_df.columns if c.startswith("C") and c[1:].isdigit()]
    coords = mds_df[coord_cols].values.astype(np.float64)

    eigvals_file = f"{prefix}.mds.eigvals"
    if os.path.isfile(eigvals_file):
        eigenvalues = np.loadtxt(eigvals_file)
    else:
        eigenvalues = np.sum(coords ** 2, axis=0)

    return coords, eigenvalues, samples


def _compute_pca_plink(bfile: str, k: int, plink_path: str,
                       threads: int, tmpdir: str):
    """Compute PCA from BED file using PLINK.

    Returns (coords, eigenvalues)."""
    prefix = os.path.join(tmpdir, "plink_pca")
    pca_cmd = [
        plink_path,
        "--bfile", bfile,
        "--allow-no-sex",
        "--pca", str(k),
        "--threads", str(threads),
        "--out", prefix,
    ]

    print(f"Running: {' '.join(pca_cmd)}", file=sys.stderr)
    r = subprocess.run(pca_cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write(r.stderr)
        raise RuntimeError(f"PLINK PCA step failed (exit {r.returncode})")

    eigenvec_df = pd.read_csv(f"{prefix}.eigenvec", sep=r"\s+", header=None)
    coords = eigenvec_df.iloc[:, 2:].values.astype(np.float64)

    eigenval_file = f"{prefix}.eigenval"
    if os.path.isfile(eigenval_file):
        eigenvalues = np.loadtxt(eigenval_file)
    else:
        eigenvalues = np.var(coords, axis=0)

    return coords, eigenvalues


def _plot_eigenvalues_array(eigenvalues: np.ndarray, file_path: str):
    """Plot a simple eigenvalue scatter and save to file_path."""
    lambdas = np.array(eigenvalues).ravel()
    lambdas = -np.sort(-lambdas)
    fig, ax = plt.subplots()
    ax.plot(lambdas, "o")
    ax.set_ylabel("Eigenvalues")
    ax.set_xlabel("Eigenvalues index")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    fig.savefig(file_path)
    plt.close(fig)


def _get_variant_samples(variant_path: str):
    """Extract sample IDs from VCF/BCF (gzipped or plain) using pysam."""
    try:
        import pysam
    except ImportError as exc:
        raise ImportError(
            "pysam is required for --vcf/--vcf-phased1/--vcf-phased2 input. "
            "Install it with 'pip install pysam'."
        ) from exc

    vf = pysam.VariantFile(variant_path)
    return list(vf.header.samples)




class SimulationInfo(object):

    def __init__(self):
        parser = argparse.ArgumentParser()

        input_group = parser.add_argument_group(
            "Input source (required, choose exactly one)"
        )
        input_source = input_group.add_mutually_exclusive_group(required=True)

        input_source.add_argument("--snps", metavar="FILE",
                                  help="SNP file in MEDEAS format (one SNP per row, space-separated integers, one column per individual).")

        input_source.add_argument("--bfile", metavar="PREFIX",
                                  help="PLINK binary file prefix (.bed/.bim/.fam). "
                                        "Conversion to pseudo-haploid individuals.",
                                  type=str, default=None)

        input_source.add_argument("--vcf", metavar="FILE",
                                  help="VCF/BCF input (gzipped or plain). "
                                       "Conversion to pseudo-haploid individuals.",
                                  type=str, default=None)

        input_source.add_argument("--vcf-phased1", metavar="FILE",
                                  help="VCF/BCF input (gzipped or plain). "
                                       "If genotype is phased, use first haplotype; "
                                       "if unphased, pseudo-haploidize.",
                                  type=str, default=None)

        input_source.add_argument("--vcf-phased2", metavar="FILE",
                                  help="VCF/BCF input (gzipped or plain). "
                                       "If genotype is phased, use second haplotype; "
                                       "if unphased, pseudo-haploidize.",
                                  type=str, default=None)

        options_group = parser.add_argument_group("options")

        options_group.add_argument("--labels", metavar="FILE", 
                                    help="File containing the labels")

        options_group.add_argument("--out-dir", metavar="DIR", default=".",
                                    help="Folder where results and temporal data should be stored (default: current working directory)")

        options_group.add_argument("-bws", "--boot-window-size", metavar="N",
                                   help="How many markers do we have in each bootstraping windows",
                                   type=int, default=100)

        options_group.add_argument("-bsn", "--bootstrap-number", metavar="N",
                                   help="How many bootstrap do we perform",
                                   type=int, default=100)

        options_group.add_argument("--topology", metavar="NEWICK",
                                   help="What is the topology of the population (newick format, following label order)",
                                   type=str, default=None)

        options_group.add_argument("--skip-calculate-matrix",
                                   help="Skip the computation of the distance matrices and the related MDS matrix",
                                   action="store_true")

        options_group.add_argument("--output-level", metavar="LEVEL",
                                   help="How many information should be printed & saved: 0 -minimal, 1 - conventional, 2 - most of it",
                                   type=int, default=0)

        options_group.add_argument("-t", "--threads", metavar="N",
                                   help="Number of parallel process to be launch. 0 (default) used all available cores",
                                   type=int, default=0)

        options_group.add_argument("--no-split",
                                   help="Skip tree and split-time estimation; compute only distance, MDS, and PCA outputs",
                                   action="store_true")

        options_group.add_argument("--detailed-output",
                                   help="Return detailed output",
                                   action="store_true")

        options_group.add_argument("--use-plink",
                                   help="Plot MDS and PCA using PLINK routines. "
                                        "Requires --bfile.",
                                   action="store_true")

        options_group.add_argument("--plink-path", metavar="FILE",
                                   help="Path to the PLINK executable (default: 'plink')",
                                   type=str, default="plink")

        options_group.add_argument("--n-dims", metavar="N",
                                   help="Number of eigenvalues to compute and plot (default: all, i.e. n_samples-1)",
                                   type=int, default=0)

        options_group.add_argument("--plot-dims", metavar="DIMS",
                                   help="Pairs of dimensions to plot, 1-based, space-separated (e.g. '1,2 3,4'). "
                                        "Default: consecutive pairs (1,2), (3,4), ...",
                                   type=str, default=None)

        options_group.add_argument("--no-mds",
                                   help="Do not generate MDS scatter plots.",
                                   action="store_true")

        options_group.add_argument("--no-pca",
                                   help="Do not generate PCA scatter plots.",
                                   action="store_true")

        args = parser.parse_args()


        ## create output folder
        self.output_folder = args.out_dir
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        ## get input file
        self.snp_file = None
        self.snp_file_type = None
        self.label_file = args.labels

        if args.bfile is not None:
            bfile_prefix = args.bfile
            if bfile_prefix.endswith(".bed") or bfile_prefix.endswith(".bim") or bfile_prefix.endswith(".fam"):
                bfile_prefix = bfile_prefix[:-4]

            bed_file = bfile_prefix + ".bed"
            bim_file = bfile_prefix + ".bim"
            fam_file = bfile_prefix + ".fam"
            for path in (bed_file, bim_file, fam_file):
                if not os.path.isfile(path):
                    sys.exit(f"Error: Missing PLINK file: {path}")

            self.snp_file = bed_file
            self.snp_file_type = "plink"
            self.n_samples = sum(1 for line in open(fam_file) if line.strip()) ## number of rows

            # Create temp label file by extracting the labels from FAM.
            if args.labels is None:             
                label_file = os.path.join(self.output_folder, "labels_from_fam.dat")
                with open(fam_file) as fin, open(label_file, "w") as fout:
                    for line in fin:
                        parts = line.strip().split()
                        if not parts:
                            continue
                        fid = parts[0]
                        fout.write(fid + "\n")
                self.label_file = label_file

        elif args.vcf is not None or args.vcf_phased1 is not None or args.vcf_phased2 is not None:
            if args.vcf is not None:
                self.snp_file = args.vcf
                self.snp_file_type = "vcf_random"
            elif args.vcf_phased1 is not None:
                self.snp_file = args.vcf_phased1
                self.snp_file_type = "vcf_phased1"
            elif args.vcf_phased2 is not None:
                self.snp_file = args.vcf_phased2
                self.snp_file_type = "vcf_phased2"

            try:
                self.n_samples = len(_get_variant_samples(self.snp_file))
            except Exception as exc:
                sys.exit(f"Error: unable to parse variant header: {exc}")

        else: ## medeas input
            self.snp_file = args.snps
            self.snp_file_type = "medeas"

            with open(args.snps) as f:
                first_line = f.readline()
            self.n_samples = len(first_line.split())    ## number of columns


        ## check that snp file exists and contains samples
        if not os.path.isfile(self.snp_file):
            sys.exit(f"Error: variant file not found: {self.snp_file}")

        if self.n_samples == 0:
            sys.exit("Error: SNP file contains no samples.")

        ## read the label file
        if self.label_file is None:
            if self.snp_file_type.split("_")[0] == "vcf":
                try:
                    self.labels = np.array(_get_variant_samples(self.snp_file))
                except Exception as exc:
                    sys.exit(f"Error: unable to parse variant header: {exc}")
            else:
                self.labels = np.array([f"ind{i+1}" for i in range(self.n_samples)])

        else:
            if not os.path.isfile(self.label_file):
                sys.exit("Error: The label file does not exist. Exiting Now.")
                    
            try:
                with open(self.label_file) as f:
                    lines = f.readlines()
            except:
                sys.exit(
                    "Error: A problem occurs when loading the label file. Please check that you use the right format.")

            labels = [l.rstrip() for l in lines]
            self.labels = np.array(labels)

        ## check that the number of samples in the snp and label file correspond
        if len(self.labels) != self.n_samples:
            sys.exit(
                f"Error: label file has {len(self.labels)} labels but input contains {self.n_samples} individuals."
            )

        _, index, numerical_labels = np.unique(self.labels,return_inverse=True, return_index = True)
        self.populations = self.labels[np.sort(index)]
        self.numerical_labels = np.sort(numerical_labels)
        self.K = len(self.populations)

        if self.K < 2 and not self.no_split:
            sys.exit("Error: At least 2 populations are required to estimate split-times. Provide a label file with at least 2 unique labels, or use --no-split to skip split-time estimation.")


        ## plink binary usage
        self.use_plink = args.use_plink
        self.plink_path = args.plink_path
        if self.use_plink:
            if args.bfile is None:
                sys.exit("Error: --use-plink is only supported with PLINK input (--bfile).")
            if args.no_split is False:
                sys.exit("Error: --use-plink requires --no-split, as plink cannot be used to estimate split-times.")

            import shutil
            if shutil.which(self.plink_path) is None:
                sys.exit(f"Error: PLINK executable not found: '{self.plink_path}'. "
                         "Install PLINK or pass the correct path with --plink-path.")


        ## dimensions to compute and plot
        self.n_dims = args.n_dims
        if self.n_dims < 0:
            sys.exit("Error: --n-dims must be >= 0 (0 means all dimensions).")
        elif self.n_dims == 0:
            self.n_dims = max(1, len(self.labels) - 1)


        # Parse --plot-dims into a list of 0-based (p, q) pairs.
        if args.plot_dims is not None:
            plot_dims = []
            for token in args.plot_dims.split():
                parts = token.split(",")
                if len(parts) != 2:
                    sys.exit(f"Error: --plot-dims token '{token}' must be two comma-separated integers (e.g. '1,2').")
                try:
                    p1, p2 = int(parts[0]) - 1, int(parts[1]) - 1
                except ValueError:
                    sys.exit(f"Error: --plot-dims token '{token}' must contain integers.")
                if p1 < 0 or p2 < 0:
                    sys.exit(f"Error: --plot-dims dimensions are 1-based; got '{token}'.")
                plot_dims.append((p1, p2))

            max_dim = max(max(p1, p2) for p1, p2 in plot_dims) + 1  # back to 1-based
            if max_dim > self.n_dims:
                sys.exit(f"Error: --plot-dims requests dimension {max_dim} but --n-dims {self.n_dims} is too small.")
            self.plot_dims = plot_dims
        else:
            self.plot_dims = None  # auto: consecutive pairs


        ## multi-threading
        self.NCORE = args.threads
        if self.NCORE < 0:
            sys.exit("Error: --threads must be >= 0 (0 means all available threads).")
        elif self.NCORE == 0:
            self.NCORE = cpu_count()

        self.bootsize = args.boot_window_size
        self.skip_calculate_matrix = args.skip_calculate_matrix
        self.bootstrap_number = args.bootstrap_number
        self.output_level = args.output_level
        self.topology = args.topology
        self.has_labels = args.labels is not None or args.bfile is not None
        self.no_split = args.no_split or not self.has_labels
        self.detailed_output = args.detailed_output
        self.no_mds = args.no_mds
        self.no_pca = args.no_pca
        self.logfile = os.path.join(self.output_folder, "simulation.log")
        self.generate_initial_output(args)

        asd_folder = "asd_matrices"
        mds_folder = "MDS_eigensystem"
        asd_full_path = os.path.join(self.output_folder,asd_folder)
        mds_full_path = os.path.join(self.output_folder,mds_folder)
        for path in [asd_full_path, mds_full_path]:
            if not os.path.exists(path):
                os.makedirs(path)
        if self.use_plink:
            self.asd_pattern = os.path.join(asd_full_path, 'p{}_plink.asd.data')
            self.vec_pattern = os.path.join(mds_full_path, 'p{}_plink.vecs.data')
        else:
            self.asd_pattern = os.path.join(asd_full_path, 'p{}.asd.data')
            self.vec_pattern = os.path.join(mds_full_path, 'p{}.vecs.data')


    def export_sfs_simple(self):
        self.sfs = np.array(self.sfs)
        max_freq = len(self.labels)
        plt.figure()
        plt.bar(self.sfs[0]/max_freq,self.sfs[1], width=0.5/max_freq)
        plt.xlabel("Mutation Frequency")
        plt.ylabel("Site count")
        plt.suptitle("Overall site frequency spectrum")

        filePath = os.path.join(self.output_folder, "SFS.pdf")
        plt.savefig(filePath)
        plt.close()
        with open(os.path.join(self.output_folder, "SFS.txt"), 'w') as f:
            np.savetxt(f, np.transpose(self.sfs).astype(int),fmt='%i')


    def export_sfs(self):
        self.sfs = np.array(self.sfs)
        max_freq = len(self.labels)
        plt.figure()
        if(len(self.sfs[0]) < 25):
            plt.bar(self.sfs[0]/max_freq,self.sfs[1], width=0.5/max_freq)
        else:
            plt.subplot(211)
            plt.bar(self.sfs[0,0:25]/max_freq, self.sfs[1,0:25],  width=0.5/max_freq)
            plt.ylabel("Site count")
            plt.subplot(212)
            plt.bar(self.sfs[0]/max_freq, self.sfs[1],  width=0.5/max_freq)
        plt.xlabel("Mutation Frequency")
        plt.ylabel("Site count")
        plt.suptitle("Overall site frequency spectrum")

        filePath = os.path.join(self.output_folder, "SFS.pdf")
        plt.savefig(filePath)
        plt.close()
        with open(os.path.join(self.output_folder, "SFS.txt"), 'w') as f:
            np.savetxt(f, np.transpose(self.sfs).astype(int),fmt='%i')


    def plot_eigenvalues_simple(self):
        with open(self.vec_pattern.format(2), 'rb') as f:
            lambdas, vecs = pickle.load(f)
        lambdas = -np.sort(-lambdas)
        fig, ax = plt.subplots()
        ax.plot(lambdas,"o")
        ax.set_ylabel("Eigenvalues")
        ax.set_xlabel("Eigenvalues index")
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        filePath = os.path.join(self.output_folder, "eigenvalues.pdf")
        fig.savefig(filePath)
        plt.close(fig)


    def plot_eigenvalues(self):
        with open(self.vec_pattern.format(2), 'rb') as f:
            lambdas, vecs = pickle.load(f)
        lambdas = -np.sort(-lambdas)
        fig, (ax1, ax2) = plt.subplots(2, 1)
        ax1.plot(lambdas,"o")
        ax1.set_ylabel("Eigenvalues")
        ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

        ax2.plot(lambdas[self.K:-2],"o")
        ax2.set_ylabel("Eigenvalues")
        ax2.set_xlabel("Eigenvalues index")
        ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
        filePath = os.path.join(self.output_folder, "eigenvalues.pdf")
        fig.savefig(filePath)
        plt.close(fig)

        plt.figure()
        plt.hist(-np.sort(-lambdas)[self.K:-2])
        plt.xlabel("Eigenvalue")
        plt.ylabel("Eigenvalues count")
        filePath = os.path.join(self.output_folder, "histogram_eigenvalues.pdf")
        plt.savefig(filePath)
        plt.close()


    def pops_contain_at_least_2_individual(self):
        _, counts = np.unique(self.labels,return_counts = True)
        return all(counts > 1)


    def plot_distance_matrix(self, delta):
        with open(self.labels_file) as f:
            lines = f.readlines()
        labels_individual = np.array([l.split()[0] for l in lines])
        if not len(delta) == len(labels_individual):
            sys.exit("Error: The number of individual in the label file is not the same as the number of individual\
in the distance matrix. Exiting Now.")
        label_pop = self.populations
        sorting_index = np.argsort(labels_individual)
        individual_per_pop = [np.sum(labels_individual == label) for label in np.sort(label_pop)]
        end_position = np.cumsum(individual_per_pop)
        start_position = np.insert(end_position, 0, 0, axis=0)
        delta_reorder = np.copy(delta)
        delta_reorder = delta_reorder[sorting_index, :]
        delta_reorder = delta_reorder[:, sorting_index]
        plt.figure()
        plt.imshow(delta_reorder)
        plt.tick_params(bottom=False, top=True, labeltop=True, labelbottom=False)

        for index_position in range(len(start_position)-1):
            plt.text((start_position[index_position] + start_position[index_position+1])/2, -2,
                     np.sort(label_pop)[index_position],
                        verticalalignment = 'bottom',
                     horizontalalignment='center',
                     rotation=90
                     )
            plt.text(-2,(start_position[index_position] + start_position[index_position+1])/2,
                     np.sort(label_pop)[index_position],
                        verticalalignment = 'center',
                     horizontalalignment='right',
                     )

        plt.xticks(start_position-1/2,"")

        plt.yticks(start_position-1/2, "")



        filePath = os.path.join(self.output_folder, "plot_distance.pdf")
        plt.savefig(filePath)
        plt.close()
        label_given = np.array(self.labels)
        if (self.K < 9):
            prop_cycle = plt.rcParams['axes.prop_cycle']
            prop_cycle = prop_cycle*(1+len(np.unique(label_given))//len(prop_cycle))
            colors = prop_cycle.by_key()['color']
        else:
            cmap = plt.get_cmap('jet')
            colors = cmap(np.linspace(0, 1.0, self.K))
        for population_index, population_label in enumerate(label_pop):
            population_position = labels_individual == population_label
            pop_mat = delta[np.ix_(population_position, population_position)]
            all_pop_value = pop_mat.flatten()
            all_pop_value = all_pop_value[all_pop_value > 0.00000001]
            plt.hist(all_pop_value, 15, label=population_label, density=1, alpha=0.75,color = colors[population_index])
        plt.legend()
        plt.xlabel("Allele sharing distance")
        plt.ylabel("# pairwise hit")
        filePath = os.path.join(self.output_folder, "Time_per_pop.pdf")
        plt.savefig(filePath)
        plt.close()
        label_given = nb_population = len(label_pop)

        if nb_population < 4:
            plt.figure()
            for pop1_index in range(nb_population):
                for pop2_index in range(pop1_index, nb_population):
                    population_position1 = labels_individual == label_pop[pop1_index]
                    population_position2 = labels_individual == label_pop[pop2_index]
                    pop_mat = delta[np.ix_(population_position1, population_position2)]
                    all_pop_value = pop_mat.flatten()
                    all_pop_value = all_pop_value[all_pop_value > 0.00000001]
                    plt.hist(all_pop_value, 20, label=label_pop[pop1_index] + "-" + label_pop[pop2_index], density=1,
                             alpha=0.5)
            plt.legend(ncol=3)
            plt.xlabel("Allele sharing distance")
            plt.ylabel("# pairwise hit")
            plt.savefig(os.path.join(self.output_folder, f"all_pop.pdf"))
            plt.close()
        elif nb_population < 9:
            for pop1_index in range(nb_population):
                plt.figure()
                for pop2_index in range(nb_population):
                    population_position1 = np.where(labels_individual == label_pop[pop1_index])[0]
                    population_position2 = np.where(labels_individual == label_pop[pop2_index])[0]
                    pop_mat = delta[np.ix_(population_position2, population_position1)]
                    all_pop_value = pop_mat.flatten()
                    all_pop_value = all_pop_value[all_pop_value > 0.00000001]
                    plt.hist(all_pop_value, 20, label=label_pop[pop1_index] + "-" + label_pop[pop2_index], density=1,
                             alpha=0.5)
                plt.legend(ncol=2)
                plt.xlabel("Allele sharing distance")
                plt.ylabel("# pairwise hit")
                plt.savefig(os.path.join(self.output_folder, f"time_pop_{label_pop[pop1_index]}.pdf"))
                plt.close()


    def set_tree(self, tree: TreeNode):
        self.tree = tree
        self.tree_with_name = copy.deepcopy(tree)
        for leave in self.tree_with_name.tips():
            leave.name = self.populations[int(leave.name)]
        # Defining the name for the population split
        split_index_matrix = -np.ones((self.K, self.K), dtype='int16')
        constraints = []
        constraints_coal_time = []
        build_split_index_matrix(tree, split_index_matrix, constraints, constraints_coal_time)
        self.split_index_matrix = split_index_matrix
        self.split_names = []
        for index_split in range(self.K - 1):
            for row in split_index_matrix:
                if index_split in row:
                    group_pop_1 = np.where(row == index_split)[0]
                    group_pop_2 = np.where(index_split == split_index_matrix[group_pop_1[0]])[0]
                    self.split_names.append((group_pop_1,group_pop_2))
                    break



    def plot_mds(self, coordinate, title: str, dir: str, labels_override=None):
        """Generic scatter plot for 2D component pairs (MDS/PCA)."""
        use_colors = self.has_labels or labels_override is not None
        label_given = np.array(self.labels if labels_override is None else labels_override)
        unique_labels = np.unique(label_given)
        nb_groups = len(unique_labels)

        if use_colors:
            if nb_groups < 9:
                prop_cycle = plt.rcParams['axes.prop_cycle']
                prop_cycle = prop_cycle*(1+nb_groups//len(prop_cycle))
                colors = prop_cycle.by_key()['color']
            else:
                cmap = plt.get_cmap('jet')
                colors = cmap(np.linspace(0, 1.0, nb_groups))

        n_dim = coordinate.shape[1]
        if self.plot_dims is not None:
            for p, q in self.plot_dims:
                if p >= n_dim or q >= n_dim:
                    sys.exit(
                        f"Error: --plot-dims requested dimension pair ({p+1},{q+1}) but only "
                        f"{n_dim} dimension(s) were computed. "
                        f"Increase --n-dims or reduce the requested dimensions."
                    )
            pairs = self.plot_dims
        else:
            pairs = [(p, p + 1) for p in range(0, n_dim - 1, 2)]
        for p, q in pairs:
            plt.rcParams.update({'font.size': 22})
            fig, ax = plt.subplots(figsize=(15, 15))
            if use_colors:
                for population_index, population_name in enumerate(unique_labels):
                    position_population = np.where(population_name == label_given)
                    color_value = colors[population_index]
                    ax.scatter(coordinate.T[p, position_population].ravel(), coordinate.T[q, position_population].ravel(), c=color_value, s=75, alpha=0.6)
                markers_color = [mlines.Line2D([], [], color=c, marker="o", linestyle='None') for c in colors]
                nb_column = nb_groups//14 + 1
                plt.legend(markers_color, unique_labels, title="Population", ncol=nb_column, bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
            else:
                ax.scatter(coordinate.T[p].ravel(), coordinate.T[q].ravel(), s=75, alpha=0.6)
            dir_plot = os.path.join(self.output_folder, dir)
            if not os.path.isdir(dir_plot):
                os.mkdir(dir_plot)
            ax.set_xlabel(f'PC. {p+1}')
            ax.set_ylabel(f'PC. {q+1}')
            fig.savefig(os.path.join(dir_plot, f'{title}{p+1}_{q+1}.pdf'),bbox_inches="tight")
            plt.close()


    def plot_diploid_mds_pca(self):
        """Compute and plot MDS and PCA using PLINK routines.

        Returns coords_mds so callers can use it for tree inference.
        """

        bfile_prefix = self.snp_file[:-4]  # remove .bed

        with tempfile.TemporaryDirectory(prefix="medeas_plink_") as tmpdir:
            out_prefix = os.path.join(tmpdir, "plink_run")
            out_prefix = self.asd_pattern.format(1)
            # --- IBS distance matrix via PLINK ---
            dist_path = _ibs_distance_plink(bfile_prefix, out_prefix, self.plink_path, "", self.NCORE)

            # --- MDS via PLINK ---
            coords_mds, eigenvalues_mds, _ = _mds_plink(
                bfile_prefix, self.n_dims, dist_path, self.plink_path, "", self.NCORE, tmpdir
            )
            _plot_eigenvalues_array(eigenvalues_mds, os.path.join(self.output_folder, "eigenvalues.pdf"))
            if not self.no_mds:
                self.plot_mds(coords_mds, "plink_MDS_", "mds_plot", labels_override=self.labels)

            # --- PCA via PLINK ---
            coords_pca, _ = _compute_pca_plink(bfile_prefix, self.n_dims, self.plink_path, self.NCORE, tmpdir)
            if not self.no_pca:
                self.plot_mds(coords_pca, "plink_PCA_", labels_override=self.labels)

        self.plink_coords_mds = coords_mds
        self.plink_coords_pca = coords_pca
        return coords_mds


    def save_tree(self):
        """Write the information about the infered tree into a file"""
        tree_filename = os.path.join(self.output_folder, "tree.txt")
        with open(tree_filename, "w") as f:
            f.write(self.tree_with_name.ascii_art())
            f.write("\n")
            f.write(str(self.tree_with_name))


    def generate_initial_output(self,args):
        with open(self.logfile, "w") as f:
            self.starting_time = datetime.datetime.now().replace(microsecond=0)
            f.write(f'starting new simulation at time: {self.starting_time} \n')
            try:
                label = subprocess.check_output(["git", "describe","--always"]).strip()
                f.write(f'you are using commit: {label}\n')
            except:
                f.write(f'No git hash tag detected \n')
            f.write("the following line was used to launch the simulation: \n")
            f.write(" ".join(sys.argv)+"\n")
            f.write("This led to the following argument being actually used: \n")
            f.write("\n \n" + "".join(100 * ["*"]) + "\n")
            for arg in vars(args):
                f.write(f'{arg}: {getattr(args, arg)} \n')
            f.write("\n" + "".join(100 * ["*"]) + "\n \n")


    def get_bootstraped_value(self, all_value):
        all_value = np.sort(np.array(all_value))
        nb_element = len(all_value)
        upper_index = int(0.975*nb_element)
        lower_index = int(0.025*nb_element)
        median_index = int(0.5*nb_element)
        return((all_value[median_index],all_value[lower_index],all_value[upper_index]))


    def get_bootstraped_value_matrix(self,matrix_to_bootstrap):
        all_bootstraped_value = []
        for column in matrix_to_bootstrap:
            all_bootstraped_value.append(self.get_bootstraped_value(column))
        return all_bootstraped_value


    def write_single_split(self, split_name, file):
        file.write("(")
        file.write("-".join(self.populations[split_name[0]]))
        file.write(")/(")
        file.write("-".join(self.populations[split_name[1]]))
        file.write(")")


    def write_header_split(self, file):
        for split_name in self.split_names:
            self.write_single_split(split_name, file)
            file.write("\t")
        file.write("\n")

    def write_header_pop(self, file):
        for population in self.populations:
            file.write(population + "\t")
        file.write("\n")


    def write_raw_computed_value(self):
        with open(os.path.join(self.output_folder, "between_population_coalescence_time.txt"), 'w') as f:
            self.write_header_split(f)
            for distances in self.all_between_pop_coalescence_time:
                for distance in distances:
                    f.write(f"{distance:.5f} \t")
                f.write("\n")
        with open(os.path.join(self.output_folder, "within_population_coalescence_time.txt"), 'w') as f:
            self.write_header_pop(f)
            for effective_sizes in self.all_within_pop_coalescence_time:
                for effective_size in effective_sizes:
                    f.write(f"{effective_size:.5f} \t")
                f.write("\n")
        with open(os.path.join(self.output_folder, "split_time.txt"), 'w') as f:
            self.write_header_split(f)
            for distances in self.all_split_time:
                for distance in distances:
                    f.write(f"{distance:.5f} \t")
                f.write("\n")
        with open(os.path.join(self.output_folder, "effective_size.txt"), 'w') as f:
            self.write_header_pop(f)
            for effective_sizes in self.all_effective_size:
                for effective_size in effective_sizes:
                    f.write(f"{effective_size:.5f} \t")
                f.write("\n")


    def write_boostraped_value(self):
        with open(os.path.join(self.output_folder, "split_bootstraped_confidence_interval.txt"), 'w') as f:
            f.write("node\t Coalescence time 50% [2.5%, 97.5%] \t split time 50% [2.5%, 97.5%] \n")
            all_bootstraped_coalescence = self.get_bootstraped_value_matrix(list(zip(*self.all_between_pop_coalescence_time)))
            all_bootstraped_split = self.get_bootstraped_value_matrix(list(zip(*self.all_split_time)))
            for coalescence_time, split_time, split_name in zip(all_bootstraped_coalescence,
                                         all_bootstraped_split,
                                         self.split_names
                                         ):
                self.write_single_split(split_name, f)
                f.write(f"\t {coalescence_time[0]:.5f} [{coalescence_time[1]:.5f},{coalescence_time[2]:.5f}] ")
                f.write(f"\t {split_time[0]:.5f} [{split_time[1]:.5f},{split_time[2]:.5f}] ")
                f.write("\n")

        with open(os.path.join(self.output_folder, "population_bootstraped_confidence_interval.txt"), 'w') as f:
            f.write("node\t Coalescence time 50% [2.5%, 97.5%] \t Effective size after split 50% [2.5%, 97.5%] \n")
            all_bootstraped_coalescence = self.get_bootstraped_value_matrix(list(zip(*self.all_within_pop_coalescence_time)))
            all_bootstraped_effective_size = self.get_bootstraped_value_matrix(list(zip(*self.all_effective_size)))
            for coalescence_time, effective_size, population_name in zip(all_bootstraped_coalescence,
                                         all_bootstraped_effective_size,
                                         self.populations
                                         ):
                f.write(population_name)
                f.write(f"\t {coalescence_time[0]:.5f} [{coalescence_time[1]:.5f},{coalescence_time[2]:.5f}] ")
                f.write(f"\t {effective_size[0]:.5f} [{effective_size[1]:.5f},{effective_size[2]:.5f}] ")
                f.write("\n")


    def generate_final_output(self):
        print("generating final output")
        self.write_raw_computed_value()
        self.write_boostraped_value()

        with open(os.path.join(self.output_folder, "all_T.txt"), 'w') as f:
            for T in self.all_T:
                    f.write(str(T)+"\n")
        with open(self.logfile, "a") as f:
            self.end_time = datetime.datetime.now().replace(microsecond=0)
            f.write(f'Simulation ended successfully at: {self.end_time} \n')
            f.write(f'job duration:  {self.end_time - self.starting_time} \n')


    def generate_mds_pca_only_output(self):
        """Write outputs available when split-time inference is disabled."""
        with open(self.logfile, "a") as f:
            self.end_time = datetime.datetime.now().replace(microsecond=0)
            f.write("Split-time inference skipped (--no-split).\n")
            f.write(f'Simulation ended successfully at: {self.end_time} \n')
            f.write(f'job duration:  {self.end_time - self.starting_time} \n')

