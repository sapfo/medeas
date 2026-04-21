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


from skbio.tree import nj, TreeNode
import pickle
from src.clustering import get_mds_coordinate
from src.clustering import build_split_index_matrix


# ---------------------------------------------------------------------------
# PLINK-based MDS/PCA helpers (adapted from mds2split.py)
# ---------------------------------------------------------------------------

def _ibs_distance_plink(bfile: str, dist_out: str, plink_path: str,
                        plink_params: str, threads: int) -> str:
    """Compute pairwise IBS distance file using PLINK --distance."""
    prefix = dist_out
    cmd = [
        plink_path,
        "--bfile", bfile,
        "--allow-no-sex",
        "--distance", "square", "gz", "flat-missing",
        "--threads", str(threads),
        "--out", prefix,
    ]
    if plink_params:
        cmd += plink_params.split()

    print(f"[mds2split] Running: {' '.join(cmd)}", file=sys.stderr)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.stderr.write(result.stderr)
        raise RuntimeError(f"PLINK distance computation failed (exit {result.returncode})")

    # flat-missing outputs are typically .dist(.gz); keep alternate stems for compatibility.
    candidates = [
        f"{prefix}.dist.gz",
        f"{prefix}.dist",
        f"{prefix}.mdist.gz",
        f"{prefix}.mdist",
        f"{prefix}.mibs.gz",
        f"{prefix}.mibs",
    ]
    for path in candidates:
        if os.path.isfile(path):
            return path

    raise RuntimeError(
        "PLINK distance output file not found after successful run. "
        f"Checked: {', '.join(candidates)}"
    )


def _write_read_dists_binary(D: np.ndarray, out_path: str):
    """Write PLINK --read-dists format: float64 lower triangle, no diagonal."""
    D = np.asarray(D, dtype=np.float64)
    n = D.shape[0]
    with open(out_path, "wb") as f:
        for i in range(1, n):
            D[i, :i].tofile(f)


def _count_snps_from_bim(bfile: str) -> int:
    """Count SNPs from the PLINK BIM file."""
    bim_path = f"{bfile}.bim"

    with open(bim_path, "r") as f:
        n_snps = sum(1 for _ in f)

    if n_snps <= 0:
        raise RuntimeError(f"No SNPs found in BIM file: {bim_path}")
    print(f"Counted {n_snps} SNPs from BIM file: {bim_path}")
    return n_snps


def _mds_plink(bfile: str, k: int, dist_path: str, plink_path: str,
               plink_params: str, threads: int, tmpdir: str):
    """Run PLINK MDS from a precomputed IBS distance matrix.

    Returns (coords, eigenvalues, samples)."""
    prefix = os.path.join(tmpdir, "plink_mds")
    read_dists_path = os.path.join(tmpdir, "read_dists.bin")

    D = np.loadtxt(dist_path)
    n_snps = _count_snps_from_bim(bfile)
    #D = D / float(n_snps)
    np.fill_diagonal(D, 0.0)
    _write_read_dists_binary(D, read_dists_path)

    mds_cmd = [
        plink_path,
        "--bfile", bfile,
        "--read-dists", read_dists_path,
        "--allow-no-sex",
        "--cluster",
        "--mds-plot", str(k), "eigvals",
        "--threads", str(threads),
        "--out", prefix,
    ]
    if plink_params:
        mds_cmd += plink_params.split()

    print(f"[mds2split] Running: {' '.join(mds_cmd)}", file=sys.stderr)
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

    print(f"[mds2split] Running: {' '.join(pca_cmd)}", file=sys.stderr)
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



class SimulationInfo(object):

    def __init__(self):
        parser = argparse.ArgumentParser()

        parser.add_argument("-sf", "--snps_file",
                            help="The name of the file from which the pattern should be read. ")
        parser.add_argument("-lf", "--labels_file", help="File containing the labels")
        parser.add_argument("-of", "--output_folder", help="Folder where results and temporal data should be store")

        parser.add_argument("-bws", "--boot_window_size",
                            help="How many markers do we have in each bootstraping windows",
                            type=int, default=100)
        parser.add_argument("-bsn","--bootstrap_number",
                            help="How many bootstrap do we perform",
                            type=int, default=100
                            )

        parser.add_argument("-t","--topology",
                            help="What is the topology of the population (newick format, following label order)",
                            type=str, default=None
                            )

        parser.add_argument("--skip_calculate_matrix",
                            help="Skip the computation of the distance matrices and the related MDS matrix",
                            action="store_true")

        parser.add_argument("--output_level", help="How many information should be printed & saved: 0 -minimal, 1 - conventional, 2 - most of it",
                            type=int, default=0)

        parser.add_argument("--threads", help="Number of parallel process to be launch. 0 (default) used all available cores",
                            type=int, default=0)

        parser.add_argument("--no_split",
                            help="Skip tree and split-time estimation; compute only distance, MDS, and PCA outputs",
                            action="store_true")

        parser.add_argument("--detailed_output",
                            help="Return detailed output",
                            action="store_true")

        parser.add_argument("--bfile",
                            help="PLINK binary file prefix (.bed/.bim/.fam). "
                                 "If provided, --snps_file and --labels_file are not needed; "
                                 "the PLINK data is converted automatically.",
                            type=str, default=None)

        parser.add_argument("--max_snps",
                           help="Maximum number of SNPs to process (for faster tests). "
                               "By default, process all SNPs.",
                           type=int, default=None)

        parser.add_argument("--diploid",
                            help="Plot MDS and PCA using PLINK routines (for diploid data). "
                                 "Requires --bfile.",
                            action="store_true")

        parser.add_argument("--plink_path",
                            help="Path to the PLINK executable (default: 'plink')",
                            type=str, default="plink")

        parser.add_argument("-k", "--k",
                    help="Number of eigenvalues to compute and plot (default: all, i.e. n_samples-1)",
                    type=int, default=0)

        args = parser.parse_args()

        if args.bfile is None and (args.snps_file is None or args.labels_file is None):
            sys.exit("Error: either --bfile or both --snps_file and --labels_file must be provided.")

        if args.diploid and args.bfile is None:
            sys.exit("Error: --diploid is only supported with PLINK input (--bfile).")

        if args.diploid and (args.snps_file is not None or args.labels_file is not None):
            sys.exit("Error: --diploid requires plink input.")

        if args.k < 0:
            sys.exit("Error: -k/--k must be >= 0 (0 means all dimensions).")

        # Create output folder first (required before bfile conversion)
        self.output_folder = args.output_folder
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

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

            self.bfile_prefix = bfile_prefix
            self.snps_pattern = None

            # Build haploid labels directly from FAM (each diploid sample -> 2 labels).
            labels_file = os.path.join(self.output_folder, "labels_from_fam.dat")
            with open(fam_file) as fin, open(labels_file, "w") as fout:
                for line in fin:
                    parts = line.strip().split()
                    if not parts:
                        continue
                    fid = parts[0]
                    fout.write(fid + "\n")
                    fout.write(fid + "\n")
            self.labels_file = labels_file
        else:
            self.bfile_prefix = None
            self.snps_pattern = args.snps_file
            if not os.path.isfile(self.snps_pattern):
                sys.exit("Error: The file containing the genotype does not exist. Exiting Now.")
            self.labels_file = args.labels_file
            if not os.path.isfile(self.labels_file):
                sys.exit("Error: The file containing the label does not exist. Exiting Now.")

        self.bootsize = args.boot_window_size
        self.skip_calculate_matrix = args.skip_calculate_matrix
        self.bootstrap_number = args.bootstrap_number
        self.output_level = args.output_level
        self.NCORE = args.threads
        if self.NCORE == 0:
            self.NCORE = cpu_count()
        self.topology = args.topology
        self.no_split = args.no_split
        self.detailed_output = args.detailed_output
        self.max_snps = args.max_snps
        if self.max_snps is not None and self.max_snps <= 0:
            sys.exit("Error: --max_snps must be a positive integer.")

        self.diploid = args.diploid
        self.plink_path = args.plink_path
        self.plink_k = args.k
        if self.diploid and self.bfile_prefix is None:
            sys.exit("Error: --diploid requires --bfile to be provided.")

        self.logfile = os.path.join(self.output_folder, "simulation.log")
        self.generate_initial_output(args)

        if self.diploid:
            # Diploid mode uses PLINK routines directly; ASD/MDS eigensystem folders are not needed.
            self.asd_pattern = None
            self.vec_pattern = None
        else:
            asd_folder = "asd_matrices"
            mds_folder = "MDS_eigensystem"
            asd_full_path = os.path.join(self.output_folder,asd_folder)
            mds_full_path = os.path.join(self.output_folder,mds_folder)
            all_path = [asd_full_path,mds_full_path]
            print(os.getcwd())
            for path in all_path:
                if not os.path.exists(path):
                    os.makedirs(path)
            self.asd_pattern = os.path.join(asd_full_path, 'p{}.asd.data')
            self.vec_pattern = os.path.join(mds_full_path, 'p{}.vecs.data')

        try:
            with open(self.labels_file) as f:
                lines = f.readlines()
        except:
            sys.exit(
                "Error: A problem occurs when loading the label file. Please check that you use the right format.")

        #labels = [l.split()[0] for l in lines]
        labels = [l.rstrip() for l in lines]
        self.labels = np.array(labels)
        _, index, numerical_labels = np.unique(self.labels,return_inverse=True, return_index = True)
        self.populations = self.labels[np.sort(index)]
        self.numerical_labels = np.sort(numerical_labels)
        self.K = len(self.populations)

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



    def plot_mds(self, coordinate, title: str, labels_override=None):
        """Generic scatter plot for 2D component pairs (MDS/PCA), haploid or diploid."""
        label_given = np.array(self.labels if labels_override is None else labels_override)
        unique_labels = np.unique(label_given)
        nb_groups = len(unique_labels)

        if nb_groups < 9:
            prop_cycle = plt.rcParams['axes.prop_cycle']
            prop_cycle = prop_cycle*(1+nb_groups//len(prop_cycle))
            colors = prop_cycle.by_key()['color']
        else:
            cmap = plt.get_cmap('jet')
            colors = cmap(np.linspace(0, 1.0, nb_groups))

        n_dim = coordinate.shape[1]
        for p in range(0, n_dim - 1, 2):
            q = p + 1
            plt.rcParams.update({'font.size': 22})
            fig, ax = plt.subplots(figsize=(15, 15))
            for population_index, population_name in enumerate(unique_labels):
                position_population = np.where(population_name == label_given)
                color_value = colors[population_index]
                ax.scatter(coordinate.T[p, position_population].ravel(), coordinate.T[q, position_population].ravel(), c=color_value, s=75, alpha = 0.6)
            plt.legend(unique_labels)
            leg = ax.get_legend()
            for point in leg.legend_handles:
                point.set_color('black')
            dir_plot = os.path.join(self.output_folder, "mds_plot")
            if not os.path.isdir(dir_plot):
                os.mkdir(dir_plot)
            markers_color = [mlines.Line2D([], [], color=marker_color, marker="o", linestyle='None') for marker_color in colors]
            nb_column = nb_groups//14 + 1
            plt.legend(markers_color, unique_labels, title="Population",ncol=nb_column,bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)

            ax.set_xlabel(f'PC. {p+1}')
            ax.set_ylabel(f'PC. {q+1}')
            fig.savefig(os.path.join(dir_plot, f'{title}{p+1}_{q+1}.pdf'),bbox_inches="tight")
            plt.close()

    def plot_diploid_mds_pca(self):
        """Plot MDS and PCA for diploid data using PLINK routines.

        Uses the same plotting function as the haploid version (plot_mds).
        """
        if self.bfile_prefix is None:
            sys.exit("Error: diploid plotting is only supported with PLINK input (--bfile).")

        bfile = self.bfile_prefix
        plink_path = self.plink_path
        threads = self.NCORE

        # Build diploid (one-per-sample) labels from FAM FID column
        fam = pd.read_csv(f"{bfile}.fam", sep=r"\s+", header=None,
                          names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
                          dtype=str)
        diploid_labels = np.array(fam["FID"].tolist())

        # k=0 means all available dimensions.
        if self.plink_k == 0:
            k = max(1, len(diploid_labels) - 1)
        else:
            k = self.plink_k

        with tempfile.TemporaryDirectory(prefix="medeas_diploid_") as tmpdir:
            out_prefix = os.path.join(tmpdir, "diploid_plink")
            # --- IBS distance via PLINK ---
            dist_path = _ibs_distance_plink(bfile, out_prefix, plink_path, "", threads)

            # --- MDS via PLINK ---
            coords_mds, eigenvalues_mds, _ = _mds_plink(
                bfile, k, dist_path, plink_path, "", threads, tmpdir
            )
            _plot_eigenvalues_array(eigenvalues_mds, os.path.join(self.output_folder, "eigenvalues.pdf"))
            self.plot_mds(coords_mds, "plink_MDS_", labels_override=diploid_labels)

            # --- PCA via PLINK ---
            coords_pca, _ = _compute_pca_plink(bfile, k, plink_path, threads, tmpdir)
            self.plot_mds(coords_pca, "plink_PCA_", labels_override=diploid_labels)

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
        with open(os.path.join(self.output_folder, "MDS_coordinate.txt"), 'w') as f:
            np.savetxt(f,  get_mds_coordinate(self, 1))
        with open(os.path.join(self.output_folder, "PCA_coordinate.txt"), 'w') as f:
            np.savetxt(f,  get_mds_coordinate(self, 2))
        with open(self.logfile, "a") as f:
            self.end_time = datetime.datetime.now().replace(microsecond=0)
            f.write(f'Simulation ended successfully at: {self.end_time} \n')
            f.write(f'job duration:  {self.end_time - self.starting_time} \n')

    def generate_mds_pca_only_output(self):
        """Write outputs available when split-time inference is disabled."""
        print("generating mds/pca-only output")
        with open(os.path.join(self.output_folder, "MDS_coordinate.txt"), 'w') as f:
            np.savetxt(f, get_mds_coordinate(self, 1))
        with open(os.path.join(self.output_folder, "PCA_coordinate.txt"), 'w') as f:
            np.savetxt(f, get_mds_coordinate(self, 2))
        with open(self.logfile, "a") as f:
            self.end_time = datetime.datetime.now().replace(microsecond=0)
            f.write("Split-time inference skipped (--no_split).\n")
            f.write(f'Simulation ended successfully at: {self.end_time} \n')
            f.write(f'job duration:  {self.end_time - self.starting_time} \n')

