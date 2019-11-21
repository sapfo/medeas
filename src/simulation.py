import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import subprocess
import datetime
import sys
from multiprocessing import cpu_count


from skbio.tree import nj, TreeNode
import pickle
from src.clustering import get_mds_coordinate
from src.clustering import build_split_index_matrix
class SimulationInfo(object):

    def __init__(self):
        parser = argparse.ArgumentParser()

        parser.add_argument("-sf", "--snps_file",
                            help="The name of the file from which the pattern should be read. ",
                            required=True)
        parser.add_argument("-lf", "--labels_file", help="File containing the labels",
                            required=True)
        parser.add_argument("-of", "--output_folder", help="Folder where results and temporal data should be store")

        parser.add_argument("-bws", "--boot_window_size",
                            help="How many markers do we have in each bootstraping windows",
                            type=int, default=100)
        parser.add_argument("-bsn","--bootstrap_number",
                            help="How many bootstrap do we perform",
                            type=int, default=100
                            )

        parser.add_argument("-t","--topology",
                            help="What is the topology of the population (newick format, following label order",
                            type=str, default=None
                            )

        parser.add_argument("--skip_calculate_matrix",
                            help="Skip the computation of the distance matrices and the related MDS matrix",
                            action="store_true")

        parser.add_argument("--output_level", help="How many information should be printed & saved: 0 -minimal, 1 - conventional, 2 - most of it",
                            type=int, default=1)

        parser.add_argument("--ncpus", help="Number of parallel process to be launch. 0 (default) used all available cores",
                            type=int, default=0)

        args = parser.parse_args()

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
        self.NCORE = args.ncpus
        if self.NCORE == 0:
            self.NCORE = cpu_count()
        self.topology = args.topology

        self.output_folder = args.output_folder
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
        self.logfile = os.path.join(self.output_folder, "simulation.log")
        self.generate_initial_output(args)

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

    def plot_eigenvalues(self):
        with open(self.vec_pattern.format(2), 'rb') as f:
            lambdas, vecs = pickle.load(f)
        lambdas = -np.sort(-lambdas)
        plt.subplot(211)
        plt.plot(lambdas,"o")

        plt.ylabel("Eigenvalues")
        plt.subplot(212)
        plt.plot(lambdas[self.K:-2],"o")
        plt.ylabel("Eigenvalues")
        plt.xlabel("Eigenvalues index")
        filePath = os.path.join(self.output_folder, "eigenvalues.pdf")
        plt.savefig(filePath)
        plt.close()

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

    def set_tree(self, tree: 'skbio.tree'):
        self.tree = tree
        self.tree_with_name = tree.deepcopy()
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



    def plot_mds(self, coordinate, title: str):
        """Plot the MDS plot
        """
        label_given = np.array(self.labels)
        label_given_index = np.copy(label_given)
        for index_label, label in enumerate(np.unique(label_given)):
            label_given_index[label_given == label] = index_label
        if (self.K < 9):
            prop_cycle = plt.rcParams['axes.prop_cycle']
            prop_cycle = prop_cycle*(1+len(np.unique(label_given))//len(prop_cycle))
            colors = prop_cycle.by_key()['color']
        else:
            cmap = plt.get_cmap('jet')
            colors = cmap(np.linspace(0, 1.0, self.K))

        for p in range(0, self.K + 2, 2):
        #for p in range(0, len(self.labels)-1, 2):
            q = p + 1
            plt.rcParams.update({'font.size': 22})
            fig, ax = plt.subplots(figsize=(15, 15))
            for population_index, population_name in enumerate(np.unique(label_given)):
                position_population = np.where(population_name == label_given)
                color_value = colors[population_index]
                ax.scatter(coordinate.T[p, position_population].ravel(), coordinate.T[q, position_population].ravel(), c=color_value, s=75, alpha = 0.6)
            plt.legend(np.unique(label_given))
            leg = ax.get_legend()
            for point in leg.legendHandles:
                point.set_color('black')
            dir_plot = os.path.join(self.output_folder, "mds_plot")
            if not os.path.isdir(dir_plot):
                os.mkdir(dir_plot)
            markers_color = [mlines.Line2D([], [], color=marker_color, marker="o", linestyle='None') for marker_color in colors]
            nb_column = self.K//14 + 1
            plt.legend(markers_color, np.unique(label_given) , title="Population",ncol=nb_column,bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)

            ax.set_xlabel(f'PC. {p+1}')
            ax.set_ylabel(f'PC. {q+1}')
            fig.savefig(os.path.join(dir_plot, f'{title}{p+1}_{q+1}.pdf'),bbox_inches="tight")
            plt.close()

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

