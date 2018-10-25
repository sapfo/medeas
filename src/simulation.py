import argparse
import os
import numpy as np
import matplotlib.pyplot as plt

class simulation_info(object):



    def __init__(self):
        parser = argparse.ArgumentParser()

        parser.add_argument("-sf", "--snps_file",
                            help="Prepare the data wihtout performing the actual analysis on them",
                            required=True)
        parser.add_argument("-lf", "--labels_file", help="File containing the labels",
                            required=True)
        parser.add_argument("-of", "--output_folder", help="Folder where results and temporal data should be store")
        parser.add_argument("-K",
                            help="Number of population. If K=0 (default), the soft detect automatically the number of population",
                            type=int, default=0)
        parser.add_argument("--n_chromosome", help="Number of chromosome If they are stored in different file",
                            type=int, default=1)
        parser.add_argument("--outgroup", help="Who is the outgroup in your data", nargs='+')

        parser.add_argument("-af", "--ancestry_file", help="File containing the ancestry of each locus")

        parser.add_argument("-bws", "--boot_window_size",
                            help="How many markers do we have in each bootstraping windows",
                            type=int, default=100)

        parser.add_argument("--simulation", help="Does the data come from a simulation",
                            action="store_true")
        parser.add_argument("--skip_calculate_matrix",
                            help="Skip the computation of the distance matrices and the related MDS matrix",
                            action="store_true")
        parser.add_argument("--skip_preprocessing", help="Directly proceed to analysis without preparing the data",
                            action="store_true")
        parser.add_argument("--skip_analysis", help="Prepare the data without performing the actual analysis on them",
                            action="store_true")
        args = parser.parse_args()
        self.ancestry_pattern = args.ancestry_file

        self.snps_pattern = args.snps_file
        self.output_folder = args.output_folder
        self.labels_file = args.labels_file
        self.chromosomes = range(1, args.n_chromosome + 1)
        self.bootsize = args.boot_window_size
        self.output_file = os.path.join(self.output_folder, 'processed', 'chr.{}.stped')
        self.K = args.K
        self.outgroups = args.outgroup
        self.skip_calculate_matrix = args.skip_calculate_matrix
        self.simulation = args.simulation
        self.skip_analysis = args.skip_analysis

        asd_folder = "asd_matrices"
        mds_folder = "MDS_eigensystem"
        asd_full_path = os.path.join(self.output_folder,asd_folder)
        mds_full_path = os.path.join(self.output_folder,mds_folder)
        all_path = [asd_full_path,mds_full_path]
        for path in all_path:
            if not os.path.exists(path):
                os.makedirs(path)
        self.asd_pattern = os.path.join(asd_full_path, 'p{}.asd.data')
        self.vec_pattern = os.path.join(mds_full_path, 'p{}.vecs.data')



    def generate_output(self):
        with open(os.path.join(self.output_folder, "all_extrapolated_distances.txt"), 'w') as f:
            np.savetxt(f, self.all_res)

    def plot_distance_matrix(self, delta):
        with open(self.labels_file) as f:
            lines = f.readlines()

        labels_individual = np.array([l.split()[0] for l in lines])
        sorted_labels_individual = np.sort(labels_individual)
        label_pop = np.unique(labels_individual)
        sorting_index = np.argsort(labels_individual)
        individual_per_pop = [np.sum(labels_individual == label) for label in np.sort(label_pop)]
        end_position = np.cumsum(individual_per_pop)
        start_position = np.insert(end_position, 0, 0, axis=0)
        print(start_position)
        delta = delta[sorting_index, :]
        delta = delta[:, sorting_index]
        plt.figure()
        plt.imshow(delta)
        plt.tick_params(bottom=False, top=True, labeltop=True, labelbottom=False)
        plt.xticks(start_position, np.sort(label_pop), rotation='vertical')
        plt.yticks(start_position, np.sort(label_pop))
        plt.savefig("plot_distance.pdf")
        plt.figure()

        for population_label in label_pop:
            population_position = sorted_labels_individual == population_label
            pop_mat = delta[np.ix_(population_position, population_position)]
            all_pop_value = pop_mat.flatten()
            all_pop_value = all_pop_value[all_pop_value > 0.00000001]
            plt.hist(all_pop_value, 15, label=population_label, density=1, alpha=0.75)
        plt.legend()
        plt.savefig("Time_per_pop.pdf")

        nb_population = len(label_pop)
        for pop1_index in range(nb_population):
            plt.figure()
            for pop2_index in range(nb_population):
                population_position1 = labels_individual == label_pop[pop1_index]
                population_position2 = labels_individual == label_pop[pop2_index]
                pop_mat = delta[np.ix_(population_position1, population_position2)]
                all_pop_value = pop_mat.flatten()
                all_pop_value = all_pop_value[all_pop_value > 0.00000001]
                plt.hist(all_pop_value, 20, label=label_pop[pop1_index] + "-" + label_pop[pop2_index], density=1,
                         alpha=0.5)
            plt.legend(ncol=2)
            plt.savefig(f"time_pop_{label_pop[pop1_index]}.pdf")
