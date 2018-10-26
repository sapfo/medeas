import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines



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
        parser.add_argument("-bsn","--bootstrap_number",
                            help="How many bootstrap do we perform",
                            type=int, default=10
                            )


        parser.add_argument("--simulation", help="Does the data come from a simulation",
                            action="store_true")
        parser.add_argument("--skip_calculate_matrix",
                            help="Skip the computation of the distance matrices and the related MDS matrix",
                            action="store_true")
        parser.add_argument("--skip_preprocessing", help="Directly proceed to analysis without preparing the data",
                            action="store_true")
        parser.add_argument("--skip_analysis", help="Prepare the data without performing the actual analysis on them",
                            action="store_true")

        parser.add_argument("--output_level", help="How many information should be printed & saved: 0 -minimal, 1 - conventional, 2 - most of it",
                            type=int, default=1)

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
        self.bootstrap_number = args.bootstrap_number
        self.output_level = args.output_level

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
        with open(self.labels_file) as f:
            lines = f.readlines()

        self.labels = [l.split()[0] for l in lines]  # + ['WCD'] * 7



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


    def plot_mds(self, arr, labels_inferred):
        """Plot the MDS plot
        """
        # TODO: autogenerate nice summary plot depending on 'npop'
        label_given = np.array(self.labels)
        label_given_index = np.copy(label_given)
        for index_label, label in enumerate(np.unique(label_given)):
            label_given_index[label_given == label] = index_label
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        markers = [".", "v", "*", "+", "x", "2", "p", "^", "s"]
        for p in range(len(np.unique(labels_inferred))+2):
            for q in range(p+1,len(np.unique(labels_inferred))+2):
                fig, ax = plt.subplots(figsize=(15,10))
                for population_index, population_name in enumerate(np.unique(label_given)):
                    position_population = np.where(population_name == label_given)
                    index_colors = labels_inferred[position_population]
                    color_value = [colors[index_color] for index_color in index_colors]
                    markers_value = markers[population_index]
                    ax.scatter(arr.T[p,position_population].ravel(), arr.T[q, position_population].ravel(), c = color_value, marker = markers_value, s=100)
                plt.legend(np.unique(label_given))
                leg = ax.get_legend()
                for point in leg.legendHandles:
                    point.set_color('black')
                dir_plot = os.path.join(self.output_folder, "mds_plot")
                if not os.path.isdir(dir_plot):
                    os.mkdir(dir_plot)
                markers_shape = [mlines.Line2D([], [], color='black', marker=marker_shape, linestyle='None') for marker_shape in markers]
                markers_color = [mlines.Line2D([], [], color=marker_color, marker="o", linestyle='None') for marker_color in colors]
                legend1 = plt.legend(markers_shape,np.unique(label_given) , loc=1,title="True population")
                legend2 = plt.legend(markers_color,np.unique(labels_inferred), loc=2,title="Inferred population")
                ax.add_artist(legend1)
                ax.add_artist(legend2)
                ax.set_xlabel(f'PC. {p+1}')
                ax.set_ylabel(f'PC. {q+1}')
                fig.savefig(os.path.join(dir_plot, f'mds_axis{p+1}_{q+1}.pdf'))