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



