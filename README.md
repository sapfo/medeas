# Documentation of Medeas
## General idea
Medeas is a python software that allows to infer coalescence time between individuals within and between populations. It also reconstructs the population tree, but requires labeled individuals, i.e. the populations to which individuals belong to need to be specified.

The goal is achieved by computing a multi-dimensional scaling (MDS) projection in a K-dimensional space and then inferring coalescence times from the eigenvalues of the MDS. The population tree is inferred using the center of mass of the various clusters in the MDS projection. 

## Usage
To run Medeas, you need python 3.6 or higher, as well as a the following packages: 
- numpy, 
- matplotlib, 
- scikit-learn, 
- scikit-bio
- scipy

Medeas requires two input files:
- A matrix of genotypes
- A file indicating the population of each individual

For now, the method works only on haploid individuals. Our recommendation is use this method on diploid individuals by considering the two sets of chromosomes independently. 
### Genotypes matrix format
The genotypes matrix is a text file, where each column is one haploid individual and each row is a position in the genome. No header nor row names should appear. Reference allele should be encoded as 1 and alternative allele as 2. Missing data are encoded with 0. 
### Label file format
The label file should contain the population of each haploid individual, one by line, without any other information. The number of lines in this file should be the same as the number of columns in the genotype matrix file. 
### Mandatory parameters

Three arguments are mandatory to launch Medeas:
- `-sf`: path to the genotype file
- `-lf`:  path to the label file
- `-of`: path to the ouput folder

### Example
An example genotypes matrix and a label file are given in the test directory ( [test/snp.dat](test/snp.dat) and  [test/pop_label.dat](test/pop_label.dat), respectively). Using them, we can launch Medeas using: 
```
python main.py -sf test/snp.dat -lf test/pop_label.dat -of test
```
### Optional parameters:
- `-bws`: The number of SNPs to include in a bootstrapping window. Notice that the largest this number is, the faster it is, but if it’s too big, the bootstrapping intervals become less reliable. Our recommendation is to use about 1% of the SNPs. 
- `-bsn`: the number of bootstrapping to perform. We recommend either to use 0 (i.e. just do one calculation with all the SNP) or 100 to have a straight way to get a 95% interval
- `-t`: The topology of the tree can be enforced by specifying this option. It should be given using the Newick format
- `--skip_calculate_matrix`: If the software was already launched, but you want to relaunch it using different parameters (for example changing the topology or the population names), you can use this option. In this case, the genotype file would not be re-read but the calculation will use some internally stored distance matrix which would make the calculation much faster if your genotype matrix is large.
- `--output_level`: This option can be changed to decide how much information is displayed to the user while running the simulation: 0 = minimal, 1 = conventional, 2 = maximum
- `-h`: Display a help message.


## Output
Medeas will generate a folder containing various information. Here is a list of the output files:
- `mds_plot`: folder containing MDS and PCA plots of the first K+2 dimension (where K is the number of populations contained in the label file)
- `between_population_coalescence_time.txt`: file giving the coalescence time between each groups of population (all bootstrap values)
- `within_population_coalescence_time.txt`: file giving the coalescence time within population. Notice that by definition (see publication) the coalescence time in the first population is 1 (all bootstrap value)
- `effective_size.txt`: file giving the effective size of each population after splitting for the main branch (all bootstrap value)
- `split_time.txt`: file giving the split time for each split (all bootstrap value)
- `population_bootstraped_confidence_interval.txt`: file giving the confidence interval for the within population coalescence time as well as the effective size after splitting from the main branch. 
- `split_bootstraped_confidence_interval.txt`: file giving the confidence interval for the between population coalescence time as well as the split time
 
- `tree.txt`: file showing the inferred tree.
- `SFS.pdf`:  file showing the folded SFS of the entire panel. 
- `eigenvalues.pdf`: plot showing the eigenvalues of the MDS. It shows both, the eigenvalues (upper panel) and the bulk eigenvalues (lower panel)
- `histogram_eigenvalues.pdf`: A histogram of the bulk eigenvalues
- `Time_per_pop.pdf` and `time_*.pdf` or `all_pop.pdf`: Several plot showing the histogram of the pairwise distance of individual between and within populations. If K <= 3 (K = number of population), all the pairwise distance are depicted in  `all_pop.pdf`. If K > 3, separeted file `time_popName.pdf` for each population are created for a better readability. 
- `plot_distance.pdf`: A graphical representation of the distance matrix
- `(MDS|PCA)_coordinate.txt`: The actual coordinate of the MDS and PCA plot. If the user wants to redo the plot with a different layout. 
- `asd_matrices` and `MDS_eigensystem`: Folders used for internal storage. They are kept and are useful if the user wants to repeat the same calculation using the option --skip-calculate-matrix but can be deleted otherwise. 

## usefull Scripts:
Two scripts are also available to change the format from a vcf file to the internal format: 
- [script/Snakefile](script/Snakefile): This is a snakemake pipeline which first convert the file [script/CEU_YRI.chr22.vcf.gz](script/CEU_YRI.chr22.vcf.gz) and [script/label.lab](script/label.lab) into Medeas internal format and then run medeas with these files.  It is the recommended option if you have `snakemake` and `conda` installed on your system. The input vcf and label file can be changed in [script/config.yaml](script/config.yaml), and the script should be launched using:
```snakemake --use-conda```
Using this command, you don’t need to have python 3.6 and plink install, conda will take care of it for you. 
- [script/launch.sh](script/launch.sh): This is a bash script which first convert the file [script/CEU_YRI.chr22.vcf.gz](script/CEU_YRI.chr22.vcf.gz) and [script/label.lab](script/label.lab) into Medeas internal format and then run medeas with thess files. This script is easy to adapt for your own data. It requires plink and common UNIX tools to be installed

## Publication
Pending.

