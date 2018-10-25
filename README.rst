======
medeas
======
listing of directories, 02.01.18
================================
main.py: main entry point
names for temp files.

recode: add spaces between alleles (numpy likes it better)

filtering for real data:
"-freqs
"-filter
"-manual
"-hard
"-setmiss
"-sparse

"-asd: makes the distance matrix
"asd_main(p=1,2, etc)

"-analyze: the actual analyses

TODO: codes to strings in hardfilt
TODO: recalculate everything at every bootstrap run


the user cannot use the p currently (marchenko pastur p=2 ONLY, for clustering its unclear but currently only p=1)



src
several important steps

#src/make_asd.py
it makes the matrix
complex part parallelization and bootstrap
(TODO: kill some old stuff)
main idea is that it splits the file or data on chunks on given size and computes the matrix of unnormalized distances for every chunks of SNPs
with this data you can reconstruct many resamplings by adding the unnormalized distances and dividing them

parallelization is done currently on one of the indices of the matrices
now i take all the columns and i spawn some processes that take the colmn number and compute the pairwise distance

(a column is an indiv.)
(a row as well)

#src/lambda_analyze.py
two main functions
- find T and L
fits the eigenvalues and finds optimal fit of total tree length and independent number os sites
- find K, pvalue hardcoded to 1e-3

#src/clustering.py
really bad. poorly organised.
used a lot for testings. so loads of plots etc. need to clean and have a spearate big file to make plots.

3 main functions. first importan function. perform clustering. every individual is assigned a label.
(TODO: remove all the old stuff so its readable)
perform_clustering: where we assign the labels
find_tree: where we find tree
find_distances: most complex and important function. constrains. etc.
==========
tracy.py: not using yet, it gives pvalues as a  function of sigma (for now using a fixed pvalue so there is a single value of the stats)
=========
transcode.py: converts a filtered output of scrm (filtered for numbers only) to the format that medeas use and saves the seed to a hardcoded filename: _tmp_seed,txt
arguments: folder, whereiscrm, splittime /currently works with two pops, one D)

within transcode.py:
create fake_labs because medeas needs it

This file is the 'medeas' main entry point.

Here is the summary of how things are organized.
On large scale, there are two folders: 'procesing' and 'src'.
The first one contains files for preprocessing of the data before
actuall analyzis, this includes various filtering steps, format
transforms, generation of auxilliary files, setting missingness etc.
(see function docstrings for details). The second folder contains modules
for actual analyzis using the eigenvalue statistics.

Here is the summary of steps in normal data pre-processing:
1. Convert the input files with SNP data and ancestry data into internal
   binary format (currently just stdlib pickle).
2. Calculate ancestry frequencies per site.
3. Manually remove some populations/individuals (if desired).
3. Filter sites depending on ancetry frequencies, this removes the sites
   where ancestry is anomaluos (i.e. to far from average).
4. Filter sites by keeping only those where amount of given ancestry is
   larger than given level.
5. Set SNP data with given ancestry as missing.
6. Filter individuals by keeping only those with given amount of non-missing
   sites.

Here is the summary of actual analysis workflow:
1. Calculate distance matrixes for p=1 and p=2 (together with several
   bootstrapped -- resampled with replacement -- matrixes).
2. Calculate the eigensystem for two matrixes.
3. Find T (total tree length) and L (effective number of markers) by fitting
   the bulk distribution of eigenvalues for p=2.
4. Find the number of populations (clusters) K using the Tracy-Widom
   statistics.
5. Perform agglomerative clustering using K and p=1 distances.
6. Find the matrix of distances between centers of mass of every cluster.
7. Use neighbour join and an outgroup to find the tree topology from
   the above matrix.
8. Use the tree to construct the constraints and equations for split times D.
9. Find the split times, repeat steps 5-9 to bootstrap
   the confidence interval.


Most of the step above are performed if an option is present in the input file
Pre-processing step:
1: -recode
2: -freqs
3: -manual
4: -filter
5 and 6 are related to -hard -miss and -sparse... not so clear what exactly

Actual analysis
1: -asd
2  -> 9: -analyse



the pattern are the file name and should be passed to it with {}
that will be replaced with the number of the chromosome



