======
medeas
======
listing of directories, 02.01.18
================================
main.py: main entry point
names for temp files.

recode: add spaces between alleles (numpy likes it better)
TODO change the number of chromosomes to argparse

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

!! currently always overriding with K=2 hard coded
TODO: if SIMULATION: for K=2"


TODO: codes to strings in hardfilt
TODO: recalculate everything at every bootstrap run
TODO: only savetxt if SIMULATION
TODO: output/saving hard coded to 2 populations


the user cannot use the p currently (marchenko pastur p=2 ONLY, for clustering its unclear but currently only p=1)

TODO: chang the output name to something else Papuan
=============
options.py

important file, file through which modules communicate.
===========
processing
filtering and recoding scripts and functions in that folder.
===========
simulate.py: script that runs many other scripts
==========
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
TODO: remove the plotting logic from this file.
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

for now, n1, n2, theta: hard coded.
=================================
installation note:
sapfo@electra:~/Projects/medeas$ python3.6 -m pip install scikit-bio
sapfo@electra:~/Projects/medeas$ sudo python3.6 -m pip install numpy
