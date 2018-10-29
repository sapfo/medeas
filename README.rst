======
medeas
======
listing of directories, 29.10.18
================================

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

If there are various chromosome, the pattern are the file name and should be passed to it with {}
that will be replaced with the number of the chromosome



TODO: recalculate everything at every bootstrap run


the user cannot use the p currently (marchenko pastur p=2 ONLY, for clustering its unclear but currently only p=2)
For diploid data (Ana's owl) it seems that p=2 do a better job at clustering. 


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

3 main functions. first important function. perform clustering. every individual is assigned a label.
perform_clustering: where we assign the labels
find_tree: where we find tree
find_distances: most complex and important function. constrains. etc.
==========




