# content of the folder

## launch.sh
This is a bash script which first convert the file script/CEU_YRI.chr22.vcf.gz and script/label.lab into Medeas internal format and then run medeas with thess files. This script is easy to adapt for your own data. It requires plink and common UNIX tools to be installed

## Snakefile 
This is a snakemake pipeline which first convert the file script/CEU_YRI.chr22.vcf.gz and script/label.lab into Medeas internal format and then run medeas with these files. It is the recommended option if you have snakemake and conda installed on your system. The input vcf and label file can be changed in script/config.yaml, and the script should be launched using: snakemake --use-conda Using this command, you don’t need to have python 3.6 and plink install, conda will take care of it for you. It can be configured with the config.yaml file.

##CEU_YRI_chr22.vcf.gz
This is a compress version of a vcf file obtained through the public data of the 1000 genome. 
We downloaded the data from the following ftp repository ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/. This repository contain vcf file and a file with the assignation of individuals to population. We subsample the vcf file from the chromosome 22 so that it contains only the CEU and YRI population, and we only kept bi-allelic SNP. 

## label.lab
This file contain the population of the individual contained in the vcf file. 