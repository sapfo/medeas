# Test dataseets

The folder contains two test datasets in different formats.

## Simulation data (simulation.*) ##
This dataset was generated using a coalescence simulation software (scrm). The data was then post-process to match medeas format.

We simulated three populations of 25 individuals with each 10'000 Loci. All the population sizes are constant. The first split happen at D = 0.2 and the second one at D = 0.4 (in unit of N) (which in unit of scrm is D = 0.1 and D = 0.2 since they work in unit of 2N).

### Files ###
- MEDEAS
    - simulation.snps
- PLINK
    - simulation.bed
    - simulation.bim
    - simulation.fam
- VCF
    - simulation.vcf.gz
    - simulation.vcf.gz.tbi
- Labels
    - simulation.label

## Real data (CEU_YRI.chr22.*) ##
This dataset contains real data from the 1000GP. Extracted are two populations CEU adn YRI.

### Files ###
- MEDEAS
    - CEU_YRI.chr22.snps
- PLINK
    - CEU_YRI.chr22.bed
    - CEU_YRI.chr22.bim
    - CEU_YRI.chr22.fam
- VCF (phased)
    - CEU_YRI.chr22.vcf.gz
    - CEU_YRI.chr22.vcf.gz.tbi
- Labels
    - CEU_YRI.chr22.label