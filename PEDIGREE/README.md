# PEDIGREE

All the code is present in this folder:

- 1_Get_F_PEDIGREE.sh: estimate the different inbreeding coefficients (F) used in this study for the simulated PEDIGREE population

- 2_creating_GRMs.sh: create the three GRMs used: Allele-sharing GRM, GCTA weighted GRM and GCTA unweighted

- 3_Simulating_ID.sh: Run inbreeding depression simulations

- 4_PEDIGREE_RANDOM_SUBSAMPLING: subsample PEDIGREE population randomly, estimate GRM matrices and run ID simulations

- 5_PEDIGREE_RANGED_SUBSAMPLING: subsample PEDIGREE population covering the entire spectrum of F values, estimate GRM matrices and run ID simulations

For the supplementary material where we we used intermediate frequencies loci as causal markers, you just need to modify one parameter in 3_Simulating_ID.sh: the maf parameter (on line 65), set it to 0.1


### functions folder

contains the code with the bash and R functions used at different steps of the analyses !


## data folder

contains data necessary to re-run the analyses from the paper, namely the bedmatrix
