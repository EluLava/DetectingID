# PEDIGREE

All the code is present in this folder:

- 0_Simulating_PED.sh: contains the code to simulate the pedigree and it's corresponding genotypes (as bed matrix)

- 1_Get_F_PEDIGREE.sh: estimate the different inbreeding coefficients (F) used in this study for the simulated PEDIGREE population

- 2_creating_GRMs.sh: create the three GRMs used: Allele-sharing GRM, GCTA weighted GRM and GCTA unweighted

- 3_Simulating_ID.sh: Run inbreeding depression simulations

- 4_PEDIGREE_RANDOM_SUBSAMPLING: subsample PEDIGREE population randomly (to 1KG data set sizes: EAS, AFR and WORLD), estimate GRM matrices and run ID simulations

- 5_PEDIGREE_RANGED_SUBSAMPLING: subsample PEDIGREE population covering the entire spectrum of F values (to 1KG data set sizes: EAS, AFR and WORLD), estimate GRM matrices and run ID simulations

- 6.1_PEDIGREE_SMALLPopSize_simPED.sh: simulate PEDIGREE populations with small sample sizes (50, 100, 250 and 500), calculate Fs and GRMs

- 6.2_PEDIGREE_SMALLPopSize_simID.sh: merge F dataset and simulate inbreeding depression on small pedigrees (simulated in 6.1)

For the supplementary material where we we used intermediate frequencies loci as causal markers, you just need to modify one parameter in 3_Simulating_ID.sh: the maf parameter (on line 65), set it to 0.1


### functions folder

contains the code with the bash and R functions used at different steps of the analyses !


## data folder

contains data necessary to re-run the analyses from the paper, namely the bedmatrix
