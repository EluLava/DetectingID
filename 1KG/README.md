# 1KG (WORLD, EAS and AFR)

All the code is present in this folder:

- 1_Get_F_1KG.sh: estimate the different inbreeding coefficients (F) used in this study for the three 1KG populations

- 2_creating_GRMs.sh: create the three GRMs used: Allele-sharing GRM, GCTA weighted GRM and GCTA unweighted

- 3_MAF_LD_Stratification.sh: contains the code to get MAF and LD stratification (LDMS) inbreeding coefficients estimation (not in the main text, presented in supplementary material BUT necessary to run 4_Simulating_ID.sh)

- 4_Simulating_ID.sh: Run inbreeding depression simulations

- 5_subsampling_INDVs.sh: Subsample individuals (50, 100, 250 and 500) (ensuring we have structure) from the WORLD dataset and resimulate ID and quantify ID with different F.


For the supplementary material where we filtered loci on MAF < 0.05, you can simply filter the bed matrix on MAF first and then rerun steps 1 to 4 with the new bed file.


For the supplementary material where we we used all segments to estimate FHBD.QUAL, you just need to modify the BCFTools part in 1_Get_F_1KG.sh and take all segments (rather that filtering on size).


For the supplementary material where we we used intermediate frequencies loci as causal markers, you just need to modify one parameter in 4_Simulating_ID.sh: the maf parameter (on line 86), set it to 0.1


### functions folder

contains the code with the bash and R functions used at different steps of the analyses !


## data folder

contains data necessary to re-run the analyses from the paper
