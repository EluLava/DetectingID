#Load the functions
source ./functions/ALL_bash_functions.sh

#Create output directories if do not exist
mkdir -p ./data/LDScores/
mkdir -p ./data/LD_MAF_bins

#Loop through populations
for CONT in {WORLD, EAS, AFR}; do

    #Create LD scores file
    LD_scores ./data/bedmatrices/${CONT}_ch_1_22 ./data/LDScores/${CONT}_LDscores_ch_1_22.score.ld

    #Run LDMD analysis
    Stratification_MAF_LD ./data/bedmatrices/${CONT}_ch_1_22 ./data/LDScores/${CONT}_LDscores_ch_1_22.score.ld \
    ./data/LD_MAF_bins/${CONT}_LD_MAF_F_per_bins

done