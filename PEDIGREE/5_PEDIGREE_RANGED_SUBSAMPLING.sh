#! /bin/bash

#SBATCH --mail-user eleonore.lavanchy@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --array 1
#SBATCH --time 03-00:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 40
#SBATCH --mem 500G
#SBATCH --job-name POLYPEDrandomSUB
#SBATCH -o STD/%x_%A_%a.stdout
#SBTACH -e STD/%x_%A_%a.stderr
#SBATCH --export=NONE

#load modules
module load gcc r

#Define function
SIMULATING_ID_RANGED_SUBSAMPLING(){

    #First argument passed to the script is continent (for # indvs subsampling)
    pop=$1

    #Open R
    R --vanilla <<EOF

    list.of.packages <- c("gaston", "Rcpp", "hierfstat", "parallel", "data.table", "devtools")
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)){install.packages(new.packages)}
    if(!("JGTeach" %in% installed.packages()[,"Package"])){library(devtools); install_github("jgx65/JGTeach")}

    #Load packages
    library(gaston)
    library(hierfstat)
    library(JGTeach)
    library(parallel)
    library(doParallel)
    library(data.table)
    library(SNPRelate)
    cl <- 80
    registerDoParallel(cl)

    source("./functions/functions.R")

    #Read the bed POLYPED (complete) matrix
    bedTOT = read.bed.matrix("./data/bedmatrices/PGPED_ch_1_2")
    #Read the associated F table
    FTOT = read.table("./data/F_datasets/PGPED_PED_allchr_allFs.txt", header = T)

    #Extract list of samples range F1: 11,924 samples
    samples.1 = FTOT[FTOT[,8] <= 0.1, 1]

    #Extract list of samples range F2: 143 samples
    samples.2 = FTOT[FTOT[,8] > 0.1 & FTOT[,8] <= 0.2, 1]

    #Extract list of samples range F3: 18 samples
    samples.3 = FTOT[FTOT[,8] > 0.2 & FTOT[,8] <= 0.3, 1]

    #Extract list of samples range F4: 3 samples
    samples.4 = FTOT[FTOT[,8] > 0.3 & FTOT[,8] <= 0.4, 1]

    #Extract list of samples range F5: 1 sample
    samples.5 = FTOT[FTOT[,8] > 0.4, 1]

    #Populations to subsample
    populations = c(2504,504,661)
    names(populations) = c("WORLD","EAS","AFR")

    #Extract "POP" subsample name
    nbindvs = populations[names(populations) == "${pop}"]

    #Loop through REPLICATES
    foreachoutput = foreach(repl=1:50, .combine=rbind) %dopar% {

        #nb idv to sumsabple from F1 and F2: numbers depends on how many INDVs we want at the end = (nbindvs we want - the mandatory inbred samples)/2 (because F1 and F2)
        #You'll need to adjust these numbers if you have a different pedigree with different ranges of inbreeding !
        nbindvsSUB = round((nbindvs - (18 + 3 + 1))/2)

        #check if bigger than possible
        if(nbindvsSUB > length(samples.2)){

            #RANDOMLY subsample INDVs range F1 
            dta_samples.1 = sort(sample(x = samples.1, size = (nbindvsSUB + (nbindvsSUB - length(samples.2))), replace = F))

            #RANDOMLY subsample INDVs range F2
            dta_samples.2 = samples.2

            #We'll take all of these range F3
            dta_samples.3 = samples.3

            #We'll take all of these range F4
            dta_samples.4 = samples.4

            #We'll take all of these range F5
            dta_samples.5 = samples.5

        } else {

            #RANDOMLY subsample INDVs range F1 
            dta_samples.1 = sort(sample(x = samples.1, size = nbindvsSUB, replace = F))

            #RANDOMLY subsample INDVs range F2
            dta_samples.2 = sort(sample(x = samples.2, size = nbindvsSUB, replace = F))

            #We'll take all of these range F3
            dta_samples.3 = samples.3

            #We'll take all of these range F4
            dta_samples.4 = samples.4

            #We'll take all of these range F5
            dta_samples.5 = samples.5
        }

        #merge list of indvs to subsample
        dta_samples = sort(c(dta_samples.1, dta_samples.2, dta_samples.3, dta_samples.4, dta_samples.5))

        #subsample bed matrix
        bedSample = select.inds(bedTOT, id %in% dta_samples)

        #filter for maf > 0 --> non monomorphic SNPs
        bedSample.snps = select.snps(bedSample, maf > 0)

        #write bedmatrix for MatchMat later on
        write.bed.matrix(bedSample.snps, paste0("./data/bedmatrices/PGPED_RandomSUB_BEDFILE_${pop}_",repl))

        #rm for mem
        rm(bedSample)

        #Extract corresponding F table
        FSample = FTOT[FTOT[,1] %in% dta_samples,]

        #Create corresponding bedmatrices
        #From bed to GDS
        snpgdsBED2GDS(bed.fn = paste0("./data/bedmatrices/PGPED_RandomSUB_BEDFILE_${pop}_",repl, ".bed"), fam.fn = paste0("./data/bedmatrices/PGPED_RandomSUB_BEDFILE_${pop}_",repl, ".fam")
        bim.fn = paste0("./data/bedmatrices/PGPED_RandomSUB_BEDFILE_${pop}_",repl, ".bim"), out.gdsfn = paste0("./data/bedmatrices/PGPED_RandomSUB_BEDFILE_${pop}_",repl, ".gds"))

        #Open GDS object
        gdsMat = snpgdsOpen(paste0("./data/bedmatrices/PGPED_RandomSUB_BEDFILE_${pop}_",repl, ".gds"))

        #From GDS to GRM gcta
        GCTAun = snpgdsGRM(gdsMat, method = "GCTA")

        #From GDS to dosage
        dos = snpgdsGetGeno(gdsMat)

        #matching Mat
        AS = matching(dos)

        #beta matching mat
        kas = beta.dosage(AS, MATCHING=TRUE)
        #pass to GRM
        grmas = kinship2grm(kas)

        #pairwise kinships (STANDARD) ratio of averages
        GCTAwe = Kc0(AS, matching=TRUE)

        #rm useless objects (for mem)
        rm(gdsMat); rm(dos); rm(AS); rm(kas)
        #rm gds and bed matrix from this subsampling
        system(paste0("./data/bedmatrices/PGPED_RandomSUB_BEDFILE_${pop}_",repl, ".*"))

        ### Set SIM parameters

        #Nb of causal loci
        ncaus = 1000
        #Filtering on MAF = 0.001 to avoid super rare alleles causing super strong effects
        maf = 0.001
        #ID (and dominance effect) = 3
        b = 3
        #Additive effects = 10
        s = 10
        #beta dist parameters, 1 & 1 = uniform distribution
        b1 = 1
        b2 = 1

        ### SIMULATION

        assign(x = paste0("PGPED_LMM_ID",b,"_selhF_selpF_demaF"), value = id.est(FSample[,4:ncol(FSample)], genos = bedSample.snps, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100, selh = FALSE, selp = FALSE, direct = FALSE, b1=b1, b2=b2))
        assign(x = paste0("PGPED_LMM_ID",b,"_selhT_selpF_demaF"), value = id.est(FSample[,4:ncol(FSample)], genos = bedSample.snps, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100, selh = TRUE, selp = FALSE, direct = FALSE, b1=b1, b2=b2))
        assign(x = paste0("PGPED_LMM_ID",b,"_selhF_selpT_demaF"), value = id.est(FSample[,4:ncol(FSample)], genos = bedSample.snps, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100, selh = FALSE, selp = TRUE, direct = FALSE, b1=b1, b2=b2))
        assign(x = paste0("PGPED_LMM_ID",b,"_selhT_selpT_demaF"), value = id.est(FSample[,4:ncol(FSample)], genos = bedSample.snps, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100, selh = TRUE, selp = TRUE, direct = FALSE, b1=b1, b2=b2))
        assign(x = paste0("PGPED_LMM_ID",b,"_selhF_selpF_demaT"), value = id.est(FSample[,4:ncol(FSample)], genos = bedSample.snps, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100, selh = FALSE, selp = FALSE, direct = TRUE, b1=b1, b2=b2))
        assign(x = paste0("PGPED_LMM_ID",b,"_selhT_selpF_demaT"), value = id.est(FSample[,4:ncol(FSample)], genos = bedSample.snps, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100, selh = TRUE, selp = FALSE, direct = TRUE, b1=b1, b2=b2))
        assign(x = paste0("PGPED_LMM_ID",b,"_selhF_selpT_demaT"), value = id.est(FSample[,4:ncol(FSample)], genos = bedSample.snps, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100, selh = FALSE, selp = TRUE, direct = TRUE, b1=b1, b2=b2))
        assign(x = paste0("PGPED_LMM_ID",b,"_selhT_selpT_demaT"), value = id.est(FSample[,4:ncol(FSample)], genos = bedSample.snps, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100, selh = TRUE, selp = TRUE, direct = TRUE, b1=b1, b2=b2))

        #Extract Simulated scenarios
        simuScen = ls(pattern = paste0("PGPED_LMM_ID",b,"_"))

        ### linear mixed effect model ID estimates (b) with allele sharing GRM

        for(sim in simuScen){

            #Loop through FSample
            for(Finb in colnames(FSample)[4:ncol(FSample)]){

                #print
                writeLines(paste0("\n", sim, " ", Finb, "\n"))

                #assign to respective object with GRM AS
                assign(x = paste0(Finb, "_", sim, "_MatGRMAS"), value = mclapply(1:100, function(z) {print(z)
                lmm.aireml(Y = get(sim)[[2]][,z], K = grmas, X = cbind(1,FSample[,colnames(FSample) == Finb]), verbose = T)[[9]][2]}
                , mc.cores = 2))

                #assign to respective object with GRM GCTA weighted
                assign(x = paste0(Finb, "_", sim, "_MatGRMkc0"), value = mclapply(1:100, function(z) {print(z)
                lmm.aireml(Y = get(sim)[[2]][,z], K = 2*GCTAwe, X = cbind(1,FSample[,colnames(FSample) == Finb]), verbose=F)[[9]][2]}
                , mc.cores = 2))

                #assign to respective object with GRM GCTA unweighted
                assign(x = paste0(Finb, "_", sim, "_MatGRMgcta"), value = mclapply(1:100, function(z) {print(z)
                lmm.aireml(Y = get(sim)[[2]][,z], K = GCTAun[[4]], X = cbind(1,FSample[,colnames(FSample) == Finb]), verbose=F)[[9]][2]}
                , mc.cores = 2))
            }

            ### After estimating all the ID est for this scenario, bringing everything together in a R object

            assign(x = paste0("ID_EST_", sim), value = data.frame(Simuscen = rep(sim, nrow(get(sim)[[6]])),
            REP = rep(repl, nrow(get(sim)[[6]])),
            get(sim)[[6]][,c(1:13)],
            Fped.lme.AS = unlist(get(paste0("Fped_", sim, "_MatGRMAS"))),
            Fas.lme.AS = unlist(get(paste0("Fas_", sim, "_MatGRMAS"))),
            Fhom.lme.AS = unlist(get(paste0("Fhom_", sim, "_MatGRMAS"))),
            FuniUN.lme.AS = unlist(get(paste0("Funi.un_", sim, "_MatGRMAS"))),
            FuniWE.lme.AS = unlist(get(paste0("Funi.we_", sim, "_MatGRMAS"))),
            Froh.BCFTools.100KB.lme.AS = unlist(get(paste0("FROHs.BCFTools.100KB_", sim, "_MatGRMAS"))),
            Froh.PLINK.100KB.lme.AS = unlist(get(paste0("FROHs.PLINK.100KB_", sim, "_MatGRMAS"))),
            Froh.BCFTools.1MB.lme.AS = unlist(get(paste0("FROHs.BCFTools.1MB_", sim, "_MatGRMAS"))),
            Froh.PLINK.1MB.lme.AS = unlist(get(paste0("FROHs.PLINK.1MB_", sim, "_MatGRMAS"))),
            Fped.lme.kc0 = unlist(get(paste0("Fped_", sim, "_MatGRMkc0"))),
            Fas.lme.kc0 = unlist(get(paste0("Fas_", sim, "_MatGRMkc0"))),
            Fhom.lme.kc0 = unlist(get(paste0("Fhom_", sim, "_MatGRMkc0"))),
            FuniUN.lme.kc0 = unlist(get(paste0("Funi.un_", sim, "_MatGRMkc0"))),
            FuniWE.lme.kc0 = unlist(get(paste0("Funi.we_", sim, "_MatGRMkc0"))),
            Froh.BCFTools.100KB.lme.kc0 = unlist(get(paste0("FROHs.BCFTools.100KB_", sim, "_MatGRMkc0"))),
            Froh.PLINK.100KB.lme.kc0 = unlist(get(paste0("FROHs.PLINK.100KB_", sim, "_MatGRMkc0"))),
            Froh.BCFTools.1MB.lme.kc0 = unlist(get(paste0("FROHs.BCFTools.1MB_", sim, "_MatGRMkc0"))),
            Froh.PLINK.1MB.lme.kc0 = unlist(get(paste0("FROHs.PLINK.1MB_", sim, "_MatGRMkc0"))),
            Fped.lme.gcta = unlist(get(paste0("Fped_", sim, "_MatGRMgcta"))),
            Fas.lme.gcta = unlist(get(paste0("Fas_", sim, "_MatGRMgcta"))),
            Fhom.lme.gcta = unlist(get(paste0("Fhom_", sim, "_MatGRMgcta"))),
            FuniUN.lme.gcta = unlist(get(paste0("Funi.un_", sim, "_MatGRMgcta"))),
            FuniWE.lme.gcta = unlist(get(paste0("Funi.we_", sim, "_MatGRMgcta"))),
            Froh.BCFTools.100KB.lme.gcta = unlist(get(paste0("FROHs.BCFTools.100KB_", sim, "_MatGRMgcta"))),
            Froh.PLINK.100KB.lme.gcta = unlist(get(paste0("FROHs.PLINK.100KB_", sim, "_MatGRMgcta"))),
            Froh.BCFTools.1MB.lme.gcta = unlist(get(paste0("FROHs.BCFTools.1MB_", sim, "_MatGRMgcta"))),
            Froh.PLINK.1MB.lme.gcta = unlist(get(paste0("FROHs.PLINK.1MB_", sim, "_MatGRMgcta")))))

        }

        ### merge all scenarios into one df for output

        SimuScen=ls(pattern = "ID_EST_")

        #create emoty df for output
        ID_EST_allscenarios = as.data.frame(matrix(ncol = ncol(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaF), nrow = 0))
        colnames(ID_EST_allscenarios) = colnames(ID_EST_PGPED_LMM_ID3_selhF_selpF_demaF)

        #Loop through scenarios
        for(sim in SimuScen){

            ID_EST_allscenarios = rbind(ID_EST_allscenarios, get(sim))

        }

        #return the output
        return(ID_EST_allscenarios)

    }

    #Then save the R object
    save(foreachoutput, file="./data/ID_EST_Datasets/PGPED_${pop}_RANGEDSUB.RData")

EOF

}

export -f SIMULATING_ID_RANGED_SUBSAMPLING

#Loop through continents
for POP in {WORLD, EAS, AFR}; do

    SIMULATING_ID_RANGED_SUBSAMPLING ${POP}

done