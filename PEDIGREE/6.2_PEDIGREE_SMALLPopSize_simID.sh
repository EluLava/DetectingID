#! /bin/bash

mkdir -p ./data/ID_EST_Datasets/NewSIMUs/

#Define function
Running_POLYPED_RANGEDSUB(){

    #First argument passed to the script is continent (for # indvs subsampling)
    SampSize=$1

    #Open R
    R --vanilla <<EOF

        #Load packages
        library(gaston)
        library(hierfstat)
        library(JGTeach)
        library(parallel)
        library(doParallel)
        library(data.table)
        library(SNPRelate)
        cl <- 3
        registerDoParallel(cl)

        source("../1KG/functions/functions.R")

        #Extract nb INDVs to subsample
        nbindvs = as.numeric("${SampSize}")

        #Read the bed matrix
        bed = read.bed.matrix("./data/bedmatrices/BED_ReducedPGPED_${SampSize}")

        #Read F tables (and merge)
        dtaFSNPs = read.table("./data/F_datasets/ReducedPGPED_${SampSize}_SNPsbasedF.txt", h = T)
        colnames(dtaFSNPs)[1] = c("sample")
        dtaFBCFTools = read.table("./data/F_datasets/CallingROHsBCFTools/ReducedPGPED_${SampSize}.txt", h = T)
        dtaF.1 = merge(dtaFSNPs, dtaFBCFTools, by = "sample")
        dtaFPLINK = read.table("./data/F_datasets/CallingROHsPLINK/ReducedPGPED_${SampSize}.txt", h = T)
        dtaF.2 = merge(dtaF.1, dtaFPLINK, by = "sample")
        FSample = as.data.frame(cbind(sample = dtaF.1[,1], pop = rep(NA, nrow(dtaF.1)), cont = rep(NA, nrow(dtaF.1)), Fped = dtaF.2[,2], Fas = dtaF.2[,3], Funi.un = dtaF.2[,4], Funi.we = dtaF.2[,5], FROHs.BCFTools.100KB = dtaF.2[,6], FROHs.PLINK.100KB = dtaF.2[,8],
	                                  FROHs.BCFTools.1MB = dtaF.2[,7], FROHs.PLINK.1MB = dtaF.2[,9]))

        #Read bed matrices
        grmas = readRDS("./data/GRMs/ReducedPGPED_${SampSize}_AS_GRM.RDS")
        kc0 = readRDS("./data/GRMs/ReducedPGPED_${SampSize}_GCTAwe_GRM.RDS")
        gcta = readRDS("./data/GRMs/ReducedPGPED_${SampSize}_GCTAun_GRM.RDS")

        ### Set SIM parameters

        #Nb of causal loci
        ncaus = 1000
        #Filtering on MAF = 0.001 to avoid super rare alleles causing super strong effects
        maf = 0.01
        #ID (and dominance effect) = 3
        b = 3
        #Additive effects = 10
        s = 10
        #beta dist parameters, 1 & 1 = uniform distribution
        b1 = 1
        b2 = 1

        ### SIMULATION
        
        PGPED_LMM_ID3_selhT_selpT_demaT = id.est(FSample[,4:ncol(FSample)], genos = bed, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100, selh = TRUE, selp = TRUE, direct = TRUE, b1=b1, b2=b2)

        #Extract Simulated scenarios
        simuScen = ls(pattern = paste0("PGPED_LMM_ID",b,"_"))

        ### linear mixed effect model ID estimates (b) with allele sharing GRM

        for(sim in simuScen){

            #Loop through FSample
            for(Finb in colnames(FSample)[4:ncol(FSample)]){

                #print
                writeLines(paste0("\n", sim, " ", Finb, "\n"))

                #assign to respective object with GRM AS Matrix
                assign(x = paste0(Finb, "_", sim, "_MatGRMAS"), value = mclapply(1:100, function(z) {print(z)
                lmm.aireml(Y = get(sim)[[2]][,z], K = grmas, X = cbind(1,FSample[,colnames(FSample) == Finb]), verbose = T)[[9]][2]}
                , mc.cores = 3))

                #assign to respective object with GRM Matrix kc0
                assign(x = paste0(Finb, "_", sim, "_MatGRMkc0"), value = mclapply(1:100, function(z) {print(z)
                lmm.aireml(Y = get(sim)[[2]][,z], K = 2*kc0, X = cbind(1,FSample[,colnames(FSample) == Finb]), verbose=F)[[9]][2]}
                , mc.cores = 3))

                #assign to respective object with GRM Matrix gcta
                assign(x = paste0(Finb, "_", sim, "_MatGRMgcta"), value = mclapply(1:100, function(z) {print(z)
                lmm.aireml(Y = get(sim)[[2]][,z], K = gcta[[4]], X = cbind(1,FSample[,colnames(FSample) == Finb]), verbose=F)[[9]][2]}
                , mc.cores = 3))
            }

            ### After estimating all the ID est for this scenario, bringing everything together in a R object

            assign(x = paste0("ID_EST_", sim), value = data.frame(get(sim)[[6]][,c(1:12)],
            Fped.lme.AS = unlist(get(paste0("Fped_", sim, "_MatGRMAS"))),
            Fas.lme.AS = unlist(get(paste0("Fas_", sim, "_MatGRMAS"))),
            FuniUN.lme.AS = unlist(get(paste0("Funi.un_", sim, "_MatGRMAS"))),
            FuniWE.lme.AS = unlist(get(paste0("Funi.we_", sim, "_MatGRMAS"))),
            Froh.BCFTools.100KB.lme.AS = unlist(get(paste0("FROHs.BCFTools.100KB_", sim, "_MatGRMAS"))),
            Froh.PLINK.100KB.lme.AS = unlist(get(paste0("FROHs.PLINK.100KB_", sim, "_MatGRMAS"))),
            Froh.BCFTools.1MB.lme.AS = unlist(get(paste0("FROHs.BCFTools.1MB_", sim, "_MatGRMAS"))),
            Froh.PLINK.1MB.lme.AS = unlist(get(paste0("FROHs.PLINK.1MB_", sim, "_MatGRMAS"))),
            Fped.lme.kc0 = unlist(get(paste0("Fped_", sim, "_MatGRMkc0"))),
            Fas.lme.kc0 = unlist(get(paste0("Fas_", sim, "_MatGRMkc0"))),
            FuniUN.lme.kc0 = unlist(get(paste0("Funi.un_", sim, "_MatGRMkc0"))),
            FuniWE.lme.kc0 = unlist(get(paste0("Funi.we_", sim, "_MatGRMkc0"))),
            Froh.BCFTools.100KB.lme.kc0 = unlist(get(paste0("FROHs.BCFTools.100KB_", sim, "_MatGRMkc0"))),
            Froh.PLINK.100KB.lme.kc0 = unlist(get(paste0("FROHs.PLINK.100KB_", sim, "_MatGRMkc0"))),
            Froh.BCFTools.1MB.lme.kc0 = unlist(get(paste0("FROHs.BCFTools.1MB_", sim, "_MatGRMkc0"))),
            Froh.PLINK.1MB.lme.kc0 = unlist(get(paste0("FROHs.PLINK.1MB_", sim, "_MatGRMkc0"))),
            Fped.lme.gcta = unlist(get(paste0("Fped_", sim, "_MatGRMgcta"))),
            Fas.lme.gcta = unlist(get(paste0("Fas_", sim, "_MatGRMgcta"))),
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
        ID_EST_allscenarios = as.data.frame(matrix(ncol = ncol(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT), nrow = 0))
        colnames(ID_EST_allscenarios) = colnames(ID_EST_PGPED_LMM_ID3_selhT_selpT_demaT)

        #Loop through scenarios
        for(sim in SimuScen){

            ID_EST_allscenarios = rbind(ID_EST_allscenarios, get(sim))

        }


        #Then save the R env.
        save(ID_EST_allscenarios, file="./data/ID_EST_Datasets/NewSIMUs/${SampSize}_SUB_POLYPED_ID3_simus_COMPLETE_ID_EST.RData")

EOF

}

export -f Running_POLYPED_RANGEDSUB

#Extract runFILE from inputFILE: one array per superpop
runFILE=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ./data/inputFILES/POLYPED_SUB_sampSize.list)

Running_POLYPED_RANGEDSUB ${runFILE}
