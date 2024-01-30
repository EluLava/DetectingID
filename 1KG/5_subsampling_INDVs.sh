#!/bin/bash

Running_WORLD_RANGEDSUB(){

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
        cl <- 5
        registerDoParallel(cl)

        source("./functions/functions.R")

        #Read the bed POLYPED (complete) matrix
        bedTOT = read.bed.matrix("./data/bedmatrices/WORLD_ch_1_22")
        #Read the associated F table
        FTOT = read.table("./data/F_datasets/WORLD_allchr_allFs.txt", header = T)
        #rm F estimated with continental all frq
        FTOT = FTOT[,c(1:3,5,7,9,11,13,15,17)]

        #read smple description file
        SampDesc = read.table("../1KG_WORLD/data/SAMPLE_DESC/integrated_call_samples_v3.20130502.ALL.panel", h = T)

        #Extract list of samples range F1: 11,924 samples
        samples.1 = FTOT[FTOT[,6] <= 0.1, c(1,3)]

        #divide by continent
        samples.1.EUR = samples.1[samples.1[,2] == "EUR",]
        samples.1.EAS = samples.1[samples.1[,2] == "EAS",]
        samples.1.AMR = samples.1[samples.1[,2] == "AMR",]
        samples.1.SAS = samples.1[samples.1[,2] == "SAS",]
        samples.1.AFR = samples.1[samples.1[,2] == "AFR",]

        #Extract list of samples range F2: 143 samples
        samples.2 = FTOT[FTOT[,6] > 0.1 & FTOT[,6] <= 0.2, c(1,3)]

        #divide by continent
        samples.2.EUR = samples.2[samples.2[,2] == "EUR",]
        samples.2.EAS = samples.2[samples.2[,2] == "EAS",]
        samples.2.AMR = samples.2[samples.2[,2] == "AMR",]
        samples.2.SAS = samples.2[samples.2[,2] == "SAS",]
        samples.2.AFR = samples.2[samples.2[,2] == "AFR",]

        #Extract list of samples range F3: 18 samples
        samples.3 = FTOT[FTOT[,6] > 0.2 & FTOT[,6] <= 0.3, c(1,3)]

        #Extract nb INDVs to subsample
        nbindvs = as.numeric("${SampSize}")

        #Loop through REPLICATES
        foreachoutput = foreach(repl=1:100, .combine=rbind, .packages=c("gaston","hierfstat","SNPRelate","JGTeach","data.table")) %dopar% {

            print(paste0("Starting replicate ", repl, " !!!"))

            #nb idv to sumsabple from F1, F2 and F3: numbers depends on how many INDVs we want at the end = (nbindvs we want - the mandatory inbred samples)/2 (because F1 and F2)
            nbindvsSUB = round((nbindvs - (0 + 0 + 8))/2)

            if(nbindvsSUB > nrow(samples.2.SAS)){

                nbsub2 = (nbindvsSUB - nrow(samples.2.SAS)) #+ (nbindvsSUB - nrow(samples.3))

                if(nbsub2/3 > nrow(samples.2.AMR)){

                    ##RANDOMLY subsample INDVs range F1
                    dta_samples.1.1 = sort(sample(x = samples.1.EUR[,1], size = nbindvsSUB/5, replace = F))
                    dta_samples.1.2 = sort(sample(x = samples.1.EAS[,1], size = nbindvsSUB/5, replace = F))
                    dta_samples.1.3 = sort(sample(x = samples.1.SAS[,1], size = nbindvsSUB/5, replace = F))
                    dta_samples.1.4 = sort(sample(x = samples.1.AMR[,1], size = nbindvsSUB/5, replace = F))
                    dta_samples.1.5 = sort(sample(x = samples.1.AFR[,1], size = nbindvsSUB/5, replace = F))

                    nbsub3 = nbsub2 - nrow(samples.2.SAS) - nrow(samples.2.AMR) - nrow(samples.2.EAS)

                    #RANDOMLY subsample INDVs range F2
                    dta_samples.2.1 = sort(sample(x = samples.2.EUR[,1], size = 1, replace = F))
                    dta_samples.2.2 = samples.2.SAS[,1]
                    dta_samples.2.3 = samples.2.EAS[,1]
                    dta_samples.2.4 = samples.2.AMR[,1]
                    dta_samples.2.5 = sort(sample(x = samples.2.AFR[,1], size = nbsub3, replace = F))

                    #We'll take all of these range F3
                    dta_samples.3 = samples.3[,1]

                } else {

                    #RANDOMLY subsample INDVs range F1
                    dta_samples.1.1 = sort(sample(x = samples.1.EUR[,1], size = nbindvsSUB/5, replace = F))
                    dta_samples.1.2 = sort(sample(x = samples.1.EAS[,1], size = nbindvsSUB/5, replace = F))
                    dta_samples.1.3 = sort(sample(x = samples.1.SAS[,1], size = nbindvsSUB/5, replace = F))
                    dta_samples.1.4 = sort(sample(x = samples.1.AMR[,1], size = nbindvsSUB/5, replace = F))
                    dta_samples.1.5 = sort(sample(x = samples.1.AFR[,1], size = nbindvsSUB/5, replace = F))

                    #RANDOMLY subsample INDVs range F2
                    dta_samples.2.1 = sort(sample(x = samples.2.EUR[,1], size = 1, replace = F))
                    dta_samples.2.2 = samples.2.SAS[,1]
                    dta_samples.2.3 = sort(sample(x = samples.2.EAS[,1], size = nbsub2/3, replace = F))
                    dta_samples.2.4 = sort(sample(x = samples.2.AMR[,1], size = nbsub2/3, replace = F))
                    dta_samples.2.5 = sort(sample(x = samples.2.AFR[,1], size = nbsub2/3, replace = F))

                    #We'll take all of these range F3
                    dta_samples.3 = samples.3[,1]
                    #We'll take all of these range F4
                    #dta_samples.4 = samples.4[,1]
                    #We'll take all of these range F5

                #dta_samples.5 = samples.5[,1]
                }

            } else {

                #RANDOMLY subsample INDVs range F1
                dta_samples.1.1 = sort(sample(x = samples.1.EUR[,1], size = nbindvsSUB/5, replace = F))
                dta_samples.1.2 = sort(sample(x = samples.1.EAS[,1], size = nbindvsSUB/5, replace = F))
                dta_samples.1.3 = sort(sample(x = samples.1.SAS[,1], size = nbindvsSUB/5, replace = F))
                dta_samples.1.4 = sort(sample(x = samples.1.AMR[,1], size = nbindvsSUB/5, replace = F))
                dta_samples.1.5 = sort(sample(x = samples.1.AFR[,1], size = nbindvsSUB/5, replace = F))

                #RANDOMLY subsample INDVs range F2
                dta_samples.2.1 = sort(sample(x = samples.2.EUR[,1], size = 1, replace = F))
                dta_samples.2.2 = sort(sample(x = samples.2.EAS[,1], size = nbindvsSUB/4, replace = F))
                dta_samples.2.3 = sort(sample(x = samples.2.SAS[,1], size = nbindvsSUB/4, replace = F))
                dta_samples.2.4 = sort(sample(x = samples.2.AMR[,1], size = nbindvsSUB/4, replace = F))
                dta_samples.2.5 = sort(sample(x = samples.2.AFR[,1], size = nbindvsSUB/4, replace = F))

                #We'll take all of these range F3
                dta_samples.3 = samples.3[,1]

            }

            #merge list of indvs to subsample
            dta_samples = sort(c(dta_samples.1.1,dta_samples.1.2,dta_samples.1.3,dta_samples.1.4,dta_samples.1.5, dta_samples.2.1,dta_samples.2.2,dta_samples.2.3,dta_samples.2.4,dta_samples.2.5, dta_samples.3))

            #check if we have correct individuals
            if(length(dta_samples) < nbindvs){

                #estimate diff
                diff = (nbindvs - length(dta_samples))
                #sample random X more indivs
                addindv = sample(x = FTOT[!(FTOT[,1] %in% dta_samples),1], size = diff, replace = F)
                dta_samples = c(dta_samples, addindv)

            } else if (length(dta_samples) > nbindvs){

                diff = (length(dta_samples) - nbindvs)
                rmind = sample(dta_samples, size = diff, replace = F)
                dta_samples = dta_samples[!(dta_samples %in% rmind)]
            }

            #subsample bed matrix
            bedSample = select.inds(bedTOT, id %in% dta_samples)

            #filter for maf > 0 --> non monomorphic SNPs
            bedSample.snps = select.snps(bedSample, maf > 0)

            #write bedmatrix for MatchMat later on
            write.bed.matrix(bedSample.snps, paste0("./data/bedmatrices/rangedWORLD_ch_1_22_${SampSize}_",repl))

            #rm for mem
            rm(bedSample)

            #Extract corresponding F table
            FSample = FTOT[FTOT[,1] %in% dta_samples,]

            #Create corresponding bedmatrices
            #From bed to GDS
            SNPRelate::snpgdsBED2GDS(bed.fn = paste0("./data/bedmatrices/rangedWORLD_ch_1_22_${SampSize}_",repl,".bed"), fam.fn = paste0("./data/bedmatrices/rangedWORLD_ch_1_22_${SampSize}_",repl,".fam"),
            bim.fn = paste0("./data/bedmatrices/rangedWORLD_ch_1_22_${SampSize}_",repl,".bim"), out.gdsfn = paste0("./data/bedmatrices/rangedWORLD_ch_1_22_${SampSize}_",repl,".gds"))

            #Open GDS object
            gdsMat = SNPRelate::snpgdsOpen(paste0("./data/bedmatrices/rangedWORLD_ch_1_22_${SampSize}_",repl,".gds"))

            #From GDS to GRM gcta
            gcta = snpgdsGRM(gdsMat, method = "GCTA")

            #From GDS to dosage
            dos = snpgdsGetGeno(gdsMat)

            #matching Mat
            AS = matching(dos)

            #beta matching mat
            kas = beta.dosage(AS, MATCHING=TRUE)
            #pass to GRM
            grmas = kinship2grm(kas)

            #pairwise kinships (STANDARD) ratio of averages
            kc0 = Kc0(AS,matching=TRUE)

            #rm useless objects (for mem)
            rm(gdsMat); rm(dos); rm(AS); rm(kas)
            #rm gds and bed matrix from this subsampling
            system(paste0("rm ./bedmatrices/rangedWORLD_ch_1_22_${SampSize}_",repl,".*"))

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

            WORLD_LMM_ID3_selhT_selpT_demaT = id.est(FSample[,4:ncol(FSample)], genos = bedSample.snps, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100, selh = TRUE, selp = TRUE, direct = TRUE, b1=b1, b2=b2)

            #Extract Simulated scenarios
            simuScen = ls(pattern = paste0("WORLD_LMM_ID",b,"_"))

            ### linear mixed effect model ID estimates (b) with allele sharing GRM

            for(sim in simuScen){

                #Loop through FSample
                for(Finb in colnames(FSample)[4:ncol(FSample)]){

                    #print
                    writeLines(paste0("\n", sim, " ", Finb, "\n"))

                    #assign to respective object with GRM AS Matrix
                    assign(x = paste0(Finb, "_", sim, "_MatGRMAS"), value = mclapply(1:100, function(z) {print(z)
                    lmm.aireml(Y = get(sim)[[2]][,z], K = grmas, X = cbind(1,FSample[,colnames(FSample) == Finb]), verbose = T)[[9]][2]}
                    , mc.cores = 1))

                    #assign to respective object with GRM Matrix kc0
                    assign(x = paste0(Finb, "_", sim, "_MatGRMkc0"), value = mclapply(1:100, function(z) {print(z)
                    lmm.aireml(Y = get(sim)[[2]][,z], K = 2*kc0, X = cbind(1,FSample[,colnames(FSample) == Finb]), verbose=F)[[9]][2]}
                    , mc.cores = 1))

                    #assign to respective object with GRM Matrix gcta
                    assign(x = paste0(Finb, "_", sim, "_MatGRMgcta"), value = mclapply(1:100, function(z) {print(z)
                    lmm.aireml(Y = get(sim)[[2]][,z], K = gcta[[4]], X = cbind(1,FSample[,colnames(FSample) == Finb]), verbose=F)[[9]][2]}
                    , mc.cores = 1))
                }

                assign(x = paste0("ID_EST_", sim), value = data.frame(Simuscen = rep(sim, nrow(get(sim)[[6]])),
                REP = rep(repl, nrow(get(sim)[[6]])),
                get(sim)[[6]][,c(1:13)],
                #Fped.lme.AS = unlist(get(paste0("Fped_", sim, "_MatGRMAS"))),
                Fas.lme.AS = unlist(get(paste0("Fas.w_", sim, "_MatGRMAS"))),
                #Fhom.lme.AS = unlist(get(paste0("Fhom_", sim, "_MatGRMAS"))),
                FuniUN.lme.AS = unlist(get(paste0("Funi.un.w_", sim, "_MatGRMAS"))),
                FuniWE.lme.AS = unlist(get(paste0("Funi.we.w_", sim, "_MatGRMAS"))),
                Froh.BCFTools.100KB.lme.AS = unlist(get(paste0("FROHs.BCFTools.100KB.w_", sim, "_MatGRMAS"))),
                Froh.PLINK.100KB.lme.AS = unlist(get(paste0("FROHs.PLINK.100KB.w_", sim, "_MatGRMAS"))),
                Froh.BCFTools.1MB.lme.AS = unlist(get(paste0("FROHs.BCFTools.1MB.w_", sim, "_MatGRMAS"))),
                Froh.PLINK.1MB.lme.AS = unlist(get(paste0("FROHs.PLINK.1MB.w_", sim, "_MatGRMAS"))),
                #Fped.lme.kc0 = unlist(get(paste0("Fped_", sim, "_MatGRMkc0"))),
                Fas.lme.kc0 = unlist(get(paste0("Fas.w_", sim, "_MatGRMkc0"))),
                #Fhom.lme.kc0 = unlist(get(paste0("Fhom_", sim, "_MatGRMkc0"))),
                FuniUN.lme.kc0 = unlist(get(paste0("Funi.un.w_", sim, "_MatGRMkc0"))),
                FuniWE.lme.kc0 = unlist(get(paste0("Funi.we.w_", sim, "_MatGRMkc0"))),
                Froh.BCFTools.100KB.lme.kc0 = unlist(get(paste0("FROHs.BCFTools.100KB.w_", sim, "_MatGRMkc0"))),
                Froh.PLINK.100KB.lme.kc0 = unlist(get(paste0("FROHs.PLINK.100KB.w_", sim, "_MatGRMkc0"))),
                Froh.BCFTools.1MB.lme.kc0 = unlist(get(paste0("FROHs.BCFTools.1MB.w_", sim, "_MatGRMkc0"))),
                Froh.PLINK.1MB.lme.kc0 = unlist(get(paste0("FROHs.PLINK.1MB.w_", sim, "_MatGRMkc0"))),
                #Fped.lme.gcta = unlist(get(paste0("Fped_", sim, "_MatGRMgcta"))),
                Fas.lme.gcta = unlist(get(paste0("Fas.w_", sim, "_MatGRMgcta"))),
                #Fhom.lme.gcta = unlist(get(paste0("Fhom.w_", sim, "_MatGRMgcta"))),
                FuniUN.lme.gcta = unlist(get(paste0("Funi.un.w_", sim, "_MatGRMgcta"))),
                FuniWE.lme.gcta = unlist(get(paste0("Funi.we.w_", sim, "_MatGRMgcta"))),
                Froh.BCFTools.100KB.lme.gcta = unlist(get(paste0("FROHs.BCFTools.100KB.w_", sim, "_MatGRMgcta"))),
                Froh.PLINK.100KB.lme.gcta = unlist(get(paste0("FROHs.PLINK.100KB.w_", sim, "_MatGRMgcta"))),
                Froh.BCFTools.1MB.lme.gcta = unlist(get(paste0("FROHs.BCFTools.1MB.w_", sim, "_MatGRMgcta"))),
                Froh.PLINK.1MB.lme.gcta = unlist(get(paste0("FROHs.PLINK.1MB.w_", sim, "_MatGRMgcta")))))

            }

            ### merge all scenarios into one df for output

            SimuScen=ls(pattern = "ID_EST_")

            #create emoty df for output
            ID_EST_allscenarios = as.data.frame(matrix(ncol = ncol(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT), nrow = 0))
            colnames(ID_EST_allscenarios) = colnames(ID_EST_WORLD_LMM_ID3_selhT_selpT_demaT)

            #Loop through scenarios
            for(sim in SimuScen){

                ID_EST_allscenarios = rbind(ID_EST_allscenarios, get(sim))

            }

            #return the output
            return(ID_EST_allscenarios)

        }

        #Then save the R env.
        save(foreachoutput, file="./data/ID_EST_Datasets/${SampSize}_SUB_WORLD_ID3_simus_COMPLETE_ID_EST_noFrecalc.RData")

EOF

}

export -f Running_WORLD_RANGEDSUB

#Extract runFILE from inputFILE: one array per superpop
runFILE=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ./data/inputFILES/WORLD_SUB_sampSize.list)

Running_WORLD_RANGEDSUB ${runFILE}
