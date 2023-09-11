#create outputdir if does not exists
mkdir -p ./data/ID_EST_Datasets/

#loop through continents
for POP {WORLD, EAS, AFR}; do

	#Describing the function/parameters
	SIMULATING_ID(){
	
		#First argument passed to the function is bedfile
		bedF=$1
	
		#Second argument passed to the function is Fdataset
		Fdata=$2
	
		#Third argument passed to the function is GRMas
		GRMas=$3
	
		#Forth argument passed to the function is GRM GCTAwe
		GCTAwe=$4
	
		#Fifth argument passed to the function is GRM GCTAun
		GCTAun=$5
	
		#Sixth argument passed to the function is F per MAF & LD bins for FAS
		FasbinsLD=$6
	
		#Seventh argument passed to the function is F per MAF & LD bins for FUNI WE
		FuniWEbinsLD=$7
	
		#Eigth argument passed to the function is F per MAF & LD bins for FUN UN
		FuniUNbinsLD=$8
	
		#Ninth argument passed to the function is output file
		output=$9
	
		R --vanilla << EOF
	
		list.of.packages <- c("gaston", "Rcpp", "hierfstat", "parallel", "data.table", "devtools")
		new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
		if(length(new.packages)){install.packages(new.packages)}
		if(!("JGTeach" %in% installed.packages()[,"Package"])){library(devtools); install_github("jgx65/JGTeach")}
	
		library(hierfstat)
		library(JGTeach)
		library(gaston)
		library(parallel)
		library(data.table)
	
		#Source the code containing the functions
		source("./functions/functions.R")
	
		#### LOAD the data ####
	
		#Read bed matrix
		bed = read.bed.matrix("${bedF}")
	
		#Read the all F df
		Fs = read.table("${Fdata}", header = T)
	
		#Read GRMas
		GRMas = readRDS("${GRMas}")
	
		#Read GCTAwe
		GCTAwe = readRDS("${GCTAwe}")
	
		#Read GCTAun
		GCTAun = readRDS("${GCTAun}")
	
		#Read the three F per bins F df
		FasbinsLD = read.table("${FasbinsLD}", header = T)
	
		#Read the three F per bins F df
		FuniWEbinsLD = read.table("${FuniWEbinsLD}", header = T)
	
		#Read the three F per bins F df
		FuniUNbinsLD = read.table("${FuniUNbinsLD}", header = T)
	
	
		##### Set SIM parameters ####
	
		#Nb of causal loci
		ncaus = 1000
	
		#Filtering on MAF = 0.001 to avoid super rare alleles causing super strong effects
		maf = 0.001
	
		#ID (and dominance effect) = 3
		b = 3
	
		#Additive effects = 10
		s = 10
	
		#beta dist parameters, 1 & 1 = uniform distribution, FYI we tried changing the distrbution and it has almost no effect !
		b1 = 1
		b2 = 1
	
		### RUN THE SIMULATION LM for all the scenarios (only the last one rpesented in main text, the rest are in supplementary material) ###
	
		assign(x = paste0("${POP}_LMM_ID",b,"_selhF_selpF_demaF"), value = id.est(Fs[,4:ncol(Fs)], genos = bed, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100, binsfas = FasbinsLD, binsfuniu = FuniUNbinsLD,
		binsfuniw = FuniWEbinsLD, selh = FALSE, selp = FALSE, direct = FALSE, b1=b1, b2=b2))
	
		assign(x = paste0("${POP}_LMM_ID",b,"_selhT_selpF_demaF"), value = id.est(Fs[,4:ncol(Fs)], genos = bed, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100, binsfas = FasbinsLD, binsfuniu = FuniUNbinsLD,
		binsfuniw = FuniWEbinsLD, selh = TRUE, selp = FALSE, direct = FALSE, b1=b1, b2=b2))
	
		assign(x = paste0("${POP}_LMM_ID",b,"_selhF_selpT_demaF"), value = id.est(Fs[,4:ncol(Fs)], genos = bed, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100, binsfas = FasbinsLD, binsfuniu = FuniUNbinsLD,
		binsfuniw = FuniWEbinsLD, selh = FALSE, selp = TRUE, direct = FALSE, b1=b1, b2=b2))
	
		assign(x = paste0("${POP}_LMM_ID",b,"_selhT_selpT_demaF"), value = id.est(Fs[,4:ncol(Fs)], genos = bed, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100, binsfas = FasbinsLD, binsfuniu = FuniUNbinsLD,
		binsfuniw = FuniWEbinsLD, selh = TRUE, selp = TRUE, direct = FALSE, b1=b1, b2=b2))
	
		assign(x = paste0("${POP}_LMM_ID",b,"_selhF_selpF_demaT"), value = id.est(Fs[,4:ncol(Fs)], genos = bed, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100, binsfas = FasbinsLD, binsfuniu = FuniUNbinsLD,
		binsfuniw = FuniWEbinsLD, selh = FALSE, selp = FALSE, direct = TRUE, b1=b1, b2=b2))
	
		assign(x = paste0("${POP}_LMM_ID",b,"_selhT_selpF_demaT"), value = id.est(Fs[,4:ncol(Fs)], genos = bed, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100, binsfas = FasbinsLD, binsfuniu = FuniUNbinsLD,
		binsfuniw = FuniWEbinsLD, selh = TRUE, selp = FALSE, direct = TRUE, b1=b1, b2=b2))
	
		assign(x = paste0("${POP}_LMM_ID",b,"_selhF_selpT_demaT"), value = id.est(Fs[,4:ncol(Fs)], genos = bed, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100, binsfas = FasbinsLD, binsfuniu = FuniUNbinsLD,
		binsfuniw = FuniWEbinsLD, selh = FALSE, selp = TRUE, direct = TRUE, b1=b1, b2=b2))
	
		assign(x = paste0("${POP}_LMM_ID",b,"_selhT_selpT_demaT"), value = id.est(Fs[,4:ncol(Fs)], genos = bed, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100, binsfas = FasbinsLD, binsfuniu = FuniUNbinsLD,
		binsfuniw = FuniWEbinsLD, selh = TRUE, selp = TRUE, direct = TRUE, b1=b1, b2=b2))
	
		#Extract Simulated scenarios
		simuScen = ls(pattern = paste0("${POP}_LMM_ID",b,"_"))
	
		### GET LM LDMS ESTIMATES ###
		#per Simulated scenario & F.SNPs (Fas, FuniUN & FuniWE)
	
		#Loop through simulations schemes (selh, selp and DEMA parameters)
		for(sim in simuScen){
	
			#Loop trough SNPs-based inbreeding coefficients (Fas, Funiu, Funiw) LDMS stratification
			for(ldmsF in c("FasbinsLD", "FuniUNbinsLD", "FuniWEbinsLD")) {
	
				#assign LDMS estimates to corresponding object
				assign(x = paste0("LDMS_",sim,"_", ldmsF), value = lapply(1:100,function(z) { print(z);
				sum(summary(lm(as.formula(paste0(sim, "[[2]][,z] ~ ", paste(colnames(get(ldmsF)), collapse= "+"))), data=get(ldmsF)))[[4]][-1,1])}))
			}
		}
	
	
		### linear mixed effect model ID estimates (b) with GRMs
	
	
		#Loop through simulations schemes (selh, selp and DEMA parameters)
		for(sim in simuScen){
	
			#Loop through Fs
			for(Finb in colnames(Fs)[4:ncol(Fs)]){
	
				#To know where we are
				writeLines(paste0("\n", sim, " ", Finb, "\n"))
	
				#assign to respective object with GRM AS Matrix
				assign(x = paste0(Finb, "_", sim, "_MatGRMAS"), value = mclapply(1:100, function(z) {print(z)
				lmm.aireml(Y = get(sim)[[2]][,z], K = GRMas, X = cbind(1,Fs[,colnames(Fs) == Finb]), verbose = T)[[9]][2]}
				, mc.cores = 30))
	
				#assign to respective object with GRM GCTAwe Matrix
				assign(x = paste0(Finb, "_", sim, "_MatGRMkc0"), value = mclapply(1:100, function(z) {print(z)
				lmm.aireml(Y = get(sim)[[2]][,z], K = 2*GCTAwe, X = cbind(1,Fs[,colnames(Fs) == Finb]), verbose=F)[[9]][2]}
				, mc.cores = 30))
	
				#assign to respective object with GRM GCTAun Matrix
				assign(x = paste0(Finb, "_", sim, "_MatGRMgcta"), value = mclapply(1:100, function(z) {print(z)
				lmm.aireml(Y = get(sim)[[2]][,z], K = GCTAun[[4]], X = cbind(1,Fs[,colnames(Fs) == Finb]), verbose=F)[[9]][2]}
				, mc.cores = 30))
	
				#If bedmatrix is for the WORLD we take Fas.w, FuniUN.w & FuniWE.w
				if("${POP}" == "WORLD"){

					#IF Finb is a F.SNPs coeff (Fas, Funiu, Funiw) do LDMS strat
					if(Finb %in% c("Fas.w", "Funi.un.w", "Funi.we.w")){
	
						#Corresponing matrix LDMS
						if(Finb == "Fas.w"){
						ldmsMat = FasbinsLD
	
						} else if (Finb == "Funi.un.w"){
						ldmsMat = FuniUNbinsLD
	
						} else if (Finb == "Funi.we.w"){
						ldmsMat = FuniWEbinsLD
	
						}
	
						#LDMS F with GRM AS
						assign(x = paste0(Finb, "_", sim, "_MatGRMAS_LDMS"), value = mclapply(1:100, function(z) {print(z)
						lmm.aireml(Y = get(sim)[[2]][,z], K = GRMas, X = cbind(1,as.matrix(ldmsMat)), verbose = T)[[9]][2]}
						, mc.cores = 30))
	
						#LDMS F with GRM GCTAwe
						assign(x = paste0(Finb, "_", sim, "_MatGRMkc0_LDMS"), value = mclapply(1:100, function(z) {print(z)
						lmm.aireml(Y = get(sim)[[2]][,z], K = 2*GCTAwe, X = cbind(1,as.matrix(ldmsMat)), verbose = T)[[9]][2]}
						, mc.cores = 30))
	
						#LDMS F with GRM GCTAun
						assign(x = paste0(Finb, "_", sim, "_MatGRMgcta_LDMS"), value = mclapply(1:100, function(z) {print(z)
						lmm.aireml(Y = get(sim)[[2]][,z], K = GCTAun[[4]], X = cbind(1,as.matrix(ldmsMat)), verbose = T)[[9]][2]}
						, mc.cores = 30))
					}
	
				#Else we take fas.c FuniUN.c & FuniWE.c
				} else {
	
					#IF Finb is a F.SNPs coeff (Fas, Funiu, Funiw) do LDMS strat
					if(Finb %in% c("Fas.c", "Funi.un.c", "Funi.we.c")){
	
						#Corresponing matrix LDMS
						if(Finb == "Fas.c"){
						ldmsMat = FasbinsLD
	
						} else if (Finb == "Funi.un.c"){
						ldmsMat = FuniUNbinsLD
	
						} else if (Finb == "Funi.we.c"){
						ldmsMat = FuniWEbinsLD
	
						}
	
						#LDMS F with GRM AS
						assign(x = paste0(Finb, "_", sim, "_MatGRMAS_LDMS"), value = mclapply(1:100, function(z) {print(z)
						lmm.aireml(Y = get(sim)[[2]][,z], K = GRMas, X = cbind(1,as.matrix(ldmsMat)), verbose = T)[[9]][2]}
						, mc.cores = 30))
	
						#LDMS F with GRM GCTAwe
						assign(x = paste0(Finb, "_", sim, "_MatGRMkc0_LDMS"), value = mclapply(1:100, function(z) {print(z)
						lmm.aireml(Y = get(sim)[[2]][,z], K = 2*GCTAwe, X = cbind(1,as.matrix(ldmsMat)), verbose = T)[[9]][2]}
						, mc.cores = 30))
	
						#LDMS F with GRM GCTAun
						assign(x = paste0(Finb, "_", sim, "_MatGRMgcta_LDMS"), value = mclapply(1:100, function(z) {print(z)
						lmm.aireml(Y = get(sim)[[2]][,z], K = GCTAun[[4]], X = cbind(1,as.matrix(ldmsMat)), verbose = T)[[9]][2]}
						, mc.cores = 30))
					}
				}
			}
	
			#After estimating all the ID est for this scenario, bringing everything together
	
			#If it is the WORLD, take the Fs.w
			if("${POP}" == "WORLD"){
	
				#Extract objects to assign
				assign(x = paste0("ID_EST_", sim), value = data.frame(get(sim)[[6]][,c(1:4,6,8,10,12,14,16,18:21)],
				Fas.lme.AS = unlist(get(paste0("Fas.w_", sim, "_MatGRMAS"))),
				FuniUN.lme.AS = unlist(get(paste0("Funi.un.w_", sim, "_MatGRMAS"))),
				FuniWE.lme.AS = unlist(get(paste0("Funi.we.w_", sim, "_MatGRMAS"))),
				Froh.BCFTools.100KB.lme.AS = unlist(get(paste0("FROHs.BCFTools.100KB.w_", sim, "_MatGRMAS"))),
				Froh.PLINK.100KB.lme.AS = unlist(get(paste0("FROHs.PLINK.100KB.w_", sim, "_MatGRMAS"))),
				Froh.BCFTools.1MB.lme.AS = unlist(get(paste0("FROHs.BCFTools.1MB.w_", sim, "_MatGRMAS"))),
				Froh.PLINK.1MB.lme.AS = unlist(get(paste0("FROHs.PLINK.1MB.w_", sim, "_MatGRMAS"))),
				Fas.lme.AS.LDMS = unlist(get(paste0("Fas.w_", sim, "_MatGRMAS_LDMS"))),
				FuniUN.lme.AS.LDMS = unlist(get(paste0("Funi.un.w_", sim, "_MatGRMAS_LDMS"))),
				FuniWE.lme.AS.LDMS = unlist(get(paste0("Funi.we.w_", sim, "_MatGRMAS_LDMS"))),
				Fas.lme.kc0 = unlist(get(paste0("Fas.w_", sim, "_MatGRMkc0"))),
				FuniUN.lme.kc0 = unlist(get(paste0("Funi.un.w_", sim, "_MatGRMkc0"))),
				FuniWE.lme.kc0 = unlist(get(paste0("Funi.we.w_", sim, "_MatGRMkc0"))),
				Froh.BCFTools.100KB.lme.kc0 = unlist(get(paste0("FROHs.BCFTools.100KB.w_", sim, "_MatGRMkc0"))),
				Froh.PLINK.100KB.lme.kc0 = unlist(get(paste0("FROHs.PLINK.100KB.w_", sim, "_MatGRMkc0"))),
				Froh.BCFTools.1MB.lme.kc0 = unlist(get(paste0("FROHs.BCFTools.1MB.w_", sim, "_MatGRMkc0"))),
				Froh.PLINK.1MB.lme.kc0 = unlist(get(paste0("FROHs.PLINK.1MB.w_", sim, "_MatGRMkc0"))),
				Fas.lme.kc0.LDMS = unlist(get(paste0("Fas.w_", sim, "_MatGRMkc0_LDMS"))),
				FuniUN.lme.kc0.LDMS = unlist(get(paste0("Funi.un.w_", sim, "_MatGRMkc0_LDMS"))),
				FuniWE.lme.kc0.LDMS = unlist(get(paste0("Funi.we.w_", sim, "_MatGRMkc0_LDMS"))),
				Fas.lme.gcta = unlist(get(paste0("Fas.w_", sim, "_MatGRMgcta"))),
				FuniUN.lme.gcta = unlist(get(paste0("Funi.un.w_", sim, "_MatGRMgcta"))),
				FuniWE.lme.gcta = unlist(get(paste0("Funi.we.w_", sim, "_MatGRMgcta"))),
				Froh.BCFTools.100KB.lme.gcta = unlist(get(paste0("FROHs.BCFTools.100KB.w_", sim, "_MatGRMgcta"))),
				Froh.PLINK.100KB.lme.gcta = unlist(get(paste0("FROHs.PLINK.100KB.w_", sim, "_MatGRMgcta"))),
				Froh.BCFTools.1MB.lme.gcta = unlist(get(paste0("FROHs.BCFTools.1MB.w_", sim, "_MatGRMgcta"))),
				Froh.PLINK.1MB.lme.gcta = unlist(get(paste0("FROHs.PLINK.1MB.w_", sim, "_MatGRMgcta"))),
				Fas.lme.gcta.LDMS = unlist(get(paste0("Fas.w_", sim, "_MatGRMgcta_LDMS"))),
				FuniUN.lme.gcta.LDMS = unlist(get(paste0("Funi.un.w_", sim, "_MatGRMgcta_LDMS"))),
				FuniWE.lme.gcta.LDMS = unlist(get(paste0("Funi.we.w_", sim, "_MatGRMgcta_LDMS")))))
	
				#Change first colnames to rm the .w for all F
				setnames(get(paste0("ID_EST_", sim))[5:11], c("Fas","Funi.un","Funi.we","Froh.BCFTools.100KB","Froh.PLINK.100KB","Froh.BCFTools.1MB","Froh.PLINK.1MB"))
	
			#If not take Fs.c
			} else {
	
				#Extract objects to assign
				assign(x = paste0("ID_EST_", sim), value = data.frame(get(sim)[[6]][,c(1:4,5,7,9,11,13,15,17,19:21)],
				Fas.lme.AS = unlist(get(paste0("Fas.c_", sim, "_MatGRMAS"))),
				FuniUN.lme.AS = unlist(get(paste0("Funi.un.c_", sim, "_MatGRMAS"))),
				FuniWE.lme.AS = unlist(get(paste0("Funi.we.c_", sim, "_MatGRMAS"))),
				Froh.BCFTools.100KB.lme.AS = unlist(get(paste0("FROHs.BCFTools.100KB.c_", sim, "_MatGRMAS"))),
				Froh.PLINK.100KB.lme.AS = unlist(get(paste0("FROHs.PLINK.100KB.c_", sim, "_MatGRMAS"))),
				Froh.BCFTools.1MB.lme.AS = unlist(get(paste0("FROHs.BCFTools.1MB.c_", sim, "_MatGRMAS"))),
				Froh.PLINK.1MB.lme.AS = unlist(get(paste0("FROHs.PLINK.1MB.c_", sim, "_MatGRMAS"))),
				Fas.lme.AS.LDMS = unlist(get(paste0("Fas.c_", sim, "_MatGRMAS_LDMS"))),
				FuniUN.lme.AS.LDMS = unlist(get(paste0("Funi.un.c_", sim, "_MatGRMAS_LDMS"))),
				FuniWE.lme.AS.LDMS = unlist(get(paste0("Funi.we.c_", sim, "_MatGRMAS_LDMS"))),
				Fas.lme.kc0 = unlist(get(paste0("Fas.c_", sim, "_MatGRMkc0"))),
				FuniUN.lme.kc0 = unlist(get(paste0("Funi.un.c_", sim, "_MatGRMkc0"))),
				FuniWE.lme.kc0 = unlist(get(paste0("Funi.we.c_", sim, "_MatGRMkc0"))),
				Froh.BCFTools.100KB.lme.kc0 = unlist(get(paste0("FROHs.BCFTools.100KB.c_", sim, "_MatGRMkc0"))),
				Froh.PLINK.100KB.lme.kc0 = unlist(get(paste0("FROHs.PLINK.100KB.c_", sim, "_MatGRMkc0"))),
				Froh.BCFTools.1MB.lme.kc0 = unlist(get(paste0("FROHs.BCFTools.1MB.c_", sim, "_MatGRMkc0"))),
				Froh.PLINK.1MB.lme.kc0 = unlist(get(paste0("FROHs.PLINK.1MB.c_", sim, "_MatGRMkc0"))),
				Fas.lme.kc0.LDMS = unlist(get(paste0("Fas.c_", sim, "_MatGRMkc0_LDMS"))),
				FuniUN.lme.kc0.LDMS = unlist(get(paste0("Funi.un.c_", sim, "_MatGRMkc0_LDMS"))),
				FuniWE.lme.kc0.LDMS = unlist(get(paste0("Funi.we.c_", sim, "_MatGRMkc0_LDMS"))),
				Fas.lme.gcta = unlist(get(paste0("Fas.c_", sim, "_MatGRMgcta"))),
				FuniUN.lme.gcta = unlist(get(paste0("Funi.un.c_", sim, "_MatGRMgcta"))),
				FuniWE.lme.gcta = unlist(get(paste0("Funi.we.c_", sim, "_MatGRMgcta"))),
				Froh.BCFTools.100KB.lme.gcta = unlist(get(paste0("FROHs.BCFTools.100KB.c_", sim, "_MatGRMgcta"))),
				Froh.PLINK.100KB.lme.gcta = unlist(get(paste0("FROHs.PLINK.100KB.c_", sim, "_MatGRMgcta"))),
				Froh.BCFTools.1MB.lme.gcta = unlist(get(paste0("FROHs.BCFTools.1MB.c_", sim, "_MatGRMgcta"))),
				Froh.PLINK.1MB.lme.gcta = unlist(get(paste0("FROHs.PLINK.1MB.c_", sim, "_MatGRMgcta"))),
				Fas.lme.gcta.LDMS = unlist(get(paste0("Fas.c_", sim, "_MatGRMgcta_LDMS"))),
				FuniUN.lme.gcta.LDMS = unlist(get(paste0("Funi.un.c_", sim, "_MatGRMgcta_LDMS"))),
				FuniWE.lme.gcta.LDMS = unlist(get(paste0("Funi.we.c_", sim, "_MatGRMgcta_LDMS")))))
	
				#Change first colnames to rm the .c for all F
				setnames(get(paste0("ID_EST_", sim))[5:11], c("Fas","Funi.un","Funi.we","Froh.BCFTools.100KB","Froh.PLINK.100KB","Froh.BCFTools.1MB","Froh.PLINK.1MB"))
			}
		}
	
		#Then save the R env. with FINAL results in /work
		save(list=ls(pattern="ID_EST_"), file = "${output}_COMPLETE_ID_EST.RData")
		#save the ENTIRE R env. in /nas
		#save(list = ls(all.names = TRUE), file = paste0("/scratch/elavanc1/ID_EST_HUMANS/Feb_2022_ID_SIM_Rdata/", strsplit("${output}", "/")[[1]][length(strsplit("${output}", "/")[[1]])], ".RData"))
	
EOF
	
	}
	
	#export the function
	export -f SIMULATING_ID
	
	SIMULATING_ID ./data/bedmatrices/${POP}_ch_1_22 \
		./data//F_datasets/${POP}_allchr_allFs.txt \
		./data/GRMs/${POP}_AS_GRM_ch_1_22.RDS \
		./data/GRMs/${POP}_GCTAun_GRM_ch_1_22.RDS \
		./data/GRMs/${POP}_GCTAwe_GRM_ch_1_22.RDS \
		./data/LD_MAF_bins/${POP}_LD_MAF_F_per_bins.Fas \
		./data/LD_MAF_bins/${POP}_LD_MAF_F_per_bins.FuniWE \
		./data/LD_MAF_bins/${POP}_LD_MAF_F_per_bins.FuniUN \
		./data/ID_EST_Datasets/${POP}_ID3_simus

done
