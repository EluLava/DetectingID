#create outputdir if does not exists
mkdir -p ./data/ID_EST_Datasets/

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

	#Sixth argument passed to the function is output file
	output=$6

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


	##### Set SIM parameters ####

	#Nb of causal loci
	ncaus = 1000

	#Filtering on MAF = 0.001 to avoid super rare alleles causing super strong effects
	maf = 0.001

	#ID (and dominance effect) = 3
	b = 3

	#Additive effects = 10
	s = 10

	#beta dist parameters, 1 & 1 = uniform distribution, FYI we tried changing the distribution and it has almost no effect !
	b1 = 1
	b2 = 1

	### RUN THE SIMULATION LM for all the scenarios (only the last one presented in main text, the rest are in supplementary material) ###

	assign(x = paste0("PGPED_LMM_ID",b,"_selhF_selpF_demaF"), value = id.est(Fs[,4:ncol(Fs)], genos = bed, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100,
	selh = FALSE, selp = FALSE, direct = FALSE, b1=b1, b2=b2))

	assign(x = paste0("PGPED_LMM_ID",b,"_selhT_selpF_demaF"), value = id.est(Fs[,4:ncol(Fs)], genos = bed, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100,
	selh = TRUE, selp = FALSE, direct = FALSE, b1=b1, b2=b2))

	assign(x = paste0("PGPED_LMM_ID",b,"_selhF_selpT_demaF"), value = id.est(Fs[,4:ncol(Fs)], genos = bed, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100,
	selh = FALSE, selp = TRUE, direct = FALSE, b1=b1, b2=b2))

	assign(x = paste0("PGPED_LMM_ID",b,"_selhT_selpT_demaF"), value = id.est(Fs[,4:ncol(Fs)], genos = bed, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100,
	selh = TRUE, selp = TRUE, direct = FALSE, b1=b1, b2=b2))

	assign(x = paste0("PGPED_LMM_ID",b,"_selhF_selpF_demaT"), value = id.est(Fs[,4:ncol(Fs)], genos = bed, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100,
	selh = FALSE, selp = FALSE, direct = TRUE, b1=b1, b2=b2))

	assign(x = paste0("PGPED_LMM_ID",b,"_selhT_selpF_demaT"), value = id.est(Fs[,4:ncol(Fs)], genos = bed, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100,
	selh = TRUE, selp = FALSE, direct = TRUE, b1=b1, b2=b2))

	assign(x = paste0("PGPED_LMM_ID",b,"_selhF_selpT_demaT"), value = id.est(Fs[,4:ncol(Fs)], genos = bed, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100,
	selh = FALSE, selp = TRUE, direct = TRUE, b1=b1, b2=b2))

	assign(x = paste0("PGPED_LMM_ID",b,"_selhT_selpT_demaT"), value = id.est(Fs[,4:ncol(Fs)], genos = bed, filter.maf = maf, n.causal = ncaus, b = b, s = s, n.rep = 100,
	selh = TRUE, selp = TRUE, direct = TRUE, b1=b1, b2=b2))

	#Extract Simulated scenarios
	simuScen = ls(pattern = paste0("PGPED_LMM_ID",b,"_"))


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
	}

		#After estimating all the ID est for this scenario, bring everything together

		#Extract objects to assign
		assign(x = paste0("ID_EST_", sim), value = data.frame(get(sim)[[6]][,c(1:13)],
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

	#Then save the R env. with FINAL results
	save(list=ls(pattern="ID_EST_"), file = "${output}_COMPLETE_ID_EST.RData")
	#save the ENTIRE R env. in /scratch
	#save(list = ls(all.names = TRUE), file = paste0("/scratch/elavanc1/ID_EST_HUMANS/Feb_2022_ID_SIM_Rdata/", strsplit("${output}", "/")[[1]][length(strsplit("${output}", "/")[[1]])], ".RData"))

EOF

}

#export the function
export -f SIMULATING_ID

SIMULATING_ID ./data/bedmatrices/PGPED_ch_1_2 \
    ./data/F_datasets/PGPED_PED_allchr_allFs.txt \
    ./data/GRMs/PGPED_AS_GRM_ch_1_22.RDS \
    ./data/GRMs/PGPED_GCTAun_GRM_ch_1_22.RDS \
	./data/GRMs/PGPED_GCTAun_GRM_ch_1_22.RDS \
    ./data/ID_EST_Datasets/PGPED_ID3_simus


