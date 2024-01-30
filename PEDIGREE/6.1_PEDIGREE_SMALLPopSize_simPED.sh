#! /bin/bash

mkdir -p ./data/msgenos/

#Define function
POLYPED_ReSim(){

    #First argument passed to the script is continent (for # indvs subsampling)
    SampSize=$1
    #Second argument is the numbers of founders
    founderSize=$2

    #GENERATE THE FOUNDERS's genotypes file (nb founders * 2 because it is haploid genomes)
    mspms $((${founderSize}*2)) 1 -t 4000 -r 4000 100000000 -p 9 > ./data/msgenos/foundersb_${SampSize}.txt

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
        cl <- 1
        registerDoParallel(cl)

	source("../1KG/functions/functions.R")

	nbindvs = as.numeric("${SampSize}")
	nbfounders = as.numeric("${founderSize}")

	#set condition to continue the while as long as we don't have a pedigree of the correct size
	condition = TRUE

	#As long as we are not ahppy, we re-simulate a pedigree
	while(condition){

	    #simulate pedigree
	    ped.pg = try(JGTeach::buildped.2sexes(founders.m = (nbfounders/2), founders.f = (nbfounders/2), fert=1.0, death.rate = 0.5, breed.prop = 0.1, n.tstep = 15))

	    #if sumulation fails (mostly because all individuals died), we re-try
	    while(inherits(ped.pg, "try-error")){

		ped.pg = try(JGTeach::buildped.2sexes(founders.m = (nbfounders/2), founders.f = (nbfounders/2), fert=1.0, death.rate = 0.5, breed.prop = 0.1, n.tstep = 15))

	    }

	    #If not the correct number of individuals in the pedigree, we re-simulate, else we move on
	    if(nrow(ped.pg) == nbindvs){

		condition = FALSE
	    }
	}

	#Read the founders file
    myfounds <- hierfstat::ms2dos("./msgenos/foundersb_${SampSize}.txt")

	#Get Fped
	Fped.pg <- diag(grm2kinship(pedARM(ped.pg[,2],ped.pg[,3])))

	#Create genotypes
	bedpg.20.ch1 <- drop.along.ped(ped = ped.pg, founders.genotypes = myfounds[['alldat']][,myfounds[['bim']][,'chr'] == 1], maplength = 20)
	bedpg.20.ch1 <- as.bed.matrix(bedpg.20.ch1[,,1] + bedpg.20.ch1[,,2])
	bedpg.20.ch1@snps[,'pos'] <- myfounds[['bim']][,'pos'][myfounds[['bim']][,'chr']==1]
	bedpg.20.ch1@snps[,'chr']<-myfounds[['bim']][,'chr'][myfounds[['bim']][,'chr']==1]
	bedpg.20.ch1@ped[,'father']<-ped.pg[,'sire']
	bedpg.20.ch1@ped[,'mother']<-ped.pg[,'dam']

	#Save the bed matrix
	write.bed.matrix(bedpg.20.ch1, "./bedmatrices/BED_ReducedPGPED_${SampSize}")

	#Merge Fped with indvs IDs
	F.data.1 = as.data.frame(cbind(ID = ped.pg[,'ind'], Fped = Fped.pg))

	#Get the different Fs
	Fas = get.fasT(bed = bedpg.20.ch1, nb.cores = 1)
	#Add indvs
	Fas.ped = as.data.frame(cbind(ID = bedpg.20.ch1@ped[,'id'], Fas = Fas))

	#merge with Fped
	F.data.2 = merge(F.data.1, Fas.ped, by = "ID")

	funi.un = get.funiuT(bed = bedpg.20.ch1, nb.cores = 1)
        #Add indvs
	Funi.un.ped = as.data.frame(cbind(ID = bedpg.20.ch1@ped[,'id'], Funi.un = funi.un))

        #merge with Fped
        F.data.3 = merge(F.data.2, Funi.un.ped, by = "ID")

        funi.we = get.funiwT(bed = bedpg.20.ch1, nb.cores = 1)
        #Add indvs
        Funi.we.ped = as.data.frame(cbind(ID = bedpg.20.ch1@ped[,'id'], Funi.we = funi.we))

        #merge with Fped
        F.data.4 = merge(F.data.3, Funi.we.ped, by = "ID")

	#Save the inbreeding coefficients file
	write.table(F.data.4, "./data/F_datasets/ReducedPGPED_${SampSize}_SNPsbasedF.txt", quote = F, col.names = T, row.names = F)

EOF

    #get 0 for second column of bim file
    awk '{print $1"\t"$2"\t"0"\t"int($4)"\t"$5"\t"$6}' ./data/bedmatrices/BED_ReducedPGPED_${SampSize}.bim >> $$; mv $$ ./data/bedmatrices/BED_ReducedPGPED_${SampSize}.bim

    #from bed to VCF for ROHs and HBD segments
    plink --bfile ./data/bedmatrices/BED_ReducedPGPED_${SampSize} --recode vcf bgz --out ./data/VCFs/VCF_ReducedPGPED_${SampSize}

    #index
    bcftools index ./data/VCFs/VCF_ReducedPGPED_${SampSize}.vcf.gz

    #Creating directory for BCFTools hbd seg
    mkdir -p ./data/F_datasets/CallingROHsBCFTools

    #Calling segments
    bcftools roh -o ./data/F_datasets/CallingROHsBCFTools/ReducedPGPED_${SampSize}.txt -O r -e "-" -G 30 ./data/VCFs/VCF_ReducedPGPED_${SampSize}.vcf.gz

    ##Calling ROHs with PLINK

    #Create directory if does not exists
    mkdir -p ./data/F_datasets/CallingROHsPLINK

    #100KB
    plink --bfile ./data/bedmatrices/BED_ReducedPGPED_${SampSize} --homozyg --homozyg-kb 100 --out ./data/F_datasets/CallingROHsPLINK/ReducedPGPED_${SampSize}_100KB
    #1MB
    plink --bfile ./data/bedmatrices/BED_ReducedPGPED_${SampSize} --homozyg --homozyg-kb 1000 --out ./data/F_datasets/CallingROHsPLINK/ReducedPGPED_${SampSize}_1MB

    #BACK to R for filtering on size
    R --vanilla << EOF

        dta = read.table("./data/F_datasets/CallingROHsBCFTools/ReducedPGPED_${SampSize}.txt")
        colnames(dta) = c("RG","Sample","Chromosome","Start","End", "Length_bp","Nbof_markers", "Quality")

        #Select fragments with the minimum size
        dta100KB = dta[dta[,6] >= 100000,]
        dta1MB = dta[dta[,6] >= 1000000,]

        #Create empty F vectors
        FROH.BCFTools.100KB = vector(mode = "numeric", length = length(unique(dta[,2])))
        names(FROH.BCFTools.100KB) = unique(dta[,2])

        FROH.BCFTools.1MB = FROH.BCFTools.100KB

        #Loop through indvs and get F
        for (indv in unique(dta[,2])) {

          ## 100KB

          #Subset the data
          dtasub = dta100KB[dta100KB[,2] == indv,]
          #Estimate F
          FROH.BCFTools.100KB[names(FROH.BCFTools.100KB) == indv] = sum(dtasub[,6])/1e8

          ## 1MB

          #Subset the data
          dtasub = dta1MB[dta1MB[,2] == indv,]
          #Estimate F
          FROH.BCFTools.1MB[names(FROH.BCFTools.1MB) == indv] = sum(dtasub[,6])/1e8

        }

        #Paste names as other column
        FROH.BCFTools = cbind(sample=names(FROH.BCFTools.100KB), FROH.BCFTools.100KB = FROH.BCFTools.100KB, FROH.BCFTools.1MB = FROH.BCFTools.1MB)

        #as numeric
        FROH.BCFTools[,2] = as.numeric(as.character(FROH.BCFTools[,2]))
        FROH.BCFTools[,3] = as.numeric(as.character(FROH.BCFTools[,3]))

        #samples' names
        FROH.BCFTools[,1] = sapply(FROH.BCFTools[,1], function(x){strsplit(x, split = "_")[[1]][2]})

        #Save the dataframe
        write.table(FROH.BCFTools, "./data/F_datasets/CallingROHsBCFTools/ReducedPGPED_${SampSize}.txt", quote = F, col.names = T, row.names = F)


    dta100KB = read.table("./data/F_datasets/CallingROHsPLINK/ReducedPGPED_${SampSize}_100KB.hom.indiv", header = T)
    dta1MB = read.table("./data/F_datasets/CallingROHsPLINK/ReducedPGPED_${SampSize}_1MB.hom.indiv", header = T)

    #Create empty F vectors
    FROH.PLINK.100KB = vector(mode = "numeric", length = length(unique(dta100KB[,2])))
    names(FROH.PLINK.100KB) = unique(dta100KB[,2])

    FROH.PLINK.1MB = FROH.PLINK.100KB

    #Loop through indvs and get F
    for (indv in unique(dta100KB[,2])) {

      ## 100KB

      #Estimate F, divided by e6 becasue in KB
      FROH.PLINK.100KB[names(FROH.PLINK.100KB) == indv] = dta100KB[,5][dta100KB[,2] == indv]/1e6

      ## 1MB

      #Estimate F
      FROH.PLINK.1MB[names(FROH.PLINK.1MB) == indv] = dta1MB[,5][dta1MB[,2] == indv]/1e6

    }

    #sanity check
    all.equal(names(FROH.PLINK.100KB),names(FROH.PLINK.1MB))

    #Paste names as other column
    FROH.PLINK = as.data.frame(cbind(sample=names(FROH.PLINK.100KB), FROH.PLINK.100KB = FROH.PLINK.100KB, FROH.PLINK.1MB = FROH.PLINK.1MB))

    #as numeric
    FROH.PLINK[,2] = as.numeric(as.character(FROH.PLINK[,2]))
    FROH.PLINK[,3] = as.numeric(as.character(FROH.PLINK[,3]))

    #Save the dataframe
    write.table(FROH.PLINK, "./data/F_datasets/CallingROHsPLINK/ReducedPGPED_${SampSize}.txt", quote = F, col.names = T, row.names = F)

    #### OK now we want to estimate the relatedness matrices

    #List packages needed
	list.of.packages <- c("SNPRelate", "hierfstat", "devtools")
	new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
	if(length(new.packages)){install.packages(new.packages)}
	if(!("JGTeach" %in% installed.packages()[,"Package"])){library(devtools); install_github("jgx65/JGTeach")}

	library(SNPRelate)
	library(hierfstat)
	library(JGTeach)

	#From bed to GDS
	snpgdsBED2GDS(bed.fn = "./data/bedmatrices/BED_ReducedPGPED_${SampSize}.bed", fam.fn = "./data/bedmatrices/BED_ReducedPGPED_${SampSize}.fam",
			   bim.fn = "./data/bedmatrices/BED_ReducedPGPED_${SampSize}.bim", out.gdsfn = "./data/bedmatrices/BED_ReducedPGPED_${SampSize}.gds")

	#Read the gds object
	gdsMat = snpgdsOpen("./data/bedmatrices/BED_ReducedPGPED_${SampSize}.gds")

	#From GDS to GRM gcta UNweighted
	GCTAun = snpgdsGRM(gdsMat, method = "GCTA")
	#save it
	saveRDS(GCTAun, "./data/GRMs/ReducedPGPED_${SampSize}_GCTAun_GRM.RDS")

	#From GDS to dosage
	dos = snpgdsGetGeno(gdsMat)

	#matching Mat
	allMatchingMat = matching(dos)

	#Estimate pairwise kinships (FROM ALLELES MATCHING)
	kas = beta.dosage(allMatchingMat, MATCHING = TRUE)
	#Pass to GRM matrix
	grmas = kinship2grm(kas)
	#write GRMas
	saveRDS(grmas, "./data/GRMs/ReducedPGPED_${SampSize}_AS_GRM.RDS")

	#pairwise kinships (STANDARD) ratio of averages: GCTA Weigthed
	GCTAwe = Kc0(allMatchingMat, matching=TRUE)
	#save it
	saveRDS(GCTAwe, "./data/GRMs/ReducedPGPED_${SampSize}_GCTAwe_GRM.RDS")

EOF

}

export -f POLYPED_ReSim

#Extract runFILE from inputFILE: one array per superpop
runFILE=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ./data/inputFILES/POLYPED_SUB_sampSize.list)

POLYPED_ReSim ${runFILE}
