#We'll estimate the F for the entire WORLD population first, then we'll loop through continents and redo the F estimation !

#Create directory if does not exists
mkdir -p VCFs

#Pass from BED to VCF for ROHs and HBD segments estimation
plink --bfile ./data/bedmatrices/WORLD_ch_1_22 --recode vcf iid bgz --out ./data/VCFs/1KG_WORLD

#index
bcftools index ./data/VCFs/1KG_WORLD.vcf.gz

##Calling HBD segments with BCFTools

#Create directory if does not exists
mkdir -p CallingROHsBCFTools

#Calling segments
bcftools roh -o ./data/CallingROHsBCFTools/1KG_WORLD.txt -O r -e "-" -G 30 ./data/VCFs/1KG_WORLD.vcf.gz

##Calling ROHs with PLINK

#Create directory if does not exists
mkdir -p CallingROHsPLINK

#100KB
plink --bfile ./data/bedmatrices/WORLD_ch_1_22 --homozyg --homozyg-kb 100 --out ./data/CallingROHsPLINK/1KG_WORLD_100KB
#1MB
plink --bfile ./data/bedmatrices/WORLD_ch_1_22 --homozyg --homozyg-kb 1000 --out ./data/CallingROHsPLINK/1KG_WORLD_1MB

#create output directory if does not exists
mkdir -p ./data/F_datasets/TMP

R --vanilla << EOF

    ###################################
    #### ROHs called with BCFTools ####
    ###################################

    dta = read.table("./data/CallingROHsBCFTools/1KG_WORLD.txt")
    colnames(dta) = c("RG","Sample","Chromosome","Start","End", "Length_bp","Nbof_markers", "Quality")
    
    #Select fragments with the minimum size
    dta100KB = dta[dta[,6] >= 100000,]
    dta1MB = dta[dta[,6] >= 1000000,]
    
    #Create empty F vectors
    FROH.BCFTools.100KB.w = vector(mode = "numeric", length = length(unique(dta[,2])))
    names(FROH.BCFTools.100KB.w) = unique(dta[,2])
    
    FROH.BCFTools.1MB.w = FROH.BCFTools.100KB.w
    
    #Loop through indvs and get F
    for (indv in unique(dta[,2])) {
      
      ## 100KB
      
      #Subset the data
      dtasub = dta100KB[dta100KB[,2] == indv,]
      #Estimate F
      FROH.BCFTools.100KB.w[names(FROH.BCFTools.100KB.w) == indv] = sum(dtasub[,6])/3.2e9
      
      ## 1MB
      
      #Subset the data
      dtasub = dta1MB[dta1MB[,2] == indv,]
      #Estimate F
      FROH.BCFTools.1MB.w[names(FROH.BCFTools.1MB.w) == indv] = sum(dtasub[,6])/3.2e9
      
    }
    
    #Paste names as other column
    FROH.BCFTools = cbind(sample=names(FROH.BCFTools.100KB.w), FROH.BCFTools.100KB.w = FROH.BCFTools.100KB.w, FROH.BCFTools.1MB.w = FROH.BCFTools.1MB.w)
    
    #as numeric
    FROH.BCFTools[,2] = as.numeric(as.character(FROH.BCFTools[,2]))
    FROH.BCFTools[,3] = as.numeric(as.character(FROH.BCFTools[,3]))
    
    #sample names
    FROH.BCFTools[,1] = sapply(FROH.BCFTools[,1], function(x){strsplit(x, split = "_")[[1]][2]})
    
    #Save the dataframe
    write.table(FROH.BCFTools, "./data/F_datasets/TMP/WORLD_allchr_both_FROH.BCFTools.size.txt", quote = F, col.names = T, row.names = F)


    ################################
    #### ROHs called with PLINK ####
    ################################

    dta100KB = read.table("./data/CallingROHsPLINK/1KG_WORLD_100KB.hom.indiv", header = T)
    dta1MB = read.table("./data/CallingROHsPLINK/1KG_WORLD_100KB.hom.indiv", header = T)
    
    #Create empty F vectors
    FROH.PLINK.100KB.w = vector(mode = "numeric", length = length(unique(dta100KB[,2])))
    names(FROH.PLINK.100KB.w) = unique(dta100KB[,2])
    
    FROH.PLINK.1MB.w = FROH.PLINK.100KB.w
    
    #Loop through indvs and get F
    for (indv in unique(dta100KB[,2])) {
      
      ## 100KB
      
      #Estimate F, divided by e6 becasue in KB
      FROH.PLINK.100KB.w[names(FROH.PLINK.100KB.w) == indv] = dta100KB[,5][dta100KB[,2] == indv]/3.2e6
      
      ## 1MB
      
      #Estimate F
      FROH.PLINK.1MB.w[names(FROH.PLINK.1MB.w) == indv] = dta1MB[,5][dta1MB[,2] == indv]/3.2e6
      
    }
    
    #sanity check
    all.equal(names(FROH.PLINK.100KB.w),names(FROH.PLINK.1MB.w))
    
    #Paste names as other column
    FROH.PLINK = as.data.frame(cbind(sample=names(FROH.PLINK.100KB.w), FROH.PLINK.100KB.w = FROH.PLINK.100KB.w, FROH.PLINK.1MB.w = FROH.PLINK.1MB.w))
    
    #as numeric
    FROH.PLINK[,2] = as.numeric(as.character(FROH.PLINK[,2]))
    FROH.PLINK[,3] = as.numeric(as.character(FROH.PLINK[,3]))
    
    #Save the dataframe
    write.table(FROH.PLINK, "./data/F_datasets/TMP/WORLD_allchr_both_FROH.PLINK.size.txt", quote = F, col.names = T, row.names = F)


    ##############################
    #### FAS, FUNIu and FUNIw ####
    ##############################

    if(!("gaston" %in% installed.packages()[,"Package"])){install.packages("gaston")}
    library(gaston)
    setThreadOptions(5)

    #Load the scripts
    source("./functions/functions.R")

    #Read bed matrix
    bed = read.bed.matrix("./data/bedmatrices/WORLD_ch_1_22")

    #Extract sample name info FROM BED
    sampleBED = bed@ped[,2]

    #Get FuniUN
    FuniUN = get.funiuT(bed, nb.cores = 5)
    #Get FuniWE
    FuniWE = get.funiwT(bed, nb.cores = 5)
    #Get fas
    Fas = get.fasT(bed, nb.cores = 5)

    #merge everyting
    SNPbasedFs = cbind(sample=names(FuniUN), Fas.w = Fas, Funi.un.w = FuniUN, Funi.we.w = FuniWE)

    #Save the table
    write.table(SNPbasedFs, "./data/F_datasets/TMP/WORLD_allchr_FasbothFuni.txt", quote = F, col.names = T, row.names = F)


    #########################################
    #### Merging and sample descriptions ####
    #########################################

    #Read sample description file
    sampleDESC = read.table("./data/SAMPLE_DESC/integrated_call_samples_v3.20130502.ALL.panel")

    #merge sample desc POP and CONT with our estimations
    dta = merge(sampleDESC[,1:3], SNPbasedFs, by = "sample")
    #merge with FROH BCFTOOLS and PLINK
    dta.2 = merge(dta, FROH.BCFTools, by = "sample")
    dta.3 = merge(dta.2, FROH.PLINK, by = "sample")

    write.table(dta.3, "./data/F_datasets/TMP/WORLD_allchr_FasbothFunibothFROHs.txt", quote = F, col.names = T, row.names = F)

    ##########################################################################################
    #### Looping through continents to subsample BEDmatrix and estimate Fas and both Funi ####
    ##########################################################################################

    #Create empty dataframe that we'll fill with fas and both funi for all continents
    SNPbasedFscont = as.data.frame(matrix(nrow = 0, ncol = 4))
    colnames(SNPbasedFscont) = c("sample", "Fas.c", "Funi.un.c", "Funi.we.c")

    #Loop through continents
    for(cont in unique(dtaDESC[,3])){

        #Extract samples names from continent
        dta_samples = unique(dtaDESC[,1][dtaDESC[,3] == cont])

        #susabmple bed matrix
        bedSample = select.inds(bed, id %in% dta_samples)

        #filter for maf > 0 --> non monomorphic SNPs
        bedSample.snps = select.snps(bedSample, maf > 0)

        #Save the bed matrix for ROHs estimation later (and GRM matrices btw)
        write.bed.matrix(bedSample.snps, paste0("./data/bedmatrices/bed_", cont, "_ch1_22"))

        #Get FuniUN
        FuniUN = get.funiuT(bedSample.snps, nb.cores = 5)
        #Get FuniWE
        FuniWE = get.funiwT(bedSample.snps, nb.cores = 5)
        #Get fas
        Fas = get.fasT(bedSample.snps, nb.cores = 5)

        #merge everyting
        SNPbasedFscont = rbind(SNPbasedFscont, cbind(sample=names(FuniUN), Fas.c = Fas, Funi.un.c = FuniUN, Funi.we.c = FuniWE))
    }

    #Save the table
    write.table(SNPbasedFscont, "./data/F_datasets/TMP/PerCONT_allchr_FasbothFuni.txt", quote = F, col.names = T, row.names = F)


EOF

#Now we'll loop through continents for both FROHs estimation

for CONT in {AFR,AMR,EAS,EUR,SAS}; do

    #Pass from BED to VCF for ROHs and HBD segments estimation
    plink --bfile ./data/bedmatrices/${CONT}_ch_1_22 --recode vcf iid bgz --out ./data/VCFs/1KG_${CONT}
    
    #index
    bcftools index ./data/VCFs/1KG_${CONT}.vcf.gz
    
    ##Calling HBD segments with BCFTools
    
    #Calling segments
    bcftools roh -o ./data/CallingROHsBCFTools/1KG_${CONT}.txt -O r -e "-" -G 30 ./data/VCFs/1KG_${CONT}.vcf.gz
    
    ##Calling ROHs with PLINK

    #100KB
    plink --bfile ./data/bedmatrices/${CONT}_ch_1_22 --homozyg --homozyg-kb 100 --out ./data/CallingROHsPLINK/1KG_${CONT}_100KB
    #1MB
    plink --bfile ./data/bedmatrices/${CONT}_ch_1_22 --homozyg --homozyg-kb 1000 --out ./data/CallingROHsPLINK/1KG_${CONT}_1MB


    #Re-open R for both FROH estimation
    R --vanilla << EOF

    ###################################
    #### ROHs called with BCFTools ####
    ###################################

    dta = read.table("./data/CallingROHsBCFTools/1KG_${CONT}.txt")
    colnames(dta) = c("RG","Sample","Chromosome","Start","End", "Length_bp","Nbof_markers", "Quality")
    
    #Subset per length
    dta100KB = dta[dta[,6] >= 100000,]
    dta1MB = dta[dta[,6] >= 1000000,]
    
    #Create empty F vectors
    FROH.BCFTools.100KB.c = vector(mode = "numeric", length = length(unique(dta[,2])))
    names(FROH.BCFTools.100KB.c) = unique(dta[,2])
    
    FROH.BCFTools.1MB.c = FROH.BCFTools.100KB.c
    
    #Loop through indvs and get F
    for (indv in unique(dta[,2])) {
      
      ## 100KB
      
      #Subset the data
      dtasub = dta100KB[dta100KB$[,2] == indv,]
      #Estimate F
      FROH.BCFTools.100KB.c[names(FROH.BCFTools.100KB.c) == indv] = sum(dtasub[,6])/3.2e9
      
      ## 1MB
      
      #Subset the data
      dtasub = dta1MB[dta1MB[,2] == indv,]
      #Estimate F
      FROH.BCFTools.1MB.c[names(FROH.BCFTools.1MB.c) == indv] = sum(dtasub[,6])/3.2e9
      
    }
    
    #Paste names as other column
    FROH.BCFTools = cbind(sample=names(FROH.BCFTools.100KB.c), FROH.BCFTools.100KB.c = FROH.BCFTools.100KB.c, FROH.BCFTools.1MB.c = FROH.BCFTools.1MB.c)
    
    #as numeric
    FROH.BCFTools[,2] = as.numeric(as.character(FROH.BCFTools[,2]))
    FROH.BCFTools[,3] = as.numeric(as.character(FROH.BCFTools[,3]))
    
    #sample names
    FROH.BCFTools[,1] = sapply(FROH.BCFTools[,1], function(x){strsplit(x, split = "_")[[1]][2]})
    
    #Save the dataframe
    write.table(FROH.BCFTools, "./data/F_datasets/TMP/${CONT}_allchr_both_FROH.BCFTools.size.txt", quote = F, col.names = T, row.names = F)


    ################################
    #### ROHs called with PLINK ####
    ################################

    dta100KB = read.table("./data/CallingROHsPLINK/1KG_${CONT}_100KB.hom.indiv", header = T)
    dta1MB = read.table("./data/CallingROHsPLINK/1KG_${CONT}_100KB.hom.indiv", header = T)
    
    #Create empty F vectors
    FROH.PLINK.100KB.c = vector(mode = "numeric", length = length(unique(dta100KB[,2])))
    names(FROH.PLINK.100KB.c) = unique(dta100KB[,2])
    
    FROH.PLINK.1MB.c = FROH.PLINK.100KB.c
    
    #Loop through indvs and get F
    for (indv in unique(dta100KB[,2])) {
      
      ## 100KB
      
      #Estimate F, divided by e6 becasue in KB
      FROH.PLINK.100KB.c[names(FROH.PLINK.100KB.c) == indv] = dta100KB[,5][dta100KB[,2] == indv]/3.2e6
      
      ## 1MB
      
      #Estimate F
      FROH.PLINK.1MB.c[names(FROH.PLINK.1MB.c) == indv] = dta1MB[,5][dta1MB[,2] == indv]/3.2e6
      
    }
    
    #sanity check
    all.equal(names(FROH.PLINK.100KB.c),names(FROH.PLINK.1MB.c))
    
    #Paste names as other column
    FROH.PLINK = as.data.frame(cbind(sample=names(FROH.PLINK.100KB.c), FROH.PLINK.100KB.c = FROH.PLINK.100KB.c, FROH.PLINK.1MB.c = FROH.PLINK.1MB.c))
    
    #as numeric
    FROH.PLINK[,2] = as.numeric(as.character(FROH.PLINK[,2]))
    FROH.PLINK[,3] = as.numeric(as.character(FROH.PLINK[,3]))
    
    #Save the dataframe
    write.table(FROH.PLINK, "./data/F_datasets/TMP/${CONT}_allchr_both_FROH.PLINK.size.txt", quote = F, col.names = T, row.names = F)

EOF

done

#Merging of Fs from all the continents and with Fs estimated with the entire dataset !
R --vanilla << EOF

    ## Read all the files

    #WORLD estimations
    dtaWORLD = read.table("./data/F_datasets/TMP/WORLD_allchr_FasbothFunibothFROHs.txt", h = T)

    #Continent SNPs-based F estimations
    dtaCONTSNPs = read.table("./data/F_datasets/TMP/PerCONT_allchr_FasbothFuni.txt", h = T)

    #We'll loop through continents for ROHs reading and merging

    #First create the empty dtaCONTBCFTools and dtaCONTPLINK
    dtaCONTBCFTools = as.data.frame(matrix(nrow = 0, ncol = 3))
    colnames(dtaCONTBCFTools) = c("sample", "FROH.BCFTools.100KB.c", "FROH.BCFTools.1MB.c")
    dtaCONTPLINK = as.data.frame(matrix(nrow = 0, ncol = 3))
    colnames(dtaCONTPLINK) = c("sample", "FROH.PLINK.100KB.c", "FROH.PLINK.1MB.c")

    for(cont in unique(dtaWORLD[,3])){

        #read the df BCFTools
        dtaBCFTools = read.table(paste0("./data/F_datasets/TMP/", cont, "_allchr_both_FROH.BCFTools.size.txt"), h = T)
        #rbind the df
        dtaCONTBCFTools = rbind(dtaCONTBCFTools, dtaBCFTools)

        #read the df PLINK
        dtaPLINK = read.table(paste0("./data/F_datasets/TMP/", cont, "_allchr_both_FROH.PLINK.size.txt"), h = T)
        #rbind the df
        dtaCONTPLINK = rbind(dtaCONTPLINK, dtaPLINK)

    }

    #Time to merge everything
    dta = merge(dtaWORLD, dtaCONTSNPs, by = "sample")
    dta.2 = merge(dta, dtaCONTBCFTools, by = "sample")
    dta.3 = merge(dta.2, dtaCONTPLINK, by = "sample")

    #Now we just order it the way I want
    dta.4 = dta.3[,c(1:4,9,5,10,6,11,7,12,8,13)]

    #save the table
    write.table(dta.4, "/data/F_datasets/WORLD_allchr_allFs.txt", quote = F, row.names = T, col.names = T)

    #We'll also save EAS and AFR tables
    write.table(dta.4[dta.4[,3] == "EAS",], "/data/F_datasets/EAS_allchr_allFs.txt", quote = F, row.names = T, col.names = T)
    write.table(dta.4[dta.4[,3] == "AFR",], "/data/F_datasets/AFR_allchr_allFs.txt", quote = F, row.names = T, col.names = T)


EOF

