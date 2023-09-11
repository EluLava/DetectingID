#Create directory if does not exists
mkdir -p VCFs

#Pass from BED to VCF for ROHs and HBD segments estimation
plink --bfile ./data/bedmatrices/PGPED_ch_1_22 --recode vcf iid bgz --out ./data/VCFs/PGPED

#index
bcftools index ./data/VCFs/PGPED.vcf.gz

##Calling HBD segments with BCFTools

#Create directory if does not exists
mkdir -p CallingROHsBCFTools

#Calling segments
bcftools roh -o ./data/CallingROHsBCFTools/PGPED.txt -O r -e "-" -G 30 ./data/VCFs/PGPED.vcf.gz

##Calling ROHs with PLINK

#Create directory if does not exists
mkdir -p CallingROHsPLINK

#100KB
plink --bfile ./data/bedmatrices/PGPED_ch_1_22 --homozyg --homozyg-kb 100 --out ./data/CallingROHsPLINK/PGPED_100KB
#1MB
plink --bfile ./data/bedmatrices/PGPED_ch_1_22 --homozyg --homozyg-kb 1000 --out ./data/CallingROHsPLINK/PGPED_1MB

#create output directory if does not exists
mkdir -p ./data/F_datasets/TMP

R --vanilla << EOF

    ###################################
    #### ROHs called with BCFTools ####
    ###################################

    dta = read.table("./data/CallingROHsBCFTools/PGPED.txt")
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
      FROH.BCFTools.100KB[names(FROH.BCFTools.100KB) == indv] = sum(dtasub[,6])/3.2e9
      
      ## 1MB
      
      #Subset the data
      dtasub = dta1MB[dta1MB[,2] == indv,]
      #Estimate F
      FROH.BCFTools.1MB[names(FROH.BCFTools.1MB) == indv] = sum(dtasub[,6])/3.2e9
      
    }
    
    #Paste names as other column
    FROH.BCFTools = cbind(sample=names(FROH.BCFTools.100KB), FROH.BCFTools.100KB = FROH.BCFTools.100KB, FROH.BCFTools.1MB = FROH.BCFTools.1MB)
    
    #as numeric
    FROH.BCFTools[,2] = as.numeric(as.character(FROH.BCFTools[,2]))
    FROH.BCFTools[,3] = as.numeric(as.character(FROH.BCFTools[,3]))
    
    #sample names
    FROH.BCFTools[,1] = sapply(FROH.BCFTools[,1], function(x){strsplit(x, split = "_")[[1]][2]})
    
    #Save the dataframe
    write.table(FROH.BCFTools, "./data/F_datasets/TMP/PGPED_allchr_both_FROH.BCFTools.size.txt", quote = F, col.names = T, row.names = F)


    ################################
    #### ROHs called with PLINK ####
    ################################

    dta100KB = read.table("./data/CallingROHsPLINK/PGPED_100KB.hom.indiv", header = T)
    dta1MB = read.table("./data/CallingROHsPLINK/PGPED_100KB.hom.indiv", header = T)
    
    #Create empty F vectors
    FROH.PLINK.100KB = vector(mode = "numeric", length = length(unique(dta100KB[,2])))
    names(FROH.PLINK.100KB) = unique(dta100KB[,2])
    
    FROH.PLINK.1MB = FROH.PLINK.100KB
    
    #Loop through indvs and get F
    for (indv in unique(dta100KB[,2])) {
      
      ## 100KB
      
      #Estimate F, divided by e6 becasue in KB
      FROH.PLINK.100KB[names(FROH.PLINK.100KB) == indv] = dta100KB[,5][dta100KB[,2] == indv]/3.2e6
      
      ## 1MB
      
      #Estimate F
      FROH.PLINK.1MB[names(FROH.PLINK.1MB) == indv] = dta1MB[,5][dta1MB[,2] == indv]/3.2e6
      
    }
    
    #sanity check
    all.equal(names(FROH.PLINK.100KB),names(FROH.PLINK.1MB))
    
    #Paste names as other column
    FROH.PLINK = as.data.frame(cbind(sample=names(FROH.PLINK.100KB), FROH.PLINK.100KB = FROH.PLINK.100KB, FROH.PLINK.1MB = FROH.PLINK.1MB))
    
    #as numeric
    FROH.PLINK[,2] = as.numeric(as.character(FROH.PLINK[,2]))
    FROH.PLINK[,3] = as.numeric(as.character(FROH.PLINK[,3]))
    
    #Save the dataframe
    write.table(FROH.PLINK, "./data/F_datasets/TMP/PGPED_allchr_both_FROH.PLINK.size.txt", quote = F, col.names = T, row.names = F)


    ##############################
    #### FAS, FUNIu and FUNIw ####
    ##############################

    if(!("gaston" %in% installed.packages()[,"Package"])){install.packages("gaston")}
    library(gaston)
    setThreadOptions(5)

    #Load the scripts
    source("./functions/functions.R")

    #Read bed matrix
    bed = read.bed.matrix("./data/bedmatrices/PGPED_ch_1_22")

    #Extract sample name info FROM BED
    sampleBED = bed@ped[,2]

    #Get FuniUN
    FuniUN = get.funiuT(bed, nb.cores = 5)
    #Get FuniWE
    FuniWE = get.funiwT(bed, nb.cores = 5)
    #Get fas
    Fas = get.fasT(bed, nb.cores = 5)

    #merge everyting
    SNPbasedFs = cbind(sample=names(FuniUN), Fas = Fas, Funi.un = FuniUN, Funi.we = FuniWE)

    #Save the table
    write.table(SNPbasedFs, "./data/F_datasets/TMP/PGPED_allchr_FasbothFuni.txt", quote = F, col.names = T, row.names = F)


    #########################################
    #### Merging and sample descriptions ####
    #########################################

    #merge sample desc POP and CONT with our estimations
    dta = merge(as.data.frame(cbind(sample = SNPbasedFs[,1], pop = rep(NA, nrow(SNPbasedFs)), cont = rep(NA, nrow(SNPbasedFs)))), SNPbasedFs, by = "sample")
    #merge with FROH BCFTOOLS and PLINK
    dta.2 = merge(dta, FROH.BCFTools, by = "sample")
    dta.3 = merge(dta.2, FROH.PLINK, by = "sample")

    write.table(dta.3, "./data/F_datasets/PGPED_PED_allchr_allFs.txt", quote = F, col.names = T, row.names = F)

EOF