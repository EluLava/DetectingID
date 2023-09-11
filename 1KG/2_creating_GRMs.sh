
#Loop through populations
for POP in {WORLD, EAS, AFR}; do

	R --vanilla << EOF

	#List packages needed
	list.of.packages <- c("SNPRelate", "hierfstat", "devtools")
	new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
	if(length(new.packages)){install.packages(new.packages)}
	if(!("JGTeach" %in% installed.packages()[,"Package"])){library(devtools); install_github("jgx65/JGTeach")}

	library(SNPRelate)
	library(hierfstat)
	library(JGTeach)
	
	#From bed to GDS
	snpgdsBED2GDS(bed.fn = "./data/bedmatrices/${POP}_ch_1_22.bed", fam.fn = "./data/bedmatrices/${POP}_ch_1_22.fam",
				   bim.fn = "./data/bedmatrices/${POP}_ch_1_22.bim", out.gdsfn = "./data/bedmatrices/${POP}_ch_1_22.gds")
	
	#Read the gds object
	gdsMat = snpgdsOpen("./data/bedmatrices/${POP}_ch_1_22.gds")
	
	#From GDS to GRM gcta UNweighted
	GCTAun = snpgdsGRM(gdsMat, method = "GCTA")
	#save it
	saveRDS(GCTAun, "./data/GRMs/${POP}_GCTAun_GRM_ch_1_22.RDS")
	
	#From GDS to dosage
	dos = snpgdsGetGeno(gdsMat)
	
	#matching Mat
	allMatchingMat = matching(dos)

	#Estimate pairwise kinships (FROM ALLELES MATCHING)
	kas = beta.dosage(allMatchingMat, MATCHING = TRUE)
	#Pass to GRM matrix
	grmas = kinship2grm(kas)
	#write GRMas
	saveRDS(grmas, "./data/GRMs/${POP}_AS_GRM_ch_1_22.RDS")

	#pairwise kinships (STANDARD) ratio of averages: GCTA Weigthed
	GCTAwe = Kc0(Mall,matching=TRUE)
	#save it
	saveRDS(GCTAwe, "./data/GRMs/${POP}_GCTAwe_GRM_ch_1_22.RDS")
	
EOF


done

