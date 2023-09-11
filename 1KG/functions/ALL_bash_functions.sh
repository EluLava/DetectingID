#!/bin/bash

#Create the function
LD_scores(){

	#First argument passed to the script is bedfile
	bedF=$1

	#Second argument passed to the script is output file
	output=$2

	#run gcta on the corresponding bedmatrix
	gcta64 --bfile bedF --ld-score-region 200 --ld-wind 200 --out output --thread-num 30

}

#Allows to get the three F  matrices (dataframes) from a Bedmatrix
Stratification_MAF_LD(){

	#First argument passed to the script is bedfile
	bedF=$1

	#ld score file is the second argument passed to the script
	ld_scores_file=$2

	#third argument passed to the script is output file
	output=$3

	R --vanilla << EOF

		#List packages needed
		list.of.packages <- c("gaston", "Rcpp")
		new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
		if(length(new.packages)) install.packages(new.packages)
		#Load gaston anyway
		library(gaston)

		#Source the code containing the function for estimating F
		source("./functions/functions.R")

		#Read the bed matrix
		bed = read.bed.matrix("${bedF}")

		#read the ld score table
		lds_seg = read.table("${ld_scores_file}", header = T)

		#Apply the function
		fs = LD.MAF.bins.getF(bed, lds_seg)

		#Set row.names as samples id
		rownames(fs[[1]]) = bed@ped[,2]
		rownames(fs[[2]]) = bed@ped[,2]
		rownames(fs[[3]]) = bed@ped[,2]

		#Write the three F.LD.MAF matrices

		#Funi Weighted
		write.table(fs[[1]], "${output}.FuniWE", quote = F, col.names = T, row.names = F)
		#Funi Unweighted
		write.table(fs[[2]], "${output}.FuniUN", quote = F, col.names = T, row.names = F)
		#Fas
		write.table(fs[[3]], "${output}.Fas", quote = F, col.names = T, row.names = F)

EOF

}


