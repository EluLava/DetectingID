#!/bin/bash

#Create genotypes (with ms) for 1000 founders

mspms 2000 20 -t 4000 -r 4000 100000000 -p 9 > foundersb.txt

R --vanilla << EOF

    library(gaston)
    library(hierfstat)
    library(JGTeach)

    #Simulate the pedigree
    ped.pg<-JGTeach::buildped.2sexes(founders.m=500,founders.f=500,fert=1.0,death.rate=0.5,breed.prop=0.1,n.tstep=30)

    #Read founder genotypes and convert into dosage
    myfounds<-hierfstat::ms2dos("./data/bedmatrices/foundersb.txt")

    #get pedigree F
    Fped.pg<-diag(grm2kinship(pedARM(ped.pg[,2],ped.pg[,3])))

    #Simulate pedigree genotypes (from founders' genotypes)
    bedpg.20.ch1<-drop.along.ped(ped=ped.pg,founders.genotypes=myfounds[,'alldat'][,myfounds$bim$chr==1],maplength=20)

    #Fill the bed matrix
    bedpg.20.ch1<-as.bed.matrix(bedpg.20.ch1[,,1]+bedpg.20.ch1[,,2])
    bedpg.20.ch1@snps[,'pos']<-myfounds$bim$pos[myfounds$bim$chr==1]
    bedpg.20.ch1@snps[,'chr']<-myfounds$bim$chr[myfounds$bim$chr==1]
    bedpg.20.ch1@ped[,'father']<-ped.pg[,'sire']
    bedpg.20.ch1@ped[,'mother']<-ped.pg[,'dam']

    #write the bed matrix
    write.bed.matrix(bedpg.20.ch1, "./data/bedmatrices/PGPED_ch_1_22")

EOF
