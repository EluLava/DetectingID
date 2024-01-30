first ped files have to be converted to bed files using plink --make-bed

## plink --file Genotypes --mak-bed --out Genotypes

#then, in an R session
#load necessary libraries
library(gaston)
library(hierfstat)
library(lme4)
library(lmerTest)

#read genotype data
bed<-read.bed.matrix("Genotypes")
#read phenotype data
pheno<-read.table("Data.txt",header=TRUE)
#check obs are odered the same way
all.equal(pheno$id,bed@ped$id)
#removes non autosomes
bedb<-bed[,bed@snps$chr <32]


#############Note##################
#we are NOT filtering for MAF or LD, 
#as it is unecessary and could lead to biased results
#FuniW and allele sharing GRM handle this well
####################################



# filtering individuals as in original paper
# from file morphology_INLA.R
data_adult_temp <- pheno[which(pheno$laflok %in% c("20","22","23","24","26","27","28","38")),]

###
###Include individuals with morphology information
###

#Only those individuals that have age1mass estimate and have been adults on one of the study islands
data_morph <- data_adult_temp[!is.na(data_adult_temp$age1mass),]

#keep only individuals with morpho infos
bedm<-bedb[match(data_morph$id,bed@ped$id),]
bedm

#Allele sharing matrix
M<-matching(bedm)

#kinship matrix
ks<-beta.dosage(M,MATCHING=TRUE)
# allele sharing F Fas
Fas<-diag(ks)

###get.funiw function ##########################
get.funiw<-function(dos){
  
  if(!class(dos)[[1]]=="bed.matrix") stop("Argument must be of class bed.matrix. Exiting")
  p<-dos@p
  het<-2*p*(1-p)
  res<-apply(gaston:::as.matrix(dos),1,function(x) {
    nas<-which(is.na(x)); 
    xs<-x[-nas];
    ps<-p[-nas];
    hets<-het[-nas];
    sum(xs^2-(1+2*ps)*xs+2*ps^2)/sum(hets)
  }
  )
  return(Funi=unlist(res))
}

##################################################

#Funiw

Funiw<-get.funiw(bedm)




myb<-outer(data_morph$all_hatchyears,data_morph$all_hatchyears,function(x,y) x==y)
#makes random factor island into a matrix for analysis with lmm.aireml
misl<-outer(data_morph$laflok,data_morph$laflok,function(x,y) x==y)

#create nested factor year:island  

yii<-factor(paste(data_morph$all_hatchyears,data_morph$laflok,sep="."))
#makes random factor year:islands into a matrix for analysis with lmm.aireml
myii<-outer(yii,yii,function(x,y) x==y)

#traits to be analysed:
#[13] "age1tarsus"     "age1wing"       "age1billD"      "age1billL"
#[17] "age1mass" 

morpho.data_morphs<-data_morph[,13:17]



TheF<-Funiw

get.ID.eval<-function(Y,TheF){
  tmp.s<-summary(lm(Y~gen_sex+TheF,data=data_morph))
  tmp.a<-anova(lm(Y~gen_sex+TheF,data=data_morph))
  tmp.lmm<-summary(lmer(Y~gen_sex+TheF+(1|laflok)+(1|yii),data=data_morph))
  #following 2 lines check lmm.aireml and lmer give the same results with no GRM
  #str(lmm.aireml(Y,cbind(rep(1,dim(data_morph)[1]),factor(data_morph$gen_sex),TheF),K=list(misl,myii),verbose=FALSE))
  #score.fixed.linear(x=cbind(TheF),Y,cbind(rep(1,dim(data_morph)[1]),factor(data_morph$gen_sex)),K=list(misl,myii),verbose=FALSE)$p
  
  tmp.lmmFull<-lmm.aireml(Y,cbind(rep(1,dim(data_morph)[1]),factor(data_morph$gen_sex),TheF),K=list(misl,myii,kinship2grm(ks)),verbose=FALSE)
  pval.TheF.lmmfull<-score.fixed.linear(x=cbind(TheF),Y,cbind(rep(1,dim(data_morph)[1]),factor(data_morph$gen_sex)),K=list(misl,myii,kinship2grm(ks)),verbose=FALSE)$p
  
  tmp.lmmGRM<-lmm.aireml(Y,cbind(rep(1,dim(data_morph)[1]),factor(data_morph$gen_sex),TheF),K=kinship2grm(ks),verbose=FALSE)
  pval.TheF.lmmGRM<-score.fixed.linear(x=cbind(TheF),Y,cbind(rep(1,dim(data_morph)[1]),factor(data_morph$gen_sex)),K=kinship2grm(ks),verbose=FALSE)$p
  res.lm<-c(tmp.s$coefficients[,1],NA,NA,NA,tmp.a[3,3],tmp.s$coefficients[3,4])
  res.lmm<-with(tmp.lmm,c(coefficients[,1],varcor$laflok,varcor$yii,NA,sigma^2,coefficients[3,5]))
  res.lmmFull<-with(tmp.lmmFull,c(BLUP_beta,tau,sigma2,pval.TheF.lmmfull))
  res.lmmGRM<-with(tmp.lmmGRM,c(BLUP_beta,NA,NA,tau,sigma2,pval.TheF.lmmGRM))
  
  res<-rbind(res.lm,res.lmmGRM,res.lmm,res.lmmFull)
  res<-data.frame(res)
  names(res)<-c("Intercept","Sex","TheF","Visl","Vyear:isl","Va","Ve","PvalTheF")
  return(res)
}


all.res.Funiw<-apply(morpho.data_morphs,2,function(x) get.ID.eval(x,Funiw))

library(xtable)
#results in a printable format
lapply(all.res.Funiw,function(x) round(x,digits=3))
lapply(all.res.Funiw,xtable,digits=3)