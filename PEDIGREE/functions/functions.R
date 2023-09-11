#################################################################################
###### functions to obtain individual Fs from large (>100k SNPs) data sets ######
#################################################################################

get.funiwT<-function(bed,nloc=100000,nb.cores=30){
#calculates overall Funiw for bed matrix with more than 100k SNPs
if (dim(bed)[2]<=nloc){ 
  get.funiw(bed)$Funi 
  }
else {
   nl<-dim(bed)[2]
   a<-split(1:nl,floor(1:nl/nloc))
   tmp_all<-parallel::mclapply(a,function(x) get.funiw(bed[,x]),mc.cores=nb.cores)
   sdens<-sum(unlist(lapply(tmp_all,function(x) x$het)))

   snums<-vector(length=dim(bed)[1])
   for(i in 1:length(tmp_all)){
     snums<-snums+tmp_all[[i]][[2]]*tmp_all[[i]][[1]]
   }
   snums/sdens
  }
}


#################################################

get.funiuT<-function(bed,nloc=100000,nb.cores=30){
#calculates overall Funiu for bed matrix with more than 100k SNPs
nl<-dim(bed)[2]
ni<-dim(bed)[1]
if (nl<=nloc){ 
  get.funiu(bed)$Funi 
  }
else {
  a<-split(1:nl,floor(1:nl/nloc))
  tmp_all<-parallel::mclapply(a,function(x) get.funiu(bed[,x]),mc.cores=nb.cores)
  snums<-vector(length=ni)
  sdens<-vector(length=ni)

  for(i in 1:length(tmp_all)) {
    sdens<-sdens+tmp_all[[i]]$nloc; 
	snums<-snums+tmp_all[[i]]$nloc*tmp_all[[i]]$Funi
	}
  snums/sdens
}
}

#################################################

get.fasT<-function(bed,nloc=100000,nb.cores=30){
#calculates overall Fas for bed matrix with more than 100k SNPs
matching<-function(dos){
    dos <- gaston::as.matrix(dos)
	NAs<-sum(is.na(dos))
    if (NAs > 0) {
        na <- matrix(rep(1, prod(dim(dos))), ncol = ncol(dos))
        ina <- which(is.na(dos))
        na[ina] <- 0
        dos[ina] <- 1
		mNAs<-tcrossprod(na)
        Mij <- 1/2 * (1 + 1/mNAs * tcrossprod(dos -
            1))
    }
    else {
        nl <- dim(dos)[2]
        Mij <- 1/2 * (1 + tcrossprod(dos - 1)/nl)
    }
	if (NAs==0) nl<-rep(nl,ni) else nl<-diag(mNAs)
    return(list(Mii=diag(Mij),MB=mean(hierfstat::mat2vec(Mij)),nl=nl))
}

if (class(bed)[[1]] != "bed.matrix") stop("bed must be of class bed.matrix. Exiting")

if (dim(bed)[2]<=nloc){ 
  diag(hierfstat::beta.dosage(bed)) 
  }
else {
   nl<-dim(bed)[2]
   ni<-dim(bed)[1]
   a<-split(1:nl,floor(1:nl/nloc))

   tmp_all<-parallel::mclapply(a,function(x) matching(bed[,x]),mc.cores=nb.cores)


   #self-allele sharing
   MiiT<-rep(0,ni)
   MbT<-0
   nlt<-rep(0,ni)
   for (i in 1:length(tmp_all)) {
     MiiT<-MiiT+tmp_all[[i]]$Mii*tmp_all[[i]]$nl
     MbT<-MbT+tmp_all[[i]]$MB*length(a[[i]])
     nlt<-nlt+tmp_all[[i]]$nl
}

   MiiT<-MiiT/nlt
   MbT<-MbT/sum(unlist(lapply(a,length)))

  #Fas inbreeding coeff
  ((MiiT*2-1)-MbT)/(1-MbT)
}
}


#################################################

#Return a list containing (in order) Funi W, Funi U, Fas
LD.MAF.bins.getF<-function(bed, lds_seg){

	#Extract quartiles from LD scores
	quartiles=summary(lds_seg$ldscore_SNP)

	#which SNP in which quartile
	lb1 = which(lds_seg$ldscore_SNP <= quartiles[2])
	lb2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3])
	lb3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP <= quartiles[5])
	lb4 = which(lds_seg$ldscore_SNP > quartiles[5])

	#Subset these SNPs
	LD1 = lds_seg$SNP[lb1]
	LD2 = lds_seg$SNP[lb2]
	LD3 = lds_seg$SNP[lb3]
	LD4 = lds_seg$SNP[lb4]

	#From the bed file, extract which snps are in MAF 7 category and in all LD categories
	MAF7LD1<-which(bed@snps$id %in% LD1 & bed@snps$maf >0.4 & bed@snps$maf <=0.5)
	MAF7LD2<-which(bed@snps$id %in% LD2 & bed@snps$maf >0.4 & bed@snps$maf <=0.5)
	MAF7LD3<-which(bed@snps$id %in% LD3 & bed@snps$maf >0.4 & bed@snps$maf <=0.5)
	MAF7LD4<-which(bed@snps$id %in% LD4 & bed@snps$maf >0.4 & bed@snps$maf <=0.5)

	#From the bed file, extract which snps are in MAF 6 category and in all LD categories
	MAF6LD1<-which(bed@snps$id %in% LD1 & bed@snps$maf >0.3 & bed@snps$maf <=0.4)
	MAF6LD2<-which(bed@snps$id %in% LD2 & bed@snps$maf >0.3 & bed@snps$maf <=0.4)
	MAF6LD3<-which(bed@snps$id %in% LD3 & bed@snps$maf >0.3 & bed@snps$maf <=0.4)
	MAF6LD4<-which(bed@snps$id %in% LD4 & bed@snps$maf >0.3 & bed@snps$maf <=0.4)

	#From the bed file, extract which snps are in MAF 5 category and in all LD categories
	MAF5LD1<-which(bed@snps$id %in% LD1 & bed@snps$maf >0.2 & bed@snps$maf <=0.3)
	MAF5LD2<-which(bed@snps$id %in% LD2 & bed@snps$maf >0.2 & bed@snps$maf <=0.3)
	MAF5LD3<-which(bed@snps$id %in% LD3 & bed@snps$maf >0.2 & bed@snps$maf <=0.3)
	MAF5LD4<-which(bed@snps$id %in% LD4 & bed@snps$maf >0.2 & bed@snps$maf <=0.3)

	#From the bed file, extract which snps are in MAF 4 category and in all LD categories
	MAF4LD1<-which(bed@snps$id %in% LD1 & bed@snps$maf >0.1 & bed@snps$maf <=0.2)
	MAF4LD2<-which(bed@snps$id %in% LD2 & bed@snps$maf >0.1 & bed@snps$maf <=0.2)
	MAF4LD3<-which(bed@snps$id %in% LD3 & bed@snps$maf >0.1 & bed@snps$maf <=0.2)
	MAF4LD4<-which(bed@snps$id %in% LD4 & bed@snps$maf >0.1 & bed@snps$maf <=0.2)

	#From the bed file, extract which snps are in MAF 2 category and in all LD categories
	MAF3LD1<-which(bed@snps$id %in% LD1 & bed@snps$maf >0.01 & bed@snps$maf <=0.1)
	MAF3LD2<-which(bed@snps$id %in% LD2 & bed@snps$maf >0.01 & bed@snps$maf <=0.1)
	MAF3LD3<-which(bed@snps$id %in% LD3 & bed@snps$maf >0.01 & bed@snps$maf <=0.1)
	MAF3LD4<-which(bed@snps$id %in% LD4 & bed@snps$maf >0.01 & bed@snps$maf <=0.1)

	#From the bed file, extract which snps are in MAF 2 category and in all LD categories
	MAF2LD1<-which(bed@snps$id %in% LD1 &  bed@snps$maf <=0.01)
	MAF2LD2<-which(bed@snps$id %in% LD2 &  bed@snps$maf <=0.01)
	MAF2LD3<-which(bed@snps$id %in% LD3 &  bed@snps$maf <=0.01)
	MAF2LD4<-which(bed@snps$id %in% LD4 &  bed@snps$maf <=0.01)

	#List with indexes for each category
	aa<-list(MAF7LD1,MAF7LD2,MAF7LD3,MAF7LD4,
	MAF6LD1,MAF6LD2,MAF6LD3,MAF6LD4,
	MAF5LD1,MAF5LD2,MAF5LD3,MAF5LD4,
	MAF4LD1,MAF4LD2,MAF4LD3,MAF4LD4,
	MAF3LD1,MAF3LD2,MAF3LD3,MAF3LD4,
	MAF2LD1,MAF2LD2,MAF2LD3,MAF2LD4)

	#Get Funi W from all variants per MAF.LD bins
	funiw_MAF <- lapply(aa,function(x) get.funiwT(bed[,x]))
	
	#Get Funi U from all variants per MAF.LD bins
	funiu_MAF <- lapply(aa,function(x) get.funiuT(bed[,x]))
		
	#Get Fas from all variants per MAF.LD bins
	fas_MAF <- lapply(aa,function(x) get.fasT(bed[,x]))

	#Pass the lists to dataframes format
	FuniwMAFxLDx<-data.frame(matrix(unlist(funiw_MAF),ncol=24))
	FuniuMAFxLDx<-data.frame(matrix(unlist(funiu_MAF),ncol=24))
	FasMAFxLDx<-data.frame(matrix(unlist(fas_MAF),ncol=24))

	#set names of matrix Funi Weigthed
	dimnames(FuniwMAFxLDx)[[2]]<-c("MAF7LD1","MAF7LD2","MAF7LD3","MAF7LD4",
	"MAF6LD1","MAF6LD2","MAF6LD3","MAF6LD4",
	"MAF5LD1","MAF5LD2","MAF5LD3","MAF5LD4",
	"MAF4LD1","MAF4LD2","MAF4LD3","MAF4LD4",
	"MAF3LD1","MAF3LD2","MAF3LD3","MAF3LD4",
	"MAF2LD1","MAF2LD2","MAF2LD3","MAF2LD4")

	#set names of matrix Funi UnWeigthed
	dimnames(FuniuMAFxLDx)[[2]]<-c("MAF7LD1","MAF7LD2","MAF7LD3","MAF7LD4",
	"MAF6LD1","MAF6LD2","MAF6LD3","MAF6LD4",
	"MAF5LD1","MAF5LD2","MAF5LD3","MAF5LD4",
	"MAF4LD1","MAF4LD2","MAF4LD3","MAF4LD4",
	"MAF3LD1","MAF3LD2","MAF3LD3","MAF3LD4",
	"MAF2LD1","MAF2LD2","MAF2LD3","MAF2LD4")
	
	#set names of matrix Fas
	dimnames(FasMAFxLDx)[[2]]<-c("MAF7LD1","MAF7LD2","MAF7LD3","MAF7LD4",
	"MAF6LD1","MAF6LD2","MAF6LD3","MAF6LD4",
	"MAF5LD1","MAF5LD2","MAF5LD3","MAF5LD4",
	"MAF4LD1","MAF4LD2","MAF4LD3","MAF4LD4",
	"MAF3LD1","MAF3LD2","MAF3LD3","MAF3LD4",
	"MAF2LD1","MAF2LD2","MAF2LD3","MAF2LD4")

	output = list(FuniwMAFxLDx, FuniuMAFxLDx, FasMAFxLDx)
	names(output) = c("Funi.we.perMAFLD.bins","Funi.un.perMAFLD.bins","Fas.perMAFLD.bins")
	
	return(output)

}

##################################################

id.est<-function(fs,binsfas=NULL,binsfuniu=NULL,binsfuniw=NULL,genos,filter.maf=0.05,n.causal=1000,b=3,s=1,H2=0.8,n.rep=1000,
                 selh=TRUE,selp=TRUE,direct=FALSE,inds=NULL,b1=1,b2=1,debug=FALSE,...){
#fs data frame containing inbreeding coefficients
#binsfas contains fas inbreeding coefficients per maf (and LD) bins
#binsfuniu contains funi_u inbreeding coefficients per maf (and LD) bins
#binsfuniw contains funi_w inbreeding coefficients per maf (and LD) bins
#genos is a bed (gaston) object
#n.causal is number of causal loci
#b is the overall directional dominance effect. 
#s is the overall directional additive effect
#H2 is the intended heritability of the trait
#n.rep is the number of traits (replicates) to simulate
#selh=TRUE implies dominance effect inversely prop to heterozygosity [2p(1-p)]
#selp=TRUE implies additive effect inversely prop to MAF
#dir=TRUE implies additive effects of Minor alleles all positives (DEMA)
#inds is a selection of individuals to use
#b1 and b2 are the parameters of the beta distribution
#from which additive and dominance effects are drawn
#b1,b2=1-> uniform dist
#b1=b2 >10 approx normal, mu=0.5
#debug opens a browser in the middle of the execution
#if(!isa(genos,"bed.matrix")) {stop("genos must be a bed.matrix object. Exiting")} #R 2.1.4
if(class(genos)[[1]]!="bed.matrix") {stop("genos must be a bed.matrix object. Exiting")} # R < 2.1.4? 

if (H2>0.999999) H2i1<-TRUE else H2i1<-FALSE
nfs<-dim(fs)[2]
#nbinfs<-dim(binfs)[2]
my.res<-matrix(numeric(n.rep*((nfs+1+6)+6)),ncol=((nfs+1+6)+6))

#subset of the data
if (!is.null(inds)){
genos<-genos[inds,]
fs<-fs[inds,] # or recalculate it?
}

#if(isa(genos,"bed.matrix")){
maf<-genos@snps$maf
#}
#else{
#af<-colMeans(genos)/2
#maf<-af
#maf[af>0.5]<-1-maf[af>0.5]
#}
candidates<-which(maf>=filter.maf)

####
#for saving phenotypes and effect sizes
#
ni<-dim(genos)[1]
phenos<-matrix(numeric(ni*n.rep),ncol=n.rep)
add.ef<-matrix(numeric(n.causal*n.rep),ncol=n.rep)
dom.ef<-matrix(numeric(n.causal*n.rep),ncol=n.rep)
causal.id<-matrix(character(n.causal*n.rep),ncol=n.rep)

for (i in 1:n.rep){
causal<-sample(candidates,size=n.causal)

X<-genos[,causal]
afc<-X@p #genos@p[causal] #X@snps$maf
mafc<-X@snps$maf
causal.id[,i]<-X@snps$id
X<-as.matrix(X)

#make sure that EF are for Minor allele count
switches<-which(mafc<afc)
X[,switches]<-2-X[,switches]
afc<-colMeans(X)/2

#get fs (fas fhom and funi) from causal genes only
causal.fs<-get.fs(X,levs=0)

#if (direct | sel) x<-order(af[causal])

#draw as from beta dist and insure \sum a=s
h<-2*afc*(1-afc)

a<-rbeta(n.causal,b1,b2);a<-a/sum(a)*s

if(selp) a<-a/2/afc
else a<-a/2/mean(afc)

#larger eff sizes go to lower allele freq
#what's below is not really necessary
#(effect sizes inversely correlated with allele freq)
# but adds a bit of realism
#sa<-sort(a,decreasing=TRUE);
#x<-order(afc);
#a[x]<-sa

#alternate plus and minus if no DEMA
if(!direct) a<-a*(-1)^(1:n.causal)



#draw ds from uniform
if (b!=0) {
d<-rbeta(n.causal,b1,b2);d<-d/sum(d)


if (selh) d<-(d/h)*b
else d<-d*b/mean(h) #d/sum(d)*b/mean(h)

         }
else {d<-rep(0,n.causal)}


#cat(sum(abs(a)),sum(d),"\n") #for debugging
ya<-(2-X)%*%a #+effect on the reference allele if dema
yb<-(X*(2-X))%*%d


if (debug) browser()


yg<-ya+yb
if (H2i1) ye<-0 else ye<-rnorm(length(yg),sd=(sd(ya)^2*(1/H2-1))^.5)

#traits sds scale to 1

y<-yg+ye

#one may want to standardize things
#ya<-ya/sd(y);yb<-yb/sd(y);yg<-yg/sd(y);y<-y/sd(y);

phenos[,i]<-y
add.ef[,i]<-a
dom.ef[,i]<-d

tfs<-data.frame(causal.fs,fs)
vars<-apply(tfs,2,var)


if (!is.null(binsfas)) {
multreg.model.binFas<-summary(lm(as.formula(paste("y ~ ", paste(colnames(binsfas), collapse= "+"))),data=binsfas))
sum.bs.binFas<-sum(multreg.model.binFas$coeff[-1,1])
adj.r.rsquared.binFas<-multreg.model.binFas$adj.r.squared
}
else {
sum.bs.binFas<-NA
adj.r.rsquared.binFas<-NA
}

if (!is.null(binsfuniu)) {
multreg.model.binFuniu<-summary(lm(as.formula(paste("y ~ ", paste(colnames(binsfuniu), collapse= "+"))),data=binsfuniu))
sum.bs.binFuniu<-sum(multreg.model.binFuniu$coeff[-1,1])
adj.r.rsquared.binFuniu<-multreg.model.binFuniu$adj.r.squared
}
else {
sum.bs.binFuniu<-NA
adj.r.rsquared.binFuniu<-NA
}

if (!is.null(binsfuniw)) {
multreg.model.binFuniw<-summary(lm(as.formula(paste("y ~ ", paste(colnames(binsfuniw), collapse= "+"))),data=binsfuniw))
sum.bs.binFuniw<-sum(multreg.model.binFuniw$coeff[-1,1])
adj.r.rsquared.binFuniw<-multreg.model.binFuniw$adj.r.squared
}
else {
sum.bs.binFuniw<-NA
adj.r.rsquared.binFuniw<-NA
}



my.res[i,]<-c(cov(cbind(y,tfs))[1,-1]/vars,sum.bs.binFas,sum.bs.binFuniu,sum.bs.binFuniw,
sd(y),sd(ya),sd(yb),sd(yg),mean(h),cor(a,d))
}
my.res<-data.frame(my.res)
names(my.res)<-c("c.Fas","c.Fhom","c.Funi_u","c.Funi_w",names(fs),
"cum.bin.Fas","cum.bin.Funi_u","cum.bin.Funi_w","sd.trait","sd.BV","sd.Dom","sd.G","h","cor.a.d")
my.res.sum<-apply(my.res,2,function(x) {mx<-mean(x);sex<-sd(x)/length(x)^.5;c(mx,sex)})
return(list(call=match.call(),phenos=phenos,causal=causal.id,add.ef=add.ef,dom.ef=dom.ef,est=my.res,summary=my.res.sum))
}


#################################################################################

