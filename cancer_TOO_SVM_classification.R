#!/usr/bin/env Rscript
#cancer classification include normals

#Get significant probes in all cancer types----------------------

rm(list=ls())


#Get significant probes in all cancer types----------------------
resfolder="/fh/fast/dai_j/CancerGenomics/prostate_methylation/result/"
#find significant probes in each cancer type
#Cancertypes=c("TCGA-PRAD","TCGA-COAD","TCGA-LUAD","TCGA-LUSC","TCGA-BRCA","Krause-EAC")
#remove EAC
Cancertypes=c("TCGA-PRAD","TCGA-COAD","TCGA-LUAD","TCGA-LUSC","TCGA-BRCA","TCGA-LIHC")
#a list to keep all sig probes
sigprobes=NULL
for (i in 1:length(Cancertypes))
{
  Cancer.type=Cancertypes[i]
  Res.file=paste0(resfolder,Cancer.type,"_output.RData")
  load(Res.file)
  bdiff <- output[,6]-output[,4]
  tmp=rownames(output)[output[,23]<0.05/nrow(output) & bdiff>0.05 & output[,4]<0.1 & output[,7]>0.05]
  tmp=tmp[grepl("^cg",tmp)] #remove non-CpG targeting probes ch.xx.xx
  sigprobes=c(sigprobes,list(probes=tmp))
  names(sigprobes)[length(sigprobes)]=Cancer.type
}

#comprobes is the vector to keep all sig probes
comprobes=NULL
#uniqprobes to keep unique sig probes in a cancer
uniqprobes=NULL
for (i in 1:length(sigprobes))
{
  comprobe=sigprobes[[i]]
  if (i==1)
  {
    comprobes=comprobe
  }
  else
  {
    comprobes=unique(c(comprobes,comprobe))
  }
  idx=1:length(sigprobes)
  idx=idx[idx!=i]
  tmp=sigprobes[[i]]
  allprobes=NULL
  for (j in idx)
  {
    allprobes=unique(c(allprobes,sigprobes[[j]]))
  }
  uniqprobes=c(uniqprobes,list(probe=tmp[!tmp %in% allprobes]))
  names(uniqprobes)[i]=Cancertypes[i]
}
length(comprobes) #39054
for (i in 1:length(sigprobes))
{
  print(paste0(Cancertypes[i],":",length(sigprobes[[i]]),",",length(uniqprobes[[i]])))
}

alluniqprobes=NULL
for (i in 1:length(sigprobes))
{
  tmp=uniqprobes[[i]]
  alluniqprobes=unique(c(alluniqprobes,tmp))
}
length(alluniqprobes) #[1] 16760

#get all methylation data and alltype factor
alltumor=tumortype=allnormal=normaltype=NULL
for (i in 1:length(Cancertypes))
{
  Cancer.type=Cancertypes[i]
  if (i>1) rm(tcga_normal_methy,tcga_tumor_methy,tcga_tumor_methy_all)
  if (Cancer.type=="TCGA-PRAD")
  {
    load("/fh/fast/dai_j/CancerGenomics/prostate_methylation/tcga_methy.RData")
    #load("tcga_clinical.RData")
  }else
  {
    RData.file=paste("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/",Cancer.type, "_", "methy.RData", sep='')
    load(RData.file)
  }
  
  if (!exists("tcga_normal_methy"))
  {
    tcga_normal_methy=normal_methy
    tcga_tumor_methy=tumor_methy
    tcga_tumor_methy_all=tumor_methy_all
  }
  tumortype=c(tumortype,rep(Cancer.type,ncol(tcga_tumor_methy_all)))
  normaltype=c(normaltype,rep(Cancer.type,ncol(tcga_normal_methy)))
  if (i==1)
  {
    alltumor=tcga_tumor_methy_all
    allnormal=tcga_normal_methy
  }
  else
  {
    tmp=rownames(tcga_tumor_methy_all)[rownames(tcga_tumor_methy_all) %in% rownames(alltumor)]
    idx1=match(tmp,rownames(alltumor))
    idx2=match(tmp,rownames(tcga_tumor_methy_all))
    alltumor=alltumor[idx1,]
    alltumor=cbind(alltumor,tcga_tumor_methy_all[idx2,])
    allnormal=allnormal[idx1,]
    allnormal=cbind(allnormal,tcga_normal_methy[idx2,])
  }
}
all(comprobes %in% rownames(alltumor)) #F
sum(alluniqprobes %in% rownames(alltumor)) #14719
sum(comprobes %in% rownames(alltumor)) #33733
tmp=rowMeans(alltumor)
idx=which(is.na(tmp))
length(idx) #8404, remove probes with missing values (some may correlated to specific cancer)
alltumor=alltumor[-idx,]
sum(comprobes %in% rownames(alltumor)) #32788
tumortype0=tumortype
tumortype=factor(tumortype,levels = Cancertypes)




library(glmnet)
x=alltumor[rownames(alltumor) %in% comprobes,]
#remove NA in x,row probe column sample
remove_NA_withmean=function(dat)
{
  row.means=rowMeans(dat,na.rm=T)
  row.means[is.na(row.means)]=0 #if a probe has all NAs, set 0
  for (i in 1:ncol(dat))
  {
    idx=which(is.na(dat[,i]))
    if (length(idx)>0)
    {
      dat[idx,i]=row.means[idx] #use mean probe-wize
    }
  }
  return(dat)
}
x=remove_NA_withmean(dat=x)
y=tumortype
#only do it for TCGA tumors
idx =as.character(y)!="Krause-EAC"
x=x[,idx]
y=y[idx]
y=droplevels(y)
#Training data
TCGAtraindat=t(x)


#remember to set #export OMP_NUM_THREADS=1 when use parallel
require(doMC)
registerDoMC(cores=12)
set.seed(1000)
Sys.time()
# library(glmnetUtils)

#Use standardize=T will select probes also with small beta values 
cvfit=cv.glmnet(TCGAtraindat, y, family="multinomial", type.multinomial = "grouped", standardize=F,parallel = T,nfolds=10, trace.it=1)
#Or try without stand
#cvfit=cv.glmnet(TCGAtraindat, y, family="multinomial", type.multinomial = "grouped", standardize=F,parallel = T,nfolds=10, trace.it=1)
plot(cvfit)
#save(cvfit,file="result/cvfit.RData")
Sys.time()
# # use customized lambda,add small lambda, slow
# lambda=10^seq(0,-5,length=100)
# lambda=10^seq(-0.3,-10,length=100) #cvfit10 Convergence for 94th lambda value not reached after maxit=100000 iterations; solutions for larger lambdas returned
# lambda=10^seq(-0.3,-20,length=100) #cvfit20 Convergence for 47th lambda value not reached after maxit=100000 iterations; solutions for larger lambdas returned 
# set.seed(1000)
# cvfit=cv.glmnet(t(x), y, family="multinomial", type.multinomial = "grouped", parallel = T,standardize=T,
#                 lambda=lambda,trace.it=1)
# png(paste0(resfolder,"Allcancers_glmprobes_cvfit10_stand.png.png"),width = 480, height = 480,type = "cairo")
# save(cvfit,file="result/cvfit20_stand.RData")
# plot(cvfit)
# dev.off()
lambda.best=cvfit$lambda.min #147 probes, 127 overlaps with 1se rule
#lambda.best=cvfit$lambda.1se #133 probes

glmcoeff1=coef(cvfit,s=lambda.best)[1]
idx=glmcoeff1$`TCGA-PRAD`@i
#the selected probes
selectprobes1=glmcoeff1$`TCGA-PRAD`@Dimnames[[1]][idx+1]
selectprobes1=selectprobes1[selectprobes1!="(Intercept)"]
glmcoeff5=coef(cvfit,s=lambda.best)[5]
idx=glmcoeff5$`TCGA-BRCA`@i
selectprobes5=glmcoeff5$`TCGA-BRCA`@Dimnames[[1]][idx+1]
selectprobes5=selectprobes5[selectprobes5!="(Intercept)"]
all(selectprobes1==selectprobes5) #T used to verify
length(selectprobes1)  #167



resps=rep(NA,nrow(TCGAtraindat))
resprob=NULL
for (i in 1:nrow(TCGAtraindat))
{
  tmp=predict(cvfit, s=lambda.best,newx = matrix(TCGAtraindat[i,], nrow=1), type = "response")
  resprob=rbind(resprob,as.numeric(tmp))
  resps[i]=dimnames(tmp)[[2]][which.max(as.numeric(tmp))]
}
colnames(resprob)=dimnames(tmp)[[2]]
table(y,resps)



#resps
#y           TCGA-BRCA TCGA-COAD TCGA-LIHC TCGA-LUAD TCGA-LUSC TCGA-PRAD
#TCGA-PRAD         0         0         0         0         0       498
#TCGA-COAD         0       296         0         0         0         0
#TCGA-LUAD         0         0         0       453         5         0
#TCGA-LUSC         0         0         0         9       361         0
#TCGA-BRCA       782         0         0         0         0         0
#TCGA-LIHC         1         0       376         0         0         0


### run svm model using lasso selected probes ###

TCGAtraindat1 <- TCGAtraindat[,colnames(TCGAtraindat) %in% selectprobes1]


library(e1071)

tune.out=tune(svm ,TCGAtraindat1,y,scale=F,  kernel ="radial", ranges =list(cost=c(0.1 ,1 ,10),gamma=c(0.01,0.1,0.5) ))

summary(tune.out)
svmfit <- svm(TCGAtraindat1,y,scale=F,kernel ="radial", cost=10,gamma = 0.1)
pred <- fitted(svmfit)
table(y,pred)


#pred
#y           TCGA-PRAD TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-BRCA TCGA-LIHC
#TCGA-PRAD       498         0         0         0         0         0
#TCGA-COAD         0       296         0         0         0         0
#TCGA-LUAD         0         0       458         0         0         0
#TCGA-LUSC         0         0         0       370         0         0
#TCGA-BRCA         0         0         0         0       782         0
#TCGA-LIHC         0         0         0         0         0       377







remove_NA_withmean=function(dat)
{
  row.means=rowMeans(dat,na.rm=T)
  row.means[is.na(row.means)]=0 #if a probe has all NAs, set 0
  for (i in 1:ncol(dat))
  {
    idx=which(is.na(dat[,i]))
    if (length(idx)>0)
    {
      dat[idx,i]=row.means[idx] #use mean probe-wize
    }
  }
  return(dat)
}
form_xpredict=function(dat1=TCGAtraindat,dat2=t(Hutch_PRAD),selectprobes=selectprobes1)
{
  res=matrix(NA,nrow=nrow(dat2),ncol=ncol(dat1))
  colnames(res)=colnames(dat1)
  rownames(res)=rownames(dat2)
  avaiprobes=intersect(colnames(dat1),colnames(dat2))
  print(paste0(sum(!selectprobes %in% avaiprobes)," probes not available!"))
  print(selectprobes[!selectprobes %in% avaiprobes])
  idx2=match(avaiprobes,colnames(dat2))
  idx1=match(avaiprobes,colnames(res))
  res[,idx1]=dat2[,idx2]
  #fit NA with mean of probe, may cause problem for linear model prediction
  res=t(remove_NA_withmean(dat=t(res)))
  return(res)
}

predit_svm_cancertype=function(dat1=TCGAtraindat,dat2=t(Hutch_PRAD),selectprobes=selectprobes1)
{
  x=form_xpredict(dat1=dat1,dat2=dat2,selectprobes=selectprobes1)
  x <- x[,colnames(x) %in% selectprobes1]
  predictcancer=predict(svmfit, x)
  return(predictcancer)
}


#Hutch_PRAD_predictprob=predict(cvfit, s=lambda.best,newx = Hutch_PRAD_x, type = "response")
#linear model to predict cancer type
predit_glmnet_cancertype=function(allprobs=Hutch_PRAD_predictprob)
{
  if (length(dimnames(allprobs))==3) #[,,1]
  {
    allprobs=allprobs[,,1,drop=T]
  }
  res=prob=rep(NA,nrow(allprobs))
  names(res)=rownames(allprobs)
  for (i in 1:nrow(allprobs))
  {
    tmp=sum(!is.na(allprobs[i,])) #all NA, check it
    if (tmp>0)
    {
      res[i]=colnames(allprobs)[which.max(allprobs[i,])]
      prob[i]=max(allprobs[i,])
    }
  }
  
  return(list(class=res,prob=prob))
}






#Hutch PRAD data
library(data.table)
Hutch_PRAD=fread("/fh/fast/stanford_j/Illumina/Data_Cleaned/Clean_Methylation/beta_main_HP.csv")
Hutch_PRAD=as.data.frame(Hutch_PRAD)
rownames(Hutch_PRAD)=Hutch_PRAD[,1]
Hutch_PRAD=Hutch_PRAD[,2:ncol(Hutch_PRAD)]
Hutch_PRAD_x=form_xpredict() 

Hutch_PRAD_predictcancer=predit_svm_cancertype()
table(Hutch_PRAD_predictcancer)

Hutch_PRAD_x=form_xpredict(dat2=t(Hutch_PRAD))
Hutch_PRAD_predictprob=predict(cvfit, s=lambda.best,newx = Hutch_PRAD_x, type = "response")
Hutch_PRAD_predictcancer=predit_glmnet_cancertype(allprobs=Hutch_PRAD_predictprob)
table(Hutch_PRAD_predictcancer$class)

#TCGA-LUSC TCGA-PRAD 
#1       524






load("/fh/fast/dai_j/CancerGenomics/Colorectal_Cancer/Data/GSE48684.RData")
CancerID <- row.names(clinicaltable)[clinicaltable$status=="cancer"]
AdenoID <- row.names(clinicaltable)[clinicaltable$status=="adenoma"]
NormalID <- row.names(clinicaltable)[clinicaltable$status=="normal-H"]
Grady_COAD_NORM <- MEall[,names(MEall)%in% NormalID]
Grady_COAD <- MEall[,names(MEall)%in% CancerID]

Grady_COAD_predictcancer=predit_svm_cancertype(dat2=t(Grady_COAD))
table(Grady_COAD_predictcancer)

Grady_COAD_x=form_xpredict(dat2=t(Grady_COAD))
Grady_COAD_predictprob=predict(cvfit, s=lambda.best,newx = Grady_COAD_x, type = "response")
Grady_COAD_predictcancer=predit_glmnet_cancertype(allprobs=Grady_COAD_predictprob)
table(Grady_COAD_predictcancer$class)

#TCGA-BRCA TCGA-COAD TCGA-LUAD 
#1        59         4 

GEO_BRCA=as.data.frame(fread("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/BRCA/GSE75067/GSE75067_series_matrix.txt",skip=58))
rownames(GEO_BRCA)=GEO_BRCA$ID_REF
if ("ID_REF" %in% colnames(GEO_BRCA)) GEO_BRCA=GEO_BRCA[,-1]
GEO_BRCA_predictcancer=predit_svm_cancertype(dat2=t(GEO_BRCA))
table(GEO_BRCA_predictcancer,useNA="ifany")



GEO_BRCA_x=form_xpredict(dat2=t(GEO_BRCA))
GEO_BRCA_predictprob=predict(cvfit, s=lambda.best,newx = GEO_BRCA_x, type = "response")
GEO_BRCA_predictcancer=predit_glmnet_cancertype(allprobs=GEO_BRCA_predictprob)
table(GEO_BRCA_predictcancer$class)

#TCGA-BRCA 
#188 



#Try lung cancers
GEO_LUNG=as.data.frame(fread("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/LUNG/GSE56044/GSE56044_series_matrix.txt",skip=69))
rownames(GEO_LUNG)=GEO_LUNG$ID_REF
if ("ID_REF" %in% colnames(GEO_LUNG)) GEO_LUNG=GEO_LUNG[,-1]
sampletable=readxl::read_xls("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/LUNG/GSE56044/GSE56044_GEO_Annotations.xls",sheet = 2)
sampleidtable=readxl::read_xls("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/LUNG/GSE56044/GSE56044_sampleID.xls")
idx=match(sampletable$SampleID,sampleidtable$ID)
sampletable$ID=sampleidtable$SapleID[idx]
GEO_LUAD=GEO_LUNG[,colnames(GEO_LUNG) %in% sampletable$ID[sampletable$Histology=="AC"]]
GEO_LUSC=GEO_LUNG[,colnames(GEO_LUNG) %in% sampletable$ID[sampletable$Histology=="SqCC"]]


GEO_LUAD_predictcancer=predit_svm_cancertype(dat2=t(GEO_LUAD))
table(GEO_LUAD_predictcancer)


GEO_LUAD_x=form_xpredict(dat2=t(GEO_LUAD))
GEO_LUAD_predictprob=predict(cvfit, s=lambda.best,newx = GEO_LUAD_x, type = "response")
GEO_LUAD_predictcancer=predit_glmnet_cancertype(allprobs=GEO_LUAD_predictprob)
table(GEO_LUAD_predictcancer$class)

#TCGA-LUAD TCGA-LUSC 
#79         4 



GEO_LUSC_predictcancer=predit_svm_cancertype(dat2=t(GEO_LUSC))
table(GEO_LUSC_predictcancer)


GEO_LUSC_x=form_xpredict(dat2=t(GEO_LUSC))
GEO_LUSC_predictprob=predict(cvfit, s=lambda.best,newx = GEO_LUSC_x, type = "response")
GEO_LUSC_predictcancer=predit_glmnet_cancertype(allprobs=GEO_LUSC_predictprob)
table(GEO_LUSC_predictcancer$class)


#TCGA-LUAD TCGA-LUSC 
#1        22












#GEO_COAD
GEO_COAD=as.data.frame(fread("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/COAD/GSE77718/GSE77718_series_matrix.txt",skip=65))
rownames(GEO_COAD)=GEO_COAD$ID_REF
GEO_COAD=GEO_COAD[,-1]
sampletable=read.table("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/COAD/GSE77718/sampletable.txt")
sampletable$type="Tumor"
idx=grepl("N",sampletable$V7)
sampletable$type[idx]="Normal"
idx=match(colnames(GEO_COAD),sampletable$V1)
GEO_COAD <- GEO_COAD[,sampletable$type[idx]=="Tumor"]


GEO_COAD_predictcancer=predit_svm_cancertype(dat2=t(GEO_COAD))
table(GEO_COAD_predictcancer)

#GEO_COAD_predictcancer
#TCGA-PRAD TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-BRCA TCGA-LIHC 
#0        94         0         2         0         0 

GEO_COAD_x=form_xpredict(dat2=t(GEO_COAD)) #0 probes not available
GEO_COAD_predictprob=predict(cvfit, s=lambda.best,newx = GEO_COAD_x, type = "response")
GEO_COAD_predictcancer=predit_glmnet_cancertype(allprobs=GEO_COAD_predictprob)
table(GEO_COAD_predictcancer$class)


#TCGA-COAD TCGA-LUSC 
#95         1






#GEO-PRAD
GEO_PRAD=as.data.frame(fread("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/PRAD/GSE76269/GSE76938_matrix_processed.txt"))
rownames(GEO_PRAD)=GEO_PRAD$TargetID
GEO_PRAD=GEO_PRAD[,-1]
idx=seq(1,ncol(GEO_PRAD),2)
GEO_PRAD=GEO_PRAD[,idx]
colnames(GEO_PRAD)=gsub(".AVG_Beta","",colnames(GEO_PRAD),fixed=T)
colnames(GEO_PRAD)=gsub(" ",".",colnames(GEO_PRAD),fixed=T)
sampletable=read.table("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/PRAD/GSE76269/sampletable.txt")
sampletable$type="Tumor"
idx=grepl("N",sampletable$V2)
sampletable$type[idx]="Normal"
table(sampletable$type)
# Normal  Tumor 
# 63     73
idx=match(colnames(GEO_PRAD),sampletable$V2)
GEO_PRAD <- GEO_PRAD[,sampletable$type[idx]=="Tumor"]



GEO_PRAD_predictcancer=predit_svm_cancertype(dat2=t(GEO_PRAD))
table(GEO_PRAD_predictcancer)

#TCGA-PRAD TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-BRCA TCGA-LIHC 
#73         0         0         0         0         0 


GEO_PRAD_x=form_xpredict(dat2=t(GEO_PRAD)) #0 probes not available
GEO_PRAD_predictprob=predict(cvfit, s=lambda.best,newx = GEO_PRAD_x, type = "response")
GEO_PRAD_predictcancer=predit_glmnet_cancertype(allprobs=GEO_PRAD_predictprob)
table(GEO_PRAD_predictcancer$class)



#GEO_Liver
GEO_LIVE=as.data.frame(fread("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/LIHC/GSE54503/GSE54503_series_matrix.txt",skip=59))
rownames(GEO_LIVE)=GEO_LIVE$ID_REF
GEO_LIVE=GEO_LIVE[,-1]
sampletable=read.table("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/LIHC/GSE54503/sampletable.txt",sep="\t")
table(sampletable$V3)
# adjacent non-tumor tissue                 HCC tumor 
# 66                        66 

idx=match(colnames(GEO_LIVE),sampletable$V1)
GEO_LIVE <- GEO_LIVE[,sampletable$V3[idx]=="HCC tumor"]


GEO_LIVE_predictcancer=predit_svm_cancertype(dat2=t(GEO_LIVE))
table(GEO_LIVE_predictcancer)




GEO_LIVE_x=form_xpredict(dat2=t(GEO_LIVE))
GEO_LIVE_predictprob=predict(cvfit, s=lambda.best,newx = GEO_LIVE_x, type = "response")
GEO_LIVE_predictcancer=predit_glmnet_cancertype(allprobs=GEO_LIVE_predictprob)
table(GEO_LIVE_predictcancer$class)

#TCGA-LIHC 
#66



#GEO_LIVE_predictcancer
#Normal Tumor
#adjacent non-tumor tissue     54    12
#HCC tumor                      1    65










