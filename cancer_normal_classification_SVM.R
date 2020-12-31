
#classify tumor vs normals using TCGA data

rm(list=ls())

load("/fh/fast/dai_j/CancerGenomics/prostate_methylation/result/TCGA-PRAD_output.RData")
bdiff <- output[,6]-output[,4]
pout <- output[output[,23]<0.05/nrow(output) & bdiff>0 & output[,4]<0.1 & output[,7]>0.05,]

load("/fh/fast/dai_j/CancerGenomics/prostate_methylation/result/TCGA-COAD_output.RData")
bdiff <- output[,6]-output[,4]
cout <- output[output[,23]<0.05/nrow(output) & bdiff>0.05 & output[,4]<0.1 & output[,7]>0.05,]

load("/fh/fast/dai_j/CancerGenomics/prostate_methylation/result/TCGA-LUSC_output.RData")
bdiff <- output[,6]-output[,4]
lsout <- output[output[,23]<0.05/nrow(output) & bdiff>0.05 & output[,4]<0.1 & output[,7]>0.05,]

load("/fh/fast/dai_j/CancerGenomics/prostate_methylation/result/TCGA-LUAD_output.RData")
bdiff <- output[,6]-output[,4]
laout <- output[output[,23]<0.05/nrow(output) & bdiff>0.05 & output[,4]<0.1 & output[,7]>0.05,]

load("/fh/fast/dai_j/CancerGenomics/prostate_methylation/result/TCGA-BRCA_output.RData")
bdiff <- output[,6]-output[,4]
bout <- output[output[,23]<0.05/nrow(output) & bdiff>0.05 & output[,4]<0.1 & output[,7]>0.05,]

load("/fh/fast/dai_j/CancerGenomics/prostate_methylation/result/TCGA-LIHC_output.RData")
bdiff <- output[,6]-output[,4]
bout <- output[output[,23]<0.05/nrow(output) & bdiff>0.05 & output[,4]<0.1 & output[,7]>0.05,]

shareprobe <- row.names(pout)[row.names(pout) %in% row.names(cout) &  row.names(pout) %in% row.names(bout) &  row.names(pout) %in% row.names(laout) &  row.names(pout) %in% row.names(lsout)]

rm(output)





rm(list=ls())

resfolder="/fh/fast/dai_j/CancerGenomics/prostate_methylation/result/"
Cancertypes=c("TCGA-PRAD","TCGA-COAD","TCGA-LUAD","TCGA-LUSC","TCGA-BRCA","TCGA-LIHC")

sigprobes=NULL

for (i in 1:length(Cancertypes))
{
  Cancer.type=Cancertypes[i]
  Res.file=paste0(resfolder,Cancer.type,"_output.RData")
  load(Res.file)
  bdiff <- output[,6]-output[,4]
  tmp=rownames(output)[output[,23]<0.05/nrow(output) & bdiff>0.05 & output[,4]<0.1 & output[,7]>0.05]
  tmp=tmp[grepl("^cg",tmp)] #remove non-CpG targeting probes ch.xx.xx
  if (i==1) shareprobes<-tmp else shareprobes<- intersect(shareprobes, tmp)
  sigprobes=c(sigprobes,list(probes=tmp))
  names(sigprobes)[length(sigprobes)]=Cancer.type
}

length(shareprobes)
#[1] 823

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
length(comprobes) #36948
for (i in 1:length(sigprobes))
{
  print(paste0(Cancertypes[i],":",length(sigprobes[[i]]),",",length(uniqprobes[[i]])))
}

#total significant probes,unique probes:
# [1] "TCGA-PRAD:12371,2855"
# [1] "TCGA-COAD:16823,5887"
# [1] "TCGA-LUAD:11165,569"
# [1] "TCGA-LUSC:16225,3954"
# [1] "TCGA-BRCA:18025,3449"
# [1] "TCGA-LIHC:12060,2106"

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
  print(dim(tcga_tumor_methy_all))
  print(dim(tcga_normal_methy))
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
  print(dim(alltumor))
  print(dim(allnormal))
}
allTCGAdat=cbind(alltumor,allnormal)
dim(allTCGAdat) #[1] 339604   3089 
tmp=rowMeans(allTCGAdat)
idx=which(is.na(tmp))
length(idx) #10973, remove probes with missing values (some may correlated to specific cancer)
allTCGAdat=allTCGAdat[-idx,]
#remove non-cpgs
idx=grepl("^cg",rownames(allTCGAdat))
allTCGAdat=allTCGAdat[idx,]
dim(allTCGAdat) #329239   3089  #all probes
#only include sig probes


allTCGAdat=allTCGAdat[rownames(allTCGAdat) %in% comprobes,]

#allTCGAdat=allTCGAdat[rownames(allTCGAdat) %in% shareprobe,]
dim(allTCGAdat) #32413  2662
alltype=c(rep("Tumor",ncol(alltumor)),rep("Normal",ncol(allnormal)))
TCGAtraindat=t(allTCGAdat)
y=alltype


require(doMC)
registerDoMC(cores=12)
set.seed(5446)
Sys.time()


cvfit=cv.glmnet(TCGAtraindat, y, family="binomial", standardize=F,parallel = T,nfolds=10, trace.it=1)
#cvfit=cv.glmnet(TCGAtraindat, y, family="binomial", standardize=F,parallel = T,nfolds=5, alpha=0.5,trace.it=1)
#save(cvfit,file="result/cvfit_tcgacancervsnormal.RData")
Sys.time()
lambda.best=cvfit$lambda.min #102 probes (all probes),100 probes (sig probes)
#lambda.best=cvfit$lambda.1se #46 probes, 47 probes
glmcoeff1=coef(cvfit,s=lambda.best)
idx=glmcoeff1@i
#the selected probes
selectprobes1=glmcoeff1@Dimnames[[1]][idx+1]
selectprobes1=selectprobes1[selectprobes1!="(Intercept)"]
length(selectprobes1)
#[1] 103
predictprob=predict(cvfit, s=lambda.best,newx = TCGAtraindat, type = "response")
predictcancer=rep("Normal",nrow(predictprob))
predictcancer[predictprob[,1]>0.5]="Tumor"
table(y,predictcancer)

#min
##predictcancer
##y        Normal Tumor
##Normal    295    13
##Tumor      15  2766

#1se
#predictcancer
#y        Normal Tumor
#Normal    291    17
#Tumor      24  2757

sum(selectprobes1%in%shareprobes)

TCGAtraindat1 <- TCGAtraindat[,colnames(TCGAtraindat) %in% selectprobes1]

y<-factor(y)
library(e1071)
tune.out=tune(svm ,TCGAtraindat1,y,scale=F,  kernel ="radial", ranges =list(cost=c(0.1 ,1 ,10,20),gamma=c(0.01,0.1,0.5,1) ))
summary(tune.out)
#bestmod =tune.out$best.model
#summary(bestmod)

svmfit <- svm(TCGAtraindat1,y,scale=F,kernel ="radial", cost=10,gamma = 0.1)
pred <- fitted(svmfit)
table(y,pred)

#pred
#y        Normal Tumor
#Normal    307     1
#Tumor       6  2775







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

predit_glmnet_cancertype=function(dat1=TCGAtraindat,dat2=t(Hutch_PRAD),selectprobes=selectprobes1)
{
  x=form_xpredict(dat1=dat1,dat2=dat2,selectprobes=selectprobes1)
  #x <- x[,colnames(x) %in% selectprobes1]
  predictprob=predict(cvfit, s=lambda.best,newx = as.matrix(x), type = "response")
  predictcancer=rep("Normal",nrow(predictprob))
  predictcancer[predictprob[,1]>0.5]="Tumor"
  return(predictcancer)
}



#Hutch PRAD data
library(data.table)
Hutch_PRAD=fread("/fh/fast/stanford_j/Illumina/Data_Cleaned/Clean_Methylation/beta_main_HP.csv")
Hutch_PRAD=as.data.frame(Hutch_PRAD)
rownames(Hutch_PRAD)=Hutch_PRAD[,1]
Hutch_PRAD=Hutch_PRAD[,2:ncol(Hutch_PRAD)]


Hutch_PRAD_predictcancer=predit_svm_cancertype()
table(Hutch_PRAD_predictcancer)
#Hutch_PRAD_predictcancer
#Normal  Tumor 
#1    524


Hutch_PRAD_predictcancer=predit_glmnet_cancertype()
table(Hutch_PRAD_predictcancer)
#Normal  Tumor 
#2    523


#Hutch_PRAD_NORM=as.data.frame(fread("/fh/fast/stanford_j/Illumina/Data_Cleaned/Clean_Methylation/beta_tumnorm.csv"))
#rownames(Hutch_PRAD_NORM)=Hutch_PRAD_NORM$V1
#if ("V1" %in% colnames(Hutch_PRAD_NORM)) Hutch_PRAD_NORM=Hutch_PRAD_NORM[,-1]
#mdat=read.csv("/fh/fast/stanford_j/Illumina/Data_Cleaned/Clean_Methylation/sample_info_clean.csv",header=T,stringsAsFactors = F)
#Hutch_PRAD_NORM=Hutch_PRAD_NORM[,colnames(Hutch_PRAD_NORM) %in% mdat$sampleid[grepl("NORMAL",mdat$comments)]]

#Hutch_PRAD_NORM_predictcancer=predit_svm_cancertype(dat2=t(Hutch_PRAD_NORM))
#table(Hutch_PRAD_NORM_predictcancer)



load("/fh/fast/dai_j/CancerGenomics/Colorectal_Cancer/Data/GSE48684.RData")
CancerID <- row.names(clinicaltable)[clinicaltable$status=="cancer"]
AdenoID <- row.names(clinicaltable)[clinicaltable$status=="adenoma"]
NormalID <- row.names(clinicaltable)[clinicaltable$status=="normal-H"]
Grady_COAD_NORM <- MEall[,names(MEall)%in% NormalID]
Grady_COAD <- MEall[,names(MEall)%in% CancerID]


Grady_COAD_predictcancer=predit_glmnet_cancertype(dat2=t(Grady_COAD))
table(Grady_COAD_predictcancer)
#Normal  Tumor 
#6     58 
Grady_COAD_predictcancer=predit_svm_cancertype(dat2=t(Grady_COAD))
table(Grady_COAD_predictcancer)

#Normal  Tumor 
#6     58

Grady_COAD_NORM_predictcancer=predit_svm_cancertype(dat2=t(Grady_COAD_NORM))
table(Grady_COAD_NORM_predictcancer)
#Normal  Tumor 
#17      0
Grady_COAD_NORM_predictcancer=predit_glmnet_cancertype(dat2=t(Grady_COAD_NORM))
table(Grady_COAD_NORM_predictcancer)
#Normal  Tumor 
#17      0



GEO_BRCA=as.data.frame(fread("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/BRCA/GSE75067/GSE75067_series_matrix.txt",skip=58))
rownames(GEO_BRCA)=GEO_BRCA$ID_REF
if ("ID_REF" %in% colnames(GEO_BRCA)) GEO_BRCA=GEO_BRCA[,-1]
GEO_BRCA_predictcancer=predit_svm_cancertype(dat2=t(GEO_BRCA))
table(GEO_BRCA_predictcancer,useNA="ifany")
#Normal  Tumor 
#3    185

GEO_BRCA_predictcancer=predit_glmnet_cancertype(dat2=t(GEO_BRCA))
table(GEO_BRCA_predictcancer,useNA="ifany")
#Normal  Tumor 
#2    186



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
GEO_LUNG_NORM=GEO_LUNG[,colnames(GEO_LUNG) %in% sampletable$ID[sampletable$Histology=="NA"]]


GEO_LUAD_predictcancer=predit_svm_cancertype(dat2=t(GEO_LUAD))
table(GEO_LUAD_predictcancer)
#Normal  Tumor 
#1     82


GEO_LUAD_predictcancer=predit_glmnet_cancertype(dat2=t(GEO_LUAD))
table(GEO_LUAD_predictcancer)
#Normal  Tumor 
#1     82

GEO_LUSC_predictcancer=predit_svm_cancertype(dat2=t(GEO_LUSC))
table(GEO_LUSC_predictcancer)
#Normal  Tumor 
#0     23


GEO_LUSC_predictcancer=predit_glmnet_cancertype(dat2=t(GEO_LUSC))
table(GEO_LUSC_predictcancer)
#Tumor 
#23 


GEO_LUNG_NORM_predictcancer=predit_svm_cancertype(dat2=t(GEO_LUNG_NORM))
table(GEO_LUNG_NORM_predictcancer)
#Normal  Tumor 
#12      0


GEO_LUNG_NORM_predictcancer=predit_glmnet_cancertype(dat2=t(GEO_LUNG_NORM))
table(GEO_LUNG_NORM_predictcancer)
#Normal 
#11   1 
















###additional data from GEO
#GEO-BRCA-NORM
GEO_BRCA_NORM=as.data.frame(fread("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/BRCA/GSE101961/GSE101961_series_matrix.txt",skip=67)) #black and white
rownames(GEO_BRCA_NORM)=GEO_BRCA_NORM$ID_REF
GEO_BRCA_NORM=GEO_BRCA_NORM[,-1]


GEO_BRCA_NORM_predictcancer=predit_svm_cancertype(dat2=t(GEO_BRCA_NORM))
table(GEO_BRCA_NORM_predictcancer)
#Normal  Tumor 
#119      2

GEO_BRCA_NORM_predictcancer=predit_glmnet_cancertype(dat2=t(GEO_BRCA_NORM))
table(GEO_BRCA_NORM_predictcancer)

# TCGA-NORM 
# 120                1



#GEO-LUNG-NORM1
GEO_LUNG_NORM1=as.data.frame(fread("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/LUNG/GSE52401/GSE52401_series_matrix.txt",skip=65))
rownames(GEO_LUNG_NORM1)=GEO_LUNG_NORM1$ID_REF
GEO_LUNG_NORM1=GEO_LUNG_NORM1[,-1]


GEO_LUNG_NORM_predictcancer=predit_svm_cancertype(dat2=t(GEO_LUNG_NORM1))
table(GEO_LUNG_NORM_predictcancer)
#Normal  Tumor 
#217     27 

GEO_LUNG_NORM_predictcancer=predit_glmnet_cancertype(dat2=t(GEO_LUNG_NORM1))
table(GEO_LUNG_NORM_predictcancer)

#Normal  Tumor 
#216     28

 

#Grady-COLON-NORM
load("/fh/fast/dai_j/CancerGenomics/Colorectal_Cancer/Data/GSE132804.RData")   
LowID <- row.names(clinicaltable)[clinicaltable$risk=="Low"& clinicaltable$platform=="HM450"] #48
MedID <- row.names(clinicaltable)[clinicaltable$risk=="Medium"& clinicaltable$platform=="HM450"] #31
HighID <- row.names(clinicaltable)[clinicaltable$risk=="High"& clinicaltable$platform=="HM450"] #49
# Lmethy <- MEall[,names(MEall)%in% LowID]
# Mmethy <-  MEall[,names(MEall)%in% MedID]
# Hmethy <- MEall[,names(MEall)%in% HighID]
Grady_COLON_NORM=MEall[,names(MEall) %in% c(LowID)]


GEO_COLON_NORM_predictcancer=predit_svm_cancertype(dat2=t(Grady_COLON_NORM))
table(GEO_COLON_NORM_predictcancer)
#Normal  Tumor 
#45      3 


GEO_COLON_NORM_predictcancer=predit_glmnet_cancertype(dat2=t(Grady_COLON_NORM))
table(GEO_COLON_NORM_predictcancer)
#Normal  Tumor 
#47      1 




#GEO_COAD
GEO_COAD=as.data.frame(fread("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/COAD/GSE77718/GSE77718_series_matrix.txt",skip=65))
rownames(GEO_COAD)=GEO_COAD$ID_REF
GEO_COAD=GEO_COAD[,-1]
sampletable=read.table("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/COAD/GSE77718/sampletable.txt")
sampletable$type="Tumor"
idx=grepl("N",sampletable$V7)
sampletable$type[idx]="Normal"

GEO_COAD_x=form_xpredict(dat2=t(GEO_COAD)) #57 probes not available


GEO_COAD_predictcancer=predit_svm_cancertype(dat2=t(GEO_COAD))
idx=match(colnames(GEO_COAD),sampletable$V1)
table(sampletable$type[idx],GEO_COAD_predictcancer)


#Normal Tumor
#Normal     92     4
#Tumor       3    93



GEO_COAD_predictcancer=predit_glmnet_cancertype(dat2=t(GEO_COAD))
idx=match(colnames(GEO_COAD),sampletable$V1)
table(sampletable$type[idx],GEO_COAD_predictcancer)

#Normal Tumor
#Normal     95     1
#Tumor      13    83


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
GEO_PRAD_x=form_xpredict(dat2=t(GEO_PRAD)) #0 probes not available


GEO_PRAD_predictcancer=predit_svm_cancertype(dat2=t(GEO_PRAD))
idx=match(colnames(GEO_PRAD),sampletable$V2)
table(sampletable$type[idx],GEO_PRAD_predictcancer)

#GEO_PRAD_predictcancer
#Normal Tumor
#Normal     54     9
#Tumor       4    69




GEO_PRAD_predictcancer=predit_glmnet_cancertype(dat2=t(GEO_PRAD))
idx=match(colnames(GEO_PRAD),sampletable$V2)
table(sampletable$type[idx],GEO_PRAD_predictcancer)

#GEO_PRAD_predictcancer
#Normal Tumor
#Normal     55    8
#Tumor       3    70





#GEO_Liver
GEO_LIVE=as.data.frame(fread("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/LIHC/GSE54503/GSE54503_series_matrix.txt",skip=59))
rownames(GEO_LIVE)=GEO_LIVE$ID_REF
GEO_LIVE=GEO_LIVE[,-1]
sampletable=read.table("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/LIHC/GSE54503/sampletable.txt",sep="\t")
table(sampletable$V3)
# adjacent non-tumor tissue                 HCC tumor 
# 66                        66 
GEO_LIVE_x=form_xpredict(dat2=t(GEO_LIVE)) 


GEO_LIVE_predictcancer=predit_svm_cancertype(dat2=t(GEO_LIVE))
idx=match(colnames(GEO_LIVE),sampletable$V1)
table(sampletable$V3[idx],GEO_LIVE_predictcancer)

#GEO_LIVE_predictcancer
#Normal Tumor
#adjacent non-tumor tissue     54    12
#HCC tumor                      1    65




GEO_LIVE_predictcancer=predit_glmnet_cancertype(dat2=t(GEO_LIVE))
idx=match(colnames(GEO_LIVE),sampletable$V1)
table(sampletable$V3[idx],GEO_LIVE_predictcancer)


#GEO_LIVE_predictcancer
#Normal Tumor
#adjacent non-tumor tissue     49    17
#HCC tumor                      1    65
