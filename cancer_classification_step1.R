#!/usr/bin/env Rscript
#classify tumor vs normals using TCGA data

Cancertypes=c("TCGA-PRAD","TCGA-COAD","TCGA-LUAD","TCGA-LUSC","TCGA-BRCA")

alltumor=tumortype=allnormal=normaltype=NULL
for (i in 1:length(Cancertypes))
{
  Cancer.type=Cancertypes[i]
  rm(tcga_normal_methy,tcga_tumor_methy,tcga_tumor_methy_all)
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

allTCGAdat=cbind(alltumor,allnormal)
dim(allTCGAdat) #348813   2662
tmp=rowMeans(allTCGAdat)
idx=which(is.na(tmp))
length(idx) #10973, remove probes with missing values (some may correlated to specific cancer)
allTCGAdat=allTCGAdat[-idx,]
#remove non-cpgs
idx=grepl("^cg",rownames(allTCGAdat))
allTCGAdat=allTCGAdat[idx,]
dim(allTCGAdat) #337005   2662
alltype=c(rep("Tumor",ncol(alltumor)),rep("Normal",ncol(allnormal)))
TCGAtraindat=t(allTCGAdat)
y=alltype
require(doMC)
registerDoMC(cores=12)
set.seed(1000)
Sys.time()
cvfit=cv.glmnet(TCGAtraindat, y, family="binomial", standardize=F,parallel = T,nfolds=10, trace.it=1)
save(cvfit,file="result/cvfit_tcgacancervsnormal.RData")

lambda.best=cvfit$lambda.min #102 probes
lambda.best=cvfit$lambda.1se #64 probes
glmcoeff1=coef(cvfit,s=lambda.best)
idx=glmcoeff1@i
#the selected probes
selectprobes1=glmcoeff1@Dimnames[[1]][idx+1]
selectprobes1=selectprobes1[selectprobes1!="(Intercept)"]
length(selectprobes1)

predictprob=predict(cvfit, s=lambda.best,newx = TCGAtraindat, type = "response")
predictcancer=rep("Normal",nrow(predictprob))
predictcancer[predictprob[,1]>0.5]="Tumor"
table(y,predictcancer)
#         predictcancer
# y        Normal Tumor
# Normal    254     4
# Tumor       5  2399
#1se
# y        Normal Tumor
# Normal    244    14
# Tumor      12  2392
quantile(predictprob[y=="Normal",1],c(0,0.5,0.75,0.9,0.95,1))
# 0%          50%          75%          90%          95%         100% 
#   0.0003157367 0.0348728103 0.1037140137 0.2178391152 0.3282303760 0.8269132166 
quantile(predictprob[y=="Tumor",1],c(0,0.01,0.05,0.1,0.5,0.75,1))
# 0%        1%        5%       10%       50%       75%      100% 
# 0.1243325 0.7955683 0.9683022 0.9903173 0.9997050 0.9999461 1.0000000 
#Tumor predict to normal
idx=which(y=="Tumor" & predictcancer=="Normal")
predictprob[idx,1]
#  0.4704941    0.3482796    0.1243325    0.3756682    0.3456957
idx=which(y=="Normal" & predictcancer=="Tumor")
#Normal predict to tumor
predictprob[idx,1]
#0.6805503    0.7076562    0.8269132    0.5569130
predictcancer=predict(cvfit, s=lambda.best,newx = TCGAtraindat, type = "class")[,1]
table(y,predictcancer)
#          predictcancer
# y        Normal Tumor
# Normal    254     4
# Tumor       5  2399

#for testing data
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
predit_cancertype=function(dat1=TCGAtraindat,dat2=t(Hutch_PRAD),selectprobes=selectprobes1)
{
  x=form_xpredict(dat1=dat1,dat2=dat2,selectprobes=selectprobes)
  predictcancer=predict(cvfit, s=lambda.best,newx = x, type = "class")[,1] 
  return(predictcancer)
}
Hutch_PRAD_predictcancer=predit_cancertype()
table(Hutch_PRAD_predictcancer)
# Normal  Tumor 
# 25    500 
#1se
# Normal  Tumor 
# 33    492 

#Use KNN for prediction
knnclust=function(dat1=TCGAtraindat,dat2=t(Hutch_PRAD),selectprobes=selectprobes1,dat1_class=as.character(alltype))
{
  #only use available probes both in training and testing
  avaiprobes=intersect(colnames(dat1),colnames(dat2))
  print(paste0(sum(!selectprobes %in% avaiprobes)," probes not available!"))
  avaiprobes=avaiprobes[avaiprobes %in% selectprobes]
  idx1=match(avaiprobes,colnames(dat1))
  idx2=match(avaiprobes,colnames(dat2))
  dat1=dat1[,idx1]
  dat2=dat2[,idx2]
  tmp=sum(is.na(dat2))
  if (tmp>0) 
  {
    warning(paste0(tmp," NA exist"))
    for (i in 1:ncol(dat2))
    {
      idx=which(is.na(dat2[,i]))
      if (length(idx)>0)
      {
        dat2[idx,i]=mean(dat2[,i],na.rm=T)
      }
    }
  }
  
  library(class)
  knnClust <- knn(train =dat1, test = dat2 , k = 5, cl = dat1_class,prob=T)
  # quantile(attributes(knnClust)$prob)
  prob=attributes(knnClust)$prob
  names(knnClust)=rownames(dat2)
  res=data.frame(class=knnClust,prob=prob)
  rownames(res)=rownames(dat2)
  return(res)
}

Hutch_PRAD_predictcancer_knn=knnclust()
table(Hutch_PRAD_predictcancer_knn$class)
# Normal  Tumor 
# 39    486
#1se
# Normal  Tumor 
# 42    483 

Hutch_PRAD_NORM_predictcancer=predit_cancertype(dat2=t(Hutch_PRAD_NORM))
table(Hutch_PRAD_NORM_predictcancer)
# Normal  Tumor 
# 4     16 
#1se
# Normal  Tumor 
# 7     13
Hutch_PRAD_NORM_predictcancer_knn=knnclust(dat2=t(Hutch_PRAD_NORM))
table(Hutch_PRAD_NORM_predictcancer_knn$class)
# Normal  Tumor 
# 6     14 
# #1se
# Normal  Tumor 
# 9     11 

Grady_COAD_predictcancer=predit_cancertype(dat2=t(Grady_COAD))
table(Grady_COAD_predictcancer)
# Normal  Tumor 
# 7     57 
#1se
# Normal  Tumor 
# 8     56 
Grady_COAD_predictcancer_knn=knnclust(dat2=t(Grady_COAD))
table(Grady_COAD_predictcancer_knn$class)
# Normal  Tumor 
# 7     57 
#1se
# Normal  Tumor 
# 6     58

Grady_COAD_NORM_predictcancer=predit_cancertype(dat2=t(Grady_COAD_NORM))
table(Grady_COAD_NORM_predictcancer)
# Normal 
# 17
#1se
# Normal 
# 17
Grady_COAD_NORM_predictcancer_knn=knnclust(dat2=t(Grady_COAD_NORM))
table(Grady_COAD_NORM_predictcancer_knn$class)
# Normal  Tumor 
# 17      0 
#1se
# Normal  Tumor 
# 17      0 

GEO_BRCA_predictcancer=predit_cancertype(dat2=t(GEO_BRCA))
table(GEO_BRCA_predictcancer,useNA="ifany")
# Normal  Tumor 
# 1    187
#1se
# Normal  Tumor 
# 1    187
GEO_BRCA_predictcancer_knn=knnclust(dat2=t(GEO_BRCA))
table(GEO_BRCA_predictcancer_knn$class)
# Normal  Tumor 
# 2    186
#1se
# Normal  Tumor 
# 1    187

GEO_LUAD_predictcancer=predit_cancertype(dat2=t(GEO_LUAD))
table(GEO_LUAD_predictcancer)
# Normal  Tumor 
# 1     82 
#1se
# Normal  Tumor 
# 1     82 
GEO_LUAD_predictcancer_knn=knnclust(dat2=t(GEO_LUAD))
table(GEO_LUAD_predictcancer_knn$class)
# Normal  Tumor 
# 2     81 
#1se
# Normal  Tumor 
# 2     81

GEO_LUSC_predictcancer=predit_cancertype(dat2=t(GEO_LUSC))
table(GEO_LUSC_predictcancer)
# Tumor 
# 23 
#1se
# Tumor 
# 23 
GEO_LUSC_predictcancer_knn=knnclust(dat2=t(GEO_LUSC))
table(GEO_LUSC_predictcancer_knn$class)
# Normal  Tumor 
# 0     23 
#1se
# Normal  Tumor 
# 0     23 

GEO_LUNG_NORM_predictcancer=predit_cancertype(dat2=t(GEO_LUNG_NORM))
table(GEO_LUNG_NORM_predictcancer)
# Normal 
# 12 
#1se
# Normal 
# 12 
GEO_LUNG_NORM_predictcancer_knn=knnclust(dat2=t(GEO_LUNG_NORM))
table(GEO_LUNG_NORM_predictcancer_knn$class)
# Normal  Tumor 
# 12      0 
#1se
# Normal  Tumor 
# 12      0 