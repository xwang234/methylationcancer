#!/usr/bin/env Rscript
#classify tumor vs normals using TCGA data

Cancertypes=c("TCGA-PRAD","TCGA-COAD","TCGA-LUAD","TCGA-LUSC","TCGA-BRCA")

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

alluniqprobes=NULL
for (i in 1:length(sigprobes))
{
  tmp=uniqprobes[[i]]
  alluniqprobes=unique(c(alluniqprobes,tmp))
}
length(alluniqprobes) #[1] 16714

#get all methylation data and alltype factor
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
dim(allTCGAdat) #337005   2662 #all probes
#only include sig probes
allTCGAdat=allTCGAdat[rownames(allTCGAdat) %in% comprobes,]
dim(allTCGAdat) #32413  2662
alltype=c(rep("Tumor",ncol(alltumor)),rep("Normal",ncol(allnormal)))
TCGAtraindat=t(allTCGAdat)
y=alltype

require(doMC)
registerDoMC(cores=12)
set.seed(1000)
Sys.time()
cvfit=cv.glmnet(TCGAtraindat, y, family="binomial", standardize=F,parallel = T,nfolds=10, trace.it=1)
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
#sig probes
# y        Normal Tumor
# Normal    247    11
# Tumor      12  2392
#1se
# y        Normal Tumor
# Normal    238    20
# Tumor      20  2384

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

#Hutch PRAD data
library(data.table)
Hutch_PRAD=fread("/fh/fast/stanford_j/Illumina/Data_Cleaned/Clean_Methylation/beta_main_HP.csv")
Hutch_PRAD=as.data.frame(Hutch_PRAD)
rownames(Hutch_PRAD)=Hutch_PRAD[,1]
Hutch_PRAD=Hutch_PRAD[,2:ncol(Hutch_PRAD)]
Hutch_PRAD_x=form_xpredict() 
Hutch_PRAD_predictcancer=predit_cancertype()
table(Hutch_PRAD_predictcancer)
# Normal  Tumor 
# 25    500 
#1se
# Normal  Tumor 
# 33    492 
#use sig probes
# Normal  Tumor 
# 2    523
#1se
# Normal  Tumor 
# 13    512 

Hutch_PRAD_predictcancer_knn=knnclust()
table(Hutch_PRAD_predictcancer_knn$class)
# Normal  Tumor 
# 39    486
#1se
# Normal  Tumor 
# 42    483 
#
# Normal  Tumor 
# 31    494 
#1se
# Normal  Tumor 
# 29    496 

#Hutch normals
Hutch_PRAD_NORM=as.data.frame(fread("/fh/fast/stanford_j/Illumina/Data_Cleaned/Clean_Methylation/beta_tumnorm.csv"))
rownames(Hutch_PRAD_NORM)=Hutch_PRAD_NORM$V1
if ("V1" %in% colnames(Hutch_PRAD_NORM)) Hutch_PRAD_NORM=Hutch_PRAD_NORM[,-1]
mdat=read.csv("/fh/fast/stanford_j/Illumina/Data_Cleaned/Clean_Methylation/sample_info_clean.csv",header=T,stringsAsFactors = F)
Hutch_PRAD_NORM=Hutch_PRAD_NORM[,colnames(Hutch_PRAD_NORM) %in% mdat$sampleid[grepl("NORMAL",mdat$comments)]]
Hutch_PRAD_NORM_x=form_xpredict(dat2=t(Hutch_PRAD_NORM))

Hutch_PRAD_NORM_predictcancer=predit_cancertype(dat2=t(Hutch_PRAD_NORM))
table(Hutch_PRAD_NORM_predictcancer)
# Normal  Tumor 
# 4     16 
#1se
# Normal  Tumor 
# 7     13
#
# Tumor 
# 20 
#1se
# Tumor 
# 20 

Hutch_PRAD_NORM_predictcancer_knn=knnclust(dat2=t(Hutch_PRAD_NORM))
table(Hutch_PRAD_NORM_predictcancer_knn$class)
# Normal  Tumor 
# 6     14 
# #1se
# Normal  Tumor 
# 9     11 
#
# Normal  Tumor 
# 5     15 
#1se
# Normal  Tumor 
# 5     15

#Try Hutch COAD
load("/fh/fast/dai_j/CancerGenomics/Colorectal_Cancer/Data/GSE48684.RData")
CancerID <- row.names(clinicaltable)[clinicaltable$status=="cancer"]
AdenoID <- row.names(clinicaltable)[clinicaltable$status=="adenoma"]
NormalID <- row.names(clinicaltable)[clinicaltable$status=="normal-H"]
Grady_COAD_NORM <- MEall[,names(MEall)%in% NormalID]
Grady_COAD <- MEall[,names(MEall)%in% CancerID]
Grady_COAD_predictcancer=predit_cancertype(dat2=t(Grady_COAD))
table(Grady_COAD_predictcancer)
# Normal  Tumor 
# 7     57 
#1se
# Normal  Tumor 
# 8     56 
#
# Normal  Tumor 
# 6     58 
#1se
# Normal  Tumor 
# 5     59 

Grady_COAD_predictcancer_knn=knnclust(dat2=t(Grady_COAD))
table(Grady_COAD_predictcancer_knn$class)
# Normal  Tumor 
# 7     57 
#1se
# Normal  Tumor 
# 6     58
#
# Normal  Tumor 
# 7     57 
#1se
# Normal  Tumor 
# 6     58

Grady_COAD_NORM_x=form_xpredict(dat2=t(Grady_COAD_NORM))
Grady_COAD_NORM_predictcancer=predit_cancertype(dat2=t(Grady_COAD_NORM))
table(Grady_COAD_NORM_predictcancer)
# Normal 
# 17
#1se
# Normal 
# 17
#
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
#
# Normal  Tumor 
# 17      0 
#1se
# Normal  Tumor 
# 17      0 

#Try GEO BRCA
GEO_BRCA=as.data.frame(fread("data/BRCA/GSE75067/GSE75067_series_matrix.txt",skip=58))
rownames(GEO_BRCA)=GEO_BRCA$ID_REF
if ("ID_REF" %in% colnames(GEO_BRCA)) GEO_BRCA=GEO_BRCA[,-1]
GEO_BRCA_predictcancer=predit_cancertype(dat2=t(GEO_BRCA))
table(GEO_BRCA_predictcancer,useNA="ifany")
# Normal  Tumor 
# 1    187
#1se
# Normal  Tumor 
# 1    187
#
# Normal  Tumor 
# 2    186 
#1se
# Normal  Tumor 
# 2    186 

GEO_BRCA_predictcancer_knn=knnclust(dat2=t(GEO_BRCA))
table(GEO_BRCA_predictcancer_knn$class)
# Normal  Tumor 
# 2    186
#1se
# Normal  Tumor 
# 1    187
#
# Normal  Tumor 
# 5    183 
#1se
# Normal  Tumor 
# 5    183 

#Try lung cancers
GEO_LUNG=as.data.frame(fread("data/LUNG/GSE56044/GSE56044_series_matrix.txt",skip=69))
rownames(GEO_LUNG)=GEO_LUNG$ID_REF
if ("ID_REF" %in% colnames(GEO_LUNG)) GEO_LUNG=GEO_LUNG[,-1]
sampletable=readxl::read_xls("data/LUNG/GSE56044/GSE56044_GEO_Annotations.xls",sheet = 2)
sampleidtable=readxl::read_xls("data/LUNG/GSE56044/GSE56044_sampleID.xls")
idx=match(sampletable$SampleID,sampleidtable$ID)
sampletable$ID=sampleidtable$SapleID[idx]
GEO_LUAD=GEO_LUNG[,colnames(GEO_LUNG) %in% sampletable$ID[sampletable$Histology=="AC"]]
GEO_LUSC=GEO_LUNG[,colnames(GEO_LUNG) %in% sampletable$ID[sampletable$Histology=="SqCC"]]
GEO_LUNG_NORM=GEO_LUNG[,colnames(GEO_LUNG) %in% sampletable$ID[sampletable$Histology=="NA"]]
GEO_LUAD_predictcancer=predit_cancertype(dat2=t(GEO_LUAD))
table(GEO_LUAD_predictcancer)
# Normal  Tumor 
# 1     82 
#1se
# Normal  Tumor 
# 1     82 
#
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
#
# Normal  Tumor 
# 1     82 

GEO_LUSC_predictcancer=predit_cancertype(dat2=t(GEO_LUSC))
table(GEO_LUSC_predictcancer)
# Tumor 
# 23 
#1se
# Tumor 
# 23 
#
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
#
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
#
# Normal  Tumor 
# 11      1 
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
#
# Normal  Tumor 
# 12      0
#1se
# Normal  Tumor 
# 12      0 

#plot some figures-----
allcolors=c("red","blue","darkorchid1","limegreen","goldenrod1","black","brown4","darkseagreen","darkgoldenrod","cadetblue4")
#allcolors=c("red","blue","darkorchid1","limegreen","black","darkcyan","goldenrod1","darkseagreen","darkseagreen1")
plot.pca=function(dat=allTCGAdat,probes=comprobes,types=alltype,yinc=1.2,pc1=1,pc2=2,main="",prefix=NULL,opt="beta")
{
  types=factor(types,levels = unique(types))
  dat=dat[rownames(dat) %in% probes,]
  row.means=rowMeans(dat,na.rm=T)
  for (i in 1:ncol(dat))
  {
    idx=which(is.na(dat[,i]))
    if (length(idx)>0)
    {
      dat[idx,i]=row.means[idx]
    }
  }
  library(matrixStats)
  dat0=dat
  if (opt!="beta")
  {
    dat[dat <= 0]=min(dat0[dat0!=0])/10
    dat[dat >= 1]=max(dat0[dat0!=1])+(1-max(dat0[dat0!=1]))*0.5
    #M values
    dat=log(dat/(1-dat),base=2)
  }
  #remove constant rows, var=0
  tmp=rowSds(as.matrix(dat))
  dat=dat[tmp>0,]
  pcadat=prcomp(t(dat),scale = T)
  expl_var <- pcadat$sdev^2/sum(pcadat$sdev^2)
  if (!is.null(prefix))
  {
    png(paste0(resfolder,prefix,"PCs_variances.png"),width = 480, height = 480,type = "cairo")
    par(mar=c(6,6,2,1))
    barplot(expl_var[1:50], ylab="EXPLAINED VARIANCE",
            names.arg=paste0("pcadat",seq(1:50)), col="darkgreen",las=2)
    dev.off()
  }
  pcadat=pcadat$x
  xmin=min(pcadat[,pc1])
  xmax=max(pcadat[,pc1])
  ymin=min(pcadat[,pc2])
  ymax=max(pcadat[,pc2])*yinc
  pch=1:length(unique(types))
  pch=pch[types]
  plot(pcadat[,pc1],pcadat[,pc2],col=alpha(allcolors[types],0.4),pch=pch,cex.axis=1.5,cex.lab=1.5,xlim=c(xmin-0.05*(xmax-xmin),xmax+0.05*(xmax-xmin)),
       ylim=c(ymin-0.05*(ymax-ymin),ymax+0.2*(ymax-ymin)),xlab=paste0("PC",pc1," (",round(expl_var[pc1]*100),"%)"),ylab=paste0("PC",pc2," (",round(expl_var[pc2]*100),"%)"),main=main,bty='l')
  #text(x=pcadat[,1],y=pcadat[,2],rownames(pcadat),col=colors[clinicaltable$`Platium response`],cex=1)
  legend("topleft",legend=c(unique(as.character(types))),col=allcolors[1:length(unique(types))],pch=unique(as.numeric(types)),cex=1,bty = "n",ncol=2)
}

#https://towardsdatascience.com/how-to-tune-hyperparameters-of-tsne-7c0596a18868
plot.tsne=function(dat=allTCGAdat,probes=comprobes,types=alltype,yinc=1.2,pc1=1,pc2=2,main="",opt="beta")
{
  types=factor(types,levels = unique(types))
  dat=dat[rownames(dat) %in% probes,]
  row.means=rowMeans(dat,na.rm=T)
  for (i in 1:ncol(dat))
  {
    idx=which(is.na(dat[,i]))
    if (length(idx)>0)
    {
      dat[idx,i]=row.means[idx]
    }
  }
  dat0=dat
  if (opt!="beta")
  {
    dat[dat <= 0]=min(dat0[dat0!=0])/10
    dat[dat >= 1]=max(dat0[dat0!=1])+(1-max(dat0[dat0!=1]))*0.5
    #M values
    dat=log(dat/(1-dat),base=2)
  }
  library(matrixStats)
  #remove constant rows, var=0
  tmp=rowSds(as.matrix(dat))
  dat=dat[tmp>0,]
  pcadat=prcomp(t(dat),scale = T)
  # expl_var <- pcadat$sdev^2/sum(pcadat$sdev^2)
  # barplot(expl_var[1:50], ylab="EXPLAINED VARIANCE",
  #         names.arg=paste0("pcadat",seq(1:50)), col="darkgreen")
  pcadat=pcadat$x
  library(Rtsne)
  set.seed(1500)
  n=as.integer(sqrt(nrow(pcadat)))
  tsnedat = Rtsne(pcadat, pca = FALSE,perplexity = n,theta=0,initial_dims = 50)
  tsnedat=tsnedat$Y
  rownames(tsnedat)=rownames(pcadat)
  xmin=min(tsnedat[,pc1])
  xmax=max(tsnedat[,pc1])
  ymin=min(tsnedat[,pc2])
  ymax=max(tsnedat[,pc2])*yinc
  par(mar=c(5,5,2,1))
  pch=1:length(unique(types))
  pch=pch[types]
  plot(tsnedat[,pc1],tsnedat[,pc2],col=alpha(allcolors[types],0.4),pch=pch,cex.axis=1.5,cex.lab=1.5,xlim=c(xmin-0.05*(xmax-xmin),xmax+0.05*(xmax-xmin)),
       ylim=c(ymin-0.05*(ymax-ymin),ymax+0.2*(ymax-ymin)),xlab=paste0("TSNE",pc1),ylab=paste0("TSNE",pc2),main=main,bty='l')
  #legend("topleft",legend=c(Cancertypes),col=allcolors[1:length(Cancertypes)],pch=1,cex=1,bty = "n",ncol=3)
  legend("topleft",legend=c(unique(as.character(types))),col=allcolors[1:length(unique(types))],pch=unique(as.numeric(types)),cex=1,bty = "n",ncol=3)
}

# png(paste0(resfolder,"AllTCGATumorvsNormal_glmprobes_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
# par(mar=c(6,6,2,1))
# plot.pca(main="Glmnet selected probes",probes=selectprobes1)
# dev.off()
alltype1=c(as.character(tumortype),rep("TCGA-NORM",ncol(allnormal))) #6 classes
png(paste0(resfolder,"AllTCGATumorvsNormal_glmprobes1_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(main="Glmnet selected probes",probes=selectprobes1,types=alltype1)
dev.off()
# png(paste0(resfolder,"AllTCGATumorvsNormal_glmprobes_PCA_PC1_PC3.png"),width = 480, height = 480,type = "cairo")
# par(mar=c(6,6,2,1))
# plot.pca(main="Glmnet selected probes",probes=selectprobes1,pc2=3)
# dev.off()
png(paste0(resfolder,"AllTCGATumorvsNormal_glmprobes1_PCA_PC1_PC3.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(main="Glmnet selected probes",probes=selectprobes1,types=alltype1,pc2=3)
dev.off()

# png(paste0(resfolder,"AllTCGATumorvsNormal_glmprobes_TSNE.png"),width = 480, height = 480,type = "cairo")
# par(mar=c(6,6,2,1))
# plot.tsne(main="Glmnet selected probes",probes=selectprobes1)
dev.off()
png(paste0(resfolder,"AllTCGATumorvsNormal_glmprobes1_TSNE.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.tsne(main="Glmnet selected probes",probes=selectprobes1,types=alltype1)
dev.off()

#include Hutch normals
png(paste0(resfolder,"AllTCGATumorvsNormal_HutchPRADNORM_glmprobes_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=t(rbind(TCGAtraindat,Hutch_PRAD_NORM_x)),types=c(as.character(alltype1),rep("Hutch_PRADN",ncol(Hutch_PRAD_NORM))),main="Glmnet selected probes",probes=selectprobes1)
dev.off()
png(paste0(resfolder,"AllTCGATumorvsNormal_HutchPRADNORM_glmprobes_TSNE.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.tsne(dat=t(rbind(TCGAtraindat,Hutch_PRAD_NORM_x)),types=c(as.character(alltype1),rep("Hutch_PRADN",ncol(Hutch_PRAD_NORM))),main="Glmnet selected probes",probes=selectprobes1)
dev.off()

#include Grady COAD normals
png(paste0(resfolder,"AllTCGATumorvsNormal_GradyCOADNORM_glmprobes_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=t(rbind(TCGAtraindat,Grady_COAD_NORM_x)),types=c(as.character(alltype1),rep("Grady_COADN",ncol(Grady_COAD_NORM))),main="Glmnet selected probes",probes=selectprobes1)
dev.off()
png(paste0(resfolder,"AllTCGATumorvsNormal_GradyCOADNORM_glmprobes_TSNE.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.tsne(dat=t(rbind(TCGAtraindat,Grady_COAD_NORM_x)),types=c(as.character(alltype1),rep("Grady_COADN",ncol(Grady_COAD_NORM))),main="Glmnet selected probes",probes=selectprobes1)
dev.off()
#study what happened to Hutch normal
normaltype1=paste0(normaltype,"-NORM")
sum(colnames(TCGAtraindat) %in% sigprobes$`TCGA-PRAD`)

plot.pca(dat=t(rbind(TCGAtraindat[alltype1 %in% c("TCGA-PRAD","TCGA-NORM"),],Hutch_PRAD_x,Hutch_PRAD_NORM_x)),types=c(rep("TCGA-PRAD",sum(alltype1=="TCGA-PRAD")),normaltype1,rep("HUTCH-PRAD",ncol(Hutch_PRAD)),rep("Hutch-NORM",ncol(Hutch_PRAD_NORM))),main="PRAD sig probes",probes=sigprobes$`TCGA-PRAD`)
plot.pca(dat=t(rbind(TCGAtraindat[alltype1 %in% c("TCGA-PRAD","TCGA-NORM"),],Hutch_PRAD_x,Hutch_PRAD_NORM_x)),types=c(rep("TCGA-PRAD",sum(alltype1=="TCGA-PRAD")),normaltype1,rep("HUTCH-PRAD",ncol(Hutch_PRAD)),rep("Hutch-NORM",ncol(Hutch_PRAD_NORM))),pc2=3,main="PRAD sig probes",probes=sigprobes$`TCGA-PRAD`)
png(paste0(resfolder,"TCGAHUTCHPRAD_NORMAL_PRADprobes_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=t(rbind(TCGAtraindat[alltype1 %in% c("TCGA-PRAD","TCGA-NORM"),],Hutch_PRAD_x,Hutch_PRAD_NORM_x)),types=c(rep("TCGA-PRAD",sum(alltype1=="TCGA-PRAD")),rep("TCGA-NORM",ncol(allnormal)),rep("HUTCH-PRAD",ncol(Hutch_PRAD)),rep("Hutch-NORM",ncol(Hutch_PRAD_NORM))),main="PRAD sig probes",probes=sigprobes$`TCGA-PRAD`)
dev.off()
plot.pca(dat=t(rbind(TCGAtraindat[alltype1 %in% c("TCGA-PRAD","TCGA-NORM"),],Hutch_PRAD_x,Hutch_PRAD_NORM_x)),types=c(rep("TCGA-PRAD",sum(alltype1=="TCGA-PRAD")),rep("TCGA-NORM",ncol(allnormal)),rep("HUTCH-PRAD",ncol(Hutch_PRAD)),rep("Hutch-NORM",pc2=3,ncol(Hutch_PRAD_NORM))),main="PRAD sig probes",probes=sigprobes$`TCGA-PRAD`)
png(paste0(resfolder,"TCGAHUTCHPRAD_NORMAL_PRADprobes_TSNE1.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.tsne(dat=t(rbind(TCGAtraindat[alltype1 %in% c("TCGA-PRAD","TCGA-NORM"),],Hutch_PRAD_x,Hutch_PRAD_NORM_x)),types=c(rep("TCGA-PRAD",sum(alltype1=="TCGA-PRAD")),normaltype1,rep("HUTCH-PRAD",ncol(Hutch_PRAD)),rep("Hutch-NORM",ncol(Hutch_PRAD_NORM))),main="PRAD sig probes",probes=sigprobes$`TCGA-PRAD`)
dev.off()
png(paste0(resfolder,"TCGAHUTCHPRAD_NORMAL_PRADprobes_TSNE2.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.tsne(dat=t(rbind(TCGAtraindat[alltype1 %in% c("TCGA-PRAD","TCGA-NORM"),],Hutch_PRAD_x,Hutch_PRAD_NORM_x)),types=c(rep("TCGA-PRAD",sum(alltype1=="TCGA-PRAD")),rep("TCGA-NORM",ncol(allnormal)),rep("HUTCH-PRAD",ncol(Hutch_PRAD)),rep("Hutch-NORM",ncol(Hutch_PRAD_NORM))),main="PRAD sig probes",probes=sigprobes$`TCGA-PRAD`)
dev.off()

png(paste0(resfolder,"TCGAPRAD_NORMAL_PRADprobes_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=t(rbind(TCGAtraindat[alltype1 %in% c("TCGA-PRAD","TCGA-NORM"),])),types=c(rep("TCGA-PRAD",sum(alltype1=="TCGA-PRAD")),rep("TCGA-NORM",ncol(allnormal))),main="PRAD sig probes",probes=sigprobes$`TCGA-PRAD`)
dev.off()
png(paste0(resfolder,"TCGAPRAD_NORMAL_PRADprobes_PCA_PC1_PC3.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=t(rbind(TCGAtraindat[alltype1 %in% c("TCGA-PRAD","TCGA-NORM"),])),types=c(rep("TCGA-PRAD",sum(alltype1=="TCGA-PRAD")),rep("TCGA-NORM",ncol(allnormal))),pc2=3,main="PRAD sig probes",probes=sigprobes$`TCGA-PRAD`)
dev.off()
#Only include TCGA PRAD data
alltype2=c(tumortype,normaltype1)
dat=t(rbind(TCGAtraindat[alltype2 %in% c("TCGA-PRAD","TCGA-PRAD-NORM"),]))
types=c(rep("TCGA-PRAD",sum(alltype1=="TCGA-PRAD")),rep("TCGA-PRAD-NORM",sum(alltype2=="TCGA-PRAD-NORM")))
png(paste0(resfolder,"TCGAPRAD_PRADNORMAL_PRADprobes_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=dat,types=types,main="PRAD sig probes",probes=sigprobes$`TCGA-PRAD`)
dev.off()
png(paste0(resfolder,"TCGAPRAD_PRADNORMAL_PRADprobes_TSNE.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.tsne(dat=dat,types=types,main="PRAD sig probes",probes=sigprobes$`TCGA-PRAD`)
dev.off()
plot.pca(dat=dat,types=types,main="PRAD sig probes",pc2=3,probes=sigprobes$`TCGA-PRAD`)

png(paste0(resfolder,"TCGAPRAD_NORMAL_PRADprobes_TSNE1.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.tsne(dat=t(rbind(TCGAtraindat[alltype1 %in% c("TCGA-PRAD","TCGA-NORM"),])),types=c(rep("TCGA-PRAD",sum(alltype1=="TCGA-PRAD")),rep("TCGA-NORM",ncol(allnormal))),main="PRAD sig probes",probes=sigprobes$`TCGA-PRAD`)
dev.off()
png(paste0(resfolder,"TCGAPRAD_NORMAL_PRADprobes_TSNE2.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.tsne(dat=t(rbind(TCGAtraindat[alltype1 %in% c("TCGA-PRAD","TCGA-NORM"),])),types=c(rep("TCGA-PRAD",sum(alltype1=="TCGA-PRAD")),normaltype1),main="PRAD sig probes",probes=sigprobes$`TCGA-PRAD`)
dev.off()
png(paste0(resfolder,"HUTCHPRAD_NORMAL_PRADprobes_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=t(rbind(Hutch_PRAD_x,Hutch_PRAD_NORM_x)),types=c(rep("HUTCH-PRAD",ncol(Hutch_PRAD)),rep("Hutch-NORM",ncol(Hutch_PRAD_NORM))),main="PRAD sig probes",probes=sigprobes$`TCGA-PRAD`)
dev.off()
png(paste0(resfolder,"HUTCHPRAD_NORMAL_PRADprobes_PCA_PC1_PC3.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=t(rbind(Hutch_PRAD_x,Hutch_PRAD_NORM_x)),types=c(rep("HUTCH-PRAD",ncol(Hutch_PRAD)),rep("Hutch-NORM",ncol(Hutch_PRAD_NORM))),pc2=3,main="PRAD sig probes",probes=sigprobes$`TCGA-PRAD`)
dev.off()
png(paste0(resfolder,"HUTCHPRAD_NORMAL_PRADprobes_TSNE.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.tsne(dat=t(rbind(Hutch_PRAD_x,Hutch_PRAD_NORM_x)),types=c(rep("HUTCH-PRAD",ncol(Hutch_PRAD)),rep("Hutch-NORM",ncol(Hutch_PRAD_NORM))),main="PRAD sig probes",probes=sigprobes$`TCGA-PRAD`)
dev.off()

#to compare sig probes in TCGA PRAD and Hutch PRAD
dat1=TCGAtraindat[alltype2 %in% "TCGA-PRAD-NORM",]
dat2=Hutch_PRAD_NORM_x
dat1=TCGAtraindat[alltype2 %in% "TCGA-PRAD",]
dat2=Hutch_PRAD_x
dat1=dat1[,colnames(dat1) %in% sigprobes$`TCGA-PRAD`]
dat2=dat2[,colnames(dat2) %in% sigprobes$`TCGA-PRAD`]
all(colnames(dat1)==colnames(dat2))
plot(colMeans(dat1),colMeans(dat2),xlab="Mean beta value in TCGA normal",ylab="Mean beta value in Hutch normal")
abline(0,1,col="red")
boxplot(t(rbind(dat1[1:50,],dat2[1:20,])),col=c(rep("red",50),rep("blue",20)),outline=F)

dat1=Hutch_PRAD_NORM
dat2=cbind(alltumor,allnormal)[,alltype2=="TCGA-PRAD-NORM"]
theprobes=intersect(rownames(dat1),rownames(dat2))
idx1=match(theprobes,rownames(dat1))
idx2=match(theprobes,rownames(dat2))
dat1=dat1[idx1,]
dat2=dat2[idx2,]
plot(rowMeans(dat2),rowMeans(dat1),xlab="Mean beta value in TCGA normal",ylab="Mean beta value in Hutch normal")
abline(0,1,col="red")
quantile(rowSds(as.matrix(dat1),na.rm=T))
# 0%         25%         50%         75%        100% 
# 0.002032003 0.026679653 0.053095487 0.092228811 0.330425150
quantile(rowSds(as.matrix(dat2),na.rm=T))
# 0%          25%          50%          75%         100% 
# 0.0007715403 0.0109566595 0.0320697411 0.0768264670 0.3494309195 
plot(rowSds(as.matrix(dat2),na.rm=T),rowSds(as.matrix(dat1),na.rm=T),xlab="SD beta value in TCGA normal",ylab="SD beta value in Hutch normal")
abline(0,1,col="red")
fit=glm(rowSds(as.matrix(dat1),na.rm=T)~rowSds(as.matrix(dat2),na.rm=T))
abline(fit,col="pink")
table(rowSds(as.matrix(dat2))>rowSds(as.matrix(dat1)))
# FALSE   TRUE 
# 269956  77365
plot(rowMeans(dat2),rowSds(as.matrix(dat2)),xlab="Mean beta value in TCGA normal",ylab="SD beta value in TCGA normal")
plot(rowMeans(dat1),rowSds(as.matrix(dat1)),xlab="Mean beta value in Hutch normal",ylab="SD beta value in Hutch normal")

tmp1=rowMeans(dat1)
quantile(tmp1)
theprobes1=rownames(dat1)[tmp1<0.1]
tmp2=rowMeans(dat2,na.rm=T)
quantile(tmp2)
theprobes2=rownames(dat2)[tmp2<0.1]
idx1=match(theprobes1,rownames(dat1))
idx2=match(theprobes1,rownames(dat2))
dat11=dat1[idx1,]
dat22=dat2[idx2,]
plot(rowMeans(dat22),rowMeans(dat11),xlab="Mean beta value in TCGA normal",ylab="Mean beta value in Hutch normal")
abline(0,1,col="red")
table(rowMeans(dat22)>rowMeans(dat11))
idx1=match(theprobes2,rownames(dat1))
idx2=match(theprobes2,rownames(dat2))
dat11=dat1[idx1,]
dat22=dat2[idx2,]
plot(rowMeans(dat22),rowMeans(dat11),xlab="Mean beta value in TCGA normal",ylab="Mean beta value in Hutch normal")
abline(0,1,col="red")

dat1=Hutch_PRAD
dat2=cbind(alltumor,allnormal)[,alltype2=="TCGA-PRAD"]
theprobes=intersect(rownames(dat1),rownames(dat2))
idx1=match(theprobes,rownames(dat1))
idx2=match(theprobes,rownames(dat2))
dat1=dat1[idx1,]
dat2=dat2[idx2,]
plot(rowMeans(dat2),rowMeans(dat1),xlab="Mean beta value in TCGA tumor",ylab="Mean beta value in Hutch tumor")
abline(0,1,col="red")
quantile(rowSds(as.matrix(dat1),na.rm=T))
# 0%         25%         50%         75%        100% 
# 0.003062375 0.025797832 0.049129575 0.087185181 0.395302498 
quantile(rowSds(as.matrix(dat2),na.rm=T))
# 0%         25%         50%         75%        100% 
# 0.001366114 0.018338349 0.056644882 0.116691981 0.357985499 
plot(rowSds(as.matrix(dat2),na.rm=T),rowSds(as.matrix(dat1),na.rm=T),xlab="SD beta value in TCGA tumor",ylab="SD beta value in Hutch tumor")
abline(0,1,col="red")
fit=glm(rowSds(as.matrix(dat1),na.rm=T)~rowSds(as.matrix(dat2),na.rm=T))
abline(fit,col="pink")
table(rowSds(as.matrix(dat2))>rowSds(as.matrix(dat1)))
# FALSE   TRUE 
# 269956  77365
plot(rowMeans(dat2),rowSds(as.matrix(dat2)),xlab="Mean beta value in TCGA tumor",ylab="SD beta value in TCGA tumor")
plot(rowMeans(dat1),rowSds(as.matrix(dat1)),xlab="Mean beta value in Hutch tumor",ylab="SD beta value in Hutch tumor")


plot(rowMeans(dat2),rowMeans(dat1),xlab="Mean beta value in TCGA tumor",ylab="Mean beta value in Hutch tumor")
abline(0,1,col="red")
tmp1=rowMeans(dat1)
quantile(tmp1)
tmp2=rowMeans(dat2,na.rm=T)
quantile(tmp2)


#do heatmap on all TCGA data
library(ComplexHeatmap)
library(GetoptLong)
library(circlize)
normaltype1=paste0(normaltype,"-NORM")
uniq_cancertype=c(Cancertypes,unique(normaltype1))
listcolor_cancertype=rep("NA",length(uniq_cancertype))
names(listcolor_cancertype)=uniq_cancertype
for (i in 1:length(uniq_cancertype))
{
  listcolor_cancertype[i]=allcolors[i]
}
listcolor_cancertype=list(Type=listcolor_cancertype)

drawheatmap=function(mat=t(TCGAtraindat[,colnames(TCGAtraindat) %in% selectprobes1]), cluster=alltype,
                     showrowname=F)
{
  
  # mat=t(scale(t(mat)))
  # mat[mat < -2]=-2
  # mat[mat > 2]=2
  ha = HeatmapAnnotation(Type = cluster, 
                         col = listcolor_cancertype,
                         show_annotation_name = T,
                         show_legend=T,
                         annotation_name_gp=gpar(fontface = "bold"),
                         # annotation_name_offset = unit(2, "cm"),
                         # annotation_name_rot = c(0, 0),
                         annotation_name_side = "right")
  
  
  ht_list = Heatmap(mat, name = "Beta",
                    col = colorRamp2(c(0,0.5,1), c("white", "pink", "red")), 
                    #column_dend_height = unit(4, "cm"),
                    cluster_rows = T,
                    cluster_columns=F,
                    column_dend_reorder=F,
                    #clustering_distance_columns="canberra",#kendall",
                    row_names_gp=gpar(fontsize = 8, fontface = "bold"),
                    column_dend_height=unit(2,"cm"),
                    row_dend_width=unit(2,"cm"),
                    top_annotation = c(ha),
                    show_row_dend=T,
                    show_heatmap_legend=T,
                    heatmap_height=unit(12,"cm"),
                    show_column_names = FALSE, show_row_names=showrowname) 
  ht_list = draw(ht_list, heatmap_legend_side = "right")
  return(ht_list)
}
png(paste0(resfolder,"AllTCGA_Heatmap_glmprobes_nostand.png"),width = 800, height = 480,type = "cairo")
drawheatmap()
dev.off()