#!/usr/bin/env Rscript
#cancer classification include normals

#Get significant probes in all cancer types----------------------
resfolder="/fh/fast/dai_j/CancerGenomics/prostate_methylation/result/"
#find significant probes in each cancer type
#Cancertypes=c("TCGA-PRAD","TCGA-COAD","TCGA-LUAD","TCGA-LUSC","TCGA-BRCA","Krause-EAC")
#remove EAC
Cancertypes=c("TCGA-PRAD","TCGA-COAD","TCGA-LUAD","TCGA-LUSC","TCGA-BRCA")
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
# [1] "Krause-EAC:21228,5152"
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
all(comprobes %in% rownames(alltumor)) #F
sum(alluniqprobes %in% rownames(alltumor)) #14719
sum(comprobes %in% rownames(alltumor)) #33733
tumortype0=tumortype
alltype=c(as.character(tumortype),rep("TCGA-NORM",ncol(allnormal)))
tumortype=factor(tumortype,levels = Cancertypes)
alltype=factor(alltype,levels=c("TCGA-NORM",Cancertypes))
allTCGAdat=cbind(alltumor,allnormal)
tmp=rowMeans(allTCGAdat)
idx=which(is.na(tmp))
length(idx) #10973, remove probes with missing values (some may correlated to specific cancer)
allTCGAdat=allTCGAdat[-idx,]
sum(comprobes %in% rownames(allTCGAdat)) #32413

allcolors=c("red","blue","darkorchid1","limegreen","goldenrod1","black","brown4","darkseagreen","darkseagreen1")
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
#png(paste0(resfolder,"Allcancers_ALLsig_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
png(paste0(resfolder,"AllTCGAsig_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(main="All significant probes")
dev.off()
#png(paste0(resfolder,"Allcancers_ALLsig_PCA_PC1_PC3.png"),width = 480, height = 480,type = "cairo")
png(paste0(resfolder,"AllTCGAsig_PCA_PC1_PC3.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(main="All significant probes",pc1=1,pc2=3)
dev.off()

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
#png(paste0(resfolder,"Allcancers_ALLsig_TSNE.png"),width = 480, height = 480,type = "cairo")
png(paste0(resfolder,"ALLsig_TSNE.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.tsne(main="All significant probes")
dev.off()

#glmnet------
library(glmnet)
x=allTCGAdat[rownames(allTCGAdat) %in% comprobes,]
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
dim(x) #32413  2662
x=remove_NA_withmean(dat=x)
dim(x) #32413  2662
y=alltype
#Training data
TCGAtraindat=t(x)
#remember to set #export OMP_NUM_THREADS=1 when use parallel
require(doMC)
registerDoMC(cores=12)
set.seed(1000)
Sys.time()

cvfit=cv.glmnet(TCGAtraindat, y, family="multinomial", type.multinomial = "grouped", standardize=F,parallel = T,nfolds=10, trace.it=1)
save(cvfit,file="result/cvfit_all.RData") #tumors and normal
#Or try without stand
#cvfit=cv.glmnet(TCGAtraindat, y, family="multinomial", type.multinomial = "grouped", standardize=F,parallel = T,nfolds=10, trace.it=1)
plot(cvfit)
Sys.time()
# use customized lambda,add small lambda, slow
# lambda=10^seq(0,-5,length=100)
# lambda=10^seq(-0.3,-10,length=100) #cvfit10 Convergence for 94th lambda value not reached after maxit=100000 iterations; solutions for larger lambdas returned
# # lambda=10^seq(-0.3,-20,length=100) #cvfit20 Convergence for 47th lambda value not reached after maxit=100000 iterations; solutions for larger lambdas returned
# set.seed(1000)
# cvfit=cv.glmnet(t(x), y, family="multinomial", type.multinomial = "grouped", parallel = T,standardize=F,
#                 lambda=lambda,trace.it=1)
# #png(paste0(resfolder,"Allcancers_glmprobes_cvfit10_stand.png.png"),width = 480, height = 480,type = "cairo")
# save(cvfit,file="result/cvfit10_all.RData")
# plot(cvfit)
# Sys.time()
#dev.off()
lambda.best=cvfit$lambda.min #196 probes, 184 overlap
lambda.best=cvfit$lambda.1se #189 probes, 297
glmcoeff1=coef(cvfit,s=lambda.best)[1]
idx=glmcoeff1$`TCGA-NORM`@i
#the selected probes
selectprobes1=glmcoeff1$`TCGA-NORM`@Dimnames[[1]][idx+1]
selectprobes1=selectprobes1[selectprobes1!="(Intercept)"]
glmcoeff5=coef(cvfit,s=lambda.best)[5]
idx=glmcoeff5$`TCGA-LUSC`@i
selectprobes5=glmcoeff5$`TCGA-LUSC`@Dimnames[[1]][idx+1]
selectprobes5=selectprobes5[selectprobes5!="(Intercept)"]
all(selectprobes1==selectprobes5) #T used to verify
length(selectprobes1)
#png(paste0(resfolder,"Allcancers_glmprobes_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
png(paste0(resfolder,"AllTCGA_glmprobes_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(main="Glmnet selected probes",probes=selectprobes1)
dev.off()
png(paste0(resfolder,"AllTCGA_glmprobes_PCA_PC1_PC3.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(main="Glmnet selected probes",probes=selectprobes1,pc1=1,pc2=3)
dev.off()
png(paste0(resfolder,"AllTCGA_glmprobes_TSNE.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.tsne(main="Glmnet selected probes",probes=selectprobes1)
dev.off()

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
# resps
# y           TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD
# TCGA-PRAD         0         0         0         0         9       489
# TCGA-COAD         0       296         0         0         0         0
# TCGA-LUAD         0         0       449         5         4         0
# TCGA-LUSC         0         0        13       355         2         0
# TCGA-BRCA       777         0         0         0         5         0
# TCGA-NORM         7         0         0         0       240        11
#1se
# resps
# y           TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD
# TCGA-PRAD         0         0         0         0        10       488
# TCGA-COAD         0       296         0         0         0         0
# TCGA-LUAD         0         0       448         5         5         0
# TCGA-LUSC         0         0        16       352         2         0
# TCGA-BRCA       776         0         0         0         6         0
# TCGA-NORM         7         0         0         0       239        12
#1se best
# resps
# y           TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD
# TCGA-NORM         3         0         0         0       249         6
# TCGA-PRAD         0         0         0         0         5       493
# TCGA-COAD         0       296         0         0         0         0
# TCGA-LUAD         0         0       454         3         1         0
# TCGA-LUSC         0         0         2       366         2         0
# TCGA-BRCA       778         0         0         0         4         0

#Testing data----------
#Hutch PRAD data
library(data.table)
Hutch_PRAD=fread("/fh/fast/stanford_j/Illumina/Data_Cleaned/Clean_Methylation/beta_main_HP.csv")
Hutch_PRAD=as.data.frame(Hutch_PRAD)
rownames(Hutch_PRAD)=Hutch_PRAD[,1]
Hutch_PRAD=Hutch_PRAD[,2:ncol(Hutch_PRAD)]
sum(!colnames(TCGAtraindat) %in% rownames(Hutch_PRAD)) #29
sum(!selectprobes1 %in% rownames(Hutch_PRAD)) #1
#form x matrix used for glmnet prediction
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
Hutch_PRAD_x=form_xpredict() #0 probes not available
Hutch_PRAD_predictprob=predict(cvfit, s=lambda.best,newx = Hutch_PRAD_x, type = "response")
#linear model to predict cancer type
predit_cancertype=function(allprobs=Hutch_PRAD_predictprob)
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
Hutch_PRAD_predictcancer=predit_cancertype()
table(Hutch_PRAD_predictcancer$class)
# TCGA-NORM TCGA-PRAD 
# 27       498 
#1se
# TCGA-NORM TCGA-PRAD 
# 30       495
#1se best
# TCGA-NORM TCGA-PRAD 
# 9       516 

png(paste0(resfolder,"Allcancers_HutchPRAD_glmprobes_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=t(rbind(TCGAtraindat,Hutch_PRAD_x)),types=c(as.character(tumortype),rep("Hutch_PRAD",nrow(Hutch_PRAD_x))),main="Glmnet selected probes",probes=selectprobes1)
dev.off()
png(paste0(resfolder,"Allcancers_HutchPRAD_glmprobes_PCA_PC1_PC3.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=t(rbind(TCGAtraindat,Hutch_PRAD_x)),types=c(as.character(tumortype),rep("Hutch_PRAD",nrow(Hutch_PRAD_x))),
         pc2=3,main="Glmnet selected probes",probes=selectprobes1)
dev.off()
png(paste0(resfolder,"Allcancers_HutchPRAD_glmprobes_PCA_PC2_PC3.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=t(rbind(TCGAtraindat,Hutch_PRAD_x)),types=c(as.character(tumortype),rep("Hutch_PRAD",nrow(Hutch_PRAD_x))),
         pc1=2,pc2=3,main="Glmnet selected probes",probes=selectprobes1)
dev.off()
png(paste0(resfolder,"Allcancers_HutchPRAD_glmprobes_TSNE.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.tsne(dat=t(rbind(TCGAtraindat,Hutch_PRAD_x)),types=c(as.character(tumortype),rep("Hutch_PRAD",nrow(Hutch_PRAD_x))),main="Glmnet selected probes",probes=selectprobes1)
dev.off()

#only use available probes both in training and testing,row:samples,col:probes,used for knnclust
form_xpredict1=function(dat1=TCGAtraindat,dat2=t(Hutch_PRAD),selectprobes=selectprobes1)
{
  
  avaiprobes=intersect(colnames(dat1),colnames(dat2))
  print(paste0(sum(!selectprobes %in% avaiprobes)," probes not available!"))
  avaiprobes=avaiprobes[avaiprobes %in% selectprobes]
  idx1=match(avaiprobes,colnames(dat1))
  idx2=match(avaiprobes,colnames(dat2))
  res1=dat1[,idx1]
  res2=dat2[,idx2]
  tmp=sum(is.na(res2))
  if (tmp>0) 
  {
    warning(paste0(tmp," NA exist"))
    for (i in 1:ncol(res2))
    {
      idx=which(is.na(res2[,i]))
      if (length(idx)>0)
      {
        res2[idx,i]=mean(res2[,i],na.rm=T)
      }
    }
  }
  return(list(dat1=res1,dat2=res2))
}

#Use KNN for prediction
alldat=form_xpredict1()
knnclust=function(dat1=alldat$dat1,dat2=alldat$dat2,dat1_class=as.character(alltype))
{
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
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 0         0         0         0        32       493 
#1se
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 0         0         0         0        32       493 
#1se best
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 0         0         0         0        36       489
quantile(Hutch_PRAD_predictcancer_knn$prob)
# 0%  25%  50%  75% 100% 
# 0.6  1.0  1.0  1.0  1.0 

#Hutch normals
Hutch_PRAD_NORM=as.data.frame(fread("/fh/fast/stanford_j/Illumina/Data_Cleaned/Clean_Methylation/beta_tumnorm.csv"))
rownames(Hutch_PRAD_NORM)=Hutch_PRAD_NORM$V1
if ("V1" %in% colnames(Hutch_PRAD_NORM)) Hutch_PRAD_NORM=Hutch_PRAD_NORM[,-1]
mdat=read.csv("/fh/fast/stanford_j/Illumina/Data_Cleaned/Clean_Methylation/sample_info_clean.csv",header=T,stringsAsFactors = F)
Hutch_PRAD_NORM=Hutch_PRAD_NORM[,colnames(Hutch_PRAD_NORM) %in% mdat$sampleid[grepl("NORMAL",mdat$comments)]]

Hutch_PRAD_NORM_x=form_xpredict(dat2=t(Hutch_PRAD_NORM)) #0 probes not available
Hutch_PRAD_NORM_predictprob=predict(cvfit, s=lambda.best,newx = Hutch_PRAD_NORM_x, type = "response")
Hutch_PRAD_NORM_predictcancer=predit_cancertype(allprobs=Hutch_PRAD_NORM_predictprob)
table(Hutch_PRAD_NORM_predictcancer$class)
# TCGA-NORM TCGA-PRAD 
# 6        14
#1se
# TCGA-NORM TCGA-PRAD 
# 9        11 
#1se best
# TCGA-NORM TCGA-PRAD 
# 1        19 

png(paste0(resfolder,"All_HutchPRADNORM_glmprobes_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=t(rbind(TCGAtraindat,Hutch_PRAD_NORM_x)),types=c(as.character(alltype),rep("Hutch_PRADNORM",nrow(Hutch_PRAD_NORM_x))),main="Glmnet selected probes",probes=selectprobes1)
dev.off()

alldat=form_xpredict1(dat2=t(Hutch_PRAD_NORM))
Hutch_PRAD_NORM_predictcancer_knn=knnclust()
table(Hutch_PRAD_NORM_predictcancer_knn$class)
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 0         0         0         0         8        12 
#1se
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 0         0         1         0         7        12
#1se best
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 0         0         1         0         8        11 

#Try Hutch COAD
load("/fh/fast/dai_j/CancerGenomics/Colorectal_Cancer/Data/GSE48684.RData")
CancerID <- row.names(clinicaltable)[clinicaltable$status=="cancer"]
AdenoID <- row.names(clinicaltable)[clinicaltable$status=="adenoma"]
NormalID <- row.names(clinicaltable)[clinicaltable$status=="normal-H"]
Grady_COAD_NORM <- MEall[,names(MEall)%in% NormalID]
Grady_COAD <- MEall[,names(MEall)%in% CancerID]
#Grady_Amethy <- MEall[,names(MEall)%in% AdenoID]
Grady_COAD_x=form_xpredict(dat2=t(Grady_COAD)) #0 probes not available
Grady_COAD_predictprob=predict(cvfit, s=lambda.best,newx = Grady_COAD_x, type = "response")
Grady_COAD_predictcancer=predit_cancertype(allprobs=Grady_COAD_predictprob)
table(Grady_COAD_predictcancer$class)
# TCGA-COAD TCGA-LUAD TCGA-NORM 
# 54         2         8 
#1se
# TCGA-COAD TCGA-LUAD TCGA-NORM 
# 54         2         8
#1se best
# TCGA-COAD TCGA-LUAD TCGA-NORM 
# 54         2         8 

png(paste0(resfolder,"Allcancers_GradyCOAD_glmprobes_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=t(rbind(TCGAtraindat,Grady_COAD_x)),types=c(as.character(tumortype),rep("Grady_COAD",nrow(Grady_COAD_x))),main="Glmnet selected probes",probes=selectprobes1)
dev.off()
png(paste0(resfolder,"Allcancers_GradyCOAD_glmprobes_PCA_PC1_PC3.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=t(rbind(TCGAtraindat,Grady_COAD_x)),types=c(as.character(tumortype),rep("Grady_COAD",nrow(Grady_COAD_x))),
         pc2=3,main="Glmnet selected probes",probes=selectprobes1)
dev.off()
png(paste0(resfolder,"Allcancers_GradyCOAD_glmprobes_PCA_PC2_PC3.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=t(rbind(TCGAtraindat,Grady_COAD_x)),types=c(as.character(tumortype),rep("Grady_COAD",nrow(Grady_COAD_x))),
         pc1=2,pc2=3,main="Glmnet selected probes",probes=selectprobes1)
dev.off()
# Grady_COAD_predictcancer_pcaknn=pcaclust(dat2=Grady_COAD_x)
# table(Grady_COAD_predictcancer_pcaknn)
# # TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-PRAD 
# # 0        56         5         3         0 

alldat=form_xpredict1(dat2=t(Grady_COAD))
Grady_COAD_predictcancer_knn=knnclust()
table(Grady_COAD_predictcancer_knn$class)
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 0        56         1         0         7         0 
#1se
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 0        54         2         0         8         0
#1se best
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 0        55         1         0         8         0 

# Hutch COAD NORMAL
Grady_COAD_NORM_x=form_xpredict(dat2=t(Grady_COAD_NORM)) #0 probes not available
Grady_COAD_NORM_predictprob=predict(cvfit, s=lambda.best,newx = Grady_COAD_NORM_x, type = "response")
Grady_COAD_NORM_predictcancer=predit_cancertype(allprobs=Grady_COAD_NORM_predictprob)
table(Grady_COAD_NORM_predictcancer$class)
# TCGA-NORM 
# 17
#1se
# TCGA-NORM 
# 17
#1se best
# TCGA-NORM 
# 17
alldat=form_xpredict1(dat2=t(Grady_COAD_NORM))
Grady_COAD_NORM_predictcancer_knn=knnclust()
table(Grady_COAD_NORM_predictcancer_knn$class)
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 0         0         0         0        17         0 
#1se
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 0         0         0         0        17         0
#1se best
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 0         0         0         0        17         0
#Try GEO BRCA

GEO_BRCA=as.data.frame(fread("data/BRCA/GSE75067/GSE75067_series_matrix.txt",skip=58))
rownames(GEO_BRCA)=GEO_BRCA$ID_REF
if ("ID_REF" %in% colnames(GEO_BRCA)) GEO_BRCA=GEO_BRCA[,-1]

GEO_BRCA_x=form_xpredict(dat2=t(GEO_BRCA)) #0 probes not available
GEO_BRCA_predictprob=predict(cvfit, s=lambda.best,newx = GEO_BRCA_x, type = "response")
GEO_BRCA_predictcancer=predit_cancertype(allprobs=GEO_BRCA_predictprob)
table(GEO_BRCA_predictcancer$class,useNA="ifany")
# TCGA-BRCA TCGA-NORM 
# 187         1 
#1se
# TCGA-BRCA TCGA-NORM 
# 187         1 
#1se best
# TCGA-BRCA TCGA-NORM 
# 186         2 
png(paste0(resfolder,"Allcancers_GEOBRCA_glmprobes_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=t(rbind(TCGAtraindat,GEO_BRCA_x)),types=c(as.character(tumortype),rep("GEO_BRCA",nrow(GEO_BRCA_x))),main="Glmnet selected probes",probes=selectprobes1)
dev.off()
png(paste0(resfolder,"Allcancers_GEOBRCA_glmprobes_PCA_PC1_PC3.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=t(rbind(TCGAtraindat,GEO_BRCA_x)),types=c(as.character(tumortype),rep("GEO_BRCA",nrow(GEO_BRCA_x))),
         pc2=3,main="Glmnet selected probes",probes=selectprobes1)
dev.off()

alldat=form_xpredict1(dat2=t(GEO_BRCA))
GEO_BRCA_predictcancer_knn=knnclust()
table(GEO_BRCA_predictcancer_knn$class)
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 180         0         0         0         8         0 
#1se
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 183         0         0         0         5         0 
#1se best
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 179         0         0         0         9         0 

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

GEO_LUAD_x=form_xpredict(dat2=t(GEO_LUAD)) #0 probes not available
GEO_LUAD_predictprob=predict(cvfit, s=lambda.best,newx = GEO_LUAD_x, type = "response")
GEO_LUAD_predictcancer=predit_cancertype(allprobs=GEO_LUAD_predictprob)
table(GEO_LUAD_predictcancer$class)
# TCGA-LUAD TCGA-LUSC TCGA-NORM 
# 80         2         1 
#1se
# TCGA-LUAD TCGA-LUSC TCGA-NORM 
# 80         2         1 
#1se best
# TCGA-LUAD TCGA-LUSC TCGA-NORM 
# 80         2         1 

png(paste0(resfolder,"Allcancers_GEOLUAD_glmprobes_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=t(rbind(TCGAtraindat,GEO_LUAD_x)),types=c(as.character(tumortype),rep("GEO_LUAD",nrow(GEO_LUAD_x))),main="Glmnet selected probes",probes=selectprobes1)
dev.off()
png(paste0(resfolder,"Allcancers_GEOLUAD_glmprobes_PCA_PC1_PC3.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=t(rbind(TCGAtraindat,GEO_LUAD_x)),types=c(as.character(tumortype),rep("GEO_LUAD",nrow(GEO_LUAD_x))),
         pc2=3,main="Glmnet selected probes",probes=selectprobes1)
dev.off()

alldat=form_xpredict1(dat2=t(GEO_LUAD))
GEO_LUAD_predictcancer_knn=knnclust()
table(GEO_LUAD_predictcancer_knn$class)
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 0         0        78         3         2         0 
#1se
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 0         0        80         1         2         0 
#1se best
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 0         0        79         2         2         0 

#LUSC
GEO_LUSC_x=form_xpredict(dat2=t(GEO_LUSC)) #0 probes not available
GEO_LUSC_predictprob=predict(cvfit, s=lambda.best,newx = GEO_LUSC_x, type = "response")
GEO_LUSC_predictcancer=predit_cancertype(allprobs=GEO_LUSC_predictprob)
table(GEO_LUSC_predictcancer$class)
# TCGA-LUAD TCGA-LUSC 
# 1        22 
#1se
# TCGA-LUAD TCGA-LUSC 
# 1        22 
#1se best
# TCGA-LUAD TCGA-LUSC 
# 1        22 

png(paste0(resfolder,"Allcancers_GEOLUSC_glmprobes_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=t(rbind(TCGAtraindat,GEO_LUSC_x)),types=c(as.character(tumortype),rep("GEO_LUSC",nrow(GEO_LUSC_x))),main="Glmnet selected probes",probes=selectprobes1)
dev.off()
png(paste0(resfolder,"Allcancers_GEOLUSC_glmprobes_PCA_PC1_PC3.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=t(rbind(TCGAtraindat,GEO_LUSC_x)),types=c(as.character(tumortype),rep("GEO_LUSC",nrow(GEO_LUSC_x))),
         pc2=3,main="Glmnet selected probes",probes=selectprobes1)
dev.off()

alldat=form_xpredict1(dat2=t(GEO_LUSC))
GEO_LUSC_predictcancer_knn=knnclust()
table(GEO_LUSC_predictcancer_knn$class)
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 0         0         3        20         0         0 
#1se
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 0         0         3        20         0         0 
#1se best
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 0         0         3        20         0         0 

#LUNG NORMAL
GEO_LUNG_NORM_x=form_xpredict(dat2=t(GEO_LUNG_NORM)) #0 probes not available
GEO_LUNG_NORM_predictprob=predict(cvfit, s=lambda.best,newx = GEO_LUNG_NORM_x, type = "response")
GEO_LUNG_NORM_predictcancer=predit_cancertype(allprobs=GEO_LUNG_NORM_predictprob)
table(GEO_LUNG_NORM_predictcancer$class)
# TCGA-NORM 
# 12 
#1se
# TCGA-NORM 
# 12
#1se best
# TCGA-NORM 
# 12 
alldat=form_xpredict1(dat2=t(GEO_LUNG_NORM))
GEO_LUNG_NORM_predictcancer_knn=knnclust()
table(GEO_LUNG_NORM_predictcancer_knn$class)
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 0         0         0         0        12         0 
#1se
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 0         0         0         0        12         0
#1se best
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-NORM TCGA-PRAD 
# 0         0         0         0        12         0
#do heatmap on TCGA data
library(ComplexHeatmap)
library(GetoptLong)
library(circlize)

uniq_cancertype=c(Cancertypes,"TCGA-NORM")
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
