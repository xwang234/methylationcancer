#!/usr/bin/env Rscript
#Breast cancer

#to run the code:
#ml fhR
#export OMP_NUM_THREADS=1
#salloc --constraint=gizmok -n 21 -t 0-4 mpirun -n 1 Rscript --no-save --no-restore ./TCGA_BRCA_methy_joint1.R
#salloc --constraint=gizmok -n 21 -t 1-5 mpirun -n 1 R --no-save --no-restore

#42 normals,370 tumors, 2 normals not have tumor
#Download data----------
setwd("/fh/fast/dai_j/CancerGenomics/prostate_methylation")
library(TCGAbiolinks) # If not installed, install from Github only!
library(SummarizedExperiment)
Cancer.type <- "TCGA-BRCA"  
RData.file <- paste("data/TCGA/",Cancer.type, "_", "methy.RData", sep='') # RData file for saving the TCGA data

print(paste0("Work on ",Cancer.type))
print("Load data------------")
if(!file.exists(RData.file)){
  
  Project <- Cancer.type
  #https://www.bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/query.html
  query.exp <- GDCquery(project = Project,
                        data.category = "DNA Methylation",
                        data.type = "Methylation Beta Value",
                        #workflow.type = workflow.type,
                        platform="Illumina Human Methylation 450",
                        legacy = FALSE
  )
  
  GDCdownload(query.exp)
  Data <- GDCprepare(query.exp, save = FALSE)
  
  Data.grangs <- Data@rowRanges
  Data.colData <- as.data.frame(Data@colData)
  Data <- assay(Data)
  
  Data.colData$sample_type <- factor(Data.colData$sample_type, levels = unique(Data.colData$sample_type))
  Data.colData <- Data.colData[order(Data.colData$sample_type),]
  Data.colData$barcode <- factor(Data.colData$barcode, levels = Data.colData$barcode)
  
  Data <- Data[,rownames(Data.colData)]
  
  
  query.clinic <- GDCquery(project = Project, 
                           data.category = "Clinical",
                           data.type = "Clinical Supplement", 
                           data.format = "BCR Biotab")
  GDCdownload(query.clinic)
  Data.clinic <- GDCprepare(query.clinic)
  
  save(Data, Data.grangs, Data.colData, Data.clinic, file = RData.file)
  #Data.grangs: annotation for all probes
}else
{
  #load data
  load(RData.file)
}

sum(is.na(Data))/nrow(Data)/ncol(Data) #0.185
#remove all NA probes (removed due to preprocessing)-----
library(matrixStats)
tmp=rowSds(Data)
idx= !is.na(tmp)
Data=Data[idx,]
dim(Data) #[1] 363780    892

pickuniqsample=function(dat=Data[, Data.colData$sample_type == "Primary Tumor"])
{
  tmp=substr(colnames(dat),1,16)
  idx=duplicated(tmp)
  dat=dat[,!idx]
  tmp=substr(colnames(dat),1,16)
  tmp1=substr(colnames(dat),1,12)
  idx1=which(!tmp1 %in% tmp1[duplicated(tmp1)])
  #samples with duplicates
  idx2=which(duplicated(tmp1))
  idxdup=NULL
  if (length(idx2)>0)
  {
    for (i in idx2)
    {
      idx=which(tmp1==tmp1[i])
      idx3=order(tmp[idx]) #pick 01A if having 01A and 01B
      idxdup=c(idxdup,idx[idx3[1]])
    }
  }
  idx=c(idx1,idxdup)
  dat=dat[,idx]
  colnames(dat)=substr(colnames(dat),1,12)
  return(dat)
}

library(cgwtools)
if (! exists("tumor_methy_all"))
{
  tumor_methy_all=pickuniqsample()
  dim(tumor_methy_all) #381223    370
  normal_methy=pickuniqsample(dat=Data[, Data.colData$sample_type == "Solid Tissue Normal"])
  dim(normal_methy) #381223     42
  all(colnames(normal_methy) %in% colnames(tumor_methy_all)) #F
  idx=match(colnames(normal_methy),colnames(tumor_methy_all))
  idx=idx[!is.na(idx)]
  tumor_methy=tumor_methy_all[,idx]
  all(colnames(normal_methy) == colnames(tumor_methy)) #F
  resave(normal_methy,tumor_methy,tumor_methy_all,file=RData.file)
}
if (! exists("clinical"))
{
  sum(is.na(Data.clinic$clinical_patient_brca$age_at_diagnosis)) #0
  sum(is.na(Data.clinic$clinical_patient_brca$gender)) #0
  sum(colnames(tumor_methy_all) %in% Data.clinic$clinical_patient_brca$bcr_patient_barcode) #782
  sum(colnames(normal_methy) %in% Data.clinic$clinical_patient_brca$bcr_patient_barcode) #96
  normalonly=colnames(normal_methy)[!colnames(normal_methy) %in% colnames(tumor_methy_all)]
  #use sample order in normalonly+tumor_methy_all
  clinical=data.frame(gender=rep(NA,length(normalonly)+ncol(tumor_methy_all)),age=NA)
  tmp=intersect(c(normalonly,colnames(tumor_methy_all)),Data.clinic$clinical_patient_brca$bcr_patient_barcode)
  idx1=match(tmp,c(normalonly,colnames(tumor_methy_all)))
  idx2=match(tmp,Data.clinic$clinical_patient_brca$bcr_patient_barcode)
  clinical$gender[idx1]=Data.clinic$clinical_patient_brca$gender[idx2]
  clinical$age[idx1]=Data.clinic$clinical_patient_brca$age_at_diagnosis[idx2]
  rownames(clinical)=c(normalonly,colnames(tumor_methy_all))
  #fill NAs in clinical
  sum(is.na(clinical$gender)) #0
  if (sum(is.na(clinical$gender))>0)
  {
    uniq_gender=unique(clinical$gender)
    num_gender=sapply(uniq_gender,function(x){sum(clinical$gender==x,na.rm=T)})
    clinical$gender[is.na(clinical$gender)]=uniq_gender[which.max(num_gender)]
  }
  clinical$age=as.integer(clinical$age)
  if (sum(is.na(clinical$age))>0)
    clinical$age[is.na(clinical$age)]=mean(as.numeric(clinical$age),na.rm=T)
  resave(clinical,file=RData.file)
}


print("Compute p-values-------------")
#function for each chunk, adjust for covariates in clinical
compute_pvalue <- function(r,n.chunk){
  #if there are extra normals, form paird data
  if (! all(colnames(normal_methy) %in% colnames(tumor_methy)))
  {
    tmp=intersect(colnames(tumor_methy),colnames(normal_methy))
    idx1=match(tmp,colnames(tumor_methy))
    idx2=match(tmp,colnames(normal_methy))
    tumor_methy1=tumor_methy[,idx1]
    normal_methy1=normal_methy[,idx2]
  }else
  {
    normal_methy1=normal_methy
    tumor_methy1=tumor_methy
  }
  #record output from each node
  #sink(paste0("outputnode",r,".txt"))
  # #to make sure same code running on each node by setting r=1
  # r=1
  # tmp=as.character(system("hostname",intern = T))
  # print(tmp)
  # tmp=as.character(system("lscpu | grep MHz",intern=T))
  # print(tmp)
  
  library(lawstat)
  
  nblock <- floor(nprobes/(n.chunk))
  ncount <- rep(nblock,n.chunk)
  if (nprobes>sum(ncount))  ncount[1:(nprobes-sum(ncount))] <- ncount[1:(nprobes-sum(ncount))]+ rep(1,nprobes-sum(ncount))
  chunks <- c(0,cumsum(ncount))
  chunkr <- seq(chunks[r]+1,chunks[r+1])
  nr <- length(chunkr)
  
  out1 <- matrix(NA,nr,25)
  rownames(out1)=rownames(normal_methy)[chunkr]
  # tmp=as.character(Sys.time())
  # print(paste0("start:",tmp))
  
  group <- c(rep(0,ncol(normal_methy1)),rep(1,ncol(tumor_methy1)))
  groupnames=c(colnames(normal_methy1),colnames(tumor_methy1))
  idx=match(groupnames,rownames(clinical))
  #covariate for tumor vs normal
  covariate=cbind(group,clinical[idx,])
  colnames(covariate)[2:ncol(covariate)]=colnames(clinical)
  covariate=as.data.frame(covariate)
  group <- c(rep(0,ncol(normal_methy)),rep(1,ncol(tumor_methy_all)))
  groupnames=c(colnames(normal_methy),colnames(tumor_methy_all))
  idx=match(groupnames,rownames(clinical))
  #covariate for alltumor vs normal
  covariateall=cbind(group,clinical[idx,])
  colnames(covariateall)[2:ncol(covariateall)]=colnames(clinical)
  covariateall=as.data.frame(covariateall)
  
  for (l in 1:nr) {
    i <- chunkr[l]
    out1[l,1] <- sd(as.numeric(normal_methy[i,]),na.rm=TRUE)
    out1[l,2] <- sd(as.numeric(tumor_methy[i,]),na.rm=TRUE)
    out1[l,3] <- sd(as.numeric(tumor_methy_all[i,]),na.rm=TRUE)
    
    
    out1[l,4] <- mean(as.numeric(normal_methy[i,]),na.rm=TRUE)
    out1[l,5] <- mean(as.numeric(tumor_methy[i,]),na.rm=TRUE)
    out1[l,6] <- mean(as.numeric(tumor_methy_all[i,]),na.rm=TRUE)
    out1[l,7] <- sd(c(as.numeric(normal_methy[i,]),as.numeric(tumor_methy_all[i,])),na.rm=T)
    
    ### compare the paired data ###
    
    y <- c(as.numeric(normal_methy1[i,]),as.numeric(tumor_methy1[i,]))  
    y <- ifelse(y<=0,0.01,y)
    y <- ifelse(y>=1,0.99,y)
    y <- log(y/(1-y),base=2)
    
    #group <- c(rep(0,ncol(normal_methy1)),rep(1,ncol(tumor_methy1)))
    idx=!is.na(y)
    covariate1=covariate[idx,]
    #make sure the left covariates have at least two levels
    num_uniqcovariate=sapply(1:ncol(covariate1), function(x) {length(unique(covariate1[,x]))})
    covariate1=covariate1[,num_uniqcovariate>1,drop=F]
    if ("group" %in% colnames(covariate1))
    {
      y <- y[idx] 
      lfit <- levene.test(y,covariate1$group,correction.method="zero.removal")
      out1[l,8] <- ifelse(is.numeric(lfit$p),lfit$p,NA)
      
      fit2 <-glm(y~.,data=covariate1,x=T,y=T)
      tmp=summary(fit2)$coef
      idx1=which(rownames(tmp)=="group")
      if (length(idx1)>0) out1[l,9] <- tmp[idx1,4]
      
      y1 <- y
      y1[covariate1$group==1] <- abs(y1[covariate1$group==1]-median(y1[covariate1$group==1]))
      y1[covariate1$group==0] <- abs(y1[covariate1$group==0]-median(y1[covariate1$group==0]))
      
      fit1 <-  glm(y1~.,data=covariate1,x=T,y=T)
      tmp=summary(fit1)$coef
      idx1=which(rownames(tmp)=="group")
      if (length(idx1)>0) out1[l,10] <- tmp[idx1,4]
      
      
      covmat <- solve(t(fit1$x) %*% fit1$x) %*% t(fit1$x*(fit1$y-fit1$fitted)) %*% (fit2$x*(fit2$y-fit2$fitted)) %*% solve(t(fit2$x) %*% fit2$x) 
      var1 <- solve(t(fit1$x) %*% fit1$x) %*% t(fit1$x*(fit1$y-fit1$fitted)) %*% (fit1$x*(fit1$y-fit1$fitted)) %*% solve(t(fit1$x) %*% fit1$x) 
      var2 <- solve(t(fit2$x) %*% fit2$x) %*% t(fit2$x*(fit2$y-fit2$fitted)) %*% (fit2$x*(fit2$y-fit2$fitted)) %*% solve(t(fit2$x) %*% fit2$x) 
      
      cov.gg <- matrix(0,2,2)
      cov.gg[1,2] <- covmat[2,2]
      cov.gg[2,1] <- covmat[2,2]
      cov.gg[1,1] <- var1[2,2]
      cov.gg[2,2] <- var2[2,2]
      
      coef.gg <- c(fit1$coef[2],fit2$coef[2])
      
      out1[l,11] <- 1-pchisq(drop(t(coef.gg) %*% solve(cov.gg) %*% coef.gg),df=2)
      out1[l,12] <- 1-pchisq(summary(fit1)$coeff[2,3]^2 + summary(fit2)$coeff[2,3]^2,df=2)
      
      zz <- max(summary(fit1)$coef[2,3],0)^2 + max(summary(fit2)$coef[2,3],0)^2
      out1[l,13] <- 0.5*pchisq(zz,df=1,lower.tail=F)+0.25*pchisq(zz,df=2,lower.tail=F)
      
      ff <- function(theta) {
        tt <- c(coef.gg[1]-theta[1], coef.gg[2]-theta[2])
        drop(t(tt) %*% solve(cov.gg) %*% tt) 
      }
      
      gr <- function(theta){
        -2*c(coef.gg[1]-theta[1], coef.gg[2]-theta[2])
      } 
      
      Amat <- matrix(c(1,0,0,1), 2, 2)
      Bmat <- c(0,0)
      
      #constrOptim can generate erros
      LRT=tryCatch(
        {
          drop(t(coef.gg) %*% solve(cov.gg) %*% coef.gg) - constrOptim(c(0.5,0.5),ff,gr,Amat,Bmat)$value
        },
        error=function(e)
        {
          return(F)
        }
      )
      #LRT <- drop(t(coef.gg) %*% solve(cov.gg) %*% coef.gg) - constrOptim(c(0.5,0.5),ff,gr,Amat,Bmat)$value
      
      pho <- cov.gg[1,2]/sqrt(cov.gg[1,1]*cov.gg[2,2])
      if (is.numeric(LRT))
      {
        out1[l,14] <- 0.5*pchisq(LRT,df=1,lower.tail=F)+0.5*(1-(1/pi)*acos(pho))*pchisq(LRT,df=2,lower.tail=F)
      }
      
      
      Amat <- matrix(c(1,0), 1, 2)
      Bmat <- c(0)
      
      #LRT <- drop(t(coef.gg) %*% solve(cov.gg) %*% coef.gg) - constrOptim(c(0.5,0),ff,gr,Amat,Bmat)$value 
      LRT=tryCatch(
        {
          drop(t(coef.gg) %*% solve(cov.gg) %*% coef.gg) - constrOptim(c(0.5,0),ff,gr,Amat,Bmat)$value
        },
        error=function(e)
        {
          return(F)
        }
      )
      if (is.numeric(LRT))
      {
        out1[l,15] <- 0.5*pchisq(LRT,df=1,lower.tail=F)+0.5*pchisq(LRT,df=2,lower.tail=F)
      }
      
      out1[l,16] <- pho
    }
    
    
    
    ### compare the all tumor vs normal data ###
    
    y <- c(as.numeric(normal_methy[i,]),as.numeric(tumor_methy_all[i,]))  
    y <- ifelse(y<=0,0.01,y)
    y <- ifelse(y>=1,0.99,y)
    y <- log(y/(1-y),base=2)
    # group <- c(rep(0,ncol(normal_methy)),rep(1,ncol(tumor_methy_all)))
    # group <- group[!is.na(y)]
    idx= !is.na(y)
    covariateall1=covariateall[idx,]
    num_uniqcovariate=sapply(1:ncol(covariateall1), function(x) {length(unique(covariateall1[,x]))})
    covariateall1=covariateall1[,num_uniqcovariate>1,drop=F]
    if ("group" %in% colnames(covariateall1))
    {
      y <- y[idx]
      lfit <- levene.test(y,covariateall1$group,correction.method="zero.removal")
      out1[l,17] <- ifelse(is.numeric(lfit$p),lfit$p,NA)
      
      fit2 <- glm(y~.,data=covariateall1,x=T,y=T)
      tmp=summary(fit2)$coef
      idx1=which(rownames(tmp)=="group")
      if (length(idx1)>0) out1[l,18] <- tmp[idx1,4]
      
      y1 <- y
      y1[covariateall1$group==1] <- abs(y1[covariateall1$group==1]-median(y1[covariateall1$group==1]))
      y1[covariateall1$group==0] <- abs(y1[covariateall1$group==0]-median(y1[covariateall1$group==0]))
      
      fit1 <- glm(y1~.,data=covariateall1,x=T,y=T)
      tmp=summary(fit1)$coef
      idx1=which(rownames(tmp)=="group")
      if (length(idx1)>0) out1[l,19] <- tmp[idx1,4]
      
      covmat <- solve(t(fit1$x) %*% fit1$x) %*% t(fit1$x*(fit1$y-fit1$fitted)) %*% (fit2$x*(fit2$y-fit2$fitted)) %*% solve(t(fit2$x) %*% fit2$x) 
      var1 <- solve(t(fit1$x) %*% fit1$x) %*% t(fit1$x*(fit1$y-fit1$fitted)) %*% (fit1$x*(fit1$y-fit1$fitted)) %*% solve(t(fit1$x) %*% fit1$x) 
      var2 <- solve(t(fit2$x) %*% fit2$x) %*% t(fit2$x*(fit2$y-fit2$fitted)) %*% (fit2$x*(fit2$y-fit2$fitted)) %*% solve(t(fit2$x) %*% fit2$x) 
      
      cov.gg <- matrix(0,2,2)
      cov.gg[1,2] <- covmat[2,2]
      cov.gg[2,1] <- covmat[2,2]
      cov.gg[1,1] <- var1[2,2]
      cov.gg[2,2] <- var2[2,2]
      
      coef.gg <- c(fit1$coef[2],fit2$coef[2])
      
      out1[l,20] <- 1-pchisq(drop(t(coef.gg) %*% solve(cov.gg) %*% coef.gg),df=2)
      out1[l,21] <- 1-pchisq(summary(fit1)$coeff[2,3]^2 + summary(fit2)$coeff[2,3]^2,df=2)
      
      zz <- max(summary(fit1)$coef[2,3],0)^2 + max(summary(fit2)$coef[2,3],0)^2
      out1[l,22] <- 0.5*pchisq(zz,df=1,lower.tail=F)+0.25*pchisq(zz,df=2,lower.tail=F)
      
      ff <- function(theta) {
        tt <- c(coef.gg[1]-theta[1], coef.gg[2]-theta[2]) 
        drop(t(tt) %*% solve(cov.gg) %*% tt)  
      }
      
      gr <- function(theta){
        -2*c(coef.gg[1]-theta[1], coef.gg[2]-theta[2])
      }  
      
      Amat <- matrix(c(1,0,0,1), 2, 2)
      Bmat <- c(0,0)
      
      #LRT <- drop(t(coef.gg) %*% solve(cov.gg) %*% coef.gg) - constrOptim(c(0.5,0.5),ff,gr,Amat,Bmat)$value
      LRT=tryCatch(
        {
          drop(t(coef.gg) %*% solve(cov.gg) %*% coef.gg) - constrOptim(c(0.5,0.5),ff,gr,Amat,Bmat)$value
        },
        error=function(e)
        {
          return(F)
        }
      )
      
      pho <- cov.gg[1,2]/sqrt(cov.gg[1,1]*cov.gg[2,2])
      if (is.numeric(LRT))
      {
        out1[l,23] <- 0.5*pchisq(LRT,df=1,lower.tail=F)+0.5*(1-(1/pi)*acos(pho))*pchisq(LRT,df=2,lower.tail=F)
      }
      
      Amat <- matrix(c(1,0), 1, 2)
      Bmat <- c(0)
      
      #LRT <- drop(t(coef.gg) %*% solve(cov.gg) %*% coef.gg) - constrOptim(c(0.5,0),ff,gr,Amat,Bmat)$value 
      LRT=tryCatch(
        {
          drop(t(coef.gg) %*% solve(cov.gg) %*% coef.gg) - constrOptim(c(0.5,0),ff,gr,Amat,Bmat)$value
        },
        error=function(e)
        {
          return(F)
        }
      )
      if (is.numeric(LRT))
      {
        out1[l,24] <- 0.5*pchisq(LRT,df=1,lower.tail=F)+0.5*pchisq(LRT,df=2,lower.tail=F)
      }
      
      out1[l,25] <- pho
    }

    #check how long time spent on 20 interations
    # if (l %% 200==0)
    # {
    #   cat(l,'..')
    #   tmp=as.character(Sys.time())
    #   print(paste0("time:",tmp))
    # }
  }
  # print("done")
  # tmp=as.character(Sys.time())
  # print(paste0("end:",tmp))
  #sink()
  list(out1)
}

#genome-wide computation--------------


library(Rmpi)

#sink("MPI_TCGA_output.txt")

### Prepare Rmpi environment, DO NOT CHANGE ###
Sys.getenv(c("SLURM_SUBMIT_DIR"))
Sys.getenv(c("HOST", "SLURM_JOB_ID", "SLURM_NODELIST","SLURM_NNODES", "SLURM_NTASKS",
             "SLUR M_CPUS_PER_TASK","SLURM_CPUS_ON_NODE","SLURM_NTASKS_PER_NODE",
             "SLURM_TASK_PID", "SLURM_ PARTITION"))
#

nWorkers <- mpi.universe.size() - 1L
njob=nWorkers
mpi.spawn.Rslaves(nslaves=nWorkers,needlog=F)
.Last <- function(){
  if (is.loaded("mpi_initialize")){
    if (mpi.comm.size(1) > 0){
      print("Please use mpi.close.Rslaves() to close slaves.")
      mpi.close.Rslaves()
    }
    print("Please use mpi.quit() to quit R")
    .Call("mpi_finalize")
  }
}


starttime<-Sys.time();starttime
mpi.bcast.Robj2slave(clinical)
mpi.remote.exec(library(lawstat))


mpi.bcast.Robj2slave(compute_pvalue)

#tumor_methy_all is too big, has to split chr

all_normal_methy=normal_methy
all_tumor_methy=tumor_methy
all_tumor_methy_all=tumor_methy_all

library(GenomicRanges)
allchrs=as.character(seqnames(Data.grangs)@values)
allprobenames=as.character(names(Data.grangs))

allout=NULL
tmp=as.character(Sys.time())
print(paste0("start:",tmp))
for (chr in unique(allchrs))
{
  cat(chr,'..')
  idx=rownames(all_normal_methy) %in% allprobenames[allchrs==chr]
  if (sum(idx)>0)
  {
    normal_methy=all_normal_methy[idx,]
    tumor_methy=all_tumor_methy[idx,]
    tumor_methy_all=all_tumor_methy_all[idx,]
    mpi.bcast.Robj2slave(normal_methy)
    mpi.bcast.Robj2slave(tumor_methy)
    mpi.bcast.Robj2slave(tumor_methy_all)
    nprobes <- nrow(normal_methy)
    mpi.bcast.Robj2slave(nprobes)
    
    #split jobs into 100 chunks
    outlist<-mpi.parSapply(1:100,FUN=compute_pvalue,n.chunk=100,job.num=njob)
    out <- do.call(rbind,outlist)
    allout=rbind(allout,out)
  }
}
save(allout,file="tmp.RData")
tmp=as.character(Sys.time())
print(paste0("end:",tmp))
output <- data.frame(allout)
idx=match(rownames(all_normal_methy),rownames(output))
output=output[idx,]
# if (nrow(output)==nrow(all_normal_methy))
#   row.names(output) <- row.names(all_normal_methy)

save(output,file=paste0("/fh/fast/dai_j/CancerGenomics/prostate_methylation/result/",Cancer.type,"_output.RData"))

## quit program 
mpi.close.Rslaves()
mpi.quit()
quit()
