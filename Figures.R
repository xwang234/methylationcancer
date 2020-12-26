#Function used to plot 4 figures for DVC paper (normals vs all tumors)

library(pROC)
library(scales) #alpha

Cancer.type="TCGA-PRAD"
Cancer.type="TCGA-COAD"
Cancer.type="TCGA-LUAD"
Cancer.type="TCGA-LUSC"
Cancer.type="TCGA-BRCA"
Cancer.type="Krause-EAC"
Cancer.type="TCGA-LIHC"
rm(tcga_normal_methy,tcga_tumor_methy,tcga_tumor_methy_all,output)


#Load result data
resfolder="/fh/fast/dai_j/CancerGenomics/prostate_methylation/result/"
Res.file=paste0(resfolder,Cancer.type,"_output.RData")
load(Res.file)

#load methylation data
if (Cancer.type=="TCGA-PRAD")
{
  load("/fh/fast/dai_j/CancerGenomics/prostate_methylation/tcga_methy.RData")
  #load("tcga_clinical.RData")
}else
{
  RData.file=paste("data/TCGA/",Cancer.type, "_", "methy.RData", sep='')
  load(RData.file)
}

if (!exists("tcga_normal_methy"))
{
  tcga_normal_methy=normal_methy
  tcga_tumor_methy=tumor_methy
  tcga_tumor_methy_all=tumor_methy_all
}

if (any(rownames(output)!=rownames(tcga_normal_methy)))
{
  stop("something is wrong!")
}

# dat.median=data.frame(tumor_median=rowMedians(as.matrix(tcga_tumor_methy_all),na.rm=T),normal_median=rowMedians(as.matrix(tcga_normal_methy),na.rm=T))
# 
# idx=which(abs(output[,3]-output[,1])<0.02)
# x=output[idx,3]-output[idx,1]
# y=-log10(output[idx,19])
# idx1=output[idx,19]<0.05/nrow(output)
# colors=c("black","red")
# plot(x,y,col=colors[as.numeric(idx1)+1],xlab="Difference of SD (beta)",ylab=expression('-log'[10]*' p-value for DVC'),main=Cancer.type)
# idx=which(abs(output[,3]-output[,1])<0.02 & output[,19]<0.05/nrow(output))
# par(mfrow=c(2,1))
# hist(dat.median$normal_median[idx],main="Normal",xlab="Beta value",breaks=20)
# hist(dat.median$tumor_median[idx],main="Tumor",xlab="Beta value",breaks=20)
# plot(dat.median$normal_median[idx],dat.median$tumor_median[idx],xlab="Median of normal beta",ylab="Median of tumor beta")


#flter probes sd (c(cancer,normal))>0.03 for figures 1 and 2
alldat=cbind.data.frame(tcga_normal_methy,tcga_tumor_methy_all)
library(matrixStats)
alldat.rowsd=rowSds(as.matrix(alldat),na.rm=T)
quantile(alldat.rowsd)
#Figure 1,volcanoplt,manually adjust ylim
output1=output[alldat.rowsd>0.05,]
idx=output1[,19]<0.05/nrow(output)
mycolors=rep("black",nrow(output1))
mycolors[idx]="red"
table(mycolors)
plot(output1[,3]-output1[,1],-log(output1[,19],base=10),xlab="Difference of SD (beta)",ylab=expression('-log'[10]*' p-value for DVC'),
     col=alpha(mycolors, 0.05),bty='l',cex.axis=1.4,cex.lab=1.4)
# #LUAD check
# idx=which(abs(output[,3]-output[,1])<0.02 & output[,19]<0.05/nrow(output))
# idx=idx[order(output[idx,19])]
# View(output[idx,])
# rownames(output)[idx[1]] #cg20592884,"cg10457846"
png(paste0(resfolder,Cancer.type,"_volcanoplot.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
#p-value came from DVC (fit1)
plot(output1[,3]-output1[,1],-log(output1[,19],base=10),xlab="Difference of SD (beta)",ylab=expression('-log'[10]*' p-value for DVC'),
     col=alpha(mycolors, 0.05),bty='l',cex.axis=1.4,cex.lab=1.4,ylim=c(0,50)) 
dev.off()

#Figure 2, Difference of mean vs Difference of SD, pick significant points,borders need to be removed.
smoothScatter1=function (x, y = NULL, nbin = 128, bandwidth, colramp = colorRampPalette(c("white", 
                                                                           blues9)), nrpoints = 100, ret.selection = FALSE, pch = ".", 
          cex = 1, col = "black", transformation = function(x) x^0.25, 
          postPlotHook = box, xlab = NULL, ylab = NULL, xlim, ylim, 
          xaxs = par("xaxs"), yaxs = par("yaxs"), ...) 
{
  if (!is.numeric(nrpoints) || nrpoints < 0 || length(nrpoints) != 1) 
    stop("'nrpoints' should be numeric scalar with value >= 0.")
  nrpoints <- round(nrpoints)
  ret.selection <- ret.selection && nrpoints > 0
  xlabel <- if (!missing(x)) deparse1(substitute(x))
  ylabel <- if (!missing(y)) deparse1(substitute(y))
  xy <- xy.coords(x, y, xlabel, ylabel)
  xlab <- if (is.null(xlab)) xy$xlab else xlab
  ylab <- if (is.null(ylab)) xy$ylab else ylab
  x <- cbind(xy$x, xy$y)[I <- is.finite(xy$x) & is.finite(xy$y), 
                         , drop = FALSE]
  if (ret.selection) iS <- which(I)
  if (!missing(xlim)) {
    stopifnot(is.numeric(xlim), length(xlim) == 2, is.finite(xlim))
    x <- x[I <- min(xlim) <= x[, 1] & x[, 1] <= max(xlim), 
           , drop = FALSE]
    if (ret.selection) iS <- iS[I]
  }else {
    xlim <- range(x[, 1])
  }
  if (!missing(ylim)) {
    stopifnot(is.numeric(ylim), length(ylim) == 2, is.finite(ylim))
    x <- x[I <- min(ylim) <= x[, 2] & x[, 2] <= max(ylim), 
           , drop = FALSE]
    if (ret.selection) iS <- iS[I]
  }else {
    ylim <- range(x[, 2])
  }
  map <- grDevices:::.smoothScatterCalcDensity(x, nbin, bandwidth)
  xm <- map$x1
  ym <- map$x2
  dens <- map$fhat
  dens[] <- transformation(dens)
  image(xm, ym, z = dens, col = colramp(256), xlab = xlab, 
        ylab = ylab, xlim = xlim, ylim = ylim, xaxs = xaxs, yaxs = yaxs, axes=F,...)
  #-0.2 is added to make sure x-axis and y-axis cross
  tmp1=round(min(xm),1)-0.2
  tmp2=round(max(xm),1)
  print(c(tmp1,tmp2))
  #seq(-0.6,0.6,0.2) generated weird series
  tmp3=tmp1
  while(tmp3[length(tmp3)]<tmp2 & tmp1 < tmp2)
  {
    tmp3=c(tmp3,round(tmp3[length(tmp3)]+0.2,1))
  }
  axis(1, at = tmp3,xpd=F,pos=0,...)
  #print(tmp3)
  axis(2, at = round(seq(0, max(ym), max(ym)/5),1),xpd = T,...)
  print(round(seq(0, max(ym), max(ym)/5),1))
  #print(seq(tmp1,tmp2,0.2))
  #print(round(seq(0, max(ym), length.out=5),1))
  #axis(1,at=NULL,xlim=c(min(xm),max(xm)),xpd = TRUE,...)
  #axis(2,at=NULL, xpd = TRUE,...)
  
  
  if (!is.null(postPlotHook)) postPlotHook()
  if (nrpoints > 0) {
    nrpoints <- min(nrow(x), ceiling(nrpoints))
    stopifnot((nx <- length(xm)) == nrow(dens), (ny <- length(ym)) == 
                ncol(dens))
    ixm <- 1L + as.integer((nx - 1) * (x[, 1] - xm[1])/(xm[nx] - 
                                                          xm[1]))
    iym <- 1L + as.integer((ny - 1) * (x[, 2] - ym[1])/(ym[ny] - 
                                                          ym[1]))
    sel <- order(dens[cbind(ixm, iym)])[seq_len(nrpoints)]
    x <- x[sel, , drop = FALSE]
    points(x, pch = pch, cex = cex, col = col)
    if (ret.selection) 
      iS[sel]
  }
}

#pick significant & y>0
idx=output1[,19]<0.05/nrow(output) & output1[,3]-output1[,1]>0
if (Cancer.type=="TCGA-LUAD")
{
  png(paste0(resfolder,Cancer.type,"_DVC_mean_SD.png"),width = 480, height = 480,type = "cairo")
  par(mar=c(6,6,2,1))
  smoothScatter1(output1[idx,6]-output1[idx,4],output1[idx,3]-output1[idx,1],postPlotHook=NULL,
                 xlab="",ylab="Difference of SD beta",
                 cex.axis=1.4,cex.lab=1.4,useRaster = TRUE,xaxs="i",yaxs="i") 
  mtext("Difference of mean beta", side=1, line=5,cex=1.4)
  dev.off()
}else
{
  png(paste0(resfolder,Cancer.type,"_DVC_mean_SD.png"),width = 480, height = 480,type = "cairo")
  par(mar=c(6,6,2,1))
  smoothScatter1(output1[idx,6]-output1[idx,4],output1[idx,3]-output1[idx,1],postPlotHook=NULL,
                 xlab="Difference of mean beta",ylab="Difference of SD beta",
                 cex.axis=1.4,cex.lab=1.4,useRaster = TRUE,xaxs="i",yaxs="i") 
  dev.off()
}
if (Cancer.type=="TCGA-LIHC") #ylim
{
  png(paste0(resfolder,Cancer.type,"_DVC_mean_SD.png"),width = 480, height = 480,type = "cairo")
  par(mar=c(6,6,2,1))
  smoothScatter1(output1[idx,6]-output1[idx,4],output1[idx,3]-output1[idx,1],postPlotHook=NULL,
                 xlab="Difference of mean beta",ylab="Difference of SD beta",
                 cex.axis=1.4,cex.lab=1.4,useRaster = TRUE,xaxs="i",yaxs="i",ylim=c(0,0.45)) 
  dev.off()
}

#Figure 3 boxplot
## now I want to find the probes with signficant DVC, but not signficant DMC, also with low methylation in normal

dat.median=data.frame(tumor_median=rowMedians(as.matrix(tcga_tumor_methy_all),na.rm=T),normal_median=rowMedians(as.matrix(tcga_normal_methy),na.rm=T),
                      tumor15=rowQuantiles(as.matrix(tcga_tumor_methy_all),probs = 0.15,na.rm=T),normal85=rowQuantiles(as.matrix(tcga_normal_methy),probs = 0.85,na.rm=T))
quantile(dat.median$normal_median)
sum(dat.median$normal_median<0.1)/nrow(output)
#PRAD
# 0%        25%        50%        75%       100% 
# 0.00638648 0.05569549 0.45380698 0.86799606 0.99429863 
#0.3375153
#COAD
# 0%         25%         50%         75%        100% 
# 0.006794376 0.059397252 0.477026244 0.851353135 0.992764957
#0.3187608
#LUAD
# 0%         25%         50%         75%        100% 
# 0.007997009 0.093078985 0.492056986 0.808729806 0.992103931 
#0.2690472
#LUSC
# 0%         25%         50%         75%        100% 
# 0.005477912 0.047120469 0.470393673 0.854058505 0.994802667
#0.3488116
#BRCA
# 0%         25%         50%         75%        100% 
# 0.005887769 0.051843180 0.486208541 0.854245675 0.994524477 
#0.3360767
#EAC
# 0%          25%          50%          75%         100% 
# 0.0001726538 0.0495650712 0.5048604375 0.8827488822 0.9897341674 
#0.317174
#LIHC
# 0%         25%         50%         75%        100% 
# 0.006878277 0.060714521 0.516979408 0.863027547 0.993677376 
# 0.3130814
clist <- which(output[,19]<0.05/nrow(output) & output[,18]>0.05/nrow(output) & dat.median$normal_median<0.1 & dat.median$tumor_median>dat.median$normal_median)

sum(output[,19]<0.05/nrow(output) & output[,18]>0.05/nrow(output))
length(clist)
#PRAD
#5068,1836
#COAD
#29192,2127
#LUAD
#15606,954
#BRCA
#43508,1562
#LIHC
#62705,3642

## here I pick the one with largest AUC, ideally with similar AUC I want a low DVC p-value, a big DMC p-value to show testing for DVC is meaningful
aucres=data.frame(i=clist,outauc= rep(0,length(clist)),p_dmc=NA,p_dvc=NA,med_norm=NA,med_tumor=NA,norm85=NA,tumor15=NA)
for (k in 1:length(clist)) {
  i <- clist[k]
  y <- c(as.numeric(tcga_normal_methy[i,]),as.numeric(tcga_tumor_methy_all[i,]))
  group <- c(rep(0,length(as.numeric(tcga_normal_methy[i,]))),rep(1,length(as.numeric(tcga_tumor_methy_all[i,]))))
  aucres$outauc[k] <- auc(roc(group~y))
  aucres$p_dmc[k]=output[i,18]
  aucres$p_dvc[k]=output[i,19]
  aucres$med_norm[k]=dat.median$normal_median[i]
  aucres$med_tumor[k]=dat.median$tumor_median[i]
  aucres$norm85[k]=dat.median$normal85[i]
  aucres$tumor15[k]=dat.median$tumor15[i]
}
aucres$ratio=aucres$p_dmc/aucres$p_dvc
aucres$med_diff=aucres$med_tumor-aucres$med_norm
# idx1=which(aucres$outauc>0.7)
# plot(-log10(aucres$p_dvc[idx1]),aucres$outauc[idx1],col=alpha("red",0.5),xlab="-log10 p-value",xlim=c(0,max(-log10(aucres$p_dvc))),ylab="AUC")
# points(-log10(aucres$p_dmc[idx1]),aucres$outauc[idx1],col=alpha("blue",0.5))
# legend("topright",legend=c("DMC","DVC"),col=c("blue","red"),pch=1)
# max(aucres$outauc)
#pick different auc cutoff to pick probe
aucres1=aucres[aucres$outauc>0.65,]
aucres1=aucres[aucres$outauc>0.7,]
idx1=order(aucres1$p_dmc,decreasing = T)
#View(aucres1[idx1,])
#plot(-log10(aucres1$p_dmc),-log10(aucres1$p_dvc),xlab="-log10 p-value (DMC)",ylab="-log10 p-value (DVC)",main="AUC>0.7")
# dat=data.frame(dmc=-log10(aucres1$p_dmc),dvc=-log10(aucres1$p_dvc),auc=aucres1$outauc)
# sp<-ggplot(dat, aes(x=dmc, y=dvc, color=auc)) + geom_point()+ theme_bw()+
#   labs(title="932 probes with AUC>0.7",x ="-log10 p-value (DMC)", y = "-log10 p-value (DVC)")+
#   theme(text = element_text(size=20)) 
# sp+scale_color_gradient(low="blue", high="red")

#plot figures3 and 4 based on different criteria
#ROC----
plotroc1=function(predict=y,ans=group,main="",xpos=0.45,ypos=0.65)
{
  idx=!is.na(predict)
  predict=predict[idx]
  ans=ans[idx]
  predict1=predict[order(predict)]
  TP=FP=rep(0,length(ans))
  if (mean(predict[ans==1]) >mean(predict[ans==0]))
  {
    for (i in 1:length(ans)) {
      TP[i] <- mean(predict[ans==1]>=predict1[i]) 
      FP[i] <- mean(predict[ans==0]>=predict1[i]) 
    }   
  }else
  {
    for (i in 1:length(ans)) {
      TP[i] <- mean(predict[ans==1]<=predict1[i]) 
      FP[i] <- mean(predict[ans==0]<=predict1[i]) 
    }
  }
  
  plot(c(FP),c(TP),type="l",bty='l',xlim=c(0,1),ylim=c(0,1),lwd=3,xlab="1-specificity",ylab="Sensitivity",main=main,cex.lab=1.4,cex.axis=1.4,cex.main=1.4)
  abline(0,1,lty=2)
  text(xpos,ypos,paste0("AUC=",round(auc(roc(ans~predict)),3)),cex=1.4)
}

plot_figure3_4=function(i,prefix,x1=0.45,y1=0.65)
{
  png(paste0(resfolder,Cancer.type,"_boxplot_DVC_",prefix,".png"),width = 480, height = 480,type = "cairo")
  par(mar=c(6,6,2,1))
  y <- c(as.numeric(tcga_normal_methy[i,]),as.numeric(tcga_tumor_methy_all[i,]))
  group <- c(rep(0,length(as.numeric(tcga_normal_methy[i,]))),rep(1,length(as.numeric(tcga_tumor_methy_all[i,]))))
  boxplot(y~group,ylim=c(0,max(y,na.rm=T)),boxlty = 0,whisklty = 0,staplelty = 0,xlab="",ylab="Methylation beta",names=rep("",2),
          outline=F,col=NULL,main=rownames(tcga_normal_methy)[i],frame.plot = FALSE,cex.axis=1.4,cex.lab=1.4,cex.main=1.4)
  axis(side = 1,at=c(0,1,2,3),labels=c("","Normal","Cancer",""),cex.axis=1.4,cex.lab=1.4)
  stripchart(y~group,vertical = TRUE, method = "jitter",jitter=0.3,pch = 1,col=c(3,2),add = TRUE)
  tmp=par("usr")
  idx=match(rownames(tcga_normal_methy)[i],rownames(output))
  tmp1=formatC(output[idx,18], format = "e", digits = 2)
  tmp2=formatC(output[idx,19], format = "e", digits = 2)
  text((tmp[1]+tmp[2])*0.4,0.92*max(y,na.rm=T),paste0("p(DMC)=",tmp1),cex=1.4)
  text((tmp[1]+tmp[2])*0.4,0.85*max(y,na.rm=T),paste0("p(DVC)=",tmp2),cex=1.4)
  dev.off()
  
  #Figure 4: ROC
  png(paste0(resfolder,Cancer.type,"_AUC_DVC_",prefix,".png"),width = 480, height = 480,type = "cairo")
  par(mar=c(6,6,2,1))
  #plot(roc(group~y),main=rownames(tcga_normal_methy)[i],print.auc=F,cex.axis=1.5,cex.main=1.5,cex.lab=1.5)
  #text(0.55,0.6,paste0("AUC=",round(auc(roc(group~y)),3)),cex=1.4)
  # text(0.8,0.4,paste0("AUC=",round(auc(roc(group~y)),3)),cex=1.4)
  plotroc1(predict=y,ans=group,main=rownames(tcga_normal_methy)[i],xpos=x1,ypos=y1)
  dev.off()
}

#Plot probe with maxAUC
plot_figure3_4(i=aucres1$i[which.max(aucres1$outauc)],prefix="maxAUC",x1=0.3,y1=0.5)
#Plot probe with max DMC p-value
plot_figure3_4(i=aucres1$i[which.max(aucres1$p_dmc)],prefix="maxDMC")
#Plot probe with max ratio of p-values DMC/DVC
plot_figure3_4(i=aucres1$i[which.max(aucres1$ratio)],prefix="maxpvaluediff")
plot_figure3_4(i=aucres1$i[which.max(aucres1$med_diff)],prefix="maxmeddiff")

#only for a subset of probes with small median normal and small norm 85 percentiles
cutoff1=quantile(aucres1$med_norm,0.5)
cutoff2=quantile(aucres1$norm85,0.15)

idx1=which(aucres1$med_norm<=cutoff1 & aucres1$norm85<=cutoff2)
length(idx1)
plot_figure3_4(i=aucres1$i[idx1][which.max(aucres1$outauc[idx1])],prefix="subset_maxAUC")
plot_figure3_4(i=aucres1$i[idx1][which.max(aucres1$p_dmc[idx1])],prefix="subset_maxDMC",x1=0.3,y1=0.5)
plot_figure3_4(i=aucres1$i[idx1][which.max(aucres1$ratio[idx1])],prefix="subset_maxpvaluediff",x1=0.3,y1=0.5)

#function including density plot, for specific cpg
plot_figure3_5=function(probeid="cg23351584",x1=0.45,y1=0.6)
{
  i=which(rownames(tcga_normal_methy)==probeid)
  png(paste0(resfolder,Cancer.type,"_boxplot_DVC_",probeid,".png"),width = 480, height = 480,type = "cairo")
  par(mar=c(6,6,2,1))
  y <- c(as.numeric(tcga_normal_methy[i,]),as.numeric(tcga_tumor_methy_all[i,]))
  group <- c(rep(0,length(as.numeric(tcga_normal_methy[i,]))),rep(1,length(as.numeric(tcga_tumor_methy_all[i,]))))
  boxplot(y~group,ylim=c(0,max(y,na.rm=T)),boxlty = 0,whisklty = 0,staplelty = 0,xlab="",ylab="Methylation beta",names=rep("",2),
          outline=F,col=NULL,main=rownames(tcga_normal_methy)[i],frame.plot = FALSE,cex.axis=1.4,cex.lab=1.4,cex.main=1.4)
  axis(side = 1,at=c(0,1,2,3),labels=c("","Normal","Cancer",""),cex.axis=1.4,cex.lab=1.4)
  stripchart(y~group,vertical = TRUE, method = "jitter",jitter=0.3,pch = 1,col=c(3,2),add = TRUE)
  tmp=par("usr")
  idx=match(rownames(tcga_normal_methy)[i],rownames(output))
  tmp1=formatC(output[idx,18], format = "e", digits = 2)
  tmp2=formatC(output[idx,19], format = "e", digits = 2)
  text((tmp[1]+tmp[2])*0.4,0.92*max(y,na.rm=T),paste0("p(DMC)=",tmp1),cex=1.4)
  text((tmp[1]+tmp[2])*0.4,0.85*max(y,na.rm=T),paste0("p(DVC)=",tmp2),cex=1.4)
  dev.off()
  
  png(paste0(resfolder,Cancer.type,"_AUC_DVC_",probeid,".png"),width = 480, height = 480,type = "cairo")
  par(mar=c(6,6,2,1))
  plotroc1(predict=y,ans=group,main=rownames(tcga_normal_methy)[i],xpos=x1,ypos=y1)
  dev.off()
  
  
  y1 <- log(as.numeric(tcga_tumor_methy_all[i,])/(1-as.numeric(tcga_tumor_methy_all[i,])),base=2)
  y0 <- log(as.numeric(tcga_normal_methy[i,]),base=2)

  
  
  
  ## you may need to adjust bandwidth (bw) to see smooth lines, also adjust ylim, possibly legend position
  png(paste0(resfolder,Cancer.type,"_density1_",probeid,".png"),width = 480, height = 480,type = "cairo")
  par(mar=c(6,6,2,1))
  bw1=bw.nrd(y1[!is.na(y1)])
  bw0=bw.nrd(y0[!is.na(y0)])
  tmp1=density(y1[!is.na(y1)],bw=bw1+0.1)
  tmp0=density(y0[!is.na(y0)],bw=bw0+0.1)

  
  xlim=c(min(tmp1$x,tmp0$x),max(c(tmp1$x,tmp0$x)))
  ylim=c(min(tmp1$y,tmp0$y),max(c(tmp1$y,tmp0$y)))
  plot(tmp1,col=2,lwd=4,main="",xlab="Methylation M value",bty='l',xlim=xlim,ylim=ylim,cex.lab=1.4,cex.axis=1.4)
  lines(tmp0,col=3,lwd=4)
  legend("topright",c("Normal","Cancer"),lwd=4,col=c(3,2),cex=1.4,bty="n")
  dev.off()
  
  png(paste0(resfolder,Cancer.type,"_density2_",probeid,".png"),width = 480, height = 480,type = "cairo")
  par(mar=c(6,6,2,1))
  bw1=bw.nrd(y1[!is.na(y1)])
  bw0=bw.nrd(y0[!is.na(y0)])
  tmp1=density(y1[!is.na(y1)],bw=bw1+0.2)
  tmp0=density(y0[!is.na(y0)],bw=bw0+0.2)
  
  
  xlim=c(min(tmp1$x,tmp0$x),max(c(tmp1$x,tmp0$x)))
  ylim=c(min(tmp1$y,tmp0$y),max(c(tmp1$y,tmp0$y)))
  plot(tmp1,col=2,lwd=4,main="",xlab="Methylation M value",bty='l',xlim=xlim,ylim=ylim,cex.lab=1.4,cex.axis=1.4)
  lines(tmp0,col=3,lwd=4)
  legend("topright",c("Normal","Cancer"),lwd=4,col=c(3,2),cex=1.4,bty="n")
  dev.off()
  
}

plot_figure3_5=function(probeid="cg23351584",xpos=0.45,ypos=0.6)
{
  i=which(rownames(tcga_normal_methy)==probeid)
  png(paste0(resfolder,Cancer.type,"_",probeid,"_boxplot_DVC",".png"),width = 480, height = 480,type = "cairo")
  par(mar=c(6,6,2,1))
  y <- c(as.numeric(tcga_normal_methy[i,]),as.numeric(tcga_tumor_methy_all[i,]))
  group <- c(rep(0,length(as.numeric(tcga_normal_methy[i,]))),rep(1,length(as.numeric(tcga_tumor_methy_all[i,]))))
  boxplot(y~group,ylim=c(0,max(y,na.rm=T)),boxlty = 0,whisklty = 0,staplelty = 0,xlab="",ylab="Methylation beta",names=rep("",2),
          outline=F,col=NULL,main=rownames(tcga_normal_methy)[i],frame.plot = FALSE,cex.axis=1.4,cex.lab=1.4,cex.main=1.4)
  axis(side = 1,at=c(0,1,2,3),labels=c("","Normal","Cancer",""),cex.axis=1.4,cex.lab=1.4)
  stripchart(y~group,vertical = TRUE, method = "jitter",jitter=0.3,pch = 1,col=c(3,2),add = TRUE)
  tmp=par("usr")
  idx=match(rownames(tcga_normal_methy)[i],rownames(output))
  tmp1=formatC(output[idx,18], format = "e", digits = 2)
  tmp2=formatC(output[idx,19], format = "e", digits = 2)
  text((tmp[1]+tmp[2])*0.4,0.92*max(y,na.rm=T),paste0("p(DMC)=",tmp1),cex=1.4)
  text((tmp[1]+tmp[2])*0.4,0.85*max(y,na.rm=T),paste0("p(DVC)=",tmp2),cex=1.4)
  dev.off()
  
  png(paste0(resfolder,Cancer.type,"_",probeid,"_AUC_DVC",".png"),width = 480, height = 480,type = "cairo")
  par(mar=c(6,6,2,1))
  plotroc1(predict=y,ans=group,main=rownames(tcga_normal_methy)[i],xpos=xpos,ypos=ypos)
  dev.off()
  
  
  ## you may need to adjust bandwidth (bw) to see smooth lines, also adjust ylim, possibly legend position
  tmp1=as.numeric(tcga_tumor_methy_all[i,])
  tmp1=ifelse(tmp1<=0,0.001,tmp1)
  tmp1 <- ifelse(tmp1>=1,0.999,tmp1)
  y1=log(tmp1/(1-tmp1),base=2)
  tmp0=as.numeric(tcga_normal_methy[i,])
  tmp0=ifelse(tmp0<=0,0.001,tmp0)
  tmp0 <- ifelse(tmp0>=1,0.999,tmp0)
  y0=log(tmp0/(1-tmp0),base=2)
  # y1 <- log(as.numeric(tcga_tumor_methy_all[i,])/(1-as.numeric(tcga_tumor_methy_all[i,])),base=2)
  # y0 <- log(as.numeric(tcga_normal_methy[i,]),base=2)
  png(paste0(resfolder,Cancer.type,"_",probeid,"_density1",".png"),width = 480, height = 480,type = "cairo")
  par(mar=c(6,6,2,1))
  bw1=bw.nrd(y1[!is.na(y1)])
  bw0=bw.nrd(y0[!is.na(y0)])
  tmp1=density(y1[!is.na(y1)],bw=bw1+0.1)
  tmp0=density(y0[!is.na(y0)],bw=bw0+0.1)
  
  xlim=c(min(tmp1$x,tmp0$x),max(c(tmp1$x,tmp0$x)))
  ylim=c(min(tmp1$y,tmp0$y),max(c(tmp1$y,tmp0$y)))
  plot(tmp1,col=2,lwd=4,xlab="Methylation M value",main=rownames(tcga_normal_methy)[i],bty='l',xlim=xlim,ylim=ylim,cex.lab=1.4,cex.axis=1.4,cex.main=1.4)
  lines(tmp0,col=3,lwd=4)
  legend("topright",c("Normal","Cancer"),lwd=4,col=c(3,2),cex=1.4,bty="n")
  dev.off()
  
  png(paste0(resfolder,Cancer.type,"_",probeid,"_density2",".png"),width = 480, height = 480,type = "cairo")
  par(mar=c(6,6,2,1))
  bw1=bw.nrd(y1[!is.na(y1)])
  bw0=bw.nrd(y0[!is.na(y0)])
  tmp1=density(y1[!is.na(y1)],bw=bw1+0.2)
  tmp0=density(y0[!is.na(y0)],bw=bw0+0.2)
  
  xlim=c(min(tmp1$x,tmp0$x),max(c(tmp1$x,tmp0$x)))
  ylim=c(min(tmp1$y,tmp0$y),max(c(tmp1$y,tmp0$y)))
  plot(tmp1,col=2,lwd=4,xlab="Methylation M value",main=rownames(tcga_normal_methy)[i],bty='l',xlim=xlim,ylim=ylim,cex.lab=1.4,cex.axis=1.4,cex.main=1.4)
  lines(tmp0,col=3,lwd=4)
  legend("topright",c("Normal","Cancer"),lwd=4,col=c(3,2),cex=1.4,bty="n")
  dev.off()
  
  png(paste0(resfolder,Cancer.type,"_",probeid,"_density3",".png"),width = 480, height = 480,type = "cairo")
  par(mar=c(6,6,2,1))
  bw1=bw.nrd(y1[!is.na(y1)])
  bw0=bw.nrd(y0[!is.na(y0)])
  tmp1=density(y1[!is.na(y1)],bw=bw1+0.3)
  tmp0=density(y0[!is.na(y0)],bw=bw0+0.3)
  
  xlim=c(min(tmp1$x,tmp0$x),max(c(tmp1$x,tmp0$x)))
  ylim=c(min(tmp1$y,tmp0$y),max(c(tmp1$y,tmp0$y)))
  plot(tmp1,col=2,lwd=4,xlab="Methylation M value",main=rownames(tcga_normal_methy)[i],bty='l',xlim=xlim,ylim=ylim,cex.lab=1.4,cex.axis=1.4,cex.main=1.4)
  lines(tmp0,col=3,lwd=4)
  legend("topright",c("Normal","Cancer"),lwd=4,col=c(3,2),cex=1.4,bty="n")
  dev.off()
  
}
#"TCGA-PRAD"
cpgs=c("cg02447380","cg05839377","cg02447380","cg00549910","cg17080504")
plot_figure3_5(probeid = cpgs[1])
plot_figure3_5(probeid = cpgs[2])
plot_figure3_5(probeid = cpgs[3])
plot_figure3_5(probeid = cpgs[4])
plot_figure3_5(probeid = cpgs[5])

#"TCGA-COAD"
cpgs=c("cg20319091","cg26090652","cg00566635")
plot_figure3_5(probeid = cpgs[1])
plot_figure3_5(probeid = cpgs[2])
plot_figure3_5(probeid = cpgs[3],xpos=0.35,ypos=0.5)

#"TCGA-LUAD"
cpgs=c("cg01696226","cg10332700","cg24578679")
plot_figure3_5(probeid = cpgs[1])
plot_figure3_5(probeid = cpgs[2])
plot_figure3_5(probeid = cpgs[3],xpos=0.35,ypos=0.5)
plot_figure3_5(probeid = cpgs[1])
plot_figure3_5(probeid = cpgs[2])

#"TCGA-LUSC"
cpgs=c("cg17147211","cg02679809")
plot_figure3_5(probeid = cpgs[1])
plot_figure3_5(probeid = cpgs[2])

#"TCGA-BRCA"
cpgs=c("cg08358907","cg23351584")
plot_figure3_5(probeid = cpgs[1])
plot_figure3_5(probeid = cpgs[2])
plot_figure3_5(probeid = "cg00030422")

#EAC
cpgs=c("cg11805669","cg07171538","cg10188897","cg06377152")
plot_figure3_5(probeid = cpgs[1])
plot_figure3_5(probeid = cpgs[2],xpos=0.35,ypos=0.5)
plot_figure3_5(probeid = cpgs[3])
plot_figure3_5(probeid = cpgs[4])

#TCGA-LIHC
cpgs=c("cg07345734","cg17641861","cg27487839")
plot_figure3_5(probeid = cpgs[1])
plot_figure3_5(probeid = cpgs[2])
plot_figure3_5(probeid = cpgs[3])


#generate new figures 11/30---

png(paste0(resfolder,Cancer.type,"_DMVCvolcanoplot.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
bdiff <- output[,6]-output[,4]
idx=output[,23]<0.05/nrow(output) & bdiff>0 & output[,4]<0.1 & output[,7]>0.05 &output[,23]>0
sum(idx)
idx1=output[,23]>=0.05/nrow(output) & bdiff>0 & output[,4]<0.1& output[,7]>0.05 &output[,23]>0
idx2= bdiff>0 & output[,4]<0.1 & output[,7]>0.05 &output[,23]>0
xmax=round(max(output[idx,6]-output[idx,4])+0.05,1)
ymax=round(max(-log(output[idx,23],base=10))+5,0)
plot(output[idx,6]-output[idx,4],-log(output[idx,23],base=10),xlim=c(0,xmax),col=alpha("red", 0.05),xlab="Difference of mean beta",ylab=expression('-log'[10]*' p-value for DMVC'),
     ylim=c(0,ymax),bty='l',cex.axis=1.4,cex.lab=1.4,cex.main=1.4)
points(output[idx1,6]-output[idx1,4],-log(output[idx1,23],base=10),col=alpha("black", 0.1))
#title("Test for hypermethylation DMVC")
dev.off()
png(paste0(resfolder,Cancer.type,"_DMCvolcanoplot.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
idx=output[,18]<0.05/nrow(output) & bdiff>0 & output[,4]<0.1 & output[,7]>0.05
sum(idx)
idx1=output[,18]>=0.05/nrow(output) & bdiff>0 & output[,4]<0.1& output[,7]>0.05
idx2= bdiff>0 & output[,4]<0.1 & output[,7]>0.05 &output[,23]>0
plot(output[idx,6]-output[idx,4],-log(output[idx,18],base=10),xlim=c(0,xmax),col=alpha("red", 0.05),xlab="Difference of mean beta",ylab=expression('-log'[10]*' p-value for DMC'),
     ylim=c(0,ymax),bty='l',cex.axis=1.4,cex.lab=1.4,cex.main=1.4)
points(output[idx1,6]-output[idx1,4],-log(output[idx1,18],base=10),col=alpha("black", 0.1))
#title("Test for hypermethylation DMC")
dev.off()



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
tmp=rowMeans(alltumor)
idx=which(is.na(tmp))
length(idx) #8404, remove probes with missing values (some may correlated to specific cancer)
alltumor=alltumor[-idx,]
sum(comprobes %in% rownames(alltumor)) #32788
tumortype0=tumortype
tumortype=factor(tumortype,levels = Cancertypes)
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
png(paste0(resfolder,"Allcancers_ALLsig_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(main="All significant probes",prefix="Allcancers_ALLsig_")
dev.off()
png(paste0(resfolder,"Allcancers_ALLsig_PCA_PC1_PC3.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(main="All significant probes",pc1=1,pc2=3)
dev.off()
png(paste0(resfolder,"Allcancers_ALLsig_PCA_PC2_PC3.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(main="All significant probes",pc1=2,pc2=3)
dev.off()
png(paste0(resfolder,"Allcancers_ALLsig_PCA_PC4_PC5.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(main="All significant probes",pc1=4,pc2=5)
dev.off()

# png(paste0(resfolder,"Allcancers_ALLsig_MPCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
# par(mar=c(6,6,2,1))
# plot.pca(main="All significant probes (M value)",prefix="Allcancers_ALLsig_M",opt="M")
# dev.off()
# png(paste0(resfolder,"Allcancers_ALLsig_MPCA_PC1_PC3.png"),width = 480, height = 480,type = "cairo")
# par(mar=c(6,6,2,1))
# plot.pca(main="All significant probes (M value)",pc1=1,pc2=3,opt="M")
# dev.off()
# png(paste0(resfolder,"Allcancers_ALLsig_MPCA_PC2_PC3.png"),width = 480, height = 480,type = "cairo")
# par(mar=c(6,6,2,1))
# plot.pca(main="All significant probes (M value)",pc1=2,pc2=3,opt="M")
# dev.off()
# png(paste0(resfolder,"Allcancers_ALLsig_MPCA_PC4_PC5.png"),width = 480, height = 480,type = "cairo")
# par(mar=c(6,6,2,1))
# plot.pca(main="All significant probes (M value)",pc1=4,pc2=5,opt="M")
# dev.off()
#plot.pca(dat=allnormal,types=normaltype)

png(paste0(resfolder,"Allcancers_ALLuniq_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=alltumor,probes=alluniqprobes,types=tumortype,main="All unique probes",prefix="Allcancers_ALLuniq_")
dev.off()
png(paste0(resfolder,"Allcancers_ALLuniq_PCA_PC1_PC3.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=alltumor,probes=alluniqprobes,types=tumortype,main="All unique probes",pc1=1,pc2=3)
dev.off()
png(paste0(resfolder,"Allcancers_ALLuniq_PCA_PC2_PC3.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=alltumor,probes=alluniqprobes,types=tumortype,main="All unique probes",pc1=2,pc2=3)
dev.off()
png(paste0(resfolder,"Allcancers_ALLuniq_PCA_PC4_PC5.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=alltumor,probes=alluniqprobes,types=tumortype,main="All unique probes",pc1=4,pc2=5)
dev.off()

# png(paste0(resfolder,"Allcancers_ALLuniq_MPCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
# par(mar=c(6,6,2,1))
# plot.pca(dat=alltumor,probes=alluniqprobes,types=tumortype,main="All unique probes (M value)",prefix="Allcancers_ALLuniq_M",opt="M")
# dev.off()
# png(paste0(resfolder,"Allcancers_ALLuniq_MPCA_PC1_PC3.png"),width = 480, height = 480,type = "cairo")
# par(mar=c(6,6,2,1))
# plot.pca(dat=alltumor,probes=alluniqprobes,types=tumortype,main="All unique probes (M value)",pc1=1,pc2=3,opt="M")
# dev.off()
# png(paste0(resfolder,"Allcancers_ALLuniq_MPCA_PC2_PC3.png"),width = 480, height = 480,type = "cairo")
# par(mar=c(6,6,2,1))
# plot.pca(dat=alltumor,probes=alluniqprobes,types=tumortype,main="All unique probes (M value)",pc1=2,pc2=3,opt="M")
# dev.off()
# png(paste0(resfolder,"Allcancers_ALLuniq_MPCA_PC4_PC5.png"),width = 480, height = 480,type = "cairo")
# par(mar=c(6,6,2,1))
# plot.pca(dat=alltumor,probes=alluniqprobes,types=tumortype,main="All unique probes (M value)",pc1=4,pc2=5,opt="M")
# dev.off()

#plot.pca(dat=alltumor,probes=alluniqprobes,types=tumortype,yinc = 1.4,pc1=2,pc2=3)

png(paste0(resfolder,"Allcancers_ALLprobes_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=alltumor,probes=rownames(alltumor),types=tumortype,main="All probes",prefix="Allcancers_ALLprobes_")
dev.off()
png(paste0(resfolder,"Allcancers_ALLprobes_PCA_PC1_PC3.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=alltumor,probes=rownames(alltumor),types=tumortype,main="All probes",pc1=1,pc2=3)
dev.off()
png(paste0(resfolder,"Allcancers_ALLprobes_PCA_PC2_PC3.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(dat=alltumor,probes=rownames(alltumor),types=tumortype,main="All probes",pc1=2,pc2=3)
dev.off()

# png(paste0(resfolder,"Allcancers_ALLprobes_MPCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
# par(mar=c(6,6,2,1))
# plot.pca(dat=alltumor,probes=rownames(alltumor),types=tumortype,main="All probes (M value)",prefix="Allcancers_ALLprobes_M",opt="M")
# dev.off()
# png(paste0(resfolder,"Allcancers_ALLprobes_MPCA_PC1_PC3.png"),width = 480, height = 480,type = "cairo")
# par(mar=c(6,6,2,1))
# plot.pca(dat=alltumor,probes=rownames(alltumor),types=tumortype,main="All probes (M value)",pc1=1,pc2=3,opt="M")
# dev.off()
# png(paste0(resfolder,"Allcancers_ALLprobes_MPCA_PC2_PC3.png"),width = 480, height = 480,type = "cairo")
# par(mar=c(6,6,2,1))
# plot.pca(dat=alltumor,probes=rownames(alltumor),types=tumortype,main="All probes (M value)",pc1=2,pc2=3,opt="M")
# dev.off()

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
png(paste0(resfolder,"Allcancers_ALLsig_TSNE.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.tsne(main="All significant probes")
dev.off()
png(paste0(resfolder,"Allcancers_ALLuniq_TSNE.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.tsne(dat=alltumor,probes=alluniqprobes,types=tumortype,main="All unique probes")
dev.off()
png(paste0(resfolder,"Allcancers_ALLprobes_TSNE.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.tsne(dat=alltumor,probes=rownames(alltumor),types=tumortype,main="All probes")
dev.off()
png(paste0(resfolder,"Allcancers_ALLprobesexptsig_TSNE.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.tsne(dat=alltumor,probes=rownames(alltumor)[! rownames(alltumor) %in% comprobes],types=tumortype,main="All probes except significat ones")
dev.off()

png(paste0(resfolder,"Allcancers_randprobes_TSNE.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
set.seed(12345)
plot.tsne(dat=alltumor,probes=rownames(alltumor)[sample(nrow(alltumor))[1:length(comprobes)]],types=tumortype,main="Random probes")
dev.off()

# png(paste0(resfolder,"Allcancers_ALLsig_MTSNE.png"),width = 480, height = 480,type = "cairo")
# par(mar=c(6,6,2,1))
# plot.tsne(main="All significant probes (M value)",opt="M")
# dev.off()
# png(paste0(resfolder,"Allcancers_ALLuniq_MTSNE.png"),width = 480, height = 480,type = "cairo")
# par(mar=c(6,6,2,1))
# plot.tsne(dat=alltumor,probes=alluniqprobes,types=tumortype,main="All unique probes (M value)",opt="M")
# dev.off()
# png(paste0(resfolder,"Allcancers_ALLprobes_MTSNE.png"),width = 480, height = 480,type = "cairo")
# par(mar=c(6,6,2,1))
# plot.tsne(dat=alltumor,probes=rownames(alltumor),types=tumortype,main="All probes (M value)",opt="M")
# dev.off()
# png(paste0(resfolder,"Allcancers_ALLprobesexptsig_MTSNE.png"),width = 480, height = 480,type = "cairo")
# par(mar=c(6,6,2,1))
# plot.tsne(dat=alltumor,probes=rownames(alltumor)[! rownames(alltumor) %in% comprobes],types=tumortype,main="All probes except significat ones (M value)",opt="M")
# dev.off()

#plot.tsne(dat=allnormal,types=normaltype,main="All significant probes")

#use a glmnet multinomial model to select probes among the pool of significant probes----
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
save(cvfit,file="result/cvfit.RData")
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
length(selectprobes1)
png(paste0(resfolder,"Allcancers_glmprobes_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(main="Glmnet selected probes",probes=selectprobes1)
dev.off()
png(paste0(resfolder,"Allcancers_glmprobes_PCA_PC1_PC3.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(main="Glmnet selected probes",probes=selectprobes1,pc1=1,pc2=3)
dev.off()
png(paste0(resfolder,"Allcancers_glmprobes_PCA_PC2_PC3.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.pca(main="Glmnet selected probes",probes=selectprobes1,pc1=2,pc2=3)
dev.off()
png(paste0(resfolder,"Allcancers_glmprobes_TSNE.png"),width = 480, height = 480,type = "cairo")
par(mar=c(6,6,2,1))
plot.tsne(main="Glmnet selected probes",probes=selectprobes1)
dev.off()

#get confusion matrix
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
# y           TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-PRAD
# TCGA-PRAD         0         0         0         0       498
# TCGA-COAD         0       296         0         0         0
# TCGA-LUAD         0         0       453         5         0
# TCGA-LUSC         0         0        14       356         0
# TCGA-BRCA       782         0         0         0         0
#1se
# resps
# y           TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-PRAD
# TCGA-PRAD         0         0         0         0       498
# TCGA-COAD         0       296         0         0         0
# TCGA-LUAD         0         0       452         6         0
# TCGA-LUSC         0         0        20       350         0
# TCGA-BRCA       782         0         0         0         0


# for (i in 1:length(Cancertypes))
# {
#   idx=which(y==Cancertypes[i])
#   print(paste0(Cancertypes[i],":"))
#   print(quantile(resprob[idx,i]))
# }

#Testing data----------
#Hutch PRAD data
library(data.table)
Hutch_PRAD=fread("/fh/fast/stanford_j/Illumina/Data_Cleaned/Clean_Methylation/beta_main_HP.csv")
Hutch_PRAD=as.data.frame(Hutch_PRAD)
rownames(Hutch_PRAD)=Hutch_PRAD[,1]
Hutch_PRAD=Hutch_PRAD[,2:ncol(Hutch_PRAD)]
sum(!colnames(TCGAtraindat) %in% rownames(Hutch_PRAD)) #23
sum(!selectprobes1 %in% rownames(Hutch_PRAD)) #0
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
# TCGA-PRAD 
# 525 
#1se
# TCGA-PRAD 
# 525 

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
# #try knn cluster based on PCA data
# pcaclust=function(dat1=TCGAtraindat,dat2=Hutch_PRAD_x,dat1_class=as.character(tumortype),selectprobes=selectprobes1)
# {
#   library(class)
#   tmp=sum(is.na(dat2))
#   if (tmp>0) 
#   {
#     warning(paste0(tmp," NA exist"))
#     for (i in 1:ncol(dat2))
#     {
#       idx=which(is.na(dat2[,i]))
#       if (length(idx)>0)
#       {
#         dat2[idx,i]=mean(dat2[,i],na.rm=T)
#       }
#     }
#   }
#   dat=rbind(dat1,dat2)
#   dat=dat[,colnames(dat) %in% selectprobes]
#   pc=prcomp(dat, scale. = T)
#   # plot(pc,type='l')
#   # tmp=pc$sdev^2
#   # plot(tmp[1:10]/sum(tmp),type='l',ylab="Variance%",xlab="")
#   # points(tmp[1:10]/sum(tmp))
#   comp <- data.frame(pc$x[,1:3])
#   knnClust <- knn(train =comp[1:nrow(dat1),], test = comp[(nrow(dat1)+1):nrow(comp),] , k = 5, cl = dat1_class,prob=T)
#   names(knnClust)=rownames(dat2)
#   table(knnClust)
#   return(knnClust)
# }
# Hutch_PRAD_predictcancer_pcaknn=pcaclust()
# table(Hutch_PRAD_predictcancer_pcaknn)
# # Hutch_PRAD_predictcancer_pcaknn
# # TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-PRAD 
# # 0         0         5         1       519 

#Use KNN for prediction
alldat=form_xpredict1()
knnclust=function(dat1=alldat$dat1,dat2=alldat$dat2,dat1_class=as.character(tumortype))
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
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-PRAD 
# 0         0         0         0       525
#1se
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-PRAD 
# 0         0         0         0       525 
quantile(Hutch_PRAD_predictcancer_knn$prob)
# 0%  25%  50%  75% 100% 
# 0.6  1.0  1.0  1.0  1.0 

#Try Hutch COAD
load("/fh/fast/dai_j/CancerGenomics/Colorectal_Cancer/Data/GSE48684.RData")
CancerID <- row.names(clinicaltable)[clinicaltable$status=="cancer"]
AdenoID <- row.names(clinicaltable)[clinicaltable$status=="adenoma"]
#NormalID <- row.names(clinicaltable)[clinicaltable$status=="normal-H"]
Grady_COAD <- MEall[,names(MEall)%in% CancerID]
#Grady_Amethy <- MEall[,names(MEall)%in% AdenoID]
Grady_COAD_x=form_xpredict(dat2=t(Grady_COAD)) #0 probes not available
Grady_COAD_predictprob=predict(cvfit, s=lambda.best,newx = Grady_COAD_x, type = "response")
Grady_COAD_predictcancer=predit_cancertype(allprobs=Grady_COAD_predictprob)
table(Grady_COAD_predictcancer$class)
# TCGA-BRCA TCGA-COAD TCGA-LUAD 
# 2        57         5
#1se
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC 
# 1        58         3         2
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
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-PRAD 
# 1        59         3         1         0 
#1se
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-PRAD 
# 1        57         5         1         0

#Try GEO BRCA
GEO_BRCA=as.data.frame(fread("data/BRCA/GSE75067/GSE75067_series_matrix.txt",skip=58))
rownames(GEO_BRCA)=GEO_BRCA$ID_REF
if ("ID_REF" %in% colnames(GEO_BRCA)) GEO_BRCA=GEO_BRCA[,-1]

GEO_BRCA_x=form_xpredict(dat2=t(GEO_BRCA)) #0 probes not available
GEO_BRCA_predictprob=predict(cvfit, s=lambda.best,newx = GEO_BRCA_x, type = "response")
GEO_BRCA_predictcancer=predit_cancertype(allprobs=GEO_BRCA_predictprob)
table(GEO_BRCA_predictcancer$class,useNA="ifany")
# TCGA-BRCA  
#188
#1se
# TCGA-BRCA  
#188

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
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-PRAD 
# 188         0         0         0         0 
#1se
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-PRAD 
# 188         0         0         0         0 

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

GEO_LUAD_x=form_xpredict(dat2=t(GEO_LUAD)) #0 probes not available
GEO_LUAD_predictprob=predict(cvfit, s=lambda.best,newx = GEO_LUAD_x, type = "response")
GEO_LUAD_predictcancer=predit_cancertype(allprobs=GEO_LUAD_predictprob)
table(GEO_LUAD_predictcancer$class)
# TCGA-COAD TCGA-LUAD TCGA-LUSC 
# 1        79         3 
#1se
# TCGA-COAD TCGA-LUAD TCGA-LUSC 
# 1        79         3 

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
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-PRAD 
# 0         0        80         3         0 
#1se
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-PRAD 
# 0         0        80         3         0 

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
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-PRAD 
# 0         0         3        20         0
#1se
# TCGA-BRCA TCGA-COAD TCGA-LUAD TCGA-LUSC TCGA-PRAD 
# 0         0         3        20         0 

#do heatmap on TCGA data
library(ComplexHeatmap)
library(GetoptLong)
library(circlize)

uniq_cancertype=Cancertypes
listcolor_cancertype=rep("NA",length(uniq_cancertype))
names(listcolor_cancertype)=uniq_cancertype
for (i in 1:length(uniq_cancertype))
{
  listcolor_cancertype[i]=allcolors[i]
}
listcolor_cancertype=list(Type=listcolor_cancertype)

#method: It can be a pre-defined character which is in ("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall").
drawheatmap=function(mat=t(TCGAtraindat[,colnames(TCGAtraindat) %in% selectprobes1]), cluster=tumortype,
                     showrowname=F,cluster_col=F,includeha=T,method="euclidean")
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
  
  if (includeha==T)
  {
    ht_list = Heatmap(mat, name = "Beta value",
                      col = colorRamp2(c(0,0.5,1), c("white", "pink", "red")), 
                      #column_dend_height = unit(4, "cm"),
                      cluster_rows = T,
                      clustering_distance_rows=method,
                      clustering_distance_columns=method,
                      cluster_columns=cluster_col,
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
  }else
  {
    ht_list = Heatmap(mat, name = "Beta value",
                      col = colorRamp2(c(0,0.5,1), c("white", "pink", "red")), 
                      #column_dend_height = unit(4, "cm"),
                      cluster_rows = T,
                      clustering_distance_rows=method,
                      clustering_distance_columns=method,
                      cluster_columns=cluster_col,
                      column_dend_reorder=F,
                      #clustering_distance_columns="canberra",#kendall",
                      row_names_gp=gpar(fontsize = 8, fontface = "bold"),
                      column_dend_height=unit(2,"cm"),
                      row_dend_width=unit(2,"cm"),
                      #top_annotation = c(ha),
                      show_row_dend=T,
                      show_heatmap_legend=T,
                      heatmap_height=unit(12,"cm"),
                      show_column_names = FALSE, show_row_names=showrowname) 
  }
  
  ht_list = draw(ht_list, heatmap_legend_side = "right")
  return(ht_list)
}
png(paste0(resfolder,"Allcancers_Heatmap_cvfit_nostand.png"),width = 800, height = 480,type = "cairo")
drawheatmap()
dev.off()
#draw on each cancer
method="spearman"
method="pearson"
for (i in 1:length(Cancertypes))
{
  probes=sigprobes[[i]]
  mat=t(TCGAtraindat[tumortype==Cancertypes[i],colnames(TCGAtraindat) %in% probes])
  #mat=mat[1:100,]
  png(paste0(resfolder,Cancertypes[i],"_",method,"_Heatmap_sigprobe.png"),width = 800, height = 480,type = "cairo")
  drawheatmap(mat=mat, cluster=as.character(tumortype[tumortype==Cancertypes[i]]),
                       showrowname=F,cluster_col=T,includeha=F,method=method)
  dev.off()
}