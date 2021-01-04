
#do figures 1 and 2------------------------------------

rm(list=ls())

library(pROC)
library(scales) #alpha
#do each cancer one by one
Cancer.type="PRAD"
Cancer.type="COAD"
Cancer.type="LUAD"
Cancer.type="LUSC"
Cancer.type="BRCA"
#Cancer.type="Krause-EAC"
Cancer.type="LIHC"
Cancertypes=c("PRAD","COAD","LUAD","LUSC","BRCA","LIHC")
resfolder="/fh/fast/dai_j/CancerGenomics/prostate_methylation/result/"
plotfolder="/fh/fast/dai_j/CancerGenomics/prostate_methylation/result/manuscript/"

for (Cancer.type in Cancertypes)
{
  
  rm(tcga_normal_methy,tcga_tumor_methy,tcga_tumor_methy_all,output)
  #result
  Res.file=paste0(resfolder,"TCGA-",Cancer.type,"_output.RData")
  load(Res.file)
  if (Cancer.type=="PRAD")
  {
    #methylation
    load("/fh/fast/dai_j/CancerGenomics/prostate_methylation/tcga_methy.RData")
    #load("tcga_clinical.RData")
  }else
  {
    RData.file=paste("data/TCGA/","TCGA-",Cancer.type, "_", "methy.RData", sep='')
    load(RData.file)
  }
  
  if (!exists("tcga_normal_methy"))
  {
    tcga_normal_methy=normal_methy
    tcga_tumor_methy=tumor_methy
    tcga_tumor_methy_all=tumor_methy_all
    rm(normal_methy,tumor_methy,tumor_methy_all)
  }
  
  if (any(rownames(output)!=rownames(tcga_normal_methy)))
  {
    stop("something is wrong!")
  }
  
  #figure1--------
  plot_figure1_2=F
  if (plot_figure1_2==T)
  {
    #flter probes sd (c(cancer,normal))>0.03 for figures 1 and 2
    alldat=cbind.data.frame(tcga_normal_methy,tcga_tumor_methy_all)
    library(matrixStats)
    alldat.rowsd=rowSds(as.matrix(alldat),na.rm=T)
    quantile(alldat.rowsd)
    #Figure 1,DVC volcanoplt,manually adjust ylim
    output1=output[alldat.rowsd>0.05,]
    idx=output1[,19]<0.05/nrow(output)
    mycolors=rep("black",nrow(output1))
    mycolors[idx]="red"
    # table(mycolors)
    plot(output1[,3]-output1[,1],-log(output1[,19],base=10),xlab="Difference of SD (beta)",ylab=expression('-log'[10]*' p-value for DVC'),
         col=alpha(mycolors, 0.05),bty='l',cex.axis=1.4,cex.lab=1.4)
    if (Cancer.type=="PRAD") ymax=30
    if (Cancer.type=="BRCA") ymax=50
    if (Cancer.type=="COAD") ymax=30
    if (Cancer.type=="LIHC") ymax=40
    if (Cancer.type=="LUAD") ymax=25
    if (Cancer.type=="LUSC") ymax=30
    png(paste0(plotfolder,"TCGA-",Cancer.type,"_volcanoplot.png"),width = 480, height = 480,type = "cairo")
    par(mar=c(6,6,2,1))
    #p-value came from DVC (fit1)
    plot(output1[,3]-output1[,1],-log(output1[,19],base=10),xlab="Difference of SD (beta)",ylab=expression('-log'[10]*' p-value for DVC'),
         col=alpha(mycolors, 0.05),bty='l',cex.axis=1.4,cex.lab=1.4,ylim=c(0,ymax))
    tmp=par("usr")
    text(0.5*(tmp[1]+tmp[2]),0.92*(tmp[4]-tmp[3])+tmp[3],Cancer.type,cex=1.5)
    dev.off()
    
    #figure2--------
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
    if (Cancer.type=="LUAD")
    {
      png(paste0(plotfolder,"TCGA-",Cancer.type,"_DVC_mean_SD.png"),width = 480, height = 480,type = "cairo")
      par(mar=c(6,6,2,1))
      smoothScatter1(output1[idx,6]-output1[idx,4],output1[idx,3]-output1[idx,1],postPlotHook=NULL,
                     xlab="",ylab="Difference of SD beta",
                     cex.axis=1.4,cex.lab=1.4,useRaster = TRUE,xaxs="i",yaxs="i") 
      tmp=par("usr")
      text(0.5*(tmp[1]+tmp[2]),0.92*(tmp[4]-tmp[3])+tmp[3],Cancer.type,cex=1.5)
      mtext("Difference of mean beta", side=1, line=5,cex=1.4)
      dev.off()
    }else
    {
      png(paste0(plotfolder,"TCGA-",Cancer.type,"_DVC_mean_SD.png"),width = 480, height = 480,type = "cairo")
      par(mar=c(6,6,2,1))
      smoothScatter1(output1[idx,6]-output1[idx,4],output1[idx,3]-output1[idx,1],postPlotHook=NULL,
                     xlab="Difference of mean beta",ylab="Difference of SD beta",
                     cex.axis=1.4,cex.lab=1.4,useRaster = TRUE,xaxs="i",yaxs="i") 
      tmp=par("usr")
      text(0.5*(tmp[1]+tmp[2]),0.92*(tmp[4]-tmp[3])+tmp[3],Cancer.type,cex=1.5)
      dev.off()
    }
    if (Cancer.type=="LIHC") #ylim
    {
      png(paste0(plotfolder,"TCGA-",Cancer.type,"_DVC_mean_SD.png"),width = 480, height = 480,type = "cairo")
      par(mar=c(6,6,2,1))
      smoothScatter1(output1[idx,6]-output1[idx,4],output1[idx,3]-output1[idx,1],postPlotHook=NULL,
                     xlab="Difference of mean beta",ylab="Difference of SD beta",
                     cex.axis=1.4,cex.lab=1.4,useRaster = TRUE,xaxs="i",yaxs="i",ylim=c(0,0.45)) 
      tmp=par("usr")
      text(0.5*(tmp[1]+tmp[2]),0.92*(tmp[4]-tmp[3])+tmp[3],Cancer.type,cex=1.5)
      dev.off()
    }
    
  }
  
  #figure3------------------------
  plotfiure3=F
  if (plotfiure3==T)
  {
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
      
      plot(c(FP),c(TP),type="l",bty='l',xlim=c(0,1),ylim=c(0,1),lwd=3,xlab="1-specificity",ylab="Sensitivity",main=main,cex.lab=1.4,cex.axis=1.4,cex.main=1.5)
      abline(0,1,lty=2)
      text(xpos,ypos,paste0("AUC=",round(auc(roc(ans~predict)),3)),cex=1.4)
    }
    plot_figure3=function(probeid="cg23351584",xpos=0.45,ypos=0.6)
    {
      i=which(rownames(tcga_normal_methy)==probeid)
      png(paste0(plotfolder,"TCGA-",Cancer.type,"_",probeid,"_boxplot_DVC",".png"),width = 480, height = 480,type = "cairo")
      par(mar=c(6,6,2,1))
      y <- c(as.numeric(tcga_normal_methy[i,]),as.numeric(tcga_tumor_methy_all[i,]))
      group <- c(rep(0,length(as.numeric(tcga_normal_methy[i,]))),rep(1,length(as.numeric(tcga_tumor_methy_all[i,]))))
      boxplot(y~group,ylim=c(0,max(y,na.rm=T)),boxlty = 0,whisklty = 0,staplelty = 0,xlab="",ylab="Methylation beta",names=rep("",2),
              outline=F,col=NULL,main=paste0(Cancer.type," ",rownames(tcga_normal_methy)[i]),frame.plot = FALSE,cex.axis=1.4,cex.lab=1.4,cex.main=1.5)
      axis(side = 1,at=c(0,1,2,3),labels=c("","Normal","Cancer",""),cex.axis=1.4,cex.lab=1.4)
      stripchart(y~group,vertical = TRUE, method = "jitter",jitter=0.3,pch = 1,col=c(3,2),add = TRUE)
      tmp=par("usr")
      idx=match(rownames(tcga_normal_methy)[i],rownames(output))
      tmp1=formatC(output[idx,18], format = "e", digits = 2)
      tmp2=formatC(output[idx,19], format = "e", digits = 2)
      text((tmp[1]+tmp[2])*0.4,0.92*max(y,na.rm=T),paste0("p(DMC)=",tmp1),cex=1.4)
      text((tmp[1]+tmp[2])*0.4,0.85*max(y,na.rm=T),paste0("p(DVC)=",tmp2),cex=1.4)
      dev.off()
      
      png(paste0(plotfolder,"TCGA-",Cancer.type,"_",probeid,"_AUC_DVC",".png"),width = 480, height = 480,type = "cairo")
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
      png(paste0(plotfolder,"TCGA-",Cancer.type,"_",probeid,"_density1",".png"),width = 480, height = 480,type = "cairo")
      par(mar=c(6,6,2,1))
      bw1=bw.nrd(y1[!is.na(y1)])
      bw0=bw.nrd(y0[!is.na(y0)])
      tmp1=density(y1[!is.na(y1)],bw=bw1+0.1)
      tmp0=density(y0[!is.na(y0)],bw=bw0+0.1)
      
      xlim=c(min(tmp1$x,tmp0$x),max(c(tmp1$x,tmp0$x)))
      ylim=c(min(tmp1$y,tmp0$y),max(c(tmp1$y,tmp0$y)))
      plot(tmp1,col=2,lwd=4,xlab="Methylation M value",main=rownames(tcga_normal_methy)[i],bty='l',xlim=xlim,ylim=ylim,cex.lab=1.4,cex.axis=1.4,cex.main=1.5)
      lines(tmp0,col=3,lwd=4)
      legend("topright",c("Normal","Cancer"),lwd=4,col=c(3,2),cex=1.4,bty="n")
      dev.off()
      
    }
    #PRAD
    if (Cancer.type=="PRAD")
    {
      cpg="cg17080504"
      plot_figure3(probeid = cpg)
    }
    if (Cancer.type=="LUSC")
    {
      cpg="cg02679809"
      plot_figure3(probeid = cpg)
    }
    if (Cancer.type=="LUAD")
    {
      cpg="cg24578679"
      plot_figure3(probeid = cpg,xpos=0.35,ypos=0.5)
    }
    if (Cancer.type=="COAD")
    {
      cpg="cg00566635"
      plot_figure3(probeid = cpg)
    }
  }
  
  
  #figure5----------------
  plotfigure5=F
  if (plotfigure5==T)
  {
    png(paste0(plotfolder,"TCGA-",Cancer.type,"_DMVCvolcanoplot.png"),width = 480, height = 480,type = "cairo")
    par(mar=c(6,6,2,1))
    bdiff <- output[,6]-output[,4]
    idx=output[,23]<0.05/nrow(output) & bdiff>0 & output[,4]<0.1 & output[,7]>0.05 &output[,23]>0
    sum(idx)
    idx1=output[,23]>=0.05/nrow(output) & bdiff>0 & output[,4]<0.1& output[,7]>0.05 &output[,23]>0
    xmax=round(max(output[idx,6]-output[idx,4])+0.05,1)
    #ymax=round(max(-log(output[idx,23],base=10))+5,0)
    ymax=200
    plot(output[idx,6]-output[idx,4],-log(output[idx,23],base=10),xlim=c(0,xmax),col=alpha("red", 0.05),xlab="Difference of mean beta",ylab=expression('-log'[10]*' p-value for DMVC'),
         ylim=c(0,ymax),bty='l',cex.axis=1.4,cex.lab=1.4,cex.main=1.4)
    points(output[idx1,6]-output[idx1,4],-log(output[idx1,23],base=10),col=alpha("black", 0.1))
    tmp=par("usr")
    text(tmp[1]+(tmp[2]-tmp[1])*0.25,tmp[3]+(tmp[4]-tmp[3])*0.92,Cancer.type,cex=1.5)
    text(tmp[1]+(tmp[2]-tmp[1])*0.25,tmp[3]+(tmp[4]-tmp[3])*0.8,paste0(sum(idx)," DMVC"),cex=1.5)
    #title("Test for hypermethylation DMVC")
    dev.off()
    png(paste0(plotfolder,"TCGA-",Cancer.type,"_DMCvolcanoplot.png"),width = 480, height = 480,type = "cairo")
    par(mar=c(6,6,2,1))
    idx=output[,18]<0.05/nrow(output) & bdiff>0 & output[,4]<0.1 & output[,7]>0.05
    sum(idx)
    idx1=output[,18]>=0.05/nrow(output) & bdiff>0 & output[,4]<0.1& output[,7]>0.05
    plot(output[idx,6]-output[idx,4],-log(output[idx,18],base=10),xlim=c(0,xmax),col=alpha("red", 0.05),xlab="Difference of mean beta",ylab=expression('-log'[10]*' p-value for DMC'),
         ylim=c(0,ymax),bty='l',cex.axis=1.4,cex.lab=1.4,cex.main=1.4)
    points(output[idx1,6]-output[idx1,4],-log(output[idx1,18],base=10),col=alpha("black", 0.1))
    tmp=par("usr")
    text(tmp[1]+(tmp[2]-tmp[1])*0.25,tmp[3]+(tmp[4]-tmp[3])*0.92,Cancer.type,cex=1.5)
    text(tmp[1]+(tmp[2]-tmp[1])*0.25,tmp[3]+(tmp[4]-tmp[3])*0.8,paste0(sum(idx)," DMC"),cex=1.5)
    #title("Test for hypermethylation DMC")
    dev.off()
  }
  
}

#figure6--------------------
load("/fh/fast/dai_j/CancerGenomics/prostate_methylation/result/cancer6_classification_data.RData") #sigprobes,comprobes,alltumor,tumortype,allnormal,normaltype

library(ComplexHeatmap)
library(GetoptLong)
library(circlize)
plotfigure6=T
if (plotfigure6==T)
{
  
  drawheatmap6=function(mat,method="euclidean",Cancer.type="PRAD")
  {
    library(dendsort)
    row_dend = dendsort(hclust(dist(mat,method="euclidean")))
    col_dend = dendsort(hclust(dist(t(mat),method="euclidean")))
    
    ht_list = Heatmap(mat, name = "Beta value",
                      col = colorRamp2(c(0,0.5,1), c("white", "pink", "red")), 
                      #column_dend_height = unit(4, "cm"),
                      cluster_rows = row_dend,
                      row_dend_reorder=T,
                      #clustering_distance_rows=method,
                      #clustering_distance_columns=method,
                      cluster_columns=col_dend,
                      column_dend_reorder=T,
                      #clustering_distance_columns="canberra",#kendall",
                      row_names_gp=gpar(fontsize = 3, fontface = "bold"),
                      column_names_gp=gpar(fontsize = 4, fontface = "bold"),
                      column_dend_height=unit(3,"cm"),
                      #row_dend_width=unit(2,"cm"),
                      #top_annotation = c(ha),
                      show_row_dend=F,
                      show_column_dend=T,
                      show_heatmap_legend=T,
                      #heatmap_height=unit(12,"cm"),
                      column_title=Cancer.type,
                      column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                      show_column_names = F, show_row_names=F) 
    
    
    ht_list = draw(ht_list, heatmap_legend_side = "right")
    return(ht_list)
  }
  sigprobes=NULL
  for (i in 1:length(Cancertypes))
  {
    Cancer.type=Cancertypes[i]
    Res.file=paste0(resfolder,"TCGA-",Cancer.type,"_output.RData")
    load(Res.file)
    bdiff <- output[,6]-output[,4]
    tmp=rownames(output)[output[,23]<0.05/nrow(output) & bdiff>0.05 & output[,4]<0.1 & output[,7]>0.05]
    tmp=tmp[grepl("^cg",tmp)] #remove non-CpG targeting probes ch.xx.xx
    sigprobes=c(sigprobes,list(probes=tmp))
    names(sigprobes)[length(sigprobes)]=Cancer.type
  }
  
  
  for (i in 1:length(Cancertypes))
  #for (i in 1:1)
  {
    cat(Cancertypes[i],'..')
    probes=sigprobes[[i]]
    Cancer.type=Cancertypes[i]
    rm(tcga_normal_methy,tcga_tumor_methy,tcga_tumor_methy_all)
    if (Cancer.type=="PRAD")
    {
      load("/fh/fast/dai_j/CancerGenomics/prostate_methylation/tcga_methy.RData")
      #load("tcga_clinical.RData")
    }else
    {
      RData.file=paste("/fh/fast/dai_j/CancerGenomics/prostate_methylation/data/TCGA/","TCGA-",Cancer.type, "_", "methy.RData", sep='')
      load(RData.file)
    }
    
    if (!exists("tcga_normal_methy"))
    {
      tcga_normal_methy=normal_methy
      tcga_tumor_methy=tumor_methy
      tcga_tumor_methy_all=tumor_methy_all
    }
    
    #mat=t(TCGAtraindat[tumortype==paste0("TCGA-",Cancertypes[i]),colnames(TCGAtraindat) %in% probes])
    mat=tcga_tumor_methy_all[rownames(tcga_tumor_methy_all) %in% sigprobes[[i]],]
    mat= as.matrix(mat)
    #mat=mat[1:1000,]
    png(paste0(plotfolder,"TCGA-",Cancertypes[i],"_Heatmap_sigprobe.png"),width = 800, height = 480,type = "cairo")
    myheatmap=drawheatmap6(mat=mat,method="euclidean",Cancer.type = Cancertypes[i])
    dev.off()
    #png(paste0(plotfolder,"TCGA-",Cancertypes[i],"_spearman_Heatmap_sigprobe.png"),width = 800, height = 480,type = "cairo")
    #myheatmap=drawheatmap6(mat=mat,method="spearman",Cancer.type = Cancertypes[i])
    #dev.off()
    #png(paste0(plotfolder,"TCGA-",Cancertypes[i],"_pearson_Heatmap_sigprobe.png"),width = 800, height = 480,type = "cairo")
    #myheatmap=drawheatmap6(mat=mat,method="pearson",Cancer.type = Cancertypes[i])
    #dev.off()
  }
  
  # colfunc <- colorRampPalette(c("white","pink","red"))
  # heatmap.2(as.matrix(mat),scale="none",col=colfunc(15),trace="none",dendrogram="column",labRow=NA)
  #compare drawheatmap6 and heatmap
  hclust.row <- hclust(get_dist(as.matrix(mat[1:100,]), method = "euclidean", stand = F))
  hclust.col <- hclust(get_dist(t(as.matrix(mat[1:100,])), method = "euclidean", stand = F))
  hclust.col1 <- hclust(dist(t(as.matrix(mat[1:100,])), method = "euclidean"))
  heatmap0=heatmap(as.matrix(mat[1:100,]),scale="none",Rowv=T,Colv=T)
  heatmap1=heatmap(as.matrix(mat[1:100,]),Rowv=as.dendrogram(hclust.row),Colv=as.dendrogram(hclust.col),scale="none")
  all(heatmap1$colInd==heatmap0$colInd) #T
  heatmap2=drawheatmap6(mat=mat[1:100,]) #reorder set to F #https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#reorder-dendrograms
  all(heatmap1$colInd==column_order(heatmap2)) #T
}


plotfigure7=T
load("/fh/fast/dai_j/CancerGenomics/prostate_methylation/result/cancer6_classification_data.RData") #sigprobes,comprobes,alltumor,tumortype,allnormal,normaltype
if (plotfigure7==T)
{
  #allcolors=c("red","blue","darkorchid1","limegreen","goldenrod1","black","brown4","darkseagreen","darkseagreen1")
  #put black first for normal
  allcolors=c("black","red","blue","darkorchid1","limegreen","goldenrod1","brown4","darkseagreen","darkseagreen1")
  
  plot.pca.norm=function(dat=allTCGAdat,probes=comprobes,types=alltype,yinc=1.2,pc1=1,pc2=2,main="",prefix=NULL,opt="beta")
  {
    #types=factor(types,levels = unique(types))
    tmp=unique(types)
    tmp=c("NORM",tmp[tmp!="NORM"])
    types=factor(types,levels = tmp)
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
    # if (!is.null(prefix))
    # {
    #   png(paste0(resfolder,prefix,"PCs_variances.png"),width = 480, height = 480,type = "cairo")
    #   par(mar=c(6,6,2,1))
    #   barplot(expl_var[1:50], ylab="EXPLAINED VARIANCE",
    #           names.arg=paste0("pcadat",seq(1:50)), col="darkgreen",las=2)
    #   dev.off()
    # }
    pcadat=pcadat$x
    xmin=min(pcadat[,pc1])
    xmax=max(pcadat[,pc1])
    ymin=min(pcadat[,pc2])
    ymax=max(pcadat[,pc2])*yinc
    pchs=c(16,0,2,3,4,5,7,8)[1:length(unique(types))]
    pch=pchs[types]
    plot(pcadat[,pc1],pcadat[,pc2],col=alpha(allcolors[types],0.8),pch=pch,cex.axis=1.5,cex.lab=1.5,xlim=c(xmin-0.05*(xmax-xmin),xmax+0.05*(xmax-xmin)),
         ylim=c(ymin-0.05*(ymax-ymin),ymax+0.2*(ymax-ymin)),xlab=paste0("PC",pc1," (",round(expl_var[pc1]*100),"%)"),ylab=paste0("PC",pc2," (",round(expl_var[pc2]*100),"%)"),main=main,bty='l')
    #text(x=pcadat[,1],y=pcadat[,2],rownames(pcadat),col=colors[clinicaltable$`Platium response`],cex=1)
    legend("topleft",legend=levels(types),col=allcolors[1:length(unique(types))],pch=pchs,cex=1.5,bty = "n",ncol=4)
  }
  
  plot.tsne.norm=function(dat=allTCGAdat,probes=comprobes,types=alltype,yinc=1.2,pc1=1,pc2=2,main="",opt="beta")
  {
    tmp=unique(types)
    tmp=c("NORM",tmp[tmp!="NORM"])
    types=factor(types,levels = tmp)
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
    pchs=c(16,0,2,3,4,5,7,8)[1:length(unique(types))]
    pch=pchs[types]
    plot(tsnedat[,pc1],tsnedat[,pc2],col=alpha(allcolors[types],0.8),pch=pch,cex.axis=1.5,cex.lab=1.5,xlim=c(xmin-0.05*(xmax-xmin),xmax+0.05*(xmax-xmin)),
         ylim=c(ymin-0.05*(ymax-ymin),ymax+0.2*(ymax-ymin)),xlab=paste0("TSNE",pc1),ylab=paste0("TSNE",pc2),main=main,bty='l')
    #legend("topleft",legend=c(Cancertypes),col=allcolors[1:length(Cancertypes)],pch=1,cex=1,bty = "n",ncol=3)
    #legend("topleft",legend=c(unique(as.character(types))),col=allcolors[1:length(unique(types))],pch=unique(as.numeric(types)),cex=1,bty = "n",ncol=3)
    legend("topleft",legend=levels(types),col=allcolors[1:length(unique(types))],pch=pchs,cex=1.5,bty = "n",ncol=4)
  }
  
  load("result/cvfit_cancer_normal_classification_SVM.RData")
  library(glmnet)
  lambda.best=cvfit$lambda.min #102 probes (all probes),100 probes (sig probes)
  #lambda.best=cvfit$lambda.1se #46 probes, 47 probes
  glmcoeff1=coef(cvfit,s=lambda.best)
  idx=glmcoeff1@i
  #the selected probes
  selectprobes1=glmcoeff1@Dimnames[[1]][idx+1]
  selectprobes1=selectprobes1[selectprobes1!="(Intercept)"]
  length(selectprobes1) #103
  alltype=c(as.character(tumortype),rep("TCGA-NORM",ncol(allnormal)))
  alltype=gsub("TCGA-","",alltype)
  #Cancers and normals
  png(paste0(plotfolder,"All6cancers_glmprobes_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
  par(mar=c(6,6,2,1))
  plot.pca.norm(main="",probes=selectprobes1)
  dev.off()
  png(paste0(plotfolder,"All6cancers_glmprobes_TSNE.png"),width = 480, height = 480,type = "cairo")
  par(mar=c(6,6,2,1))
  plot.tsne.norm(main="",probes=selectprobes1)
  dev.off()
  png(paste0(plotfolder,"All6cancers_glmprobes_PCA_PC1_PC3.png"),width = 480, height = 480,type = "cairo")
  par(mar=c(6,6,2,1))
  plot.pca.norm(main="",probes=selectprobes1,pc2=3)
  dev.off()
  png(paste0(plotfolder,"All6cancers_glmprobes_PCA_PC1_PC4.png"),width = 480, height = 480,type = "cairo")
  par(mar=c(6,6,2,1))
  plot.pca.norm(main="",probes=selectprobes1,pc2=4)
  dev.off()
  
  ###
  load("result/cvfit_cancer_TOO_SVM_classification.RData")
  lambda.best=cvfit$lambda.min 
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
  
  #use the same pattern to draw
  allcolors1=allcolors[2:length(allcolors)]
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
    pchs=c(0,2,3,4,5,7,8)[1:length(unique(types))]
    pch=pchs[types]
    plot(pcadat[,pc1],pcadat[,pc2],col=alpha(allcolors1[types],0.8),pch=pch,cex.axis=1.5,cex.lab=1.5,xlim=c(xmin-0.05*(xmax-xmin),xmax+0.05*(xmax-xmin)),
         ylim=c(ymin-0.05*(ymax-ymin),ymax+0.2*(ymax-ymin)),xlab=paste0("PC",pc1," (",round(expl_var[pc1]*100),"%)"),ylab=paste0("PC",pc2," (",round(expl_var[pc2]*100),"%)"),main=main,bty='l')
    #text(x=pcadat[,1],y=pcadat[,2],rownames(pcadat),col=colors[clinicaltable$`Platium response`],cex=1)
    #legend("topleft",legend=c(unique(as.character(types))),col=allcolors[1:length(unique(types))],pch=unique(as.numeric(types)),cex=1.5,bty = "n",ncol=4)
    legend("topleft",legend=levels(types),col=allcolors1[1:length(unique(types))],pch=pchs,cex=1.5,bty = "n",ncol=4)
    
  }
  
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
    pchs=c(0,2,3,4,5,7,8)[1:length(unique(types))]
    pch=pchs[types]
    plot(tsnedat[,pc1],tsnedat[,pc2],col=alpha(allcolors1[types],0.8),pch=pch,cex.axis=1.5,cex.lab=1.5,xlim=c(xmin-0.05*(xmax-xmin),xmax+0.05*(xmax-xmin)),
         ylim=c(ymin-0.05*(ymax-ymin),ymax+0.2*(ymax-ymin)),xlab=paste0("TSNE",pc1),ylab=paste0("TSNE",pc2),main=main,bty='l')
    #legend("topleft",legend=c(Cancertypes),col=allcolors[1:length(Cancertypes)],pch=1,cex=1,bty = "n",ncol=3)
    #legend("topleft",legend=c(unique(as.character(types))),col=allcolors[1:length(unique(types))],pch=unique(as.numeric(types)),cex=1,bty = "n",ncol=3)
    legend("topleft",legend=levels(types),col=allcolors1[1:length(unique(types))],pch=pchs,cex=1.5,bty = "n",ncol=4)
  }
  
  alltype=c(as.character(tumortype))
  alltype=gsub("TCGA-","",alltype)
  #only cancers
  png(paste0(plotfolder,"Allcancers_glmprobes_PCA_PC1_PC2.png"),width = 480, height = 480,type = "cairo")
  par(mar=c(6,6,2,1))
  plot.pca(dat=alltumor,main="",probes=selectprobes1)
  dev.off()
  png(paste0(plotfolder,"Allcancers_TSNE.png"),width = 480, height = 480,type = "cairo")
  par(mar=c(6,6,2,1))
  plot.tsne(dat=alltumor,main="",probes=selectprobes1)
  dev.off()
  png(paste0(plotfolder,"Allcancers_glmprobes_PCA_PC1_PC3.png"),width = 480, height = 480,type = "cairo")
  par(mar=c(6,6,2,1))
  plot.pca(dat=alltumor,main="",probes=selectprobes1,pc2=3)
  dev.off()
  png(paste0(plotfolder,"Allcancers_glmprobes_PCA_PC1_PC4.png"),width = 480, height = 480,type = "cairo")
  par(mar=c(6,6,2,1))
  plot.pca(dat=alltumor,main="",probes=selectprobes1,pc2=4)
  dev.off()
  png(paste0(plotfolder,"Allcancers_glmprobes_PCA_PC1_PC5.png"),width = 480, height = 480,type = "cairo")
  par(mar=c(6,6,2,1))
  plot.pca(dat=alltumor,main="",probes=selectprobes1,pc2=5)
  dev.off()
  plot.pca(dat=alltumor,main="",probes=selectprobes1,pc1=2,pc2=4)
  
  ###
  load("result/cvfit_cancer_TOO_SVM_classification.RData")
  lambda.best=cvfit$lambda.min 
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
  
  uniq_cancertype=Cancertypes
  listcolor_cancertype=rep("NA",length(uniq_cancertype))
  names(listcolor_cancertype)=uniq_cancertype
  for (i in 1:length(uniq_cancertype))
  {
    listcolor_cancertype[i]=allcolors1[i]
  }
  listcolor_cancertype=list(Cancer=listcolor_cancertype)
  
  drawheatmap=function(mat, cluster=gsub("TCGA-","",tumortype),
                       showrowname=F,cluster_col=F,columnorder=NULL,cluster_row=T,roworder=NULL,method="euclidean")
  {
    
    # mat=t(scale(t(mat)))
    # mat[mat < -2]=-2
    # mat[mat > 2]=2
    ha = HeatmapAnnotation(Cancer = cluster, 
                           col = listcolor_cancertype,
                           show_annotation_name = T,
                           show_legend=T,
                           annotation_name_gp=gpar(fontface = "bold"),
                           # annotation_name_offset = unit(2, "cm"),
                           # annotation_name_rot = c(0, 0),
                           annotation_name_side = "right")
    

      ht_list = Heatmap(mat, name = "Beta value",
                        col = colorRamp2(c(0,0.5,1), c("white", "pink", "red")), 
                        #column_dend_height = unit(4, "cm"),
                        cluster_rows = cluster_row,
                        row_order=roworder,
                        row_dend_reorder=F,
                        clustering_distance_rows=method,
                        clustering_distance_columns=method,
                        cluster_columns=cluster_col,
                        column_order=columnorder,
                        column_dend_reorder=F,
                        #clustering_distance_columns="canberra",#kendall",
                        row_names_gp=gpar(fontsize = 8, fontface = "bold"),
                        column_dend_height=unit(2,"cm"),
                        row_dend_width=unit(2,"cm"),
                        top_annotation = c(ha),
                        show_row_dend=F,
                        show_column_dend=F,
                        show_heatmap_legend=T,
                        heatmap_height=unit(12,"cm"),
                        show_column_names = FALSE, show_row_names=showrowname) 

    
    ht_list = draw(ht_list, heatmap_legend_side = "right")
    return(ht_list)
  }
  
  mat=alltumor[rownames(alltumor) %in% selectprobes1,]
  hclust.row <- hclust(get_dist(as.matrix(mat), method = "euclidean", stand = F))
  hclust.row1=dendsort(hclust.row)
  roworder=hclust.row1$order
  neworder=NULL
  #reorder samples within each cancer type
  for (i in 1:length(Cancertypes))
  {
    idx=which(tumortype==paste0("TCGA-",Cancertypes[i]))
    hclust.col <- hclust(get_dist(t(as.matrix(mat[,idx])), method = "euclidean", stand = F))
    hclust.col1=dendsort(hclust.col)
    neworder=c(neworder,idx[hclust.col1$order])
  }
  png(paste0(plotfolder,"All6cancers_Heatmap_cvfit_nostand_dendsort.png"),width = 800, height = 480,type = "cairo")
  drawheatmap(mat=mat,columnorder=neworder,roworder = roworder,cluster_row=F)
  dev.off()
}
#