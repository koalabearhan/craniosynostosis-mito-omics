rm(list = ls())
library(circlize)
library(ggsci)
library(parallel)
library(tidyverse)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 


genecor.parallel <- function(data,gene,cl){
  cl <- makeCluster(cl)
  y <- as.numeric(data[gene,])
  rownames <- rownames(data)
  dataframes <- do.call(rbind, parLapply(cl=cl,rownames, function(x){
    dd  <- cor.test(as.numeric(data[x,]), y, type="spearman")
    data.frame(Gene_1=gene, Gene_2=x, cor=dd$estimate, p.value=dd$p.value)
  }))
  stopCluster(cl)
  return(dataframes)
}


genecor_circleplot <- function(x){
  Corr <- data.frame(rbind(data.frame(Gene=x[,1], Correlation=x[,3]), 
                           data.frame(Gene=x[,2], Correlation=x[,3])), stringsAsFactors = F)      
  Corr$Index <- seq(1,nrow(Corr),1) 
  Corr <- Corr[order(Corr[,1]),] 
  corrsp <- split(Corr,Corr$Gene)
  corrspe <- lapply(corrsp, function(x){x$Gene_Start<-0
  
  
  if (nrow(x)==1){x$Gene_End<-1}else{
    x$Gene_End<-sum(abs(x$Correlation))} 
  x})
  GeneID <- do.call(rbind,corrspe)
  GeneID <- GeneID[!duplicated(GeneID$Gene),]
  
  
  mycol <- pal_d3("category20c")(20)
  n <- nrow(GeneID)
  GeneID$Color <- mycol[1:n]
  

  Corr[,2] <- abs(Corr[,2]) 
  corrsl <- split(Corr,Corr$Gene)
  aaaaa <- c()
  corrspl <- lapply(corrsl,function(x){nn<-nrow(x)
  for (i in 1:nn){
    aaaaa[1] <- 0
    aaaaa[i+1] <- x$Correlation[i]+aaaaa[i]}
  bbbbb <- data.frame(V4=aaaaa[1:nn],V5=aaaaa[2:(nn+1)])
  bbbbbb <- cbind(x,bbbbb)
  bbbbbb
  })
  Corr <- do.call(rbind,corrspl)
  
  
  Corr <- Corr[order(Corr$Index),]
  
 
  x$start_1 <- Corr$V4[1:(nrow(Corr)/2)]
  x$end_1 <- Corr$V5[1:(nrow(Corr)/2)]
  x$start_2 <- Corr$V4[(nrow(Corr)/2 + 1):nrow(Corr)]
  x$end_2 <- Corr$V5[(nrow(Corr)/2 + 1):nrow(Corr)]
  
  
  color <- data.frame(colorRampPalette(c("#67BE54", "#FFFFFF", "#F82C2B"))(201))
  
  for (i in 1:nrow(x)){
    x[i,8] <- substring(color[x[i,3] * 100 + 101, 1], 1, 7)
  }
  names(x)[8] <- "color"
  
  
  circos.clear()
  circos.par(start.degree = 90, 
             gap.degree = 5, 
             track.margin = c(0,0.23), 
             cell.padding = c(0,0,0,0)
  )
  circos.initialize(factors = GeneID$Gene,
                    xlim = cbind(GeneID$Gene_Start, GeneID$Gene_End))
  
  
  circos.trackPlotRegion(ylim = c(0, 1), factors = GeneID$Gene, 
                         track.height = 0.05, 
                         panel.fun = function(x, y) {
                           name = get.cell.meta.data("sector.index") 
                           i = get.cell.meta.data("sector.numeric.index") 
                           xlim = get.cell.meta.data("xlim")
                           ylim = get.cell.meta.data("ylim")
                           circos.text(x = mean(xlim), y = 1,
                                       labels = name,
                                       cex = 1, 
                                       niceFacing = TRUE, 
                                       facing = "bending", 
                                       adj = c(0.5, -2.8), 
                                       font = 2 
                           
                           circos.rect(xleft = xlim[1], 
                                       ybottom = ylim[1],
                                       xright = xlim[2], 
                                       ytop = ylim[2],
                                       col = GeneID$Color[i],
                                       border = GeneID$Color[i])
                           
                           circos.axis(labels.cex = 0.7, 
                                       direction = "outside"
                           )})
  
  
  for(i in 1:nrow(x)){
    circos.link(sector.index1 = x$Gene_1[i], 
                point1 = c(x[i, 4], x[i, 5]),
                sector.index2 = x$Gene_2[i], 
                point2 = c(x[i, 6], x[i, 7]),
                col = paste(x$color[i], "C9", sep = ""), 
                border = FALSE, 
                rou = 0.7
    )}
  
  
  i <- seq(0,0.995,0.005)
  rect(-1+i/2, #xleft
       -1, #ybottom
       -0.9975+i/2, #xright
       -0.96, #ytop
       col = paste(as.character(color[,1]), "FF", sep = ""),
       border = paste(as.character(color[,1]), "FF", sep = ""))
  text(-0.97, -1.03, "-1")
  text(-0.51, -1.03, "1")
}

inputtemp <- read.csv("merge_matrix_new.csv", header = T,sep = "," ,row.names = 1)
aaa<-read.table("machine_learning_intersection.txt",header = T)
aaa<-aaa[,1]
inputtemp<-as.data.frame(t(inputtemp))
inputtemp<-dplyr::select(inputtemp,aaa)
clin<-read.table("clinicaldata.txt",header = T)
inputtemp=inputtemp[clin$Accession,]

ccc<-inputtemp
ddd<-clin
colnames(ddd)[1]="ID"
ccc$ID<-rownames(ccc)
eee<-inner_join(ddd,ccc,by="ID")
write.table(eee,file = "mergeROC.txt",row.names = F,sep = "\t")

genecorl <- lapply(colnames(inputtemp),function(x){
  ddd <- genecor.parallel(data = t(inputtemp), cl=1, gene=x) 
  ddd  
})
genecor <- do.call(rbind, genecorl)


genecorr <- genecor[-which(genecor$p.value==0),]


genecorrr<-genecorr[!duplicated(genecorr$cor),]


genecorrr$p.value <- NULL
genecor_circleplot(genecorrr)

pdf("mergecorrelation.pdf", width = 5, height = 5)
genecor_circleplot(genecorrr)
dev.off()



