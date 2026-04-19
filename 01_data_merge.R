rm(list = ls())

library(limma)
library(sva)
library(tidyverse)
library(limma)
library(ggrepel)
library(ggthemes)
library(tidyverse)
library(pheatmap)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 



express <- read.csv("GSE50796matrix.csv", row.names = 1)
group_list <- read.table("GSE50796sample.txt", sep = "\t", header = T)
express<-dplyr::select(express,group_list$Accession)

# express$ID=NA
# for (i in 1:nrow(express)) {
#   express$ID[i]<-unlist(str_split_fixed(rownames(express)[i],"//",3))[2]
# }
# express$ID<-str_replace_all(express$ID," ","")
# express<-distinct(express,ID,.keep_all = T)
# rownames(express)<-express$ID
# express<-dplyr::select(express,-ID)



express<-express[which(unlist(str_split_fixed(rownames(express),"///",2))[,2]==""),]
express=na.omit(express)
range(express)
# if yes no need to log transfer, if above this range, have to do log transfer.
# express <- log2(express+1)
# range(express)
express=rbind(geneNames=colnames(express), express)
write.table(express, file="GSE50796.txt", sep="\t", quote=F, col.names=F)

rm(list = ls())










files=c("GSE27976.txt", "GSE50796.txt")       


geneList=list()
for(i in 1:length(files)){
  inputFile=files[i]
  rt=read.table(inputFile, header=T, sep="\t",check.names=F)
  header=unlist(strsplit(inputFile, "\\.|\\-"))
  geneList[[header[1]]]=as.vector(rt[,1])
}
intersectGenes=Reduce(intersect, geneList)


allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
  inputFile=files[i]
  header=unlist(strsplit(inputFile, "\\.|\\-"))
 
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  rt=avereps(data)
  colnames(rt)=paste0(header[1], "_", colnames(rt))

  
  qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
  if(LogC){
    rt[rt<0]=0
    rt=log2(rt+1)}
  if(header[1] != "TCGA"){
    rt=normalizeBetweenArrays(rt)
  }
  
  if(i==1){
    allTab=rt[intersectGenes,]
  }else{
    allTab=cbind(allTab, rt[intersectGenes,])
  }
  batchType=c(batchType, rep(i,ncol(rt)))
}


outTab=ComBat(allTab, batchType, par.prior=TRUE)

library(FactoMineR)
library(factoextra)
library(pca3d) 
pca <- prcomp(t(allTab),scale. = TRUE) 
pheno<-data.frame(ID=colnames(allTab))
pheno$cancer<-pheno$ID
pheno[1:249,2]<-"GSE279767"
pheno[250:263,2]<-"GSE50796"

gr <- factor(pheno$cancer) 

pca3d(pca, group = gr, show.centroids = T, legend = 'topleft') 
snapshotPCA3d(file="before_batch_effect_correction.png")
library(FactoMineR)
library(factoextra)
pca <- prcomp(t(outTab),scale. = TRUE) 
pca3d(pca, group = gr, show.centroids = T, legend = 'topleft') 
snapshotPCA3d(file="after_batch_effect_correction.png")



clin=read.table("clinicaldata.txt",header = T)
clin=dplyr::filter(clin,condition == "Craniosynostosis")
clin=clin$Accession
colnames(outTab)=str_replace_all(colnames(outTab),"GSE27976_","")
colnames(outTab)=str_replace_all(colnames(outTab),"GSE50796_","")

outTab1=outTab[,clin]
outTab1=as.data.frame(outTab1)
outTab=as.data.frame(outTab)

save(outTab,file = "merge_all.RDATA")
outTab=rbind(geneNames=colnames(outTab), outTab)
write.table(outTab, file="merge_all.txt", sep="\t", quote=F, col.names=F)

save(outTab1,file = "merge_Craniosynostosis.RDATA")
outTab1=rbind(geneNames=colnames(outTab1), outTab1)
write.table(outTab1, file="merge_Craniosynostosis.txt", sep="\t", quote=F, col.names=F)

