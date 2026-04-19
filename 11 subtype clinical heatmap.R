
rm(list = ls())
library(pheatmap)         
expFile="geneexpr.txt"         
clusterFile="Cluster.txt"     
cliFile="clinicaldata.txt"           



exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
exp=t(exp)
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)


sameSample=intersect(row.names(exp), row.names(cluster))
exp=exp[sameSample, , drop=F]
cluster=cluster[sameSample, , drop=F]
expCluster=cbind(exp, cluster)

Project=gsub("(.*?)\\_.*", "\\1", rownames(expCluster))
library(tidyverse)
expCluster=cbind(expCluster, Project)


cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(expCluster), row.names(cli))
expCluster=expCluster[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
data=cbind(expCluster, cli)


data=data[order(data$cluster),]
Type=data[,((ncol(exp)+1):ncol(data))]
Type=Type[,-2]
data=t(data[,1:ncol(exp)])


bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
CluCol=bioCol[1:length(levels(factor(Type$cluster)))]
names(CluCol)=levels(factor(Type$cluster))
ann_colors[["cluster"]]=CluCol


pdf("heatmap.pdf", height=2, width=6)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         cluster_rows =F,
         scale="row",
         show_colnames=F,
         fontsize=6,
         fontsize_row=6,
         fontsize_col=6)
dev.off()










