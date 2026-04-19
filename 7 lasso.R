library(glmnet)
library(tidyverse)
rm(list = ls())
set.seed(1234)

data=read.csv("merge_matrix_new.csv",header = T,sep = ",",row.names = 1)
gene=read.table("machine_learning_intersection.txt",header = T,sep = "\t")
gene=gene[,1]
data=data[gene,]
data=as.data.frame(t(data))
data$ID=rownames(data)
clin=read.table("clinicaldata.txt",header = T,sep = "\t")
colnames(clin)[1]="ID"
aaa=intersect(clin$ID,data$ID)
clin=dplyr::filter(clin, ID %in% aaa)
data=dplyr::filter(data, ID %in% aaa)
data1=inner_join(clin,data,by = "ID")


rt=dplyr::select(data1,gene)
rownames(rt)=data1$ID



x=as.matrix(rt)
y=data1$condition

fit=glmnet(x, y, family = "binomial", alpha=1)
pdf(file="fit.pdf",width=6,height=5.5)
plot(fit)
dev.off()

cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)
pdf(file="cvfit.pdf",width=6,height=5.5)
plot(cvfit)
dev.off()


coef=coef(fit, s = cvfit$lambda.min)
coef



active.min = which(coef != 0)-1
active.min = active.min[-1]
geneids <- colnames(x)[active.min]
geneids



index.min = coef[active.min+1]
index.min

combine <- cbind(geneids, index.min)
colnames(combine)=c("ID","risk")
combine=as.data.frame(combine)
write.csv(combine,"risk.csv",row.names = F)
