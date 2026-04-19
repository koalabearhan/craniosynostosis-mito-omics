rm(list = ls())
load("merge_Craniosynostosis.RDATA")
Merge<-outTab1
genename<-read.table("machine_learning_intersection.txt",header = T,sep = "\t")
AAA<-intersect(genename[,1],rownames(Merge))
geneexpr<-Merge[AAA,]

write.table(geneexpr,"geneexpr.txt",quote = F,row.names =T ,sep = "\t")
save(geneexpr,file = "geneexpr.RDATA")
library(ConsensusClusterPlus)



data=geneexpr
data=as.matrix(data)


maxK=9
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=50,
                             pItem=0.8,
                             pFeature=1,
                             clusterAlg="pam",
                             distance="euclidean",
                             seed=123456,
                             plot="png")



clusterNum=2      
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("cluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$cluster))
cluster$cluster=letter[match(cluster$cluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="Cluster.txt", sep="\t", quote=F, col.names=F)

