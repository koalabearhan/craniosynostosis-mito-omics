

rm(list = ls())
library(limma)
library(ggrepel)
library(ggthemes)
library(tidyverse)
library(pheatmap)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 



load("merge_all.RDATA")
express <- outTab
group_list <- read.table("clinicaldata.txt", sep = "\t", header = T)
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
#express <- log2(express+1)

cols=rainbow(ncol(express)) 
pdf(file="before_normalization.PDF", width=10, height=6)
boxplot(express,outline = F,col =cols) 
dev.off()


# Normalise data quantiles
library(preprocessCore)
express.norm <- normalize.quantiles(as.matrix(express))

colnames(express.norm) <- colnames(express)
rownames(express.norm) <- rownames(express)

cols=rainbow(ncol(express.norm)) 
pdf(file="after_normalization.PDF", width=10, height=6)
boxplot(express.norm,outline = F,col =cols) 
dev.off()

#express.norm <- express
write.csv(express.norm,file = "merge_matrix_new.csv",row.names = T)



design <- model.matrix(~ 0 + factor(group_list$condition))
colnames(design) <- levels(factor(group_list$condition))
rownames(design) <- colnames(express.norm)
design


contrast.matrix <- makeContrasts(Craniosynostosis-Control, levels = design)


fit <- lmFit(express.norm,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)

# b vs. a
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")


write.csv(x, "differential_expression_results.csv", quote = F)


p.cut<-0.05
logFC.cut<-0


volcano<-x
volcano$type[(volcano$P.Value > p.cut|volcano$P.Value=="NA")|(volcano$logFC < logFC.cut)& volcano$logFC > -logFC.cut] <- "none significant"
volcano$type[volcano$P.Value <= p.cut & volcano$logFC >= logFC.cut] <- "up-regulated"
volcano$type[volcano$P.Value <= p.cut & volcano$logFC <= -logFC.cut] <- "down-regulated"  

p = ggplot(volcano,aes(logFC,-1*log10(P.Value),color=type))       
p + geom_point()

x_lim <- max(volcano$logFC,-x$logFC) 
gg=p + geom_point( aes(size = abs(logFC)),alpha = 0.4) + xlim(-x_lim,x_lim) +labs(x="log2(FC)",y="-log10(P.Value)")+
  scale_color_manual(values =c("blue","grey","red"))+
  geom_hline(aes(yintercept=-1*log10(p.cut)),colour="black", linetype="dashed") +
  geom_vline(xintercept=c(-logFC.cut,logFC.cut),colour="black", linetype="dashed")
print(gg) 


pdf(file="volcano_plot.PDF", width=6, height=6)
gg
dev.off()



volcano1<-dplyr::filter(volcano,P.Value < 0.05)
volcano1<-volcano1[order(volcano1$logFC,decreasing = T),]
aaa<-length(rownames(volcano1))
genename<-c(rownames(volcano1)[1:20],rownames(volcano1)[(aaa-19):aaa])

express1<-as.data.frame(express.norm)
express1$ID<-rownames(express1)
phonedata<-dplyr::filter(express1,ID %in% genename)
phonedata<-dplyr::select(phonedata,-ID)
group_list<-group_list[order(group_list$condition,decreasing = T),]

group_list[,2]

annotation_col<-data.frame(group_list[,2])   
rownames(annotation_col)<-group_list[,1]
colnames(annotation_col)<-"Group"
phonedata<-select(phonedata,rownames(annotation_col))

pheatmap(phonedata,scale = "row",cluster_row=T,cluster_col=F, annotation_col = annotation_col,show_colnames = F
         ,border_color = F,color = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()
pdf(file="heatmap.PDF", width=6, height=6)
pheatmap(phonedata,scale = "row",cluster_row=T,cluster_col=F, annotation_col = annotation_col,show_colnames = F
         ,border_color = F,color = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()




rm(list = ls())
options(stringsAsFactors = FALSE)
library(tidyverse)

logFC_filter=0
P.Value_filter=0.05
merge<-read.csv("differential_expression_results.csv",header = T,sep = ",")
colnames(merge)[1]<-"ID"
merge<-filter(merge,P.Value < P.Value_filter)
dat<-filter(merge,logFC >logFC_filter)
merge_UP<-dat$ID
dat<-filter(merge,logFC < -(logFC_filter))
merge_DOWN<-dat$ID


Mitochondrion=read.table("Mitochondrion.txt",header=T)
Mitochondrion=Mitochondrion[,1]

library(VennDiagram)
library(RColorBrewer)

aaa<-list(merge_UP=merge_UP , Mitochondrion=Mitochondrion )


p = venn.diagram(aaa,fill = c(brewer.pal(7,"Set1")[1:2]),
                 alpha = c(0.5, 0.5), cex = 2,cat.dist = c(0,0.025),
                 cat.cex=1.5,lty =2, fontfamily ="sans",fontface = "bold",cat.fontfamily = "sans",cat.fontface = "bold", 
                 resolution =300,filename = NULL)

pdf("Intersection_with_upregulated_genes.pdf")
grid.draw(p)
dev.off()



aaa<-list(merge_DOWN=merge_DOWN , Mitochondrion=Mitochondrion)


p = venn.diagram(aaa,fill = c(brewer.pal(7,"Set1")[1:2]),
                 alpha = c(0.5,0.5), cex = 2,cat.dist = c(0, 0.025),
                 cat.cex=1.5,lty =2, fontfamily ="sans",fontface = "bold",cat.fontfamily = "sans",cat.fontface = "bold", 
                 resolution =300,filename = NULL)

pdf("Intersection_with_downregulated_genes.pdf")
grid.draw(p)
dev.off()



a<-intersect(merge_UP,Mitochondrion)


upgene<-unique(c(a))

d<-intersect(merge_DOWN,Mitochondrion)

downgene<-unique(c(d))
difgene<-c(upgene,downgene)

upgene<-data.frame(genename=upgene)

downgene<-data.frame(genename=downgene)

difgene<-data.frame(genename=difgene)
write.table(difgene,"difgene.txt",row.names = F,quote = F)
write.table(upgene,"upgene.txt",row.names = F,quote = F)
write.table(downgene,"downgene.txt",row.names = F,quote = F)
difgeneall<-c(merge_UP,merge_DOWN)
write.table(difgeneall,"ALLdifgene.txt",row.names = F,quote = F)


