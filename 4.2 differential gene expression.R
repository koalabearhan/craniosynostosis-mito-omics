rm(list = ls())
library(tidyverse)
library(ggpubr)
rt<-read.csv("merge_matrix_new.csv",header = T,sep = ",",row.names = 1)
clindata<-read.table("clinicaldata.txt",header = T,sep = "\t")
ccc=intersect(clindata$Accession,colnames(rt))
rt=rt[,ccc]
range(rt)


fegene<-read.table("degs.txt",header = T,sep = "\t")
fegene<-fegene[,1]
aaa<-intersect(fegene,rownames(rt))
rt=rt[aaa,]
rt<-as.data.frame(t(rt))
rt$ID<-rownames(rt)


colnames(clindata)[1]<-"ID"

drawdara<-inner_join(clindata,rt,by="ID")
gene=fegene[1]

pdata<-drawdara
colnames(pdata)[2]<-"group"
pdata_melt <- reshape2::melt(pdata,
                             id.vars = c("ID","group"),
                             variable.name ="genename",
                             value.name = "expression")
pdata_melt$expression<-log2(pdata_melt$expression)

c <- ggplot(pdata_melt,
            aes(x=genename, y=expression, 
                fill = group, 
                color = group)) + 
  geom_boxplot(notch = F, alpha = 0.95, 
               outlier.shape = 16,
               outlier.colour = "black", 
               outlier.size = 0.65) +
  
  scale_fill_manual(values= c("#E7B800","#E31A1C","#2E9FDF")) +
  #ggtitle(paste0(Index[[i]],"_signature score")) + 
  labs("gene expression")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5, size = 10), 
        axis.text.y = element_text(angle = 90, size = 12),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position = "top")+
  theme(plot.title = element_text(hjust = 0.5)) 


# `****` = 1e-04, `***` = 0.001, `**` = 0.01, `*` = 0.05, ns = 1
p <- c + stat_compare_means(label = "p.signif",method = "t.test")


ggsave(p, filename =paste0("merge","differential_expression_plot.pdf"), width = 10, height = 6)

drawdata=drawdara
drawdata$condition<-factor(drawdata$condition,levels = c('Craniosynostosis','Control'))
my_comparisons <- list( c("Craniosynostosis", "Control") )

Index=fegene
Index=str_replace_all(Index,"-",".")
colnames(drawdata)=str_replace_all(colnames(drawdata),"-",".")
for (gene in Index) {
  p <- ggboxplot(drawdata, x = "condition", y = gene,
                 fill = "condition",legend=F,palette =c("#E7B800", "#00AFBB"),bxp.errorbar=T)+
    theme(legend.position='none')+
    ylab(label = paste(gene,'expression' ,sep = ' ' ))+
    stat_compare_means(comparisons = my_comparisons,method = 't.test',aes(label = ..p.signif..))
  ggsave(p,filename = paste(gene,'differential_expression.pdf',sep=' '),width = 5.6,height = 4.22)
  print(gene)
}


































