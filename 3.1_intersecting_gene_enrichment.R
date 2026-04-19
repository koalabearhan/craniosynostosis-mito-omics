rm(list=ls())

library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(tidyverse)

genelist_input <- fread(file="degs.txt", header = T, sep='\t', data.table = F)
genename <- as.character(genelist_input[,1]) 

gene_map <- biomaRt::select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID"))
gene_map
write.csv(as.data.frame(gene_map),"gene_conversion.csv",row.names =F)
genelist_input<-gene_map[,2]

head(genelist_input)
genelist_input<-na.omit(genelist_input)

Go_result_BP <- enrichGO(genelist_input, 'org.Hs.eg.db', ont="BP", pvalueCutoff=1) 

goplot(Go_result_BP, showCategory=5) 

p1<-dotplot(Go_result_BP, showCategory=20) 
p1<-p1 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(p1,filename ='GO_BP_dot_plot.pdf',width = 7.23,height = 8)
ggsave(p1,filename ='GO_BP_dot_plot.TIFF',width = 7.23,height = 8)

barplot(Go_result_BP, showCategory=20) 
y=as.data.frame(Go_result_BP)
y$geneID=as.character(sapply(y$geneID,function(x)paste(gene_map$SYMBOL[match(strsplit(x,"/")[[1]],as.character(gene_map$ENTREZID))],collapse="/")))
write.csv(y,"GO-BP.csv",row.names =F)
save(y,file = 'GO-BP.RDATA')

Go_result_CC <- enrichGO(genelist_input, 'org.Hs.eg.db', ont="CC", pvalueCutoff=10) 

goplot(Go_result_CC, showCategory=5) 

p1<-dotplot(Go_result_CC, showCategory=20) 
p1<-p1 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(p1,filename ='GO_CC_dot_plot.pdf',width = 7.23,height = 8)
ggsave(p1,filename ='GO_CC_dot_plot.TIFF',width = 7.23,height = 8)

barplot(Go_result_CC, showCategory=20) 
y=as.data.frame(Go_result_CC)
y$geneID=as.character(sapply(y$geneID,function(x)paste(gene_map$SYMBOL[match(strsplit(x,"/")[[1]],as.character(gene_map$ENTREZID))],collapse="/")))
write.csv(y,"GO-CC.csv",row.names =F)
save(y,file = 'GO-CC.RDATA')

Go_result_MF <- enrichGO(genelist_input, 'org.Hs.eg.db',ont="MF", pvalueCutoff=1000) 

goplot(Go_result_MF, showCategory=5) 

p1<-dotplot(Go_result_MF, showCategory=20) 
p1<-p1 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(p1,filename ='GO_MF_dot_plot.pdf',width = 7.23,height = 8)
ggsave(p1,filename ='GO_MF_dot_plot.TIFF',width = 7.23,height = 8)

barplot(Go_result_MF, showCategory=20) 
y=as.data.frame(Go_result_MF)
y$geneID=as.character(sapply(y$geneID,function(x)paste(gene_map$SYMBOL[match(strsplit(x,"/")[[1]],as.character(gene_map$ENTREZID))],collapse="/")))

write.csv(y,"GO-MF.csv",row.names =F)
save(y,file = 'GO-MF.RDATA')



go <- enrichGO(genelist_input, OrgDb = "org.Hs.eg.db", ont="all")
library(ggplot2)
p <- dotplot(go, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale="free")
p <-p + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(p,filename ='GO_dot_plot.pdf',width = 7.23,height = 8)
ggsave(p,filename ='GO_dot_plot.TIFF',width = 7.23,height = 8)
write.csv(as.data.frame(go),"GO.csv",row.names =F)




KEGG_result <- enrichKEGG(genelist_input, keyType = "kegg",pvalueCutoff=1,qvalueCutoff=1,pAdjustMethod = "BH", minGSSize = 5, maxGSSize = 500,organism = "hsa", use_internal_data=T)  #KEGG富集分析


barplot(KEGG_result, showCategory=20)

p1<-dotplot(KEGG_result, showCategory=20) 
p1<-p1 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(p1,filename ='KEGG_dot_plot.pdf',width = 7.23,height = 8)
ggsave(p1,filename ='KEGG_dot_plot.TIFF',width = 7.23,height = 8)



pdf(file="KEGG_circos.pdf",width = 10,height = 7)
kkx=setReadable(KEGG_result, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(kkx, showCategory = 5, circular = TRUE, colorEdge = TRUE,node_label="all")
dev.off()
x=as.data.frame(KEGG_result)
x$geneID=as.character(sapply(x$geneID,function(x)paste(gene_map$SYMBOL[match(strsplit(x,"/")[[1]],as.character(gene_map$ENTREZID))],collapse="/")))


write.csv(as.data.frame(x),"KEGG.csv",row.names =F)
save(x,file = 'KEGG.RDATA')



y =KEGG_result@result


forward <- as.numeric(sub("/\\d+$", "", y$GeneRatio))
backward <- as.numeric(sub("^\\d+/", "", y$GeneRatio))

y$GeneRatio = forward/backward

showCategory =20

font.size =12

library(ggplot2)
library(forcats)
library(dplyr)

AAA=y %>% 
 
  arrange(pvalue) %>% 
  slice(1:showCategory)
  
  p1=ggplot(AAA,aes(GeneRatio,forcats::fct_reorder(Description,Count)))+ 
 
  geom_point(aes(color=pvalue, size = Count)) +
 
  scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
 
  scale_size_continuous(range=c(3, 8))+
  
  labs(y=NULL) +
  
  ggtitle("")+
  
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = font.size, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = font.size, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = font.size),
        axis.title.y = element_text(angle=90))

p1<-p1 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(filename ='KEGG_dot_plot2.pdf',width = 8,height = 8)
ggsave(filename ='KEGG_dot_plot2.TIFF',width = 8,height = 8)


y =Go_result_MF@result


forward <- as.numeric(sub("/\\d+$", "", y$GeneRatio))
backward <- as.numeric(sub("^\\d+/", "", y$GeneRatio))

y$GeneRatio = forward/backward

showCategory =20

font.size =12

library(ggplot2)
library(forcats)
library(dplyr)

AAA=y %>% 
  
  arrange(pvalue) %>% 
  slice(1:showCategory)

p1=ggplot(AAA,aes(GeneRatio,forcats::fct_reorder(Description,Count)))+ 
  
  geom_point(aes(color=pvalue, size = Count)) +
  
  scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
  
  scale_size_continuous(range=c(3, 8))+
  
  labs(y=NULL) +
 
  ggtitle("")+
  
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = font.size, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = font.size, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = font.size),
        axis.title.y = element_text(angle=90))

p1<-p1 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(filename ='Go_result_MF2.pdf',width = 8,height = 8)
ggsave(filename ='Go_result_MF2.TIFF',width = 8,height = 8)



y =Go_result_BP@result


forward <- as.numeric(sub("/\\d+$", "", y$GeneRatio))
backward <- as.numeric(sub("^\\d+/", "", y$GeneRatio))

y$GeneRatio = forward/backward

showCategory =20

font.size =12

library(ggplot2)
library(forcats)
library(dplyr)

AAA=y %>% 
 
  arrange(pvalue) %>% 
  slice(1:showCategory)

p1=ggplot(AAA,aes(GeneRatio,forcats::fct_reorder(Description,Count)))+ 
 
  geom_point(aes(color=pvalue, size = Count)) +

  scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
 
  scale_size_continuous(range=c(3, 8))+

  labs(y=NULL) +
  
  ggtitle("")+
  
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = font.size, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = font.size, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = font.size),
        axis.title.y = element_text(angle=90))

p1<-p1 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(filename ='Go_result_BP2.pdf',width = 8,height = 8)
ggsave(filename ='Go_result_BP2.TIFF',width = 8,height = 8)


y =Go_result_CC@result


forward <- as.numeric(sub("/\\d+$", "", y$GeneRatio))
backward <- as.numeric(sub("^\\d+/", "", y$GeneRatio))

y$GeneRatio = forward/backward

showCategory =20

font.size =12

library(ggplot2)
library(forcats)
library(dplyr)

AAA=y %>% 
  
  arrange(pvalue) %>% 
  slice(1:showCategory)

p1=ggplot(AAA,aes(GeneRatio,forcats::fct_reorder(Description,Count)))+ 
  
  geom_point(aes(color=pvalue, size = Count)) +
 
  scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
 
  scale_size_continuous(range=c(3, 8))+

  labs(y=NULL) +
 
  ggtitle("")+
 
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = font.size, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = font.size, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = font.size),
        axis.title.y = element_text(angle=90))

p1<-p1 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(filename ='Go_result_CC2.pdf',width = 8,height = 8)
ggsave(filename ='Go_result_CC2.TIFF',width = 8,height = 8)




library(ggpubr)
library(tidyverse)
GODATA<-go@result
GODATA$"-log10(Pvalue)"<-  -log10(GODATA$pvalue)
GODATA$yyy<-  -log10(GODATA$pvalue)
colnames(GODATA)


write.csv(GODATA,"GODATA.csv",row.names =F)


aaa<-filter(GODATA,ONTOLOGY=='BP')
aaaa<-aaa[1:10,]

bbb<-filter(GODATA,ONTOLOGY=='CC')
bbbb<-bbb[1:10,]

ccc<-filter(GODATA,ONTOLOGY=='MF')
cccc<-ccc[1:10,]


drawdata<-rbind(aaaa,bbbb,cccc)


ggbarplot(drawdata, x = "Description", y = 'yyy',
          fill = "ONTOLOGY",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          palette = c('#5FB404','#01DFD7','#C238E5'),            # jco journal color palett. see ?ggpar
          sort.val = "desc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90           # Rotate vertically x axis texts
)
ggbarplot(drawdata, x = "Description", y = "yyy",
          fill = "ONTOLOGY",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          palette = c('#5FB404','#12B5EC','#C238E5'),            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in dscending order
          sort.by.groups = TRUE,      # Sort inside each group
          x.text.angle = 75, # Rotate vertically x axis texts
          ylab = '-log10(P-Value)',
          xlab = 'Pathway'
)
ggsave(filename ='GO2.pdf',width = 7.23,height = 8)
ggsave(filename ='GO2.TIFF',width = 7.23,height = 8)
