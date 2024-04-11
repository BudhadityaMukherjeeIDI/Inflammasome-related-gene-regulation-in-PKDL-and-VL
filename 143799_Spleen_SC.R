spleen=read.csv("spleen_counts.csv",header=TRUE)
head(spleen)
dim(spleen)
spleen=na.omit(spleen)
write.csv(spleen,"spleen_nona.csv")
spleen1=read.csv("spleen_nona.csv",header=TRUE,row.names = 1)
head(spleen1)
dim(spleen1)

sampleinfo=read.csv("sampleinfo.csv",header = TRUE,row.names = 1)
head(sampleinfo)
dim(sampleinfo)

sampleinfo$Levels=factor(sampleinfo$Levels)
class(sampleinfo$Levels)
library("DESeq2")

dds=DESeqDataSetFromMatrix(countData = spleen1, colData = sampleinfo, design = ~Levels )
keep= rowSums(counts(dds)) > 0
dds=dds[keep,]
nrow(dds)

dds$Levels <- relevel(dds$Levels, ref = "Uninfected_Spleen")
levels(dds$Levels)
dds=DESeq(dds)
library(ggplot2)

normcount=counts(dds,normalized=T)
class(normcount)
normcounts=as.data.frame(normcount)
nrow(normcounts)
write.csv(normcounts,"normcounts_143799.csv")
pca <- prcomp(t(normcounts),scale. = T)
pca3=cbind(sampleinfo, pca$x)
library("RColorBrewer")
ggplot(data=pca3,aes(x = PC1, y= PC2, col=Levels))+ geom_point()+theme_classic()+scale_color_brewer(palette="Dark2")

resultsNames(dds)
#"Levels_Infected_Donovani_vs_Uninfected_Spleen"
donovani=results(dds,name="Levels_Infected_Donovani_vs_Uninfected_Spleen", alpha = 0.05)
donovani=donovani[order(donovani$padj),]
class(donovani)
donovani=as.data.frame(donovani)
write.csv(donovani,"143799_donovani_spleen_ordered_padj.csv")
dim(donovani)

res=results(dds,name="Levels_Infected_Donovani_vs_Uninfected_Spleen", alpha = 0.05)
##MA PLOT
DESeq2::plotMA(res, ylim=c(-5,5), cex=0.45,colSig="blue",colLine="grey39", colNonSig="gray39", main= "Murine_Infected_Donovani_vs_Uninfected_Spleen")
abline(h = 0, col = "gray39", lwd = 3)

library(dbplyr)

filtered=donovani %>% dplyr::filter(donovani$padj < 0.05)
write.csv(filtered,"143799_spleen_statistically_sig_genes.csv")
dim(filtered)

reslfc=lfcShrink(dds,coef="Levels_Infected_Donovani_vs_Uninfected_Spleen", type="apeglm")
reslfc=as.data.frame(reslfc)
dim(reslfc)
#reslfc$log2FoldChange < -1
#reslfc=na.omit(reslfc)
reslfc$Diffexpressed= 'NO'
reslfc$Diffexpressed[reslfc$log2FoldChange>1 & reslfc$padj<0.05]="UP"
reslfc$Diffexpressed[reslfc$log2FoldChange< -1 & reslfc$padj<0.05]="DOWN"
reslfc$delabel=NA

table(reslfc$Diffexpressed)

write.csv(reslfc,"143799_donovani_SPLEEN_after_shrinkage_overall expression values.csv")
ji2=read.csv("143799_donovani_SPLEEN_after_shrinkage_overall expression values.csv")
stat1_row <- ji2[ji2$X == "Stat1", ]
gbp5_row<-ji2[ji2$X == "Gbp5", ]
Fpr1_row=ji2[ji2$X == "Fpr1", ]
ggplot(data=ji2,aes(x=log2FoldChange,y=-log10(pvalue),col=Diffexpressed,label=delabel))+
  geom_point(size=1)+
  theme_minimal()+ylim(0,150)+
  scale_color_manual(values = c("red","black","blue"))+
  theme(text=element_text(size=15))+
  #geom_text_repel(data = gbp5_row, aes(label = X),color="black", nudge_x = 0.1, nudge_y = 0.1, size = 4)+
  #geom_text_repel(data = stat1_row, aes(label = X),color="black", nudge_x = 0.1, nudge_y = 0.1, size = 4)+
  #geom_text_repel(data = Fpr1_row, aes(label = X),color="black", nudge_x = 0.1, nudge_y = 0.1, size = 4)+
  theme_classic()
ggplot(data=reslfc,aes(x=log2FoldChange,y=-log10(pvalue),col=Diffexpressed,label=delabel))+
  geom_point(size=1)+
  theme_minimal()+ylim(0,200)+
  scale_color_manual(values = c("red","black","green"))+
  theme(text=element_text(size=15))+theme_classic()

ggplot(data=reslfc,aes(x=log2FoldChange,y=-log10(pvalue),col=Diffexpressed,label=delabel))+
  geom_point(size=1)+
  theme_minimal()+ylim(0,20)+
  scale_color_manual(values = c("red","black","green"))+
  theme(text=element_text(size=15))+theme_classic()


#Levels_Infected_Infantum_vs_Uninfected_Spleen
infantum=results(dds,name="Levels_Infected_Infantum_vs_Uninfected_Spleen", alpha = 0.05)
infantum=infantum[order(infantum$padj),]
class(infantum)
infantum=as.data.frame(infantum)
write.csv(infantum,"143799_infantum_ordered_padj.csv")

res2=results(dds,name="Levels_Infected_Infantum_vs_Uninfected_Spleen", alpha = 0.05)
##MA PLOT
DESeq2::plotMA(res2, ylim=c(-5,5), cex=0.45,colSig="blue",colLine="grey39", colNonSig="gray39", main= "Murine_Infected_Infantum_vs_Uninfected_Spleen")
abline(h = 0, col = "gray39", lwd = 3)

library(dbplyr)

filtered3=infantum %>% dplyr::filter(infantum$padj < 0.05)
write.csv(filtered3,"143799_infantum_statistically_sig_genes.csv")
dim(filtered3)


#Volcano Plot
BiocManager::install("apeglm")
library("apeglm")
reslfc1=lfcShrink(dds,coef="Levels_Infected_Infantum_vs_Uninfected_Spleen", type="apeglm")
reslfc1=as.data.frame(reslfc1)
dim(reslfc1)
#reslfc$log2FoldChange < -1
#reslfc=na.omit(reslfc)
reslfc1$Diffexpressed= 'NO'
reslfc1$Diffexpressed[reslfc1$log2FoldChange>1 & reslfc1$padj<0.05]="UP"
reslfc1$Diffexpressed[reslfc1$log2FoldChange< -1 & reslfc1$padj<0.05]="DOWN"
reslfc1$delabel=NA

table(reslfc1$Diffexpressed)

write.csv(reslfc1,"143799_INFANTUM_spleen_after_shrinkage_overall expression values NEW.csv")
ki2=read.csv("143799_INFANTUM_spleen_after_shrinkage_overall expression values NEW.csv")
stat1_row <- ki2[ki2$X == "Stat1", ]
gbp5_row<-ki2[ki2$X == "Gbp5", ]
Fpr1_row=ki2[ki2$X == "Fpr1", ]
ggplot(data=ki2,aes(x=log2FoldChange,y=-log10(pvalue),col=Diffexpressed,label=delabel))+
  geom_point(size=1)+
  theme_minimal()+ylim(0,150)+
  scale_color_manual(values = c("red","black","blue"))+
  theme(text=element_text(size=15))+
  #geom_text_repel(data = gbp5_row, aes(label = X),color="black", nudge_x = 0.1, nudge_y = 0.1, size = 4)+
  #geom_text_repel(data = stat1_row, aes(label = X),color="black", nudge_x = 0.1, nudge_y = 0.1, size = 4)+
  #geom_text_repel(data = Fpr1_row, aes(label = X),color="black", nudge_x = 0.1, nudge_y = 0.1, size = 4)+
  theme_classic()
ggplot(data=reslfc1,aes(x=log2FoldChange,y=-log10(pvalue),col=Diffexpressed,label=delabel))+
  geom_point(size=1)+
  theme_minimal()+ylim(0,200)+
  scale_color_manual(values = c("red","black","green"))+
  theme(text=element_text(size=15))+theme_classic()


ggplot(data=reslfc1,aes(x=log2FoldChange,y=-log10(pvalue),col=Diffexpressed,label=delabel))+
  geom_point(size=1)+
  theme_minimal()+ylim(0,20)+
  scale_color_manual(values = c("red","black","green"))+
  theme(text=element_text(size=15))+theme_classic()


