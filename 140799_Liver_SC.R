library(GEOquery)
my.gse<-c("GSE140799")
geo.gse <- getGEO(GEO=my.gse, filename=NULL, 
                  GSElimits=NULL, GSEMatrix=TRUE, AnnotGPL=FALSE, 
                  getGPL=FALSE)
geo.gse


my.geo.gse <- geo.gse[[1]]
pdata=pData(my.geo.gse)
nrow(pdata)

gsms="XXXXXXXXXXXXXXXXXXXXXXXXXXXX11111222223333344444555556666"
nchar(gsms)
sml= c()
for (i in 1:nchar(gsms)){sml[i]= substr(gsms, i, i)}
sel= which(sml!= "X")

sml= sml[sel]
my.geo.gse= my.geo.gse[, sel ]

pData(my.geo.gse)$title
library("dplyr")
# make data frame of phenotypic Data
pdata <- as.data.frame(pData(my.geo.gse), stringsAsFactors=F)
colnames(pdata)
my.pdata=select(pdata,"title","geo_accession","source_name_ch1","organism_ch1","characteristics_ch1.1","characteristics_ch1.3","characteristics_ch1.4")

group_num=c("11111222223333344444555556666")
group_names<-paste("Mouse",c("liver_d42_Inf_Unreated","liver_d42_Inf_treated",  
                             "liver_d36_Inf_Unreated", "liver_d36_naive", "liver_d36_Inf_treated",
                             "liver_d42_naive"),sep="_")

group_by_numbering<-function(group_numbering, vector, group_names=NULL){
  #     gsub(pattern, replacement, string)
  group_numbering<- gsub(" ","",group_numbering)
  group <- rep(NA,length(vector))
  for (i in 1:nchar(group_numbering)){
    num<-as.numeric(substr(group_numbering,i,i))
    if(is.null(group_names)) 
      group[i] <- num
    else group[i] <- group_names[num]
  }
  return(group)
}

# allocate groups
my.pdata[["group"]]<-group_by_numbering(group_num, pdata$title, group_names)

table(my.pdata$group)
assayData<-exprs(my.geo.gse)

# check if data need to be normalized 
need.log.transformation<- function(data){
  qx <- as.numeric(quantile(data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  ##This line calculates the quantiles of the input data using the quantile 
  #function. The quantiles at specific percentiles (0%, 25%, 50%, 75%, 99%, and 100%) are calculated.
  #remove any NA (missing) values in the data when calculating the quantiles.
  # check if data need to be normalized 
  LogC <- (qx[5] > 100) ||   ##This condition checks if the 99th percentile of the data is greater than 100.
    #If so, it suggests that there may be very large values in the data
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  return(LogC)
}
LogC<-need.log.transformation(assayData)
LogC

# Normalisation
if (LogC) { 
  assayData[which(assayData<= 0)] <- NaN
  # take log2 transform
  assayData<-na.omit(log2(assayData))
}

my.rma<-ExpressionSet(assayData=assayData,
                      phenoData=AnnotatedDataFrame(my.pdata))
head(exprs(my.rma))
class(exprs(my.rma))
write.csv(exprs(my.rma),"GSE140799_NORMEXPRDATA.csv")

library(limma)
library(ggrepel)
library(RColorBrewer)
pca <- prcomp(t(exprs(my.rma)),scale. = TRUE)
#pca2=pca$x
pca3=cbind(my.pdata, pca$x)
ggplot(data=pca3,aes(x = PC1, y= PC2, col=group))+geom_point(pch=19,cex=2)+theme_classic()+scale_color_brewer(palette="Dark2")

# prepare design matrix
# check if rows of phenotype is sorted same way as columns of exprs matrix
table(rownames(pData(my.rma))==colnames(exprs(my.rma)))
pData(my.rma)$group <- as.factor(pData(my.rma)$group)
category=pData(my.rma)$group
class(category)
category=relevel(category,ref = "Mouse_liver_d36_naive")
my.design=model.matrix(~0+category,pData(my.rma))
rownames(my.design) <-pData(my.rma)$group
colnames(my.design)=levels(category)
# fit linear model
fit <- lmFit(my.rma, my.design)

my.contrasts <- makeContrasts(Day36n_Vs_Day36=Mouse_liver_d36_Inf_Unreated-Mouse_liver_d36_naive,
                              Day42n_Vs_Day42=Mouse_liver_d42_Inf_Unreated-Mouse_liver_d42_naive,
                              levels = colnames(my.design))

my.contrasts
write.csv(my.contrasts,"GSE140799_LIVER_Comparison.csv")

fits.fun<-function(x)(contrasts.fit(fit, contrasts=my.contrasts[, x]));

contrast.fits <- lapply(colnames(my.contrasts),fits.fun)
names(contrast.fits)<-colnames(my.contrasts)
ebs.fun<-function(x)(eBayes(x, proportion=0.1, trend=FALSE, robust=FALSE))

contrast.ebs <- lapply(contrast.fits, ebs.fun)
tts.fun<-function(x)(topTable(x, adjust.method="BH", number=length(x$coefficients), sort.by="none"));

contrast.tts <- lapply(contrast.ebs,tts.fun)

test.fun<-function(x){
  store<-decideTests(x, method="separate", 
                     adjust.method="BH", p.value=0.05, 
                     lfc=1)
  return(store)
};


contrast.tests <- lapply(contrast.ebs, test.fun)

contrast.tests.df<-lapply(contrast.tests,function(x){
  store<-as.data.frame(x)
  colnames(store)<-"test"
  return(store)
})

head(contrast.tests.df[[1]])

# show in one matrix whether up or down regulated in diff comparisons
tests.mat <- do.call(cbind, contrast.tests)
colnames(tests.mat) <- names(contrast.tests)
head(tests.mat,3)
table(tests.mat)#-1 is downregulated, +1 is upregulated,0 is no change

gpl <- getGEO('GPL13912', destdir=".")
View(gpl)
tab= Table(gpl)
colnames(Table(gpl))
gene= tab$GENE_SYMBOL
ID= tab$ID
gene= as.data.frame(gene)
dim(gene)
ID= as.data.frame(ID)
dim(ID)
e= cbind(gene, ID)

#DAY36NAIVE VS DAY36 LIVER INFECTED(untreated)
ji= contrast.tts[[1]]
ji= as.data.frame(ji)
ID= rownames(ji)
ji= cbind(ID, ji)
rownames(ji)= NULL
class(ji$ID)
ji$ID=as.integer(ji$ID)
ji$log2FoldChange= ji$logFC
ji$pvalue= ji$adj.P.Val
ji= inner_join(e, ji, by= "ID")

ji$Diffexpressed='NO'
ji$Diffexpressed[ji$log2FoldChange>1 & ji$pvalue<0.05]="UP"
ji$Diffexpressed[ji$log2FoldChange< -1 & ji$pvalue<0.05]="DOWN"
ji$delabel=NA
View(ji)

write.csv(ji,"140799_Day36naive_LIver vs Day36 infected untreated Liver overall expr values NEW.csv")
ji2=read.csv("140799_Day36naive_LIver vs Day36 infected untreated Liver overall expr values NEW.csv")
stat1_row <- ji2[ji2$gene == "Stat1", ]
gbp5_row<-ji2[ji2$gene == "Gbp5", ]
Fpr1_row=ji2[ji2$gene == "Fpr1", ]

table(ji$Diffexpressed)
library(ggplot2)
library(ggrepel)
ggplot(data=ji2,aes(x=log2FoldChange,y=-log10(pvalue),col=Diffexpressed,label=delabel))+
  geom_point(size=1)+
  theme_minimal()+
  scale_color_manual(values = c("red","black","blue"))+
  theme(text=element_text(size=15))+
  #geom_text_repel(data = gbp5_row, aes(label = gene),color="black", nudge_x = 0.1, nudge_y = 0.1, size = 4)+
  #geom_text_repel(data = stat1_row, aes(label = gene),color="black", nudge_x = 0.1, nudge_y = 0.1, size = 4)+
  #geom_text_repel(data = Fpr1_row, aes(label = gene),color="black", nudge_x = 0.1, nudge_y = 0.1, size = 4)+
  theme_classic()

ggplot(data=ji,aes(x=log2FoldChange,y=-log10(pvalue),col=Diffexpressed,label=delabel))+
  geom_point(size=1)+
  theme_minimal()+
  scale_color_manual(values = c("red","black","green"))+
  theme(text=element_text(size=15))+theme_classic()

ggplot(data=ji,aes(x=log2FoldChange,y=-log10(pvalue),col=Diffexpressed,label=delabel))+
  geom_point(size=1)+
  theme_minimal()+ylim(0,20)+
  scale_color_manual(values = c("red","black","green"))+
  theme(text=element_text(size=15))+theme_classic()

#DAY42NAIVE VS DAY42 INFECTED UNTREATED LIVER
ki= contrast.tts[[2]]
ki= as.data.frame(ki)
ID= rownames(ki)
ki= cbind(ID, ki)
rownames(ki)= NULL
class(ki$ID)
ki$ID=as.integer(ki$ID)
ki$log2FoldChange= ki$logFC
ki$pvalue= ki$adj.P.Val
ki= inner_join(e, ki, by= "ID")

ki$Diffexpressed='NO'
ki$Diffexpressed[ki$log2FoldChange>1 & ki$pvalue<0.05]="UP"
ki$Diffexpressed[ki$log2FoldChange< -1 & ki$pvalue<0.05]="DOWN"
ki$delabel=NA
View(ki)
write.csv(ki,"140799_Day42 naive_LIver vs Day42 infected untreated Liver overall expr values NEW.csv")
ji2=read.csv("140799_Day42 naive_LIver vs Day42 infected untreated Liver overall expr values NEW.csv")

stat1_row <- ji2[ji2$gene == "Stat1", ]
gbp5_row<-ji2[ji2$gene == "Gbp5", ]
Fpr1_row=ji2[ji2$gene == "Fpr1", ]

table(ji$Diffexpressed)
library(ggplot2)
library(ggrepel)
ggplot(data=ji2,aes(x=log2FoldChange,y=-log10(pvalue),col=Diffexpressed,label=delabel))+
  geom_point(size=1)+
  theme_minimal()+
  scale_color_manual(values = c("red","black","blue"))+
  theme(text=element_text(size=15))+
  geom_text_repel(data = gbp5_row, aes(label = gene),color="black", nudge_x = 0.1, nudge_y = 0.1, size = 4)+
  geom_text_repel(data = stat1_row, aes(label = gene),color="black", nudge_x = 0.1, nudge_y = 0.1, size = 4)+
  geom_text_repel(data = Fpr1_row, aes(label = gene),color="black", nudge_x = 0.1, nudge_y = 0.1, size = 4)+
  theme_classic()


table(ki$Diffexpressed)
ggplot(data=ki,aes(x=log2FoldChange,y=-log10(pvalue),col=Diffexpressed,label=delabel))+
  geom_point(size=1)+
  theme_minimal()+
  scale_color_manual(values = c("red","black","green"))+
  theme(text=element_text(size=15))+theme_classic()

ggplot(data=ki,aes(x=log2FoldChange,y=-log10(pvalue),col=Diffexpressed,label=delabel))+
  geom_point(size=1)+
  theme_minimal()+ylim(0,20)+
  scale_color_manual(values = c("red","black","green"))+
  theme(text=element_text(size=15))+theme_classic()
