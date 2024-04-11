
library(GEOquery)
my.gse<-c("GSE113376")
geo.gse <- getGEO(GEO=my.gse, filename=NULL, 
                  GSElimits=NULL, GSEMatrix=TRUE, AnnotGPL=FALSE, 
                  getGPL=FALSE)
geo.gse


my.geo.gse <- geo.gse[[1]]

my.geo.gse <- geo.gse[[1]]
gsms= "1111222223333344445555566666777788888XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
nchar(gsms)

sml= c()
for (i in 1:nchar(gsms)){sml[i]= substr(gsms, i, i)}

# eliminate ones we don't want ##Eliminated liver and BLOOD for now from here
##The architecture of sample arrngement was found from the GSE113375 GEO DATABASE VIEWER

sel= which(sml!= "X")
sml= sml[sel]
my.geo.gse= my.geo.gse[, sel ]

# The phenotype data is contained inside pData slot. 
head(pData(my.geo.gse),3)
dim(pData(my.geo.gse))
# the expression matrix is contained inside exprs slot
head(exprs(my.geo.gse),3)
dim(exprs(my.geo.gse))

pData(my.geo.gse)$title
class(pData(my.geo.gse)$title)
dim(exprs(my.geo.gse))

library("dplyr")
# make data frame of phenotypic Data
pdata <- as.data.frame(pData(my.geo.gse), stringsAsFactors=F)
colnames(pdata)

my.pdata=select(pdata,"title","geo_accession","source_name_ch1","organism_ch1","characteristics_ch1.1","characteristics_ch1.3","characteristics_ch1.4")

group_num<-c("1111222223333344445555566666777788888")
nchar(group_num)
group_names<-paste("Mouse",c("spleen_naive_d0","spleen_inf_d15",  
                             "spleen_inf_d17", "spleen_inf_d21", "spleen_inf_d36",
                             "spleen_naive_d36","spleen_inf_d42", "spleen_naive_d42"),sep="_")

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
write.csv(exprs(my.rma),"GSE113376_NORMEXPRDATA.csv")

library(limma)
library(ggrepel)
library(RColorBrewer)
pca <- prcomp(t(exprs(my.rma)),scale. = TRUE)
#pca2=pca$x
pca3=cbind(my.pdata, pca$x)
ggplot(data=pca3,aes(x = PC1, y= PC2, col=group))+geom_point(pch=19,cex=2)+theme_classic()+scale_color_brewer(palette="Dark2")

#+geom_text_repel(data = pca3, aes(label = group), hjust = 0, vjust = 0)


# prepare design matrix

pData(my.rma)$group <- as.factor(pData(my.rma)$group)
category=pData(my.rma)$group
class(category)
category=relevel(category,ref = "Mouse_spleen_naive_d0")
my.design=model.matrix(~0+category,pData(my.rma))
rownames(my.design) <-pData(my.rma)$group
colnames(my.design)=levels(category)

# fit linear model
fit <- lmFit(my.rma, my.design)

my.contrasts <- makeContrasts(Naive_vs_Day36=Mouse_spleen_inf_d36-Mouse_spleen_naive_d36,
                              Naive_Vs_Day42=Mouse_spleen_inf_d42-Mouse_spleen_naive_d42,
                              levels=colnames(my.design))
my.contrasts

write.csv(my.contrasts,"GSE113376_Spleen_Comparison.csv")

#EXPLANATION OF THIS PART OF CODE
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


library(GEOquery)
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


#day36vsDAY36naive
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
table(ji$Diffexpressed)

write.csv(ji,"113376_Day36naivespleen vs Day36 infected spleen overall expr values NEW.csv")
ji2=read.csv("113376_Day36naivespleen vs Day36 infected spleen overall expr values NEW.csv")
stat1_row <- ji2[ji2$gene == "Stat1", ]
gbp5_row<-ji2[ji2$gene == "Gbp5", ]
Fpr1_row=ji2[ji2$gene == "Fpr1", ]
table(ji$Diffexpressed)
ggplot(data=ji2,aes(x=log2FoldChange,y=-log10(pvalue),col=Diffexpressed,label=delabel))+
  geom_point(size=1)+
  theme_minimal()+
  scale_color_manual(values = c("red","black","blue"))+
  theme(text=element_text(size=15))+
  #geom_text_repel(data = gbp5_row, aes(label = gene),color="red", nudge_x = 0.1, nudge_y = 0.1, size = 3)+
  #geom_text_repel(data = stat1_row, aes(label = gene),color="red", nudge_x = 0.1, nudge_y = 0.1, size = 3)+
  #geom_text_repel(data = Fpr1_row, aes(label = gene),color="red", nudge_x = 0.1, nudge_y = 0.1, size = 3)+
  theme_classic()




ggplot(data=ji,aes(x=log2FoldChange,y=-log10(pvalue),col=Diffexpressed,label=delabel))+
  geom_point(size=1)+
  theme_minimal()+ylim(0,20)+
  scale_color_manual(values = c("red","black","green"))+
  theme(text=element_text(size=15))+theme_classic()

results= ji[,c(4,9,10)]
head(results)
results= results %>% mutate(new_p= ifelse(pvalue<0.05, TRUE, FALSE))
results$pvalue= NULL
str(results)
class(results)
#colnames(results)=c("baseMean","log2FoldChange","padj")
#results=as.data.frame(results)
#write.csv(results, "results.csv")

##MA PLOT
##Gene plotter function was not working
colors <- ifelse(results$new_p, "blue", "gray39")
plot(results$AveExpr, results$log2FoldChange, col = colors, pch = 19,cex=0.45,
     ylim = c(-5, 5), xlim =c(-1 , 1), xlab = "Mean Expression", ylab = "Log Fold Change", main = "Day_36_Spleen_Naive_Vs_Day_36_Spleen_infected", abline(h = 0, col = "gray39", lwd = 3))

#Day42vs42infected
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
write.csv(ki,"113376_Day42naivespleen vs Day42 infected spleen overall expr values NEW.csv")
ki2=read.csv("113376_Day42naivespleen vs Day42 infected spleen overall expr values NEW.csv")
stat1_row <- ki2[ki2$gene == "Stat1", ]
gbp5_row<-ki2[ki2$gene == "Gbp5", ]
Fpr1_row=ki2[ki2$gene == "Fpr1", ]
table(ki$Diffexpressed)
ggplot(data=ki2,aes(x=log2FoldChange,y=-log10(pvalue),col=Diffexpressed,label=delabel))+
  geom_point(size=1)+
  theme_minimal()+
  scale_color_manual(values = c("red","black","blue"))+
  theme(text=element_text(size=15))+
  geom_text_repel(data = gbp5_row, aes(label = gene),color="red", nudge_x = 0.1, nudge_y = 0.1, size = 3)+
  geom_text_repel(data = stat1_row, aes(label = gene),color="red", nudge_x = 0.1, nudge_y = 0.1, size = 3)+
  geom_text_repel(data = Fpr1_row, aes(label = gene),color="red", nudge_x = 0.1, nudge_y = 0.1, size = 3)+
  theme_classic()

ggplot(data=ki,aes(x=log2FoldChange,y=-log10(pvalue),col=Diffexpressed,label=delabel))+
  geom_point(size=1)+
  theme_minimal()+ylim(0,20)+
  scale_color_manual(values = c("red","black","green"))+
  theme(text=element_text(size=15))+theme_classic()

results= ki[,c(4,9,10)]
head(results)
results= results %>% mutate(new_p= ifelse(pvalue<0.05, TRUE, FALSE))
results$pvalue= NULL
str(results)
class(results)
#colnames(results)=c("baseMean","log2FoldChange","padj")
#results=as.data.frame(results)
#write.csv(results, "results.csv")


##MA PLOT
##Gene plotter function was not working
colors <- ifelse(results$new_p, "blue", "gray39")
plot(results$AveExpr, results$log2FoldChange, col = colors, pch = 19,cex=0.45,
     ylim = c(-5, 5), xlim =c(-1 , 1), xlab = "Mean Expression", ylab = "Log Fold Change", main = "Day_42_Spleen_Naive_Vs_Day_42_Spleen_infected", abline(h = 0, col = "gray39", lwd = 3))

