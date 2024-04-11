library(GEOquery)

my.gse<-c("GSE134661")


#get published processed data and metadata from GEO
geo.gse <- getGEO(GEO=my.gse, filename=NULL, 
                  GSElimits=NULL, GSEMatrix=TRUE, AnnotGPL=FALSE, 
                  getGPL=FALSE)
geo.gse



# data available from only one platform GPL16025
my.geo.gse <- geo.gse[[1]]
# object is now an ExpressionSet
class(my.geo.gse)
View(my.geo.gse)

#The phenotype data is contained inside pData slot. 
head(pData(my.geo.gse),3)
dim(pData(my.geo.gse))
# the expression matrix is contained inside exprs slot
head(exprs(my.geo.gse),3)
dim(exprs(my.geo.gse))

pData(my.geo.gse)$title

library("dplyr")
# make data frame of phenotypic Data
pdata <- as.data.frame(pData(my.geo.gse), stringsAsFactors=F)
colnames(pdata)
class(pdata)

##SUBSETTING THE COLUMNS OF OUR INTEREST
my.pdata=dplyr::select(pdata,"title","geo_accession","source_name_ch1","organism_ch1","characteristics_ch1.1")
head(my.pdata)
my.pdata
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

group_num<-c("1122334455")
group_names<-paste("Mouse_spleen",c("naive","v_4wks","v_16wks", "nv_4wks", "nv_16wks"),sep="_")
# allocate groups
my.pdata[["group"]]<-group_by_numbering(group_num, pdata$title, group_names)
# num of replicates in each group
table(my.pdata$group)

##store processed phenotype data
write.table(pdata, "GSE134661_PhenoData_1.txt", 
            sep="\t", quote=F,col.names=NA)

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

if (LogC) { 
  assayData[which(assayData<= 0)] <- NaN
  # take log2 transform
  assayData<-na.omit(log2(assayData))
}

my.rma<-ExpressionSet(assayData=assayData,
                      phenoData=AnnotatedDataFrame(my.pdata))
head(exprs(my.rma))
class(exprs(my.rma))
write.csv(exprs(my.rma),"134661_NORMEXPRDATA_1.csv")

library(limma)
library(ggrepel)
pca <- prcomp(t(exprs(my.rma)))
dim(exprs(my.rma))
dim(my.pdata)
pca3=cbind(my.pdata, pca$x)
ggplot(data=pca3,aes(x = PC1, y= PC2, col=group)) + geom_point()+theme_classic()



table(exprs(my.rma)<0)

as.numeric(quantile(exprs(my.rma), c(0., 0.025, 0.05, 0.075, 0.099, .10), na.rm=T))


# check if rows of phenotype is sorted same way as columns of exprs matrix
table(rownames(pData(my.rma))==colnames(exprs(my.rma)))

# prepare design matrix
sample.lst<-pData(my.rma)$group
group.lst<-factor(sample.lst)

my.design <- model.matrix(~0 + group, pData(my.rma))

rownames(my.design) <- sample.lst
colnames(my.design) <- levels(group.lst)

my.design
fit <- lmFit(my.rma, my.design)

# specify the comparison btw group of interest
my.contrasts <- makeContrasts(n_vs_nv_4wks= Mouse_spleen_nv_4wks-Mouse_spleen_naive,
                              n_vs_nv_16wks= Mouse_spleen_nv_16wks-Mouse_spleen_naive,
                              levels = my.design)
my.contrasts

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
                     lfc=0)
  return(store)
};

contrast.tests <- lapply(contrast.ebs, test.fun)

contrast.tests.df<-lapply(contrast.tests,function(x){
  store<-as.data.frame(x)
  colnames(store)<-"test"
  return(store)
})

head(contrast.tests.df[[1]])
tests.mat <- do.call(cbind, contrast.tests)
colnames(tests.mat) <- names(contrast.tests)

library(dplyr)

y= contrast.tts
aa= as.data.frame(contrast.tts)
dim(aa)
head(aa)

gpl <- getGEO('GPL16570', destdir=".")
View(gpl)
tab= Table(gpl)
gene= tab$gene_assignment
Probe_ID= tab$ID
gene= as.data.frame(gene)
dim(gene)
Probe_ID= as.data.frame(Probe_ID)
dim(Probe_ID)

gene_probe= cbind(gene, Probe_ID)
View(gene_probe)

install.packages("reshape")
library(reshape)

#If we see in the table of gpl we can see that the gene assignment contains
#values other than that of the gene symbol hence we need to split that
##delimitter is //

e= tidyr::separate(gene_probe, gene , into = c("Col1", "Col2", "Col3", "Col4"), sep="\\//")
e$Col4= NULL
e$Col1= NULL

e= e[complete.cases(e[ , 1]),] #To remove the genes whose symbol is not available
head(e)

dim(aa)
dim(e)

colnames(e)=c("Symbol","Description","Probe_ID")
vi=contrast.tts[[2]]
vi=as.data.frame(vi)
Probe_ID= rownames(vi)
dim(vi)
vi= cbind(Probe_ID, vi)
rownames(vi)= NULL
vi$log2FoldChange= vi$logFC
vi$pvalue= vi$adj.P.Val
vi= inner_join(e, vi, by= "Probe_ID")
dim(vi)

##VOLCANOPLOT
vi$diffexpressed='NO'
vi$diffexpressed[vi$log2FoldChange>1 & vi$pvalue<0.05]="UP"
vi$diffexpressed[vi$log2FoldChange< -1 & vi$pvalue<0.05]="DOWN"
vi$delabel=NA
View(vi)
write.csv(vi,"Overall expression values_DEG combined.csv")
vi3=read.csv("Overall expression values_DEG combined.csv")
stat1_row <- vi3[vi3$Symbol == " Stat1 ", ]
gbp5_row<-vi3[vi3$Symbol == " Gbp5 ", ]
Fpr1_row=vi3[vi3$Symbol == " Fpr1 ", ]
library("ggrepel")
ggplot(data=vi3,aes(x=log2FoldChange,y=-log10(pvalue),col=diffexpressed,label=delabel))+
  geom_point(size=1)+
  theme_minimal()+
  scale_color_manual(values = c("red","black","blue"))+theme(text=element_text(size=15))+
  #geom_text_repel(data = gbp5_row, aes(label = Symbol),color="black", nudge_x = 0.1, nudge_y = 0.1, size = 4)+
  #geom_text_repel(data = stat1_row, aes(label = Symbol),color="black", nudge_x = 0.1, nudge_y = 0.1, size = 4)+
  #geom_text_repel(data = Fpr1_row, aes(label = Symbol),color="black", nudge_x = 0.1, nudge_y = 0.1, size = 4)+
  theme_classic()

ggplot(data = vi3, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label = delabel)) +
  geom_point(size = 1) +
  theme_minimal() +
  scale_color_manual(values = c("red", "black", "blue")) +
  theme(text = element_text(size = 15)) +
  geom_text_repel(data = stat1_row, aes(label = Symbol), color = "black", nudge_x = 0.1, nudge_y = 0.1, size = 4) +
  geom_text_repel(data = gbp5_row, aes(label = Symbol), color = "black", nudge_x = 0.1, nudge_y = 0.1, size = 4) +
  geom_text_repel(data = Fpr1_row, aes(label = Symbol), color = "black", nudge_x = 0.1, nudge_y = 0.1, size = 4) +
  theme_classic()


##MA plot
results= vi[,c(5,10,11)]
head(results)
results= results %>% mutate(new_p= ifelse(pvalue<0.05, TRUE, FALSE))
results$pvalue= NULL
str(results)
class(results)
#colnames(results)=c("baseMean","log2FoldChange","padj")
#results=as.data.frame(results)
write.csv(results, "results.csv")


##MA PLOT
##Gene plotter function was not working
colors <- ifelse(results$new_p, "blue", "gray")
plot(results$AveExpr, results$log2FoldChange, col = colors, pch = 19,cex=0.45,
     ylim = c(-5, 5), xlab = "Mean Expression", ylab = "Log Fold Change", main = "Murine_16_weeks_VL_infected_vs_Naive", abline(h = 0, col = "gray", lwd = 3))




