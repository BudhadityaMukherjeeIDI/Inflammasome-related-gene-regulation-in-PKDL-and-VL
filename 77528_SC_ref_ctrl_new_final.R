library(GEOquery)

my.gse<-c("GSE77528")

# create directory to store files downloaded from GEO and results anlyzed
dir.create(c("geo_downloads","results"))

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

# The phenotype data is contained inside pData slot. 
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

##SUBSETTING THE COLUMNS OF OUR INTEREST
my.pdata=select(pdata,"title","geo_accession","source_name_ch1","organism_ch1")
head(my.pdata)
my.pdata

##This is to make a function to group the samples
#alternatively can be done manually by Excel

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

## Grouping 
group_num<-c("111111112222222222222233333333333333344444444")
group_names<-paste("Human",c("VL_infected","asymptomatic","control", "treated"),sep="_")
# allocate groups
my.pdata[["group"]]<-group_by_numbering(group_num, pdata$title, group_names)
# num of replicates in each group
table(my.pdata$group)

# store processed phenotype data
write.table(pdata, "GSE77528_PhenoData.txt", 
            sep="\t", quote=F,col.names=NA)

assayData<-exprs(my.geo.gse)#exprs object in expression set contains the expression values wrt genes
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
write.csv(exprs(my.rma),"GSE77528_NORMEXPRDATA_new.csv")

library(limma)
library(ggrepel)
pca <- prcomp(t(exprs(my.rma)))
pca3=cbind(my.pdata, pca$x) #Binding the dataset of first two PCA components with phenotypic data since we need to color it wrt sample
ggplot(data=pca3,aes(x = PC1, y= PC2, col=group)) + geom_point()+theme_classic()

#Quantile values
as.numeric(quantile(exprs(my.rma), c(0., 0.025, 0.05, 0.075, 0.099, .10), na.rm=T))

# check if rows of phenotype is sorted same way as columns of exprs matrix
table(rownames(pData(my.rma))==colnames(exprs(my.rma)))

pData(my.rma)$group <- as.factor(pData(my.rma)$group)
category=pData(my.rma)$group
class(category)
category=relevel(category,ref = "Human_control")
my.design=model.matrix(~0+category,pData(my.rma))
rownames(my.design) <-pData(my.rma)$group
colnames(my.design)=levels(category)


# fit linear model
fit <- lmFit(my.rma, my.design)

# specify the comparison btw group of interest, Here when we are specifying 2 contrast how come the Volcano
#is generated for one
my.contrasts <- makeContrasts(VL_vs_control= Human_VL_infected-Human_control,
                               levels = colnames(my.design))
my.contrasts

write.csv(my.contrasts,"GSE77528_Comparison.csv")

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

##This is for annotation file,since fData is not working for this...
gpl <- getGEO('GPL10558', destdir=".")
colnames(Table(gpl))
Gene= Table(gpl)$ILMN_Gene
Probe_ID= Table(gpl)$ID
Description= Table(gpl)$Definition
e= cbind(Description, Probe_ID)
e= cbind(Gene, e)
#head(f, 12)
f= as.data.frame(e)
f= f[complete.cases(f[ , 1]),]
f= f[f$Gene!="", ]
e=f

View(contrast.tts)

vi= contrast.tts[[1]] ##Only taking the 1st index containing comparison between control vs VL_infected
vi= as.data.frame(vi)
Probe_ID= rownames(vi)
vi= cbind(Probe_ID, vi)
rownames(vi)= NULL
vi$log2FoldChange= vi$logFC
vi$pvalue= vi$adj.P.Val
vi= inner_join(e, vi, by= "Probe_ID")

##VOLCANOPLOT

vi$Diffexpressed='NO'
vi$Diffexpressed[vi$log2FoldChange>1 & vi$pvalue<0.05]="UP"
vi$Diffexpressed[vi$log2FoldChange< -1 & vi$pvalue<0.05]="DOWN"
vi$delabel=NA
View(vi)
table(vi$Diffexpressed)
write.csv(vi,"77528_FINAL3_EXCEL.csv")
vi3=read.csv("77528_FINAL3_EXCEL.csv")
library(ggplot2)
library(ggrepel)
##ORIGINAL Volcano axis
stat1_row <- vi3[vi3$Gene == "STAT1", ]
ggplot(data=vi3,aes(x=log2FoldChange,y=-log10(pvalue),col=Diffexpressed,label=delabel))+
  geom_point(size=1)+
  theme_minimal()+
  scale_color_manual(values = c("red","black","blue"))+
  #geom_text_repel(data = stat1_row, aes(label = Gene),color="black", nudge_x = 0.1, nudge_y = 0.1, size = 4)+
  theme(text=element_text(size=15))+
  theme_classic()

install.packages("ggalt")
library(ggalt)
#Volcano yaxis 20
ggplot(data=vi,aes(x=log2FoldChange,y=-log10(pvalue),col=Diffexpressed,label=delabel))+
  geom_point(size=1)+
  theme_minimal()+ylim(0,20)+
  scale_color_manual(values = c("red","black","green"))+
  theme(text=element_text(size=15))+theme_classic()
  



results= vi[,c(5,10,11)]
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
    ylim = c(-5, 5), xlab = "Mean Expression", ylab = "Log Fold Change", main = "Control_Vs_Human_VL_infected", abline(h = 0, col = "gray39", lwd = 3))



