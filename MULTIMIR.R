BiocManager::install("multiMiR")
library("multiMiR")
mRNA=read.csv("PKDL.csv")
head(mRNA)
mRNA=mRNA[,1]
mRNA=as.vector(mRNA)
mRNA=as.data.frame(mRNA)

mirna=read.csv("miRNA_PKDL.csv")
mirna=as.vector(mirna)
mirna=as.data.frame(mirna)

PKDL= get_multimir(org     = "hsa",
                         mirna   = mirna,
                         target  = mRNA,
                         table   = "predicted",
                         summary = TRUE,
                         predicted.cutoff.type = "p",
                         predicted.cutoff      = 10,
                         use.tibble = TRUE)

table(PKDL@data$type)

miRNAgenes=PKDL@data
write.csv(miRNAgenes,"PKDL_PREDICTED_TOP10.csv")

PKDL1= get_multimir(org     = "hsa",
                   mirna   = mirna,
                   target  = mRNA,
                   table   = "validated",
                   summary = TRUE,
                   use.tibble = TRUE)

table(PKDL1@data$type)

miRNAvalidated=PKDL1@data
write.csv(miRNAvalidated,"PKD_VALIDATED.csv")

##VL
vlmrna=read.csv("mrna_VL.csv")
head(vlmrna)
dim(vlmrna)
vlmrna=as.data.frame(vlmrna)

vlmiRNA=read.csv("miRNA_VL.csv")
vlmiRNA=as.data.frame(vlmiRNA)


VL= get_multimir(org     = "hsa",
                   mirna   = vlmiRNA,
                   target  = vlmrna,
                   table   = "predicted",
                   summary = TRUE,
                   predicted.cutoff.type = "p",
                   predicted.cutoff      = 5,
                   use.tibble = TRUE)
table(VL@data$type)
vlpredicted=VL@data
write.csv(vlpredicted, "VL_predicted_05.csv")

VL1= get_multimir(org     = "hsa",
                 mirna   = vlmiRNA,
                 target  = vlmrna,
                 table   = "validated",
                 summary = TRUE,
                 use.tibble = TRUE)
table(VL1@data$type)
vlvalidated=VL1@data
write.csv(vlvalidated, "VL_validated_NEW.csv")


