library(PCAtools)
library(dplyr)


# PCA of A matrix

logtpmData<-read.csv(file="../../../data/processed_RNA-seq_data/log_tpm.csv",row.names = 1, check.names = FALSE)
MetaData<-read.csv(file = "../../../data/processed_RNA-seq_data/metadata.tsv", row.names = 1, sep = "\t", check.names = FALSE)

logtpmData<- logtpmData %>%mutate_all(as.numeric)
logtpmData <- t(scale(t(logtpmData)))

p <- pca(logtpmData,metadata = MetaData, removeVar = 0.1)

biplot(p,
       x = 'PC1',
       y = 'PC4',
       ntopLoadings = 5,fillBoxedLoadings = "white",widthLoadingsArrows = 1.0,drawConnectorsLoadings = TRUE,
       sizeLoadingsNames = 5.0,
       borderWidth = 1.0,gridlines.minor = FALSE,gridlines.major = FALSE,vlineWidth = 1,hlineWidth = 1,
       colby = "oscillation",
       shape = "sample_ID",
       colLegendTitle = "Oscillation",
       shapeLegendTitle = "Timepoint",
       hline = 0,
       vline = 0,
       legendPosition = "right",borderColour = "black",
       lab = NULL,
       showLoadings = TRUE)

pairsplot(p,components = c("PC1","PC2","PC3","PC4"))

write.csv(p["loadings"],"../../../data/processed_RNA-seq_data/log_tpm_PCA_loadings.csv")
write.csv(p["rotated"],"../../../data/processed_RNA-seq_data/log_tpm_PCA_rotated.csv")
write.csv(p["metadata"],"../../../data/processed_RNA-seq_data/log_tpm_PCA_metadata.csv")
write.csv(p["variance"],"../../../data/processed_RNA-seq_data/log_tpm_PCA_variances.csv")

screeplot(p)