library("DESeq2")
library(dplyr)
library(EnhancedVolcano)

metadata<-read.csv(file = "../../../data/processed_RNA-seq_data/metadata.tsv", sep= "\t", row.names = 1,check.names = FALSE)
countdata<-read.csv(file = "../../../data/processed_RNA-seq_data/FeatureCounts.tsv", sep= "\t", row.names = 1,check.names = FALSE)

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = metadata,
                              design= ~ oscillation)

dds <- DESeq(dds)

# SO_osc_vs_Control

res <- results(dds, contrast=c("oscillation","S+O","none"))
write.csv(res,"../../../data/processed_RNA-seq_data/DESeq_res_O2_sd_SO_osc_vs_Control.csv")

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 10e-5,
                FCcutoff = 1.0
)

# O_osc_vs_Control

res <- results(dds, contrast=c("oscillation","O","none"))
write.csv(res,"../../../data/processed_RNA-seq_data/DESeq_res_O2_sd_O_osc_vs_Control.csv")

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 10e-5,
                FCcutoff = 1.0
)

# S_osc_vs_Control

res <- results(dds, contrast=c("oscillation","S","none"))
write.csv(res,"../../../data/processed_RNA-seq_data/DESeq_res_O2_sd_S_osc_vs_Control.csv")

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 10e-5,
                FCcutoff = 1.0
)

### Do only with new metadata

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = metadata,
                              design= ~ condition)

dds <- DESeq(dds)

# S_famine_vs_Control

res <- results(dds, contrast=c("condition","WT_fed-batch_S_oscillation_famine","WT_fed-batch_control"))
write.csv(res,"../../../data/processed_RNA-seq_data/DESeq_res_O2_sd_S_famine_vs_Control.csv")

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 10e-5,
                FCcutoff = 1.0
)

# S_feast_vs_Control

res <- results(dds, contrast=c("condition","WT_fed-batch_S_oscillation_feast","WT_fed-batch_control"))
write.csv(res,"../../../data/processed_RNA-seq_data/DESeq_res_O2_sd_S_feast_vs_Control.csv")

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 10e-5,
                FCcutoff = 1.0
)


# S_osc_famine_vs_feast

res <- results(dds, contrast=c("condition","WT_fed-batch_S_oscillation_famine","WT_fed-batch_S_oscillation_feast"))
write.csv(res,"../../../data/processed_RNA-seq_data/DESeq_res_O2_sd_S_osc_famine_vs_feast.csv")

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 10e-5,
                FCcutoff = 1.0
)
