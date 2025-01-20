# This is Project 2 on Enrichment Analysis

# Project: Pathway Enrichment Analysis of Glioblastoma
# Cancer Subtype: Glioblastoma Multiforme (GBM)

# Setting working directory
getwd()
setwd("D:/HackBio/Advanced_Bioinformatics_Cancer")
getwd()

# Enable libraries
library("TCGAbiolinks")
library(SummarizedExperiment) # Install alongside TCGAbiolinks
library(biomaRt)
library("edgeR")
library("limma")
library("EDASeq")
library(gplots)

# Project Information
gb_proj <- getProjectSummary("TCGA-GBM")
print(gb_proj)

# Step 1: Data Acquisition
# Download and preprocess RNA-Seq data for GBM from TCGA
gbmQ <- GDCquery(project = "TCGA-GBM",
                 data.category = "Transcriptome Profiling",
                 data.type = "Gene Expression Quantification")

GDCdownload(gbmQ)
gbm.data <- GDCprepare(gbmQ)

# Explore metadata
gbm.meta <- colData(gbm.data)
head(gbm.meta)

# Create a simplified metadata dataframe
simpleMeta <- data.frame("barcode" = gbm.meta$barcode,
                         "race" = gbm.meta$race,
                         "tumor_type" = gbm.meta$tumor_descriptor)

# Step 2: Data Preprocessing
# Select unstranded dataset and downsize to 5 primary and 5 recurrent samples
gbm.raw.data <- assays(gbm.data)$unstranded
selectedBarcodes <- c(subset(simpleMeta, tumor_type == "Recurrence")$barcode[1:5],
                      subset(simpleMeta, tumor_type == "Primary")$barcode[1:5])
selectedData <- gbm.raw.data[, selectedBarcodes]

# Normalize and filter the data
normData <- TCGAanalyze_Normalization(tabDF = selectedData, 
                                      geneInfo = geneInfoHT, 
                                      method = "geneLength")
filtData <- TCGAanalyze_Filtering(tabDF = normData, 
                                  method = "quantile",
                                  qnt.cut = 0.25)

# Step 3: Differential Expression Analysis
selectResults <- TCGAanalyze_DEA(mat1 = filtData[, 1:5],
                                 mat2 = filtData[, 6:10],
                                 Cond1type = "Recurrence",
                                 Cond2type = "Primary",
                                 pipeline = 'edgeR',
                                 fdr.cut = 0.01,
                                 logFC.cut = 2)

selectResults.level <- TCGAanalyze_LevelTab(selectResults, "Recurrence", "Primary",
                                            filtData[, 1:5], filtData[, 6:10])

# Step 4: Pathway Enrichment Analysis
# Convert Ensembl IDs to gene IDs
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

upreg.genes <- rownames(subset(selectResults.level, logFC > 2))
dnreg.genes <- rownames(subset(selectResults.level, logFC < -2))

upreg.genes <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                     filters = 'ensembl_gene_id',
                     values = upreg.genes, mart = mart)$hgnc_symbol
dnreg.genes <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                     filters = 'ensembl_gene_id',
                     values = dnreg.genes, mart = mart)$hgnc_symbol

up.EA <- TCGAanalyze_EAcomplete(TFname = "Upregulated", upreg.genes)
dn.EA <- TCGAanalyze_EAcomplete(TFname = "Downregulated", dnreg.genes)

# Step 5: Visualization
# Heatmap of Differentially Expressed Genes
heat.data <- filtData[rownames(selectResults.level),]

cancer.type <- c(rep("Recurrence", 5), rep("Primary", 5))
ccodes <- ifelse(cancer.type == "Recurrence", "red", "blue")

heatmap.2(x = as.matrix(heat.data),
          col = hcl.colors(10, palette = 'Blue-Red 2'),
          Rowv = F, Colv = T,
          scale = 'row',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.5, cexCol = 1,
          main = "Heatmap of Recurrence and Primary",
          na.color = 'black',
          ColSideColors = ccodes)

# Volcano Plot
plot(x = selectResults.level$logFC, y = -log10(selectResults.level$FDR),
     pch = 20, col = ifelse(selectResults.level$logFC > 2, "red",
                            ifelse(selectResults.level$logFC < -2, "blue", "black")),
     xlab = "Log Fold Change", ylab = "-log10(FDR)",
     main = "Volcano Plot of Differential Expression")

# Pathway Enrichment Barplots
TCGAvisualize_EAbarplot(tf = rownames(up.EA$ResBP),
                        GOBPTab = up.EA$ResBP,
                        GOCCTab = up.EA$ResCC,
                        GOMFTab = up.EA$ResMF,
                        PathTab = up.EA$ResPat,
                        nRGTab = upreg.genes,
                        nBar = 10,
                        text.size = 2,
                        fig.width = 30,
                        fig.height = 15)

TCGAvisualize_EAbarplot(tf = rownames(dn.EA$ResBP),
                        GOBPTab = dn.EA$ResBP,
                        GOCCTab = dn.EA$ResCC,
                        GOMFTab = dn.EA$ResMF,
                        PathTab = dn.EA$ResPat,
                        nRGTab = dnreg.genes,
                        nBar = 10,
                        text.size = 2,
                        fig.width = 30,
                        fig.height = 15)

