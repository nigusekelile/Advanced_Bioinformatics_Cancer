# this is Project 1 on Cancer biomarker discovery

# Project: Biomarker Discovery for Early Cancer Detection
# Cancer Subtype: Pancreatic Adenocarcinoma (PAAD)

# setting working directory
getwd()
setwd("C:/Users/nigus/Documents/HackBio")
getwd()

# enable libraries
library("TCGAbiolinks")
library(SummarizedExperiment) #install alongside with TCGAbiolinks
library(biomaRt)
library("edgeR")
library("limma")
library("EDASeq")
library(gplots)

# project information
getProjectSummary("TCGA-PAAD")
?GDCquery

# download and reprocess data
paadQ <- GDCquery(project = "TCGA-PAAD",
                 data.category = "Transcriptome Profiling",
                 data.type = "Gene Expression Quantification")

GDCdownload(paadQ)
paad.data <- GDCprepare(paadQ)
head(paad.data)

# Explore some metadata information
paad.data$race
paad.data$barcode
paad.data$tumor_descriptor
paad.data$ajcc_pathologic_stage


# we can create a simple metadata for our use case
simpleMeta <- data.frame("barcode" = paad.data$barcode,
                         "race" = paad.data$race,
                         "tumor_type" = paad.data$ajcc_pathologic_stage)

# Preprocessing, normalization and Filtration

# select unstranded dataset
paad.raw.data <- assays(paad.data)
dim(paad.raw.data$unstranded) #data size

# let's downsize our data to 5 primary and 5 recurring data alone
selectedBarcodes <- c(subset(simpleMeta, is.na(tumor_type) )$barcode[c(1:3)],
                      subset(simpleMeta, tumor_type == "Stage IA")$barcode[c(1:5)])

selectedData <- paad.raw.data$unstranded[, c(selectedBarcodes)]


dim(selectedData)



# Data normalization and filtering
normData <- TCGAanalyze_Normalization(tabDF = selectedData, 
                                      geneInfo = geneInfoHT, 
                                      method = "geneLength")
dim(normData)

# then filter 
filtData <- TCGAanalyze_Filtering(tabDF = normData, 
                                  method = "quantile",
                                  qnt.cut = 0.25)
dim(filtData)


# differential expression analysis
selectResults <- 
  TCGAanalyze_DEA(mat1 = filtData[, c(selectedBarcodes)[1:3]],
                  mat2 = filtData[, c(selectedBarcodes)[4:8]],
                  Cond1type = "Missing Tumor Type",
                  Cond2type = "Stage IA",
                  pipeline = 'edgeR', # it can be "limma"
                  fdr.cut = 0.01,
                  logFC.cut = 2)

# Differential expression analysis with Treatment levels

selectResults.level <- 
  TCGAanalyze_LevelTab(selectResults, "Missing Tumor Type", "Stage IA",
                       filtData[, c(selectedBarcodes)[1:3]],
                       filtData[, c(selectedBarcodes)[4:8]])


# now we can visualize with a heatmap
head(selectResults.level)
dim(selectResults.level)

heat.data <- filtData[rownames(selectResults.level),]

# color the plot by the kind of tumor
cancer.type <- c(rep("Missing Tumor Type", 5), rep("Stage IA", 5))

ccodes <- c()

for (i in cancer.type) {
  if (i == "Missing Tumor Type") {
    ccodes <- c(ccodes, "red")
  }else{
    ccodes <- c(ccodes, "blue")
  }
}

# Generate the heatmap
heatmap.2(x = as.matrix(heat.data),
          col = hcl.colors(10, palette = 'Blue-Red 2'), # search hcl colors in R
          Rowv = FALSE, Colv = TRUE,
          scale = 'row',
          sepcolor = 'black',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.5, cexCol = 1,
          main = "Missing Tumor Type vs Stage IA",
          na.color = 'black',
          ColSideColors = ccodes)


# Extract data for the two conditions
condition1_data <- filtData[, cancer.type == "Missing Tumor Type"]
condition2_data <- filtData[, cancer.type == "Stage IA"]

# Normalize data if required (e.g., log transformation)
condition1_data <- log2(condition1_data + 1)
condition2_data <- log2(condition2_data + 1)

# Function to calculate correlation matrix
calculate_coexpression <- function(data) {
  cor(data, method = "pearson")
}

# Calculate coexpression for both conditions
coexpression_condition1 <- calculate_coexpression(condition1_data)
coexpression_condition2 <- calculate_coexpression(condition2_data)

# Difference in coexpression matrices
coexpression_difference <- coexpression_condition2 - coexpression_condition1

# Identify biomarkers (genes with significant changes in coexpression)
# This could involve thresholding or statistical tests depending on your analysis
# Plot coexpression matrices
heatmap.2(coexpression_condition1,
          main = "Coexpression Matrix - Missing Tumor Type",
          col = hcl.colors(10, palette = 'Blue-Red 2'))

heatmap.2(coexpression_condition2,
          main = "Coexpression Matrix - Stage IA",
          col = hcl.colors(10, palette = 'Blue-Red 2'))

# Plot difference in coexpression
heatmap.2(coexpression_difference,
          main = "Difference in Coexpression",
          col = hcl.colors(10, palette = 'Blue-Red 2'))

# Example: Identify genes with high coexpression changes
threshold <- 0.03
biomarkers <- which(abs(coexpression_difference) > threshold, arr.ind = TRUE)

# Print potential biomarkers
print(biomarkers)


