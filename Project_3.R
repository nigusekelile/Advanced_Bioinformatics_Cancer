# this is Project 3 on cancer Drug Design

# Project: Population Sensitive Drug Design for Cancer Treatment
# Cancer Subtype: Ovarian Serous Cystadenocarcinoma (OV)

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

# project infomation
getProjectSummary("TCGA-OV")
?GDCquery

# download and reprocess data
ovQ <- GDCquery(project = "TCGA-OV",
                 data.category = "Transcriptome Profiling",
                 data.type = "Gene Expression Quantification")

GDCdownload(ovQ)
ov.data <- GDCprepare(ovQ)
head(ov.data)

# Explore some metadata information
ov.data$race
ov.data$tumor_descriptor
ov.data$barcode

# we can create a simple metadata for our use case
simpleMeta <- data.frame("barcode" = ov.data$barcode,
                         "race" = ov.data$race,
                         "tumor_type" = ov.data$tumor_descriptor)

