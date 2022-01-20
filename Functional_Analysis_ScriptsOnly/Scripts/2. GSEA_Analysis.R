# Load libraries
library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(tidyverse)
library(stringr)
library(dplyr)

set.seed(888888)

## Note: replace diversity order, e.g. Q2.00, with desired diversity order ##
## throughout this R code                                                  ##

data <- read.csv("Output/DiffExp_age_plus_div_Q3.00.csv")

colnames(data)[1] <- "gene"

ortho <- read.csv("Raw_data/killifish_human_ortho.csv")

## mapping killifish transcripts to human orthologues ## 

data$gene <- plyr::mapvalues(data$gene, from=c(ortho$Gene.stable.ID), to=c(ortho$Gene.stable.ID.1))

data$gene[!startsWith(data$gene , "ENSG0")] <- NA
data <- na.omit(data)
old_data <- data

## Summarise those transcripts that have duplicate matches ## 
# Some genes are duplicated, due to many:one mapping. Take mean of these and remove duplicates #

data$FC <- 2^data$log2FoldChange

data <- data %>% group_by(gene) %>% summarise_all(funs(mean))

data <- as.data.frame(data)

data$log2FoldChange <- log2(data$FC)

#data[-1] <- sapply(data[-1], as.integer)

##foldchanges for GSEA analysis##

foldchanges <- data$log2FoldChange

names(foldchanges) <- data$gene

## Sort fold changes in decreasing order
foldchanges <- sort(foldchanges, decreasing = TRUE)

head(foldchanges)

# GSEA using gene sets associated with BP Gene Ontology terms
gseaGO <-gseGO(
  foldchanges,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  exponent = 1,
  minGSSize = 10,
  eps = 1e-10,
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  verbose = TRUE,
  seed = TRUE,
  by = "fgsea",
)

gseaGO_results <- gseaGO@result
save(gseaGO, file = "Output/GSEA_analysis/GSEA_Q3.00.rda")

gsea_simple <- clusterProfiler::simplify(
  gseaGO,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang"
)

gsea_simple_results <- gsea_simple@result
write.csv(gsea_simple_results, "Output/GSEA_analysis/GSEA_Q3.00.csv")