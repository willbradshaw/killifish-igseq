## Analysis using DEseq2, controlling for AB diversity ##

library("DESeq2")
library("tidyverse")

## Loading Data ##

data <- read.csv("Raw_data/rna_data_raw_counts.csv", sep=',')
meta_data <- read.csv("Raw_data/meta_data_with_Qvalues.csv")

## Running Analysis ##

## Note: replace diversity order e.g. 'Q4.00' with desired diversity order ##
## in DESeqDataSetFromMatrix function, and name output file suitably ##

row.names(data) <- data$ensgene
data <- data[-1]

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = meta_data,
                              design = ~age +Q4.00)

dds <- estimateSizeFactors(dds)
normalised <- counts(dds, normalized=TRUE)
write.csv(normalised, "Output/Normalised_RNAseq_Counts.csv")

dds <- DESeq(dds)
#contrast <- c('age', 'young', 'old')
#results <- results(dds, contrast = contrast)
results <- results(dds)
#results %>% data.frame() %>% View()
coeffs <- coef(dds)

write.csv(results, "Output/DiffExp_age_plus_div_Q4.00.csv")
write.csv(coeffs, "Output/DiffExp_age_plus_div_Q4.00_coefficients.csv")

