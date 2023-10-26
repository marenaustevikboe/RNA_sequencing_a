#Loading and installing packages
library(tidyverse)
library(here)
library(dplyr)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(ggplot2)
library(magrittr)
library(enrichplot)

BiocManager::install("pathview")
library(pathview)
#Do not update
n

# Reading dataset
RNAseq_raw <- read_csv2(here("Data", "RNAseq.csv"))
glimpse(RNAseq_raw)

# Filter rows with p-value less than 0.05 (significant)
RNAseq_significant <- filter(RNAseq_raw, `P-adj` < 0.05)
glimpse(RNAseq_significant)

#Tidying colums
dge_df<- RNAseq_significant[, c("GeneID",	"Gene", "Base mean", 
                                "log2(FC)", "StdErr", "Wald-Stats",
                                "P-value",	"P-adj")]

# we want the log2 fold change 
original_gene_list <- dge_df$"log2(FC)"

# name the vector
names(original_gene_list) <- dge_df$"GeneID"

# omit any NA values 
gene_list<-na.omit(original_gene_list)
glimpse(gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db", 
             pAdjustMethod = "none")

glimpse(gse)

require(DOSE)
dotplot <- dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
# Adjust the size of the y-axis text
dotplot <- dotplot + theme(axis.text.y = element_text(size = 5))
dotplot
