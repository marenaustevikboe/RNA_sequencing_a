#In this script I do gene set enrichment analysis
#Prior to analysis, I filtered the RNAseq data with p-adj<0.05, and used NO filter on log2FC
#I compare my filtered data to downloaded hallmark gene sets

#loading libraries
library(tidyverse)
library(here)
library(dplyr)

# Reading dataset
RNAseq_raw <- read_csv2(here("Data", "RNAseq.csv"))
glimpse(RNAseq_raw)

# Filter rows with p-value less than 0.05 (significant)
RNAseq_significant <- filter(RNAseq_raw, `P-adj` < 0.05)
glimpse(RNAseq_significant)

#Choosing colums I want to work with
dge_df<- RNAseq_significant[, c("GeneID",	"Gene", "Base mean", 
                         "log2(FC)", "StdErr", "Wald-Stats",
                         "P-value",	"P-adj")]

#installing packages and loading libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!("clusterProfiler" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("clusterProfiler", update = FALSE)
}
library(clusterProfiler)

if (!("msigdbr" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("msigdbr", update = FALSE)
}
library(msigdbr)

if (!("org.Hs.eg.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}
library(org.Hs.eg.db)

# Attach the DESeq2 library
if (!require("BiocManager", quietly = TRUE))
  BiocManager::install("DESeq2")
library(DESeq2)

# Attach the ggplot2 library for plotting
library(ggplot2)

# We will need this so we can use the pipe: %>%
library(magrittr)

#Define gene set
mm_hallmark_sets <- msigdbr(
  species = "Homo sapiens", 
  category = "H")
head(mm_hallmark_sets)

# Check if duplicates
any(duplicated(dge_df$Gene))

dup_Genes <- dge_df %>%
  dplyr::filter(duplicated(Gene)) %>%
  dplyr::pull(Gene)

dge_df %>%
  dplyr::filter(Gene %in% dup_Genes) %>%
  dplyr::arrange(Gene)

#remove duplicates, keep rows with highest absolute value of Log2FC
filtered_dge_df <- dge_df %>%
  # Sort so that the highest absolute values of the log2 fold change are at the
  # top
  dplyr::arrange(dplyr::desc(abs(`log2(FC)`))) %>%
  # Filter out the duplicated rows using `dplyr::distinct()`
  dplyr::distinct(Gene, .keep_all = TRUE)

#confirm that duplicates are removed
any(duplicated(filtered_dge_df$Gene))

# Let's create a named vector ranked based on the log2 fold change values
lfc_vector <- filtered_dge_df$`log2(FC)`
names(lfc_vector) <- filtered_dge_df$Gene

# We need to sort the log2 fold change values in descending order here
lfc_vector <- sort(lfc_vector, decreasing = TRUE)

# Look at first entries of the ranked log2 fold change vector
head(lfc_vector)

# Set the seed so our results are reproducible:
set.seed(2020)

library(stats)

#Performing Gene set enrichment analysis
gsea_results <- GSEA(
  geneList = lfc_vector, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    mm_hallmark_sets,
    gs_name,
    gene_symbol
  ))

# We can access the results from our `gsea_results` object using `@result`
head(gsea_results@result)
glimpse(gsea_results@result)

#Save results as dataframe
gsea_result_df <- data.frame(gsea_results@result)
glimpse(gsea_result_df)

#Write results to file
readr::write_tsv(
  gsea_result_df,
  file.path(
    "Results",
    "SRP123625_gsea_results_hallmark.tsv"
  )
)

#Making a plot with the most positive NES value
gsea_result_df %>%
  # This returns the 3 rows with the largest NES values
  dplyr::slice_max(NES, n = 2)

most_positive_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "HALLMARK_ALLOGRAFT_REJECTION",
  title = "HALLMARK_ALLOGRAFT_REJECTION",
  color.line = "#0d76ff"
)
most_positive_nes_plot

ggplot2::ggsave(file.path("Plots", "SRP123625_gsea_enrich_positive_plot.png"),
                plot = most_positive_nes_plot
)

#Making a plot with the most negative NES value
gsea_result_df %>%
  # Return the 3 rows with the smallest (most negative) NES values
  dplyr::slice_min(NES, n = 3)

most_negative_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  title = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  color.line = "#0d76ff"
)
most_negative_nes_plot

ggplot2::ggsave(file.path("Plots", "SRP123625_gsea_enrich_negative_plot.png"),
                plot = most_negative_nes_plot
)

#Making dotplot
require(DOSE)
dotplot <- dotplot(gsea_results, showCategory=10, split=".sign") + facet_grid(.~.sign)
dotplot

ggplot2::ggsave(file.path("Plots", "SRP123625_gsea_enrich_dotplot.png"),
                plot =dotplot
)

#Making a category net plot
# categorySize can be either 'pvalue' or 'geneNum'
cnetplot <- cnetplot(gsea_results, categorySize="pvalue", color.params = list(foldChange=dge_df$"log2(FC)"), showCategory = 3)
ggplot2::ggsave(file.path("Plots", "SRP123625_gsea_enrich_cnetplot.png"),
                plot =cnetplot
)
cnetplot

#Making a rigdde plot
install.packages("ggridges")
library(ggridges)
ridgeplot <- ridgeplot(gsea_results) + labs(x = "enrichment distribution")
ggplot2::ggsave(file.path("Plots", "SRP123625_gsea_enrich_ridgeplot.png"),
                plot = ridgeplot
)
ridgeplot
