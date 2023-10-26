#Loading and installing packages
library(tidyverse)
library(here)
library(dplyr)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(ggplot2)
library(magrittr)

BiocManager::install("enrichplot")
library(enrichplot)
#Do update
a

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

#Perform Go term analysis
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

#Extracting column with GO descriptions
GO_descriptions <- gse$Description
# Convert list to a data frame (if needed)
GO_descriptions_df <- as.data.frame(GO_descriptions)
# Save the data frame to a text file
write.table(GO_descriptions_df, file = "regulated_GO_terms.txt")

#Selecting results related with GO-term including fats, lipids or cholesterol
Selected_GO_terms <- gse %>% filter(str_detect(Description, "lipid|PIP3|fat|fatty|cholesterol"))
# Example: View only the "column1" and "column2" from your dataset
subset_df <- Selected_GO_terms[, c("Description")]
subset_df

#removing GO terms containing word "fate" or "Fate"
Selected_GO_terms_filtered <- Selected_GO_terms %>% filter(!str_detect(Description, "fate"))
subset_df_filtered <- Selected_GO_terms_filtered[, c("Description")]
subset_df_filtered

#Creating dotplot
require(DOSE)
dotplot <- dotplot(Selected_GO_terms_filtered, showCategory=10, split=".sign") + facet_grid(.~.sign)
dotplot

# Adjust the size of the y-axis text
dotplot <- dotplot + theme(axis.text.y = element_text(size = 6))
dotplot

#Saving dotplot
ggplot2::ggsave(file.path("Plots", "SRP123625_GO_term_Analysis_Lipid_Fat_Cholesterol.png"),
                plot = dotplot)

#Making rigdeplot
ridgeplot_2 <- ridgeplot(Selected_GO_terms_filtered) + labs(x = "enrichment distribution")
ridgeplot_2
# Adjust the size of the y-axis text
ridgeplot_2 <- ridgeplot_2 + theme(axis.text.y = element_text(size = 7))
ridgeplot_2

#Saving plot
ggplot2::ggsave(file.path("Plots", "SRP123625_GO_term_lipid_fat_cholesterol_rigdeplot.png"),
                plot = ridgeplot_2)

#Making a category net plot
#categorySize can be either 'pvalue' or 'geneNum'
cnetplot(Selected_GO_terms_filtered, categorySize="pvalue", foldChange=gene_list, showCategory = 3)

#Selecting results related with GO-term including insulin or PIP3
Selected_GO_terms_2 <- gse %>% filter(str_detect(Description, "PIP3|insulin"))
dotplot2 <- dotplot(Selected_GO_terms_2, showCategory=10, split=".sign") + facet_grid(.~.sign)
dotplot2
#Saving dotplot
ggplot2::ggsave(file.path("Plots", "SRP123625_GO_term_insulin.png"),
                plot = dotplot2)


#Selecting results related with GO-term including stem cell
Selected_GO_terms_3 <- gse %>% filter(str_detect(Description, "stem cell"))
#making dotplot
dotplot3 <- dotplot(Selected_GO_terms_3, showCategory=10, split=".sign") + facet_grid(.~.sign)
dotplot3
#Saving dotplot
ggplot2::ggsave(file.path("Plots", "SRP123625_GO_term_stem_cell.png"),
                plot = dotplot3)

#Making rigdeplot
ridgeplot <- ridgeplot(Selected_GO_terms_3) + labs(x = "enrichment distribution")
#Saving plot
ggplot2::ggsave(file.path("Plots", "SRP123625_GO_term_stem_cell_rigdeplot.png"),
                plot = ridgeplot)