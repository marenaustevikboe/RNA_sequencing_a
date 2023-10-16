library(tidyverse)
library(here)

# Reading data
RNAseq = read_csv2(here("Data", "RNAseq.csv"))
glimpse(RNAseq)

# Load necessary libraries 
library(dplyr)

# Filter rows with p-value less than 0.05 (significant)
significant_rows <- RNAseq %>%
  filter("P-value" < 0.05)
glimpse(significant_rows)

#Filter great log2FC values
tidyData <- significant_rows %>%
  filter("Log2(FC)"<-1) %>%
  filter("Log2(FC)">1)

#Choosing colums I want to work with
extracted_columns <- tidyData[, c("GeneID", "log2(FC)")]
glimpse(extracted_columns)

