#loading libraries
library(tidyverse)
library(here)

# Reading dataset
RNAseq = read_csv2(here("Data", "RNAseq.csv"))
glimpse(RNAseq)

# Load necessary libraries 
library(dplyr)

# Filter rows with p-value less than 0.05 (significant)
significant_rows <- filter(RNAseq, `P-value` < 0.05)
glimpse(significant_rows)

#Filter log2FC values
tidyData<-significant_rows %>%
  filter(`log2(FC)`<(-1)|`log2(FC)`>1)
glimpse(tidyData)

#Choosing colums I want to work with
extracted_columns <- tidyData[, c("GeneID", "log2(FC)")]
glimpse(extracted_columns)

library(ggplot2)





