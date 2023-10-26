GSEA_analysis_hallmarks.R:
#In this script I do gene set enrichment analysis
#Prior to analysis, I filtered the RNAseq data with p-adj<0.05, and NO filtering of log2FC
#I compare my filtered data to downloaded hallmark gene sets
#Results saved as SRP123625_gsea_hallmarks_results.tsv

GSEA_analysis_GO.R:
#In this script I do gene set enrichment analysis
#Prior to analysis, I filtered the RNAseq data with p-adj<0.05, and NO filtering of log2FC
#I compare my filtered data to downloaded ontology gene sets
#Results saved as SRP123625_gsea_GO_results.tsv

GSEA_analysis_curated.R:
#In this script I do gene set enrichment analysis
#Prior to analysis, I filtered the RNAseq data with p-adj<0.05, and NO filtering of log2FC
#I compare my filtered data to downloaded curated gene sets
#Results saved as SRP123625_gsea_curated_results.tsv

SRP123625_gsea_GO_dotplot_PIK3_regulated.png:
#dotplot demonstrating deregulated pathways affected by PIK3

Excel RNAseq file under "Data" folder:
#Contains 2 sheets of RNAseq data, one un-filtered, and one filtered by p-adjusted value

CSV RNAseq file under "Data" folder:
#Comma-delimited file of RNAseq data, not filtered by p-value