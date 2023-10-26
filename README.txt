RNAseq: project folder
Created: 20.10.2023
User: Maren Austevik Bø

- Data: folder containing data
- Scripts: folder containing scripts
- Files: folder containing files

GSEA_analysis_hallmarks.R:
#In this script I do gene set enrichment analysis
#I compare my filtered data to downloaded hallmark gene sets
#Results saved as SRP123625_gsea_hallmarks_results.tsv

GSEA_analysis_GO.R:
#In this script I do gene set enrichment analysis
#I compare my filtered data to downloaded ontology gene sets
#Results saved as SRP123625_gsea_GO_results.tsv

GSEA_analysis_GO_small_geneSetSize.R:
#In this script I do gene set enrichment analysis, filtering for small gene set size
#I compare my filtered data to downloaded ontology gene sets
#Results saved as SRP123625_gsea_GO_small_geneSetSize_results.tsv

GSEA_analysis_curated.R:
#In this script I do gene set enrichment analysis
#I compare my filtered data to downloaded curated gene sets
#Results saved as SRP123625_gsea_curated_results.tsv

GSEA_analysis_curated_small_geneSetSize.R:
##In this script I do gene set enrichment analysis, filtering for small gene set size
#I compare my filtered data to downloaded curated gene sets
#Results saved as SRP123625_gsea_curated_small_geneSetSize_results.tsv

GSEA_analysis_immunologic.R:
#In this script I do gene set enrichment analysis
#I compare my filtered data to downloaded immunologic gene sets
#Results saved as SRP123625_gsea_immunologic_results.tsv

GSEA_analysis_oncogenic.R:
#In this script I do gene set enrichment analysis
#I compare my filtered data to downloaded oncogenic gene sets
#Results saved as SRP123625_gsea_C6_results.tsv

GSEA_analysis_C4.R:
#In this script I do gene set enrichment analysis
#I compare my filtered data to downloaded C4 gene sets
#Results saved as SRP123625_gsea_C4_results.tsv

GO_term_analysis.R
#In this script I do gene ontology term analysis
#Prior to analysis, I filtered the RNAseq data with p-adj<0.05

regulated_GO_terms text document:
#containing description of all regulated pathways discovered with the GO analysis in script GO_term_analysis.R

SRP123625_gsea_GO_dotplot_PIK3_regulated.png:
#dotplot from GO gsea demonstrating deregulated pathways affected by PIK3

SRP123625_gsea_GO_dotplot_PIK3_regulated.png:
#Dotplot from GO gsea demonstrating deregulated gene sets containing PIK3

SRP123625_GO_term_Analysis_Lipid_Fat_Cholesterol.png:
#Dotplot from GO analysis showing deregulated fat- lipid- and cholesterol-pathways

SRP123625_GO_term_lipid_fat_cholesterol_rigdeplot.png:
#Rigdeplot from GO analysis showing deregulated fat- lipid- and cholesterol-pathways

SRP123625_GO_term_stem_cell.png:
#Dotplot from GO analysis showing deregulated stem cell pathways

SRP123625_GO_term_stem_cell_rigdeplot.png:
#Rigdeplot from GO analysis showing deregulated stem cell pathways

SRP123625_gsea_C6_dotplot_PTEN_regulated.png:
#Dotplot from oncogenic gsea showing deregulated gene sets assosiated with knockdown of PTEN

SRP123625_gsea_curated_dotplot_PIP5K1B_regulated_stemCellGeneSets.png:
#Dotplot from curated gsea demonstrating deregulated gene sets assosiated with stem cell properties, all of them including PIP5K1B

SRP123625_gsea_enrich_cnetplot.png:
#cnetplot from hallmark gsea analysis
#just a test of creating a cnetplot

SRP123625_gsea_enrich_negative_plot.png:
#gsea plot from hallmark gsea, demonstrating hallmark with most negative NES value
#just a test of creating a gsea plot

SRP123625_gsea_enrich_positive_plot.png:
#gsea plot from hallmark gsea, demonstrating hallmark with most positive NES value
#just a test of creating a gsea plot

Excel RNAseq file under "Data" folder:
#Contains 2 sheets of RNAseq data, one un-filtered, and one filtered by p-adjusted value

CSV RNAseq file under "Data" folder:
#Comma-delimited file of RNAseq data, not filtered by p-value


