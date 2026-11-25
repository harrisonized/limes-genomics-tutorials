## LIMES Genomics Tutorials

This follows the R tutorials for the "Genomics and Epigenetics" course.

## Installation

Install the following packages in R:

```R
install.packages("reshape2")
install.packages("pheatmap")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("ggpubr")
install.packages("DESeq2")
install.packages("BiocManager")
BiocManager::install(c(
    "DESeq2",
    "biomaRt",
    "org.Hs.eg.db", 
    "clusterProfiler")
)

```

### How to Run

After setting up the files to be analyzed, run your desired script from the command line. For example:

```bash
cd path/to/limes-genomics-tutorials
Rscript R/bulkrnaseq.R
```
