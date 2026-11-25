suppressPackageStartupMessages(library('dplyr'))
library("ggplot2")
suppressPackageStartupMessages(library("clusterProfiler"))

# Functions
# GeneExprBox
# plotVolcano
# GOEA

#' Gene Expression Box Plot'
#' 
GeneExprBox <- function(gene, norm, gene_anno) {

    # Create a pattern to find an exact match for the gene symbol
    exactgene <- paste("^", gene, "$", sep = "") 

    # Use the pattern to search for the gene symbol in the 'gene_anno' data frame
    ID <- gene_anno[grep(exactgene, gene_anno$symbol),]

    # Check if the gene is present, unique, or absent based on the search results
    if (nrow(ID) == 1) {
        "Gene present"
    } else {
        if (nrow(ID) > 1) {
            stop("ID not unique")
        } else {
            stop("Gene not present")
        }
    }

    # Subset 'norm' matrix for the identified gene's expression values
    GeneExpr <- norm[rownames(norm) %in% ID$ENTREZ,]

    # Create a data frame 'GenePlotTable' for plotting gene expression
    GenePlotTable <- as.data.frame(t(GeneExpr))
    colnames(GenePlotTable) <- "value"
    GenePlotTable$cond <- sampletable$cond

    # Create a boxplot to visualize gene expression across different conditions
    ggplot(GenePlotTable, aes(x = cond, y = value)) +
        geom_boxplot() +
        expand_limits(y = 0) +
        xlab("Condition") +
        ylab("Normalized Expression") +
        ggtitle(paste("Gene Expression of ", gene, sep = "")) + 
        theme_bw() +
        theme(plot.title = element_text(size = 8, face = "bold"))
}


#' Plot Volcano
#' 
plotVolcano <- function(data, gene_anno) {
    data <- as.data.frame(data)

    # Match gene symbols to ENTREZ IDs and create a 'symbol' column
    idx <- match(rownames(data), gene_anno$ENTREZ)
    data$symbol <- paste(gene_anno$ENTREZ[idx], gene_anno$symbol[idx], sep = ": ")

    # Determine genes meeting the threshold criteria
    data$threshold <- ifelse(data$padj < 0.05 & abs(data$log2FoldChange) > log(2, 2), "TRUE", "FALSE")

    # Select top 50 genes for labeling on the volcano plot based on fold change and adjusted p-value
    label.FC <- data[data$threshold == "TRUE", ]
    label.FC <- label.FC[order(abs(label.FC$log2FoldChange), decreasing = TRUE), ]
    if (nrow(label.FC) > 50) {
        label.FC <- as.character(label.FC[c(1:50), "symbol"])
    } else {
        label.FC <- as.character(label.FC$symbol)
    }
    label.p <- data[data$threshold == "TRUE", ]
    label.p <- label.p[order(abs(label.p$padj), decreasing = FALSE), ]
    if (nrow(label.p) > 50) {
        label.p <- as.character(label.p[c(1:50), "symbol"])
    } else {
        label.p <- as.character(label.p$symbol)
    }
    label <- unique(c(label.FC, label.p))

    # Assign labels to the data based on the threshold and label criteria
    data$label <- ifelse(data$symbol %in% label & data$threshold == "TRUE", as.character(data$symbol), "")

    # Create a volcano plot with adjusted p-values and fold changes
    ggplot(data = na.omit(data), aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
        geom_point(alpha = 0.4, size = 1.75) +
        scale_color_manual(values = c("grey", "firebrick")) +
        scale_x_continuous() +
        scale_y_continuous() +
        xlab("log2(FoldChange)") +
        ylab("-log10(p-value)") +
        geom_vline(xintercept = 0, colour = "black") +
        geom_vline(xintercept = c(-log(2, 2), log(2, 2)), colour = "red") +
        geom_hline(yintercept = -log(0.05, 10), colour = "red") +
        geom_text_repel(data = na.omit(data[!data$label == "", ]), aes(label = label), size = 3) +
        guides(colour = FALSE) +
        theme_bw()
}


#' Functional Analysis
#'
GOEA <- function(data) {
    # Convert input data to a data frame if it's not already in that format
    data <- as.data.frame(data)

    # Remove rows with NA values in the 'pvalue' column
    data <- data[!is.na(data$pvalue),]

    # Filter genes based on adjusted p-value and fold change criteria
    genes <- rownames(data[data$padj < 0.05 & abs(data$log2FoldChange) > log(2, 2),])

    # Perform gene ontology enrichment analysis using 'enrichGO' from clusterProfiler
    eGO <- clusterProfiler::enrichGO(
        gene = genes,
        OrgDb = org.Mm.eg.db,
        ont = "BP",
        pAdjustMethod = "bonferroni",
        pvalueCutoff = 0.1,
        qvalueCutoff = 0.1,
        readable = TRUE
    )

    # Generate a dot plot to visualize enriched gene ontology terms
    dotplot(eGO, showCategory = 20, font.size = 8)
}
