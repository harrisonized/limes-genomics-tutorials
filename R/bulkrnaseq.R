wd = dirname(this.path::here())  # wd = '~/github/R/limes-genomics-tutorials'

suppressPackageStartupMessages(library('logr'))
library('optparse')

suppressPackageStartupMessages(library('dplyr'))
library("reshape2")
library("ggplot2")
library("ggrepel")
library("ggpubr")
library("pheatmap")

suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("clusterProfiler"))
library("biomaRt")
library("org.Mm.eg.db")

import::from(file.path(wd, 'R', 'functions', 'plotting.R'),
    'GeneExprBox', 'plotVolcano', 'GOEA',
    .character_only=TRUE)

# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='data/GSE63310_RAW',
                metavar='data/GSE63310_RAW', type="character",
                help="input directory"),

    make_option(c("-m", "--metadata-file"), default='data/GSE63310_series_matrix.txt',
                metavar='data/GSE63310_series_matrix.txt', type="character",
                help="metadata file"),

    make_option(c("-o", "--output-dir"), default="output",
                metavar="output", type="character",
                help="output directory for the data"),

    make_option(c("-b", "--biomart-table"), default='data/biomart_table.txt',
                metavar='data/biomart_table.txt', type="character",
                help="biomart table"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]

# Start Log
start_time = Sys.time()
log <- log_open(paste0("bulkrnaseq-", strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read Data

# Read counts matrix
files <- dir(
    file.path(wd, opt[['input-dir']]),
    full.names=TRUE
)
for(i in 1:length(files)) {
    tmp <- read.delim(files[i])
    colnames(tmp)[3] <- unlist(lapply(
        strsplit(tools::file_path_sans_ext(basename(files[i])), "_"),
        function(x) x[1])
    )  # clean column names
    if(i == 1){
        counts <- as.data.frame(tmp[,c(1,3)])
    } else {
        counts <- merge(counts,tmp[,-2], by="EntrezID")
    }
}
rownames(counts) <- as.character(counts$EntrezID)  # move EntrezID to rowname
counts$EntrezID <- NULL

# Read metadata
sampleinfo_raw <- read.delim(
    file.path(wd, opt[['metadata-file']]),
    skip=30, header=FALSE, stringsAsFactors = FALSE
) 

# Preprocess the metadata file
sampleinfo <- t(sampleinfo_raw[
    grep("!Sample_title|!Sample_geo_accession|!Sample_characteristics_ch1", sampleinfo_raw$V1), ]
)
sampleinfo <- sampleinfo[-c(4,10),]  # Removes rows at positions 4 and 10 from sampleinfo.

sampletable <- as.data.frame(sampleinfo[-1,1] )
colnames(sampletable)[1] <- "description"

for (i in 1: ncol(sampleinfo)) {
    if(i==1){
        temp1 <- unlist(lapply(strsplit(sampleinfo[,i],", ")[-1],function(x) x[1]))
        sampletable[,2] <- temp1
        colnames(sampletable)[2] <- "cond"
        temp2 <- unlist(lapply(strsplit(sampleinfo[,i],", ")[-1],function(x) x[2]))
        sampletable[,3] <- temp2
        colnames(sampletable)[3] <- "cell"
        sampletable[,1] <- paste(temp1, temp2, sep="_")
    } else if(i==2){
        sampletable[,4] <- sampleinfo[-1,2]
        colnames(sampletable)[4] <- "GSE"
    } else{
        name <- unlist(lapply(strsplit(sampleinfo[,i],": ")[-1],function(x) x[1]))[1]
        temp <- unlist(lapply(strsplit(sampleinfo[,i],": ")[-1],function(x) x[2]))
        sampletable[,i+2] <- temp 
        colnames(sampletable)[i+2] <- name
        }
}
rownames(sampletable) <- sampletable$description


counts <- counts[,as.character(sampletable$GSE)] # filter columns
# identical(colnames(counts),sampletable$GSE) # test if sample annotation and expression table are compatible
colnames(counts) <- sampletable$description # rename columns of the counts matrix


# ----------------------------------------------------------------------
# Gene Annotation

# Get BioMart
if (file.exists(file.path(wd, opt[['biomart-table']]))) {
    bm <- read.delim(file.path(wd, opt[['biomart-table']]), header= TRUE, stringsAsFactors=FALSE)
} else {
    
    # Access ensembl data base and download relevant information:
    ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")  # needs internet connection!
    bm <- getBM(
        attributes = c(
            "entrezgene_id",
            "ensembl_gene_id", 
            "external_gene_name", 
            "description", 
            "gene_biotype",
            "chromosome_name"),
        mart = ensembl)

    # Export for future reuse
    write.table(bm, file.path(wd, opt[['biomart-table']]), sep="\t")
}

# create gene annotation dataframe
gene_anno <- as.data.frame(rownames(counts))
colnames(gene_anno) <- "ENTREZ"
IDindex <- match(gene_anno[['ENTREZ']], bm$entrezgene)
gene_anno[['ENSID']] <- bm$ensembl_gene_id[IDindex]          # Store Ensembl gene IDs
gene_anno[['symbol']] <- bm$external_gene_name[IDindex]      # Store external gene symbols/names
gene_anno[['description']] <- bm$description[IDindex]        # Store gene descriptions
gene_anno[['type']] <- bm$gene_biotype[IDindex]              # Store gene biotypes/types
gene_anno[['chr']] <- bm$chromosome_name[IDindex]            # Store chromosome names

# how many are missing?
# sum(is.na(gene_anno$ENSID))


# ----------------------------------------------------------------------
# Gene Type Distribution

# Create a data frame 'TypeCounts' containing the counts of unique values in the 'gene_anno$type' column
TypeCounts <- as.data.frame(table(gene_anno$type))

# Rename the columns of 'TypeCounts' to "Type" and "Frequency"
colnames(TypeCounts) <- c("Type", "Frequency")

# Filter 'TypeCounts' to retain rows where 'Frequency' is greater than 0
TypeCounts <- subset(TypeCounts, Frequency > 0)

# Create a bar plot using 'ggplot' to visualize the frequency of different types
ggplot(TypeCounts, aes(x = Type, y = Frequency, label = Frequency)) +
  geom_bar(stat = "identity", fill = "lightblue", colour = "grey") +  # Create a bar plot
  theme_bw() +  # Set a black-and-white theme
  geom_text(size = 3, position = position_stack(vjust = 1)) +  # Add text labels on top of bars
  guides(fill = FALSE) +  # Remove the legend for fill
  theme(
    text = element_text(size = 10),  # Set text size
    axis.text.x = element_text(angle = 90, hjust = 1)  # Rotate x-axis labels
  ) +
  xlab("")  # Set an empty label for the x-axis


# ----------------------------------------------------------------------
# Preprocessing

# Create a DESeqDataSet object from count data, sample metadata, and a design formula
dds <- DESeqDataSetFromMatrix(
  countData = counts,    # The count data matrix (gene counts for samples)
  colData = sampletable, # The sample metadata including conditions, batches, etc.
  design = ~ batch + cond # The design formula specifying the model (batch effect and condition of interest)
)

# How many genes in this data set are not expressed at all?
# sum(rowSums(counts)==0)

# Calculate row sums for the count matrix 'counts' and create a data frame 'rsums'
rsums <- as.data.frame(rowSums(counts)[rowSums(counts) > 0])

# Rename the column of 'rsums' as "rowsums"
colnames(rsums) <- "rowsums"

# Create a histogram using 'ggplot' with data from 'rsums' and plot the distribution of row sums
ggplot(rsums, aes(x = rowsums)) + 
  geom_histogram(bins = 100) +  # Create a histogram with 100 bins
  scale_x_log10(limits = c(1, 10000000)) +  # Set x-axis to log scale within specified limits
  theme_bw()  # Apply a black-and-white theme


# Create a logical vector 'genes_to_keep' based on a condition for row sums in the 'counts' matrix
genes_to_keep <- rowSums(counts) >= 1000

# Calculate the sum of TRUE values in the 'genes_to_keep' vector
sum(genes_to_keep)

rsums_filt <- as.data.frame(rowSums(counts[genes_to_keep,]))
colnames(rsums_filt) <- "rowsums"

ggplot(rsums_filt, aes(x = rowsums)) + 
  geom_histogram(bins=100) + 
  scale_x_log10(limits=c(1,10000000)) +
  theme_bw()

dds <- dds[genes_to_keep,]


# ----------------------------------------------------------------------
# Calculate Statistics

dds <- DESeq(dds)


# Reshape the subset of 'counts' matrix (filtered by 'genes_to_keep') into a long format using 'melt' function
counts_melt <- melt(counts[genes_to_keep,])

# Rename the columns of 'counts_melt' as "Sample" and "Counts"
colnames(counts_melt) <- c("Sample", "Counts")

# Create a violin plot using 'ggplot' to visualize the distribution of counts across samples
ggplot(counts_melt, aes(x = Sample, y = log(Counts + 1, 10))) + 
  geom_violin() +  # Create a violin plot
  scale_y_continuous() +  # Set the y-axis to be continuous
  ylab("raw Counts") +  # Label for the y-axis
  theme_bw() +  # Apply a black-and-white theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels

dds$sizeFactor

norm <- as.data.frame(counts(dds, normalized=T))

norm_melt <- melt(norm)
colnames(norm_melt) <- c("Sample","normCounts")

ggplot(norm_melt, aes(x=Sample, y=log(normCounts+1,10))) + 
  geom_violin() +
  scale_y_continuous()+
  ylab("normalized Counts")+
  theme_bw()+
  theme(axis.text.x = element_text(angle =90, hjust = 1))


# ----------------------------------------------------------------------
# Variance Stabilization

rlogCounts <- rlog(dds, blind = TRUE)


# ----------------------------------------------------------------------
# EDA

# Sort the 'norm' matrix rows by row means in decreasing order and select the top 10 rows
tmp <- norm[order(rowMeans(norm), decreasing = TRUE), ]
tmp <- tmp[1:10, ]

# Reshape the transposed 'tmp' matrix into a long format using 'melt'
tmp <- melt(t(tmp))
colnames(tmp) <- c("sample", "gene", "value")  # Rename columns in the melted data frame

# Match genes in 'tmp' with 'bm$entrezgene' and create a new column 'symbol' based on the matching indices
idx <- match(tmp$gene, bm$entrezgene)
tmp$symbol <- as.factor(paste(tmp$gene, bm$external_gene_name[idx], sep = ": "))
tmp$symbol <- factor(tmp$symbol, levels = rev(unique(tmp$symbol)))

# Create a boxplot using 'ggplot' to visualize the expression of top 10 highly expressed genes
ggplot(tmp, aes(x = tmp$symbol, y = value)) +
    geom_boxplot() +  # Create a boxplot
    xlab("Gene") +  # Label for the x-axis
    ylab("Normalized Expression") +  # Label for the y-axis
    ggtitle("Expression of 10 highest expressed genes") +  # Title for the plot
    theme_bw() +  # Apply a black-and-white theme
    coord_flip() +  # Flip the coordinates (swap x and y axes)
    theme(
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1),  # Adjust x-axis text
        plot.title = element_text(size = 8, face = "bold")  # Adjust plot title appearance
    )


# ----------------------------------------------------------------------
# Boxplot of a selected gene

# Call the function 'GeneExprBox' with the gene symbol "Cd82" to visualize its expression
GeneExprBox("Cd82", norm, gene_anno)


# ----------------------------------------------------------------------
# Hierarchical Clustering of normalized counts


# Transpose the 'norm' matrix to get data suitable for heatmap visualization
HeatmapData <- t(norm)

# Scale and center the transposed data for heatmap visualization
HeatmapData_scale <- t(scale(HeatmapData, scale = TRUE, center = TRUE))

# The `scale` function standardizes the columns (genes) of the data matrix:
# - `scale = TRUE` scales each column (gene) to have mean 0 and standard deviation 1.
# - `center = TRUE` centers the data by subtracting the mean of each column.

# Generate a heatmap using the processed and scaled data with 'pheatmap'
pheatmap(HeatmapData_scale,
         main = "Hierarchical Clustering of all present genes",
         show_rownames = FALSE,
         show_colnames = TRUE,
         fontsize = 7, 
         annotation_col = sampletable[, c("cond", "batch")])


# ----------------------------------------------------------------------
# Sample-to-Sample distance

# Calculate sample distances based on normalized counts using rlog transformation
sampleDist <- dist(t(assay(rlogCounts)))

# Convert the sample distance object to a matrix for heatmap visualization
sampleDistmatrix <- as.matrix(sampleDist)

# Set row and column names for the distance matrix using sample descriptions from 'sampletable'
rownames(sampleDistmatrix) <- sampletable$description
colnames(sampleDistmatrix) <- sampletable$description

# Generate a heatmap showing sample distances based on rlog-transformed counts
pheatmap(sampleDistmatrix,
         clustering_distance_rows = sampleDist,
         clustering_distance_cols = sampleDist, 
         main = "Sample Distances based on rlog-transformed counts",
         annotation_row = sampletable[, c("cond", "batch")], 
         annotation_col = sampletable[, c("cond", "batch")],
         show_rownames = FALSE,
         fontsize = 10)


# ----------------------------------------------------------------------
# PCA

# Perform PCA (Principal Component Analysis) on rlog-transformed counts
pcaData <- plotPCA(
  rlogCounts,
  intgroup = c("cond", "batch"),
  ntop = nrow(assay(rlogCounts)),
  returnData = TRUE
)

# Calculate the percentage of variance explained by each principal component
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Create a scatter plot of the first two principal components (PC1 and PC2)
# Color points by 'cond' (condition) and shape by 'batch'
ggplot(pcaData, aes(x = PC1, y = PC2, color = cond, shape = batch)) +
    geom_point(size = 3) +  # Plot points with size 3
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +  # Label x-axis with PC1 variance
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +  # Label y-axis with PC2 variance
    theme_bw() +  # Apply a black-and-white theme
    ggtitle("Principal Component Analysis based on rlog-transformed counts")  # Set plot title


# Loadings

# Perform Principal Component Analysis (PCA) on rlog-transformed counts
pca <- prcomp(t(assay(rlogCounts)))

# Specify the principal component (PC) of interest
PC <- 1 

# Extract loadings for the specified principal component and select top 20 loadings
loadings <- pca$rotation[, PC]
loadings <- loadings[order(loadings, decreasing = TRUE)]
loadings_20 <- loadings[c(1:20, (length(loadings) - 19):length(loadings))]
loadings <- names(loadings_20)

# Subset the normalized counts data using the selected top loadings
HC <- norm[loadings,]

# Match row names of HC with 'bm$entrezgene' to obtain gene names for annotation
idx <- match(rownames(HC), bm$entrezgene)
rownames(HC) <- paste(rownames(HC), bm$external_gene_name[idx], sep=": ")

# Scale the transposed HC data for heatmap visualization
HC_scale <- t(scale(t(HC)))

# Create a heatmap to visualize hierarchical clustering of top 20 loadings of PC1
pheatmap(HC_scale,
         main = paste("Hierarchical Clustering of top 20 loadings of PC1 in both directions"),
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = sampletable[, c("cond", "batch")],
         cluster_cols = TRUE,
         fontsize = 8)


# ----------------------------------------------------------------------
# DEG analysis

# Perform differential expression analysis comparing conditions 'Basal' and 'LP' using DESeq2
Results_BasalvsLP <- results(dds, contrast = c("cond", "Basal", "LP"), alpha = 0.05)
Results_BasalvsLP

# Convert the results of differential expression (Results_BasalvsLP) into a data frame
# and create a histogram to visualize the distribution of p-values
ggplot(as.data.frame(Results_BasalvsLP), aes(x = pvalue)) + 
  geom_histogram(aes(y = ..count..), binwidth = 0.01) +
  theme_bw()

summary(Results_BasalvsLP)

df <- as.data.frame(Results_BasalvsLP)

# Upregulated genes
nrow(na.omit(df[df$padj < 0.05 & df$log2FoldChange > log(2,2),]))
# Downregulated genes
nrow(na.omit(df[df$padj < 0.05 & df$log2FoldChange < -log(2,2),]))

# Call the plotVolcano function with Results_BasalvsLP as input data
plotVolcano(Results_BasalvsLP, gene_anno)


# ----------------------------------------------------------------------
# Functional analysis

GOEA(Results_BasalvsLP)


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
