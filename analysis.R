install.packages("BiocManager")

library(BiocManager)

BiocManager::install(c("edgeR", "DESeq2"))

install("Biobase")

library("edgeR")

BiocManager::install("DESeq2")

BiocManager::install("org.Mm.eg.db")

BiocManager::install("gplots")

BiocManager::install("Glimma")

BiocManager::install("sva")

install.packages(c("assertthat", "dplyr", "ggplot2"))

library("RColorBrewer")

library("edgeR")

library("gplots")

library("DESeq2")

library("assertthat")

library("dplyr")

library("ggplot2")

library("Glimma")

library("org.Mm.eg.db")

library("limma")

#descriptive statistics of genes expression

#boxplot of different experiments

boxplot(JQ1E2Rep2_FPKM, JQ1E2Rep1_FPKM, JQ1vehRep2_FPKM, JQ1vehRep1_FPKM,
        E2Rep2_FPKM, E2Rep1_FPKM, VehRep2_FPKM, VehRep1_FPKM)

boxplot(Initial_data_R$MCF7_JQ1_E2_rep2, Initial_data_R$MCF7_JQ1_E2_rep1,
        Initial_data_R$MCF7_JQ1_Veh_rep2, Initial_data_R$MCF7_JQ1_Veh_rep1,
        Initial_data_R$MCF7_E2_rep2, Initial_data_R$MCF7_E2_rep1,
        Initial_data_R$MCF7_Veh_rep2, Initial_data_R$MCF7_Veh_rep1,
        main = "boxplots for raw data",
        xlab = "number of experiment",
        ylab = "expression levels")

boxplot(Initial_data_R$MCF7_JQ1_E2_rep2_FPKM, Initial_data_R$MCF7_JQ1_E2_rep1_FPKM,
        Initial_data_R$MCF7_JQ1_Veh_rep2_FPKM, Initial_data_R$MCF7_JQ1_Veh_rep1_FPKM,
        Initial_data_R$MCF7_E2_rep2_FPKM, Initial_data_R$MCF7_E2_rep1_FPKM,
        Initial_data_R$MCF7_Veh_rep2_FPKM, Initial_data_R$MCF7_Veh_rep1_FPKM,
        main = "boxplots for FPKM",
        xlab = "experiment number",
        ylab = "expression levels")

summary(Initial_data_R$MCF7_JQ1_E2_rep2_FPKM, Initial_data_R$MCF7_JQ1_E2_rep1_FPKM,
        Initial_data_R$MCF7_JQ1_Veh_rep2_FPKM, Initial_data_R$MCF7_JQ1_Veh_rep1_FPKM,
        Initial_data_R$MCF7_E2_rep2_FPKM, Initial_data_R$MCF7_E2_rep1_FPKM,
        Initial_data_R$MCF7_Veh_rep2_FPKM, Initial_data_R$MCF7_Veh_rep1_FPKM)

# histogram

#PartB data normalization with EdgeR

# creating dgelist

Counts <- (Data_loading_RAW_R2)

colnames(Counts) <- c("MCF7_JQ1_E2_rep2", "MCF7_JQ1_E2_rep1",
                      "MCF7_JQ1_Veh_rep2", "MCF7_JQ1_Veh_rep1",
                      "Control", "Control.1", "COntrol.2", "Control.3")

dim(Counts)

head(Counts)

dgList <- DGEList(counts = Counts, genes = rownames(Counts))

dgList

dgList$samples

head(dgList$counts)

head(dgList$genes)

# Filtering process

countsPerMillion <- cpm(dgList) # cpm computes the counts per million

summary(countsPerMillion)

countCheck <- countsPerMillion > 1 # Check which genes > 1 cpm

head(countCheck) # This produces a logical matrix with TRUEs and FALSEs

table(rowSums(countCheck)) # we see How many genes have TRUEs in all the samples

keep <- which(rowSums(countCheck) >= 2) # Keep genes with > 1 cpm in at least 2 samples

dgList <- dgList[keep, ] # How many genes do we have left?

summary(cpm(dgList)) # Compare this to the original summary

head(dgList$counts)

countsPerMillion <- cpm(dgList) # Calculate cpm for the filtered dgList

#use the first sample

plot(countsPerMillion[, 1], dgList$counts[, 1], ylim = c(0, 50), xlim = c(0, 2))

# Add lines for 5 and 10 read counts

abline(v = 0.67, h = 10, col = "red")

abline(v = 1.34, h = 20, col = "blue")

# TMM normalisation

# How does the function work?

?calcNormFactors

# Check the "norm.factors" column

dgList$samples

# Apply TMM normalisation to the DGEList object

dgList <- calcNormFactors(dgList, method = "TMM") ## Normalisation for composition bias

#Quality control, library size and distribution plots

dgList$samples$lib.size #First, we can check how many reads we have for each sample in the dglist

#We can also plot the library sizes as a barplot to see whether there are any major discrepancies between the samples more easily.

barplot(dgList$samples$lib.size, names = colnames(dgList), las = 2)

title("Barplot of library sizes")

# the plot does not look good need to fix naming

#creating boxplot

logcounts <- cpm(dgList, log = TRUE)

boxplot(logcounts, xlab = "", ylab = "Log2 counts per million", las = 2)

abline(h = median(logcounts), col = "blue")

title("Boxplots of logCPMs (normalised)")

# Multidimensional scaling plots

plotMDS(logcounts)

#Hierarchical clustering with heatmaps

# We estimate the variance for each row in the logcounts matrix

var_genes <- apply(logcounts, 1, var)

head(var_genes)

# Get the gene names for the top 500 most variable genes

select_var <- names(sort(var_genes, decreasing = TRUE))[1:500]

head(select_var)

# Subset logcounts matrix

highly_variable_lcpm <- logcounts[select_var, ]

dim(highly_variable_lcpm)

head(highly_variable_lcpm)

## creating heatmap

?heatmap.2

x <- highly_variable_lcpm

y <- data.matrix(x)

heatmap(as.matrix(highly_variable_lcpm),
        main = "Top 500 most variable genes across samples")

##Batch effect correction with sva

# Take the log of the data

logBulkTPM <- log2(Counts + 1)

# Find the highest and lowest expression per row

maxTPM <- apply(logBulkTPM, 1, max)

minTPM <- apply(logBulkTPM, 1, min)

# Filter out genes according to our criteria

logBulkTPM <- logBulkTPM[(maxTPM >= 2) & ((maxTPM - minTPM) >= 2), ]

# Select the top 500 genes by variance

numTopGenes <- 500

idxOfTop <- order(apply(logBulkTPM, 1, var), decreasing = TRUE)[1:numTopGenes]

logBulkTPM <- logBulkTPM[idxOfTop, ]

# different method for PCA analysis

data <- Counts[, c(1:8)]

pc <- princomp(data, cor = TRUE, scores = TRUE)

pc

plot(pc, type = "l") ### elbow curve

biplot(pc)

### estimate dispersions :biological coefficient of variation (BCV)

# Common dispersion

## differential gene expression analysis

# We can use the column names of the DGEList object to extract the sample data

## before creating the the matrix create table with names control and tested

colnames(dgList)

sampleType <- rep("Control", ncol(dgList))

sampleType[grep("JQ1", colnames(dgList))] <- "JQ1"

# Finally, we can build the design matrix

designMat <- model.matrix(~sampleType)

View(designMat)

# Common dispersion

dgList <- estimateGLMCommonDisp(dgList, design = designMat)

plotBCV(dgList)

# Trended dispersion

dgList <- estimateGLMTrendedDisp(dgList, design = designMat)

plotBCV(dgList)

# Tagwise dispersion

dgList <- estimateGLMTagwiseDisp(dgList, design = designMat)

plotBCV(dgList)

### p-value(Colm)

apply(highly_variable_lcpm, 1, function(x) {
  t.test(x[1:4], x[5:8])$p.value
})

### the genes according to the p-value and choose most higly significant

## find the gene that is most significan and use it for biological anlysis

### gene enrichment anlysis with clusterprofiler

table(highly_variable_lcpm)

# Fit the GLM

fit <- glmFit(dgList, designMat)

# Determine which coefficient to select from the design matrix

View(designMat)

# Calculate likelihood ratio tests

lrt <- glmLRT(fit, coef = 2)

# Table the result of the DGE analysis

edgeR_result <- topTags(lrt, n = Inf, adjust.method = "BH", sort.by = "PValue")$table

View(edgeR_result)

# Filter out non-DEGs according to a threshold of adjusted P < 0.05

degs <- edgeR_result[edgeR_result$FDR < 0.05, ]

# Extract up- and down-regulated DEGs using the fold change

upreg_degs <- degs[degs$logFC > 0, ]

downreg_degs <- degs[degs$logFC < 0, ]

# Store the DEG information

write.csv(upreg_degs, file = "./upregulated_genes.csv",
          row.names = FALSE, quote = FALSE)

write.csv(downreg_degs, file = "./downregulated_genes.csv",
          row.names = FALSE, quote = FALSE)

# Obtain a list of DEG names

deGenes <- decideTestsDGE(lrt, p = 0.05)

deGenes <- rownames(lrt)[as.logical(deGenes)]

# Display a smear plot and add lines for 1 and -1 logFCs

plotSmear(lrt, de.tags = deGenes)

abline(h = c(-1, 1), col = "blue")

# Install and load the required package

BiocManager::install("clusterProfiler")

library(clusterProfiler)

# KEGG pathway enrichment for top 1000 DEGs

pathEnrich <- enrichKEGG(gene = degs$genes[1:1000],
                         organism = "hsa", # Homo sapiens
                         keyType = "kegg", # KEGG database
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH",
                         universe = edgeR_result$genes) # Background genes

View(pathEnrich@result)

# Store the significantly enriched pathways (P < 0.05)

signif_pathEnrich <- pathEnrich@result[pathEnrich@result$p.adjust < 0.05, ]

write.csv(signif_pathEnrich, file = "./signif_pathEnrich.csv",
          row.names = FALSE, quote = FALSE)

# Show the 10 most significant pathways as a dot plot

dotplot(pathEnrich, x = "Count", showCategory = 10, orderBy = "Count")
