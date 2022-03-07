library('DESeq2')
library(edgeR)
library(tidyverse)
library("pheatmap")
library(grid)
library(AnnotationHub)

# load count matrix
cts <- as.matrix(read_csv("count_matrix_filtered.csv"))
tags <- cts[,1]
cts <- as.matrix(cts[,-1], row.names = cts[,1])
mode(cts) <- "integer"
rownames(cts) <- tags
head(cts,2)
# these are raw counts so technical variation is expected.

# # Let's see if TMM normalization (using function from edgeR) addresses this
# DGEList <- DGEList(cts)
# DGEList.norm <- calcNormFactors(DGEList,method='TMM')
# cpm.norm <- cpm(DGEList.norm)
# head(cpm.norm,2)
# # much better

# loading metadata
colData <- data.frame(read_csv("metaData.csv"))
rownames(colData) <- colData[,1] 
all(rownames(colData) == colnames(cts))

# construct DESeqDataSet Object
dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=colData,
                              design=~condition)

# set factor levels manually; the first factor is the reference to which other levels are compared

# subset relevant samples
comp_samples = c('LB0','LB1','M9LQ','M9SF') # first is control, second is treatment
# for this dataset, want to keep  all conditions
# for NAND, need to rethink this since there are so many conditions
NAME = 'condition_M9LQ_vs_LB0' # 'condition_M9LQ_vs_LB0', 'condition_M9SF_vs_LB0'
# note that using all samples instead of downsampling is going to induce possible overestimation of dispersion 
# also  mean is way off
dds <- dds[,dds$condition %in% comp_samples]

dds$condition <- factor(dds$condition, levels = comp_samples)
dds$condition 



# dds <- estimateSizeFactors(dds)
# dds <- estimateDispersions(dds)

# get all means 
# baseMeanPerLvl <- sapply( levels(dds$condition), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$condition == lvl, drop=F] ) )

dds <- DESeq(dds)
resultsNames(dds)
# write.csv(counts(dds,normalized=TRUE),file='count_matrix_filtered_normalized.csv')
# write.csv(dispersions(dds),file='dispersions_M9LQ_vs_LB0')

res <- results(dds,name=NAME)

# write.csv(as.data.frame(res), file="condition_M9SF_vs_LB1_results.csv")

# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# resLFC <- lfcShrink(dds, coef='Condition_LB1_vs_LB0',type='apeglm')
# resLFC

# resOrdered <- res[order(res$log2FoldChange),]
# head(resOrdered,100)

# informative plots
plotCounts(dds,gene=which.max(res$log2FoldChange),intgroup='condition')
plotMA(res,colSig='blue',colLine='red')
plotDispEsts(dds, ylim = c(1e-6, 1e1) )

# volcano plot
topT <- as.data.frame(res)
#Adjusted P values (FDR Q values)
with(topT, plot(res$log2FoldChange, -log10(res$padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(topT$pvalue[topT$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)

# rlog transform
cts_rlog = rlogTransformation(dds)

# PCA
plotPCA(cts_rlog, intgroup=c( "condition","replicate")) + coord_fixed()

# clustermap
# select = order(rowMeans(assay(cts_rlog)), decreasing = TRUE)[1:50]
# hm <- pheatmap( assay(cts_rlog)[select, ],
#           scale = "column",
#           annotation_col = as.data.frame(colData(cts_rlog)[, "condition"]))

normalized_counts <- counts(dds,normalized=TRUE)

select = order(rowMeans(normalized_counts), decreasing = TRUE)[1:30]
hm <- pheatmap( normalized_counts[select, ],
                scale = "row",
                annotation_col = as.data.frame(colData(cts_rlog)[, c("condition","replicate")]))

grid.draw(rectGrob(gp=gpar(fill="black", lwd=0)))
grid.draw(hm)
grid.gedit("layout", gp = gpar(col = "white", text = ""))

# quartz.save("clustermap30.pdf", type="pdf")


# enrichment analysis ----------
# Create background dataset for hypergeometric testing using all genes tested for
# significance in the results

all_genes <- as.character(rownames(res))
# extract significant results
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_genes <- as.character(rownames(signif_res))


# run GO enrichment analysis
ego <- enrichGO


ah <- AnnotationHub()
query(ah, c("OrgDb", "Takifugu"))










