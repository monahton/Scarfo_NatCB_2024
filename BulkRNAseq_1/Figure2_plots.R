suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(RColorBrewer))

## Input featureCounts file
fcount.file <- "Scarfo_HEC2023/BulkRNAseq_1/featureCounts_results_corrected.txt"
## Contrasts to perform
contr <- "ACE_pos:ACE_neg"
contr.list <- unlist(strsplit(contr, ","))
## Design Formula
formula.str <- "donor,condition"
design.str <- gsub(",", "+", formula.str)
num_cond <- str_count(formula.str, ",") + 1
design.formula <- as.formula(paste("~", design.str, sep = ""))
## Read featureCounts results
if (!file.exists(fcount.file)) {
  stop("Error: featureCounts file not found.")
} else {
  read.counts <- read.table(fcount.file, 
                            header = TRUE, 
                            check.names = FALSE) %>%
    column_to_rownames(var = "Geneid") %>%
    select(-c(1:5)) 
}
orig.names <- names(read.counts)
ds <- list()
cond <- list()
for (i in seq(orig.names[1:length(orig.names)])) {
  ds[i] <- unlist(strsplit(orig.names[i], ":"))[2]
  cond[i] <- unlist(strsplit(orig.names[i], ":"))[3]
}
ds <- unlist(ds)
cond <- unlist(cond)
## Samples - Conditions
sample_info <- as.data.frame(str_split_fixed(cond, ",", num_cond))
colnames(sample_info) <- str_split_fixed(formula.str, ",", num_cond)
rownames(sample_info) <- ds
colnames(read.counts) <- ds
sample_info$condition <- factor(sample_info$condition)
last_condition <- str_split_fixed(formula.str, ",", num_cond)[num_cond]
batch_cond <- str_split_fixed(formula.str, ",", num_cond)[num_cond - 1]
## thresholds
adj.pval.thr <- 0.05
logfc.thr <- 0
## DESeq2
DESeq.ds <- DESeqDataSetFromMatrix(countData = read.counts,
                                   colData = sample_info,
                                   design = design.formula)
## Remove genes without any counts
DESeq.ds <- DESeq.ds[rowSums(counts(DESeq.ds)) > 0, ]
## RunDGE
DESeq.ds <- DESeq(DESeq.ds)
## obtain regularized log-transformed values
DESeq.rlog <- rlogTransformation(DESeq.ds, blind = TRUE)
assay(DESeq.rlog) <- limma::removeBatchEffect(assay(DESeq.rlog), DESeq.rlog[[batch_cond]])
ctrs <- unlist(strsplit(contr.list[1], ":"))
c1 <- ctrs[1]
c2 <- ctrs[2]
DGE.results <- results(DESeq.ds,
                       pAdjustMethod = "BH",
                       independentFiltering = TRUE,
                       contrast = c(last_condition, c1, c2),
                       alpha = adj.pval.thr)
# Change format
DGE.results.annot <- as.data.frame(DGE.results) %>% 
  dplyr::select(log2FoldChange, lfcSE, baseMean, pvalue, padj)
colnames(DGE.results.annot) <- c("logFC", "lfcSE", "baseMean", "PValue", "FDR")
## sort the results according to the adjusted p-value
DGE.results.sorted <- DGE.results.annot[order(DGE.results.annot$FDR), ]
## Subset for only significant genes
DGE.results.sig <- subset(DGE.results.sorted, 
                          FDR < adj.pval.thr & abs(logFC) > logfc.thr)
DGEgenes <- rownames(DGE.results.sig)


# Fig2A ----------------------------------------------------------------

## surface genes
if(!file.exists("Scarfo_HEC2023/BulkRNAseq_1/table_S3_surfaceome.xlsx")){
  download.file("http://wlab.ethz.ch/surfaceome/table_S3_surfaceome.xlsx", 
                destfile = "Scarfo_HEC2023/BulkRNAseq_1/table_S3_surfaceome.xlsx")  
}
sg <- read.xlsx("Scarfo_HEC2023/BulkRNAseq_1/table_S3_surfaceome.xlsx", 
                sheet = "in silico surfaceome only",
                colNames = T,
                startRow = 2)
sig_sg <- DGE.results.sig[c(row.names(DGE.results.sig) %in% sg$UniProt.gene), ]
top10_sig_sg <- sig_sg[c(1:10),]

s_info <- sample_info[sample_info[[last_condition]] == c1 | sample_info[[last_condition]] == c2,, drop = FALSE]
hm.mat_DGEgenes <- assay(DESeq.rlog)[row.names(top10_sig_sg), ]
hm.mat_DGEgenes.filt <- hm.mat_DGEgenes[, row.names(s_info)]

pheatmap::pheatmap(hm.mat_DGEgenes.filt,
                   scale = "row",
                   color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(100),
                   annotation_col = s_info,
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   show_rownames = TRUE,
                   show_colnames = TRUE,
                   angle_col = 90,
                   fontsize = 15)

