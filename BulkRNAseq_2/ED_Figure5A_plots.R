suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))

## Input featureCounts file
fcount.file <- "Scarfo_HEC2023/BulkRNAseq_2/featureCounts_results_corrected.txt"
## Contrasts to perform
contr <- "DLL4_plus:CD32_plus"
contr.list <- unlist(strsplit(contr, ","))
## Design Formula
formula.str <- "condition"
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


# ED_Fig5A -----------------------------------------------------------------

## scorecard
EHT_scorecard = c("CDH5","RUNX1","PTPRC", "SPN","CD44","ITGA2B","HLF","GFI1","MLLT3","HOXA9","MKI67","NRP2","NR2F2","GJA4","HEY1","DLL4","CXCR4","SOX17","SMAD6","GJA5","HES4","MECOM","NID2","PRND","ESM1","PDGFB","PALMD","TMEM100","EDN1","LTBP4","HEY2","SULF1","IL33","CYP26B1","ADGRG6","COL23A1","GATA6","BMX","TMCC3","DKK1","AGTR2","FBN2","ELN","ALDH1A1","NKX2-3","PROCR","GATA3","GBP4","MYCN","KCNK17","MYB","STAT5A","SMIM24","RAB27B","SPINK2")
hm.mat_DGEgenes <- assay(DESeq.rlog)[DGEgenes, ]
counts_eht <- hm.mat_DGEgenes[c(rownames(hm.mat_DGEgenes) %in% EHT_scorecard), ]
eht_sig <- counts_eht[c(rownames(counts_eht) %in% row.names(DGE.results.sig)), ]
eht_sig <- eht_sig[,c(1,3,6,2,5,8)]
annot <- data.frame(row.names = colnames(eht_sig),
                    Sample = c(rep("DLL4+", 3), rep("CD32+", 3)),
                    Donor = rep(c("RS227","RS242","RS243"), 2))

annot_colors <- list(Sample=c("DLL4+" ="blueviolet",
                              "CD32+" = "gold1"),
                     Donor=c("RS227"="darkslategray3",
                             "RS242"="brown3",
                             "RS243"="limegreen"))

pheatmap::pheatmap(eht_sig,
                   scale = "row",
                   annotation_col = annot,
                   annotation_colors = annot_colors,
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   show_rownames = TRUE,
                   show_colnames = TRUE,
                   angle_col = 90,
                   legend = TRUE,
                   legend_breaks = c(-1.5,0,1.5), # legend breaks to be shown
                   legend_labels = c(-1.5,0,1.5), # legend labels to be shown
                   treeheight_row = 25*96/72, # size of row dendogram
                   treeheight_col = 15*96/72, # size of column dendogram
                   fontsize = 10*96/72,
                   fontsize_row = 8*96/72,
                   fontsize_col = 10*96/72,
                   cellwidth = 20*96/72,
                   cellheight = 8*96/72,
                   border_color = NA
)

