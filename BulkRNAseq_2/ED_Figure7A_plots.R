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


# ED_Fig7A -----------------------------------------------------------------

gmt.obj <- read.gmt("Scarfo_HEC2023/BulkRNAseq_1/c5.all.v7.2.symbols.gmt")
organism.db <- org.Hs.eg.db
# Remove genes with no FDR
DGE.results.annot <- DGE.results.annot[!is.na(DGE.results.annot$FDR), ]
gene_univ <- row.names(DGE.results.annot)
# get positive and negative gene names
gene.res.top <- subset(DGE.results.annot, FDR < adj.pval.thr  & abs(logFC) > logfc.thr)
gene.res.top.pos <- subset(gene.res.top, logFC > 0)
gene.res.top.neg <- subset(gene.res.top, logFC < 0)
contr.enricher.neg <- enricher(row.names(gene.res.top.neg),
                           universe = gene_univ,
                           pvalueCutoff = adj.pval.thr,
                           TERM2GENE = gmt.obj)
contr.enricher.res <- as.data.frame(contr.enricher.neg)

toi <- c("GO_RESPONSE_TO_BMP", "GO_REGULATION_OF_BMP_SIGNALING_PATHWAY")
GO_terms_int <- contr.enricher.res[contr.enricher.res$Description %in% toi, ]
# reformatting the ID column
tmp <- as.data.frame(str_split_fixed(GO_terms_int$ID, "GO_", 2))
tmp2 <- as.data.frame(str_replace_all(tmp$V2, "_", " "))
colnames(tmp2) <- "id"
tmp2$Terms <- str_to_title(tmp2$id)
GO_terms_final <- cbind(GO_terms_int,tmp2)
GO_terms_final <- GO_terms_final[,c(11,9,6)]
colnames(GO_terms_final)[3] <- "Adjusted p-value"

ggplot(data=GO_terms_final, 
            aes(x=Terms, 
                y=Count, 
                fill = `Adjusted p-value`)) +
  geom_bar(stat="identity", 
           width = 0.7) +
  scale_x_discrete(limits = GO_terms_final$Terms) +
  scale_fill_gradient(low="red",high="blue") +
  coord_flip()+
  ggtitle("Downregulated in DLL4+ vs CD32+") +
  theme_bw() +
  theme(aspect.ratio = 0.8)




