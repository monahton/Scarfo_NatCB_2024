suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(clusterProfiler))

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


# Fig1B ---------------------------------------------------------------------

plotPCA(DESeq.rlog, intgroup = last_condition, ntop = 500) +
  theme_bw() +
  ggtitle(label = "Principal Component Analysis (PCA)",
          subtitle = "Top 500 DEG (with batch correction)") +
  geom_point(size = 2) +
  #coord_fixed(ratio=3) +
  geom_text(aes(label = rownames(DESeq.rlog@colData)), vjust = 0, nudge_y = 0.5,
            show.legend = FALSE)

# Fig1C ----------------------------------------------------------------

## DEGs
s_info <- sample_info[sample_info[[last_condition]] == c1 | sample_info[[last_condition]] == c2,, drop = FALSE]
hm.mat_DGEgenes <- assay(DESeq.rlog)[DGEgenes, ]
hm.mat_DGEgenes.filt <- hm.mat_DGEgenes[, row.names(s_info)]

pheatmap::pheatmap(hm.mat_DGEgenes.filt,
                   color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(100),
                   fontsize_row = 5,
                   show_rownames = FALSE,
                   scale = "row",
                   angle_col = 90,
                   cluster_rows = TRUE,
                   annotation_col = s_info,
                   cluster_cols = TRUE)

# Fig1D ----------------------------------------------------------------

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
markers <- c("CD34","CDH5","PECAM1","TEK","NR2F2","FLRT2","BMX","GJA5","DLL4","CXCR4","HEY2","MYB","GFI1","CD44")
DESeq.rlog_mat <- assay(DESeq.rlog)
DESeq.rlog_mat <- DESeq.rlog_mat[c(row.names(DESeq.rlog_mat) %in% markers), ]
DESeq.rlog_mat <- DESeq.rlog_mat[markers,]
# Specify colors
ann_colors = list(
  condition = c(ACE_neg = "gold3",
                ACE_pos = "turquoise3"),
  donor = c(E1 = "#E76BF3",
            E2 = "#F8766D",
            E3 = "#529EFF",
            E4 = "springgreen3")
)
pheatmap::pheatmap(DESeq.rlog_mat,
                   color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(100),
                   annotation_col = s_info,
                   annotation_names_row = FALSE,
                   annotation_colors = ann_colors,
                   show_rownames = TRUE,
                   show_colnames = TRUE,
                   cluster_rows = FALSE,
                   angle_col = 45,
                   fontsize = 14,
                   fontsize_row = 14,
                   fontsize_col = 18,
                   cellheight = 20,
                   cellwidth = 30)


# Fig1E --------------------------------------------

gmt.obj <- read.gmt("Scarfo_HEC2023/BulkRNAseq_1/c5.all.v7.2.symbols.gmt")
organism.db <- org.Hs.eg.db
# Remove genes with no FDR
DGE.results.annot <- DGE.results.annot[!is.na(DGE.results.annot$FDR), ]
gene_univ <- row.names(DGE.results.annot)
# get positive and negative gene names
gene.res.top <- subset(DGE.results.annot, FDR < adj.pval.thr  & abs(logFC) > logfc.thr)
gene.res.top.pos <- subset(gene.res.top, logFC > 0)
gene.res.top.neg <- subset(gene.res.top, logFC < 0)

## Over Representation Analysis
comp <- c("gene.res.top.pos","gene.res.top.neg")
for (i in comp) {
  print(paste0("Analyzing...",i," genes"))
  obj_i <- get(i)
  contr.enricher <- enricher(row.names(obj_i),
                             universe = gene_univ,
                             pvalueCutoff = adj.pval.thr,
                             TERM2GENE = gmt.obj)
  contr.enricher.res <- as.data.frame(contr.enricher)
  assign(paste("contr.enricher",i,sep = "."), contr.enricher.res)
}

## barplots
terms.pos <- data.frame(
  ID = c("GO_LEUKOCYTE_MIGRATION","GO_AMEBOIDAL_TYPE_CELL_MIGRATION","GO_TISSUE_MIGRATION","GO_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_MIGRATION","GO_REGULATION_OF_EPITHELIAL_CELL_MIGRATION","GO_ENDOTHELIAL_CELL_MIGRATION","GO_REGULATION_OF_ENDOTHELIAL_CELL_MIGRATION","GO_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX","GO_EXTRACELLULAR_MATRIX","GO_CELL_CELL_JUNCTION","GO_NEGATIVE_REGULATION_OF_CELL_ADHESION","GO_REGULATION_OF_CELL_CELL_ADHESION","GO_LEUKOCYTE_CELL_CELL_ADHESION","GO_EXTRACELLULAR_MATRIX_BINDING"),
  Type = c(rep("Migration", 7), rep("Adhesion", 7)))
terms.pos <- left_join(terms.pos,contr.enricher.gene.res.top.pos,by="ID")

terms.neg <- data.frame(
  ID = c("GO_POSITIVE_REGULATION_OF_CELL_CYCLE_PROCESS","GO_POSITIVE_REGULATION_OF_CELL_CYCLE","GO_POSITIVE_REGULATION_OF_CELL_CYCLE_G2_M_PHASE_TRANSITION","GO_POSITIVE_REGULATION_OF_CELL_CYCLE_PHASE_TRANSITION","GO_MITOTIC_CELL_CYCLE_CHECKPOINT","GO_MITOTIC_NUCLEAR_DIVISION","GO_POSITIVE_REGULATION_OF_MITOTIC_CELL_CYCLE"),
  Type = "Cell Cycle")
terms.neg <- left_join(terms.neg,contr.enricher.gene.res.top.neg,by="ID")
terms.to.plot <- rbind(terms.pos, terms.neg)
terms.to.plot <- na.omit(terms.to.plot)
parsed.terms <- as.data.frame(str_split_fixed(terms.to.plot$ID, "GO_", 2))
renamed.terms <- as.data.frame(str_replace_all(parsed.terms$V2, "_", " "))
colnames(renamed.terms) <- "id"
renamed.terms$Terms <- str_to_title(renamed.terms$id)
terms.to.plot <- cbind(terms.to.plot,renamed.terms)
terms.to.plot <- terms.to.plot[,c(12,10,7,2)]
terms.to.plot$Type <- factor(terms.to.plot$Type, levels = c("Cell Cycle", "Migration", "Adhesion"))
terms.to.plot$Terms <- factor(terms.to.plot$Terms, levels = terms.to.plot$Terms)
terms.to.plot <- terms.to.plot %>%
  mutate(color = ifelse(Type == "Cell Cycle", "gold3", "turquoise3"),
         regulated = ifelse(Type == "Cell Cycle", "Upregulated in ACEneg", "Upregulated in ACE+"))
terms.to.plot$regulated <- factor(terms.to.plot$regulated, 
                                  levels = c("Upregulated in ACEneg","Upregulated in ACE+"))
terms.to.plot$log_pvalue <- -log10(terms.to.plot$p.adjust)
ggplot(data=terms.to.plot, aes(x=Terms, 
                               y=log_pvalue,
                               fill=regulated,
                               xmin = min(log_pvalue),
                               xmax = max(log_pvalue))) +
  geom_bar(stat="identity", width = 0.7) +
  scale_fill_manual(values = c("gold3", "turquoise3"))+
  theme_bw() +
  coord_flip() +
  ylab(expression(-log[10]("adjusted p-value"))) +
  xlab("GO Terms") +
  #scale_x_discrete(limits = GO_terms_all$Term) +
  facet_grid(Type ~., scales = "free_y", space = "free_y") +
  theme(
    axis.title.x = element_text(colour="black",family = "Arial", size=16),
    axis.text.x = element_text(colour="black",family = "Arial", size=12),
    axis.title.y = element_text(colour="black",family = "Arial", size=16),
    axis.text.y = element_text(color="black",family = "Arial", size=12),
    legend.text = element_text(colour="black",family = "Arial", size=12),
    legend.title = element_text(colour="white",family = "Arial", size=12),
    strip.background = element_rect(colour="black",
                                    fill="white"),
    strip.text = element_text(family = "Arial", size=12))

## GSEA
DGE.results.annot.sorted <- DGE.results.annot[order(DGE.results.annot$logFC, decreasing = TRUE), ]
# LogFC with Symbols
gene.res.logfc.symbol <- as.vector(DGE.results.annot.sorted$logFC)
names(gene.res.logfc.symbol) <- row.names(DGE.results.annot.sorted)
contr.gsea <- GSEA(gene.res.logfc.symbol, 
                   TERM2GENE = gmt.obj,
                   nPerm = 10000, 
                   pvalueCutoff = adj.pval.thr)
contr.gsea.res <- as.data.frame(contr.gsea)




