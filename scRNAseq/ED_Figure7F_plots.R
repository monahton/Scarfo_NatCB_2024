suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(clusterProfiler))

## load Seurat data
sample <- readRDS("Scarfo_HEC2023/scRNAseq/WNTd_HE.rds")
## Seurat markers
markers <- FindAllMarkers(sample, 
                          only.pos = T, 
                          min.pct = 0.1, 
                          logfc.threshold = 0.25)
## define RUNX1 clusters
runx1.clusters <- subset(markers, markers$gene=="RUNX1")
runx1.clusters <- droplevels(runx1.clusters$cluster) 
## subset RUNX1 clusters into new object
runx1 <- subset(sample, idents = runx1.clusters)
## differential markers
degs_clu11_vs_clu0_1_2 <- FindMarkers(object = runx1, 
                                      ident.1 = 11, 
                                      ident.2 = c(0,1,2), 
                                      only.pos = FALSE, 
                                      min.pct = 0, 
                                      test.use = "wilcox", 
                                      logfc.threshold = 0, 
                                      min.cells.group = 0)
## thresholds
adj.pval.thr <- 0.05
logfc.thr <- 0


# Remove genes with no FDR
markers <- degs_clu11_vs_clu0_1_2[!is.na(degs_clu11_vs_clu0_1_2$p_val_adj), ]
gene_univ <- row.names(markers)
# get positive and negative gene names
gene.res.top <- subset(markers, p_val_adj < adj.pval.thr  & abs(avg_logFC) > logfc.thr)
gene.res.top.pos <- subset(gene.res.top, avg_logFC > 0)
gene.res.top.neg <- subset(gene.res.top, avg_logFC < 0)
organism.db <- org.Hs.eg.db
contr.enrich_path <- enrichPathway(gene = gene.res.top.neg,
                                   organism = organism.db,
                                   universe = gene_univ,
                                   pvalueCutoff = adj.pval.thr)
contr.enrich_path.symb <- setReadable(contr.enrich_path, org_db, keyType = "ENTREZID")
toi <- c("RHO GTPases activate IQGAPs", "RHO GTPases Activate ROCKs")
GO_terms_int <- contr.enrich_path.symb[contr.enrich_path.symb$Description %in% toi, ]
colnames(GO_terms_int)[2] <- "Terms"
GO_terms_int$`Adjusted p-value` <- GO_terms_int$qvalue

ggplot(data=GO_terms_int, 
            aes(x=Terms, 
                y=Count, 
                fill = `Adjusted p-value`)) +
  geom_bar(stat="identity", 
           width = 0.7) +
  scale_fill_gradient(low="red",high="blue") +
  coord_flip()+
  xlab("Terms\n") +
  ggtitle("Downregulated in Cluster 11 vs 0,1,2") +
  theme_bw() +
  theme(aspect.ratio = 0.8)










