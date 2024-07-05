suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(clusterProfiler))

## load Seurat data
sample <- readRDS("Scarfo_HEC2023/scRNAseq/WNTd_HE.rds")

# Fig4A ------------------------------------------------------------------

DimPlot(object = sample,
        reduction = "umap",
        group.by = "RNA_snn_res.0.6",
        label = T,
        repel = T,
        pt.size = 0.8,
        label.size = 8,
        cols = c("0" = "#F8766D",
                 "1" = "#B79F00",
                 "2" = "#00BA38",
                 "3" = "#EC8239",
                 "4"= "#00B81B",
                 "5"= "#7C96FF",
                 "6"= "#00C085",
                 "7" = "#00C1A7",
                 "8" = "#00BFC4",
                 "9" = "#EF67EB",
                 "10" = "#00B2F3",
                 "11"= "#00BFC4",
                 "12"= "#00A6FF",
                 "13"= "#DB8E00",
                 "14"= "#FD61D3",
                 "15"= "#FF63B6",
                 "16"= "#619CFF",
                 "17"= "#F564E3",
                 "18"= "#00BADE",
                 "19"= "#B385FF",
                 "20"= "#D874FD",
                 "21"= "#FF6B94")) +
  ggtitle("") +
  xlab("UMAP1") +
  ylab("UMAP2") +
  theme(legend.text = element_text(size = 16*96/72),
        axis.text.x = element_text(color = "black", family = "Arial",size = 14*96/72),
        axis.text.y = element_text(color = "black", family = "Arial",size = 14*96/72),
        axis.title.x = element_text(color = "black", family = "Arial",size = 16*96/72),
        axis.title.y = element_text(color = "black", family = "Arial",size = 16*96/72),
        aspect.ratio = 1,
        panel.background = element_blank()) +
  guides(color = guide_legend(override.aes = list(size=4)))


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
DimPlot(object = runx1,
        reduction = "umap",
        group.by = "RNA_snn_res.0.6",
        label = T,
        repel = T,
        pt.size = 0.8,
        label.size = 8)


# Fig4B ----------------------------------------------------------------

## Monocle3
Idents(runx1) <- "RNA_snn_res.0.6"
start = WhichCells(runx1,idents = "0")
expression_matrix <- runx1@assays[["RNA"]]@counts
cell_metadata <- runx1@meta.data
cell_metadata$orig.ident <- factor(cell_metadata$orig.ident, levels = unique(cell_metadata$orig.ident))
gene_annotation <- data.frame("gene_short_name" = rownames(runx1))
rownames(gene_annotation) <- gene_annotation$gene_short_name
cds <- new_cell_data_set(expression_data = expression_matrix, 
                         cell_metadata = cell_metadata, 
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds)
## Construct and assign the made up partition
recreate.partition <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition
## Assign the cluster info
list_cluster <- runx1@meta.data[["RNA_snn_res.0.6"]]
names(list_cluster) <- runx1@assays[["RNA"]]@data@Dimnames[[2]]
cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
cds@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
## Assign UMAP coordinate
reducedDims(cds)[["UMAP"]] <- runx1@reductions[["umap"]]@cell.embeddings
### Reset colnames and rownames (consequence of UMAP replacement)
rownames(cds@principal_graph_aux[['UMAP']]$dp_mst) <- NULL
## Learn Graph
cds <- learn_graph(cds = cds,use_partition = T,learn_graph_control=list(ncenter=220,minimal_branch_len=15),verbose = T)

plot_cells(cds,
           color_cells_by = "RNA_snn_res.0.6",
           label_groups_by_cluster=FALSE,
           label_leaves=F,
           label_branch_points=F,
           graph_label_size=1.5,
           group_label_size = 6,
           cell_size = 1)


# Fig4C -------------------------------------------------------------------

## differential markers
degs_clu11_vs_clu0_1_2 <- FindMarkers(object = runx1, ident.1 = 11, ident.2 = c(0,1,2), only.pos = FALSE, min.pct = 0, test.use = "wilcox", logfc.threshold = 0, min.cells.group = 0)
degs_clu11_vs_clu16_17 <- FindMarkers(object = runx1, ident.1 = 11, ident.2 = c(16,17), only.pos = FALSE, min.pct = 0, test.use = "wilcox", logfc.threshold = 0, min.cells.group = 0)

## GSEA
adj.pval.thr <- 0.05
gmt.obj <- read.gmt("Scarfo_HEC2023/scRNAseq/c5.all.v7.2.symbols.gmt")
comp <- c("clu11_vs_clu0_1_2", "clu11_vs_clu16_17")
for (i in comp) {
  degs <- get(paste("degs",i,sep="_"))
  gene.res.sorted <- degs[order(degs$avg_logFC, decreasing = TRUE), ]
  # LogFC with Symbols
  gene.res.logfc.symbol <- as.vector(gene.res.sorted$avg_logFC)
  names(gene.res.logfc.symbol) <- row.names(gene.res.sorted)
  contr.gsea <- GSEA(gene.res.logfc.symbol, 
                     TERM2GENE = gmt.obj,
                     nPerm = 10000, 
                     pvalueCutoff = adj.pval.thr)
  contr.gsea.res <- as.data.frame(contr.gsea)
  assign(paste0("contr.gsea.res_",i), contr.gsea.res)
  if (i == "clu11_vs_clu0_1_2") {
    terms <- c("GO_CELL_CELL_JUNCTION_ORGANIZATION",
               "GO_ACTIN_CYTOSKELETON",
               "GO_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT",
               "GO_STRUCTURAL_CONSTITUENT_OF_RIBOSOME",
               "GO_RIBOSOME",
               "GO_TRANSLATIONAL_INITIATION")
  } else {
    terms <- c("GO_MITOTIC_NUCLEAR_DIVISION",
               "GO_CELL_DIVISION",
               "GO_MITOTIC_CELL_CYCLE",
               "GO_REGULATION_OF_NOTCH_SIGNALING_PATHWAY",
               "GO_CYTOSOLIC_RIBOSOME",
               "GO_ENDOTHELIAL_CELL_MIGRATION")
  }
  ti <- contr.gsea.res[c(contr.gsea.res$Description %in% terms), ]
  ti <- mutate(ti, color = case_when(ti$NES > 0 ~ "Upregulated",
                                     ti$NES < 0 ~ "Downregulated"))
  ti <- arrange(ti, NES)
  ti$Description <- factor(ti$Description, levels = ti$Description)
  p <- ggplot(data=ti, aes(x=NES,  y=Description, fill = color)) +
    geom_bar(stat="identity", width = 0.7) +
    scale_fill_manual(values = c("Upregulated" = "red","Downregulated" = "blue"),name = "") +
    ylab("Terms") +
    ggtitle(paste(i)) +
    theme_bw() +
    theme(
      axis.title.x = element_text(colour="black",family = "Arial", size=12*96/72),
      axis.text.x = element_text(colour="black",family = "Arial", size=10*96/72),
      axis.title.y = element_text(colour="black",family = "Arial", size=12*96/72),
      axis.text.y = element_text(colour="black",family = "Arial", size=8*96/72),
      aspect.ratio = 0.8)
  assign(paste0("p_",i), p)
}

# Fig4D -------------------------------------------------------------------

runx1@meta.data[["Phase"]] <- recode_factor(runx1@meta.data[["Phase"]],
                                            "G1" = "G1", 
                                            "S" = "S", 
                                            "G2M" = "G2/M")
## grouping clusters 0,1,2 and cluster 16,17
runx1@meta.data[["RNA_snn_res.0.6"]] <- recode_factor(runx1@meta.data[["RNA_snn_res.0.6"]],
                                                          "0" = "0,1,2", 
                                                          "1" = "0,1,2", 
                                                          "2" = "0,1,2",
                                                          "11" = "11",
                                                          "16" = "16,17",
                                                          "17" = "16,17")
cc_umap <- DimPlot(object = runx1,
                   reduction = "umap",
                   group.by = "Phase",
                   split.by = "RNA_snn_res.0.6",
                   label = F,
                   repel = T,
                   pt.size = 0.8,
                   label.size = 8,
                   cols = c("G1" = "#E6AB02", 
                            "S" = "#7570B3", 
                            "G2/M" = "#66A61E")) +
  ggtitle("") +
  xlab("UMAP1") +
  ylab("UMAP2") +
  theme(legend.text = element_text(size = 12*96/72),
        axis.text.x = element_text(color = "black", family = "Arial",size = 10*96/72),
        axis.text.y = element_text(color = "black", family = "Arial",size = 10*96/72),
        axis.title.x = element_text(color = "black", family = "Arial",size = 12*96/72),
        axis.title.y = element_text(color = "black", family = "Arial",size = 12*96/72),
        aspect.ratio = 1,
        plot.margin = margin(0, 0, 0, 0, "pt"),
        panel.background = element_blank()) +
  guides(color = guide_legend(override.aes = list(size=4)))

## average cell counts per cluster
ggp_obj_clu <- melt(table(runx1$Phase, runx1$RNA_snn_res.0.6))
colnames(ggp_obj_clu) <- c("Phase","Cluster", "Count")
cc.clusters <- c("0,1,2","11","16,17")
ggp_obj_clu <- ggp_obj_clu[c(ggp_obj_clu$Cluster %in% cc.clusters), ]
for (i in cc.clusters) {
  totalB <- subset(ggp_obj_clu, ggp_obj_clu$Cluster==i)
  totalB <- totalB %>% 
    mutate(percentage = Count/sum(Count))
  totalB$percentage <- totalB$percentage*100
  totalB$fraction = totalB$Count / sum(totalB$Count)
  # Compute the cumulative percentages
  totalB$ymax <- cumsum(totalB$fraction)
  totalB$ymin <- c(0, head(totalB$ymax, n=-1))
  totalB$labelPosition <- (totalB$ymax + totalB$ymin) / 2
  # Compute a good label
  totalB$label <- paste0(totalB$Phase, ":","\n",round(totalB$percentage, 2), "%")
  assign(paste0("totalB_clu_", i), totalB)
  p <- ggplot(totalB, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Phase)) +
    geom_rect() +
    geom_text(x=5, aes(y=labelPosition, label=label), size=5) +
    #scale_fill_brewer(palette = "Set1") +
    scale_fill_manual(values = c("G1" = "#E6AB02", 
                                 "S" = "#7570B3", 
                                 "G2/M" = "#66A61E")) +
    coord_polar(theta="y") +
    xlim(c(2, 5)) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0, "pt")) +
    theme(legend.position = "none")
  assign(paste0("p_clu_",i), p)
}
p_final <- cc_umap / (`p_clu_0,1,2` | p_clu_11 | `p_clu_16,17`)

# Fig4E -------------------------------------------------------------------

cds <- order_cells(cds, root_cells = start)
cds <- estimate_size_factors(cds)
gi <- c("FCGR2B","HES1","HEY1", "HEY2")
cds_gi <- cds[rowData(cds)$gene_short_name %in% gi,]
plot_genes_in_pseudotime(cds_subset = cds_gi, 
                         ncol = 2,
                         color_cells_by = "RNA_snn_res.0.6", 
                         min_expr = NULL,
                         cell_size = 1)







