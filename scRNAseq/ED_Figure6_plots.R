suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(dyno))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(patchwork))

# ED_Fig6A ------------------------------------------------------------------

## load Seurat data
sample <- readRDS("Scarfo_HEC2023/scRNAseq/WNTd_HE.rds")
sample@meta.data[["Cell.Annotations"]] <- recode_factor(sample@meta.data[["RNA_snn_res.0.6"]],
                                                        "0" = "HECs", 
                                                        "1" = "HECs", 
                                                        "2" = "HECs", 
                                                        "3" = "Venous cells", 
                                                        "4" = "Venous cells", 
                                                        "5" = "Arterial ECs", 
                                                        "6" = "Arterial ECs", 
                                                        "7" = "Venous cells", 
                                                        "8" = "ECs (M phase)",
                                                        "9" = "Venous cells",
                                                        "10" = "EndoMT cells",
                                                        "11" = "HECs",
                                                        "12" = "Arterial ECs",
                                                        "13" = "Lymphatic ECs",
                                                        "14" = "Arterial ECs",
                                                        "15" = "Allantois/Placenta",
                                                        "16" = "HECs",
                                                        "17" = "HECs",
                                                        "18" = "EndoMT cells",
                                                        "19" = "SM cells",
                                                        "20" = "EndoMT cells",
                                                        "21" = "EndoMT cells")

DimPlot(object = sample,
        reduction = "umap",
        group.by = "Cell.Annotations",
        label = F,
        repel = T,
        pt.size = 0.8,
        label.size = 6,
        cols = c('HECs' = "#377EB8", 
                 'Lymphatic ECs' = "#4DAF4A",
                 'Arterial ECs' = "#E41A1C",
                 'Venous cells' = "#FF7F00",
                 'EndoMT cells' = "#E6AB02",
                 'Allantois/Placenta' = "#984EA3",
                 'SM cells' = "#A6761D",
                 "ECs (M phase)" = "#A50026"))



# ED_Fig6B - ED_Fig6C ------------------------------------------------------------------
FeaturePlot(object = sample,
            features = c("RUNX1", "FCGR2B"),
            cols = c("grey90","blue"),
            order = T,
            pt.size = 0.1)


# ED_Fig6D ------------------------------------------------------------------

runx1 <- readRDS("Scarfo_HEC2023/scRNAseq/RUNX1_clusters.rds")

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
cds <- order_cells(cds, root_cells = start)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=T,
           label_branch_points=T,
           group_label_size = 8,
           graph_label_size = 3,
           cell_size = 1,
           trajectory_graph_segment_size = 1)



# ED_Fig6E ------------------------------------------------------------------

cds <- estimate_size_factors(cds)
gi <- c("H19","KCNK17","RUNX1", "MYB", "SPN")
cds_gi <- cds[rowData(cds)$gene_short_name %in% gi,]
plot_genes_in_pseudotime(cds_subset = cds_gi, 
                         ncol = 3,
                         color_cells_by = "RNA_snn_res.0.6", 
                         min_expr = NULL,
                         cell_size = 1)



# ED_Fig6F - ED_Fig6G ------------------------------------------------------------------

## load Seurat data
sample <- readRDS("Scarfo_HEC2023/scRNAseq/WNTd_HE.rds")
runx1 <- subset(sample, idents = c("0","1","2","11","16","17"))
Idents(runx1) <- "RNA_snn_new_res.0.6"
# preparing data for dyno
df <- as.matrix(runx1[["RNA"]]@data)
counts <- Matrix::t(as(as.matrix(runx1@assays$RNA@counts), 'sparseMatrix'))
expression <- Matrix::t(as(as.matrix(runx1@assays$RNA@data), 'sparseMatrix'))
dataset_runx1 <- wrap_expression(expression = expression,
                                 counts = counts)
dataset_runx1 <- dataset_runx1 %>% 
  add_prior_information(start_id = "Sample1_ATCTCTATCCAATCCC-1", # from cluster 0
                        start_n = 0) %>%
  add_grouping(grouping = runx1$RNA_snn_new_res.0.6) %>%
  add_dimred(Embeddings(runx1, "umap"))
# selecting method
guidelines <- guidelines_shiny(dataset_runx1)
methods_selected <- guidelines$methods_selected
methods_selected
# running top model
model_paga_t <- infer_trajectory(dataset_runx1, 
                                 method = ti_paga_tree(), 
                                 verbose = T,
                                 give_priors = c("start_id","start_n"))
# plotting trajectory and pseudotime scores
patchwork::wrap_plots(
  plot_dimred(model_paga_t, size_cells = 0.5, grouping = dataset_runx1$grouping) + 
    ggtitle(" RUNX1 Clusters") +
    scale_color_manual(values= c("#F8766D","#B79F00","#00BFC4","#619CFF",
                                 "#F564E3","#00BA38")) +
    labs(color = "Clusters") +
    guides(color = guide_legend(override.aes = list(size = 2))),
  plot_dimred(model_paga_t, size_cells = 0.5, "pseudotime", pseudotime = calculate_pseudotime(model_paga_t)) + 
    ggtitle("Pseudotime")
)














