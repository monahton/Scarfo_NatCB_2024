suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(outliers))
suppressPackageStartupMessages(library(R.utils))


# ED_Fig5B ----------------------------------------------------------------

## load scRNAseq data from public dataset GSE162950
load("Scarfo_HEC2023/scRNAseq/sample_final.Rdata")
sample <- SetIdent(sample, value = sample@meta.data$seurat_clusters)
DimPlot(sample, label=TRUE, group.by = "seurat_clusters")
sample@meta.data[["orig.ident"]] <- factor(sample@meta.data[["orig.ident"]], 
                                           levels = c("AGM_CS10_Liu",
                                                      "AGM_CS11_Liu",
                                                      "AGM_CS13_Liu",
                                                      "AGM_4wk_658",
                                                      "AGM_5wk_555",
                                                      "AGM_5wk_575",
                                                      "AGM_6wk_563",
                                                      "YolkSac_CS11_Liu"))
DimPlot(sample, 
        label=F,
        repel = T,
        reduction="umap", 
        group.by="orig.ident",
        pt.size=0.1,
        cols = c("AGM_CS10_Liu" = "#00BFC4",
                 "AGM_CS11_Liu" = "#00A9FF",
                 "AGM_CS13_Liu" = "#C77CFF",
                 "AGM_4wk_658" = "#F8766D",
                 "AGM_5wk_555" = "#CD9600",
                 "AGM_5wk_575" = "#7CAE00",
                 "AGM_6wk_563" = "#00BE67",
                 "YolkSac_CS11_Liu" = "#FF61CC"))


# ED_Fig5C ----------------------------------------------------------------

pal <- viridis(n = 10, option = "D")
FeaturePlot(sample,
            features = c("HOXA9", "HOXA10"),
            ncol = 2,
            repel = T,
            pt.size = 0.05,
            order = T,
            cols = pal)


# ED_Fig5E ----------------------------------------------------------------

devtools::install_github("sturgeonlab/Luff-etal-2021/SingleR/versioncontrolscripts")
library(versioncontrolscripts)

## load reference dataset
cd32_dll4 <- read.table("Scarfo_HEC2023/BulkRNAseq_2/featureCounts_results_tpm.txt", 
                        header=T, sep="\t")

## load Seurat data
load("Scarfo_HEC2023/scRNAseq/sample_final.Rdata")
sample <- SetIdent(sample, value = sample@meta.data$seurat_clusters)
sample_noYS <- subset(sample, orig.ident != "YolkSac_CS11_Liu")
sample_noYS@meta.data[["orig.ident"]] <-
  factor(
    sample_noYS@meta.data[["orig.ident"]],
    levels = c(
      "AGM_CS10_Liu",
      "AGM_CS11_Liu",
      "AGM_CS13_Liu",
      "AGM_4wk_658",
      "AGM_5wk_555",
      "AGM_5wk_575",
      "AGM_6wk_563"
    )
  )

## defining HEC and AEC
aec <- WhichCells(sample_noYS, 
                 expression = CDH5 > 0.2 & 
                   CXCR4 > 0.3 & 
                   GJA5 > 0.2 & 
                   DLL4 > 0.1 & 
                   PTPRC < 0.1 & 
                   SPN < 0.05 & 
                   HEY2 > 0.15)
hec <- WhichCells(sample_noYS, 
                  expression =  RUNX1 > 0.1 & 
                    CDH5 > 0.2 & 
                    HOXA9 > 0.05 & 
                    PTPRC < 0.1 & 
                    ITGA2B < 0.1 & 
                    SPN < 0.05)
Idents(sample_noYS, cells=aec) = "AEC"
Idents(sample_noYS, cells=hec) = "HEC"
sample_noYS$CellType <- Idents(sample_noYS)
sample_noYS <- subset(sample_noYS, CellType=="AEC" | CellType=="HEC")
Idents(sample_noYS) <- "CellType"
## setting up reference dataset
bulk <- cd32_dll4[,-c(2:6)]
bulk <- bulk[,c(1,2,3,4,6,7,9)]
colnames(bulk) <- c("genes",
                    "RS227_CD184","RS227_CD32_plus",
                    "RS242_CD184","RS242_CD32_plus",
                    "RS243_CD184","RS243_CD32_plus")
bulk$CD184 <- rowMeans(bulk[,c(2,4,6)])
bulk$CD32_plus <- rowMeans(bulk[,c(3,5,7)])
bulk <- bulk[,c(1,9,8)]
types = list("CD32_plus","CD184")
main_types = list("CD32_plus","CD184")
x = bulk[,1]
loss = bulk[,-1]
rownames(loss) = make.names(x, unique = TRUE)
bulk <- loss
bulk = as.matrix(bulk)
types = as.character(types)
main_types = as.character(main_types)
myref = list(data = bulk, main_types = main_types, types = types)
myref[["de.genes"]] = CreateVariableGeneSet(bulk,types,200)
myref[["de.genes.main"]] = CreateVariableGeneSet(bulk,main_types,300)
labels <- sample_noYS@meta.data$orig.ident
stage <- sample_noYS@meta.data$CellType
stage2 <- data.frame(stage)
stage2$stage <- as.character(stage2$stage)
stage2 <- stage2$stage
# create singlecell dataset for comparison to bulk RNAseq
data <- GetAssayData(sample_noYS, assay = "RNA", slot = "data")
data <- as.matrix(data)
# run main SingleR function & add names for cluster
obj <- SingleR.CreateObject(data, 
                            myref, 
                            clusters = NULL, 
                            variable.genes = "de", 
                            fine.tune = F)
obj[["names"]] <- stage2
obj[["names2"]] <- labels
data <- t(obj$SingleR.single.main$r)
cell_types = as.data.frame(obj$names)
colnames(cell_types) = 'CellTypes'
rownames(cell_types) = colnames(data)
test <- as.data.frame(obj$names2)
colnames(test) <- "Samples"
rownames(test) <- colnames(data)
clusters <- cbind(cell_types, test)
breaksList = seq(0, 1, by = 0.01)
# export relative correlation scores and names for each cell in cs dataset
scores = obj$SingleR.single.main$r
datascores <- as.matrix(scores)
scores <- data.frame(datascores)
type <- data.frame(obj$names)
cells <- data.frame(obj$names2)
scores$type <- type$obj.names
scores$cells <- cells$obj.names2

#successive sorting orders cells in aesthetically-pleasing manner
for (k in unique(scores$cells)) {
  clu <- filter(scores, cells == k)
  clu <- clu[order(clu$CD184), ]
  clu <- clu[order(clu$CD32_plus), ]
  clu_cells <- rownames(clu)
  assign(paste0("clu",k,"_cells"),clu_cells)
}

## reorder both clusters and HEC and AEC cells
agm_4 <- filter(scores, cells == "AGM_4wk_658")
agm_4_hec <- filter(agm_4, type == "HEC")
agm_4_hec <- agm_4_hec[order(agm_4_hec$CD184),]
agm_4_hec <- agm_4_hec[order(agm_4_hec$CD32_plus),]
agm_4_aec <- filter(agm_4, type == "AEC")
agm_4_aec <- agm_4_aec[order(agm_4_aec$CD184),]
agm_4_aec <- agm_4_aec[order(agm_4_aec$CD32_plus),]

agm_5_1 <- filter(scores, cells == "AGM_5wk_555")
agm_5_1_hec <- filter(agm_5_1, type == "HEC")
agm_5_1_hec <- agm_5_1_hec[order(agm_5_1_hec$CD184),]
agm_5_1_hec <- agm_5_1_hec[order(agm_5_1_hec$CD32_plus),]
agm_5_1_aec <- filter(agm_5_1, type == "AEC")
agm_5_1_aec <- agm_5_1_aec[order(agm_5_1_aec$CD184),]
agm_5_1_aec <- agm_5_1_aec[order(agm_5_1_aec$CD32_plus),]

agm_5_2 <- filter(scores, cells == "AGM_5wk_575")
agm_5_2_hec <- filter(agm_5_2, type == "HEC")
agm_5_2_hec <- agm_5_2_hec[order(agm_5_2_hec$CD184),]
agm_5_2_hec <- agm_5_2_hec[order(agm_5_2_hec$CD32_plus),]
agm_5_2_aec <- filter(agm_5_2, type == "AEC")
agm_5_2_aec <- agm_5_2_aec[order(agm_5_2_aec$CD184),]
agm_5_2_aec <- agm_5_2_aec[order(agm_5_2_aec$CD32_plus),]

agm_6 <- filter(scores, cells == "AGM_6wk_563")
agm_6_hec <- filter(agm_6, type == "HEC")
agm_6_hec <- agm_6_hec[order(agm_6_hec$CD184),]
agm_6_hec <- agm_6_hec[order(agm_6_hec$CD32_plus),]
agm_6_aec <- filter(agm_6, type == "AEC")
agm_6_aec <- agm_6_aec[order(agm_6_aec$CD184),]
agm_6_aec <- agm_6_aec[order(agm_6_aec$CD32_plus),]

agm_cs10 <- filter(scores, cells == "AGM_CS10_Liu")
agm_cs10_hec <- filter(agm_cs10, type == "HEC")
agm_cs10_hec <- agm_cs10_hec[order(agm_cs10_hec$CD184),]
agm_cs10_hec <- agm_cs10_hec[order(agm_cs10_hec$CD32_plus),]
agm_cs10_aec <- filter(agm_cs10, type == "AEC")
agm_cs10_aec <- agm_cs10_aec[order(agm_cs10_aec$CD184),]
agm_cs10_aec <- agm_cs10_aec[order(agm_cs10_aec$CD32_plus),]

agm_cs11 <- filter(scores, cells == "AGM_CS11_Liu")
agm_cs11_hec <- filter(agm_cs11, type == "HEC")
agm_cs11_hec <- agm_cs11_hec[order(agm_cs11_hec$CD184),]
agm_cs11_hec <- agm_cs11_hec[order(agm_cs11_hec$CD32_plus),]
agm_cs11_aec <- filter(agm_cs11, type == "AEC")
agm_cs11_aec <- agm_cs11_aec[order(agm_cs11_aec$CD184),]
agm_cs11_aec <- agm_cs11_aec[order(agm_cs11_aec$CD32_plus),]

agm_cs13 <- filter(scores, cells == "AGM_CS13_Liu")
agm_cs13_hec <- filter(agm_cs13, type == "HEC")
agm_cs13_hec <- agm_cs13_hec[order(agm_cs13_hec$CD184),]
agm_cs13_hec <- agm_cs13_hec[order(agm_cs13_hec$CD32_plus),]
agm_cs13_aec <- filter(agm_cs13, type == "AEC")
agm_cs13_aec <- agm_cs13_aec[order(agm_cs13_aec$CD184),]
agm_cs13_aec <- agm_cs13_aec[order(agm_cs13_aec$CD32_plus),]

agm_4_hec.cells <- row.names(agm_4_hec)
agm_4_aec.cells <- row.names(agm_4_aec)
agm_5_1_hec.cells <- row.names(agm_5_1_hec)
agm_5_1_aec.cells <- row.names(agm_5_1_aec)
agm_5_2_hec.cells <- row.names(agm_5_2_hec)
agm_5_2_aec.cells <- row.names(agm_5_2_aec)
agm_6_hec.cells <- row.names(agm_6_hec)
agm_6_aec.cells <- row.names(agm_6_aec)
agm_cs10_hec.cells <- row.names(agm_cs10_hec)
agm_cs10_aec.cells <- row.names(agm_cs10_aec)
agm_cs11_hec.cells <- row.names(agm_cs11_hec)
agm_cs11_aec.cells <- row.names(agm_cs11_aec)
agm_cs13_hec.cells <- row.names(agm_cs13_hec)
agm_cs13_aec.cells <- row.names(agm_cs13_aec)

order <- c(
  agm_cs10_aec.cells,
  agm_cs11_aec.cells,
  agm_cs13_aec.cells,
  agm_4_aec.cells,
  agm_5_1_aec.cells,
  agm_5_2_aec.cells,
  agm_6_aec.cells,
  agm_cs10_hec.cells,
  agm_cs11_hec.cells,
  agm_cs13_hec.cells,
  agm_4_hec.cells,
  agm_5_1_hec.cells,
  agm_5_2_hec.cells,
  agm_6_hec.cells
)
data_reorder <- data[,match(order, colnames(data))]
annot_colors <- list(CellTypes=c(HEC = "red",
                                 AEC = "royalblue3"),
                     Samples = c(AGM_CS10_Liu = "#00BFC4",
                                 AGM_CS11_Liu = "#00A9FF",
                                 AGM_CS13_Liu = "#C77CFF",
                                 AGM_4wk_658 = "#F8766D",
                                 AGM_5wk_555 = "#CD9600",
                                 AGM_5wk_575 = "#7CAE00",
                                 AGM_6wk_563 = "#00BE67"))
annot <- clusters
colours <- annot_colors
colAnn <- HeatmapAnnotation(df = clusters,
                            which = 'col',
                            show_annotation_name = T,
                            col = annot_colors,
                            annotation_legend_param = list(
                              CellTypes = list(
                                title = "Cell Type",
                                title_gp = gpar(fontsize = 14*96/72,
                                                fontface = "bold"), 
                                labels_gp = gpar(fontsize = 14*96/72),
                                grid_height = unit(0.7, "cm")),
                              Samples = list(
                                title = "Samples",
                                title_gp = gpar(fontsize = 14*96/72,
                                                fontface = "bold"), 
                                labels_gp = gpar(fontsize = 14*96/72),
                                grid_height = unit(0.7, "cm"))),
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))
counts_scaled <- t(scale(t(data_reorder)))
counts_scaled_ordered <- counts_scaled[, match(rownames(clusters),colnames(counts_scaled))]
clusters$DF <- paste(clusters$CellTypes, clusters$Samples, sep = "_")
dd <- cluster_within_group(counts_scaled_ordered, clusters$DF)
hmap <- Heatmap(counts_scaled_ordered,
                name = "Z-score",
                show_column_names = F,
                show_column_dend = F,
                show_parent_dend_line = F,
                show_row_dend = F,
                cluster_rows = F,
                cluster_columns = dd,
                column_gap = unit(c(2,2,5,2,5,2,5,2,
                                    2,5,5,2,5), "mm"),
                top_annotation=colAnn,
                column_split = 14,
                clustering_method_columns = "ward.D2",
                heatmap_legend_param = list(
                  title = "Z-Scores",
                  title_gp = gpar(fontsize = 14*96/72,
                                  fontface = "bold"), 
                  labels_gp = gpar(fontsize = 12*96/72),
                  legend_height = unit(2, "cm")
                ),
                width = unit(15*96/72, "cm"), 
                height = unit(5*96/72, "cm"),
                use_raster = TRUE, # for raster quality, needs to be tested
                raster_quality = 5)
draw(hmap, heatmap_legend_side="right", annotation_legend_side="right")



# ED_Fig5F ----------------------------------------------------------------

## load scRNAseq data from public dataset GSE162950
load("Scarfo_HEC2023/scRNAseq/seurat_object.Rdata")
EHT_scorecard <- c("CXCR4","SOX17","GJA5","MECOM","HOXA9","SPN","CD44","ITGA2B","HLF","GFI1","MLLT3","KCNK17","MYB","STAT5A","SMIM24","RAB27B","SPINK2")
## Select HE cells
sample <- SetIdent(sample, value = sample@meta.data$seurat_clusters)
index1 = GetAssayData(sample, slot="counts")["RUNX1",]>0
index2 = GetAssayData(sample, slot="counts")["CDH5",]>0
index3 = GetAssayData(sample, slot="counts")["PTPRC",]==0
index4 = GetAssayData(sample, slot="counts")["FCGR2B",]==0
bars = colnames(sample)[Idents(sample) %in% c(0,3)]
bars1 = colnames(sample)[index1 & index2 & index3 & index4]
bars2 = intersect(bars, bars1)
## Establish identity for FCGR2B- HE cells
Idents(sample, cells=bars2) = "FCGR2B=0"
index5 = GetAssayData(sample, slot="counts")["FCGR2B",]>0
bars3 = colnames(sample)[index1 & index2 & index3 & index5]
bars4 = intersect(bars, bars3)
## Establish identity for FCGR2B+ HE cells
Idents(sample, cells=bars4) = "FCGR2B_exp"
sample$CellType <- Idents(sample)
sample_sub <- subset(sample, CellType=="FCGR2B_exp" | CellType=="FCGR2B=0")
DotPlot(sample_sub, 
        features=EHT_scorecard, 
        cols=c("grey90","red3"),
        dot.scale = 10) +
  coord_flip()








