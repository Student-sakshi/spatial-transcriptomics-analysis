
# Sakshi Parate - 7073022

setwd("/Users/sam/Documents/UdS/Sem3/SCB/Assignments/Project3/project_3_dataset") #set working directry

Sys.setenv(GITHUB_PAT = "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")
Sys.getenv("GITHUB_PAT")

install.packages("BiocManager")
BiocManager::install(c("Biobase", "BiocGenerics"))
install.packages("devtools")
install.packages("hdf5r")
devtools::install_github("satijalab/seurat")
devtools::install_github("sqjin/CellChat")
library(CellChat)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

#Week 1

#1 Spatial Transcriptomics Technology

list.files("/Users/sam/Documents/UdS/Sem3/SCB/Assignments/Project3/project_3_dataset")
list.files("Section_1")

#2 Spatial Transcriptomics Data in Seurat 

#2.1 Data loading and structure

brain <- Load10X_Spatial(data.dir = "Section_1", filename = "V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix.h5") #load data

#2.2 Object inspection and modalities

GetAssayData(brain, assay = "Spatial", layer = "counts") #row counts stored here
GetAssayData(brain, assay = "Spatial", layer = "data") #normalizd data stored here or SCTransform
GetTissueCoordinates(brain) #spatial coordinates stored her
brain@images$slice1@image #image data stored here

#2.3 Spatial context of gene expression

brain <- NormalizeData(brain, assay = "Spatial") #normalize data
counts <- GetAssayData(brain, assay = "Spatial", layer = "counts") #extract raw count matrx
gene_detection_rate <- rowSums(counts > 0) / ncol(counts) 
broad_gene <- names(sort(gene_detection_rate, decreasing = TRUE))[1] #select broad and restricted
restricted_gene <- names(sort(gene_detection_rate))[1]
SpatialFeaturePlot(brain,features = c(broad_gene, restricted_gene)) #spatial exp

#2.4 Tissue coverage assessment

tissue_pos <- read.csv("Section_1/spatial/tissue_positions_list.csv",header = FALSE) #read spaceranger file
colnames(tissue_pos) <- c("barcode", "tissue", "row", "col", "imagerow", "imagecol") #assign col names
brain$tissue <- tissue_pos$tissue[match(colnames(brain), tissue_pos$barcode)] #atach tissue labels to seurat obj
SpatialDimPlot(brain, group.by = "tissue") #tissue vs non tissue 
total_spots <- length(brain$tissue)
non_tissue_spots <- sum(brain$tissue == 0)
proportion_removed <- non_tissue_spots / total_spots
proportion_removed #filtering

#3 Data Preprocessing

#3.1 Filtering strategy and sensitivity

head(brain@meta.data)
brain_relaxed <- subset(brain, subset = nFeature_Spatial > 200) #relaxed
brain_stringent <- subset(brain, subset = nFeature_Spatial > 1500) #stringent
ncol(brain_relaxed) #num of retained spots
ncol(brain_stringent)
nrow(brain_relaxed) #num of retained genes
nrow(brain_stringent)
SpatialDimPlot(brain_relaxed) + ggtitle("Relaxed filtering (>200 genes)")
SpatialDimPlot(brain_stringent) + ggtitle("Stringent filtering (>1500 genes)")

#3.2 Normalization with SCTransform

BiocManager::install('glmGamPoi')
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE) #normalisation by sctranform

#3.3 Normalization comparison 

brain_log <- NormalizeData(brain, assay = "Spatial", normalization.method = "LogNormalize") #log normal
gene <- broad_gene #select broad gene
broad_gene #Bc1
log_expr <- GetAssayData(brain_log, assay = "Spatial", layer = "data")[gene, ] #retrieve normal exp values
sct_expr <- GetAssayData(brain, assay = "SCT", layer = "data")[gene, ]
var(log_expr) #compare var
var(sct_expr)

#4 Dimensionality Reduction and Clustering 

#4.1 Dimensionality reduction

brain <- RunPCA(brain, assay = "SCT", verbose = FALSE) #pca on sct nor data
ElbowPlot(brain, ndims = 50) #var of first 50 PC
brain <- RunUMAP(brain, dims = 1:30)
DimPlot(brain, reduction = "umap")

#4.2 Clustering

brain <- FindNeighbors(brain, dims = 1:30) #nearest neighbor graph using pca
brain <- FindClusters(brain, resolution = 0.5) #graph based clustering
DimPlot(brain, reduction = "umap", label = TRUE) #clusters on umap
SpatialDimPlot(brain, group.by = "seurat_clusters", label = TRUE) #clusters on tisue slide

#Week 2: Differential Expression & Comparative Integration

#5 Differential Expression Analysis

#5.1 Cluster marker robustness

markers_wilcox <- FindAllMarkers(brain,assay = "SCT",test.use = "wilcox",only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25) #wilcoxon ranksum test
markers_negbinom <- FindAllMarkers(brain,assay = "SCT",test.use = "negbinom",only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25) #negative binomial test
library(dplyr)
top_wilcox <- markers_wilcox %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 20) #select top markers
top_negbinom <- markers_negbinom %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 20)
overlap_summary <- top_wilcox %>% select(cluster, gene) %>% inner_join(top_negbinom %>% select(cluster, gene), by = c("cluster", "gene")) %>% count(cluster, name = "overlap_genes") #quantify marker overlap
overlap_summary

#5.2 Spatial coherence of gene expression

DefaultAssay(brain) <- "Spatial"
brain <- NormalizeData(brain, assay = "Spatial") #log norm
brain <- ScaleData(brain, assay = "Spatial")
install.packages("Rfast2")
brain <- FindSpatiallyVariableFeatures( brain,assay = "Spatial",layer = "scale.data",selection.method = "moransi") #compute moransi
spatial_genes <- SpatiallyVariableFeatures(brain) #spatially variable genes
head(spatial_genes) #top spatially variable genes
top_spatial_genes <- spatial_genes[1:3] #select top 3 
top_spatial_genes
SpatialFeaturePlot(brain, features = top_spatial_genes) #spatial exp

#6 Integration of Multiple Samples

#6.1 Merging without batch correction

brain_1 <- Load10X_Spatial( data.dir = "Section_1",filename = "V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix.h5") #load sec 1
counts_2 <- Read10X_h5("Section_2/V1_Mouse_Brain_Sagittal_Posterior_Section_2_filtered_feature_bc_matrix.h5")
brain_2 <- CreateSeuratObject(counts = counts_2, project = "Section_2") #load sec2
brain_merged <- merge(brain_1, y = brain_2, add.cell.ids = c("Section1", "Section2"), project = "Merged_noBatch") #merge datasets without batch corection
brain_merged <- SCTransform(brain_merged, verbose = FALSE) #sctransform normalisation
brain_merged <- RunPCA(brain_merged, verbose = FALSE) #dim reduction
brain_merged <- RunUMAP(brain_merged, dims = 1:30)
brain_merged <- FindNeighbors(brain_merged, dims = 1:30) #clustering
brain_merged <- FindClusters(brain_merged, resolution = 0.5)
DimPlot(brain_merged, group.by = "orig.ident") + ggtitle("Merged samples without batch correction") 
DimPlot(brain_merged, label = TRUE) + ggtitle("Clusters in merged dataset") #umap coloured by clusters
SpatialDimPlot(brain_1, label = TRUE) + ggtitle("Spatial clusters – Section 1") #tissue plot

#6.2 Data integration and batch correction

brain_2[["Spatial"]] <- brain_2[["RNA"]]
DefaultAssay(brain_2) <- "Spatial"
brain_2[["RNA"]] <- NULL
brain_list <- list(Section1 = brain_1,Section2 = brain_2) #put datasets into list
brain_list <- lapply(brain_list, function(obj) {
  obj <- SCTransform(obj, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE);
  DefaultAssay(obj) <- "SCT";
  obj
}) #sct norm by spatial assay
features <- SelectIntegrationFeatures(object.list = brain_list,nfeatures = 3000) #select integration features only from shared genes
brain_list <- PrepSCTIntegration(object.list = brain_list,anchor.features = features) #prep for integration
anchors <- FindIntegrationAnchors(object.list = brain_list,normalization.method = "SCT",anchor.features = features) #find anchors
brain_integrated <- IntegrateData(anchorset = anchors,normalization.method = "SCT") #integrate data
brain_integrated <- RunPCA(brain_integrated, verbose = FALSE) #dim reduction
brain_integrated <- RunUMAP(brain_integrated, dims = 1:30)
brain_integrated <- FindNeighbors(brain_integrated, dims = 1:30) #clustering
brain_integrated <- FindClusters(brain_integrated, resolution = 0.5)
DimPlot(brain_integrated, group.by = "orig.ident") + ggtitle("Integrated data with batch correction") 
DimPlot(brain_integrated, label = TRUE) + ggtitle("Clusters after SCT integration")

#6.3 Assessment of batch correction effects

DimPlot(brain_merged,group.by = "orig.ident") + ggtitle("Before batch correction: merged samples") #6.1
DimPlot(brain_merged,label = TRUE) + ggtitle("Clusters before batch corrction")
DimPlot(brain_integrated,group.by = "orig.ident") + ggtitle("After batch correction: integrated samples") #6.2
DimPlot(brain_integrated,label = TRUE) + ggtitle("Clusters after batch correction")

#Week 3

#7 Cell-Type Annotation and Validation

#7.1 Label transfer with confidence assessment

ref <- readRDS("allen_cortex.rds")  #load ref
ref <- SCTransform(ref, verbose = FALSE) #sct norm
ref <- RunPCA(ref, verbose = FALSE) #dim reduction
ref <- RunUMAP(ref, dims = 1:30)
DefaultAssay(brain) <- "Spatial"
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE) #apply sct to spatial dataset
brain <- RunPCA(brain, verbose = FALSE) #pca
head(colnames(ref@meta.data)) 
colnames(ref@meta.data) 
anchors <- FindTransferAnchors(reference = ref,query = brain,normalization.method = "SCT",dims = 1:30) #find anchors betw ref n spatial
predictions <- TransferData(anchorset = anchors,refdata = ref$subclass,dims = 1:30) #transfer n compute scores
brain <- AddMetaData(brain, metadata = predictions)
head(brain@meta.data[, c("predicted.id", "prediction.score.max")])
SpatialFeaturePlot(brain,features = "prediction.score.max")
SpatialDimPlot(brain,group.by = "predicted.id",label = TRUE)

#7.2 Marker-based validation of annotations

DefaultAssay(brain) <- "SCT"
#Astrocytes
p1 <- SpatialFeaturePlot(brain,features = "prediction.score.max") + ggtitle("Astro: prediction confidence") 
p2 <- SpatialFeaturePlot(brain,features = "Aldh1l1") + ggtitle("Astro marker: Aldh1l1") 
p3 <- SpatialFeaturePlot(brain,features = "Gfap") + ggtitle("Astro weak marker: Gfap") 
#Oligodendrocytes
p4 <- SpatialFeaturePlot(brain,features = "prediction.score.max") + ggtitle("Oligo: prediction confidence")
p5 <- SpatialFeaturePlot(brain,features = "Mbp") + ggtitle("Oligo marker: Mbp")
p6 <- SpatialFeaturePlot(brain,features = "Plp1") + ggtitle("Oligo weak marker: Plp1")
#L2/3 IT neurons
p7 <- SpatialFeaturePlot(brain,features = "prediction.score.max") + ggtitle("L2/3 IT: prediction confidence")
p8 <- SpatialFeaturePlot(brain,features = "Reln") + ggtitle("L2/3 IT marker: Reln")
p9 <- SpatialFeaturePlot(brain,features = "Cux2") + ggtitle("L2/3 IT weak marker: Cux2")
p1 + p2
p1 + p3
p4 + p5
p4 + p6
p7 + p8
p7 + p9

#8 Spot Deconvolution and Consistency Analysis 

#8.1 Conceptual motivation

#8.2 Reference sensitivity analysis 

ref_full <- ref
ref_down <- subset(ref,cells = unlist(lapply(split(colnames(ref), ref$subclass),function(x) sample(x, min(250, length(x))))))
anchors_full <- FindTransferAnchors(reference = ref_full,query = brain,normalization.method = "SCT",dims = 1:30)
pred_full <- TransferData(anchorset = anchors_full,refdata = ref_full$subclass,dims = 1:30)
anchors_down <- FindTransferAnchors(reference = ref_down,query = brain,normalization.method = "SCT",dims = 1:30)
pred_down <- TransferData(anchorset = anchors_down,refdata = ref_down$subclass,dims = 1:30)
cor(pred_full$prediction.score.max,pred_down$prediction.score.max,use = "complete.obs")

#8.3 Concordance with annotation results

colnames(predictions)
head(colnames(predictions))
head(brain@meta.data$prediction.score.max)
oligo_deconv <- predictions[["prediction.score.Oligo"]] #deconvolution
astro_deconv <- predictions[["prediction.score.Astro"]]
label_conf <- brain@meta.data$prediction.score.max #label transfer confidence
cor_oligo <- cor(oligo_deconv,label_conf,method = "pearson",use = "complete.obs") #pearson corelation for oligo
cor_astro <- cor(astro_deconv,label_conf,method = "pearson",use = "complete.obs") #pearson correlation for astro
cor_oligo
cor_astro

#8.4 Reference perturbation test

ref_no_oligo <- subset(ref, subset = subclass != "Oligo") #rem oligo
anchors_no_oligo <- FindTransferAnchors(reference = ref_no_oligo,query = brain,normalization.method = "SCT",dims = 1:30)
pred_no_oligo <- TransferData(anchorset = anchors_no_oligo,refdata = ref_no_oligo$subclass,dims = 1:30)
cor(pred_full$prediction.score.max,pred_no_oligo$prediction.score.max,use = "complete.obs")

# Week 4

#9 Cell-cell communication 

DefaultAssay(brain) <- "SCT"
Idents(brain) <- brain$predicted.id
cellchat <- createCellChat(object = GetAssayData(brain, assay = "SCT", layer = "data"),meta = brain@meta.data,group.by = "predicted.id") #creat cellchat obj
CellChatDB <- CellChatDB.mouse
cellchat@DB <- CellChatDB #add db
cellchat <- subsetData(cellchat) #preprocessing
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
graphics.off()
options(device = "RStudioGD")
par(mar = c(1, 1, 1, 1))
netVisual_circle(cellchat@net$count,vertex.weight = as.numeric(table(cellchat@idents)),weight.scale = TRUE,label.edge = FALSE,title.name = "Number of interactions") #circle plot
par(mar = c(1, 1, 1, 1))
netVisual_circle(cellchat@net$weight,vertex.weight = as.numeric(table(cellchat@idents)),weight.scale = TRUE,label.edge = FALSE,title.name = "Interaction strength") #spatial plot

#10 Summary



