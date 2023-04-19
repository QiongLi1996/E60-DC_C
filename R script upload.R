# Single-cell data reading, pre-processing, normalization, dimensionality reduction, clustering, and definition.
library(Seurat)
library(Matrix)
library(ggplot2)
# read scRNA-seq data
data <- Read10X("../data/E60/E60 scRNA-seq/")
# create Seurat object
data <- CreateSeuratObject(
  data,
  project = "e60", 
  min.cells = 5,
  min.features = 5)
# subset cells based on quality control
data <- subset(data, subset = nFeature_RNA > 500)
# calculate mitochondrial genes percent
mt.genes=c('ND1','ND2','COX1',
           'COX2','ATP8','ATP6',
           'COX3','ND3','ND4L',
           'ND4','ND5','ND6',
           'CYTB')
head(rownames(data))
kp=mt.genes %in% rownames(data)
table(kp)
mt.genes=mt.genes[kp]
C<-GetAssayData(object = data, slot = "counts")
percent.mito <- Matrix::colSums(C[mt.genes,])/Matrix::colSums(C)*100
data <- AddMetaData(data, percent.mito, col.name = "percent_mito")
# subset cells based on mitochondrial genes percent
selected_mito <- WhichCells(data, expression = percent_mito < 10)
data <- subset(data, cells = selected_mito)
# normalize data
data <- SCTransform(data, ncells = 3000, verbose = FALSE)
# run PCA
data <- RunPCA(object = data, verbose = F,npcs = 50)
# find neighbors
dim.use <- 1:20
data <- FindNeighbors(data, dims = dim.use)
# cluster cells
data <- FindClusters(data, resolution = 1.5)
# run UMAP
data <- RunUMAP(data, dims = dim.use)
# UMAP plot
DimPlot(data, group.by = 'seurat_clusters',label = T)
# find cluster marker gene
library(dplyr) 
sce.markers <- FindAllMarkers(object = data,only.pos = TRUE,
                              min.pct = 0.25, 
                              thresh.use = 0.25)
# rename cell types
data <- RenameIdents(data, '0'='Oral Epithelial',
                     '1'='Oral Epithelial',
                     '2'='Oral Epithelial',
                     '3'='Oral Epithelial',
                     '4'='SDL',
                     '5'='Oral Epithelial',
                     '6'='Oral Epithelial',
                     '7'='Oral Epithelial',
                     '8'='Oral Epithelial',
                     '9'='Oral Epithelial',
                     '10'='Oral Epithelial',
                     '11'='Oral Epithelial',
                     '12'='Fibroblast',
                     '13'='Fibroblast',
                     '14'='Oral Epithelial',
                     '15'='Oral Epithelial',
                     '16'='Oral Epithelial',
                     '17'='Oral Epithelial',
                     '18'='Fibroblast',
                     '19'='SDL',
                     '20'='Endothelial',
                     '21'='Fibroblast',
                     '22'='Oral Epithelial',
                     '23'='Oral Epithelial',
                     '24'='Oral Epithelial',
                     '25'='Fibroblast',
                     '26'='SDL',
                     '27'='Fibroblast',
                     '28'='Oral Epithelial',
                     '29'='Oral Epithelial',
                     '30'='Myeloid',
                     '31'='SDL',
                     '32'='Fibroblast',
                     '33'='SDL')
data@meta.data$celltype <- Idents(data)
# cell types UMAP plot
library(ggsci)
DimPlot(data, reduction = 'umap', group.by = 'celltype',
        label = F, pt.size = 0.5)+scale_color_npg()
ggsave('all celltype umap.pdf',width = 7,height = 6)
# cell types marker gene heatmap
library(dplyr) 
# find cell types marker gene
sce.markers <- FindAllMarkers(object = data,only.pos = TRUE,
                              min.pct = 0.25, 
                              thresh.use = 0.25)
# Define color
color <- pal_npg(palette = c("nrc"), alpha = 1)(5)
# Define celltype order
ord = c('Oral Epithelial','SDL', 'Fibroblast', 'Endothelial','Myeloid' )
data$celltype = factor(data$celltype ,levels = ord)
# select top6 marker genes for heatmap
top6 <- sce.markers %>% group_by(cluster) %>% top_n(6, avg_log2FC)
ll = split(top6$gene,top6$cluster)
ll = ll[ord]
rmg=names(table(unlist(ll))[table(unlist(ll))>1])
ll = lapply(ll, function(x) x[!x %in% rmg])
library(ggplot2)
DoHeatmap(data,
          features = unlist(ll),
          group.by = 'celltype',
          assay = 'SCT',
          group.colors = color,label = F)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))
ggsave(filename = "all celltype marker_pheatmap.pdf",width = 8,height = 5)

# read spatial transcriptomics data
cortex = Load10X_Spatial(data.dir = "../data/E60/E60 DC+C 10X ST/E60DC_C/2.Basic_analysis/outs/",slice = "E60")
# spatial transcriptomics data normalization
cortex <- SCTransform(cortex,assay = "Spatial")
# spatial transcriptomics data reduction
cortex <- RunPCA(cortex)
cortex <- RunUMAP(cortex, dims = 1:30, label = T)
# find clusters
cortex <- FindNeighbors(cortex,dims = 1:30)
cortex <- FindClusters(cortex,resolution = 0.5)

# Integrated analysis of single-cell transcriptomics and spatial transcriptomics data.
# scRNA-seq all cell types map to spatial transcriptomics
anchors <- FindTransferAnchors(reference = data, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = data$celltype, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]], dims = 1:50)
cortex[["predictions"]] <- predictions.assay
DefaultAssay(cortex) <- "predictions"

# plot mapping results
table(data$celltype)
SpatialFeaturePlot(cortex, features = c('Oral Epithelial','SDL', 'Fibroblast', 
                                        'Endothelial','Myeloid'), 
                   pt.size.factor = 1.6,min.cutoff = 0.3,ncol = 2, 
                   crop = TRUE, alpha = c(0,1)) + scale_alpha()
ggsave('all celltypes map.mix.pdf',width = 8, height = 15)


# normalization and analysis of fibroblast subpopulation data.
cells_sub <- subset(data@meta.data, 
                    celltype %in%'Fibroblast')
scRNA_sub <- subset(data, 
                    cells=row.names(cells_sub))
scRNA_sub <- SCTransform(scRNA_sub, ncells = 3000, verbose = FALSE)%>% RunPCA(
  verbose = F,npcs = 50)%>% FindNeighbors(dims = dim.use)%>% FindClusters(resolution = 0.5)%>% RunUMAP(dims = dim.use)
# fibroblast subclusters UMAP plot
DimPlot(scRNA_sub, group.by = 'seurat_clusters', label = T)
# find cluster marker gene of fibroblast subclusters
library(dplyr) 
sce.markers <- FindAllMarkers(object = scRNA_sub,only.pos = TRUE,
                              min.pct = 0.25, 
                              thresh.use = 0.25)
# Fibroblast subclusters map to spatial transcriptomics
DefaultAssay(cortex) <- "SCT"
anchors <- FindTransferAnchors(reference = scRNA_sub, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = scRNA_sub$SCT_snn_res.0.5, 
                                  prediction.assay = TRUE, weight.reduction = cortex[["pca"]], dims = 1:50)
cortex[["predictions"]] <- predictions.assay
DefaultAssay(cortex) <- "predictions"
# plot Fibroblast subclusters mapping results
table(scRNA_sub$SCT_snn_res.0.5)
SpatialFeaturePlot(cortex, features = c("0", '1','2','3','4','5','6','7','8','9',
                                        '10','11'), pt.size.factor = 1.6,
         min.cutoff = 0,ncol = 2, crop = TRUE, alpha = c(0,1)) + scale_alpha()
ggsave('Fibroblast subclusters map.mix.pdf',width = 10, height = 20)
# We found that Fibroblast subcluster 2 maps to dental follicle area, than we extracted subcluster 2 to identify DFDP
cells_sub <- subset(scRNA_sub@meta.data, 
                    SCT_snn_res.0.5 %in%'2')
Fibro_c2 <- subset(scRNA_sub, 
                    cells=row.names(cells_sub))
Fibro_c2 <- SCTransform(Fibro_c2, ncells = 3000, verbose = FALSE)%>% RunPCA(
  verbose = F,npcs = 50)%>% FindNeighbors(dims = dim.use)%>% FindClusters(resolution = 0.6)%>% RunUMAP(dims = dim.use)
# Fibroblast subcluster 2 UMAP plot
DimPlot(Fibro_c2, group.by = 'SCT_snn_res.0.6', label = T)
table(Fibro_c2@meta.data$SCT_snn_res.0.6)
# read the gene list of mesenchymal between DC and SDL based on spatial transcriptomics loupe browser
geneList <- read.csv('../gene list.csv')
geneList <- as.character(geneList)
#Calculate scores for each cell based on the gene list
Fibro_c2 <- AddModuleScore(
  object = Fibro_c2,
  features = geneList,
  ctrl = 100, 
  name = 'score')
#Generate a feature plot based on the scores calculated above
FeaturePlot(Fibro_c2, features = "score")

# rename Fibro_c2
Fibro_c2 <- RenameIdents(Fibro_c2, '0'='c2-0',
                          '1'='c2-1',
                          '2'='c2-2')
Fibro_c2@meta.data$sub_cluster <- Idents(Fibro_c2)

## Mapping Fibro_c2 cluster to all Fibroblast subcluster
metadata <- scRNA_sub@meta.data
metadata$cell_type <- NULL
metadata$cell_type <- Fibro_c2@meta.data[match(rownames(scRNA_sub@meta.data), 
                                                rownames(Fibro_c2@meta.data)), 
                                          'sub_cluster']
metadata$cell_type <- as.character(metadata$cell_type)
metadata$SCT_snn_res.0.5 <- as.character(metadata$SCT_snn_res.0.5)
print(table(metadata$cell_type))
celltype_names <- NULL
for(i in 1:dim(metadata)[1]){
  print(i)
  sub_data <- metadata[i,]
  if(is.na(sub_data$cell_type)){
    print('Change value')
    sub_data$cell_type <- sub_data$SCT_snn_res.0.5
    celltype_names <- c(celltype_names, sub_data$cell_type)
  }else{
    celltype_names <- c(celltype_names, sub_data$cell_type)
  }
}
print(table(celltype_names))

metadata$cell_type <- celltype_names
metadata$cell_type <- factor(metadata$cell_type)
metadata$cell_type <- factor(metadata$cell_type, 
                             levels = c(names(table(Fibro_c2$sub_cluster)), 0,1,3:11))
print(table(metadata$cell_type))

scRNA_sub@meta.data <- metadata
DimPlot(scRNA_sub, reduction='umap', label=T, label.size=5, pt.size=1, group.by='SCT_snn_res.0.5') + ggtitle('original celltype')
DimPlot(scRNA_sub, reduction='umap', label=T, label.size=5, pt.size=1, group.by='cell_type') + ggtitle('new celltype')

# Fibroblast subclusters map to spatial transcriptomics
DefaultAssay(cortex) <- "SCT"
anchors <- FindTransferAnchors(reference = scRNA_sub, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = scRNA_sub$cell_type, 
                                  prediction.assay = TRUE, weight.reduction = cortex[["pca"]], dims = 1:50)
cortex[["predictions"]] <- predictions.assay
DefaultAssay(cortex) <- "predictions"
# plot Fibroblast subclusters mapping results
table(scRNA_sub$cell_type)
SpatialFeaturePlot(cortex, features = c("0", '1',"c2-0", "c2-1", "c2-2",'3','4','5','6','7','8','9',
                                        '10','11'), pt.size.factor = 1.6,
                   min.cutoff = 0,ncol = 2, crop = TRUE, alpha = c(0,1)) + scale_alpha()
ggsave('Fibroblast subclusters c2 map.mix.pdf',width = 10, height = 20)

# rename Fibroblast sub-cell types
Idents(scRNA_sub) <- scRNA_sub$cell_type
scRNA_sub <- RenameIdents(scRNA_sub, "c2-0"='DFDP',
                          "c2-1"='TDFCs',
                          "c2-2"='TDFCs',
                          '0'='BMSCs',
                          '1'='Dental mes',
                          '3'='Dental mes',
                          '4'='Dental mes',
                          '5'='Subcutaneous mes',
                          '6'='Subcutaneous mes',
                          '7'='Subcutaneous mes',
                          '8'='Subcutaneous mes',
                          '9'='Subcutaneous mes',
                          '10'='BMSCs',
                          '11'='Subcutaneous mes')
scRNA_sub@meta.data$defined <- Idents(scRNA_sub)
# Fibroblast sub-cell types UMAP plot
library(ggsci)
DimPlot(scRNA_sub, reduction = 'umap', group.by = 'defined',
        label = F, pt.size = 0.5) + scale_color_d3()
ggsave('Fibroblast sub-celltypes umap.pdf',width = 7,height = 6)

# Fibroblast sub-cell types marker gene heatmap
sce.markers <- FindAllMarkers(object = scRNA_sub, only.pos = TRUE,
                              min.pct = 0.25, 
                              thresh.use = 0.25)
color <- pal_d3( alpha = 1)(5)
ord = c('DFDP','TDFCs','BMSCs','Subcutaneous mes','Dental mes')
scRNA_sub$defined = factor(scRNA_sub$defined ,levels = ord)
top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
ll = split(top10$gene,top10$cluster)
ll = ll[ord]
rmg=names(table(unlist(ll))[table(unlist(ll))>1])
ll = lapply(ll, function(x) x[!x %in% rmg])
library(ggplot2)
DoHeatmap(scRNA_sub,
          features = unlist(ll),
          group.by = 'defined',
          assay = 'SCT',
          group.colors = color,label = F)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))
ggsave(filename = "Fibroblast sub-celltypes marker_pheatmap.pdf",width = 8,height = 5.5)

# FeaturePlot marker genes of DFDP, TDFCs and BMSCs
FeaturePlot(scRNA_sub,features = c('LAMC3'),min.cutoff = 0.5,
            cols = c("lightgrey","red"))
ggsave('LAMC3 FeaturePlot.pdf',width = 6,height = 6)
FeaturePlot(scRNA_sub,features = c('FBLN1'),min.cutoff = 1.5,
            cols = c("lightgrey","red"))
ggsave('FBLN1 FeaturePlot.pdf',width = 6,height = 6)
FeaturePlot(scRNA_sub,features = c('COL14A1'),min.cutoff = 0.5,
            cols = c("lightgrey","red"))
ggsave('COL14A1 FeaturePlot.pdf',width = 6,height = 6)
FeaturePlot(scRNA_sub,features = c('HHIP'),min.cutoff = 0.5,
            cols = c("lightgrey","red"))
ggsave('HHIP FeaturePlot.pdf',width = 6,height = 6)
FeaturePlot(scRNA_sub,features = c('TNC'),min.cutoff = 0.5,
            cols = c("lightgrey","red"))
ggsave('TNC FeaturePlot.pdf',width = 6,height = 6)
FeaturePlot(scRNA_sub,features = c('IBSP'),min.cutoff = 0.5,
            cols = c("lightgrey","red"))
ggsave('IBSP FeaturePlot.pdf',width = 6,height = 6)
FeaturePlot(scRNA_sub,features = c('RUNX2'),min.cutoff = 0.5,
            cols = c("lightgrey","red"))
ggsave('RUNX2 FeaturePlot.pdf',width = 6,height = 6)

# Fibroblast sub-cell types map to spatial transcriptomics
DefaultAssay(cortex) <- "SCT"
anchors <- FindTransferAnchors(reference = scRNA_sub, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = scRNA_sub$defined, 
                                  prediction.assay = TRUE, weight.reduction = cortex[["pca"]], dims = 1:50)
cortex[["predictions"]] <- predictions.assay
DefaultAssay(cortex) <- "predictions"

# plot Fibroblast sub-cell types mapping results
table(scRNA_sub$defined)
SpatialFeaturePlot(cortex, features = c("BMSCs", "Dental mes",
                                        "Subcutaneous mes",  "TDFCs",   
                                        'DFDP'), pt.size.factor = 1.6,
                   min.cutoff = 0.3,ncol = 2, crop = TRUE, alpha = c(0,1)) + 
  scale_alpha()
ggsave('Fibroblast sub-celltypes map.mix.pdf',width = 10, height = 20)

# DFDP vs TDFCs Differential analysis
deg = FindMarkers(scRNA_sub,ident.1 = 'DFDP',
                  ident.2 = 'TDFCs')
write.csv(deg, file = "deg.csv", 
          append = FALSE, quote = TRUE, sep = ",",
          eol = "\n", na = "NA", dec = ".", row.names = TRUE,
          col.names = TRUE, qmethod = c("escape", "double"),
          fileEncoding = "")
head(deg[order(deg$p_val),])
# Differential genes volcano plot
library(EnhancedVolcano) 
EnhancedVolcano(deg,
                pCutoff = 5e-02,
                FCcutoff = 0.58,
                lab = rownames(deg),
                x = 'avg_log2FC',
                y = 'p_val_adj')
ggsave(filename = 'DFDP vs TDFCs EnhancedVolcano.pdf',width = 8,height = 7)

# KEGG pathway analysis
up.gene = deg[deg$avg_log2FC > 0,]
up.gene$gene_symbol = rownames(up.gene)
KEGG_database <- 'ssc'
library(clusterProfiler)
GO_database <- 'org.Ss.eg.db'
gene <- bitr(up.gene$gene_symbol,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
KEGG<-enrichKEGG(gene$ENTREZID,
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
barplot(KEGG,title = 'KEGG Pathway')+
  scale_fill_gradient(low = 'Chocolate',high = 'Wheat')
ggsave(filename = 'KEGG pathway.pdf',height = 5,width = 7)

# StackedVlnPlot
library(ggplot2)
modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
  p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          #axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          plot.margin = plot.margin )
  return(p)
}
## main function
StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

color <- pal_d3( alpha = 1)(5)
StackedVlnPlot(scRNA_sub, c('LAMC3','COL6A3','TNC','NPNT','COL1A2','COL1A1',
                          'COL4A2','LAMB1','CTNNB1','IGF1','CCND1',
                          'NKD1','WNT5A','TEAD1','CCND1','CCN2'), pt.size=0, cols=color)
ggsave('StackedVlnPlot.pdf',width = 6,height = 20)


