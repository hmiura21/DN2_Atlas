#load packages
library(Seurat)
library(dittoSeq)
library(ggpubr)
library(cowplot)
library(ensembldb)
library(AnnotationHub)
library(RCurl)
library(dplyr)
library(patchwork)

setwd("/Users/honokamiura/Downloads/Research/DN2")

HIV_Holla <-readRDS('/Users/honokamiura/Downloads/Research/DN2/HIV_GSE157966_B_cells.rds')

META<-HIV_Holla[[]]

DefaultAssay(HIV_Holla) <- "RNA"

#QC %MT^
HIV_Holla[["percent.mt"]] <- PercentageFeatureSet(HIV_Holla, pattern = "^MT.")

#violin plot to look at percent.mt, nFeature_RNA
VlnPlot(HIV_Holla, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
HIV_Holla <- subset(HIV_Holla, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
HIV_Holla <- NormalizeData(HIV_Holla)

#identify highly variably features (genes) total of 2000
HIV_Holla <- FindVariableFeatures(HIV_Holla, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
HIV_Holla <- ScaleData(HIV_Holla)

# Perform PCA and color by cell cycle phase
HIV_Holla <- RunPCA(HIV_Holla)

#Cluster the cells
HIV_Holla <- FindNeighbors(HIV_Holla, dims = 1:10)
HIV_Holla <- FindClusters(HIV_Holla, resolution = 0.9)

#Run UMAP
HIV_Holla <- RunUMAP(HIV_Holla, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(HIV_Holla)

saveRDS(HIV_Holla, "HIV_Holla_QC.RDS")

#Start here: 
HIV_Holla<-readRDS('HIV_Holla_QC.RDS')


#Looking for Tbet+, 
DotPlot(object = HIV_Holla, #best=3
        features = c('TBX21'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Look for cytotoxic markers best=0
DotPlot(object = HIV_Holla, #group.by = 'Clusters',
        features = c('GZMB', 'GZMA','GZMM','PRF1','NKG7','KLRK1','KLRD1','GNLY'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#dn2 markers (1)
DotPlot(object = HIV_Holla, group.by = 'seurat_clusters', 
        features = c("SOX5",'HLA-DRB5','TBX21',"ITGAX","BATF",'JAZF1','NFATC2',
                     'ZEB2','NR4A2','MS4A1','TOX','FGR','IL10RA','FCRL5',
                     'SREBF2','TFEB','TFEC','ZBTB32','MEF2A','POU2F2','SPI1',
                     #NEGATIVE
                     'LTB','VPREB3','SELL','JCHAIN','CD27','IGHD'),
        cols = c('RdBu') ) + RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 7,face = "bold"))



library(UCell)
library(Nebulosa)

# Markers to Show for DN2: Density Plot
markers <- list()
markers$cytotoxic<- c('GZMB', 'GZMA','GZMM','PRF1','NKG7','KLRK1','KLRD1','GNLY')
markers$DN2_pos<- c("SOX5",'HLA-DRB5','TBX21',"ITGAX","BATF",'JAZF1','NFATC2',
                    'ZEB2','NR4A2','MS4A1','TOX','FGR','IL10RA','FCRL5',
                    'SREBF2','TFEB','TFEC','ZBTB32','MEF2A','POU2F2','SPI1'
)
markers$DN2_neg<- c('LTB','VPREB3','SELL','JCHAIN','CD27','IGHD')
markers$tbet<- c('TBX21')


HIV_Holla <- AddModuleScore_UCell(HIV_Holla,features = markers)

#Marker cytotoxic Score
HIV_Holla$cytotoxic_UCell
plot_density(
  HIV_Holla,
  c('cytotoxic_UCell'),
  slot = NULL,
  joint = F,
  reduction = 'umap',
  dims = c(1, 2),
  method = c("ks", "wkde"),
  adjust = 1,
  size = 1,
  shape = 16,
  combine = TRUE,
  pal = "plasma"
)& NoAxes() & NoLegend()
#Marker DN2 Pos Score
HIV_Holla$DN2_pos_UCell
plot_density(
  HIV_Holla,
  c('DN2_pos_UCell'),
  slot = NULL,
  joint = F,
  reduction = 'umap',
  dims = c(1, 2),
  method = c("ks", "wkde"),
  adjust = 1,
  size = 1,
  shape = 16,
  combine = TRUE,
  pal = "plasma"
)& NoAxes() & NoLegend()
#Marker DN2 Neg Score
HIV_Holla$DN2_neg_UCell
plot_density(
  HIV_Holla,
  c('DN2_neg_UCell'),
  slot = NULL,
  joint = F,
  reduction = 'umap',
  dims = c(1, 2),
  method = c("ks", "wkde"),
  adjust = 1,
  size = 1,
  shape = 16,
  combine = TRUE,
  pal = "plasma"
)& NoAxes() & NoLegend()
#Marker tbet
HIV_Holla$tbet_UCell
plot_density(
  HIV_Holla,
  c('tbet_UCell'),
  slot = NULL,
  joint = F,
  reduction = 'umap',
  dims = c(1, 2),
  method = c("ks", "wkde"),
  adjust = 1,
  size = 1,
  shape = 16,
  combine = TRUE,
  pal = "plasma"
)& NoAxes() & NoLegend()


#ADD METADATA for DN2 cluster and other clusters
sort(unique(HIV_Holla$seurat_clusters))

unique(HIV_Holla$seurat_clusters)
Idents(HIV_Holla)<-'seurat_clusters'
current.cluster.ids<-levels(HIV_Holla@active.ident)
current.cluster.ids<-c( '0','1','2','3','4','5','6','7','8','9')
new.cluster.ids<-c('DN2','other','other','other','other','other','other','DN2','other','other')

HIV_Holla@meta.data$Atlas_DN2_Cluster<-HIV_Holla@meta.data$seurat_clusters
HIV_Holla@meta.data$Atlas_DN2_Cluster<-plyr::mapvalues(x=HIV_Holla@meta.data$Atlas_DN2_Cluster, from=current.cluster.ids, to=new.cluster.ids)

unique(HIV_Holla$Atlas_DN2_Cluster)


#ADD METADATA for cytotoxic cluster and other clusters
unique(HIV_Holla$seurat_clusters)
Idents(HIV_Holla)<-'seurat_clusters'
current.cluster.ids<-levels(HIV_Holla@active.ident)
current.cluster.ids<-c( '0','1','2','3','4','5','6','7','8','9')
new.cluster.ids<-c('other','other','other','other','other','cytotoxic','other','other','other','other')


HIV_Holla@meta.data$Atlas_Cytotoxic_Cluster<-HIV_Holla@meta.data$seurat_clusters
HIV_Holla@meta.data$Atlas_Cytotoxic_Cluster<-plyr::mapvalues(x=HIV_Holla@meta.data$Atlas_Cytotoxic_Cluster, from=current.cluster.ids, to=new.cluster.ids)

unique(HIV_Holla$Atlas_Cytotoxic_Cluster)



#DN2 proportion by disease
dittoBarPlot(
  object = HIV_Holla,
  var = "Atlas_DN2_Cluster",    #colors stacked
  group.by = "Atlas_Disease",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

dittoBarPlot(
  object = HIV_Holla,
  var = "Atlas_Cytotoxic_Cluster",    #colors stacked
  group.by = "Atlas_Disease",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))





# Subset Cluster 0,7 (DN2) 
clusters_to_combine <- c(0,7)
HIV_Holla_Subset07<-subset(x = HIV_Holla, idents =clusters_to_combine)


#QC %MT^
HIV_Holla_Subset07[["percent.mt"]] <- PercentageFeatureSet(HIV_Holla_Subset07, pattern = "^MT-")

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
HIV_Holla_Subset07 <- subset(HIV_Holla_Subset07, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
HIV_Holla_Subset07 <- NormalizeData(HIV_Holla_Subset07)

#identify highly variably features (genes) total of 2000
HIV_Holla_Subset07 <- FindVariableFeatures(HIV_Holla_Subset07, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
HIV_Holla_Subset07 <- ScaleData(HIV_Holla_Subset07)

# Perform PCA 
HIV_Holla_Subset07 <- RunPCA(HIV_Holla_Subset07)

#Cluster the cells
HIV_Holla_Subset07 <- FindNeighbors(HIV_Holla_Subset07, dims = 1:10)
HIV_Holla_Subset07 <- FindClusters(HIV_Holla_Subset07, resolution = 0.5)

#Run UMAP
HIV_Holla_Subset07 <- RunUMAP(HIV_Holla_Subset07, dims = 1:10)


#Dim Plot showing Clusters
DimPlot(HIV_Holla_Subset07, reduction = "umap")

saveRDS(HIV_Holla_Subset07, "HIV_Holla_Subset07_DN2_QC.RDS")

HIV_Holla_Subset07<-readRDS('HIV_Holla_Subset07_DN2_QC.RDS')





# Subset Cluster 5 (cytotoxic) 
HIV_Holla_Subset5<-subset(x = HIV_Holla, idents ='5')

#QC %MT^
HIV_Holla_Subset5[["percent.mt"]] <- PercentageFeatureSet(HIV_Holla_Subset5, pattern = "^MT-")

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
HIV_Holla_Subset5 <- subset(HIV_Holla_Subset5, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
HIV_Holla_Subset5 <- NormalizeData(HIV_Holla_Subset5)

#identify highly variably features (genes) total of 2000
HIV_Holla_Subset5 <- FindVariableFeatures(HIV_Holla_Subset5, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
HIV_Holla_Subset5 <- ScaleData(HIV_Holla_Subset5)

# Perform PCA 
HIV_Holla_Subset5 <- RunPCA(HIV_Holla_Subset5)

#Cluster the cells
HIV_Holla_Subset5 <- FindNeighbors(HIV_Holla_Subset5, dims = 1:10)
HIV_Holla_Subset5 <- FindClusters(HIV_Holla_Subset5, resolution = 0.5)

#Run UMAP
HIV_Holla_Subset5 <- RunUMAP(HIV_Holla_Subset5, dims = 1:10)


#Dim Plot showing Clusters
DimPlot(HIV_Holla_Subset5, reduction = "umap")

saveRDS(HIV_Holla_Subset5, "HIV_Holla_cytotoxic_QC.RDS")

HIV_Holla_Subset5<-readRDS('HIV_Holla_cytotoxic_QC.RDS')
