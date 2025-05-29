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

Malaria_Holla <-readRDS('/Users/honokamiura/Downloads/Research/DN2/Malaria_GSE149729_B_cells.rds')

META<-Malaria_Holla[[]]

DefaultAssay(Malaria_Holla) <- "RNA"

#QC %MT^
Malaria_Holla[["percent.mt"]] <- PercentageFeatureSet(Malaria_Holla, pattern = "^MT.")

#violin plot to look at percent.mt, nFeature_RNA
VlnPlot(Malaria_Holla, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
Malaria_Holla <- subset(Malaria_Holla, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
Malaria_Holla <- NormalizeData(Malaria_Holla)

#identify highly variably features (genes) total of 2000
Malaria_Holla <- FindVariableFeatures(Malaria_Holla, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
Malaria_Holla <- ScaleData(Malaria_Holla)

# Perform PCA and color by cell cycle phase
Malaria_Holla <- RunPCA(Malaria_Holla)

#Cluster the cells
Malaria_Holla <- FindNeighbors(Malaria_Holla, dims = 1:10)
Malaria_Holla <- FindClusters(Malaria_Holla, resolution = 0.9)

#Run UMAP
Malaria_Holla <- RunUMAP(Malaria_Holla, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(Malaria_Holla)

saveRDS(Malaria_Holla, "Malaria_Holla_QC.RDS")

#Start here: 
Malaria_Holla<-readRDS('Malaria_Holla_QC.RDS')


#Looking for Tbet+, 
DotPlot(object = Malaria_Holla, #best=3
        features = c('TBX21'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Look for cytotoxic markers best=0
DotPlot(object = Malaria_Holla, #group.by = 'Clusters',
        features = c('GZMB', 'GZMA','GZMM','PRF1','NKG7','KLRK1','KLRD1','GNLY'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#dn2 markers (1)
DotPlot(object = Malaria_Holla, group.by = 'seurat_clusters', 
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


Malaria_Holla <- AddModuleScore_UCell(Malaria_Holla,features = markers)

#Marker cytotoxic Score
Malaria_Holla$cytotoxic_UCell
plot_density(
  Malaria_Holla,
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
Malaria_Holla$DN2_pos_UCell
plot_density(
  Malaria_Holla,
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
Malaria_Holla$DN2_neg_UCell
plot_density(
  Malaria_Holla,
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
Malaria_Holla$tbet_UCell
plot_density(
  Malaria_Holla,
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
sort(unique(Malaria_Holla$seurat_clusters))

unique(Malaria_Holla$seurat_clusters)
Idents(Malaria_Holla)<-'seurat_clusters'
current.cluster.ids<-levels(Malaria_Holla@active.ident)
current.cluster.ids<-c( '0','1','2','3','4','5','6','7')
new.cluster.ids<-c('other','other','other','DN2','other','other','other','other' )

Malaria_Holla@meta.data$Atlas_DN2_Cluster<-Malaria_Holla@meta.data$seurat_clusters
Malaria_Holla@meta.data$Atlas_DN2_Cluster<-plyr::mapvalues(x=Malaria_Holla@meta.data$Atlas_DN2_Cluster, from=current.cluster.ids, to=new.cluster.ids)

unique(Malaria_Holla$Atlas_DN2_Cluster)


#ADD METADATA for cytotoxic cluster and other clusters
unique(Malaria_Holla$seurat_clusters)
Idents(Malaria_Holla)<-'seurat_clusters'
current.cluster.ids<-levels(Malaria_Holla@active.ident)
current.cluster.ids<-c( '0','1','2','3','4','5','6','7')
new.cluster.ids<-c('other','other','other','other','other','other','other','cytotoxic' )

Malaria_Holla@meta.data$Atlas_Cytotoxic_Cluster<-Malaria_Holla@meta.data$seurat_clusters
Malaria_Holla@meta.data$Atlas_Cytotoxic_Cluster<-plyr::mapvalues(x=Malaria_Holla@meta.data$Atlas_Cytotoxic_Cluster, from=current.cluster.ids, to=new.cluster.ids)

unique(Malaria_Holla$Atlas_Cytotoxic_Cluster)

unique(Malaria_Holla$Atlas_Disease)


#DN2 proportion by disease
dittoBarPlot(
  object = Malaria_Holla,
  var = "Atlas_DN2_Cluster",    #colors stacked
  group.by = "Atlas_Disease",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

dittoBarPlot(
  object = Malaria_Holla,
  var = "Atlas_Cytotoxic_Cluster",    #colors stacked
  group.by = "Atlas_Disease",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))




# Subset Cluster 3 (DN2) 
Malaria_Holla_Subset3<-subset(x = Malaria_Holla, idents ='3')

#QC %MT^
Malaria_Holla_Subset3[["percent.mt"]] <- PercentageFeatureSet(Malaria_Holla_Subset3, pattern = "^MT-")

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
Malaria_Holla_Subset3 <- subset(Malaria_Holla_Subset3, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
Malaria_Holla_Subset3 <- NormalizeData(Malaria_Holla_Subset3)

#identify highly variably features (genes) total of 2000
Malaria_Holla_Subset3 <- FindVariableFeatures(Malaria_Holla_Subset3, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
Malaria_Holla_Subset3 <- ScaleData(Malaria_Holla_Subset3)

# Perform PCA 
Malaria_Holla_Subset3 <- RunPCA(Malaria_Holla_Subset3)

#Cluster the cells
Malaria_Holla_Subset3 <- FindNeighbors(Malaria_Holla_Subset3, dims = 1:10)
Malaria_Holla_Subset3 <- FindClusters(Malaria_Holla_Subset3, resolution = 0.5)

#Run UMAP
Malaria_Holla_Subset3 <- RunUMAP(Malaria_Holla_Subset3, dims = 1:10)


#Dim Plot showing Clusters
DimPlot(Malaria_Holla_Subset3, reduction = "umap")

saveRDS(Malaria_Holla_Subset3, "Malaria_Holla_DN2_QC.RDS")

Malaria_Holla_Subset3<-readRDS('Malaria_Holla_DN2_QC.RDS')





# Subset Cluster 7 (cytotoxic) 
Malaria_Holla_Subset7<-subset(x = Malaria_Holla, idents ='7')

#QC %MT^
Malaria_Holla_Subset7[["percent.mt"]] <- PercentageFeatureSet(Malaria_Holla_Subset7, pattern = "^MT-")

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
Malaria_Holla_Subset7 <- subset(Malaria_Holla_Subset7, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
Malaria_Holla_Subset7 <- NormalizeData(Malaria_Holla_Subset7)

#identify highly variably features (genes) total of 2000
Malaria_Holla_Subset7 <- FindVariableFeatures(Malaria_Holla_Subset7, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
Malaria_Holla_Subset7 <- ScaleData(Malaria_Holla_Subset7)

# Perform PCA 
Malaria_Holla_Subset7 <- RunPCA(Malaria_Holla_Subset7)

#Cluster the cells
Malaria_Holla_Subset7 <- FindNeighbors(Malaria_Holla_Subset7, dims = 1:10)
Malaria_Holla_Subset7 <- FindClusters(Malaria_Holla_Subset7, resolution = 0.5)

#Run UMAP
Malaria_Holla_Subset7 <- RunUMAP(Malaria_Holla_Subset7, dims = 1:10)


#Dim Plot showing Clusters
DimPlot(Malaria_Holla_Subset7, reduction = "umap")

saveRDS(Malaria_Holla_Subset7, "Malaria_Holla_cytotoxic_QC.RDS")

Malaria_Holla_Subset7<-readRDS('Malaria_Holla_cytotoxic_QC.RDS')
