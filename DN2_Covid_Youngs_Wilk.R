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

Covid_Youngs_Wilk <-readRDS('/Users/honokamiura/Downloads/Research/DN2/SEVERE_COVID_GSE178404_GSE150728_B_cells.rds')

META<-Covid_Youngs_Wilk[[]]

DefaultAssay(Covid_Youngs_Wilk) <- "RNA"

#QC %MT^
Covid_Youngs_Wilk[["percent.mt"]] <- PercentageFeatureSet(Covid_Youngs_Wilk, pattern = "^MT.")

#violin plot to look at percent.mt, nFeature_RNA
VlnPlot(Covid_Youngs_Wilk, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
Covid_Youngs_Wilk <- subset(Covid_Youngs_Wilk, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
Covid_Youngs_Wilk <- NormalizeData(Covid_Youngs_Wilk)

#identify highly variably features (genes) total of 2000
Covid_Youngs_Wilk <- FindVariableFeatures(Covid_Youngs_Wilk, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
Covid_Youngs_Wilk <- ScaleData(Covid_Youngs_Wilk)

# Perform PCA and color by cell cycle phase
Covid_Youngs_Wilk <- RunPCA(Covid_Youngs_Wilk)

#Cluster the cells
Covid_Youngs_Wilk <- FindNeighbors(Covid_Youngs_Wilk, dims = 1:10)
Covid_Youngs_Wilk <- FindClusters(Covid_Youngs_Wilk, resolution = 0.9)

#Run UMAP
Covid_Youngs_Wilk <- RunUMAP(Covid_Youngs_Wilk, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(Covid_Youngs_Wilk)

saveRDS(Covid_Youngs_Wilk, "Covid_Youngs_Wilk_QC.RDS")

#Start here: 
Covid_Youngs_Wilk<-readRDS('Covid_Youngs_Wilk_QC.RDS')


#Looking for Tbet+, 
DotPlot(object = Covid_Youngs_Wilk, #best=3
        features = c('TBX21'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Look for cytotoxic markers best=0
DotPlot(object = Covid_Youngs_Wilk, #group.by = 'Clusters',
        features = c('GZMB','GZMA','GZMM','PRF1','NKG7','KLRK1','KLRD1','GNLY'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#dn2 markers (1)
DotPlot(object = Covid_Youngs_Wilk, group.by = 'seurat_clusters', 
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


Covid_Youngs_Wilk <- AddModuleScore_UCell(Covid_Youngs_Wilk,features = markers)

#Marker cytotoxic Score
Covid_Youngs_Wilk$cytotoxic_UCell
plot_density(
  Covid_Youngs_Wilk,
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
Covid_Youngs_Wilk$DN2_pos_UCell
plot_density(
  Covid_Youngs_Wilk,
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
Covid_Youngs_Wilk$DN2_neg_UCell
plot_density(
  Covid_Youngs_Wilk,
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
Covid_Youngs_Wilk$tbet_UCell
plot_density(
  Covid_Youngs_Wilk,
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
sort(unique(Covid_Youngs_Wilk$seurat_clusters))

unique(Covid_Youngs_Wilk$seurat_clusters)
Idents(Covid_Youngs_Wilk)<-'seurat_clusters'
current.cluster.ids<-levels(Covid_Youngs_Wilk@active.ident)
current.cluster.ids<-c( '0','1','2','3','4','5','6','7','8')
new.cluster.ids<-c('other','other','other','other','other','other','DN2','other','other')

Covid_Youngs_Wilk@meta.data$Atlas_DN2_Cluster<-Covid_Youngs_Wilk@meta.data$seurat_clusters
Covid_Youngs_Wilk@meta.data$Atlas_DN2_Cluster<-plyr::mapvalues(x=Covid_Youngs_Wilk@meta.data$Atlas_DN2_Cluster, from=current.cluster.ids, to=new.cluster.ids)



#ADD METADATA for cytotoxic cluster and other clusters
unique(Covid_Youngs_Wilk$seurat_clusters)
Idents(Covid_Youngs_Wilk)<-'seurat_clusters'
current.cluster.ids<-levels(Covid_Youngs_Wilk@active.ident)
current.cluster.ids<-c( '0','1','2','3','4','5','6','7','8')
new.cluster.ids<-c('other','cytotoxic','other','other','other','other','other','other','other')

Covid_Youngs_Wilk@meta.data$Atlas_Cytotoxic_Cluster<-Covid_Youngs_Wilk@meta.data$seurat_clusters
Covid_Youngs_Wilk@meta.data$Atlas_Cytotoxic_Cluster<-plyr::mapvalues(x=Covid_Youngs_Wilk@meta.data$Atlas_Cytotoxic_Cluster, from=current.cluster.ids, to=new.cluster.ids)

unique(Covid_Youngs_Wilk$Atlas_Cytotoxic_Cluster)


#DN2 proportion 
dittoBarPlot(
  object = Covid_Youngs_Wilk,
  var = "Atlas_Cytotoxic_Cluster",    #colors stacked
  group.by = "Atlas_Disease",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

#DN2 proportion 
dittoBarPlot(
  object = Covid_Youngs_Wilk,
  var = "Atlas_DN2_Cluster",    #colors stacked
  group.by = "Atlas_Disease",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))



# Subset Cluster 6 (DN2) 
Covid_Youngs_Wilk_Subset6<-subset(x = Covid_Youngs_Wilk, idents ='6')

#QC %MT^
Covid_Youngs_Wilk_Subset6[["percent.mt"]] <- PercentageFeatureSet(Covid_Youngs_Wilk_Subset6, pattern = "^MT-")

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
Covid_Youngs_Wilk_Subset6 <- subset(Covid_Youngs_Wilk_Subset6, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
Covid_Youngs_Wilk_Subset6 <- NormalizeData(Covid_Youngs_Wilk_Subset6)

#identify highly variably features (genes) total of 2000
Covid_Youngs_Wilk_Subset6 <- FindVariableFeatures(Covid_Youngs_Wilk_Subset6, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
Covid_Youngs_Wilk_Subset6 <- ScaleData(Covid_Youngs_Wilk_Subset6)

# Perform PCA 
Covid_Youngs_Wilk_Subset6 <- RunPCA(Covid_Youngs_Wilk_Subset6)

#Cluster the cells
Covid_Youngs_Wilk_Subset6 <- FindNeighbors(Covid_Youngs_Wilk_Subset6, dims = 1:10)
Covid_Youngs_Wilk_Subset6 <- FindClusters(Covid_Youngs_Wilk_Subset6, resolution = 0.5)

#Run UMAP
Covid_Youngs_Wilk_Subset6 <- RunUMAP(Covid_Youngs_Wilk_Subset6, dims = 1:10)


#Dim Plot showing Clusters
DimPlot(Covid_Youngs_Wilk_Subset6, reduction = "umap")

saveRDS(Covid_Youngs_Wilk_Subset6, "Covid_Youngs_Wilk_DN2_QC.RDS")

Covid_Youngs_Wilk_Subset6<-readRDS('Covid_Youngs_Wilk_DN2_QC.RDS')





# Subset Cluster 1 (cytotoxic) 
Covid_Youngs_Wilk_Subset1<-subset(x = Covid_Youngs_Wilk, idents ='1')

#QC %MT^
Covid_Youngs_Wilk_Subset1[["percent.mt"]] <- PercentageFeatureSet(Covid_Youngs_Wilk_Subset1, pattern = "^MT-")

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
Covid_Youngs_Wilk_Subset1 <- subset(Covid_Youngs_Wilk_Subset1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
Covid_Youngs_Wilk_Subset1 <- NormalizeData(Covid_Youngs_Wilk_Subset1)

#identify highly variably features (genes) total of 2000
Covid_Youngs_Wilk_Subset1 <- FindVariableFeatures(Covid_Youngs_Wilk_Subset1, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
Covid_Youngs_Wilk_Subset1 <- ScaleData(Covid_Youngs_Wilk_Subset1)

# Perform PCA 
Covid_Youngs_Wilk_Subset1 <- RunPCA(Covid_Youngs_Wilk_Subset1)

#Cluster the cells
Covid_Youngs_Wilk_Subset1 <- FindNeighbors(Covid_Youngs_Wilk_Subset1, dims = 1:10)
Covid_Youngs_Wilk_Subset1 <- FindClusters(Covid_Youngs_Wilk_Subset1, resolution = 0.5)

#Run UMAP
Covid_Youngs_Wilk_Subset1 <- RunUMAP(Covid_Youngs_Wilk_Subset1, dims = 1:10)


#Dim Plot showing Clusters
DimPlot(Covid_Youngs_Wilk_Subset1, reduction = "umap")

saveRDS(Covid_Youngs_Wilk_Subset1, "Covid_Youngs_Wilk_cytotoxic_QC.RDS")

Covid_Youngs_Wilk_Subset1<-readRDS('Covid_Youngs_Wilk_cytotoxic_QC.RDS')
