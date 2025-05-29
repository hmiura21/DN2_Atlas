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

CD_Martin <-readRDS('/Users/honokamiura/Downloads/Research/DN2/CD_GSE134809_Martin_Bcells.rds')

META<-CD_Martin[[]]

DefaultAssay(CD_Martin) <- "RNA"

#QC %MT^
CD_Martin[["percent.mt"]] <- PercentageFeatureSet(CD_Martin, pattern = "^MT.")

#violin plot to look at percent.mt, nFeature_RNA
VlnPlot(CD_Martin, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
CD_Martin <- subset(CD_Martin, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
CD_Martin <- NormalizeData(CD_Martin)

#identify highly variably features (genes) total of 2000
CD_Martin <- FindVariableFeatures(CD_Martin, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
CD_Martin <- ScaleData(CD_Martin)

# Perform PCA and color by cell cycle phase
CD_Martin <- RunPCA(CD_Martin)

#Cluster the cells
CD_Martin <- FindNeighbors(CD_Martin, dims = 1:10)
CD_Martin <- FindClusters(CD_Martin, resolution = 0.9)

#Run UMAP
CD_Martin <- RunUMAP(CD_Martin, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(CD_Martin)

saveRDS(CD_Martin, "CD_Martin_QC.RDS")

#Start here: 
CD_Martin<-readRDS('CD_Martin_QC.RDS')


#Looking for Tbet+, 
DotPlot(object = CD_Martin, #best=2,0 (maybe 1,3)
        features = c('TBX21'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Look for cytotoxic markers 
DotPlot(object = CD_Martin, #group.by = 'Clusters',
        features = c('GZMB', 'GZMA','GZMM','PRF1','NKG7','KLRK1','KLRD1','GNLY'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#dn2 markers (3)
DotPlot(object = CD_Martin, group.by = 'seurat_clusters', 
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


CD_Martin <- AddModuleScore_UCell(CD_Martin,features = markers)

#Marker cytotoxic Score
CD_Martin$cytotoxic_UCell
plot_density(
  CD_Martin,
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
CD_Martin$DN2_pos_UCell
plot_density(
  CD_Martin,
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
CD_Martin$DN2_neg_UCell
plot_density(
  CD_Martin,
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
CD_Martin$tbet_UCell
plot_density(
  CD_Martin,
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



dittoBarPlot(
  object = CD_Martin,
  var = "clusters_annotated",    #colors stacked
  group.by = "seurat_clusters",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))


unique(CD_Martin$Atlas_Donor)

#ethnicty, study id, disease, donor, 

#ADD METADATA

#FOR WEBB B CELLS
#add study id
CD_Martin@meta.data$Atlas_Study_ID<-'Martin_CD'
CD_Martin@meta.data$Atlas_Disease<-'CD'

META<-CD_Martin[[]]

unique(CD_Martin$atlas)


saveRDS(CD_Martin, "CD_Martin_updated_meta.rds")
CD_Martin<-readRDS('CD_Martin_updated_meta.RDS')


#PROP ANALYSIS


#ADD METADATA for DN2 cluster and other clusters
unique(CD_Martin$seurat_clusters)
Idents(CD_Martin)<-'seurat_clusters'
current.cluster.ids<-levels(CD_Martin@active.ident)
current.cluster.ids<-c( '0','1','2','3','4','5','6','7','8','9','10','11','12','13','14')
new.cluster.ids<-c('other','other','other','other','other','other','other','other','other','other','other','DN2','other','other','other')

CD_Martin@meta.data$Atlas_DN2_Cluster<-CD_Martin@meta.data$seurat_clusters
CD_Martin@meta.data$Atlas_DN2_Cluster<-plyr::mapvalues(x=CD_Martin@meta.data$Atlas_DN2_Cluster, from=current.cluster.ids, to=new.cluster.ids)

unique(CD_Martin$Atlas_DN2_Cluster)


#ADD METADATA for cytotoxic cluster and other clusters
unique(CD_Martin$seurat_clusters)
Idents(CD_Martin)<-'seurat_clusters'
current.cluster.ids<-levels(CD_Martin@active.ident)
current.cluster.ids<-c('0','1','2','3','4','5','6','7','8','9','10','11','12','13','14')
new.cluster.ids<-c('other','other','other','other','other','other','other','other','other','other','cytotoxic','other','other','other','other')

CD_Martin@meta.data$Atlas_Cytotoxic_Cluster<-CD_Martin@meta.data$seurat_clusters
CD_Martin@meta.data$Atlas_Cytotoxic_Cluster<-plyr::mapvalues(x=CD_Martin@meta.data$Atlas_Cytotoxic_Cluster, from=current.cluster.ids, to=new.cluster.ids)

unique(CD_Martin$Atlas_Cytotoxic_Cluster)



#DN2 proportion 
dittoBarPlot(
  object = CD_Martin,
  var = "Atlas_Cytotoxic_Cluster",    #colors stacked
  group.by = "Atlas_Disease",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

#DN2 proportion by age
dittoBarPlot(
  object = CD_Martin,
  var = "Atlas_DN2_Cluster",    #colors stacked
  group.by = "Age",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

#DN2 proportion by ethnicity
dittoBarPlot(
  object = CD_Martin,
  var = "Atlas_DN2_Cluster",    #colors stacked
  group.by = "Atlas_Ethnicity_Broad",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))


#Statistical Significance
library("scProportionTest")

prop_test <- sc_utils(CD_Martin)

#prop test for DN2 cluster vs age 
prop_test_1_vs_2 <- permutation_test(
  prop_test, cluster_identity = "Atlas_DN2_Cluster",
  sample_1 = "adult", sample_2 = "child",
  sample_identity = "Age"
)
permutation_plot(prop_test_1_vs_2)


# Subset Cluster 11 (DN2) 
CD_Martin_Subset11<-subset(x = CD_Martin, idents ='11')

#QC %MT^
CD_Martin_Subset11[["percent.mt"]] <- PercentageFeatureSet(CD_Martin_Subset11, pattern = "^MT-")

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
CD_Martin_Subset11 <- subset(CD_Martin_Subset11, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
CD_Martin_Subset11 <- NormalizeData(CD_Martin_Subset11)

#identify highly variably features (genes) total of 2000
CD_Martin_Subset11 <- FindVariableFeatures(CD_Martin_Subset11, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
CD_Martin_Subset11 <- ScaleData(CD_Martin_Subset11)

# Perform PCA 
CD_Martin_Subset11 <- RunPCA(CD_Martin_Subset11)

#Cluster the cells
CD_Martin_Subset11 <- FindNeighbors(CD_Martin_Subset11, dims = 1:10)
CD_Martin_Subset11 <- FindClusters(CD_Martin_Subset11, resolution = 0.5)

#Run UMAP
CD_Martin_Subset11 <- RunUMAP(CD_Martin_Subset11, dims = 1:10)


#Dim Plot showing Clusters
DimPlot(CD_Martin_Subset11, reduction = "umap")

saveRDS(CD_Martin_Subset11, "CD_Martin_DN2_QC.RDS")

CD_Martin_Subset11<-readRDS('CD_Martin_DN2_QC.RDS')





# Subset Cluster 10 (cytotoxic) 
CD_Martin_Subset10<-subset(x = CD_Martin, idents ='4')

#QC %MT^
CD_Martin_Subset10[["percent.mt"]] <- PercentageFeatureSet(CD_Martin_Subset10, pattern = "^MT-")

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
CD_Martin_Subset10 <- subset(CD_Martin_Subset10, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
CD_Martin_Subset10 <- NormalizeData(CD_Martin_Subset10)

#identify highly variably features (genes) total of 2000
CD_Martin_Subset10 <- FindVariableFeatures(CD_Martin_Subset10, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
CD_Martin_Subset10 <- ScaleData(CD_Martin_Subset10)

# Perform PCA 
CD_Martin_Subset10 <- RunPCA(CD_Martin_Subset10)

#Cluster the cells
CD_Martin_Subset10 <- FindNeighbors(CD_Martin_Subset10, dims = 1:10)
CD_Martin_Subset10 <- FindClusters(CD_Martin_Subset10, resolution = 0.5)

#Run UMAP
CD_Martin_Subset10 <- RunUMAP(CD_Martin_Subset10, dims = 1:10)


#Dim Plot showing Clusters
DimPlot(CD_Martin_Subset10, reduction = "umap")

saveRDS(CD_Martin_Subset10, "CD_Martin_cytotoxic_QC.RDS")

CD_Martin_Subset10<-readRDS('CD_Martin_cytotoxic_QC.RDS')

