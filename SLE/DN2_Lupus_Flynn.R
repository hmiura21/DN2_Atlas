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

Lupus_Flynn <-readRDS('/Users/honokamiura/Downloads/Research/DN2/GSE135779_SLE_Flynn_Bcells.rds')

META<-Lupus_Flynn[[]]

DefaultAssay(merged_vaccines_python) <- "RNA"

#QC %MT^
Lupus_Flynn[["percent.mt"]] <- PercentageFeatureSet(Lupus_Flynn, pattern = "^MT.")

#violin plot to look at percent.mt, nFeature_RNA
VlnPlot(Lupus_Flynn, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
Lupus_Flynn <- subset(Lupus_Flynn, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
Lupus_Flynn <- NormalizeData(Lupus_Flynn)

#identify highly variably features (genes) total of 2000
Lupus_Flynn <- FindVariableFeatures(Lupus_Flynn, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
Lupus_Flynn <- ScaleData(Lupus_Flynn)

# Perform PCA and color by cell cycle phase
Lupus_Flynn <- RunPCA(Lupus_Flynn)

#Cluster the cells
Lupus_Flynn <- FindNeighbors(Lupus_Flynn, dims = 1:10)
Lupus_Flynn <- FindClusters(Lupus_Flynn, resolution = 0.5)

#Run UMAP
Lupus_Flynn <- RunUMAP(Lupus_Flynn, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(Lupus_Flynn)

saveRDS(Lupus_Flynn, "Lupus_Flynn_QC.RDS")

#Start here: 
Lupus_Flynn<-readRDS('Lupus_Flynn_QC.RDS')



#Looking for Tbet+, (8,5?)
DotPlot(object = Lupus_Flynn, #best=2,0 (maybe 1,3)
        features = c('TBX21'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Look for cytotoxic markers (4,5,8)
DotPlot(object = Lupus_Flynn, #group.by = 'Clusters',
        features = c('GZMB', 'GZMA','GZMM','PRF1','NKG7','KLRK1','KLRD1','GNLY'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#dn2 markers (10,8)
DotPlot(object = Lupus_Flynn, group.by = 'seurat_clusters', 
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
markers$DN2_neg<- c('VPREB3','SELL','CD27','IGHD')
markers$tbet<- c('TBX21')


Lupus_Flynn <- AddModuleScore_UCell(Lupus_Flynn,features = markers)

#Marker cytotoxic Score
Lupus_Flynn$cytotoxic_UCell
plot_density(
  Lupus_Flynn,
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
Lupus_Flynn$DN2_pos_UCell
plot_density(
  Lupus_Flynn,
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
Lupus_Flynn$DN2_neg_UCell
plot_density(
  Lupus_Flynn,
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
Lupus_Flynn$tbet_UCell
plot_density(
  Lupus_Flynn,
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
  object = Lupus_Flynn,
  var = "clusters_annotated",    #colors stacked
  group.by = "seurat_clusters",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))




#ADD METADATA

#FOR WEBB B CELLS
#add study id
Lupus_Flynn@meta.data$Atlas_Study_ID<-'Flynn_Lupus'
Lupus_Flynn@meta.data$Atlas_Disease<-'SLE Active'

META<-Lupus_Flynn[[]]

unique(Lupus_Flynn$Atlas_Study_ID)


saveRDS(Lupus_Flynn, "Lupus_Flynn_updated_meta.rds")
Lupus_Flynn<-readRDS('Lupus_Flynn_updated_meta.RDS')



#PROP ANALYSIS


#ADD METADATA for DN2 cluster and other clusters
unique(Lupus_Flynn$seurat_clusters)
Idents(Lupus_Flynn)<-'seurat_clusters'
current.cluster.ids<-levels(Lupus_Flynn@active.ident)
current.cluster.ids<-c( '0','1','2','3','4','5','6','7','8','9','10','11')
new.cluster.ids<-c('other','other','other','other','other','other','other','other','DN2','other','other','other' )

Lupus_Flynn@meta.data$Atlas_DN2_Cluster<-Lupus_Flynn@meta.data$seurat_clusters
Lupus_Flynn@meta.data$Atlas_DN2_Cluster<-plyr::mapvalues(x=Lupus_Flynn@meta.data$Atlas_DN2_Cluster, from=current.cluster.ids, to=new.cluster.ids)

unique(Lupus_Flynn$Atlas_DN2_Cluster)


#ADD METADATA for cytotoxic cluster and other clusters
unique(Lupus_Flynn$seurat_clusters)
Idents(Lupus_Flynn)<-'seurat_clusters'
current.cluster.ids<-levels(Lupus_Flynn@active.ident)
current.cluster.ids<-c('0','1','2','3','4','5','6','7','8','9','10','11')
new.cluster.ids<-c('other','other','other','other','other','cytotoxic','other','other','other','other','other','other')

Lupus_Flynn@meta.data$Atlas_Cytotoxic_Cluster<-Lupus_Flynn@meta.data$seurat_clusters
Lupus_Flynn@meta.data$Atlas_Cytotoxic_Cluster<-plyr::mapvalues(x=Lupus_Flynn@meta.data$Atlas_Cytotoxic_Cluster, from=current.cluster.ids, to=new.cluster.ids)

unique(Lupus_Flynn$Atlas_Cytotoxic_Cluster)



#DN2 proportion 
dittoBarPlot(
  object = Lupus_Flynn,
  var = "Atlas_Cytotoxic_Cluster",    #colors stacked
  group.by = "Atlas_Disease",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

#DN2 proportion by age
dittoBarPlot(
  object = Lupus_Flynn,
  var = "Atlas_Cytotoxic_Cluster",    #colors stacked
  group.by = "Age",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

#DN2 proportion by ethnicity
dittoBarPlot(
  object = Lupus_Flynn,
  var = "Atlas_Cytotoxic_Cluster",    #colors stacked
  group.by = "Atlas_Ethnicity_Broad",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))


#Statistical Significance
library("scProportionTest")

prop_test <- sc_utils(Lupus_Flynn)

#prop test for DN2 cluster vs age 
prop_test_1_vs_2 <- permutation_test(
  prop_test, cluster_identity = "Atlas_DN2_Cluster",
  sample_1 = "adult", sample_2 = "child",
  sample_identity = "Age"
)
permutation_plot(prop_test_1_vs_2)



unique(Lupus_Flynn$Atlas_Donor)

#ADD METADATA

#add study id
Lupus_Flynn@meta.data$Atlas_Study_ID<-'Flynn_Lupus'
Lupus_Flynn@meta.data$Atlas_Disease<-'SLE Active'
META<-Lupus_Flynn[[]]

saveRDS(Lupus_Flynn, "Lupus_Flynn_updated_meta.rds")
Lupus_Flynn<-readRDS('Lupus_Flynn_updated_meta.RDS')






# Subset Cluster 8 (DN2/cytotoxic) 
Lupus_Flynn_Subset8<-subset(x = Lupus_Flynn, idents ='8')

#QC %MT^
Lupus_Flynn_Subset8[["percent.mt"]] <- PercentageFeatureSet(Lupus_Flynn_Subset8, pattern = "^MT-")

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
Lupus_Flynn_Subset8 <- subset(Lupus_Flynn_Subset8, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
Lupus_Flynn_Subset8 <- NormalizeData(Lupus_Flynn_Subset8)

#identify highly variably features (genes) total of 2000
Lupus_Flynn_Subset8 <- FindVariableFeatures(Lupus_Flynn_Subset8, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
Lupus_Flynn_Subset8 <- ScaleData(Lupus_Flynn_Subset8)

# Perform PCA 
Lupus_Flynn_Subset8 <- RunPCA(Lupus_Flynn_Subset8)

#Cluster the cells
Lupus_Flynn_Subset8 <- FindNeighbors(Lupus_Flynn_Subset8, dims = 1:10)
Lupus_Flynn_Subset8 <- FindClusters(Lupus_Flynn_Subset8, resolution = 0.3)

#Run UMAP
Lupus_Flynn_Subset8 <- RunUMAP(Lupus_Flynn_Subset8, dims = 1:10)


#Dim Plot showing Clusters
DimPlot(Lupus_Flynn_Subset8, reduction = "umap")

saveRDS(Lupus_Flynn_Subset8, "Lupus_Flynn_DN2_cytotoxic_QC.RDS")

Lupus_Flynn_Subset8<-readRDS('Lupus_Flynn_DN2_cytotoxic_QC.RDS')

