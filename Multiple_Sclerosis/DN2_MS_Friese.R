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

MS_Friese <-readRDS('/Users/honokamiura/Downloads/Research/DN2/GSE144744_MS_B_cells.rds')

META<-MS_Friese[[]]

DefaultAssay(MS_Friese) <- "RNA"

#QC %MT^
MS_Friese[["percent.mt"]] <- PercentageFeatureSet(MS_Friese, pattern = "^MT.")

#violin plot to look at percent.mt, nFeature_RNA
VlnPlot(MS_Friese, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
MS_Friese <- subset(MS_Friese, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
MS_Friese <- NormalizeData(MS_Friese)

#identify highly variably features (genes) total of 2000
MS_Friese <- FindVariableFeatures(MS_Friese, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
MS_Friese <- ScaleData(MS_Friese)

# Perform PCA and color by cell cycle phase
MS_Friese <- RunPCA(MS_Friese)

#Cluster the cells
MS_Friese <- FindNeighbors(MS_Friese, dims = 1:10)
MS_Friese <- FindClusters(MS_Friese, resolution = 0.9)

#Run UMAP
MS_Friese <- RunUMAP(MS_Friese, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(MS_Friese)

saveRDS(MS_Friese, "MS_Friese_QC.RDS")

#Start here: 
MS_Friese<-readRDS('MS_Friese_QC.RDS')


#Looking for Tbet+, 
DotPlot(object = MS_Friese, #best=3
        features = c('TBX21'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Look for cytotoxic markers best=0
DotPlot(object = MS_Friese, #group.by = 'Clusters',
        features = c('GZMB', 'GZMA','GZMM','PRF1','NKG7','KLRK1','KLRD1','GNLY'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#dn2 markers (1)
DotPlot(object = MS_Friese, group.by = 'seurat_clusters', 
        features = c("SOX5",'HLA-DRB5','TBX21',"ITGAX","BATF",'JAZF1','NFATC2',
                     'ZEB2','NR4A2','MS4A1','TOX','FGR','IL10RA','FCRL5',
                     'SREBF2','TFEB','TFEC','ZBTB32','MEF2A','POU2F2','SPI1',
                     #NEGATIVE
                     'LTB','VPREB3','SELL','JCHAIN','CD27','IGHD'),
        cols = c('RdBu') ) + RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 7,face = "bold"))


unique(MS_Friese$Atlas_Donor)

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


MS_Friese <- AddModuleScore_UCell(MS_Friese,features = markers)

#Marker cytotoxic Score
MS_Friese$cytotoxic_UCell
plot_density(
  MS_Friese,
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
MS_Friese$DN2_pos_UCell
plot_density(
  MS_Friese,
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
MS_Friese$DN2_neg_UCell
plot_density(
  MS_Friese,
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
MS_Friese$tbet_UCell
plot_density(
  MS_Friese,
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
  object = MS_Friese,
  var = "clusters_annotated",    #colors stacked
  group.by = "seurat_clusters",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))



#ADD METADATA

#FOR WEBB B CELLS
#add study id
MS_Friese@meta.data$Atlas_Study_ID<-'Friese_MS'

META<-MS_Friese[[]]

#sex, donor, study, disease ;;;; need ethnicity, age range
unique(MS_Friese$Atlas_Study_ID)

saveRDS(MS_Friese, "MS_Friese_updated_meta.rds")
MS_Friese<-readRDS('MS_Friese_updated_meta.RDS')

#Age Updated
META<-MS_Friese[[]]
write.csv(META,'MS_Friese_cell_id_meta.csv')
MS_Friese_cell_id_meta<-read.csv('MS_Friese_cell_id_meta.csv', row.names = NULL)

donor_id_meta<-(META%>%distinct(Atlas_Donor, .keep_all=TRUE))
write.csv(donor_id_meta,'MS_Friese_donor_id_meta.csv')
donor_id_meta<-read.csv('MS_Friese_donor_id_meta.csv')

merged_df<- merge(donor_id_meta, MS_Friese_cell_id_meta, by= 'Atlas_Donor', all=TRUE)
#check to see which column has the cell ids/barcodes
write.csv(merged_df,'merged_df.csv')
merged_df2<-read.csv('merged_df.csv', row.names = 5) #row.names = column with cell barcodes/ids
drop <- c("X.1")
merged_df2 = merged_df2[,!(names(merged_df2) %in% drop)]


MS_Friese<-AddMetaData(MS_Friese,merged_df2 )
unique(MS_Friese$Atlas_Age_Category)
saveRDS(MS_Friese, "MS_Friese_updated_meta.rds")
MS_Friese<-readRDS('MS_Friese_updated_meta.RDS')



#PROP ANALYSIS


#ADD METADATA for DN2 cluster and other clusters
unique(MS_Friese$seurat_clusters)
Idents(MS_Friese)<-'seurat_clusters'
current.cluster.ids<-levels(MS_Friese@active.ident)
current.cluster.ids<-c( '0','1','2','3','4')
new.cluster.ids<-c('DN2','other','other','other','other')

MS_Friese@meta.data$Atlas_DN2_Cluster<-MS_Friese@meta.data$seurat_clusters
MS_Friese@meta.data$Atlas_DN2_Cluster<-plyr::mapvalues(x=MS_Friese@meta.data$Atlas_DN2_Cluster, from=current.cluster.ids, to=new.cluster.ids)

unique(MS_Friese$Atlas_DN2_Cluster)


#ADD METADATA for cytotoxic cluster and other clusters
unique(MS_Friese$seurat_clusters)
Idents(MS_Friese)<-'seurat_clusters'
current.cluster.ids<-levels(MS_Friese@active.ident)
current.cluster.ids<-c('0','1','2','3','4')
new.cluster.ids<-c('other','other','other','other','cytotoxic')

MS_Friese@meta.data$Atlas_Cytotoxic_Cluster<-MS_Friese@meta.data$seurat_clusters
MS_Friese@meta.data$Atlas_Cytotoxic_Cluster<-plyr::mapvalues(x=MS_Friese@meta.data$Atlas_Cytotoxic_Cluster, from=current.cluster.ids, to=new.cluster.ids)

unique(MS_Friese$Atlas_Cytotoxic_Cluster)


#reorder age x axis labeling
MS_Friese$Atlas_Age_Category <- factor(MS_Friese$Atlas_Age_Category, levels = c("YA", "MA"))
MS_Friese$Atlas_Disease <- factor(MS_Friese$Atlas_Disease, levels = c("Healthy",'MS'))



#DN2 proportion 
dittoBarPlot(
  object = MS_Friese,
  var = "Atlas_DN2_Cluster",    #colors stacked
  group.by = "Atlas_Disease",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

#DN2 proportion by age
dittoBarPlot(
  object = MS_Friese,
  var = "Atlas_DN2_Cluster",    #colors stacked
  group.by = "Atlas_Age_Category",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

#DN2 proportion by age
dittoBarPlot(
  object = MS_Friese,
  var = "Atlas_Age_Category",    #colors stacked
  group.by = "Atlas_DN2_Cluster",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

#Cytotoxic proportion by age
dittoBarPlot(
  object = MS_Friese,
  var = "Atlas_Cytotoxic_Cluster",    #colors stacked
  group.by = "Atlas_Age_Category",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

#Cytotoxic proportion by disease
dittoBarPlot(
  object = MS_Friese,
  var = "Atlas_Cytotoxic_Cluster",    #colors stacked
  group.by = "Atlas_Disease",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))


#Statistical Significance
library("scProportionTest")

prop_test <- sc_utils(MS_Friese)

#prop test for DN2 cluster vs age 
prop_test_1_vs_2 <- permutation_test(
  prop_test, cluster_identity = "Atlas_DN2_Cluster",
  sample_1 = "adult", sample_2 = "child",
  sample_identity = "Age"
)
permutation_plot(prop_test_1_vs_2)






# Subset Cluster 0 (DN2) 
MS_Friese_Subset0<-subset(x = MS_Friese, idents ='0')

#QC %MT^
MS_Friese_Subset0[["percent.mt"]] <- PercentageFeatureSet(MS_Friese_Subset0, pattern = "^MT-")

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
MS_Friese_Subset0 <- subset(MS_Friese_Subset0, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
MS_Friese_Subset0 <- NormalizeData(MS_Friese_Subset0)

#identify highly variably features (genes) total of 2000
MS_Friese_Subset0 <- FindVariableFeatures(MS_Friese_Subset0, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
MS_Friese_Subset0 <- ScaleData(MS_Friese_Subset0)

# Perform PCA 
MS_Friese_Subset0 <- RunPCA(MS_Friese_Subset0)

#Cluster the cells
MS_Friese_Subset0 <- FindNeighbors(MS_Friese_Subset0, dims = 1:10)
MS_Friese_Subset0 <- FindClusters(MS_Friese_Subset0, resolution = 0.5)

#Run UMAP
MS_Friese_Subset0 <- RunUMAP(MS_Friese_Subset0, dims = 1:10)


#Dim Plot showing Clusters
DimPlot(MS_Friese_Subset0, reduction = "umap")

saveRDS(MS_Friese_Subset0, "MS_Friese_DN2_QC.RDS")

MS_Friese_Subset0<-readRDS('MS_Friese_DN2_QC.RDS')





# Subset Cluster 4 (cytotoxic) 
MS_Friese_Subset4<-subset(x = MS_Friese, idents ='4')

#QC %MT^
MS_Friese_Subset4[["percent.mt"]] <- PercentageFeatureSet(MS_Friese_Subset4, pattern = "^MT-")

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
MS_Friese_Subset4 <- subset(MS_Friese_Subset4, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
MS_Friese_Subset4 <- NormalizeData(MS_Friese_Subset4)

#identify highly variably features (genes) total of 2000
MS_Friese_Subset4 <- FindVariableFeatures(MS_Friese_Subset4, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
MS_Friese_Subset4 <- ScaleData(MS_Friese_Subset4)

# Perform PCA 
MS_Friese_Subset4 <- RunPCA(MS_Friese_Subset4)

#Cluster the cells
MS_Friese_Subset4 <- FindNeighbors(MS_Friese_Subset4, dims = 1:10)
MS_Friese_Subset4 <- FindClusters(MS_Friese_Subset4, resolution = 0.2)

#Run UMAP
MS_Friese_Subset4 <- RunUMAP(MS_Friese_Subset4, dims = 1:10)


#Dim Plot showing Clusters
DimPlot(MS_Friese_Subset4, reduction = "umap")

saveRDS(MS_Friese_Subset4, "MS_Friese_Subset4_cytotoxic_QC.RDS")

MS_Friese_Subset4<-readRDS('MS_Friese_Subset4_cytotoxic_QC.RDS')
