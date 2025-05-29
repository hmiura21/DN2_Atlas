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
library(Azimuth)

setwd("/Users/honokamiura/Downloads/Research/DN2")

Lupus_Webb <-readRDS('/Users/honokamiura/Downloads/Research/DN2/GSE189050_final_seurat.rds')

META<-Lupus_Webb[[]]

unique(Lupus_Webb$classification)

DefaultAssay(Lupus_Webb) <- "RNA"
DefaultAssay(Bcell_Lupus_Webb) <- "RNA"



#violin plot to look at percent.mt, nFeature_RNA
VlnPlot(Lupus_Webb, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3)

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
Lupus_Webb <- subset(Lupus_Webb, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent_mt < 5)

#Normalization at default
Lupus_Webb <- NormalizeData(Lupus_Webb)

#identify highly variably features (genes) total of 2000
Lupus_Webb <- FindVariableFeatures(Lupus_Webb, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
Lupus_Webb <- ScaleData(Lupus_Webb)

# Perform PCA and color by cell cycle phase
Lupus_Webb <- RunPCA(Lupus_Webb)

#Cluster the cells
Lupus_Webb <- FindNeighbors(Lupus_Webb, dims = 1:10)
Lupus_Webb <- FindClusters(Lupus_Webb, resolution = 0.5)

#Run UMAP
Lupus_Webb <- RunUMAP(Lupus_Webb, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(Lupus_Webb)

saveRDS(Lupus_Webb, "Lupus_Webb_QC.RDS")

#Start here: 
Lupus_Webb<-readRDS('Lupus_Webb_QC.RDS')


#Look for only Significant DN2 markers
DotPlot(object = Lupus_Webb, #group.by = 'Clusters',
        features = c('CD19','CD20','CD21','CD22'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

BiocManager::install("TFBSTools", type = "source", force = TRUE)

#Looking for Tbet+,
DotPlot(object = Lupus_Webb, #best=2,0 (maybe 1,3)
        features = c('TBX21'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))


library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)


# The RunAzimuth function can take a Seurat object as input
Lupus_Webb <- RunAzimuth(Lupus_Webb, reference = "pbmcref")
DimPlot(Lupus_Webb, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()


DimPlot(object = Lupus_Webb, reduction = "umap", shuffle = F,
        cols= c('0'= 'snow2',
                '1'='snow2',
                '2'='snow2',
                '3'='pink',
                '4'='snow2',
                '5'='snow2', 
                '6'='snow2',
                '7'='snow2', 
                '8'='snow2',
                '9'='snow2',
                '10'='snow2',
                '11'='red'))

clusters_to_combine <- c(1,6)
Bcell_Lupus_Webb<-subset(x = Lupus_Webb, idents =clusters_to_combine)

DimPlot(Bcell_Lupus_Webb)

#QC %MT^
Bcell_Lupus_Webb[["percent.mt"]] <- PercentageFeatureSet(Bcell_Lupus_Webb, pattern = "^MT-")

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
Bcell_Lupus_Webb <- subset(Bcell_Lupus_Webb, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
Bcell_Lupus_Webb <- NormalizeData(Bcell_Lupus_Webb)

#identify highly variably features (genes) total of 2000
Bcell_Lupus_Webb <- FindVariableFeatures(Bcell_Lupus_Webb, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
Bcell_Lupus_Webb <- ScaleData(Bcell_Lupus_Webb)

# Perform PCA 
Bcell_Lupus_Webb <- RunPCA(Bcell_Lupus_Webb)

#Cluster the cells
Bcell_Lupus_Webb <- FindNeighbors(Bcell_Lupus_Webb, dims = 1:10)
Bcell_Lupus_Webb <- FindClusters(Bcell_Lupus_Webb, resolution = 0.5)

#Run UMAP
Bcell_Lupus_Webb <- RunUMAP(Bcell_Lupus_Webb, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(Bcell_Lupus_Webb, reduction = "umap")

saveRDS(Bcell_Lupus_Webb, "Bcell_Lupus_Webb_QC.RDS")

Bcell_Lupus_Webb<-readRDS('Bcell_Lupus_Webb_QC.RDS')


#dn2 markers
DotPlot(object = Bcell_Lupus_Webb, group.by = 'seurat_clusters', 
        features = c("SOX5",'HLA-DRB5','TBX21',"ITGAX","BATF",'JAZF1','NFATC2',
                     'ZEB2','NR4A2','MS4A1','TOX','FGR','IL10RA','FCRL5',
                     'SREBF2','TFEB','TFEC','ZBTB32','MEF2A','POU2F2','SPI1',
                     #NEGATIVE
                     'VPREB3','SELL','CD27','IGHD'),
        cols = c('RdBu') ) + RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 7,face = "bold"))
#Looking for Tbet+,
DotPlot(object = Bcell_Lupus_Webb, 
        features = c('TBX21'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Looking for Key characteristics+,
DotPlot(object = Bcell_Lupus_Webb, #best= 0,1,2,3,4,5
        features = c( 'CD27'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))
DotPlot(object = Bcell_Lupus_Webb, #best=2,5
        features = c( 'IGHD'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))


#Look for cytotoxic markers 
DotPlot(object = Bcell_Lupus_Webb, #best= 5 (maybe 3)
        features = c('GZMB', 'GZMA','GZMM','PRF1','NKG7','KLRK1','KLRD1','GNLY'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))



# Markers to Show for DN2: Density Plot
markers <- list()
markers$cytotoxic<- c('GZMB', 'GZMA','GZMM','PRF1','NKG7','KLRK1','KLRD1','GNLY')
markers$DN2_pos<- c("SOX5",'HLA-DRB5','TBX21',"ITGAX","BATF",'JAZF1','NFATC2',
                'ZEB2','NR4A2','MS4A1','TOX','FGR','IL10RA','FCRL5',
                'SREBF2','TFEB','TFEC','ZBTB32','MEF2A','POU2F2','SPI1'
                )
markers$DN2_neg<- c('VPREB3','SELL','CD27','IGHD')
markers$tbet<- c('TBX21')


Bcell_Lupus_Webb <- AddModuleScore_UCell(Bcell_Lupus_Webb,features = markers)

#Marker cytotoxic Score
Bcell_Lupus_Webb$cytotoxic_UCell
plot_density(
  Bcell_Lupus_Webb,
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
Bcell_Lupus_Webb$DN2_pos_UCell
plot_density(
  Bcell_Lupus_Webb,
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
Bcell_Lupus_Webb$DN2_neg_UCell
plot_density(
  Bcell_Lupus_Webb,
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
Bcell_Lupus_Webb$tbet_UCell
plot_density(
  Bcell_Lupus_Webb,
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


#ADD METADATA

#FOR WEBB B CELLS
#add study id
Bcell_Lupus_Webb@meta.data$Atlas_Study_ID<-'Webb_Lupus'
Bcell_Lupus_Webb@meta.data$Atlas_Sex<-'female'
META<-Bcell_Lupus_Webb[[]]

#Disease
unique(Bcell_Lupus_Webb$classification)
Idents(Bcell_Lupus_Webb)<-'classification'
current.cluster.ids<-levels(Bcell_Lupus_Webb@active.ident)
current.cluster.ids<-c( "SLE INACT"   ,                "SLE ACT"  ,                   
                        "Control")
new.cluster.ids<-c( "SLE Inactive"   ,                "SLE Active"  ,                   
                    "Healthy")

Bcell_Lupus_Webb@meta.data$Atlas_Disease<-Bcell_Lupus_Webb@meta.data$classification
Bcell_Lupus_Webb@meta.data$Atlas_Disease<-plyr::mapvalues(x=Bcell_Lupus_Webb@meta.data$Atlas_Disease, from=current.cluster.ids, to=new.cluster.ids)

unique(Bcell_Lupus_Webb$Atlas_Disease)


#Ethnicity- Broad
unique(Bcell_Lupus_Webb$ancestry)
Idents(Bcell_Lupus_Webb)<-'ancestry'
current.cluster.ids<-levels(Bcell_Lupus_Webb@active.ident)
current.cluster.ids<-c( "EA"   ,                "AA")
new.cluster.ids<-c( "EUR"   ,                "AFR")

Bcell_Lupus_Webb@meta.data$Atlas_Ethnicity_Broad<-Bcell_Lupus_Webb@meta.data$ancestry
Bcell_Lupus_Webb@meta.data$Atlas_Ethnicity_Broad<-plyr::mapvalues(x=Bcell_Lupus_Webb@meta.data$Atlas_Ethnicity_Broad, from=current.cluster.ids, to=new.cluster.ids)

unique(Bcell_Lupus_Webb$Atlas_Ethnicity_Broad)

#Ethnicity- Fine
Bcell_Lupus_Webb@meta.data$Atlas_Ethnicity_Fine<-Bcell_Lupus_Webb@meta.data$ancestry
Bcell_Lupus_Webb@meta.data$Atlas_Ethnicity_Fine<-plyr::mapvalues(x=Bcell_Lupus_Webb@meta.data$Atlas_Ethnicity_Fine, from=current.cluster.ids, to=new.cluster.ids)

unique(Bcell_Lupus_Webb$Atlas_Ethnicity_Fine)

#Atlas Donor
unique(Bcell_Lupus_Webb$subject_id)
Idents(Bcell_Lupus_Webb)<-'subject_id'
current.cluster.ids<-levels(Bcell_Lupus_Webb@active.ident)
current.cluster.ids
new.cluster.ids<-levels(Bcell_Lupus_Webb@active.ident)

Bcell_Lupus_Webb@meta.data$Atlas_Donor<-Bcell_Lupus_Webb@meta.data$'subject_id'
Bcell_Lupus_Webb@meta.data$Atlas_Donor<-plyr::mapvalues(x=Bcell_Lupus_Webb@meta.data$Atlas_Donor, from=current.cluster.ids, to=new.cluster.ids)

unique(Bcell_Lupus_Webb$Atlas_Donor)

#Age Updated
META<-Bcell_Lupus_Webb[[]]
write.csv(META,'Bcell_Lupus_Webb_cell_id_meta.csv')
Bcell_Lupus_Webb_cell_id_meta<-read.csv('Bcell_Lupus_Webb_cell_id_meta.csv', row.names = NULL)

donor_id_meta<-(META%>%distinct(Atlas_Donor, .keep_all=TRUE))
write.csv(donor_id_meta,'Bcell_Lupus_Webb_donor_id_meta.csv')
donor_id_meta<-read.csv('Bcell_Lupus_Webb_donor_id_meta.csv')

merged_df<- merge(donor_id_meta, Bcell_Lupus_Webb_cell_id_meta, by= 'Atlas_Donor', all=TRUE)
#check to see which column has the cell ids/barcodes
write.csv(merged_df,'merged_df.csv')
merged_df2<-read.csv('merged_df.csv', row.names = 5) #row.names = column with cell barcodes/ids
drop <- c("X.1")
merged_df2 = merged_df2[,!(names(merged_df2) %in% drop)]


Bcell_Lupus_Webb<-AddMetaData(Bcell_Lupus_Webb,merged_df2 )
unique(Bcell_Lupus_Webb$Atlas_Age_Category)
saveRDS(Bcell_Lupus_Webb, "Bcell_Lupus_Webb_updated_meta.rds")
Bcell_Lupus_Webb<-readRDS('Bcell_Lupus_Webb_updated_meta.rds')



#PROP ANALYSIS

# DN2 Ratio
cluster_2_count <- sum(Bcell_Lupus_Webb@meta.data$seurat_clusters == "2")
all_clusters_count <- table(Bcell_Lupus_Webb@meta.data$seurat_clusters)
total_cells <- sum(all_clusters_count)

unique(Bcell_Lupus_Webb$Atlas_Age_Category)

dittoBarPlot(
  object = Bcell_Lupus_Webb,
  var = "Atlas_Disease",    #colors stacked
  group.by = "seurat_clusters",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

dittoBarPlot(
  object = Bcell_Lupus_Webb,
  var = "Atlas_Disease",    #colors stacked
  group.by = "Atlas_Cytotoxic_Cluster",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))


#ADD METADATA for DN2 cluster and other clusters
sort(unique(Bcell_Lupus_Webb$seurat_clusters))

unique(Bcell_Lupus_Webb$seurat_clusters)
Idents(Bcell_Lupus_Webb)<-'seurat_clusters'
current.cluster.ids<-levels(Bcell_Lupus_Webb@active.ident)
current.cluster.ids<-c( '0','2','3','1','4')
new.cluster.ids<-c('other','other','DN2','other','other' )

Bcell_Lupus_Webb@meta.data$Atlas_DN2_Cluster<-Bcell_Lupus_Webb@meta.data$seurat_clusters
Bcell_Lupus_Webb@meta.data$Atlas_DN2_Cluster<-plyr::mapvalues(x=Bcell_Lupus_Webb@meta.data$Atlas_DN2_Cluster, from=current.cluster.ids, to=new.cluster.ids)

unique(Bcell_Lupus_Webb$Atlas_DN2_Cluster)


#ADD METADATA for cytotoxic cluster and other clusters
unique(Bcell_Lupus_Webb$seurat_clusters)
Idents(Bcell_Lupus_Webb)<-'seurat_clusters'
current.cluster.ids<-levels(Bcell_Lupus_Webb@active.ident)
current.cluster.ids<-c( '0','2','3','1','4')
new.cluster.ids<-c('other','other','other','other','cytotoxic')

Bcell_Lupus_Webb@meta.data$Atlas_Cytotoxic_Cluster<-Bcell_Lupus_Webb@meta.data$seurat_clusters
Bcell_Lupus_Webb@meta.data$Atlas_Cytotoxic_Cluster<-plyr::mapvalues(x=Bcell_Lupus_Webb@meta.data$Atlas_Cytotoxic_Cluster, from=current.cluster.ids, to=new.cluster.ids)

unique(Bcell_Lupus_Webb$Atlas_Cytotoxic_Cluster)




#reorder age x axis labeling
Bcell_Lupus_Webb$Atlas_Age_Category <- factor(Bcell_Lupus_Webb$Atlas_Age_Category, levels = c("YA", "MA", 'OA'))
Bcell_Lupus_Webb$Atlas_Disease <- factor(Bcell_Lupus_Webb$Atlas_Disease, levels = c("Healthy", "SLE Inactive", 'SLE Active'))


#DN2 proportion by age
dittoBarPlot(
  object = Bcell_Lupus_Webb,
  var = "Atlas_DN2_Cluster",    #colors stacked
  group.by = "Atlas_Age_Category",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

#DN2 proportion by disease
dittoBarPlot(
  object = Bcell_Lupus_Webb,
  var = "Atlas_DN2_Cluster",    #colors stacked
  group.by = "Atlas_Disease",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

#DN2 proportion by ethnicity
dittoBarPlot(
  object = Bcell_Lupus_Webb,
  var = "Atlas_DN2_Cluster",    #colors stacked
  group.by = "Atlas_Ethnicity_Broad",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))


#Cytotoxic proportion by age
dittoBarPlot(
  object = Bcell_Lupus_Webb,
  var = "Atlas_Cytotoxic_Cluster",    #colors stacked
  group.by = "Atlas_Age_Category",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

#Cytotoxic proportion by disease
dittoBarPlot(
  object = Bcell_Lupus_Webb,
  var = "Atlas_Cytotoxic_Cluster",    #colors stacked
  group.by = "Atlas_Disease",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

#Cytotoxic proportion by ethnicity
dittoBarPlot(
  object = Bcell_Lupus_Webb,
  var = "Atlas_Cytotoxic_Cluster",    #colors stacked
  group.by = "Atlas_Ethnicity_Broad",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

unique(Bcell_Lupus_Webb$clusters_annotated)


#DN2 proportion by clusters annotated
dittoBarPlot(
  object = Bcell_Lupus_Webb,
  var = "clusters_annotated",    #colors stacked
  group.by = "seurat_clusters",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

dittoBarPlot(
  object = Bcell_Lupus_Webb,
  var = "seurat_clusters",    #colors stacked
  group.by = "clusters_annotated",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))




#Statistical Significance
library("scProportionTest")

prop_test <- sc_utils(Bcell_Lupus_Webb)


#prop test for DN2 cluster vs age 
prop_test_1_vs_2 <- permutation_test(
  prop_test, cluster_identity = "Atlas_DN2_Cluster",
  sample_1 = "YA", sample_2 = "MA",
  sample_identity = "Atlas_Age_Category"
)
permutation_plot(prop_test_1_vs_2)

prop_test_1_vs_3 <- permutation_test(
  prop_test, cluster_identity = "Atlas_DN2_Cluster",
  sample_1 = "YA", sample_2 = "OA",
  sample_identity = "Atlas_Age_Category"
)
permutation_plot(prop_test_1_vs_3)

prop_test_2_vs_3 <- permutation_test(
  prop_test, cluster_identity = "Atlas_DN2_Cluster",
  sample_1 = "MA", sample_2 = "OA",
  sample_identity = "Atlas_Age_Category"
)
permutation_plot(prop_test_2_vs_3)

#prop test for DN2 cluster vs disease 
prop_test_1_vs_2 <- permutation_test(
  prop_test, cluster_identity = "Atlas_DN2_Cluster",
  sample_1 = "Healthy", sample_2 = "SLE Inactive",
  sample_identity = "Atlas_Disease"
)
permutation_plot(prop_test_1_vs_2)

prop_test_1_vs_3 <- permutation_test(
  prop_test, cluster_identity = "Atlas_DN2_Cluster",
  sample_1 = "Healthy", sample_2 = "SLE Active",
  sample_identity = "Atlas_Disease"
)
permutation_plot(prop_test_1_vs_3)

prop_test_2_vs_3 <- permutation_test(
  prop_test, cluster_identity = "Atlas_DN2_Cluster",
  sample_1 = "SLE Inactive", sample_2 = "SLE Active",
  sample_identity = "Atlas_Disease"
)
permutation_plot(prop_test_2_vs_3)

#Prop test for DN2 cluster vs ethnicity
prop_test_1_vs_2 <- permutation_test(
  prop_test, cluster_identity = "Atlas_DN2_Cluster",
  sample_1 = "EUR", sample_2 = "AFR",
  sample_identity = "Atlas_Ethnicity_Broad"
)
permutation_plot(prop_test_1_vs_2)



#prop test for cytotoxic cluster vs age 
prop_test_1_vs_2 <- permutation_test(
  prop_test, cluster_identity = "Atlas_Cytotoxic_Cluster",
  sample_1 = "YA", sample_2 = "MA",
  sample_identity = "Atlas_Age_Category"
)
permutation_plot(prop_test_1_vs_2)

prop_test_1_vs_3 <- permutation_test(
  prop_test, cluster_identity = "Atlas_Cytotoxic_Cluster",
  sample_1 = "YA", sample_2 = "OA",
  sample_identity = "Atlas_Age_Category"
)
permutation_plot(prop_test_1_vs_3)

prop_test_2_vs_3 <- permutation_test(
  prop_test, cluster_identity = "Atlas_Cytotoxic_Cluster",
  sample_1 = "MA", sample_2 = "OA",
  sample_identity = "Atlas_Age_Category"
)
permutation_plot(prop_test_2_vs_3)

#prop test for cytotoxic cluster vs disease 
prop_test_1_vs_2 <- permutation_test(
  prop_test, cluster_identity = "Atlas_Cytotoxic_Cluster",
  sample_1 = "Healthy", sample_2 = "SLE Inactive",
  sample_identity = "Atlas_Disease"
)
permutation_plot(prop_test_1_vs_2)

prop_test_1_vs_3 <- permutation_test(
  prop_test, cluster_identity = "Atlas_Cytotoxic_Cluster",
  sample_1 = "Healthy", sample_2 = "SLE Active",
  sample_identity = "Atlas_Disease"
)
permutation_plot(prop_test_1_vs_3)

prop_test_2_vs_3 <- permutation_test(
  prop_test, cluster_identity = "Atlas_Cytotoxic_Cluster",
  sample_1 = "SLE Inactive", sample_2 = "SLE Active",
  sample_identity = "Atlas_Disease"
)
permutation_plot(prop_test_2_vs_3)

#Prop test for cytotoxic cluster vs ethnicity
prop_test_1_vs_2 <- permutation_test(
  prop_test, cluster_identity = "Atlas_Cytotoxic_Cluster",
  sample_1 = "EUR", sample_2 = "AFR",
  sample_identity = "Atlas_Ethnicity_Broad"
)
permutation_plot(prop_test_1_vs_2)




# Subset Cluster 3 (DN2) 
Bcell_Lupus_Webb_Subset3<-subset(x = Bcell_Lupus_Webb, idents ='3')

#QC %MT^
Bcell_Lupus_Webb_Subset3[["percent.mt"]] <- PercentageFeatureSet(Bcell_Lupus_Webb_Subset3, pattern = "^MT-")

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
Bcell_Lupus_Webb_Subset3 <- subset(Bcell_Lupus_Webb_Subset3, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
Bcell_Lupus_Webb_Subset3 <- NormalizeData(Bcell_Lupus_Webb_Subset3)

#identify highly variably features (genes) total of 2000
Bcell_Lupus_Webb_Subset3 <- FindVariableFeatures(Bcell_Lupus_Webb_Subset3, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
Bcell_Lupus_Webb_Subset3 <- ScaleData(Bcell_Lupus_Webb_Subset3)

# Perform PCA 
Bcell_Lupus_Webb_Subset3 <- RunPCA(Bcell_Lupus_Webb_Subset3)

#Cluster the cells
Bcell_Lupus_Webb_Subset3 <- FindNeighbors(Bcell_Lupus_Webb_Subset3, dims = 1:10)
Bcell_Lupus_Webb_Subset3 <- FindClusters(Bcell_Lupus_Webb_Subset3, resolution = 0.5)

#Run UMAP
Bcell_Lupus_Webb_Subset3 <- RunUMAP(Bcell_Lupus_Webb_Subset3, dims = 1:10)


#Dim Plot showing Clusters
DimPlot(Bcell_Lupus_Webb_Subset3, reduction = "umap")

saveRDS(Bcell_Lupus_Webb_Subset3, "Lupus_Webb_DN2_QC.RDS")

Bcell_Lupus_Webb_Subset3<-readRDS('Lupus_Webb_DN2_QC.RDS')





# Subset Cluster 4 (cytotoxic) 
Bcell_Lupus_Webb_Subset4<-subset(x = Bcell_Lupus_Webb, idents ='4')

#QC %MT^
Bcell_Lupus_Webb_Subset4[["percent.mt"]] <- PercentageFeatureSet(Bcell_Lupus_Webb_Subset4, pattern = "^MT-")

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
Bcell_Lupus_Webb_Subset4 <- subset(Bcell_Lupus_Webb_Subset4, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
Bcell_Lupus_Webb_Subset4 <- NormalizeData(Bcell_Lupus_Webb_Subset4)

#identify highly variably features (genes) total of 2000
Bcell_Lupus_Webb_Subset4 <- FindVariableFeatures(Bcell_Lupus_Webb_Subset4, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
Bcell_Lupus_Webb_Subset4 <- ScaleData(Bcell_Lupus_Webb_Subset4)

# Perform PCA 
Bcell_Lupus_Webb_Subset4 <- RunPCA(Bcell_Lupus_Webb_Subset4)

#Cluster the cells
Bcell_Lupus_Webb_Subset4 <- FindNeighbors(Bcell_Lupus_Webb_Subset4, dims = 1:10)
Bcell_Lupus_Webb_Subset4 <- FindClusters(Bcell_Lupus_Webb_Subset4, resolution = 0.5)

#Run UMAP
Bcell_Lupus_Webb_Subset4 <- RunUMAP(Bcell_Lupus_Webb_Subset4, dims = 1:10)


#Dim Plot showing Clusters
DimPlot(Bcell_Lupus_Webb_Subset4, reduction = "umap")

saveRDS(Bcell_Lupus_Webb_Subset4, "Lupus_Webb_cytotoxic_QC.RDS")

Bcell_Lupus_Webb_Subset4<-readRDS('Lupus_Webb_cytotoxic_QC.RDS')
