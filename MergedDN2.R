#load packages
library(Seurat)
library(dittoSeq)
library(ggpubr)
library(cowplot)
library(ensembldb)
library(AnnotationHub)
library(RCurl)
library(dplyr)
library(SeuratData)
library(patchwork)

setwd("/Users/honokamiura/Downloads/Research/DN2")


DN2_Jimmie_forMerge<-readRDS('DN2_Jimmie_Healthy_Subset11_updated_meta.RDS')
DN2_Powell_Female_forMerge<-readRDS('DN2_Powell_Female_clean_updated_meta.RDS')
DN2_Powell_Male_forMerge<-readRDS('DN2_Powell_Male_clean_updated_meta.RDS')
DN2_AIDA_Female_forMerge<-readRDS('DN2_AIDA_Female_Subset6_updated_meta.RDS')
DN2_AIDA_Male_forMerge<-readRDS('DN2_AIDA_Male_Subset7_updated_meta.RDS')
DN2_Barreiro_forMerge<-readRDS('DN2_Barreiro_Subset8_QC_updated_meta.RDS')



merged_DN2 <- merge(DN2_Jimmie_forMerge, y=c(DN2_Powell_Female_forMerge,DN2_Powell_Male_forMerge,
                                                     DN2_AIDA_Female_forMerge,DN2_AIDA_Male_forMerge,DN2_Barreiro_forMerge))


unique(merged_DN2$Atlas_Study_ID)

#RPCA PIPELINE

# split the dataset into a list of two seurat objects (stim and CTRL)
merged_DN2.list <- SplitObject(merged_DN2, split.by = "Atlas_Study_ID")

# normalize and identify variable features for each dataset independently
merged_DN2.list <- lapply(X = merged_DN2.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = merged_DN2.list)
merged_DN2.list <- lapply(X = merged_DN2.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

immune.anchors <- FindIntegrationAnchors(object.list = merged_DN2.list, anchor.features = features, reduction = "cca")

# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.2)

# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "Atlas_Study_ID")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

DefaultAssay(immune.combined) <- "RNA"

immune.combined<-NormalizeData(immune.combined)

saveRDS(immune.combined, 'DN2_integrated_healthy_11.25.23.rds')
immune.combined <-readRDS('DN2_integrated_healthy_11.25.23.rds')



# Subset Out Cluster 2 (cytotoxic)
clusters_to_combine <- c(0,1,2)
merged_DN2_clean<-subset(x = immune.combined, idents =clusters_to_combine)

#RPCA PIPELINE

# split the dataset into a list of two seurat objects (stim and CTRL)
merged_DN2_clean.list <- SplitObject(merged_DN2_clean, split.by = "Atlas_Study_ID")

# normalize and identify variable features for each dataset independently
merged_DN2_clean.list <- lapply(X = merged_DN2_clean.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = merged_DN2_clean.list)
merged_DN2_clean.list <- lapply(X = merged_DN2_clean.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

immune.anchors.clean <- FindIntegrationAnchors(object.list = merged_DN2_clean.list, anchor.features = features, reduction = "cca")

# this command creates an 'integrated' data assay
immune.combined.clean <- IntegrateData(anchorset = immune.anchors.clean)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined.clean) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined.clean <- ScaleData(immune.combined.clean, verbose = FALSE)
immune.combined.clean <- RunPCA(immune.combined.clean, npcs = 30, verbose = FALSE)
immune.combined.clean <- RunUMAP(immune.combined.clean, reduction = "pca", dims = 1:30)
immune.combined.clean <- FindNeighbors(immune.combined.clean, reduction = "pca", dims = 1:30)
immune.combined.clean <- FindClusters(immune.combined.clean, resolution = 0.2)

# Visualization
p1 <- DimPlot(immune.combined.clean, reduction = "umap", group.by = "Atlas_Study_ID")
p2 <- DimPlot(immune.combined.clean, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

DimPlot(immune.combined.clean, reduction = "umap", group.by = "Atlas_Study_ID")
DimPlot(immune.combined.clean, reduction = "umap", label = TRUE, repel = TRUE)


DefaultAssay(immune.combined.clean) <- "RNA"

immune.combined.clean<-NormalizeData(immune.combined.clean)

saveRDS(immune.combined.clean, 'DN2_integrated_healthy_clean_11.25.23.rds')
immune.combined.clean <-readRDS('DN2_integrated_healthy_clean_11.25.23.rds')



#dn2 markers
DotPlot(object = immune.combined.clean, group.by = 'seurat_clusters', 
        features = c("SOX5",'HLA-DRB5','TBX21',"ITGAX","BATF",'JAZF1','NFATC2',
                     'ZEB2','NR4A2','MS4A1','TOX','FGR','IL10RA','FCRL5',
                     'SREBF2','TFEB','TFEC','ZBTB32','MEF2A','POU2F2','SPI1',
                     #NEGATIVE
                     'LTB','VPREB3','SELL','JCHAIN','CD27','IGHD'),
        cols = c('RdBu') ) + RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 7,face = "bold"))

#suggested
DotPlot(object = immune.combined.clean, group.by = 'seurat_clusters', 
        features = c('ZEB2','CD24','CXCR5'),
        cols = c('RdBu') ) + RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 7,face = "bold"))


#Look for cytotoxic markers 
DotPlot(object = immune.combined, #group.by = 'Clusters',
        features = c('GZMB', 'GZMA','GZMM','PRF1','NKG7','KLRK1','KLRD1','GNLY'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))



DimPlot(object = immune.combined.clean, reduction = "umap", shuffle = F, group.by="Atlas_Sex", 
        cols=c('female'='red',
               'male'='blue'))

unique(immune.combined.clean$Atlas_Age_Category)
DimPlot(object = immune.combined.clean, reduction = "umap", shuffle = F, group.by="Atlas_Age_Category", 
        cols=c('YA'='red',
               'MA'='blue',
               'OA'='green'))

unique(immune.combined$Atlas_Ethnicity_Broad)
DimPlot(object = immune.combined.clean, reduction = "umap", shuffle = F, group.by="Atlas_Ethnicity_Broad", 
        cols=c('EUR'='green',
               'ASI'='red'))

unique(immune.combined$Atlas_Study_ID)
DimPlot(object = immune.combined, reduction = "umap", shuffle = F, group.by="Atlas_Study_ID", 
        cols=c('Jimmie_YE_2022'='snow',
               'Powell_2022'='snow',
               'AIDA'='snow',
               'Barreiro'='red'))

unique(immune.combined$Atlas_Ethnicity_Fine)
DimPlot(object = immune.combined, reduction = "umap", shuffle = F, group.by="Atlas_Ethnicity_Fine", 
        cols=c('EUR'='snow',
               'ASI'='snow',
               'HIS'='snow',
               'Indian'='snow',
               'Singaporean Chinese'='snow',
               'Malaysian'='snow',
               'European'='snow',
               'Japanese'='snow',
               'Korean'='snow',
               'AFR'='black'))


#left=red , right= blue
unique(immune.combined.clean$DN2_Subtypes)
DimPlot(object = immune.combined.clean, reduction = "umap", shuffle = F, group.by="DN2_Subtypes", 
        cols=c('DN2.B'='pink',
               'DN2.A'='cyan',
               'DN2.B.male'='pink',
               'DN2.A.male'='cyan',
               'DN2.AIDA.Female.A'='magenta',
               'DN2.AIDA.Female.B'='navy',
               'DN2.AIDA.Male.B'='navy',
               'DN2.AIDA.Male.A'='magenta',
               'DN2.Jimmie.A'='blue',
               'DN2.Jimmie.B'='red'))

#individual effect
dittoBarPlot(
  object = immune.combined.clean,
  var = "Atlas_Donor", 
  group.by = "seurat_clusters", 
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

#barplot by sex
dittoBarPlot(
  object = immune.combined.clean,
  var = "Atlas_Sex", 
  group.by = "seurat_clusters", 
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

#barplot by ethnicity broad
dittoBarPlot(
  object = Korean_DN2s,
  var = "Atlas_Donor", 
  group.by = "seurat_clusters", 
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))&NoLegend()

#barplot by ethnicity broad
dittoBarPlot(
  object = immune.combined.clean,
  var = "seurat_clusters", 
  group.by = "Atlas_Ethnicity_Broad", 
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))


#barplot by ethnicity fine
dittoBarPlot(
  object = immune.combined.clean,
  var = "Atlas_Ethnicity_Fine", 
  group.by = "seurat_clusters", 
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

#barplot by study id
dittoBarPlot(
  object = immune.combined.clean,
  var = "Atlas_Study_ID", 
  group.by = "seurat_clusters", 
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

#barplot by age
dittoBarPlot(
  object = immune.combined.clean,
  var = "Atlas_Age_Category", 
  group.by = "seurat_clusters", 
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

#barplot by age
dittoBarPlot(
  object = immune.combined.clean,
  var = "seurat_clusters", 
  group.by = "Atlas_Age_Category", 
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))


#transcriptional differences
#compare 0 to another
DN2_merged_Cluster0_Markers <- FindMarkers(immune.combined.clean, only.pos = F,
                                                ident.1 = ('0'),
                                                ident.2 = c('1','2','3'),
                                                min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_merged_Cluster0_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_merged_Cluster0_Markers.csv")

#compare 1 to another
DN2_merged_Cluster1_Markers <- FindMarkers(immune.combined.clean, only.pos = F,
                                           ident.1 = ('1'),
                                           ident.2 = c('0','2','3'),
                                           min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_merged_Cluster1_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_merged_Cluster1_Markers.csv")

#compare 2 to another
DN2_merged_Cluster2_Markers <- FindMarkers(immune.combined.clean, only.pos = F,
                                           ident.1 = ('2'),
                                           ident.2 = c('0','1','3'),
                                           min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_merged_Cluster2_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_merged_Cluster2_Markers.csv")


#Age Updated
META<-immune.combined.clean[[]]
write.csv(META,'Powell_Male_cell_id_meta.csv')
Powell_Male_cell_id_meta<-read.csv('Powell_Male_cell_id_meta.csv', row.names = NULL)

donor_id_meta<-(META%>%distinct(Atlas_Donor, .keep_all=TRUE))
write.csv(donor_id_meta,'integratedDN2s_donor_id_meta.csv')
donor_id_meta<-read.csv('Powell_Male_donor_id_meta.csv')


unique(immune.combined.clean$Atlas_Ethnicity_Fine)


#Korean population analysis- ethnicity fine

Idents(immune.combined.clean)<-'Atlas_Ethnicity_Fine'

Korean_DN2s<-subset(immune.combined.clean, idents=c('Korean'))








