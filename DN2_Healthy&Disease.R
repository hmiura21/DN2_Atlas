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


install.packages("devtools")
remotes::install_github("bnprks/BPCells", quiet = TRUE)
library(BPCells)
library(HDF5Array)

#DN2 Clusters
DN2_Healthy<-readRDS("MergedB_DN2_QC.RDS")
DN2_Lupus_Webb<-readRDS("Lupus_Webb_DN2_QC.RDS")
DN2_Lupus_Flynn<-readRDS("Lupus_Flynn_DN2_cytotoxic_QC.RDS") #same as cytotoxic
DN2_MS_Friese<-readRDS("MS_Friese_DN2_QC.RDS")
DN2_Crohns_Martin<-readRDS("CD_Martin_DN2_QC.RDS")
DN2_Covid_Youngs_Wilk<-readRDS("Covid_Youngs_Wilk_DN2_QC.RDS")
DN2_Malaria_Holla<-readRDS("Malaria_Holla_DN2_QC.RDS") 
DN2_HIV_Holla<-readRDS("HIV_Holla_Subset07_DN2_QC.RDS") #2 clusters merged together

#Cytotoxic Clusters
cytotoxic_Healthy<-readRDS("MergedB_cytotoxic_QC.RDS")
cytotoxic_Lupus_Webb<-readRDS("Lupus_Webb_cytotoxic_QC.RDS")
cytotoxic_Lupus_Flynn<-readRDS("Lupus_Flynn_DN2_cytotoxic_QC.RDS") #same as cytotoxic
cytotoxic_MS_Friese<-readRDS("MS_Friese_Subset4_cytotoxic_QC.RDS")
cytotoxic_Crohns_Martin<-readRDS("CD_Martin_cytotoxic_QC.RDS")
#cytotoxic_Covid_Youngs_Wilk<-readRDS("Covid_Youngs_Wilk_cytotoxic_QC.RDS")
cytotoxic_Malaria_Holla<-readRDS("Malaria_Holla_cytotoxic_QC.RDS") 
cytotoxic_HIV_Holla<-readRDS("HIV_Holla_cytotoxic_QC.RDS")

DefaultAssay(DN2_Healthy) <- "RNA"
DefaultAssay(DN2_Lupus_Webb) <- "RNA"
DefaultAssay(DN2_Lupus_Flynn) <- "RNA"
DefaultAssay(DN2_MS_Friese) <- "RNA"
DefaultAssay(DN2_Crohns_Martin) <- "RNA"
DefaultAssay(DN2_Covid_Youngs_Wilk) <- "RNA"
DefaultAssay(DN2_Malaria_Holla) <- "RNA"
DefaultAssay(DN2_HIV_Holla) <- "RNA"

DN2_Healthy[['SCT']]<-NULL
DN2_Lupus_Webb[['SCT']]<-NULL
DN2_Lupus_Flynn[['SCT']]<-NULL
DN2_MS_Friese[['SCT']]<-NULL
DN2_Crohns_Martin[['SCT']]<-NULL
DN2_Covid_Youngs_Wilk[['SCT']]<-NULL
DN2_Malaria_Holla[['SCT']]<-NULL
DN2_HIV_Holla[['SCT']]<-NULL


DefaultAssay(cytotoxic_Healthy) <- "RNA"
DefaultAssay(cytotoxic_Lupus_Webb) <- "RNA"
DefaultAssay(cytotoxic_Lupus_Flynn) <- "RNA"
DefaultAssay(cytotoxic_MS_Friese) <- "RNA"
DefaultAssay(cytotoxic_Crohns_Martin) <- "RNA"
#DefaultAssay(cytotoxic_Covid_Youngs_Wilk) <- "RNA"
DefaultAssay(cytotoxic_Malaria_Holla) <- "RNA"
DefaultAssay(cytotoxic_HIV_Holla) <- "RNA"

cytotoxic_Healthy[['SCT']]<-NULL
cytotoxic_Lupus_Webb[['SCT']]<-NULL
cytotoxic_Lupus_Flynn[['SCT']]<-NULL
cytotoxic_MS_Friese[['SCT']]<-NULL
cytotoxic_Crohns_Martin[['SCT']]<-NULL
#cytotoxic_Covid_Youngs_Wilk[['SCT']]<-NULL
cytotoxic_Malaria_Holla[['SCT']]<-NULL
cytotoxic_HIV_Holla[['SCT']]<-NULL


#add study id
DN2_MS_Friese@meta.data$Atlas_Study_ID<-'Friese_MS'
DN2_Malaria_Holla@meta.data$Atlas_Study_ID<-'Holla_Malaria'
DN2_HIV_Holla@meta.data$Atlas_Study_ID<-'Holla_HIV'
DN2_Covid_Youngs_Wilk@meta.data$Atlas_Study_ID<-'Youngs_Wilk_Covid'

#add study id
cytotoxic_Healthy@meta.data$Atlas_Study_ID<-'Friese_MS'
cytotoxic_Malaria_Holla@meta.data$Atlas_Study_ID<-'Holla_Malaria'
cytotoxic_HIV_Holla@meta.data$Atlas_Study_ID<-'Holla_HIV'
#cytotoxic_Covid_Youngs_Wilk@meta.data$Atlas_Study_ID<-'Youngs_Wilk_Covid'


merged_DN2 <- merge(DN2_Healthy, y=c(DN2_Lupus_Webb,DN2_Lupus_Flynn,DN2_MS_Friese,DN2_Crohns_Martin,DN2_Covid_Youngs_Wilk,
                    DN2_Malaria_Holla,DN2_HIV_Holla))       

merged_cytotoxic <- merge(cytotoxic_Healthy, y=c(cytotoxic_Lupus_Webb,cytotoxic_Lupus_Flynn,cytotoxic_MS_Friese,cytotoxic_Crohns_Martin,
                          cytotoxic_Malaria_Holla,cytotoxic_HIV_Holla))

merged_all<-merge(DN2_Healthy, y=c(DN2_Lupus_Webb,DN2_Lupus_Flynn,DN2_MS_Friese,DN2_Crohns_Martin,DN2_Covid_Youngs_Wilk,
                                   DN2_Malaria_Holla,DN2_HIV_Holla,cytotoxic_Healthy,cytotoxic_Lupus_Webb,cytotoxic_MS_Friese,
                                   cytotoxic_Crohns_Martin,cytotoxic_Malaria_Holla,cytotoxic_HIV_Holla,
                                   #cytotoxic_Covid_Youngs_Wilk
))

merged_Healthy<-merge(DN2_Healthy, y=c(cytotoxic_Healthy))  


unique(merged_cytotoxic$Atlas_Study_ID)



#RPCA PIPELINE

# split the dataset into a list of two seurat objects (stim and CTRL)
merged_cytotoxic.list <- SplitObject(merged_cytotoxic, split.by = "Atlas_Study_ID")

# normalize and identify variable features for each dataset independently
merged_cytotoxic.list <- lapply(X = merged_cytotoxic.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = merged_cytotoxic.list)
merged_cytotoxic.list <- lapply(X = merged_cytotoxic.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})


#ALTERED METHOD
npcs <- 28  # Set the number of principal components to 10, adjust as needed

merged_cytotoxic.list <- lapply(X = merged_cytotoxic.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, npcs = npcs, features = features, verbose = FALSE)  # Specify npcs here
})



immune.anchors <- FindIntegrationAnchors(object.list = merged_cytotoxic.list, anchor.features = features, reduction = "rpca", dims=1:28)


# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors, k.weight = 66) #k.weight must be lower than 67

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.7)

# Visualization
DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "Atlas_Study_ID")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

DefaultAssay(immune.combined) <- "RNA"

immune.combined<-NormalizeData(immune.combined)

saveRDS(immune.combined, 'merged_DN2_RPCA_Healthy&Disease.rds')
immune.combined <-readRDS('merged_DN2_RPCA_Healthy.rds')




DimPlot(object = immune.combined, reduction = "umap", shuffle = F, group.by="Atlas_Sex", 
        cols=c('female'='red',
               'male'='blue'))

unique(immune.combined.clean$Atlas_Age_Category)
DimPlot(object = immune.combined, reduction = "umap", shuffle = F, group.by="Atlas_Age_Category", 
        cols=c('YA'='green',
               'MA'='blue',
               'OA'='red'))

unique(immune.combined$Atlas_Ethnicity_Broad)
DimPlot(object = immune.combined, reduction = "umap", shuffle = F, group.by="Atlas_Ethnicity_Broad", 
        cols=c('EUR'='blue',
               'ASI'='green',
               'AFR'='orange',
               'Ashkenazi_Jewish'='red'))

unique(immune.combined$Atlas_Ethnicity_Fine)
DimPlot(object = immune.combined, reduction = "umap", shuffle = F, group.by="Atlas_Ethnicity_Fine", 
        cols=c('EUR'='cyan',
               'ASI'='green',
               'AFR'='orange',
               'Malaysian'='pink',
               'Indian'='purple',
               'Singaporean Chinese'='yellow',
               'European'='cyan',
               'Japanese'='blue',
               'Korean'='black',
               'Ashkenazi_Jewish'='red'))

unique(immune.combined$Atlas_Study_ID)
DimPlot(object = immune.combined, reduction = "umap", shuffle = F, group.by="Atlas_Study_ID", 
        cols=c('Jimmie_YE_2022'='red',
               'Powell_2022'='cyan',
               'AIDA'='green',
               'Barreiro'='pink',
               'Webb_Lupus'='purple',
               'Flynn_Lupus'='yellow',
               'Friese_MS'='blue',
               'Martin_CD'='black',
               'Youngs_Wilk_Covid'='tan',
               'Holla_Malaria'='darkgreen',
               'Holla_HIV'='orange'))

unique(immune.combined$Atlas_Disease)
DimPlot(object = immune.combined, reduction = "umap", shuffle = F, group.by="Atlas_Disease", 
        cols=c('Healthy'='snow',
               'SLE Inactive'='snow',
               'SLE Active'='snow',
               'MS'='snow',
               'CD'='snow',
               'Severe_COVID'='snow',
               'Malaria'='snow',
               'HIV'='snow'))

unique(immune.combined$seurat_clusters)
DimPlot(object = immune.combined, reduction = "umap", shuffle = F, group.by="seurat_clusters", 
        cols=c('0'='snow',
               '1'='blue',
               '2'='snow', #?
               '3'='green',
               '4'='purple',
               '5'='pink', 
               '6'='black', #?
               '7'='orange',
               '8'='brown',
               '9'='green',
               '10'='red')) #?


#DN2 Prop by age
dittoBarPlot(
  object = immune.combined,
  var = "Atlas_Disease",    #colors stacked
  group.by = "seurat_clusters",     #x axis
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))
