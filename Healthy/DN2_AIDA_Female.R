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
library(UCell)
library(Nebulosa)


getwd()
setwd("/Users/honokamiura/Downloads/Research/DN2")

DN2_AIDA_Female <-readRDS('/Users/honokamiura/Downloads/Research/DN2/AIDA_female_B.rds')

#QC %MT^
DN2_AIDA_Female[["percent.mt"]] <- PercentageFeatureSet(DN2_AIDA_Female, pattern = "^MT.")

#violin plot to look at percent.mt, nFeature_RNA
VlnPlot(DN2_AIDA_Female, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
DN2_AIDA_Female <- subset(DN2_AIDA_Female, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
DN2_AIDA_Female <- NormalizeData(DN2_AIDA_Female)

#identify highly variably features (genes) total of 2000
DN2_AIDA_Female <- FindVariableFeatures(DN2_AIDA_Female, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
all.genes <- rownames(DN2_AIDA_Female)
DN2_AIDA_Female <- ScaleData(DN2_AIDA_Female, features = all.genes)

# Perform PCA and color by cell cycle phase
DN2_AIDA_Female <- RunPCA(DN2_AIDA_Female)

#Cluster the cells
DN2_AIDA_Female <- FindNeighbors(DN2_AIDA_Female, dims = 1:10)
DN2_AIDA_Female <- FindClusters(DN2_AIDA_Female, resolution = 0.5)

#Run UMAP
DN2_AIDA_Female <- RunUMAP(DN2_AIDA_Female, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(DN2_AIDA_Female, reduction = "umap")

saveRDS(DN2_AIDA_Female, "DN2_AIDA_Female_QC.RDS")

#Start here: 
DN2_AIDA_Female<-readRDS('DN2_AIDA_Female_QC.RDS')

Idents(DN2_AIDA_Female)<-'seurat_clusters'

dittoBarPlot(
  object = DN2_AIDA_Female,
  var = "sample_uuid", 
  group.by = "seurat_clusters", 
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))&NoLegend()



#Looking for CD11c+(ITGAX), Tbet+ (TBX21), CD21-(CR2), Vreb3-, Ltb-,
#CD32b=FCGR2B, HLA-DR=HLA-DRA and HLA-DRB, IGD=IGHD, CD62L=SELL
DotPlot(object = DN2_AIDA_Female, #group.by = 'Clusters',
        features = c(
          'CD19','ITGAX','TBX21','ZEB2', 'TLR7','FCGR2B','CD69','HLA-DRA', 'HLA-DRB','CD86','FCRL5','CD22',
          'VPREB3','LTB','CXCR5','CR2','FCRL4','CD27','IGHD','TRAF5','CD24','CD38','SELL'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Look for only Significant DN2 markers
DotPlot(object = DN2_AIDA_Female, #group.by = 'Clusters',
        features = c('ITGAX','TBX21','CR2','VPREB3','LTB'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Looking for Tbet+,
DotPlot(object = DN2_AIDA_Female, #group.by = 'Clusters',
        features = c('TBX21'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Look for cytotoxic markers
DotPlot(object = DN2_AIDA_Female, #group.by = 'Clusters',
        features = c('GZMB', 'GZMA','GZMM','PRF1','NKG7','KLRK1','KLRD1','GNLY'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))


DimPlot(object = DN2_AIDA_Female, reduction = "umap", shuffle = T,
        cols= c('0'= 'snow2',
                '1'='snow2',
                '2'='snow2',
                '3'='snow2',
                '4'='snow2',
                '5'='snow2', 
                '6'='red',
                '7'='snow2', 
                '8'='snow2',
                '9'='snow2'))



# Subset Cluster 6
DN2_AIDA_Female_Subset6<-subset(x = DN2_AIDA_Female, idents ='6')

#QC %MT^
DN2_AIDA_Female_Subset6[["percent.mt"]] <- PercentageFeatureSet(DN2_AIDA_Female_Subset6, pattern = "^MT-")

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
DN2_AIDA_Female_Subset6 <- subset(DN2_AIDA_Female_Subset6, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
DN2_AIDA_Female_Subset6 <- NormalizeData(DN2_AIDA_Female_Subset6)

#identify highly variably features (genes) total of 2000
DN2_AIDA_Female_Subset6 <- FindVariableFeatures(DN2_AIDA_Female_Subset6, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
all.genes.DN2 <- rownames(DN2_AIDA_Female_Subset6)
DN2_AIDA_Female_Subset6 <- ScaleData(DN2_AIDA_Female_Subset6, features = all.genes.DN2)

# Perform PCA 
DN2_AIDA_Female_Subset6 <- RunPCA(DN2_AIDA_Female_Subset6)

#Cluster the cells
DN2_AIDA_Female_Subset6 <- FindNeighbors(DN2_AIDA_Female_Subset6, dims = 1:10)
DN2_AIDA_Female_Subset6 <- FindClusters(DN2_AIDA_Female_Subset6, resolution = 0.1)

#Run UMAP
DN2_AIDA_Female_Subset6 <- RunUMAP(DN2_AIDA_Female_Subset6, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(DN2_AIDA_Female_Subset6, reduction = "umap")

saveRDS(DN2_AIDA_Female_Subset6, "DN2_AIDA_Female_Subset6_QC.RDS")

DN2_AIDA_Female_Subset6<-readRDS('DN2_AIDA_Female_Subset6_QC.RDS')


#Rename Cluster
Idents(DN2_AIDA_Female_Subset6) <- "seurat_clusters"
current.cluster.ids <- levels(DN2_AIDA_Female_Subset6)
current.cluster.ids
new.cluster.ids<-c("DN2.AIDA.Female.A", "DN2.AIDA.Female.B")

DN2_AIDA_Female_Subset6@meta.data$DN2_Subtypes<-DN2_AIDA_Female_Subset6@meta.data$seurat_clusters
DN2_AIDA_Female_Subset6@meta.data$DN2_Subtypes<-plyr::mapvalues(x=DN2_AIDA_Female_Subset6@meta.data$DN2_Subtypes, from=current.cluster.ids, to=new.cluster.ids)

DimPlot(DN2_AIDA_Female_Subset6, group.by="DN2_Subtypes")

saveRDS(DN2_AIDA_Female_Subset6, "DN2_AIDA_Female_Subtypes.RDS")

DN2_AIDA_Female_Subset6<-readRDS('DN2_AIDA_Female_Subtypes.RDS')



#add study id
DN2_AIDA_Female@meta.data$Atlas_Study_ID<-'AIDA'
META<-DN2_AIDA_Female[[]]

DN2_AIDA_Female@meta.data$Atlas_Disease<-'Healthy'

#see what column labels are stored for individual info
unique(DN2_AIDA_Female$sex)
Idents(DN2_AIDA_Female)<-'sex'
current.cluster.ids<-levels(DN2_AIDA_Female@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_AIDA_Female@active.ident)

DN2_AIDA_Female@meta.data$Atlas_Sex<-DN2_AIDA_Female@meta.data$sex
DN2_AIDA_Female@meta.data$Atlas_Sex<-plyr::mapvalues(x=DN2_AIDA_Female@meta.data$Atlas_Sex, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_AIDA_Female$Atlas_Sex)

#Ethnicity- Broad
unique(DN2_AIDA_Female$self_reported_ethnicity)
Idents(DN2_AIDA_Female)<-'self_reported_ethnicity'
current.cluster.ids<-levels(DN2_AIDA_Female@active.ident)
current.cluster.ids<-c("Malaysian" ,          "Indian"        ,      "Singaporean Chinese", "European"  ,         
                         "Japanese"    ,        "Korean" )
new.cluster.ids<-c("ASI" ,          "ASI"        ,      "ASI", "EUR"  ,         
                   "ASI"    ,        "ASI" )

DN2_AIDA_Female@meta.data$Atlas_Ethnicity_Broad<-DN2_AIDA_Female@meta.data$self_reported_ethnicity
DN2_AIDA_Female@meta.data$Atlas_Ethnicity_Broad<-plyr::mapvalues(x=DN2_AIDA_Female@meta.data$Atlas_Ethnicity_Broad, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_AIDA_Female$Atlas_Ethnicity_Broad)

#Ethnicity- Fine
unique(DN2_AIDA_Female$self_reported_ethnicity)
Idents(DN2_AIDA_Female)<-'self_reported_ethnicity'
current.cluster.ids<-levels(DN2_AIDA_Female@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_AIDA_Female@active.ident)

DN2_AIDA_Female@meta.data$Atlas_Ethnicity_Fine<-DN2_AIDA_Female@meta.data$self_reported_ethnicity
DN2_AIDA_Female@meta.data$Atlas_Ethnicity_Fine<-plyr::mapvalues(x=DN2_AIDA_Female@meta.data$Atlas_Ethnicity_Fine, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_AIDA_Female$Atlas_Ethnicity_Fine)

#Atlas Donor
unique(DN2_AIDA_Female$donor_id)
Idents(DN2_AIDA_Female)<-'donor_id'
current.cluster.ids<-levels(DN2_AIDA_Female@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_AIDA_Female@active.ident)

DN2_AIDA_Female@meta.data$Atlas_Donor<-DN2_AIDA_Female@meta.data$'donor_id'
DN2_AIDA_Female@meta.data$Atlas_Donor<-plyr::mapvalues(x=DN2_AIDA_Female@meta.data$Atlas_Donor, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_AIDA_Female$Atlas_Donor)

#Age Updated
META<-DN2_AIDA_Female[[]]
write.csv(META,'AIDA_Female_cell_id_meta.csv')
AIDA_Female_cell_id_meta<-read.csv('AIDA_Female_cell_id_meta.csv', row.names = NULL)

donor_id_meta<-(META%>%distinct(Atlas_Donor, .keep_all=TRUE))
write.csv(donor_id_meta,'AIDA_Female_donor_id_meta.csv')
donor_id_meta<-read.csv('AIDA_Female_donor_id_meta.csv')

merged_df<- merge(donor_id_meta, AIDA_Female_cell_id_meta, by= 'Atlas_Donor', all=TRUE)
#check to see which column has the cell ids/barcodes
write.csv(merged_df,'merged_df.csv')
merged_df2<-read.csv('merged_df.csv', row.names = 5) #row.names = column with cell barcodes/ids
drop <- c("X.1")
merged_df2 = merged_df2[,!(names(merged_df2) %in% drop)]


DN2_AIDA_Female<-AddMetaData(DN2_AIDA_Female,merged_df2 )
unique(DN2_AIDA_Female$Atlas_Age_Category)
saveRDS(DN2_AIDA_Female, "DN2_AIDA_Female_updated_meta.rds")
DN2_AIDA_Female<-readRDS('DN2_AIDA_Female_updated_meta.RDS')




#add study id
DN2_AIDA_Female_Subset6@meta.data$Atlas_Study_ID<-'AIDA'

DN2_AIDA_Female_Subset6@meta.data$Atlas_Disease<-'Healthy'

#see what column labels are stored for individual info
unique(DN2_AIDA_Female_Subset6$sex)
Idents(DN2_AIDA_Female_Subset6)<-'sex'
current.cluster.ids<-levels(DN2_AIDA_Female_Subset6@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_AIDA_Female_Subset6@active.ident)

DN2_AIDA_Female_Subset6@meta.data$Atlas_Sex<-DN2_AIDA_Female_Subset6@meta.data$sex
DN2_AIDA_Female_Subset6@meta.data$Atlas_Sex<-plyr::mapvalues(x=DN2_AIDA_Female_Subset6@meta.data$Atlas_Sex, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_AIDA_Female_Subset6$Atlas_Sex)

#Ethnicity- Broad
unique(DN2_AIDA_Female_Subset6$self_reported_ethnicity)
Idents(DN2_AIDA_Female_Subset6)<-'self_reported_ethnicity'
current.cluster.ids<-levels(DN2_AIDA_Female_Subset6@active.ident)
current.cluster.ids<-c("Malaysian" ,          "Indian"        ,      "Singaporean Chinese", "European"  ,         
                       "Japanese"    ,        "Korean" )
new.cluster.ids<-c("ASI" ,          "ASI"        ,      "ASI", "EUR"  ,         
                   "ASI"    ,        "ASI" )

DN2_AIDA_Female_Subset6@meta.data$Atlas_Ethnicity_Broad<-DN2_AIDA_Female_Subset6@meta.data$self_reported_ethnicity
DN2_AIDA_Female_Subset6@meta.data$Atlas_Ethnicity_Broad<-plyr::mapvalues(x=DN2_AIDA_Female_Subset6@meta.data$Atlas_Ethnicity_Broad, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_AIDA_Female_Subset6$Atlas_Ethnicity_Broad)

#Ethnicity- Fine
unique(DN2_AIDA_Female_Subset6$self_reported_ethnicity)
Idents(DN2_AIDA_Female_Subset6)<-'self_reported_ethnicity'
current.cluster.ids<-levels(DN2_AIDA_Female_Subset6@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_AIDA_Female_Subset6@active.ident)

DN2_AIDA_Female_Subset6@meta.data$Atlas_Ethnicity_Fine<-DN2_AIDA_Female_Subset6@meta.data$self_reported_ethnicity
DN2_AIDA_Female_Subset6@meta.data$Atlas_Ethnicity_Fine<-plyr::mapvalues(x=DN2_AIDA_Female_Subset6@meta.data$Atlas_Ethnicity_Fine, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_AIDA_Female_Subset6$Atlas_Ethnicity_Fine)

#Atlas Donor
unique(DN2_AIDA_Female_Subset6$donor_id)
Idents(DN2_AIDA_Female_Subset6)<-'donor_id'
current.cluster.ids<-levels(DN2_AIDA_Female_Subset6@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_AIDA_Female_Subset6@active.ident)

DN2_AIDA_Female_Subset6@meta.data$Atlas_Donor<-DN2_AIDA_Female_Subset6@meta.data$'donor_id'
DN2_AIDA_Female_Subset6@meta.data$Atlas_Donor<-plyr::mapvalues(x=DN2_AIDA_Female_Subset6@meta.data$Atlas_Donor, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_AIDA_Female_Subset6$Atlas_Donor)

#Age Updated
META<-DN2_AIDA_Female_Subset6[[]]
write.csv(META,'AIDA_Female_Subset6_cell_id_meta.csv')
AIDA_Female_Subset6_cell_id_meta<-read.csv('AIDA_Female_Subset6_cell_id_meta.csv', row.names = NULL)

donor_id_meta<-(META%>%distinct(Atlas_Donor, .keep_all=TRUE))
write.csv(donor_id_meta,'AIDA_Female_Subset6_donor_id_meta.csv')
donor_id_meta<-read.csv('AIDA_Female_Subset6_donor_id_meta.csv')

merged_df<- merge(donor_id_meta, AIDA_Female_Subset6_cell_id_meta, by= 'Atlas_Donor', all=TRUE)
#check to see which column has the cell ids/barcodes
write.csv(merged_df,'merged_df.csv')
merged_df2<-read.csv('merged_df.csv', row.names = 5) #row.names = column with cell barcodes/ids
drop <- c("X.1")
merged_df2 = merged_df2[,!(names(merged_df2) %in% drop)]


DN2_AIDA_Female_Subset6<-AddMetaData(DN2_AIDA_Female_Subset6,merged_df2 )
unique(DN2_AIDA_Female_Subset6$Atlas_Age_Category)
saveRDS(DN2_AIDA_Female_Subset6, "DN2_AIDA_Female_Subset6_updated_meta.rds")
DN2_AIDA_Female_Subset6<-readRDS('DN2_AIDA_Female_Subset6_updated_meta.RDS')
