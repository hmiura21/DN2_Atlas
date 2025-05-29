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

DN2_AIDA_Male <-readRDS('/Users/honokamiura/Downloads/Research/DN2/AIDA_Male_B.rds')

#QC %MT^
DN2_AIDA_Male[["percent.mt"]] <- PercentageFeatureSet(DN2_AIDA_Male, pattern = "^MT.")

#violin plot to look at percent.mt, nFeature_RNA
VlnPlot(DN2_AIDA_Male, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
DN2_AIDA_Male <- subset(DN2_AIDA_Male, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
DN2_AIDA_Male <- NormalizeData(DN2_AIDA_Male)

#identify highly variably features (genes) total of 2000
DN2_AIDA_Male <- FindVariableFeatures(DN2_AIDA_Male, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
DN2_AIDA_Male <- ScaleData(DN2_AIDA_Male)

# Perform PCA and color by cell cycle phase
DN2_AIDA_Male <- RunPCA(DN2_AIDA_Male)

#Cluster the cells
DN2_AIDA_Male <- FindNeighbors(DN2_AIDA_Male, dims = 1:10)
DN2_AIDA_Male <- FindClusters(DN2_AIDA_Male, resolution = 0.5)

#Run UMAP
DN2_AIDA_Male <- RunUMAP(DN2_AIDA_Male, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(DN2_AIDA_Male, reduction = "umap")

saveRDS(DN2_AIDA_Male, "DN2_AIDA_Male_QC.RDS")

#Start here: 
DN2_AIDA_Male<-readRDS('DN2_AIDA_Male_QC.RDS')

Idents(DN2_AIDA_Male)<-'seurat_clusters'

dittoBarPlot(
  object = DN2_AIDA_Male,
  var = "sample_uuid", 
  group.by = "seurat_clusters", 
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))&NoLegend()



#Looking for CD11c+(ITGAX), Tbet+ (TBX21), CD21-(CR2), Vreb3-, Ltb-,
#CD32b=FCGR2B, HLA-DR=HLA-DRA and HLA-DRB, IGD=IGHD, CD62L=SELL
DotPlot(object = DN2_AIDA_Male, #group.by = 'Clusters',
        features = c(
          'CD19','ITGAX','TBX21','ZEB2', 'TLR7','FCGR2B','CD69','HLA-DRA', 'HLA-DRB','CD86','FCRL5','CD22',
          'VPREB3','LTB','CXCR5','CR2','FCRL4','CD27','IGHD','TRAF5','CD24','CD38','SELL'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Look for only Significant DN2 markers
DotPlot(object = DN2_AIDA_Male, #group.by = 'Clusters',
        features = c('ITGAX','TBX21','CR2','VPREB3','LTB'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Looking for Tbet+,
DotPlot(object = DN2_AIDA_Male, #group.by = 'Clusters',
        features = c('TBX21'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Look for cytotoxic markers
DotPlot(object = DN2_AIDA_Male, #group.by = 'Clusters',
        features = c('GZMB', 'GZMA','GZMM','PRF1','NKG7','KLRK1','KLRD1','GNLY'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))


DimPlot(object = DN2_AIDA_Male, reduction = "umap", shuffle = T,
        cols= c('0'= 'snow2',
                '1'='snow2',
                '2'='snow2',
                '3'='snow2',
                '4'='snow2',
                '5'='snow2', 
                '6'='snow2',
                '7'='red', 
                '8'='snow2',
                '9'='snow2',
                '10'='snow2'))



# Subset Cluster 7
DN2_AIDA_Male_Subset7<-subset(x = DN2_AIDA_Male, idents ='7')

#QC %MT^
DN2_AIDA_Male_Subset7[["percent.mt"]] <- PercentageFeatureSet(DN2_AIDA_Male_Subset7, pattern = "^MT-")

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
DN2_AIDA_Male_Subset7 <- subset(DN2_AIDA_Male_Subset7, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
DN2_AIDA_Male_Subset7 <- NormalizeData(DN2_AIDA_Male_Subset7)

#identify highly variably features (genes) total of 2000
DN2_AIDA_Male_Subset7 <- FindVariableFeatures(DN2_AIDA_Male_Subset7, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
all.genes.DN2 <- rownames(DN2_AIDA_Male_Subset7)
DN2_AIDA_Male_Subset7 <- ScaleData(DN2_AIDA_Male_Subset7, features = all.genes.DN2)

# Perform PCA 
DN2_AIDA_Male_Subset7 <- RunPCA(DN2_AIDA_Male_Subset7)

#Cluster the cells
DN2_AIDA_Male_Subset7 <- FindNeighbors(DN2_AIDA_Male_Subset7, dims = 1:10)
DN2_AIDA_Male_Subset7 <- FindClusters(DN2_AIDA_Male_Subset7, resolution = 0.2)

#Run UMAP
DN2_AIDA_Male_Subset7 <- RunUMAP(DN2_AIDA_Male_Subset7, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(DN2_AIDA_Male_Subset7, reduction = "umap")

saveRDS(DN2_AIDA_Male_Subset7, "DN2_AIDA_Male_Subset7_QC.RDS")

DN2_AIDA_Male_Subset7<-readRDS('DN2_AIDA_Male_Subset7_QC.RDS')



#Rename Cluster
Idents(DN2_AIDA_Male_Subset7) <- "seurat_clusters"
current.cluster.ids <- levels(DN2_AIDA_Male_Subset7)
current.cluster.ids
new.cluster.ids<-c("DN2.AIDA.Male.A", "DN2.AIDA.Male.B")

DN2_AIDA_Male_Subset7@meta.data$DN2_Subtypes<-DN2_AIDA_Male_Subset7@meta.data$seurat_clusters
DN2_AIDA_Male_Subset7@meta.data$DN2_Subtypes<-plyr::mapvalues(x=DN2_AIDA_Male_Subset7@meta.data$DN2_Subtypes, from=current.cluster.ids, to=new.cluster.ids)

DimPlot(DN2_AIDA_Male_Subset7, group.by="DN2_Subtypes")

saveRDS(DN2_AIDA_Male_Subset7, "DN2_AIDA_Male_Subtypes.RDS")

DN2_AIDA_Male_Subset7<-readRDS('DN2_AIDA_Male_Subtypes.RDS')

unique(DN2_AIDA_Male_Subset7@meta.data$DN2_Subtypes)




#add study id
DN2_AIDA_Male@meta.data$Atlas_Study_ID<-'AIDA'
META<-DN2_AIDA_Male[[]]

DN2_AIDA_Male@meta.data$Atlas_Disease<-'Healthy'

#see what column labels are stored for individual info
unique(DN2_AIDA_Male$sex)
Idents(DN2_AIDA_Male)<-'sex'
current.cluster.ids<-levels(DN2_AIDA_Male@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_AIDA_Male@active.ident)

DN2_AIDA_Male@meta.data$Atlas_Sex<-DN2_AIDA_Male@meta.data$sex
DN2_AIDA_Male@meta.data$Atlas_Sex<-plyr::mapvalues(x=DN2_AIDA_Male@meta.data$Atlas_Sex, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_AIDA_Male$Atlas_Sex)

#Ethnicity- Broad
unique(DN2_AIDA_Male$self_reported_ethnicity)
Idents(DN2_AIDA_Male)<-'self_reported_ethnicity'
current.cluster.ids<-levels(DN2_AIDA_Male@active.ident)
current.cluster.ids<-c("Singaporean Chinese" ,"Malaysian"    ,       "European"    ,        "Indian"       ,       "Japanese" ,          
                       "Korean"  )
new.cluster.ids<-c("ASI" ,"ASI"    ,       "EUR"    ,        "ASI"       ,       "ASI" ,          
                   "ASI"  )

DN2_AIDA_Male@meta.data$Atlas_Ethnicity_Broad<-DN2_AIDA_Male@meta.data$self_reported_ethnicity
DN2_AIDA_Male@meta.data$Atlas_Ethnicity_Broad<-plyr::mapvalues(x=DN2_AIDA_Male@meta.data$Atlas_Ethnicity_Broad, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_AIDA_Male$Atlas_Ethnicity_Broad)

#Ethnicity- Fine
unique(DN2_AIDA_Male$self_reported_ethnicity)
Idents(DN2_AIDA_Male)<-'self_reported_ethnicity'
current.cluster.ids<-levels(DN2_AIDA_Male@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_AIDA_Male@active.ident)

DN2_AIDA_Male@meta.data$Atlas_Ethnicity_Fine<-DN2_AIDA_Male@meta.data$self_reported_ethnicity
DN2_AIDA_Male@meta.data$Atlas_Ethnicity_Fine<-plyr::mapvalues(x=DN2_AIDA_Male@meta.data$Atlas_Ethnicity_Fine, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_AIDA_Male$Atlas_Ethnicity_Fine)

#Atlas Donor
unique(DN2_AIDA_Male$donor_id)
Idents(DN2_AIDA_Male)<-'donor_id'
current.cluster.ids<-levels(DN2_AIDA_Male@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_AIDA_Male@active.ident)

DN2_AIDA_Male@meta.data$Atlas_Donor<-DN2_AIDA_Male@meta.data$'donor_id'
DN2_AIDA_Male@meta.data$Atlas_Donor<-plyr::mapvalues(x=DN2_AIDA_Male@meta.data$Atlas_Donor, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_AIDA_Male$Atlas_Donor)

#Age Updated
META<-DN2_AIDA_Male[[]]
write.csv(META,'AIDA_Male_cell_id_meta.csv')
AIDA_Male_cell_id_meta<-read.csv('AIDA_Male_cell_id_meta.csv', row.names = NULL)

donor_id_meta<-(META%>%distinct(Atlas_Donor, .keep_all=TRUE))
write.csv(donor_id_meta,'AIDA_Male_donor_id_meta.csv')
donor_id_meta<-read.csv('AIDA_Male_donor_id_meta.csv')

merged_df<- merge(donor_id_meta, AIDA_Male_cell_id_meta, by= 'Atlas_Donor', all=TRUE)
#check to see which column has the cell ids/barcodes
write.csv(merged_df,'merged_df.csv')
merged_df2<-read.csv('merged_df.csv', row.names = 5) #row.names = column with cell barcodes/ids
drop <- c("X.1")
merged_df2 = merged_df2[,!(names(merged_df2) %in% drop)]


DN2_AIDA_Male<-AddMetaData(DN2_AIDA_Male,merged_df2 )
unique(DN2_AIDA_Male$Atlas_Age_Category)
saveRDS(DN2_AIDA_Male, "DN2_AIDA_Male_updated_meta.rds")
DN2_AIDA_Male<-readRDS('DN2_AIDA_Male_updated_meta.RDS')





#add study id
DN2_AIDA_Male_Subset7@meta.data$Atlas_Study_ID<-'AIDA'

DN2_AIDA_Male_Subset7@meta.data$Atlas_Disease<-'Healthy'

#see what column labels are stored for individual info
unique(DN2_AIDA_Male_Subset7$sex)
Idents(DN2_AIDA_Male_Subset7)<-'sex'
current.cluster.ids<-levels(DN2_AIDA_Male_Subset7@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_AIDA_Male_Subset7@active.ident)

DN2_AIDA_Male_Subset7@meta.data$Atlas_Sex<-DN2_AIDA_Male_Subset7@meta.data$sex
DN2_AIDA_Male_Subset7@meta.data$Atlas_Sex<-plyr::mapvalues(x=DN2_AIDA_Male_Subset7@meta.data$Atlas_Sex, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_AIDA_Male_Subset7$Atlas_Sex)

#Ethnicity- Broad
unique(DN2_AIDA_Male_Subset7$self_reported_ethnicity)
Idents(DN2_AIDA_Male_Subset7)<-'self_reported_ethnicity'
current.cluster.ids<-levels(DN2_AIDA_Male_Subset7@active.ident)
current.cluster.ids<-c("Singaporean Chinese" ,"Malaysian"    ,       "European"    ,        "Indian"       ,       "Japanese" ,          
                       "Korean"  )
new.cluster.ids<-c("ASI" ,"ASI"    ,       "EUR"    ,        "ASI"       ,       "ASI" ,          
                   "ASI"  )

DN2_AIDA_Male_Subset7@meta.data$Atlas_Ethnicity_Broad<-DN2_AIDA_Male_Subset7@meta.data$self_reported_ethnicity
DN2_AIDA_Male_Subset7@meta.data$Atlas_Ethnicity_Broad<-plyr::mapvalues(x=DN2_AIDA_Male_Subset7@meta.data$Atlas_Ethnicity_Broad, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_AIDA_Male_Subset7$Atlas_Ethnicity_Broad)

#Ethnicity- Fine
unique(DN2_AIDA_Male_Subset7$self_reported_ethnicity)
Idents(DN2_AIDA_Male_Subset7)<-'self_reported_ethnicity'
current.cluster.ids<-levels(DN2_AIDA_Male_Subset7@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_AIDA_Male_Subset7@active.ident)

DN2_AIDA_Male_Subset7@meta.data$Atlas_Ethnicity_Fine<-DN2_AIDA_Male_Subset7@meta.data$self_reported_ethnicity
DN2_AIDA_Male_Subset7@meta.data$Atlas_Ethnicity_Fine<-plyr::mapvalues(x=DN2_AIDA_Male_Subset7@meta.data$Atlas_Ethnicity_Fine, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_AIDA_Male_Subset7$Atlas_Ethnicity_Fine)

#Atlas Donor
unique(DN2_AIDA_Male_Subset7$donor_id)
Idents(DN2_AIDA_Male_Subset7)<-'donor_id'
current.cluster.ids<-levels(DN2_AIDA_Male_Subset7@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_AIDA_Male_Subset7@active.ident)

DN2_AIDA_Male_Subset7@meta.data$Atlas_Donor<-DN2_AIDA_Male_Subset7@meta.data$'donor_id'
DN2_AIDA_Male_Subset7@meta.data$Atlas_Donor<-plyr::mapvalues(x=DN2_AIDA_Male_Subset7@meta.data$Atlas_Donor, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_AIDA_Male_Subset7$Atlas_Donor)

#Age Updated
META<-DN2_AIDA_Male_Subset7[[]]
write.csv(META,'AIDA_Male_Subset7_cell_id_meta.csv')
AIDA_Male_Subset7_cell_id_meta<-read.csv('AIDA_Male_Subset7_cell_id_meta.csv', row.names = NULL)

donor_id_meta<-(META%>%distinct(Atlas_Donor, .keep_all=TRUE))
write.csv(donor_id_meta,'AIDA_Male_Subset7_donor_id_meta.csv')
donor_id_meta<-read.csv('AIDA_Male_Subset7_donor_id_meta.csv')

merged_df<- merge(donor_id_meta, AIDA_Male_Subset7_cell_id_meta, by= 'Atlas_Donor', all=TRUE)
#check to see which column has the cell ids/barcodes
write.csv(merged_df,'merged_df.csv')
merged_df2<-read.csv('merged_df.csv', row.names = 5) #row.names = column with cell barcodes/ids
drop <- c("X.1")
merged_df2 = merged_df2[,!(names(merged_df2) %in% drop)]


DN2_AIDA_Male_Subset7<-AddMetaData(DN2_AIDA_Male_Subset7,merged_df2 )
unique(DN2_AIDA_Male_Subset7$Atlas_Age_Category)
saveRDS(DN2_AIDA_Male_Subset7, "DN2_AIDA_Male_Subset7_updated_meta.rds")
DN2_AIDA_Male_Subset7<-readRDS('DN2_AIDA_Male_Subset7_updated_meta.RDS')
