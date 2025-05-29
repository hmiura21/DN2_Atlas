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

setwd("/Users/honokamiura/Downloads/Research/DN2")

DN2_Barreiro <-readRDS('/Users/honokamiura/Downloads/Research/DN2/Barreiro_B_cells.rds')


#QC %MT^
DN2_Barreiro[["percent.mt"]] <- PercentageFeatureSet(DN2_Barreiro, pattern = "^MT.")

#violin plot to look at percent.mt, nFeature_RNA
VlnPlot(DN2_Barreiro, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
DN2_Barreiro <- subset(DN2_Barreiro, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
DN2_Barreiro <- NormalizeData(DN2_Barreiro)

#identify highly variably features (genes) total of 2000
DN2_Barreiro <- FindVariableFeatures(DN2_Barreiro, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
DN2_Barreiro <- ScaleData(DN2_Barreiro)

# Perform PCA and color by cell cycle phase
DN2_Barreiro <- RunPCA(DN2_Barreiro)

#Cluster the cells
DN2_Barreiro <- FindNeighbors(DN2_Barreiro, dims = 1:10)
DN2_Barreiro <- FindClusters(DN2_Barreiro, resolution = 0.6)

#Run UMAP
DN2_Barreiro <- RunUMAP(DN2_Barreiro, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(DN2_Barreiro, reduction = "umap")

saveRDS(DN2_Barreiro, "DN2_Barreiro_QC.RDS")

#Start here: 
DN2_Barreiro<-readRDS('DN2_Barreiro_QC.RDS')

Idents(DN2_Barreiro)<-'seurat_clusters'

#Looking for CD11c+(ITGAX), Tbet+ (TBX21), CD21-(CR2), Vreb3-, Ltb-,
#CD32b=FCGR2B, HLA-DR=HLA-DRA and HLA-DRB, IGD=IGHD, CD62L=SELL
DotPlot(object = DN2_Barreiro, #group.by = 'Clusters',
        features = c(
          'CD19','ITGAX','TBX21','ZEB2', 'TLR7','FCGR2B','CD69', 'HLA-DRB','CD86','FCRL5','CD22',
          'VPREB3','LTB','CXCR5','CR2','FCRL4','CD27','IGHD','TRAF5','CD24','CD38','SELL'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#dn2 markers
DotPlot(object = DN2_Barreiro, group.by = 'seurat_clusters', 
        features = c("SOX5",'HLA-DRB5','TBX21',"ITGAX","BATF",'JAZF1','NFATC2',
                     'ZEB2','NR4A2','MS4A1','TOX','FGR','IL10RA','FCRL5',
                     'SREBF2','TFEB','TFEC','ZBTB32','MEF2A','POU2F2','SPI1',
                     #NEGATIVE
                     'LTB','VPREB3','SELL','JCHAIN','CD27','IGHD'),
        cols = c('RdBu') ) + RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 7,face = "bold"))


#Look for only Significant DN2 markers
DotPlot(object = DN2_Barreiro, #group.by = 'Clusters',
        features = c('ITGAX','TBX21','CR2','VPREB3','LTB'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Looking for Tbet+,
DotPlot(object = DN2_Barreiro, #group.by = 'Clusters',
        features = c('TBX21'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Look for cytotoxic markers 
DotPlot(object = DN2_Barreiro, #group.by = 'Clusters',
        features = c('GZMB', 'GZMA','GZMM','PRF1','NKG7','KLRK1','KLRD1','GNLY'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))


#violin plot to look at cytotoxicity by cluster
VlnPlot(DN2_Barreiro, features= c('GZMB', 'GZMA','GZMM','PRF1','NKG7','KLRK1','KLRD1','GNLY'), idents=c('0','1','2','3','4','5','6','7','8','9','10','11'))

VlnPlot(DN2_Barreiro, features= 'TBX21', idents=c('0','1','2','3','4','5','6','7','8','9','10','11'))

#individual effect
dittoBarPlot(
  object = DN2_Barreiro,
  var = "orig.ident", 
  group.by = "seurat_clusters", 
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))&NoLegend()


# Subset Cluster 8
DN2_Barreiro<-subset(x = DN2_Barreiro, idents ='8')

#QC %MT^
DN2_Barreiro[["percent.mt"]] <- PercentageFeatureSet(DN2_Barreiro, pattern = "^MT-")

#violin plot to look at percent.mt, nFeature_RNA
VlnPlot(DN2_Barreiro, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
DN2_Barreiro <- subset(DN2_Barreiro, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
DN2_Barreiro <- NormalizeData(DN2_Barreiro)

#identify highly variably features (genes) total of 2000
DN2_Barreiro <- FindVariableFeatures(DN2_Barreiro, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
DN2_Barreiro <- ScaleData(DN2_Barreiro)

# Perform PCA 
DN2_Barreiro <- RunPCA(DN2_Barreiro)

#Cluster the cells
DN2_Barreiro <- FindNeighbors(DN2_Barreiro, dims = 1:10)
DN2_Barreiro <- FindClusters(DN2_Barreiro, resolution = 0.5)

#Run UMAP
DN2_Barreiro <- RunUMAP(DN2_Barreiro, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(DN2_Barreiro_Subset8_QC, reduction = "umap")

saveRDS(DN2_Barreiro, "DN2_Barreiro_Subset8_QC.RDS")

DN2_Barreiro_Subset8_QC<-readRDS('DN2_Barreiro_Subset8_QC.RDS')






#Rename Cluster
Idents(DN2_Barreiro_Subset8_QC) <- "seurat_clusters"
current.cluster.ids <- levels(DN2_Barreiro_Subset8_QC)
current.cluster.ids
new.cluster.ids<-c("DN2.Barriero.A", "DN2.Barriero.B")

DN2_Barreiro_Subset8_QC@meta.data$DN2_Subtypes<-DN2_Barreiro_Subset8_QC@meta.data$seurat_clusters
DN2_Barreiro_Subset8_QC@meta.data$DN2_Subtypes<-plyr::mapvalues(x=DN2_Barreiro_Subset8_QC@meta.data$DN2_Subtypes, from=current.cluster.ids, to=new.cluster.ids)

DimPlot(DN2_Barreiro_Subset8_QC, group.by="DN2_Subtypes")

saveRDS(DN2_Barreiro_Subset8_QC, "DN2_Barreiro_Subtypes.RDS")

DN2_Barreiro_Subset8_QC<-readRDS('DN2_Barreiro_Subtypes.RDS')





#add study id
DN2_Barreiro@meta.data$Atlas_Study_ID<-'Barreiro'
META<-DN2_Barreiro[[]]

DN2_Barreiro@meta.data$Atlas_Disease<-'Healthy'

#see what column labels are stored for individual info
unique(DN2_Barreiro$Atlas_Sex)
Idents(DN2_Barreiro)<-'Atlas_Sex'
current.cluster.ids<-levels(DN2_Barreiro@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_Barreiro@active.ident)

DN2_Barreiro@meta.data$Atlas_Sex<-DN2_Barreiro@meta.data$sex
DN2_Barreiro@meta.data$Atlas_Sex<-plyr::mapvalues(x=DN2_Barreiro@meta.data$Atlas_Sex, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Barreiro$Atlas_Sex)

#Ethnicity- Broad
unique(DN2_Barreiro$Atlas_Ancestry)
Idents(DN2_Barreiro)<-'Atlas_Ancestry'
current.cluster.ids<-levels(DN2_Barreiro@active.ident)
current.cluster.ids<-c( "AFR"   , "EUR")
new.cluster.ids<-c( "AFR"   ,                "EUR"  )

DN2_Barreiro@meta.data$Atlas_Ethnicity_Broad<-DN2_Barreiro@meta.data$Atlas_Ancestry
DN2_Barreiro@meta.data$Atlas_Ethnicity_Broad<-plyr::mapvalues(x=DN2_Barreiro@meta.data$Atlas_Ethnicity_Broad, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Barreiro$Atlas_Ethnicity_Broad)

#Ethnicity- Fine
DN2_Barreiro@meta.data$Atlas_Ethnicity_Fine<-DN2_Barreiro@meta.data$Atlas_Ancestry
DN2_Barreiro@meta.data$Atlas_Ethnicity_Fine<-plyr::mapvalues(x=DN2_Barreiro@meta.data$Atlas_Ethnicity_Fine, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Barreiro$Atlas_Ethnicity_Fine)

#Atlas Donor
unique(DN2_Barreiro$SOC_indiv_ID)
Idents(DN2_Barreiro)<-'SOC_indiv_ID'
current.cluster.ids<-levels(DN2_Barreiro@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_Barreiro@active.ident)

DN2_Barreiro@meta.data$Atlas_Donor<-DN2_Barreiro@meta.data$'SOC_indiv_ID'
DN2_Barreiro@meta.data$Atlas_Donor<-plyr::mapvalues(x=DN2_Barreiro@meta.data$Atlas_Donor, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Barreiro$Atlas_Donor)

#Age (Dont use?)
META<-DN2_Barreiro[[]]
donor_id_meta<-(META%>%distinct(Atlas_Donor, .keep_all=TRUE))

write.csv(donor_id_meta,'Barreiro_donor_id_meta.csv')
donor_id_meta<-read.csv('Barreiro_donor_id_meta.csv')
merged_df<- merge(donor_id_meta, META, by= 'Atlas_Donor', all=TRUE)

DN2_Barreiro<-AddMetaData(DN2_Barreiro,merged_df )
unique(DN2_Barreiro$Atlas_Age_Category)

saveRDS(DN2_Barreiro, "DN2_Barreiro_updatedmeta.RDS")
DN2_Barreiro<-readRDS('DN2_Barreiro_updatedmeta.RDS')

#Age Updated
META<-DN2_Barreiro[[]]
write.csv(META,'Barreiro_cell_id_meta.csv')
Barreiro_cell_id_meta<-read.csv('Barreiro_cell_id_meta.csv', row.names = NULL)

donor_id_meta<-(META%>%distinct(Atlas_Donor, .keep_all=TRUE))
write.csv(donor_id_meta,'Barreiro_donor_id_meta.csv') #correct age here
donor_id_meta<-read.csv('Barreiro_donor_id_meta.csv')

merged_df<- merge(donor_id_meta, Barreiro_cell_id_meta, by= 'Atlas_Donor', all=TRUE)
#check to see which column has the cell ids/barcodes
write.csv(merged_df,'merged_df.csv')
merged_df2<-read.csv('merged_df.csv', row.names = 5) #row.names = column with cell barcodes/ids
drop <- c("X.1")
merged_df2 = merged_df2[,!(names(merged_df2) %in% drop)]


DN2_Barreiro<-AddMetaData(DN2_Barreiro,merged_df2 )
unique(DN2_Barreiro$Atlas_Age_Category)
DimPlot(DN2_Barreiro, group.by = 'seurat_clusters')
saveRDS(DN2_Barreiro, "DN2_Barreiro_updated_meta.rds")
DN2_Barreiro<-readRDS('DN2_Barreiro_updated_meta.RDS')







#add study id
DN2_Barreiro_Subset8_QC@meta.data$Atlas_Study_ID<-'Barreiro'
META<-DN2_Barreiro_Subset8_QC[[]]

DN2_Barreiro_Subset8_QC@meta.data$Atlas_Disease<-'Healthy'

#see what column labels are stored for individual info
unique(DN2_Barreiro_Subset8_QC$Atlas_Sex)
Idents(DN2_Barreiro_Subset8_QC)<-'Atlas_Sex'
current.cluster.ids<-levels(DN2_Barreiro_Subset8_QC@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_Barreiro_Subset8_QC@active.ident)

DN2_Barreiro_Subset8_QC@meta.data$Atlas_Sex<-DN2_Barreiro_Subset8_QC@meta.data$Atlas_Sex
DN2_Barreiro_Subset8_QC@meta.data$Atlas_Sex<-plyr::mapvalues(x=DN2_Barreiro_Subset8_QC@meta.data$Atlas_Sex, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Barreiro_Subset8_QC$Atlas_Sex)

#Ethnicity- Broad
unique(DN2_Barreiro_Subset8_QC$Atlas_Ancestry)
Idents(DN2_Barreiro_Subset8_QC)<-'Atlas_Ancestry'
current.cluster.ids<-levels(DN2_Barreiro_Subset8_QC@active.ident)
current.cluster.ids<-c( "EUR"   , "AFR")
new.cluster.ids<-c( "EUR"   , "AFR")

DN2_Barreiro_Subset8_QC@meta.data$Atlas_Ethnicity_Broad<-DN2_Barreiro_Subset8_QC@meta.data$Atlas_Ancestry
DN2_Barreiro_Subset8_QC@meta.data$Atlas_Ethnicity_Broad<-plyr::mapvalues(x=DN2_Barreiro_Subset8_QC@meta.data$Atlas_Ethnicity_Broad, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Barreiro_Subset8_QC$Atlas_Ethnicity_Broad)

#Ethnicity- Fine
DN2_Barreiro_Subset8_QC@meta.data$Atlas_Ethnicity_Fine<-DN2_Barreiro_Subset8_QC@meta.data$Atlas_Ancestry
DN2_Barreiro_Subset8_QC@meta.data$Atlas_Ethnicity_Fine<-plyr::mapvalues(x=DN2_Barreiro_Subset8_QC@meta.data$Atlas_Ethnicity_Fine, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Barreiro_Subset8_QC$Atlas_Ethnicity_Fine)

#Atlas Donor
unique(DN2_Barreiro_Subset8_QC$SOC_indiv_ID)
Idents(DN2_Barreiro_Subset8_QC)<-'SOC_indiv_ID'
current.cluster.ids<-levels(DN2_Barreiro_Subset8_QC@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_Barreiro_Subset8_QC@active.ident)

DN2_Barreiro_Subset8_QC@meta.data$Atlas_Donor<-DN2_Barreiro_Subset8_QC@meta.data$'SOC_indiv_ID'
DN2_Barreiro_Subset8_QC@meta.data$Atlas_Donor<-plyr::mapvalues(x=DN2_Barreiro_Subset8_QC@meta.data$Atlas_Donor, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Barreiro_Subset8_QC$Atlas_Donor)

#Age Updated
META<-DN2_Barreiro_Subset8_QC[[]]
write.csv(META,'Barreiro_Subset8_cell_id_meta.csv')
Barreiro_Subset8_cell_id_meta<-read.csv('Barreiro_Subset8_cell_id_meta.csv', row.names = NULL)

donor_id_meta<-(META%>%distinct(Atlas_Donor, .keep_all=TRUE))
write.csv(donor_id_meta,'Barreiro_Subset8_donor_id_meta.csv')
donor_id_meta<-read.csv('Barreiro_Subset8_donor_id_meta.csv')

merged_df<- merge(donor_id_meta, Barreiro_Subset8_cell_id_meta, by= 'Atlas_Donor', all=TRUE)
#check to see which column has the cell ids/barcodes
write.csv(merged_df,'merged_df.csv')
merged_df2<-read.csv('merged_df.csv', row.names = 5) #row.names = column with cell barcodes/ids
drop <- c("X.1")
merged_df2 = merged_df2[,!(names(merged_df2) %in% drop)]


DN2_Barreiro_Subset8_QC<-AddMetaData(DN2_Barreiro_Subset8_QC,merged_df2 )
unique(DN2_Barreiro_Subset8_QC$Atlas_Age_Category)
saveRDS(DN2_Barreiro_Subset8_QC, "DN2_Barreiro_Subset8_QC_updated_meta.rds")
DN2_Barreiro_Subset8_QC<-readRDS('DN2_Barreiro_Subset8_QC_updated_meta.RDS')



