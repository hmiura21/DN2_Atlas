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
setwd("/Users/honokamiura/Downloads/Olink_Vaccine_Research")

DN2_Powell_Female <-readRDS('/Users/honokamiura/Downloads/Research/DN2/Powell_female_atlas_Bcells.rds')
head(rownames(DN2_Powell_Female))


#QC %MT^
DN2_Powell_Female[["percent.mt"]] <- PercentageFeatureSet(DN2_Powell_Female, pattern = "^MT.")

#violin plot to look at percent.mt, nFeature_RNA
VlnPlot(DN2_Powell_Female, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
DN2_Powell_Female <- subset(DN2_Powell_Female, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
DN2_Powell_Female <- NormalizeData(DN2_Powell_Female)

#identify highly variably features (genes) total of 2000
DN2_Powell_Female <- FindVariableFeatures(DN2_Powell_Female, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
all.genes <- rownames(DN2_Powell_Female)
DN2_Powell_Female <- ScaleData(DN2_Powell_Female, features = all.genes)

# Perform PCA and color by cell cycle phase
DN2_Powell_Female <- RunPCA(DN2_Powell_Female)

#Cluster the cells
DN2_Powell_Female <- FindNeighbors(DN2_Powell_Female, dims = 1:10)
DN2_Powell_Female <- FindClusters(DN2_Powell_Female, resolution = 0.5)

#Run UMAP
DN2_Powell_Female <- RunUMAP(DN2_Powell_Female, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(DN2_Powell_Female, reduction = "umap")

saveRDS(DN2_Powell_Female, "DN2_Female_QC.RDS")

#Start here: 
DN2_Powell_Female<-readRDS('DN2_Female_QC.RDS')

Idents(DN2_Powell_Female)<-'seurat_clusters'

#closer look at female 702_703
dittoBarPlot(
  object = DN2_Powell_Female,
  var = "individual", 
  group.by = "seurat_clusters", 
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))&NoLegend()

DimPlot(object = DN2_Powell_Female, reduction = "umap", shuffle = F, group.by="individual", order=T, pt.size = .5,
        cols=c('702_703'='red'))


DimPlot(object = DN2_Powell_Female, reduction = "umap", shuffle = F, group.by="individual", order=T, pt.size = .5,
        cols=c('685_686'='red'))

META<-DN2_Powell_Female[[]]

#compare one group (0) to another
DN2_Powell_Female_Cluster0_Markers <- FindMarkers(DN2_Powell_Female, only.pos = F,
                                                  ident.1 = ('0'),
                                                  ident.2 = c('1','2','3','4','5','6','7','8','9','10'),
                                                  min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_Powell_Female_Cluster0_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_Powell_Female_Cluster0_Markers.csv")

#compare one group (1) to another
DN2_Powell_Female_Cluster1_Markers <- FindMarkers(DN2_Powell_Female, only.pos = F,
                                                  ident.1 = ('1'),
                                                  ident.2 = c('0','2','3','4','5','6','7','8','9','10'),
                                                  min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_Powell_Female_Cluster1_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_Powell_Female_Cluster1_Markers.csv")

#compare one group (2) to another
DN2_Powell_Female_Cluster2_Markers <- FindMarkers(DN2_Powell_Female, only.pos = F,
                                                  ident.1 = ('2'),
                                                  ident.2 = c('0','1','3','4','5','6','7','8','9','10'),
                                                  min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_Powell_Female_Cluster2_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_Powell_Female_Cluster2_Markers.csv")

#compare one group (3) to another
DN2_Powell_Female_Cluster3_Markers <- FindMarkers(DN2_Powell_Female, only.pos = F,
                                                  ident.1 = ('3'),
                                                  ident.2 = c('0','1','2','4','5','6','7','8','9','10'),
                                                  min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_Powell_Female_Cluster3_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_Powell_Female_Cluster3_Markers.csv")

#compare one group (4) to another
DN2_Powell_Female_Cluster4_Markers <- FindMarkers(DN2_Powell_Female, only.pos = F,
                                                  ident.1 = ('4'),
                                                  ident.2 = c('0','1','2','3','5','6','7','8','9','10'),
                                                  min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_Powell_Female_Cluster4_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_Powell_Female_Cluster4_Markers.csv")


#compare one group (5) to another
DN2_Powell_Female_Cluster5_Markers <- FindMarkers(DN2_Powell_Female, only.pos = F,
                                                  ident.1 = ('5'),
                                                  ident.2 = c('0','1','2','3','4','6','7','8','9','10'),
                                                  min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_Powell_Female_Cluster5_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_Powell_Female_Cluster5_Markers.csv")


#compare one group (6) to another
DN2_Powell_Female_Cluster6_Markers <- FindMarkers(DN2_Powell_Female, only.pos = F,
                                                  ident.1 = ('6'),
                                                  ident.2 = c('0','1','2','3','4','5','7','8','9','10'),
                                                  min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_Powell_Female_Cluster6_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_Powell_Female_Cluster6_Markers.csv")

#compare one group (7) to another
DN2_Powell_Female_Cluster7_Markers <- FindMarkers(DN2_Powell_Female, only.pos = F,
                                                  ident.1 = ('7'),
                                                  ident.2 = c('0','1','2','3','4','5','6','8','9','10'),
                                                  min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_Powell_Female_Cluster7_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_Powell_Female_Cluster7_Markers.csv")

#compare one group (8) to another
DN2_Powell_Female_Cluster8_Markers <- FindMarkers(DN2_Powell_Female, only.pos = F,
                                                  ident.1 = ('8'),
                                                  ident.2 = c('0','1','2','3','4','5','6','7','9','10'),
                                                  min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_Powell_Female_Cluster8_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_Powell_Female_Cluster8_Markers.csv")

#compare one group (9) to another
DN2_Powell_Female_Cluster9_Markers <- FindMarkers(DN2_Powell_Female, only.pos = F,
                                                ident.1 = ('9'),
                                                ident.2 = c('0','1','2','3','4','5','6','7','8','10'),
                                                min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_Powell_Female_Cluster9_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_Powell_Female_Cluster9_Markers.csv")


#compare one group (10) to another
DN2_Powell_Female_Cluster10_Markers <- FindMarkers(DN2_Powell_Female, only.pos = F,
                                                  ident.1 = ('10'),
                                                  ident.2 = c('0','1','2','3','4','5','6','7','8','9'),
                                                  min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_Powell_Female_Cluster10_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_Powell_Female_Cluster10_Markers.csv")


#Looking for CD11c+(ITGAX), Tbet+ (TBX21), CD21-(CR2), Vreb3-, Ltb-,
#CD32b=FCGR2B, HLA-DR=HLA-DRA and HLA-DRB, IGD=IGHD, CD62L=SELL
DotPlot(object = DN2_Powell_Female, #group.by = 'Clusters',
        features = c(
                     'CD19','ITGAX','TBX21','ZEB2', 'TLR7','FCGR2B','CD69','HLA-DRA', 'HLA-DRB','CD86','FCRL5','CD22',
                     'VPREB3','LTB','CXCR5','CR2','FCRL4','CD27','IGHD','TRAF5','CD24','CD38','SELL'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Look for only Significant DN2 markers
DotPlot(object = DN2_Powell_Female, #group.by = 'Clusters',
        features = c('ITGAX','TBX21','CR2','VPREB3','LTB'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Looking for Tbet+,
DotPlot(object = DN2_Powell_Female, #group.by = 'Clusters',
        features = c('TBX21'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Look for cytotoxic markers (cluster 10)
DotPlot(object = DN2_Powell_Female, #group.by = 'Clusters',
        features = c('GZMB', 'GZMA','GZMM','PRF1','NKG7','KLRK1','KLRD1','GNLY'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))


#Look for OTHER cytotoxic markers (cluster 10)
DotPlot(object = DN2_Powell_Female, #group.by = 'Clusters',
        features = c('TNF', 'FAS', 'TNFSF10', 'LTA', 'GZMH', 'GZMK', 'KLRG1', 'CCL5', 'ITGAL'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))



#Dim Plot showing Clusters
DimPlot(DN2_Powell_Female, reduction = "umap")

DimPlot(object = DN2_Powell_Female, reduction = "umap", shuffle = T, order = c('6'),
        cols= c('0'= 'snow2',
                '1'='snow2',
                '2'='snow2',
                '3'='snow2',
                '4'='snow2',
                '5'='snow2', 
                '6'='red',
                '7'='snow2', 
                '8'='snow2',
                '9'='snow2',
                '10'='snow2'))

DimPlot(object = DN2_Powell_Female, reduction = "umap", shuffle = T, order = c('10'),
        cols= c('0'= 'snow2',
                '1'='snow2',
                '2'='snow2',
                '3'='snow2',
                '4'='snow2',
                '5'='snow2', 
                '6'='snow2',
                '7'='snow2', 
                '8'='snow2',
                '9'='snow2',
                '10'='blue'))

# Markers to Show for DN2: Density Plot
markers <- list()
markers$positive<- c('CD19','ITGAX','TBX21','ZEB2', 'TLR7','FCGR2B','CD69','HLA-DRA','CD86','FCRL5','CD22')
markers$negative<- c('VPREB3','LTB','CXCR5','CR2','FCRL4','CD27','IGHD','TRAF5','CD24','CD38','SELL')
markers$positive_import<- c('ITGAX','TBX21')
markers$negative_import<-c('CR2','VPREB3','LTB')

DN2_Powell_Female <- AddModuleScore_UCell(DN2_Powell_Female,features = markers)

#Marker Positive Score
DN2_Powell_Female$positive_UCell
plot_density(
  DN2_Powell_Female,
  c('positive_UCell'),
  slot = NULL,
  joint = F,
  reduction = NULL,
  dims = c(1, 2),
  method = c("ks", "wkde"),
  adjust = 1,
  size = 1,
  shape = 16,
  combine = TRUE,
  pal = "plasma"
)& NoAxes() & NoLegend()

#Marker Negative Score
DN2_Powell_Female$negative_UCell
plot_density(
  DN2_Powell_Female,
  c('negative_UCell'),
  slot = NULL,
  joint = F,
  reduction = NULL,
  dims = c(1, 2),
  method = c("ks", "wkde"),
  adjust = 1,
  size = 1,
  shape = 16,
  combine = TRUE,
  pal = "plasma"
)& NoAxes() & NoLegend()

#Marker Positive Import Score
DN2_Powell_Female$positive_import_UCell
plot_density(
  DN2_Powell_Female,
  c('positive_import_UCell'),
  slot = NULL,
  joint = F,
  reduction = NULL,
  dims = c(1, 2),
  method = c("ks", "wkde"),
  adjust = 1,
  size = 1,
  shape = 16,
  combine = TRUE,
  pal = "plasma"
)& NoAxes() & NoLegend()

#Marker Negative Import Score
DN2_Powell_Female$negative_import_UCell
plot_density(
  DN2_Powell_Female,
  c('negative_import_UCell'),
  slot = NULL,
  joint = F,
  reduction = NULL,
  dims = c(1, 2),
  method = c("ks", "wkde"),
  adjust = 1,
  size = 1,
  shape = 16,
  combine = TRUE,
  pal = "plasma"
)& NoAxes() & NoLegend()





# Subset Cluster 6 and 10 (DN2)
clusters_to_combine <- c(6, 10)
DN2_Subset_Female<-subset(x = DN2_Powell_Female, idents =clusters_to_combine)

#QC %MT^
DN2_Subset_Female[["percent.mt"]] <- PercentageFeatureSet(DN2_Subset_Female, pattern = "^MT-")

#violin plot to look at percent.mt, nFeature_RNA
VlnPlot(DN2_Powell_Female, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
DN2_Subset_Female <- subset(DN2_Subset_Female, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
DN2_Subset_Female <- NormalizeData(DN2_Subset_Female)

#identify highly variably features (genes) total of 2000
DN2_Subset_Female <- FindVariableFeatures(DN2_Subset_Female, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
all.genes.DN2 <- rownames(DN2_Subset_Female)
DN2_Subset_Female <- ScaleData(DN2_Subset_Female, features = all.genes.DN2)

# Perform PCA 
DN2_Subset_Female <- RunPCA(DN2_Subset_Female)

#Cluster the cells
DN2_Subset_Female <- FindNeighbors(DN2_Subset_Female, dims = 1:10)
DN2_Subset_Female <- FindClusters(DN2_Subset_Female, resolution = 0.2)

#Run UMAP
DN2_Subset_Female <- RunUMAP(DN2_Subset_Female, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(DN2_Subset_Female, reduction = "umap")

saveRDS(DN2_Subset_Female, "DN2_Subset_Female_QC.RDS")

DN2_Subset_Female<-readRDS('DN2_Subset_Female_QC.RDS')

#Markers for all DN2
DotPlot(object = DN2_Subset_Female, #group.by = 'Clusters',
        features = c(
          'CD19','ITGAX','TBX21','ZEB2', 'TLR7','FCGR2B','CD69','HLA-DRA', 'HLA-DRB','CD86','FCRL5','CD22',
          'VPREB3','LTB','CXCR5','CR2','FCRL4','CD27','IGHD','TRAF5','CD24','CD38','SELL'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Look for only Significant DN2 markers
DotPlot(object = DN2_Subset_Female, #group.by = 'Clusters',
        features = c('ITGAX','TBX21','CR2','VPREB3','LTB'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Looking for Tbet+,
DotPlot(object = DN2_Subset_Female, #group.by = 'Clusters',
        features = c('TBX21'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Look for cytotoxic markers (to identify original cluster 10)
DotPlot(object = DN2_Subset_Female, #group.by = 'Clusters',
        features = c('GZMB', 'GZMA','GZMM','PRF1','NKG7','KLRK1','KLRD1','GNLY'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))




# Subset Cluster 0,1,2 (remove 3 from batch effect)
clusters_to_combine <- c(0,1,2)
DN2_Subset_Female_clean<-subset(x = DN2_Subset_Female, idents =clusters_to_combine)

#QC %MT^
DN2_Subset_Female_clean[["percent.mt"]] <- PercentageFeatureSet(DN2_Subset_Female_clean, pattern = "^MT-")

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
DN2_Subset_Female_clean <- subset(DN2_Subset_Female_clean, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
DN2_Subset_Female_clean <- NormalizeData(DN2_Subset_Female_clean)

#identify highly variably features (genes) total of 2000
DN2_Subset_Female_clean <- FindVariableFeatures(DN2_Subset_Female_clean, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
all.genes.DN2 <- rownames(DN2_Subset_Female_clean)
DN2_Subset_Female_clean <- ScaleData(DN2_Subset_Female_clean, features = all.genes.DN2)

# Perform PCA 
DN2_Subset_Female_clean <- RunPCA(DN2_Subset_Female_clean)

#Cluster the cells
DN2_Subset_Female_clean <- FindNeighbors(DN2_Subset_Female_clean, dims = 1:10)
DN2_Subset_Female_clean <- FindClusters(DN2_Subset_Female_clean, resolution = 0.2)

#Run UMAP
DN2_Subset_Female_clean <- RunUMAP(DN2_Subset_Female_clean, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(DN2_Subset_Female_clean, reduction = "umap")

saveRDS(DN2_Subset_Female_clean, "DN2_Subset_Female_clean_QC.RDS")

DN2_Subset_Female_clean<-readRDS('DN2_Subset_Female_clean_QC.RDS')

#Rename Cluster
Idents(DN2_Subset_Female_clean) <- "seurat_clusters"
current.cluster.ids <- levels(DN2_Subset_Female_clean)
current.cluster.ids
new.cluster.ids<-c("DN2.A", "DN2.B", "DN2.Cytotoxic")

DN2_Subset_Female_clean@meta.data$DN2_Subtypes<-DN2_Subset_Female_clean@meta.data$seurat_clusters
DN2_Subset_Female_clean@meta.data$DN2_Subtypes<-plyr::mapvalues(x=DN2_Subset_Female_clean@meta.data$DN2_Subtypes, from=current.cluster.ids, to=new.cluster.ids)

DimPlot(DN2_Subset_Female_clean, group.by="DN2_Subtypes")

saveRDS(DN2_Subset_Female_clean, "DN2_Female_Subtypes.RDS")

DN2_Subset_Female_clean<-readRDS('DN2_Female_Subtypes.RDS')


  
#add study id
DN2_Powell_Female@meta.data$Atlas_Study_ID<-'Powell_2022'
DN2_Powell_Female@meta.data$Atlas_Disease<-'Healthy'
META<-DN2_Powell_Female[[]]

#see what column labels are stored for individual info
unique(DN2_Powell_Female$sex)
Idents(DN2_Powell_Female)<-'sex'
current.cluster.ids<-levels(DN2_Powell_Female@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_Powell_Female@active.ident)

DN2_Powell_Female@meta.data$Atlas_Sex<-DN2_Powell_Female@meta.data$sex
DN2_Powell_Female@meta.data$Atlas_Sex<-plyr::mapvalues(x=DN2_Powell_Female@meta.data$Atlas_Sex, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Powell_Female$Atlas_Sex)

#Ethnicity- Broad
unique(DN2_Powell_Female$ethnicity)
Idents(DN2_Powell_Female)<-'ethnicity'
current.cluster.ids<-levels(DN2_Powell_Female@active.ident)
current.cluster.ids<-c("European" )
new.cluster.ids<-c("EUR" )

DN2_Powell_Female@meta.data$Atlas_Ethnicity_Broad<-DN2_Powell_Female@meta.data$ethnicity
DN2_Powell_Female@meta.data$Atlas_Ethnicity_Broad<-plyr::mapvalues(x=DN2_Powell_Female@meta.data$Atlas_Ethnicity_Broad, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Powell_Female$Atlas_Ethnicity_Broad)

#Ethnicity- Fine
DN2_Powell_Female@meta.data$Atlas_Ethnicity_Fine<-DN2_Powell_Female@meta.data$ethnicity
DN2_Powell_Female@meta.data$Atlas_Ethnicity_Fine<-plyr::mapvalues(x=DN2_Powell_Female@meta.data$Atlas_Ethnicity_Fine, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Powell_Female$Atlas_Ethnicity_Fine)

#Atlas Donor
unique(DN2_Powell_Female$individual)
Idents(DN2_Powell_Female)<-'individual'
current.cluster.ids<-levels(DN2_Powell_Female@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_Powell_Female@active.ident)

DN2_Powell_Female@meta.data$Atlas_Donor<-DN2_Powell_Female@meta.data$'individual'
DN2_Powell_Female@meta.data$Atlas_Donor<-plyr::mapvalues(x=DN2_Powell_Female@meta.data$Atlas_Donor, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Powell_Female$Atlas_Donor)

#Age Updated
META<-DN2_Powell_Female[[]]
write.csv(META,'Powell_Female_cell_id_meta.csv')
Powell_Female_cell_id_meta<-read.csv('Powell_Female_cell_id_meta.csv', row.names = NULL)

donor_id_meta<-(META%>%distinct(Atlas_Donor, .keep_all=TRUE))
write.csv(donor_id_meta,'Powell_Female_donor_id_meta.csv')
donor_id_meta<-read.csv('Powell_Female_donor_id_meta.csv')

merged_df<- merge(donor_id_meta, Powell_Female_cell_id_meta, by= 'Atlas_Donor', all=TRUE)
#check to see which column has the cell ids/barcodes
write.csv(merged_df,'merged_df.csv')
merged_df2<-read.csv('merged_df.csv', row.names = 5) #row.names = column with cell barcodes/ids
drop <- c("X.1")
merged_df2 = merged_df2[,!(names(merged_df2) %in% drop)]


DN2_Powell_Female<-AddMetaData(DN2_Powell_Female,merged_df2 )
unique(DN2_Powell_Female$Atlas_Age_Category)
saveRDS(DN2_Powell_Female, "DN2_Powell_Female_updated_meta.rds")
DN2_Powell_Female<-readRDS('DN2_Powell_Female_updated_meta.RDS')




#add study id
DN2_Subset_Female_clean@meta.data$Atlas_Study_ID<-'Powell_2022'
DN2_Subset_Female_clean@meta.data$Atlas_Disease<-'Healthy'
META<-DN2_Subset_Female_clean[[]]


#see what column labels are stored for individual info
unique(DN2_Subset_Female_clean$sex)
Idents(DN2_Subset_Female_clean)<-'sex'
current.cluster.ids<-levels(DN2_Subset_Female_clean@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_Subset_Female_clean@active.ident)

DN2_Subset_Female_clean@meta.data$Atlas_Sex<-DN2_Subset_Female_clean@meta.data$sex
DN2_Subset_Female_clean@meta.data$Atlas_Sex<-plyr::mapvalues(x=DN2_Subset_Female_clean@meta.data$Atlas_Sex, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Subset_Female_clean$Atlas_Sex)

#Ethnicity- Broad
unique(DN2_Subset_Female_clean$ethnicity)
Idents(DN2_Subset_Female_clean)<-'ethnicity'
current.cluster.ids<-levels(DN2_Subset_Female_clean@active.ident)
current.cluster.ids<-c("European" )
new.cluster.ids<-c("EUR" )

DN2_Subset_Female_clean@meta.data$Atlas_Ethnicity_Broad<-DN2_Subset_Female_clean@meta.data$ethnicity
DN2_Subset_Female_clean@meta.data$Atlas_Ethnicity_Broad<-plyr::mapvalues(x=DN2_Subset_Female_clean@meta.data$Atlas_Ethnicity_Broad, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Subset_Female_clean$Atlas_Ethnicity_Broad)

#Ethnicity- Fine
DN2_Subset_Female_clean@meta.data$Atlas_Ethnicity_Fine<-DN2_Subset_Female_clean@meta.data$ethnicity
DN2_Subset_Female_clean@meta.data$Atlas_Ethnicity_Fine<-plyr::mapvalues(x=DN2_Subset_Female_clean@meta.data$Atlas_Ethnicity_Fine, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Subset_Female_clean$Atlas_Ethnicity_Fine)

#Atlas Donor
unique(DN2_Subset_Female_clean$individual)
Idents(DN2_Subset_Female_clean)<-'individual'
current.cluster.ids<-levels(DN2_Subset_Female_clean@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_Subset_Female_clean@active.ident)

DN2_Subset_Female_clean@meta.data$Atlas_Donor<-DN2_Subset_Female_clean@meta.data$'individual'
DN2_Subset_Female_clean@meta.data$Atlas_Donor<-plyr::mapvalues(x=DN2_Subset_Female_clean@meta.data$Atlas_Donor, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Subset_Female_clean$Atlas_Donor)

#Age Updated
META<-DN2_Subset_Female_clean[[]]
write.csv(META,'Powell_Female_clean_cell_id_meta.csv')
Powell_Female_clean_cell_id_meta<-read.csv('Powell_Female_cell_id_meta.csv', row.names = NULL)

donor_id_meta<-(META%>%distinct(Atlas_Donor, .keep_all=TRUE))
write.csv(donor_id_meta,'Powell_Female_clean_donor_id_meta.csv')
donor_id_meta<-read.csv('Powell_Female_clean_donor_id_meta.csv')

merged_df<- merge(donor_id_meta, Powell_Female_clean_cell_id_meta, by= 'Atlas_Donor', all=TRUE)
#check to see which column has the cell ids/barcodes
write.csv(merged_df,'merged_df.csv')
merged_df2<-read.csv('merged_df.csv', row.names = 5) #row.names = column with cell barcodes/ids
drop <- c("X.1")
merged_df2 = merged_df2[,!(names(merged_df2) %in% drop)]


DN2_Subset_Female_clean<-AddMetaData(DN2_Subset_Female_clean,merged_df2 )
unique(DN2_Subset_Female_clean$Atlas_Age_Category)
saveRDS(DN2_Subset_Female_clean, "DN2_Powell_Female_clean_updated_meta.rds")
DN2_Subset_Female_clean<-readRDS('DN2_Powell_Female_clean_updated_meta.RDS')






# Subset Cluster 6  (DN2)
DN2_Subset_Female<-subset(x = DN2_Powell_Female, idents ='6')


#QC %MT^
DN2_Subset_Female[["percent.mt"]] <- PercentageFeatureSet(DN2_Subset_Female, pattern = "^MT-")

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
DN2_Subset_Female <- subset(DN2_Subset_Female, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
DN2_Subset_Female <- NormalizeData(DN2_Subset_Female)

#identify highly variably features (genes) total of 2000
DN2_Subset_Female <- FindVariableFeatures(DN2_Subset_Female, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
DN2_Subset_Female <- ScaleData(DN2_Subset_Female)

# Perform PCA 
DN2_Subset_Female <- RunPCA(DN2_Subset_Female)

#Cluster the cells
DN2_Subset_Female <- FindNeighbors(DN2_Subset_Female, dims = 1:10)
DN2_Subset_Female <- FindClusters(DN2_Subset_Female, resolution = 0.2)

#Run UMAP
DN2_Subset_Female <- RunUMAP(DN2_Subset_Female, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(DN2_Subset_Female, reduction = "umap")

saveRDS(DN2_Subset_Female, "DN2_Subset_Female_QC.RDS")

DN2_Subset_Female<-readRDS('DN2_Subset_Female_QC.RDS')

#closer look at female 702_703
dittoBarPlot(
  object = DN2_Subset_Female,
  var = "individual", 
  group.by = "seurat_clusters", 
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))&NoLegend()

DimPlot(object = DN2_Subset_Female, reduction = "umap", shuffle = F, group.by="individual", order=T,# pt.size = .5,
        cols=c('702_703'='blue'))



# Subset Cluster 0,1 (remove individual effect)
clusters_to_combine <- c(0,1)
DN2_Subset_Female_noI<-subset(x = DN2_Subset_Female, idents =clusters_to_combine)

#QC %MT^
DN2_Subset_Female_noI[["percent.mt"]] <- PercentageFeatureSet(DN2_Subset_Female_noI, pattern = "^MT-")

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
DN2_Subset_Female_noI <- subset(DN2_Subset_Female_noI, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
DN2_Subset_Female_noI <- NormalizeData(DN2_Subset_Female_noI)

#identify highly variably features (genes) total of 2000
DN2_Subset_Female_noI <- FindVariableFeatures(DN2_Subset_Female_noI, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
DN2_Subset_Female_noI <- ScaleData(DN2_Subset_Female_noI)

# Perform PCA 
DN2_Subset_Female_noI <- RunPCA(DN2_Subset_Female_noI)

#Cluster the cells
DN2_Subset_Female_noI <- FindNeighbors(DN2_Subset_Female_noI, dims = 1:10)
DN2_Subset_Female_noI <- FindClusters(DN2_Subset_Female_noI, resolution = 0.1)

#Run UMAP
DN2_Subset_Female_noI <- RunUMAP(DN2_Subset_Female_noI, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(DN2_Subset_Female_noI, reduction = "umap")

saveRDS(DN2_Subset_Female_noI, "DN2_Subset_Female_noI.RDS")

DN2_Subset_Female_noI<-readRDS('DN2_Subset_Female_noI.RDS')


#Rename Cluster
Idents(DN2_Subset_Female_noI) <- "seurat_clusters"
current.cluster.ids <- levels(DN2_Subset_Female_noI)
current.cluster.ids
new.cluster.ids<-c("DN2.A", "DN2.B")

DN2_Subset_Female_noI@meta.data$DN2_Subtypes<-DN2_Subset_Female_noI@meta.data$seurat_clusters
DN2_Subset_Female_noI@meta.data$DN2_Subtypes<-plyr::mapvalues(x=DN2_Subset_Female_noI@meta.data$DN2_Subtypes, from=current.cluster.ids, to=new.cluster.ids)

DimPlot(DN2_Subset_Female_noI, group.by="DN2_Subtypes")

saveRDS(DN2_Subset_Female_noI, "DN2_Female_noI_Subtypes.RDS")

DN2_Subset_Female_noI<-readRDS('DN2_Female_noI_Subtypes.RDS')





#add study id
DN2_Subset_Female_noI@meta.data$Atlas_Study_ID<-'Powell_2022'
DN2_Subset_Female_noI@meta.data$Atlas_Disease<-'Healthy'
META<-DN2_Subset_Female_noI[[]]


#see what column labels are stored for individual info
unique(DN2_Subset_Female_noI$sex)
Idents(DN2_Subset_Female_noI)<-'sex'
current.cluster.ids<-levels(DN2_Subset_Female_noI@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_Subset_Female_noI@active.ident)

DN2_Subset_Female_noI@meta.data$Atlas_Sex<-DN2_Subset_Female_noI@meta.data$sex
DN2_Subset_Female_noI@meta.data$Atlas_Sex<-plyr::mapvalues(x=DN2_Subset_Female_noI@meta.data$Atlas_Sex, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Subset_Female_clean$Atlas_Sex)

#Ethnicity- Broad
unique(DN2_Subset_Female_clean$ethnicity)
Idents(DN2_Subset_Female_clean)<-'ethnicity'
current.cluster.ids<-levels(DN2_Subset_Female_clean@active.ident)
current.cluster.ids<-c("European" )
new.cluster.ids<-c("EUR" )

DN2_Subset_Female_clean@meta.data$Atlas_Ethnicity_Broad<-DN2_Subset_Female_clean@meta.data$ethnicity
DN2_Subset_Female_clean@meta.data$Atlas_Ethnicity_Broad<-plyr::mapvalues(x=DN2_Subset_Female_clean@meta.data$Atlas_Ethnicity_Broad, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Subset_Female_clean$Atlas_Ethnicity_Broad)

#Ethnicity- Fine
DN2_Subset_Female_clean@meta.data$Atlas_Ethnicity_Fine<-DN2_Subset_Female_clean@meta.data$ethnicity
DN2_Subset_Female_clean@meta.data$Atlas_Ethnicity_Fine<-plyr::mapvalues(x=DN2_Subset_Female_clean@meta.data$Atlas_Ethnicity_Fine, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Subset_Female_clean$Atlas_Ethnicity_Fine)

#Atlas Donor
unique(DN2_Subset_Female_clean$individual)
Idents(DN2_Subset_Female_clean)<-'individual'
current.cluster.ids<-levels(DN2_Subset_Female_clean@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_Subset_Female_clean@active.ident)

DN2_Subset_Female_clean@meta.data$Atlas_Donor<-DN2_Subset_Female_clean@meta.data$'individual'
DN2_Subset_Female_clean@meta.data$Atlas_Donor<-plyr::mapvalues(x=DN2_Subset_Female_clean@meta.data$Atlas_Donor, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Subset_Female_clean$Atlas_Donor)

#Age Updated
META<-DN2_Subset_Female_clean[[]]
write.csv(META,'Powell_Female_clean_cell_id_meta.csv')
Powell_Female_clean_cell_id_meta<-read.csv('Powell_Female_cell_id_meta.csv', row.names = NULL)

donor_id_meta<-(META%>%distinct(Atlas_Donor, .keep_all=TRUE))
write.csv(donor_id_meta,'Powell_Female_clean_donor_id_meta.csv')
donor_id_meta<-read.csv('Powell_Female_clean_donor_id_meta.csv')

merged_df<- merge(donor_id_meta, Powell_Female_clean_cell_id_meta, by= 'Atlas_Donor', all=TRUE)
#check to see which column has the cell ids/barcodes
write.csv(merged_df,'merged_df.csv')
merged_df2<-read.csv('merged_df.csv', row.names = 5) #row.names = column with cell barcodes/ids
drop <- c("X.1")
merged_df2 = merged_df2[,!(names(merged_df2) %in% drop)]


DN2_Subset_Female_clean<-AddMetaData(DN2_Subset_Female_clean,merged_df2 )
unique(DN2_Subset_Female_clean$Atlas_Age_Category)
saveRDS(DN2_Subset_Female_clean, "DN2_Powell_Female_clean_updated_meta.rds")
DN2_Subset_Female_clean<-readRDS('DN2_Powell_Female_clean_updated_meta.RDS')

META<-DN2_Subset_Female_clean[[]]



#closer look at female 702_703
dittoBarPlot(
  object = DN2_Subset_Female,
  var = "individual", 
  group.by = "seurat_clusters", 
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))&NoLegend()

DimPlot(object = DN2_Subset_Female, reduction = "umap", shuffle = F, group.by="individual", order=T, pt.size = .5,
        cols=c('702_703'='red'))


DimPlot(object = DN2_Subset_Female, reduction = "umap", shuffle = F, group.by="individual", order=T, pt.size = .5,
        cols=c('685_686'='red'))

#compare blood cancer patient to other
DN2_Powell_Female_BloodCancer_Markers <- FindMarkers(DN2_Subset_Female, only.pos = F,
                                                  ident.1 = ('3'),
                                                  ident.2 = c('0','1','2'),
                                                  min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_Powell_Female_BloodCancer_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_Powell_Female_BloodCancer_Markers.csv")
