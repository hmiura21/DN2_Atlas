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

DN2_Powell_Male <-readRDS('/Users/honokamiura/Downloads/Research/DN2/powell_male_bcells.rds')

#QC %MT^
DN2_Powell_Male[["percent.mt"]] <- PercentageFeatureSet(DN2_Powell_Male, pattern = "^MT.")

#violin plot to look at percent.mt, nFeature_RNA
VlnPlot(DN2_Powell_Male, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
DN2_Powell_Male <- subset(DN2_Powell_Male, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
DN2_Powell_Male <- NormalizeData(DN2_Powell_Male)

#identify highly variably features (genes) total of 2000
DN2_Powell_Male <- FindVariableFeatures(DN2_Powell_Male, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
DN2_Powell_Male <- ScaleData(DN2_Powell_Male)

# Perform PCA and color by cell cycle phase
DN2_Powell_Male <- RunPCA(DN2_Powell_Male)

#Cluster the cells
DN2_Powell_Male <- FindNeighbors(DN2_Powell_Male, dims = 1:10)
DN2_Powell_Male <- FindClusters(DN2_Powell_Male, resolution = 0.8)

#Run UMAP
DN2_Powell_Male <- RunUMAP(DN2_Powell_Male, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(DN2_Powell_Male, reduction = "umap")

saveRDS(DN2_Powell_Male, "DN2_Male_QC.RDS")

#Start here: 
DN2_Powell_Male<-readRDS('DN2_Male_QC.RDS')

Idents(DN2_Powell_Male)<-'seurat_clusters'

#individual effect
dittoBarPlot(
  object = DN2_Powell_Male,
  var = "individual", 
  group.by = "seurat_clusters", 
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))&NoLegend()



#compare one group (0) to another
DN2_Powell_Male_Cluster0_Markers <- FindMarkers(DN2_Powell_Male, only.pos = F,
                                                  ident.1 = ('0'),
                                                  ident.2 = c('1','2','3','4','5','6','7','8'),
                                                  min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_Powell_Male_Cluster0_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_Powell_Male_Cluster0_Markers.csv")

#compare one group (1) to another
DN2_Powell_Male_Cluster1_Markers <- FindMarkers(DN2_Powell_Male, only.pos = F,
                                                  ident.1 = ('1'),
                                                  ident.2 = c('0','2','3','4','5','6','7','8'),
                                                  min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_Powell_Male_Cluster1_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_Powell_Male_Cluster1_Markers.csv")

#compare one group (2) to another
DN2_Powell_Male_Cluster2_Markers <- FindMarkers(DN2_Powell_Male, only.pos = F,
                                                  ident.1 = ('2'),
                                                  ident.2 = c('0','1','3','4','5','6','7','8'),
                                                  min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_Powell_Male_Cluster2_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_Powell_Male_Cluster2_Markers.csv")

#compare one group (3) to another
DN2_Powell_Male_Cluster3_Markers <- FindMarkers(DN2_Powell_Male, only.pos = F,
                                                  ident.1 = ('3'),
                                                  ident.2 = c('0','1','2','4','5','6','7','8'),
                                                  min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_Powell_Male_Cluster3_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_Powell_Male_Cluster3_Markers.csv")

#compare one group (4) to another
DN2_Powell_Male_Cluster4_Markers <- FindMarkers(DN2_Powell_Male, only.pos = F,
                                                  ident.1 = ('4'),
                                                  ident.2 = c('0','1','2','3','5','6','7','8'),
                                                  min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_Powell_Male_Cluster4_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_Powell_Male_Cluster4_Markers.csv")


#compare one group (5) to another
DN2_Powell_Male_Cluster5_Markers <- FindMarkers(DN2_Powell_Male, only.pos = F,
                                                  ident.1 = ('5'),
                                                  ident.2 = c('0','1','2','3','4','6','7','8'),
                                                  min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_Powell_Male_Cluster5_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_Powell_Male_Cluster5_Markers.csv")


#compare one group (6) to another
DN2_Powell_Male_Cluster6_Markers <- FindMarkers(DN2_Powell_Male, only.pos = F,
                                                  ident.1 = ('6'),
                                                  ident.2 = c('0','1','2','3','4','5','7','8'),
                                                  min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_Powell_Male_Cluster6_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_Powell_Male_Cluster6_Markers.csv")

#compare one group (7) to another
DN2_Powell_Male_Cluster7_Markers <- FindMarkers(DN2_Powell_Male, only.pos = F,
                                                  ident.1 = ('7'),
                                                  ident.2 = c('0','1','2','3','4','5','6','8'),
                                                  min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_Powell_Male_Cluster7_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_Powell_Male_Cluster7_Markers.csv")

#compare one group (8) to another
DN2_Powell_Male_Cluster8_Markers <- FindMarkers(DN2_Powell_Male, only.pos = F,
                                                  ident.1 = ('8'),
                                                  ident.2 = c('0','1','2','3','4','5','6','7'),
                                                  min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_Powell_Male_Cluster8_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_Powell_Male_Cluster8_Markers.csv")

#Looking for CD11c+(ITGAX), Tbet+ (TBX21), CD21-(CR2), Vreb3-, Ltb-,
#CD32b=FCGR2B, HLA-DR=HLA-DRA and HLA-DRB, IGD=IGHD, CD62L=SELL
DotPlot(object = DN2_Powell_Male, #group.by = 'Clusters',
        features = c(
          'CD19','ITGAX','TBX21','ZEB2', 'TLR7','FCGR2B','CD69','HLA-DRA', 'HLA-DRB','CD86','FCRL5','CD22',
          'VPREB3','LTB','CXCR5','CR2','FCRL4','CD27','IGHD','TRAF5','CD24','CD38','SELL'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#dn2 markers
DotPlot(object = DN2_Powell_Male, group.by = 'seurat_clusters', 
        features = c("SOX5",'HLA-DRB5','TBX21',"ITGAX","BATF",'JAZF1','NFATC2',
                     'ZEB2','NR4A2','MS4A1','TOX','FGR','IL10RA','FCRL5',
                     'SREBF2','TFEB','TFEC','ZBTB32','MEF2A','POU2F2','SPI1',
                     #NEGATIVE
                     'LTB','VPREB3','SELL','JCHAIN','CD27','IGHD'),
        cols = c('RdBu') ) + RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 7,face = "bold"))


#Look for only Significant DN2 markers
DotPlot(object = DN2_Powell_Male, #group.by = 'Clusters',
        features = c('ITGAX','TBX21','CR2','VPREB3','LTB'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Looking for Tbet+,
DotPlot(object = DN2_Powell_Male, #group.by = 'Clusters',
        features = c('TBX21'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Look for cytotoxic markers 
DotPlot(object = DN2_Powell_Male, #group.by = 'Clusters',
        features = c('GZMB', 'GZMA','GZMM','PRF1','NKG7','KLRK1','KLRD1','GNLY'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Look for OTHER cytotoxic markers 
DotPlot(object = DN2_Powell_Male, #group.by = 'Clusters',
        features = c('TNF', 'FAS', 'TNFSF10', 'LTA', 'GZMH', 'GZMK', 'KLRG1', 'CCL5', 'ITGAL'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Look for granzyme cytotoxic markers 
DotPlot(object = DN2_Powell_Male, #group.by = 'Clusters',
        features = c('GZMB'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))



#Dim Plot showing Clusters
DimPlot(DN2_Powell_Male, reduction = "umap")

DimPlot(object = DN2_Powell_Male, reduction = "umap", shuffle = FALSE, order = 0:11,
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
                '10'='snow2',
                '11'='snow2'
                ))

DimPlot(object = DN2_Powell_Male, reduction = "umap", shuffle = FALSE, order = 0:11,
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
                '10'='snow2',
                '11'='blue'
        ))


# Markers to Show for DN2: Density Plot
markers <- list()
markers$positive<- c('CD19','ITGAX','TBX21','ZEB2', 'TLR7','FCGR2B','CD69','HLA-DRA','CD86','FCRL5','CD22')
markers$negative<- c('VPREB3','LTB','CXCR5','CR2','FCRL4','CD27','IGHD','TRAF5','CD24','CD38','SELL')
markers$positive_import<- c('ITGAX','TBX21')
markers$negative_import<-c('CR2','VPREB3','LTB')

DN2_Powell_Male <- AddModuleScore_UCell(DN2_Powell_Male,features = markers)

#Marker Positive Score
DN2_Powell_Male$positive_UCell
plot_density(
  DN2_COMBINED,
  c('NKG7'),
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
DN2_Powell_Male$negative_UCell
plot_density(
  DN2_Powell_Male,
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
DN2_Powell_Male$positive_import_UCell
plot_density(
  DN2_Powell_Male,
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
DN2_Powell_Male$negative_import_UCell
plot_density(
  DN2_Powell_Male,
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





# Subset Cluster 11 and 7 (DN2)
clusters_to_combine <- c(11, 7)
DN2_Subset_Male<-subset(x = DN2_Powell_Male, idents =clusters_to_combine)


#QC %MT^
DN2_Subset_Male[["percent.mt"]] <- PercentageFeatureSet(DN2_Subset_Male, pattern = "^MT-")

#violin plot to look at percent.mt, nFeature_RNA
VlnPlot(DN2_Powell_Male, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
DN2_Subset_Male <- subset(DN2_Subset_Male, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
DN2_Subset_Male <- NormalizeData(DN2_Subset_Male)

#identify highly variably features (genes) total of 2000
DN2_Subset_Male <- FindVariableFeatures(DN2_Subset_Male, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
all.genes.DN2 <- rownames(DN2_Subset_Male)
DN2_Subset_Male <- ScaleData(DN2_Subset_Male, features = all.genes.DN2)

# Perform PCA 
DN2_Subset_Male <- RunPCA(DN2_Subset_Male)

#Cluster the cells
DN2_Subset_Male <- FindNeighbors(DN2_Subset_Male, dims = 1:10)
DN2_Subset_Male <- FindClusters(DN2_Subset_Male, resolution = 0.2)

#Run UMAP
DN2_Subset_Male <- RunUMAP(DN2_Subset_Male, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(DN2_Subset_Male, reduction = "umap")

saveRDS(DN2_Subset_Male, "DN2_Subset_Male_QC.RDS")

DN2_Subset_Male<-readRDS('DN2_Subset_Male_QC.RDS')

#Markers for all DN2
DotPlot(object = DN2_Subset_Male, #group.by = 'Clusters',
        features = c(
          'CD19','ITGAX','TBX21','ZEB2', 'TLR7','FCGR2B','CD69','HLA-DRA', 'HLA-DRB','CD86','FCRL5','CD22',
          'VPREB3','LTB','CXCR5','CR2','FCRL4','CD27','IGHD','TRAF5','CD24','CD38','SELL'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Look for only Significant DN2 markers
DotPlot(object = DN2_Subset_Male, #group.by = 'Clusters',
        features = c('ITGAX','TBX21','CR2','VPREB3','LTB'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Looking for Tbet+,
DotPlot(object = DN2_Subset_Male, #group.by = 'Clusters',
        features = c('TBX21'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Look for cytotoxic markers 
DotPlot(object = DN2_Subset_Male, #group.by = 'Clusters',
        features = c('GZMB', 'GZMA','GZMM','PRF1','NKG7','KLRK1','KLRD1','GNLY'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))






# Subset Cluster 0,1,2 (remove 3 from batch effect)
clusters_to_combine <- c(0,1,2)
DN2_Subset_Male_clean<-subset(x = DN2_Subset_Male, idents =clusters_to_combine)

#QC %MT^
DN2_Subset_Male_clean[["percent.mt"]] <- PercentageFeatureSet(DN2_Subset_Male_clean, pattern = "^MT-")

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
DN2_Subset_Male_clean <- subset(DN2_Subset_Male_clean, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
DN2_Subset_Male_clean <- NormalizeData(DN2_Subset_Male_clean)

#identify highly variably features (genes) total of 2000
DN2_Subset_Male_clean <- FindVariableFeatures(DN2_Subset_Male_clean, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
all.genes.DN2 <- rownames(DN2_Subset_Male_clean)
DN2_Subset_Male_clean <- ScaleData(DN2_Subset_Male_clean, features = all.genes.DN2)

# Perform PCA 
DN2_Subset_Male_clean <- RunPCA(DN2_Subset_Male_clean)

#Cluster the cells
DN2_Subset_Male_clean <- FindNeighbors(DN2_Subset_Male_clean, dims = 1:10)
DN2_Subset_Male_clean <- FindClusters(DN2_Subset_Male_clean, resolution = 0.2)

#Run UMAP
DN2_Subset_Male_clean <- RunUMAP(DN2_Subset_Male_clean, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(DN2_Subset_Male_clean, reduction = "umap")

saveRDS(DN2_Subset_Male_clean, "DN2_Subset_Male_clean_QC.RDS")

DN2_Subset_Male_clean<-readRDS('DN2_Subset_Male_clean_QC.RDS')

#Rename Cluster
Idents(DN2_Subset_Male_clean) <- "seurat_clusters"
current.cluster.ids <- levels(DN2_Subset_Male_clean)
current.cluster.ids
new.cluster.ids<-c("DN2.B.male", "DN2.A.male", "DN2.Cytotoxic.male")

DN2_Subset_Male_clean@meta.data$DN2_Subtypes<-DN2_Subset_Male_clean@meta.data$seurat_clusters
DN2_Subset_Male_clean@meta.data$DN2_Subtypes<-plyr::mapvalues(x=DN2_Subset_Male_clean@meta.data$DN2_Subtypes, from=current.cluster.ids, to=new.cluster.ids)

DimPlot(DN2_Subset_Male_clean, group.by="DN2_Subtypes")

saveRDS(DN2_Subset_Male_clean, "DN2_Male_Subtypes.RDS")

DN2_Subset_Male_clean<-readRDS('DN2_Male_Subtypes.RDS')




#add study id
DN2_Powell_Male@meta.data$Atlas_Study_ID<-'Powell_2022'
DN2_Powell_Male@meta.data$Atlas_Disease<-'Healthy'
META<-DN2_Powell_Male[[]]

#see what column labels are stored for individual info
unique(DN2_Powell_Male$sex)
Idents(DN2_Powell_Male)<-'sex'
current.cluster.ids<-levels(DN2_Powell_Male@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_Powell_Male@active.ident)

DN2_Powell_Male@meta.data$Atlas_Sex<-DN2_Powell_Male@meta.data$sex
DN2_Powell_Male@meta.data$Atlas_Sex<-plyr::mapvalues(x=DN2_Powell_Male@meta.data$Atlas_Sex, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Powell_Male$Atlas_Sex)

#Ethnicity- Broad
unique(DN2_Powell_Male$ethnicity)
Idents(DN2_Powell_Male)<-'ethnicity'
current.cluster.ids<-levels(DN2_Powell_Male@active.ident)
current.cluster.ids<-c("European" )
new.cluster.ids<-c("EUR" )

DN2_Powell_Male@meta.data$Atlas_Ethnicity_Broad<-DN2_Powell_Male@meta.data$ethnicity
DN2_Powell_Male@meta.data$Atlas_Ethnicity_Broad<-plyr::mapvalues(x=DN2_Powell_Male@meta.data$Atlas_Ethnicity_Broad, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Powell_Male$Atlas_Ethnicity_Broad)

#Ethnicity- Fine
DN2_Powell_Male@meta.data$Atlas_Ethnicity_Fine<-DN2_Powell_Male@meta.data$ethnicity
DN2_Powell_Male@meta.data$Atlas_Ethnicity_Fine<-plyr::mapvalues(x=DN2_Powell_Male@meta.data$Atlas_Ethnicity_Fine, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Powell_Male$Atlas_Ethnicity_Fine)

#Atlas Donor
unique(DN2_Powell_Male$individual)
Idents(DN2_Powell_Male)<-'individual'
current.cluster.ids<-levels(DN2_Powell_Male@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_Powell_Male@active.ident)

DN2_Powell_Male@meta.data$Atlas_Donor<-DN2_Powell_Male@meta.data$'individual'
DN2_Powell_Male@meta.data$Atlas_Donor<-plyr::mapvalues(x=DN2_Powell_Male@meta.data$Atlas_Donor, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Powell_Male$Atlas_Donor)

#Age Updated
META<-DN2_Powell_Male[[]]
write.csv(META,'Powell_Male_cell_id_meta.csv')
Powell_Male_cell_id_meta<-read.csv('Powell_Male_cell_id_meta.csv', row.names = NULL)

donor_id_meta<-(META%>%distinct(Atlas_Donor, .keep_all=TRUE))
write.csv(donor_id_meta,'Powell_Male_donor_id_meta.csv')
donor_id_meta<-read.csv('Powell_Male_donor_id_meta.csv')

merged_df<- merge(donor_id_meta, Powell_Male_cell_id_meta, by= 'Atlas_Donor', all=TRUE)
#check to see which column has the cell ids/barcodes
write.csv(merged_df,'merged_df.csv')
merged_df2<-read.csv('merged_df.csv', row.names = 5) #row.names = column with cell barcodes/ids
drop <- c("X.1")
merged_df2 = merged_df2[,!(names(merged_df2) %in% drop)]


DN2_Powell_Male<-AddMetaData(DN2_Powell_Male,merged_df2 )
unique(DN2_Powell_Male$Atlas_Age_Category)
saveRDS(DN2_Powell_Male, "DN2_Powell_Male_updated_meta.rds")
DN2_Powell_Male<-readRDS('DN2_Powell_Male_updated_meta.RDS')





#add study id
DN2_Subset_Male_clean@meta.data$Atlas_Study_ID<-'Powell_2022'
DN2_Subset_Male_clean@meta.data$Atlas_Disease<-'Healthy'
META<-DN2_Subset_Male_clean[[]]

#see what column labels are stored for individual info
unique(DN2_Subset_Male_clean$sex)
Idents(DN2_Subset_Male_clean)<-'sex'
current.cluster.ids<-levels(DN2_Subset_Male_clean@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_Subset_Male_clean@active.ident)

DN2_Subset_Male_clean@meta.data$Atlas_Sex<-DN2_Subset_Male_clean@meta.data$sex
DN2_Subset_Male_clean@meta.data$Atlas_Sex<-plyr::mapvalues(x=DN2_Subset_Male_clean@meta.data$Atlas_Sex, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Subset_Male_clean$Atlas_Sex)

#Ethnicity- Broad
unique(DN2_Subset_Male_clean$ethnicity)
Idents(DN2_Subset_Male_clean)<-'ethnicity'
current.cluster.ids<-levels(DN2_Subset_Male_clean@active.ident)
current.cluster.ids<-c("European" )
new.cluster.ids<-c("EUR" )

DN2_Subset_Male_clean@meta.data$Atlas_Ethnicity_Broad<-DN2_Subset_Male_clean@meta.data$ethnicity
DN2_Subset_Male_clean@meta.data$Atlas_Ethnicity_Broad<-plyr::mapvalues(x=DN2_Subset_Male_clean@meta.data$Atlas_Ethnicity_Broad, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Subset_Male_clean$Atlas_Ethnicity_Broad)

#Ethnicity- Fine
DN2_Subset_Male_clean@meta.data$Atlas_Ethnicity_Fine<-DN2_Subset_Male_clean@meta.data$ethnicity
DN2_Subset_Male_clean@meta.data$Atlas_Ethnicity_Fine<-plyr::mapvalues(x=DN2_Subset_Male_clean@meta.data$Atlas_Ethnicity_Fine, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Subset_Male_clean$Atlas_Ethnicity_Fine)

#Atlas Donor
unique(DN2_Subset_Male_clean$individual)
Idents(DN2_Subset_Male_clean)<-'individual'
current.cluster.ids<-levels(DN2_Subset_Male_clean@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_Subset_Male_clean@active.ident)

DN2_Subset_Male_clean@meta.data$Atlas_Donor<-DN2_Subset_Male_clean@meta.data$'individual'
DN2_Subset_Male_clean@meta.data$Atlas_Donor<-plyr::mapvalues(x=DN2_Subset_Male_clean@meta.data$Atlas_Donor, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Subset_Male_clean$Atlas_Donor)

#Age Updated
META<-DN2_Subset_Male_clean[[]]
write.csv(META,'Powell_Male_clean_cell_id_meta.csv')
Powell_Male_clean_cell_id_meta<-read.csv('Powell_Male_cell_id_meta.csv', row.names = NULL)

donor_id_meta<-(META%>%distinct(Atlas_Donor, .keep_all=TRUE))
write.csv(donor_id_meta,'Powell_Male_clean_donor_id_meta.csv')
donor_id_meta<-read.csv('Powell_Male_clean_donor_id_meta.csv')

merged_df<- merge(donor_id_meta, Powell_Male_clean_cell_id_meta, by= 'Atlas_Donor', all=TRUE)
#check to see which column has the cell ids/barcodes
write.csv(merged_df,'merged_df.csv')
merged_df2<-read.csv('merged_df.csv', row.names = 5) #row.names = column with cell barcodes/ids
drop <- c("X.1")
merged_df2 = merged_df2[,!(names(merged_df2) %in% drop)]


DN2_Subset_Male_clean<-AddMetaData(DN2_Subset_Male_clean,merged_df2 )
unique(DN2_Subset_Male_clean$Atlas_Age_Category)
saveRDS(DN2_Subset_Male_clean, "DN2_Powell_Male_clean_updated_meta.rds")
DN2_Subset_Male_clean<-readRDS('DN2_Powell_Male_clean_updated_meta.RDS')



#no cytotoxic

# Subset Cluster 0,1,2 (remove 3 from batch effect)
clusters_to_combine <- c(0,1)
DN2_Subset_Male_clean<-subset(x = DN2_Subset_Male, idents =clusters_to_combine)

#QC %MT^
DN2_Subset_Male_clean[["percent.mt"]] <- PercentageFeatureSet(DN2_Subset_Male_clean, pattern = "^MT-")

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
DN2_Subset_Male_clean <- subset(DN2_Subset_Male_clean, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
DN2_Subset_Male_clean <- NormalizeData(DN2_Subset_Male_clean)

#identify highly variably features (genes) total of 2000
DN2_Subset_Male_clean <- FindVariableFeatures(DN2_Subset_Male_clean, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
DN2_Subset_Male_clean <- ScaleData(DN2_Subset_Male_clean)

# Perform PCA 
DN2_Subset_Male_clean <- RunPCA(DN2_Subset_Male_clean)

#Cluster the cells
DN2_Subset_Male_clean <- FindNeighbors(DN2_Subset_Male_clean, dims = 1:10)
DN2_Subset_Male_clean <- FindClusters(DN2_Subset_Male_clean, resolution = 0.2)

#Run UMAP
DN2_Subset_Male_clean <- RunUMAP(DN2_Subset_Male_clean, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(DN2_Subset_Male_clean, reduction = "umap")

saveRDS(DN2_Subset_Male_clean, "DN2_Subset_Male_clean_QC.RDS")

DN2_Subset_Male_clean<-readRDS('DN2_Subset_Male_clean_QC.RDS')






# Subset Cluster 7
clusters_to_combine <- c(7)
DN2_Subset_Male<-subset(x = DN2_Powell_Male, idents =clusters_to_combine)


#QC %MT^
DN2_Subset_Male[["percent.mt"]] <- PercentageFeatureSet(DN2_Subset_Male, pattern = "^MT-")

#violin plot to look at percent.mt, nFeature_RNA
VlnPlot(DN2_Powell_Male, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
DN2_Subset_Male <- subset(DN2_Subset_Male, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
DN2_Subset_Male <- NormalizeData(DN2_Subset_Male)

#identify highly variably features (genes) total of 2000
DN2_Subset_Male <- FindVariableFeatures(DN2_Subset_Male, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
DN2_Subset_Male <- ScaleData(DN2_Subset_Male)

# Perform PCA 
DN2_Subset_Male <- RunPCA(DN2_Subset_Male)

#Cluster the cells
DN2_Subset_Male <- FindNeighbors(DN2_Subset_Male, dims = 1:10)
DN2_Subset_Male <- FindClusters(DN2_Subset_Male, resolution = 0.2)

#Run UMAP
DN2_Subset_Male <- RunUMAP(DN2_Subset_Male, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(DN2_Subset_Male, reduction = "umap")

saveRDS(DN2_Subset_Male, "DN2_Subset_Male_QC.RDS")

DN2_Subset_Male<-readRDS('DN2_Subset_Male_QC.RDS')