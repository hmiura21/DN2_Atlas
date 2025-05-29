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

DN2_Jimmie_Healthy <-readRDS('/Users/honokamiura/Downloads/Research/DN2/Jimmie_Ye_healthy_B_cells.rds')

#QC %MT^
DN2_Jimmie_Healthy[["percent.mt"]] <- PercentageFeatureSet(DN2_Jimmie_Healthy, pattern = "^MT.")

#violin plot to look at percent.mt, nFeature_RNA
VlnPlot(DN2_Jimmie_Healthy, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
DN2_Jimmie_Healthy <- subset(DN2_Jimmie_Healthy, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
DN2_Jimmie_Healthy <- NormalizeData(DN2_Jimmie_Healthy)

#identify highly variably features (genes) total of 2000
DN2_Jimmie_Healthy <- FindVariableFeatures(DN2_Jimmie_Healthy, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
DN2_Jimmie_Healthy <- ScaleData(DN2_Jimmie_Healthy)

# Perform PCA and color by cell cycle phase
DN2_Jimmie_Healthy <- RunPCA(DN2_Jimmie_Healthy)

#Cluster the cells
DN2_Jimmie_Healthy <- FindNeighbors(DN2_Jimmie_Healthy, dims = 1:10)
DN2_Jimmie_Healthy <- FindClusters(DN2_Jimmie_Healthy, resolution = 0.8)

#Run UMAP
DN2_Jimmie_Healthy <- RunUMAP(DN2_Jimmie_Healthy, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(DN2_Jimmie_Healthy, reduction = "umap")

saveRDS(DN2_Jimmie_Healthy, "DN2_Jimmie_Healthy_QC.RDS")

#Start here: 
DN2_Jimmie_Healthy<-readRDS('DN2_Jimmie_Healthy_QC.RDS')

Idents(DN2_Jimmie_Healthy)<-'seurat_clusters'

META<-DN2_Jimmie_Healthy[[]]
#Looking for CD11c+(ITGAX), Tbet+ (TBX21), CD21-(CR2), Vreb3-, Ltb-,
#CD32b=FCGR2B, HLA-DR=HLA-DRA and HLA-DRB, IGD=IGHD, CD62L=SELL
DotPlot(object = DN2_Jimmie_Healthy, #group.by = 'Clusters',
        features = c(
          'CD19','ITGAX','TBX21','ZEB2', 'TLR7','FCGR2B','CD69','HLA-DRA', 'HLA-DRB','CD86','FCRL5','CD22',
          'VPREB3','LTB','CXCR5','CR2','FCRL4','CD27','IGHD','TRAF5','CD24','CD38','SELL'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#dn2 markers
DotPlot(object = DN2_Jimmie_Healthy, group.by = 'seurat_clusters', 
        features = c("SOX5",'HLA-DRB5','TBX21',"ITGAX","BATF",'JAZF1','NFATC2',
                     'ZEB2','NR4A2','MS4A1','TOX','FGR','IL10RA','FCRL5',
                     'SREBF2','TFEB','TFEC','ZBTB32','MEF2A','POU2F2','SPI1',
                     #NEGATIVE
                     'LTB','VPREB3','SELL','JCHAIN','CD27','IGHD'),
        cols = c('RdBu') ) + RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 7,face = "bold"))


#Look for only Significant DN2 markers
DotPlot(object = DN2_Jimmie_Healthy, #group.by = 'Clusters',
        features = c('ITGAX','TBX21','CR2','VPREB3','LTB'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Looking for Tbet+,
DotPlot(object = DN2_Jimmie_Healthy, #group.by = 'Clusters',
        features = c('TBX21'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Look for cytotoxic markers 
DotPlot(object = DN2_Jimmie_Healthy, #group.by = 'Clusters',
        features = c('GZMB', 'GZMA','GZMM','PRF1','NKG7','KLRK1','KLRD1','GNLY'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#violin plot to look at cytotoxicity by cluster
VlnPlot(DN2_Jimmie_Healthy, features= c('GZMB', 'GZMA','GZMM','PRF1','NKG7','KLRK1','KLRD1','GNLY'), idents=c('0','1','2','3','4','5','6','7','8','9','10','11'))





# Subset Cluster 11
DN2_Jimmie_Healthy_Subset11<-subset(x = DN2_Jimmie_Healthy, idents ='11')

#QC %MT^
DN2_Jimmie_Healthy_Subset11[["percent.mt"]] <- PercentageFeatureSet(DN2_Jimmie_Healthy_Subset11, pattern = "^MT-")

#QC to eliminate anything below 200 and above 2500 features and %mt of more than 5%
DN2_Jimmie_Healthy_Subset11 <- subset(DN2_Jimmie_Healthy_Subset11, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization at default
DN2_Jimmie_Healthy_Subset11 <- NormalizeData(DN2_Jimmie_Healthy_Subset11)

#identify highly variably features (genes) total of 2000
DN2_Jimmie_Healthy_Subset11 <- FindVariableFeatures(DN2_Jimmie_Healthy_Subset11, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
DN2_Jimmie_Healthy_Subset11 <- ScaleData(DN2_Jimmie_Healthy_Subset11)

# Perform PCA 
DN2_Jimmie_Healthy_Subset11 <- RunPCA(DN2_Jimmie_Healthy_Subset11)

#Cluster the cells
DN2_Jimmie_Healthy_Subset11 <- FindNeighbors(DN2_Jimmie_Healthy_Subset11, dims = 1:10)
DN2_Jimmie_Healthy_Subset11 <- FindClusters(DN2_Jimmie_Healthy_Subset11, resolution = 0.2)

#Run UMAP
DN2_Jimmie_Healthy_Subset11 <- RunUMAP(DN2_Jimmie_Healthy_Subset11, dims = 1:10)

#Dim Plot showing Clusters
DimPlot(DN2_Jimmie_Healthy_Subset11, reduction = "umap")

saveRDS(DN2_Jimmie_Healthy_Subset11, "DN2_Jimmie_Healthy_Subset11_QC.RDS")

DN2_Jimmie_Healthy_Subset11<-readRDS('DN2_Jimmie_Healthy_Subset11_QC.RDS')

#compare cluster 0 to cluster 1
DN2_Jimmie_Healthy_0v1 <- FindMarkers(DN2_Jimmie_Healthy_Subset11, only.pos = F,
                                                ident.1 = ('0'),
                                                ident.2 = c('1'),
                                                min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_Jimmie_Healthy_0v1,"/Users/honokamiura/Downloads/Research/DN2/DN2_Jimmie_Healthy_0v1.csv")

#Stacked Barplot for sex vs cluster
?dittoBarPlot()
dittoBarPlot(
  object = DN2_Jimmie_Healthy_Subset11,
  var = "sex", 
  group.by = "seurat_clusters", 
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

DimPlot(object = DN2_Jimmie_Healthy_Subset11, reduction = "umap", shuffle = F, group.by="sex", 
        cols=c('female'='snow',
               'male'='blue'))

#Rename Cluster
Idents(DN2_Jimmie_Healthy_Subset11) <- "seurat_clusters"
current.cluster.ids <- levels(DN2_Jimmie_Healthy_Subset11)
current.cluster.ids
new.cluster.ids<-c("DN2.Jimmie.A", "DN2.Jimmie.B")

DN2_Jimmie_Healthy_Subset11@meta.data$DN2_Subtypes<-DN2_Jimmie_Healthy_Subset11@meta.data$seurat_clusters
DN2_Jimmie_Healthy_Subset11@meta.data$DN2_Subtypes<-plyr::mapvalues(x=DN2_Jimmie_Healthy_Subset11@meta.data$DN2_Subtypes, from=current.cluster.ids, to=new.cluster.ids)

DimPlot(DN2_Jimmie_Healthy_Subset11, group.by="DN2_Subtypes")

saveRDS(DN2_Jimmie_Healthy_Subset11, "DN2_Jimmie_Subtypes.RDS")

DN2_Jimmie_Healthy_Subset11<-readRDS('DN2_Jimmie_Subtypes.RDS')



#add study id
DN2_Jimmie_Healthy@meta.data$Atlas_Study_ID<-'Jimmie_YE_2022'
META<-DN2_Jimmie_Healthy[[]]

DN2_Jimmie_Healthy@meta.data$Atlas_Disease<-'Healthy'

#see what column labels are stored for individual info
unique(DN2_Jimmie_Healthy$sex)
Idents(DN2_Jimmie_Healthy)<-'sex'
current.cluster.ids<-levels(DN2_Jimmie_Healthy@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_Jimmie_Healthy@active.ident)

DN2_Jimmie_Healthy@meta.data$Atlas_Sex<-DN2_Jimmie_Healthy@meta.data$sex
DN2_Jimmie_Healthy@meta.data$Atlas_Sex<-plyr::mapvalues(x=DN2_Jimmie_Healthy@meta.data$Atlas_Sex, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Jimmie_Healthy$Atlas_Sex)

#Ethnicity- Broad
unique(DN2_Jimmie_Healthy$ethnicity)
Idents(DN2_Jimmie_Healthy)<-'ethnicity'
current.cluster.ids<-levels(DN2_Jimmie_Healthy@active.ident)
current.cluster.ids<-c( "European"   ,                "Asian"  ,                   
                       "Hispanic or Latin American")
new.cluster.ids<-c( "EUR"   ,                "ASI"  ,                   
                    "HIS")

DN2_Jimmie_Healthy@meta.data$Atlas_Ethnicity_Broad<-DN2_Jimmie_Healthy@meta.data$ethnicity
DN2_Jimmie_Healthy@meta.data$Atlas_Ethnicity_Broad<-plyr::mapvalues(x=DN2_Jimmie_Healthy@meta.data$Atlas_Ethnicity_Broad, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Jimmie_Healthy$Atlas_Ethnicity_Broad)

#Ethnicity- Fine
DN2_Jimmie_Healthy@meta.data$Atlas_Ethnicity_Fine<-DN2_Jimmie_Healthy@meta.data$ethnicity
DN2_Jimmie_Healthy@meta.data$Atlas_Ethnicity_Fine<-plyr::mapvalues(x=DN2_Jimmie_Healthy@meta.data$Atlas_Ethnicity_Fine, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Jimmie_Healthy$Atlas_Ethnicity_Fine)

#Atlas Donor
unique(DN2_Jimmie_Healthy$donor_uuid)
Idents(DN2_Jimmie_Healthy)<-'donor_uuid'
current.cluster.ids<-levels(DN2_Jimmie_Healthy@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_Jimmie_Healthy@active.ident)

DN2_Jimmie_Healthy@meta.data$Atlas_Donor<-DN2_Jimmie_Healthy@meta.data$'donor_uuid'
DN2_Jimmie_Healthy@meta.data$Atlas_Donor<-plyr::mapvalues(x=DN2_Jimmie_Healthy@meta.data$Atlas_Donor, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Jimmie_Healthy$Atlas_Donor)

#Age 
META<-DN2_Jimmie_Healthy[[]]
donor_id_meta<-(META%>%distinct(Atlas_Donor, .keep_all=TRUE))

write.csv(donor_id_meta,'Jimmie_donor_id_meta.csv')
donor_id_meta<-read.csv('Jimmie_donor_id_meta.csv')
merged_df<- merge(donor_id_meta, META, by= 'Atlas_Donor', all=TRUE)

DN2_Jimmie_Healthy<-AddMetaData(DN2_Jimmie_Healthy,merged_df )
unique(DN2_Jimmie_Healthy$Atlas_Age_Category)

saveRDS(DN2_Jimmie_Healthy, "DN2_Jimmie_Healthy_updatedmeta.RDS")
DN2_Jimmie_Healthy<-readRDS('DN2_Jimmie_Healthy_updatedmeta.RDS')

#Age Updated
META<-DN2_Jimmie_Healthy[[]]
write.csv(META,'Jimmie_cell_id_meta.csv')
Jimmie_cell_id_meta<-read.csv('Jimmie_cell_id_meta.csv', row.names = NULL)

donor_id_meta<-(META%>%distinct(Atlas_Donor, .keep_all=TRUE))
write.csv(donor_id_meta,'Jimmie_donor_id_meta.csv')
donor_id_meta<-read.csv('Jimmie_donor_id_meta.csv')

merged_df<- merge(donor_id_meta, Jimmie_cell_id_meta, by= 'Atlas_Donor', all=TRUE)
#check to see which column has the cell ids/barcodes
write.csv(merged_df,'merged_df.csv')
merged_df2<-read.csv('merged_df.csv', row.names = 5) #row.names = column with cell barcodes/ids
drop <- c("X.1")
merged_df2 = merged_df2[,!(names(merged_df2) %in% drop)]


DN2_Jimmie_Healthy<-AddMetaData(DN2_Jimmie_Healthy,merged_df2 )
unique(DN2_Jimmie_Healthy$Atlas_Age_Category)
DimPlot(DN2_Jimmie_Healthy, group.by = 'seurat_clusters')
saveRDS(DN2_Jimmie_Healthy, "DN2_Jimmie_Healthy_updated_meta.rds")
DN2_Jimmie_Healthy<-readRDS('DN2_Jimmie_Healthy_updated_meta.RDS')



#add study id
DN2_Jimmie_Healthy_Subset11@meta.data$Atlas_Study_ID<-'Jimmie_YE_2022'
META<-DN2_Jimmie_Healthy_Subset11[[]]

DN2_Jimmie_Healthy_Subset11@meta.data$Atlas_Disease<-'Healthy'

#see what column labels are stored for individual info
unique(DN2_Jimmie_Healthy_Subset11$sex)
Idents(DN2_Jimmie_Healthy_Subset11)<-'sex'
current.cluster.ids<-levels(DN2_Jimmie_Healthy_Subset11@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_Jimmie_Healthy_Subset11@active.ident)

DN2_Jimmie_Healthy_Subset11@meta.data$Atlas_Sex<-DN2_Jimmie_Healthy_Subset11@meta.data$sex
DN2_Jimmie_Healthy_Subset11@meta.data$Atlas_Sex<-plyr::mapvalues(x=DN2_Jimmie_Healthy_Subset11@meta.data$Atlas_Sex, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Jimmie_Healthy_Subset11$Atlas_Sex)

#Ethnicity- Broad
unique(DN2_Jimmie_Healthy_Subset11$ethnicity)
Idents(DN2_Jimmie_Healthy_Subset11)<-'ethnicity'
current.cluster.ids<-levels(DN2_Jimmie_Healthy_Subset11@active.ident)
current.cluster.ids<-c( "European"   ,                "Asian"  ,                   
                        "Hispanic or Latin American")
new.cluster.ids<-c( "EUR"   ,                "ASI"  ,                   
                    "HIS")

DN2_Jimmie_Healthy_Subset11@meta.data$Atlas_Ethnicity_Broad<-DN2_Jimmie_Healthy_Subset11@meta.data$ethnicity
DN2_Jimmie_Healthy_Subset11@meta.data$Atlas_Ethnicity_Broad<-plyr::mapvalues(x=DN2_Jimmie_Healthy_Subset11@meta.data$Atlas_Ethnicity_Broad, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Jimmie_Healthy_Subset11$Atlas_Ethnicity_Broad)

#Ethnicity- Fine
DN2_Jimmie_Healthy_Subset11@meta.data$Atlas_Ethnicity_Fine<-DN2_Jimmie_Healthy_Subset11@meta.data$ethnicity
DN2_Jimmie_Healthy_Subset11@meta.data$Atlas_Ethnicity_Fine<-plyr::mapvalues(x=DN2_Jimmie_Healthy_Subset11@meta.data$Atlas_Ethnicity_Fine, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Jimmie_Healthy_Subset11$Atlas_Ethnicity_Fine)

#Atlas Donor
unique(DN2_Jimmie_Healthy_Subset11$donor_uuid)
Idents(DN2_Jimmie_Healthy_Subset11)<-'donor_uuid'
current.cluster.ids<-levels(DN2_Jimmie_Healthy_Subset11@active.ident)
current.cluster.ids
new.cluster.ids<-levels(DN2_Jimmie_Healthy_Subset11@active.ident)

DN2_Jimmie_Healthy_Subset11@meta.data$Atlas_Donor<-DN2_Jimmie_Healthy_Subset11@meta.data$'donor_uuid'
DN2_Jimmie_Healthy_Subset11@meta.data$Atlas_Donor<-plyr::mapvalues(x=DN2_Jimmie_Healthy_Subset11@meta.data$Atlas_Donor, from=current.cluster.ids, to=new.cluster.ids)

unique(DN2_Jimmie_Healthy_Subset11$Atlas_Donor)

#Age Updated
META<-DN2_Jimmie_Healthy_Subset11[[]]
write.csv(META,'Jimmie_Subset11_cell_id_meta.csv')
Jimmie_Subset11_cell_id_meta<-read.csv('Jimmie_Subset11_cell_id_meta.csv', row.names = NULL)

donor_id_meta<-(META%>%distinct(Atlas_Donor, .keep_all=TRUE))
write.csv(donor_id_meta,'Jimmie_Subset11_donor_id_meta.csv')
donor_id_meta<-read.csv('Jimmie_Subset11_donor_id_meta.csv')

merged_df<- merge(donor_id_meta, Jimmie_Subset11_cell_id_meta, by= 'Atlas_Donor', all=TRUE)
#check to see which column has the cell ids/barcodes
write.csv(merged_df,'merged_df.csv')
merged_df2<-read.csv('merged_df.csv', row.names = 5) #row.names = column with cell barcodes/ids
drop <- c("X.1")
merged_df2 = merged_df2[,!(names(merged_df2) %in% drop)]


DN2_Jimmie_Healthy_Subset11<-AddMetaData(DN2_Jimmie_Healthy_Subset11,merged_df2 )
unique(DN2_Jimmie_Healthy_Subset11$Atlas_Age_Category)
saveRDS(DN2_Jimmie_Healthy_Subset11, "DN2_Jimmie_Healthy_Subset11_updated_meta.rds")
DN2_Jimmie_Healthy_Subset11<-readRDS('DN2_Jimmie_Healthy_Subset11_updated_meta.RDS')




#SHINYAPP
#get packages ready
library(rsconnect)
library(Seurat)
library(ShinyCell)

#connect to account
rsconnect::setAccountInfo(name='hmiura1221',
                          token='93EAD3895AF84DA902AAD43FDAF30183',
                          secret='JhJweWwnaU166CRJNUbLcMacrHjXuRBeFkr0dz5r')

scConf = createConfig(DN2_Jimmie_Healthy)
write.csv(scConf,'scConf_meta_dn2_shiny.csv')
#scConf= modColours(scConf,  meta.to.mod = 'sex',
#                   new.colours = c('red3','blue3'))
makeShinyApp(DN2_Combined, scConf, gene.mapping = TRUE,
             shiny.title = "DN2 B Cell Characterization for Jimme Ye Data Set", #title for overall page
             default.gene1 = "TBX21", #default first two data you see (dn2, and cytotoxicity)
             shiny.dir = "shinyApp_DN2_JimmieYe/") 

rsconnect::deployApp("shinyApp_DN2_JimmieYe/") #done; only when ready!

