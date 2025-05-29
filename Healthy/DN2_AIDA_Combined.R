
#COMBINE MALE AND FEMALE
#merge male dn2s with female dn2s
DN2_AIDA_Combined<-merge(DN2_AIDA_Female_Subset6,DN2_AIDA_Male_Subset7)

#normalize, scale, pca, etc

#Normalization at default
DN2_AIDA_Combined <- NormalizeData(DN2_AIDA_Combined)

#identify highly variably features (genes) total of 2000
DN2_AIDA_Combined <- FindVariableFeatures(DN2_AIDA_Combined, selection.method = "vst", nfeatures = 2000)

#Scaling data by making mean expression 0 and variance is 1. allow highly expressed genes to dominate
all.genes <- rownames(DN2_AIDA_Combined)
DN2_AIDA_Combined <- ScaleData(DN2_AIDA_Combined, features = all.genes)

# Perform PCA and color by cell cycle phase
DN2_AIDA_Combined <- RunPCA(DN2_AIDA_Combined)

#Cluster the cells
DN2_AIDA_Combined <- FindNeighbors(DN2_AIDA_Combined, dims = 1:10)
DN2_AIDA_Combined <- FindClusters(DN2_AIDA_Combined, resolution = 0.1)

#Run UMAP
DN2_AIDA_Combined <- RunUMAP(DN2_AIDA_Combined, dims = 1:10)

#DimPlot showing Clusters
DimPlot(DN2_AIDA_Combined, reduction = "umap")

#save and read
saveRDS(DN2_AIDA_Combined, "DN2_AIDA_Combined.RDS")
DN2_AIDA_Combined<-readRDS('DN2_AIDA_Combined.RDS')


#dimplot dn2 subtypes (split.by)
DimPlot(DN2_AIDA_Combined, reduction = "umap", split.by='sex')

DimPlot(object = DN2_AIDA_Combined, reduction = "umap", shuffle = F, group.by="sex", 
        cols=c('female'='red',
               'male'='blue'))

DN2_AIDA_Combined@meta.data

#A
DimPlot(object = DN2_AIDA_Combined, reduction = "umap", shuffle = F, group.by="DN2_Subtypes",
        cols=c('DN2.AIDA.Male.A'='blue',
               'DN2.AIDA.Male.B'='snow2',
               'DN2.AIDA.Female.A'='red',
               'DN2.AIDA.Female.B'='snow2'))
#B
DimPlot(object = DN2_AIDA_Combined, reduction = "umap", shuffle = F, group.by="DN2_Subtypes", 
        cols=c('DN2.AIDA.Male.A'='snow2',
               'DN2.AIDA.Male.B'='blue',
               'DN2.AIDA.Female.A'='snow2',
               'DN2.AIDA.Female.B'='red'))

#UP UNTIL HERE

DotPlot(object = DN2_AIDA_Combined, #group.by = 'Clusters',
        features = c(
          'CD19','ITGAX','TBX21','ZEB2', 'TLR7','FCGR2B','CD69','HLA-DRA', 'HLA-DRB','CD86','FCRL5','CD22',
          'VPREB3','LTB','CXCR5','CR2','FCRL4','CD27','IGHD','TRAF5','CD24','CD38','SELL'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Look for only Significant DN2 markers
DotPlot(object = DN2_AIDA_Combined, #group.by = 'Clusters',
        features = c('ITGAX','TBX21','CR2','VPREB3','LTB'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#Looking for Tbet+,
DotPlot(object = DN2_AIDA_Combined, group.by = 'DN2_Subtypes',
        features = c('TBX21'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

DotPlot(object=DN2_AIDA_Combined, group.by='justCytotoxic',
        features = c('GZMB', 'GZMA','GZMM','PRF1','NKG7','KLRK1','KLRD1','GNLY'
        ),
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#violin plot to look at cytotoxicity by sex
?VlnPlot
VlnPlot(DN2_AIDA_Combined, features= c('GZMB', 'GZMA','GZMM','PRF1','NKG7','KLRK1','KLRD1','GNLY'), idents='2',split.by = 'sex') & theme(legend.position = 'right')

?VlnPlot()

#Find Gene Markers
#compare one group to another
DN2_AIDA_Combined_DN2.A_Markers <- FindMarkers(DN2_AIDA_Combined, only.pos = F,
                                          ident.1 = ('0'),
                                          ident.2 = c('1','2'),
                                          min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_AIDA_Combined_DN2.A_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_AIDA_Combined_DN2.A_Markers.csv")


DN2_AIDA_Combined_DN2.B_Markers <- FindMarkers(DN2_AIDA_Combined, only.pos = F,
                                          ident.1 = ('1'),
                                          ident.2 = c('0','2'),
                                          min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_AIDA_Combined_DN2.B_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_AIDA_Combined_DN2.B_Markers.csv")


DN2_AIDA_Combined_DN2.cytotoxic_Markers <- FindMarkers(DN2_AIDA_Combined, only.pos = F,
                                                  ident.1 = ('2'),
                                                  ident.2 = c('0','1'),
                                                  min.pct = 0.1, logfc.threshold = 0.5)
write.csv(DN2_AIDA_Combined_DN2.cytotoxic_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_AIDA_Combined_DN2.cytotoxic_Markers.csv")


#differentiate by sex for DN2. B
DN2_AIDA_Combined_DN2.B_sex_Markers <- FindMarkers(DN2_AIDA_Combined, only.pos = F,
                                              subset.ident = ('1'), group.by = 'sex',
                                              ident.1 = c('male'), ident.2 = 'female',
                                              min.pct = 0.1, logfc.threshold = 0.25)


write.csv(DN2_AIDA_Combined_DN2.B_sex_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_AIDA_Combined_DN2.B_sex_Markers.csv")

#differentiate by sex for ALL
DN2_AIDA_Combined_DN2.B_sex_all_Markers <- FindMarkers(DN2_AIDA_Combined, only.pos = F,
                                                  group.by = 'sex',
                                                  ident.1 = c('male'), ident.2 = 'female',
                                                  min.pct = 0.1, logfc.threshold = 0.25)


write.csv(DN2_AIDA_Combined_DN2.B_sex_all_Markers,"/Users/honokamiura/Downloads/Research/DN2/DN2_AIDA_Combined_DN2.B_sex_all_Markers.csv")


# Markers to Show: Density Plot
markers <- list()
markers$DN2.A<- c('CLEC2D', 'IGHD', 'CD79A', 'FCRL5', 'FCRL3', 'FGR', 'LAPTM5')
markers$DN2.B<- c("CD79A", "CIB1", "PSAP", "CD19", "HLA-DPA1", "ITGB2", "HLA-DPB1", "ZEB2", "HCK", "IGKC")
markers$DN2.C.female<-c("TSPO", "TAX1BP3", "RHOA", "IFI30", "GSTP1", "S100A11", "JUN", "CIB1", "SPINT2", "LAPTM5")
markers$DN2.Cytotoxic<-c('GZMB', 'GZMA','GZMM','PRF1','NKG7','KLRK1','KLRD1','GNLY')

DN2_AIDA_Combined <- AddModuleScore_UCell(DN2_AIDA_Combined,features = markers)

#Marker DN2.A
DN2_AIDA_Combined$DN2.A_UCell
plot_density(
  DN2_AIDA_Combined,
  c('FCRL5'),
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

#Marker DN2.B
DN2_AIDA_Combined$DN2.B_UCell
plot_density(
  DN2_AIDA_Combined,
  c('DN2.B_UCell'),
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

#Marker DN2.A
DN2_AIDA_Combined$DN2.C.female_UCell
plot_density(
  DN2_AIDA_Combined,
  c('DN2.C.female_UCell'),
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

#Marker DN2.B. Female
DN2_AIDA_Combined$DN2.B_UCell
plot_density(
  DN2_AIDA_Combined,
  c('DN2.B_UCell'),
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


library(dittoSeq)
library(ggpubr)
library(cowplot)

#Stacked Barplot for Individuals
?dittoBarPlot()
dittoBarPlot(
  object = DN2_AIDA_Combined,
  var = "individual", 
  group.by = "seurat_clusters", 
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))&NoLegend()

#Stacked Barplot for sex vs cluster
?dittoBarPlot()
dittoBarPlot(
  object = DN2_AIDA_Combined,
  var = "sex", 
  group.by = "seurat_clusters", 
  retain.factor.levels = T, 
  theme=theme_half_open(), main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold"))

#Statistical Significance
library("scProportionTest")

prop_test <- sc_utils(DN2_AIDA_Combined)

#unique(Batch1_3_Tumor_only_UPDATED$Treatment)

#prop test for cluster vs sex 
prop_test <- permutation_test(
  prop_test, cluster_identity = "seurat_clusters",
  sample_1 = "male", sample_2 = "female",
  sample_identity = "sex"
)
permutation_plot(prop_test)


#get packages ready
library(rsconnect)
library(Seurat)
library(ShinyCell)

#connect to account
rsconnect::setAccountInfo(name='hmiura1221',
                          token='93EAD3895AF84DA902AAD43FDAF30183',
                          secret='JhJweWwnaU166CRJNUbLcMacrHjXuRBeFkr0dz5r')

scConf = createConfig(DN2_AIDA_Combined)
write.csv(scConf,'scConf_meta_dn2_shiny.csv')
scConf= modColours(scConf,  meta.to.mod = 'sex',
                   new.colours = c('red3','blue3'))
makeShinyApp(DN2_AIDA_Combined, scConf, gene.mapping = TRUE,
             shiny.title = "DN2 B Cell Characterization for Males and Females", #title for overall page
             default.gene1 = "TBX21", default.gene2 = "NKG7", #default first two data you see (dn2, and cytotoxicity)
             shiny.dir = "shinyApp_DN2/") 

rsconnect::deployApp("shinyApp_DN2/") #done; only when ready!



# Markers to Show for DN2: Density Plot
markers <- list()
markers$DN2.A<- c('CLEC2D', 'EZR', 'CD79A', 'FCRL5', 'FCRL3', 'FGR', 'IGHD', 'IGHG3', 'LAPTM5')
markers$DN2.B<- c('TAGLN2', 'COTL1', 'ANXA2', 'ARPC1B', 'MARCKS', 'CAPG')
markers$DN2.Cytotoxic<- c('GNLY', 'PRF1', 'SRGN', 'GZMH', 'GZMB', 'NKG7')
markers$DN2.B.Male<-c('CXCR4')
markers$DN2.B.Female<-c('FGR', 'IGHA1', 'IGHA2', 'IGHG3')

DN2_AIDA_Combined <- AddModuleScore_UCell(DN2_AIDA_Combined,features = markers)

#Marker DN2.A-ACTIVATING CELL SURFACE RECEPTORS
DN2_AIDA_Combined$DN2.A_UCell
plot_density(
  DN2_AIDA_Combined,
  c('DN2.A_UCell'),
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

#Marker DN2.B- ACTIN 
DN2_AIDA_Combined$DN2.B_UCell
plot_density(
  DN2_AIDA_Combined,
  c('DN2.B_UCell'),
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

#Marker DN2.Cytotoxic- Cytolytic granule
DN2_AIDA_Combined$DN2.Cytotoxic_UCell
plot_density(
  DN2_AIDA_Combined,
  c('DN2.Cytotoxic_UCell'),
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

#Marker DN2.B.Male- Immunoglobulin Receptor binding
DN2_AIDA_Combined$DN2.B.Male_UCell
plot_density(
  DN2_AIDA_Combined,
  c('DN2.B.Male_UCell'),
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

#Marker DN2.B.Female- Immunoglobulin Receptor binding
DN2_AIDA_Combined$DN2.B.Female_UCell
plot_density(
  DN2_AIDA_Combined,
  c('DN2.B.Female_UCell'),
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



DN2_AIDA_Combined@meta.data$Study_ID<-'Powell_2022'
META<-DN2_AIDA_Combined[[]]

DN2_AIDA_Combined@meta.data$Health_Status<-'Healthy'

