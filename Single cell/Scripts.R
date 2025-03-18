# Save 
> saveRDS(pBAT.data, file='id.rds')
> save.image("myfile_image_D3_pBAT") 
> savehistory("myfile_history_D3_pBAT") 
> q() 

> RDS<-readRDS('/bin/id.rds') 
> RDS 
> library(dplyr) > library(Seurat) 
> setwd("/bin/test ") 
> data_dir <- "/bin/test" 
> library(dplyr) 
> library(Seurat) 

> data < Read10X(data.dir = test)
> dim(data) #check
> pBAT.data<- CreateSeuratObject(counts = pBAT.data, project = "pBAT.data", min.cells = 3, min features = 200)
> pBAT.data
> pBAT.data[["percent.mt"|] <- PercentageFeatureSet(BAT.data, pattern = "MT-")
> head(pBAT.data@meta.data, 5)

# plot
> p=VlnPlot(pBAT.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
> pdf("data.pdf") 
>p
> dev.off() 
> getwd() 

> p=plot1 <- FeatureScatter(pBAT.data, feature1 = "nCount_RNA", feature2 = "percent.mt") 
> plot2 <- FeatureScatter(pBAT.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
> p=CombinePlots(plots = list(plot1, plot2))
> pdf("qualitycontrol-after2-D3_2_pBAT.pdf")
> p
> dev.off()
> getwd() 

#data normalization
> N.data <- NormalizeData(data)
> N.data@assays$RNA[30:34,1:3]
> N.data <- FindVariableFeatures(N.data,selection.method = "vst",nfeatures = 2000) 
> top10 <- head(VariableFeatures(N.data), 10);top10 

> p=plot1 <- VariableFeaturePlot(N.data)
> pdf("N.data.pdf")
> p
> dev.off)

> p=plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
> pdf("N.data.pdf")
> p
> all.genes <- rownames(N.data)
> N.data <- ScaleData(N.data, features = all.genes)
> N.data[["RNA"I]@scale.data[30:34,1:3]
> N.data <- RunPCA(N.data, features = VariableFeatures(object = N.data))
> p=VizDimLoadings(N.data, dims = 1:2, reduction = "pca")
>pdf("VizDimLoadings-gene.pdf")
> p
> dev.off

#DimPlot(cds2, reduction = "umap",split.by = 'orig.ident' label=TRUE, ncol=2,group.by = "seurat_clusters")
> p= DimPlot(N.data, reduction = "pca"‚split.by = 'ident")
> pdf("VizDimLoadings-gene-2PCA-display.pdf")
> p
> dev.off()
> pdf("test.pdf")
> DimHeatmap(N.data, dims = 1:20, cells = 500, balanced = TRUE)
> pdf("DimHeatmap-20PCA-display.pdf")
> p
> dev.off()

> N.data <- JackStraw(N.data, num.replicate = 100) 
> N.data <- ScoreJackStraw(N.data, dims = 1:20) 
> p= JackStrawPlot(N.data, dims = 1:20)
> pdf("JackStrawPlot.pdf")
>p 

> p= ElbowPlot(N.data) 
> pdf("ElbowPlot.pdf") 
> p
> dev.off() 

> p= DimPlot(N.data, reduction = "pca") + NoLegend() 
> pdf("DimPlot.pdf")
>p
> dev.off()
