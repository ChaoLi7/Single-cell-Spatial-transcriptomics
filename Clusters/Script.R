# Data Visualization with Nonlinear Dimensionality Reduction Methods
> N.data <- FindNeighbors(BAT.data, dims = 1:20)
> N.data <- FindClusters (N.data, resolution = 0.5)
> N.data <- FindClustersN.data, resolution = 0.2)
> table(N.data@active.ident)
> table(N.data@active.ident)
> head(subset(as.data.frame(BAT.data@active.ident),N.data@active.ident=="2"))
> N.data <- RunTSNE(N.data, dims = 1:20)
> p=DimPlot(N.data,reduction = "tsne",pt.size = 1,label = TRUE, cols = c(#99c9fb', 'green3', '#FF7F00', 'pink', orchid', 're d','dodgerblue', 'yellow', 'grey60', '#FB9A99','black','#c51b8a'))
> pdf("D3-N_TSNE-0.5.pdf")
> p
> dev.off()

> N.data<-RunUMAP(N.data, dims = 1:20)
> p=DimPlot(BAT_UAMP, reduction = "umap", label = TRUE)
> pdf("D3-N_UAMP-0.5.pdf")
> p
> dev.off()

> N.data.makers <- FindAl|Markers (N.data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5<-N.data.makers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
> p=DoHeatmap(N.data, features = top5$gene) + NoLegend()
> pdf("D3-N_tsne_DoHeatmap-0.5-5gene.pdf")
> p
> dev.offl
>top100<-N.data.makers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_ log2FC)
> write.table(top100,file="min.pct0.1_log2fc0.25_markers.txt",quote=F,sep="It", row.names=F,col.names=T)
> N.data.makers <- FindAllMarkers (N.data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
> head(N.data.makers)
> N.data.makers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
> top5<-N.data.makers %>% group by(cluster) %>% top n(n = 5, wt =avg.1og2FC) > top5
> p=DoHeatmap(N.data, features = top5$gene) + NoLegend)
> pdf("D3-N_TSNE_DoHeatmap-0.5-5gene.pdf")
> p
> dev.off()

> top100<-N.data.makers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
> write.table(top100,file="top100_markers.txt", quote=F,sep="|t", row.names=F,col.names=T)

