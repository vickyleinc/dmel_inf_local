
ob.list <- list(wtF_seurat_variable, wtJ_seurat_variable, ninesixh_seurat_variable, adult_seurat_variable)
genes.use <- c()
for (i in 1:length(ob.list)) {
  genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 1000))
}
genes.use <- names(which(table(genes.use) > 1))
for (i in 1:length(ob.list)) {
  genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
}

# Run multi-set CCA
combined <- RunMultiCCA(ob.list, genes.use = genes.use, num.ccs = 15)

#adding metadata for MetageneBicorPlot to separate 4 groups
combined@meta.data$identity <- rownames(combined@meta.data)
head(combined@meta.data)
combined@meta.data$identity <- substr(combined@meta.data$identity, 1, nchar(combined@meta.data$identity)-17) 
combined@meta.data$identity <- substr(combined@meta.data$identity,6, nchar(combined@meta.data$identity)) 

#Determine how many CC dimensions to use
MetageneBicorPlot(combined, grouping.var = "identity", dims.eval = 1:15)

tail(combined@meta.data)

##CCA alignment

combined_integrated <- AlignSubspace(combined,
                                     reduction.type = "cca",
                                     grouping.var = "orig.ident",
                                     dims.align = 1:10)

# t-SNE and Clustering
combined_clusters <- FindClusters(combined_integrated, reduction.type = "cca.aligned",
                                    dims.use = 1:10, save.SNN = T, resolution = 0.8)
combined_tsne <- RunTSNE(combined_clusters,
                               reduction.use = "cca.aligned",
                               dims.use = 1:10)

# Visualization
TSNEPlot(combined_tsne, do.label = T)

FeaturePlot(object = combined_tsne, features.plot = c("vas", "hh", "drm", "tj", "bam", "Trim9", "Dh44-R2", "Fas3", "CG43693"), min.cutoff = "q9", cols.use = c("lightgrey","blue"), pt.size = 0.5)



