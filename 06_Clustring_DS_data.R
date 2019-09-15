library(scater)
library(ggplot2)
library(biomaRt)
library(tibble)
library(data.table)
library(kBET)
library(cluster)

# Gene ID Mapping Biomart ====
# Functions ####
load('/Volumes/Ati-Archive/HCA/SC_Protocols/Version3/R_analysis/Biomart_mmus_mapping_table.RData')
load('/Volumes/Ati-Archive/HCA/SC_Protocols/Version3/R_analysis/Biomart_hsap_mapping_table.RData')
mapIDs <- function (count.mat, species){
  ens_ids <- rownames(count.mat)
  mat.rownames <- as.character(lapply(ens_ids, function (x) unlist(strsplit(x, split = ".", fixed = T))[1]))
  mat.dup.rows <- which(duplicated(mat.rownames))
  if (length(mat.dup.rows) != 0){
    mat.rownames <- mat.rownames[-mat.dup.rows]
    count.mat <- count.mat[-mat.dup.rows,]
  }
  if (species == "hsap"){
    hsap.genes.with.ids <- intersect(mat.rownames, rownames(hsap.GID.mapping.final))
    hsap.rownames.mapped <- mat.rownames
    names(hsap.rownames.mapped) <- hsap.rownames.mapped
    hsap.rownames.mapped[hsap.genes.with.ids] <- hsap.GID.mapping.final[hsap.genes.with.ids, "hgnc_symbol"]
    rownames(count.mat) <- hsap.rownames.mapped
  }
  else if (species == "mmus"){
    mmus.genes.with.ids <- intersect(mat.rownames, rownames(mmus.GID.mapping.final))
    mmus.rownames.mapped <- mat.rownames
    names(mmus.rownames.mapped) <- mmus.rownames.mapped
    mmus.rownames.mapped[mmus.genes.with.ids] <- mmus.GID.mapping.final[mmus.genes.with.ids, "mgi_symbol"]
    rownames(count.mat) <- mmus.rownames.mapped
  }
  return(count.mat)
  
}

load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/col.clusters.RData")
col


# CELseq2 Clustering ####
# loading Annotated seurat objects #
load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/celseq_seu.obj.RData")
CELseq2.hsap.obj <- celseq
rm(celseq)
CELseq2.hsap.obj@ident <- CELseq2.hsap.obj@meta.data$nnet2
names(CELseq2.hsap.obj@ident) <- rownames(CELseq2.hsap.obj@meta.data)
CELseq2.hsap.metadata <- CELseq2.hsap.obj@meta.data

load("/Volumes/Ati-Archive/HCA/SCE_Robjects/Downsampled/DS_secondRound/Separated_DS_datasets/CELseq2_DS_20K.RData")
CELseq2.20K <- mapIDs(CELseq2.20K, "hsap")
colnames(CELseq2.20K) <- gsub("\\.", "_", colnames(CELseq2.20K))
common.cells <- intersect(colnames(CELseq2.20K), colnames(CELseq2.hsap.obj@scale.data))
common.genes <- intersect(rownames(CELseq2.20K), rownames(CELseq2.hsap.obj@scale.data))
CELseq2.20K.clean <- CELseq2.20K[common.genes, common.cells]

CELseq2.DS.obj <- CreateSeuratObject(raw.data = as.matrix(CELseq2.20K.clean), project = "CELseq2")
CELseq2.DS.obj <- AddMetaData(object = CELseq2.DS.obj, metadata = CELseq2.hsap.obj@ident, col.name = "clusters")

CELseq2.DS.obj <- NormalizeData(object = CELseq2.DS.obj, normalization.method = "LogNormalize", scale.factor = 10000)
dim(CELseq2.DS.obj@data)
CELseq2.DS.obj <- FindVariableGenes(object = CELseq2.DS.obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.2, x.high.cutoff = 5, y.cutoff = 0.49)
length(CELseq2.DS.obj@var.genes)
CELseq2.DS.obj <- ScaleData(object = CELseq2.DS.obj)#, vars.to.regress = c("nUMI", "percent.mito"))
CELseq2.DS.obj <- RunPCA(object = CELseq2.DS.obj, pc.genes = CELseq2.DS.obj@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = CELseq2.DS.obj, dim.1 = 1, dim.2 = 2)
#CELseq2.DS.obj <- JackStraw(object = CELseq2.DS.obj, num.replicate = 100, display.progress = FALSE)
#JackStrawPlot(object = CELseq2.DS.obj, PCs = 1:12)
PCElbowPlot(object = CELseq2.DS.obj)
CELseq2.DS.obj <- FindClusters(object = CELseq2.DS.obj, reduction.type = "pca", dims.use = 1:8, resolution = 0.6, print.output = 0, save.SNN = TRUE)
CELseq2.DS.obj <- RunTSNE(object = CELseq2.DS.obj, dims.use = 1:8, do.fast = TRUE)

CELseq2.colors <- c("blueviolet","aquamarine","green4","maroon","orange","black","red","deepskyblue3", "grey")
#pdf("CELseq2_clustering_DS20K_clusters.pdf")
png(file="CELseq2_clustering_DS20K_clusters.png", width = 8, height = 6, units= "in",res = 600)
TSNEPlot(object = CELseq2.DS.obj, group.by = "clusters", colors.use = CELseq2.colors, pt.size = 2.5)
dev.off()
save(CELseq2.DS.obj, file="CELseq2_clustering_DS20K_dim1to8_res06_OBJ.RData")
png(file="/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/elbowPlots/CELseq2_clustering_DS20K_clusters_elbowplot.png", width = 6, height = 6, units= "in",res = 600)
PCElbowPlot(object = CELseq2.DS.obj)
dev.off()

library(scran)
library(igraph)
library(cluster)
CELseq2.graph <- graph.adjacency(CELseq2.DS.obj@snn, weighted = T, mode="undirected")
CELseq2.clusters <- CELseq2.DS.obj@ident
CELseq2.graph.mod.score <- modularity(CELseq2.graph, CELseq2.clusters)
CELseq2.cluster.mod.score <- clusterModularity(CELseq2.graph, CELseq2.clusters)#, get.values=FALSE) 
rm(CELseq2.hsap.obj)
rm(CELseq2.20K)
png(file="CELseq2_clustering_DS20K_graph.png", width = 8, height = 6, units= "in",res = 600)
plot.igraph(CELseq2.graph, layout= as.matrix(CELseq2.DS.obj@dr$tsne@cell.embeddings), edge.width = E(graph = CELseq2.graph)$weight, vertex.label = NA,vertex.size = 0 )
dev.off()
#Silhouette
CELseq2.PCA.data <- CELseq2.DS.obj@dr$pca@cell.embeddings[,1:8]
CELseq2.dd <- as.matrix(dist(CELseq2.PCA.data))
CELseq2.sil <- summary(silhouette(as.numeric(CELseq2.clusters), CELseq2.dd))$avg.width



# MARSseq Clustering ####
# loading Annotated seurat objects #
load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/marsseq_seu.obj.RData")
MARSseq.hsap.obj <- marsseq
rm(marsseq)
MARSseq.hsap.metadata <- MARSseq.hsap.obj@meta.data
MARSseq.hsap.obj@ident <- MARSseq.hsap.obj@meta.data$nnet2
names(MARSseq.hsap.obj@ident) <- rownames(MARSseq.hsap.obj@meta.data)

load("/Volumes/Ati-Archive/HCA/SCE_Robjects/Downsampled/DS_secondRound/Separated_DS_datasets/MARSseq_DS_20K.RData")
MARSseq.20K <- mapIDs(MARSseq.20K, "hsap")
colnames(MARSseq.20K) <- gsub("\\.", "_", colnames(MARSseq.20K))
common.cells <- intersect(colnames(MARSseq.20K), colnames(MARSseq.hsap.obj@scale.data))
common.genes <- intersect(rownames(MARSseq.20K), rownames(MARSseq.hsap.obj@scale.data))
MARSseq.20K.clean <- MARSseq.20K[common.genes, common.cells]

MARSseq.DS.obj <- CreateSeuratObject(raw.data = as.matrix(MARSseq.20K.clean), project = "MARSseq")
MARSseq.DS.obj <- AddMetaData(object = MARSseq.DS.obj, metadata = MARSseq.hsap.obj@ident, col.name = "clusters")

MARSseq.DS.obj <- NormalizeData(object = MARSseq.DS.obj, normalization.method = "LogNormalize", scale.factor = 10000)
dim(MARSseq.DS.obj@data)
MARSseq.DS.obj <- FindVariableGenes(object = MARSseq.DS.obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.25, x.high.cutoff = 5, y.cutoff = 0.5)
length(MARSseq.DS.obj@var.genes)
MARSseq.DS.obj <- ScaleData(object = MARSseq.DS.obj)#, vars.to.regress = c("nUMI", "percent.mito"))
MARSseq.DS.obj <- RunPCA(object = MARSseq.DS.obj, pc.genes = MARSseq.DS.obj@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = MARSseq.DS.obj, dim.1 = 1, dim.2 = 2)
#MARSseq.DS.obj <- JackStraw(object = MARSseq.DS.obj, num.replicate = 100, display.progress = FALSE)
#JackStrawPlot(object = MARSseq.DS.obj, PCs = 1:12)
PCElbowPlot(object = MARSseq.DS.obj)
MARSseq.DS.obj <- FindClusters(object = MARSseq.DS.obj, reduction.type = "pca", dims.use = 1:8, resolution = 0.6, print.output = 0, save.SNN = TRUE)
MARSseq.DS.obj <- RunTSNE(object = MARSseq.DS.obj, dims.use = 1:8, do.fast = TRUE)

MARSseq.colors <- c("blueviolet","aquamarine","green4","maroon","orange","black","red","coral2","deepskyblue3", "grey")
png("MARSseq_clustering_DS20K_clusters.png", width = 8, height = 6, units= "in",res = 600)
TSNEPlot(object = MARSseq.DS.obj, group.by = "clusters", colors.use = MARSseq.colors, pt.size = 2.5)
dev.off()
save(MARSseq.DS.obj, file="MARSseq_clustering_DS20K_dim1to8_res06_OBJ.RData")
png(file="/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/elbowPlots/MARSseq_clustering_DS20K_clusters_elbowplot.png", width = 6, height = 6, units= "in",res = 600)
PCElbowPlot(object = MARSseq.DS.obj)
dev.off()

library(scran)
library(igraph)
MARSseq.graph <- graph.adjacency(MARSseq.DS.obj@snn, weighted = T, mode="undirected")
MARSseq.clusters <- MARSseq.DS.obj@ident
MARSseq.graph.mod.score <- modularity(MARSseq.graph, MARSseq.clusters)
MARSseq.cluster.mod.score <- clusterModularity(MARSseq.graph, MARSseq.clusters)#, get.values=FALSE) 
rm(MARSseq.hsap.obj)
png(file="MARSseq_clustering_DS20K_graph.png", width = 8, height = 6, units= "in",res = 600)
plot.igraph(MARSseq.graph, layout= as.matrix(MARSseq.DS.obj@dr$tsne@cell.embeddings), edge.width = E(graph = MARSseq.graph)$weight, vertex.label = NA,vertex.size = 0 )
dev.off()

#Silhouette
MARSseq.PCA.data <- MARSseq.DS.obj@dr$pca@cell.embeddings[,1:8]
MARSseq.dd <- as.matrix(dist(MARSseq.PCA.data))
MARSseq.sil <- summary(silhouette(as.numeric(MARSseq.clusters), MARSseq.dd))$avg.width


# QUARTZseq Clustering ####
# loading Annotated seurat objects #
load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/quartzseq_seu.obj.RData")
QUARTZseq.hsap.obj <- quartzseq
rm(quartzseq)
QUARTZseq.hsap.obj@ident <- QUARTZseq.hsap.obj@meta.data$nnet2
names(QUARTZseq.hsap.obj@ident) <- rownames(QUARTZseq.hsap.obj@meta.data)
QUARTZseq.hsap.metadata <- QUARTZseq.hsap.obj@meta.data

load("/Volumes/Ati-Archive/HCA/SCE_Robjects/Downsampled/DS_secondRound/Separated_DS_datasets/QUARTZseq_DS_20K.RData")
QUARTZseq.20K <- mapIDs(QUARTZseq.20K, "hsap")
colnames(QUARTZseq.20K) <- gsub("\\.", "_", colnames(QUARTZseq.20K))
common.cells <- intersect(colnames(QUARTZseq.20K), colnames(QUARTZseq.hsap.obj@scale.data))
common.genes <- intersect(rownames(QUARTZseq.20K), rownames(QUARTZseq.hsap.obj@scale.data))
QUARTZseq.20K.clean <- QUARTZseq.20K[common.genes, common.cells]

QUARTZseq.DS.obj <- CreateSeuratObject(raw.data = as.matrix(QUARTZseq.20K.clean), project = "QUARTZseq")
QUARTZseq.DS.obj <- AddMetaData(object = QUARTZseq.DS.obj, metadata = QUARTZseq.hsap.obj@ident, col.name = "clusters")

QUARTZseq.DS.obj <- NormalizeData(object = QUARTZseq.DS.obj, normalization.method = "LogNormalize", scale.factor = 10000)
dim(QUARTZseq.DS.obj@data)
QUARTZseq.DS.obj <- FindVariableGenes(object = QUARTZseq.DS.obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.22, x.high.cutoff = 5, y.cutoff = 0.5)
length(QUARTZseq.DS.obj@var.genes)
QUARTZseq.DS.obj <- ScaleData(object = QUARTZseq.DS.obj)#, vars.to.regress = c("nUMI", "percent.mito"))
QUARTZseq.DS.obj <- RunPCA(object = QUARTZseq.DS.obj, pc.genes = QUARTZseq.DS.obj@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = QUARTZseq.DS.obj, dim.1 = 1, dim.2 = 2)
#QUARTZseq.DS.obj <- JackStraw(object = QUARTZseq.DS.obj, num.replicate = 100, display.progress = FALSE)
#JackStrawPlot(object = QUARTZseq.DS.obj, PCs = 1:12)
PCElbowPlot(object = QUARTZseq.DS.obj)
QUARTZseq.DS.obj <- FindClusters(object = QUARTZseq.DS.obj, reduction.type = "pca", dims.use = 1:8, resolution = 0.6, print.output = 0, save.SNN = TRUE)
QUARTZseq.DS.obj <- RunTSNE(object = QUARTZseq.DS.obj, dims.use = 1:8, do.fast = TRUE)

QUARTZseq.colors <- c("blueviolet","aquamarine","green4","maroon","orange","black","red","coral2","deepskyblue3", "grey")
png("QUARTZseq_clustering_DS20K_clusters.png", width = 8, height = 6, units= "in",res = 600)
TSNEPlot(object = QUARTZseq.DS.obj, group.by = "clusters", colors.use = QUARTZseq.colors, pt.size = 2.5)
dev.off()
save(QUARTZseq.DS.obj, file="QUARTZseq_clustering_DS20K_dim1to8_res06_OBJ.RData")
png(file="/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/elbowPlots/QUARTZseq_clustering_DS20K_clusters_elbowplot.png", width = 6, height = 6, units= "in",res = 600)
PCElbowPlot(object = QUARTZseq.DS.obj)
dev.off()


library(scran)
library(igraph)
QUARTZseq.graph <- graph.adjacency(QUARTZseq.DS.obj@snn, weighted = T, mode="undirected")
QUARTZseq.clusters <- QUARTZseq.DS.obj@ident
QUARTZseq.graph.mod.score <- modularity(QUARTZseq.graph, QUARTZseq.clusters)
QUARTZseq.cluster.mod.score <- clusterModularity(QUARTZseq.graph, QUARTZseq.clusters)#, get.values=FALSE) 
rm(QUARTZseq.hsap.obj)
rm(QUARTZseq.20K)
png(file="QUARTZseq_clustering_DS20K_graph.png", width = 8, height = 6, units= "in",res = 600)
plot.igraph(QUARTZseq.graph, layout= as.matrix(QUARTZseq.DS.obj@dr$tsne@cell.embeddings), edge.width = E(graph = QUARTZseq.graph)$weight, vertex.label = NA,vertex.size = 0 )
dev.off()

#Silhouette
QUARTZseq.PCA.data <- QUARTZseq.DS.obj@dr$pca@cell.embeddings[,1:8]
QUARTZseq.dd <- as.matrix(dist(QUARTZseq.PCA.data))
QUARTZseq.sil <- summary(silhouette(as.numeric(QUARTZseq.clusters), QUARTZseq.dd))$avg.width


# SCRBseq Clustering ####
# loading Annotated seurat objects #
load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/scrbseq_seu.obj.RData")
SCRBseq.hsap.obj <- scrbseq
rm(scrbseq)
SCRBseq.hsap.obj@ident <- SCRBseq.hsap.obj@meta.data$nnet2
names(SCRBseq.hsap.obj@ident) <- rownames(SCRBseq.hsap.obj@meta.data)
SCRBseq.hsap.metadata <- SCRBseq.hsap.obj@meta.data

load("/Volumes/Ati-Archive/HCA/SCE_Robjects/Downsampled/DS_secondRound/Separated_DS_datasets/SCRBseq_DS_20K.RData")
SCRBseq.20K <- mapIDs(SCRBseq.20K, "hsap")
colnames(SCRBseq.20K) <- gsub("\\.", "_", colnames(SCRBseq.20K))
common.cells <- intersect(colnames(SCRBseq.20K), colnames(SCRBseq.hsap.obj@scale.data))
common.genes <- intersect(rownames(SCRBseq.20K), rownames(SCRBseq.hsap.obj@scale.data))
SCRBseq.20K.clean <- SCRBseq.20K[common.genes, common.cells]

SCRBseq.DS.obj <- CreateSeuratObject(raw.data = as.matrix(SCRBseq.20K.clean), project = "SCRBseq")
SCRBseq.DS.obj <- AddMetaData(object = SCRBseq.DS.obj, metadata = SCRBseq.hsap.obj@ident, col.name = "clusters")

SCRBseq.DS.obj <- NormalizeData(object = SCRBseq.DS.obj, normalization.method = "LogNormalize", scale.factor = 10000)
dim(SCRBseq.DS.obj@data)
SCRBseq.DS.obj <- FindVariableGenes(object = SCRBseq.DS.obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.3, x.high.cutoff = 5, y.cutoff = 0.5)
length(SCRBseq.DS.obj@var.genes)
SCRBseq.DS.obj <- ScaleData(object = SCRBseq.DS.obj)#, vars.to.regress = c("nUMI", "percent.mito"))
SCRBseq.DS.obj <- RunPCA(object = SCRBseq.DS.obj, pc.genes = SCRBseq.DS.obj@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = SCRBseq.DS.obj, dim.1 = 1, dim.2 = 2)
#SCRBseq.DS.obj <- JackStraw(object = SCRBseq.DS.obj, num.replicate = 100, display.progress = FALSE)
#JackStrawPlot(object = SCRBseq.DS.obj, PCs = 1:12)
PCElbowPlot(object = SCRBseq.DS.obj)
SCRBseq.DS.obj <- FindClusters(object = SCRBseq.DS.obj, reduction.type = "pca", dims.use = 1:8, resolution = 0.6, print.output = 0, save.SNN = TRUE)
SCRBseq.DS.obj <- RunTSNE(object = SCRBseq.DS.obj, dims.use = 1:8, do.fast = TRUE)

SCRBseq.colors <- c("blueviolet","aquamarine","green4","maroon","orange","black","red","coral2","deepskyblue3", "grey")
png("SCRBseq_clustering_DS20K_clusters.png", width = 8, height = 6, units= "in",res = 600)
TSNEPlot(object = SCRBseq.DS.obj, group.by = "clusters", colors.use = SCRBseq.colors, pt.size = 2.5)
dev.off()
save(SCRBseq.DS.obj, file="SCRBseq_clustering_DS20K_dim1to8_res06_OBJ.RData")
png(file="/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/elbowPlots/SCRBseq_clustering_DS20K_clusters_elbowplot.png", width = 6, height = 6, units= "in",res = 600)
PCElbowPlot(object = SCRBseq.DS.obj)
dev.off()


library(scran)
library(igraph)
SCRBseq.graph <- graph.adjacency(SCRBseq.DS.obj@snn, weighted = T, mode="undirected")
SCRBseq.clusters <- SCRBseq.DS.obj@ident
SCRBseq.graph.mod.score <- modularity(SCRBseq.graph, SCRBseq.clusters)
SCRBseq.cluster.mod.score <- clusterModularity(SCRBseq.graph, SCRBseq.clusters)#, get.values=FALSE) 
rm(SCRBseq.hsap.obj)
rm(SCRBseq.20K)
png(file="SCRBseq_clustering_DS20K_graph.png", width = 8, height = 6, units= "in",res = 600)
plot.igraph(SCRBseq.graph, layout= as.matrix(SCRBseq.DS.obj@dr$tsne@cell.embeddings), edge.width = E(graph = SCRBseq.graph)$weight, vertex.label = NA,vertex.size = 0 )
dev.off()

#Silhouette
SCRBseq.PCA.data <- SCRBseq.DS.obj@dr$pca@cell.embeddings[,1:8]
SCRBseq.dd <- as.matrix(dist(SCRBseq.PCA.data))
SCRBseq.sil <- summary(silhouette(as.numeric(SCRBseq.clusters), SCRBseq.dd))$avg.width




# SMARTseq2 Clustering ####
# loading Annotated seurat objects #
load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/smartseq_seu.obj.RData")
SMARTseq2.hsap.obj <- smartseq
rm(smartseq)
SMARTseq2.hsap.obj@ident <- SMARTseq2.hsap.obj@meta.data$nnet2
names(SMARTseq2.hsap.obj@ident) <- rownames(SMARTseq2.hsap.obj@meta.data)
SMARTseq2.hsap.metadata <- SMARTseq2.hsap.obj@meta.data

load("/Volumes/Ati-Archive/HCA/SCE_Robjects/Downsampled/DS_secondRound/Human/Separated_DS_hsap_datasets/SMARTseq2_DS_hsap_20K.RData")
SMARTseq2.20K <- mapIDs(SMARTseq2.20K, "hsap")
colnames(SMARTseq2.20K) <- gsub("\\.", "_", colnames(SMARTseq2.20K))
common.cells <- intersect(colnames(SMARTseq2.20K), colnames(SMARTseq2.hsap.obj@scale.data))
common.genes <- intersect(rownames(SMARTseq2.20K), rownames(SMARTseq2.hsap.obj@scale.data))
SMARTseq2.20K.clean <- SMARTseq2.20K[common.genes, common.cells]

SMARTseq2.DS.obj <- CreateSeuratObject(raw.data = as.matrix(SMARTseq2.20K.clean), project = "SMARTseq2")
SMARTseq2.DS.obj <- AddMetaData(object = SMARTseq2.DS.obj, metadata = SMARTseq2.hsap.obj@ident, col.name = "clusters")

SMARTseq2.DS.obj <- NormalizeData(object = SMARTseq2.DS.obj, normalization.method = "LogNormalize", scale.factor = 10000)
dim(SMARTseq2.DS.obj@data)
SMARTseq2.DS.obj <- FindVariableGenes(object = SMARTseq2.DS.obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.3, x.high.cutoff = 5, y.cutoff = 0.55)
length(SMARTseq2.DS.obj@var.genes)
SMARTseq2.DS.obj <- ScaleData(object = SMARTseq2.DS.obj)#, vars.to.regress = c("nUMI", "percent.mito"))
SMARTseq2.DS.obj <- RunPCA(object = SMARTseq2.DS.obj, pc.genes = SMARTseq2.DS.obj@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = SMARTseq2.DS.obj, dim.1 = 1, dim.2 = 2)
#SMARTseq2.DS.obj <- JackStraw(object = SMARTseq2.DS.obj, num.replicate = 100, display.progress = FALSE)
#JackStrawPlot(object = SMARTseq2.DS.obj, PCs = 1:12)
PCElbowPlot(object = SMARTseq2.DS.obj)
SMARTseq2.DS.obj <- FindClusters(object = SMARTseq2.DS.obj, reduction.type = "pca", dims.use = 1:8, resolution = 0.6, print.output = 0, save.SNN = TRUE)
SMARTseq2.DS.obj <- RunTSNE(object = SMARTseq2.DS.obj, dims.use = 1:8, do.fast = TRUE)

SMARTseq2.colors <- c("blueviolet","aquamarine","green4","maroon","orange","black","red","deepskyblue3", "grey")
png("SMARTseq2_clustering_DS20K_clusters.png", width = 8, height = 6, units= "in",res = 600)
TSNEPlot(object = SMARTseq2.DS.obj, group.by = "clusters", colors.use = SMARTseq2.colors, pt.size = 2.5)
dev.off()
save(SMARTseq2.DS.obj, file="SMARTseq2_clustering_DS20K_dim1to8_res06_OBJ.RData")
png(file="/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/elbowPlots/SMARTseq2_clustering_DS20K_clusters_elbowplot.png", width = 6, height = 6, units= "in",res = 600)
PCElbowPlot(object = SMARTseq2.DS.obj)
dev.off()


library(scran)
library(igraph)
SMARTseq2.graph <- graph.adjacency(SMARTseq2.DS.obj@snn, weighted = T, mode="undirected")
SMARTseq2.clusters <- SMARTseq2.DS.obj@ident
SMARTseq2.graph.mod.score <- modularity(SMARTseq2.graph, SMARTseq2.clusters)
SMARTseq2.cluster.mod.score <- clusterModularity(SMARTseq2.graph, SMARTseq2.clusters)#, get.values=FALSE) 
rm(SMARTseq2.hsap.obj)
rm(SMARTseq2.DS.UMI)
png(file="SMARTseq2_clustering_DS20K_graph.png", width = 8, height = 6, units= "in",res = 600)
plot.igraph(SMARTseq2.graph, layout= as.matrix(SMARTseq2.DS.obj@dr$tsne@cell.embeddings), edge.width = E(graph = SMARTseq2.graph)$weight, vertex.label = NA,vertex.size = 0 )
dev.off()

#Silhouette
SMARTseq2.PCA.data <- SMARTseq2.DS.obj@dr$pca@cell.embeddings[,1:8]
SMARTseq2.dd <- as.matrix(dist(SMARTseq2.PCA.data))
SMARTseq2.sil <- summary(silhouette(as.numeric(SMARTseq2.clusters), SMARTseq2.dd))$avg.width


# C1HTsmall Clustering ####
# loading Annotated seurat objects #
load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/c1ht.s_seu.obj.RData")
C1HTsmall.hsap.obj <- c1ht.s
rm(c1ht.s)
C1HTsmall.hsap.obj@ident <- C1HTsmall.hsap.obj@meta.data$nnet2
names(C1HTsmall.hsap.obj@ident) <- rownames(C1HTsmall.hsap.obj@meta.data)
C1HTsmall.hsap.metadata <- C1HTsmall.hsap.obj@meta.data

load("/Volumes/Ati-Archive/HCA/SCE_Robjects/Downsampled/DS_secondRound/Separated_DS_datasets/C1HTsmall_DS_20K.RData")
C1HTsmall.20K <- mapIDs(C1HTsmall.20K, "hsap")
colnames(C1HTsmall.20K) <- gsub("\\.", "_", colnames(C1HTsmall.20K))
common.cells <- intersect(colnames(C1HTsmall.20K), colnames(C1HTsmall.hsap.obj@scale.data))
common.genes <- intersect(rownames(C1HTsmall.20K), rownames(C1HTsmall.hsap.obj@scale.data))
C1HTsmall.20K.clean <- C1HTsmall.20K[common.genes, common.cells]

C1HTsmall.DS.obj <- CreateSeuratObject(raw.data = as.matrix(C1HTsmall.20K.clean), project = "C1HTsmall")
C1HTsmall.DS.obj <- AddMetaData(object = C1HTsmall.DS.obj, metadata = C1HTsmall.hsap.obj@ident, col.name = "clusters")

C1HTsmall.DS.obj <- NormalizeData(object = C1HTsmall.DS.obj, normalization.method = "LogNormalize", scale.factor = 10000)
dim(C1HTsmall.DS.obj@data)
C1HTsmall.DS.obj <- FindVariableGenes(object = C1HTsmall.DS.obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.24, x.high.cutoff = 5, y.cutoff = 0.5)
length(C1HTsmall.DS.obj@var.genes)
C1HTsmall.DS.obj <- ScaleData(object = C1HTsmall.DS.obj)#, vars.to.regress = c("nUMI", "percent.mito"))
C1HTsmall.DS.obj <- RunPCA(object = C1HTsmall.DS.obj, pc.genes = C1HTsmall.DS.obj@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = C1HTsmall.DS.obj, dim.1 = 1, dim.2 = 2)
#C1HTsmall.DS.obj <- JackStraw(object = C1HTsmall.DS.obj, num.replicate = 100, display.progress = FALSE)
#JackStrawPlot(object = C1HTsmall.DS.obj, PCs = 1:12)
PCElbowPlot(object = C1HTsmall.DS.obj)
C1HTsmall.DS.obj <- FindClusters(object = C1HTsmall.DS.obj, reduction.type = "pca", dims.use = 1:8, resolution = 0.6, print.output = 0, save.SNN = TRUE)
C1HTsmall.DS.obj <- RunTSNE(object = C1HTsmall.DS.obj, dims.use = 1:8, do.fast = TRUE)

C1HTsmall.colors <- c("blueviolet","aquamarine","green4","maroon","orange","black","red","coral2","deepskyblue3", "grey")
png("C1HTsmall_clustering_DS20K_clusters.png", width = 8, height = 6, units= "in",res = 600)
TSNEPlot(object = C1HTsmall.DS.obj, group.by = "clusters", colors.use = C1HTsmall.colors, pt.size = 2.5)
dev.off()
save(C1HTsmall.DS.obj, file="C1HTsmall_clustering_DS20K_dim1to8_res06_OBJ.RData")
png(file="/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/elbowPlots/C1HTsmall_clustering_DS20K_clusters_elbowplot.png", width = 6, height = 6, units= "in",res = 600)
PCElbowPlot(object = C1HTsmall.DS.obj)
dev.off()

library(scran)
library(igraph)
C1HTsmall.graph <- graph.adjacency(C1HTsmall.DS.obj@snn, weighted = T, mode="undirected")
C1HTsmall.clusters <- C1HTsmall.DS.obj@ident
C1HTsmall.graph.mod.score <- modularity(C1HTsmall.graph, C1HTsmall.clusters)
C1HTsmall.cluster.mod.score <- clusterModularity(C1HTsmall.graph, C1HTsmall.clusters)#, get.values=FALSE) 
rm(C1HTsmall.hsap.obj)
rm(C1HTsmall.20K)
png(file="C1HTsmall_clustering_DS20K_graph.png", width = 8, height = 6, units= "in",res = 600)
plot.igraph(C1HTsmall.graph, layout= as.matrix(C1HTsmall.DS.obj@dr$tsne@cell.embeddings), edge.width = E(graph = C1HTsmall.graph)$weight, vertex.label = NA,vertex.size = 0 )
dev.off()

#Silhouette
C1HTsmall.PCA.data <- C1HTsmall.DS.obj@dr$pca@cell.embeddings[,1:8]
C1HTsmall.dd <- as.matrix(dist(C1HTsmall.PCA.data))
C1HTsmall.sil <- summary(silhouette(as.numeric(C1HTsmall.clusters), C1HTsmall.dd))$avg.width


# C1HTmedium Clustering ####
# loading Annotated seurat objects #
load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/c1ht.m_seu.obj.RData")
C1HTmedium.hsap.obj <- c1ht.m
rm(c1ht.m)
C1HTmedium.hsap.obj@ident <- C1HTmedium.hsap.obj@meta.data$nnet2
names(C1HTmedium.hsap.obj@ident) <- rownames(C1HTmedium.hsap.obj@meta.data)
C1HTmedium.hsap.metadata <- C1HTmedium.hsap.obj@meta.data

load("/Volumes/Ati-Archive/HCA/SCE_Robjects/Downsampled/DS_secondRound/Separated_DS_datasets/C1HTmedium_DS_20K.RData")
C1HTmedium.20K <- mapIDs(C1HTmedium.20K, "hsap")
colnames(C1HTmedium.20K) <- gsub("\\.", "_", colnames(C1HTmedium.20K))
common.cells <- intersect(colnames(C1HTmedium.20K), colnames(C1HTmedium.hsap.obj@scale.data))
common.genes <- intersect(rownames(C1HTmedium.20K), rownames(C1HTmedium.hsap.obj@scale.data))
C1HTmedium.20K.clean <- C1HTmedium.20K[common.genes, common.cells]

C1HTmedium.DS.obj <- CreateSeuratObject(raw.data = as.matrix(C1HTmedium.20K.clean), project = "C1HTmedium")
C1HTmedium.DS.obj <- AddMetaData(object = C1HTmedium.DS.obj, metadata = C1HTmedium.hsap.obj@ident, col.name = "clusters")

C1HTmedium.DS.obj <- NormalizeData(object = C1HTmedium.DS.obj, normalization.method = "LogNormalize", scale.factor = 10000)
dim(C1HTmedium.DS.obj@data)
C1HTmedium.DS.obj <- FindVariableGenes(object = C1HTmedium.DS.obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.28, x.high.cutoff = 5, y.cutoff = 0.5)
length(C1HTmedium.DS.obj@var.genes)
C1HTmedium.DS.obj <- ScaleData(object = C1HTmedium.DS.obj)#, vars.to.regress = c("nUMI", "percent.mito"))
C1HTmedium.DS.obj <- RunPCA(object = C1HTmedium.DS.obj, pc.genes = C1HTmedium.DS.obj@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = C1HTmedium.DS.obj, dim.1 = 1, dim.2 = 2)
#C1HTmedium.DS.obj <- JackStraw(object = C1HTmedium.DS.obj, num.replicate = 100, display.progress = FALSE)
#JackStrawPlot(object = C1HTmedium.DS.obj, PCs = 1:12)
PCElbowPlot(object = C1HTmedium.DS.obj)
C1HTmedium.DS.obj <- FindClusters(object = C1HTmedium.DS.obj, reduction.type = "pca", dims.use = 1:8, resolution = 0.6, print.output = 0, save.SNN = TRUE)
C1HTmedium.DS.obj <- RunTSNE(object = C1HTmedium.DS.obj, dims.use = 1:8, do.fast = TRUE)

C1HTmedium.colors <- c("blueviolet","aquamarine","green4","maroon","orange","black","red","coral2","deepskyblue3", "grey")
png("C1HTmedium_clustering_DS20K_clusters.png", width = 8, height = 6, units= "in",res = 600)
TSNEPlot(object = C1HTmedium.DS.obj, group.by = "clusters", colors.use = C1HTmedium.colors, pt.size = 2.5)
dev.off()
save(C1HTmedium.DS.obj, file="C1HTmedium_clustering_DS20K_dim1to8_res06_OBJ.RData")
length(which(C1HTsmall.DS.obj@meta.data$clusters == "Megakaryocytes"))
png(file="/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/elbowPlots/C1HTmedium_clustering_DS20K_clusters_elbowplot.png", width = 6, height = 6, units= "in",res = 600)
PCElbowPlot(object = C1HTmedium.DS.obj)
dev.off()

library(scran)
library(igraph)
C1HTmedium.graph <- graph.adjacency(C1HTmedium.DS.obj@snn, weighted = T, mode="undirected")
C1HTmedium.clusters <- C1HTmedium.DS.obj@ident
C1HTmedium.graph.mod.score <- modularity(C1HTmedium.graph, C1HTmedium.clusters)
C1HTmedium.cluster.mod.score <- clusterModularity(C1HTmedium.graph, C1HTmedium.clusters)#, get.values=FALSE) 
rm(C1HTmedium.hsap.obj)
rm(C1HTmedium.20K)

png(file="C1HTmedium_clustering_DS20K_graph.png", width = 8, height = 6, units= "in",res = 600)
plot.igraph(C1HTmedium.graph, layout= as.matrix(C1HTmedium.DS.obj@dr$tsne@cell.embeddings), edge.width = E(graph = C1HTmedium.graph)$weight, vertex.label = NA,vertex.size = 0 )
dev.off()
#Silhouette
C1HTmedium.PCA.data <- C1HTmedium.DS.obj@dr$pca@cell.embeddings[,1:8]
C1HTmedium.dd <- as.matrix(dist(C1HTmedium.PCA.data))
C1HTmedium.sil <- summary(silhouette(as.numeric(C1HTmedium.clusters), C1HTmedium.dd))$avg.width


# Chromium Clustering ####
# loading Annotated seurat objects #
load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/chromium_seu.obj.RData")
Chromium.hsap.obj <- chromium
rm(chromium)
Chromium.hsap.obj@ident <- Chromium.hsap.obj@meta.data$nnet2
names(Chromium.hsap.obj@ident) <- rownames(Chromium.hsap.obj@meta.data)
Chromium.hsap.metadata <- Chromium.hsap.obj@meta.data

load("/Volumes/Ati-Archive/HCA/SCE_Robjects/Downsampled/DS_secondRound/Separated_DS_datasets/Chromium_DS_20K.RData")
Chromium.20K <- mapIDs(Chromium.20K, "hsap")
colnames(Chromium.20K) <- gsub("\\.", "_", colnames(Chromium.20K))
common.cells <- intersect(colnames(Chromium.20K), colnames(Chromium.hsap.obj@scale.data))
common.genes <- intersect(rownames(Chromium.20K), rownames(Chromium.hsap.obj@scale.data))
Chromium.20K.clean <- Chromium.20K[common.genes, common.cells]

Chromium.DS.obj <- CreateSeuratObject(raw.data = as.matrix(Chromium.20K.clean), project = "Chromium")
Chromium.DS.obj <- AddMetaData(object = Chromium.DS.obj, metadata = Chromium.hsap.obj@ident, col.name = "clusters")

Chromium.DS.obj <- NormalizeData(object = Chromium.DS.obj, normalization.method = "LogNormalize", scale.factor = 10000)
dim(Chromium.DS.obj@data)
Chromium.DS.obj <- FindVariableGenes(object = Chromium.DS.obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.20, x.high.cutoff = 5, y.cutoff = 0.45)
length(Chromium.DS.obj@var.genes)
Chromium.DS.obj <- ScaleData(object = Chromium.DS.obj)#, vars.to.regress = c("nUMI", "percent.mito"))
Chromium.DS.obj <- RunPCA(object = Chromium.DS.obj, pc.genes = Chromium.DS.obj@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = Chromium.DS.obj, dim.1 = 1, dim.2 = 2)
#Chromium.DS.obj <- JackStraw(object = Chromium.DS.obj, num.replicate = 100, display.progress = FALSE)
#JackStrawPlot(object = Chromium.DS.obj, PCs = 1:12)
PCElbowPlot(object = Chromium.DS.obj)
Chromium.DS.obj <- FindClusters(object = Chromium.DS.obj, reduction.type = "pca", dims.use = 1:8, resolution = 0.6, print.output = 0, save.SNN = TRUE)
Chromium.DS.obj <- RunTSNE(object = Chromium.DS.obj, dims.use = 1:8, do.fast = TRUE)

Chromium.colors <- c("blueviolet","aquamarine","green4","maroon","orange","black","red","deepskyblue3", "grey")
png("Chromium_clustering_DS20K_clusters.png", width = 8, height = 6, units= "in",res = 600)
TSNEPlot(object = Chromium.DS.obj, group.by = "clusters", colors.use = Chromium.colors, pt.size = 2.5)
dev.off()
save(Chromium.DS.obj, file="Chromium_clustering_DS20K_dim1to8_res06_OBJ.RData")
png(file="/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/elbowPlots/Chromium_clustering_DS20K_clusters_elbowplot.png", width = 6, height = 6, units= "in",res = 600)
PCElbowPlot(object = Chromium.DS.obj)
dev.off()


library(scran)
library(igraph)
Chromium.graph <- graph.adjacency(Chromium.DS.obj@snn, weighted = T, mode="undirected")
Chromium.clusters <- Chromium.DS.obj@ident
Chromium.graph.mod.score <- modularity(Chromium.graph, Chromium.clusters)
Chromium.cluster.mod.score <- clusterModularity(Chromium.graph, Chromium.clusters)#, get.values=FALSE) 
rm(Chromium.hsap.obj)
rm(Chromium.20K)
png(file="Chromium_clustering_DS20K_graph.png", width = 8, height = 6, units= "in",res = 600)
plot.igraph(Chromium.graph, layout= as.matrix(Chromium.DS.obj@dr$tsne@cell.embeddings), edge.width = E(graph = Chromium.graph)$weight, vertex.label = NA,vertex.size = 0 )
dev.off()
#Silhouette
Chromium.PCA.data <- Chromium.DS.obj@dr$pca@cell.embeddings[,1:8]
Chromium.dd <- as.matrix(dist(Chromium.PCA.data))
Chromium.sil <- summary(silhouette(as.numeric(Chromium.clusters), Chromium.dd))$avg.width


# ChromiumNuclei Clustering ####
# loading Annotated seurat objects #
load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/nuclei_seu.obj.RData")
ChromiumNuclei.hsap.obj <- nuclei
rm(nuclei)
ChromiumNuclei.hsap.obj@ident <- ChromiumNuclei.hsap.obj@meta.data$nnet2
names(ChromiumNuclei.hsap.obj@ident) <- rownames(ChromiumNuclei.hsap.obj@meta.data)
ChromiumNuclei.hsap.metadata <- ChromiumNuclei.hsap.obj@meta.data

load("/Volumes/Ati-Archive/HCA/SCE_Robjects/Downsampled/DS_secondRound/Separated_DS_datasets/ChromiumNuclei_DS_20K.RData")
ChromiumNuclei.20K <- mapIDs(ChromiumNuclei.20K, "hsap")
colnames(ChromiumNuclei.20K) <- gsub("\\.", "_", colnames(ChromiumNuclei.20K))
common.cells <- intersect(colnames(ChromiumNuclei.20K), colnames(ChromiumNuclei.hsap.obj@scale.data))
common.genes <- intersect(rownames(ChromiumNuclei.20K), rownames(ChromiumNuclei.hsap.obj@scale.data))
ChromiumNuclei.20K.clean <- ChromiumNuclei.20K[common.genes, common.cells]

ChromiumNuclei.DS.obj <- CreateSeuratObject(raw.data = as.matrix(ChromiumNuclei.20K.clean), project = "ChromiumNuclei")
ChromiumNuclei.DS.obj <- AddMetaData(object = ChromiumNuclei.DS.obj, metadata = ChromiumNuclei.hsap.obj@ident, col.name = "clusters")

ChromiumNuclei.DS.obj <- NormalizeData(object = ChromiumNuclei.DS.obj, normalization.method = "LogNormalize", scale.factor = 10000)
dim(ChromiumNuclei.DS.obj@data)
ChromiumNuclei.DS.obj <- FindVariableGenes(object = ChromiumNuclei.DS.obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.28, x.high.cutoff = 5, y.cutoff = 0.5)
length(ChromiumNuclei.DS.obj@var.genes)
ChromiumNuclei.DS.obj <- ScaleData(object = ChromiumNuclei.DS.obj)#, vars.to.regress = c("nUMI", "percent.mito"))
ChromiumNuclei.DS.obj <- RunPCA(object = ChromiumNuclei.DS.obj, pc.genes = ChromiumNuclei.DS.obj@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = ChromiumNuclei.DS.obj, dim.1 = 1, dim.2 = 2)
#ChromiumNuclei.DS.obj <- JackStraw(object = ChromiumNuclei.DS.obj, num.replicate = 100, display.progress = FALSE)
#JackStrawPlot(object = ChromiumNuclei.DS.obj, PCs = 1:12)
PCElbowPlot(object = ChromiumNuclei.DS.obj)
ChromiumNuclei.DS.obj <- FindClusters(object = ChromiumNuclei.DS.obj, reduction.type = "pca", dims.use = 1:8, resolution = 0.6, print.output = 0, save.SNN = TRUE)
ChromiumNuclei.DS.obj <- RunTSNE(object = ChromiumNuclei.DS.obj, dims.use = 1:8, do.fast = TRUE)

ChromiumNuclei.colors <- c("blueviolet","aquamarine","green4","maroon","orange","black","red","coral2","deepskyblue3", "grey")
png("ChromiumNuclei_clustering_DS20K_clusters.png", width = 8, height = 6, units= "in",res = 600)
TSNEPlot(object = ChromiumNuclei.DS.obj, group.by = "clusters", colors.use = ChromiumNuclei.colors, pt.size = 2.5)
dev.off()
save(ChromiumNuclei.DS.obj, file="ChromiumNuclei_clustering_DS20K_dim1to8_res06_OBJ.RData")
length(which(ChromiumNuclei.DS.obj@meta.data$clusters == "Megakaryocytes"))
png(file="/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/elbowPlots/ChromiumNuclei_clustering_DS20K_clusters_elbowplot.png", width = 6, height = 6, units= "in",res = 600)
PCElbowPlot(object = ChromiumNuclei.DS.obj)
dev.off()


library(scran)
library(igraph)
ChromiumNuclei.graph <- graph.adjacency(ChromiumNuclei.DS.obj@snn, weighted = T, mode="undirected")
ChromiumNuclei.clusters <- ChromiumNuclei.DS.obj@ident
ChromiumNuclei.graph.mod.score <- modularity(ChromiumNuclei.graph, ChromiumNuclei.clusters)
ChromiumNuclei.cluster.mod.score <- clusterModularity(ChromiumNuclei.graph, ChromiumNuclei.clusters)#, get.values=FALSE) 
rm(ChromiumNuclei.hsap.obj)
rm(ChromiumNuclei.20K)
png(file="ChromiumNuclei_clustering_DS20K_graph.png", width = 8, height = 6, units= "in",res = 600)
plot.igraph(ChromiumNuclei.graph, layout= as.matrix(ChromiumNuclei.DS.obj@dr$tsne@cell.embeddings), edge.width = E(graph = ChromiumNuclei.graph)$weight, vertex.label = NA,vertex.size = 0 )
dev.off()
#Silhouette
ChromiumNuclei.PCA.data <- ChromiumNuclei.DS.obj@dr$pca@cell.embeddings[,1:8]
ChromiumNuclei.dd <- as.matrix(dist(ChromiumNuclei.PCA.data))
ChromiumNuclei.sil <- summary(silhouette(as.numeric(ChromiumNuclei.clusters), ChromiumNuclei.dd))$avg.width


# ddSEQ Clustering ####
# loading Annotated seurat objects #
load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/ddseq_seu.obj.RData")
ddSEQ.hsap.obj <- ddseq
rm(ddseq)
ddSEQ.hsap.obj@ident <- ddSEQ.hsap.obj@meta.data$nnet2
names(ddSEQ.hsap.obj@ident) <- rownames(ddSEQ.hsap.obj@meta.data)
ddSEQ.hsap.metadata <- ddSEQ.hsap.obj@meta.data

load("/Volumes/Ati-Archive/HCA/SCE_Robjects/Downsampled/DS_secondRound/Separated_DS_datasets/ddSEQ_DS_20K.RData")
ddSEQ.20K <- mapIDs(ddSEQ.20K, "hsap")
colnames(ddSEQ.20K) <- gsub("\\.", "_", colnames(ddSEQ.20K))
common.cells <- intersect(colnames(ddSEQ.20K), colnames(ddSEQ.hsap.obj@scale.data))
common.genes <- intersect(rownames(ddSEQ.20K), rownames(ddSEQ.hsap.obj@scale.data))
ddSEQ.20K.clean <- ddSEQ.20K[common.genes, common.cells]

ddSEQ.DS.obj <- CreateSeuratObject(raw.data = as.matrix(ddSEQ.20K.clean), project = "ddSEQ")
ddSEQ.DS.obj <- AddMetaData(object = ddSEQ.DS.obj, metadata = ddSEQ.hsap.obj@ident, col.name = "clusters")

ddSEQ.DS.obj <- NormalizeData(object = ddSEQ.DS.obj, normalization.method = "LogNormalize", scale.factor = 10000)
dim(ddSEQ.DS.obj@data)
ddSEQ.DS.obj <- FindVariableGenes(object = ddSEQ.DS.obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.23, x.high.cutoff = 5, y.cutoff = 0.5)
length(ddSEQ.DS.obj@var.genes)
ddSEQ.DS.obj <- ScaleData(object = ddSEQ.DS.obj)#, vars.to.regress = c("nUMI", "percent.mito"))
ddSEQ.DS.obj <- RunPCA(object = ddSEQ.DS.obj, pc.genes = ddSEQ.DS.obj@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = ddSEQ.DS.obj, dim.1 = 1, dim.2 = 2)
#ddSEQ.DS.obj <- JackStraw(object = ddSEQ.DS.obj, num.replicate = 100, display.progress = FALSE)
#JackStrawPlot(object = ddSEQ.DS.obj, PCs = 1:12)
PCElbowPlot(object = ddSEQ.DS.obj)
ddSEQ.DS.obj <- FindClusters(object = ddSEQ.DS.obj, reduction.type = "pca", dims.use = 1:8, resolution = 0.6, print.output = 0, save.SNN = TRUE)
ddSEQ.DS.obj <- RunTSNE(object = ddSEQ.DS.obj, dims.use = 1:8, do.fast = TRUE)

ddSEQ.colors <- c("blueviolet","aquamarine","green4","maroon","orange","black","red","coral2","deepskyblue3", "grey")
png("ddSEQ_clustering_DS20K_clusters.png", width = 8, height = 6, units= "in",res = 600)
TSNEPlot(object = ddSEQ.DS.obj, group.by = "clusters", colors.use = ddSEQ.colors, pt.size = 2)
dev.off()
save(ddSEQ.DS.obj, file="ddSEQ_clustering_DS20K_dim1to8_res06_OBJ.RData")
png(file="/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/elbowPlots/ddSEQ_clustering_DS20K_clusters_elbowplot.png", width = 6, height = 6, units= "in",res = 600)
PCElbowPlot(object = ddSEQ.DS.obj)
dev.off()

library(scran)
library(igraph)
ddSEQ.graph <- graph.adjacency(ddSEQ.DS.obj@snn, weighted = T, mode="undirected")
ddSEQ.clusters <- ddSEQ.DS.obj@ident
ddSEQ.graph.mod.score <- modularity(ddSEQ.graph, ddSEQ.clusters)
ddSEQ.cluster.mod.score <- clusterModularity(ddSEQ.graph, ddSEQ.clusters)#, get.values=FALSE) 
rm(ddSEQ.hsap.obj)
rm(ddSEQ.20K)
png(file="ddSEQ_clustering_DS20K_graph.png", width = 8, height = 6, units= "in",res = 600)
plot.igraph(ddSEQ.graph, layout= as.matrix(ddSEQ.DS.obj@dr$tsne@cell.embeddings), edge.width = E(graph = ddSEQ.graph)$weight, vertex.label = NA,vertex.size = 0 )
dev.off()
#Silhouette
ddSEQ.PCA.data <- ddSEQ.DS.obj@dr$pca@cell.embeddings[,1:8]
ddSEQ.dd <- as.matrix(dist(ddSEQ.PCA.data))
ddSEQ.sil <- summary(silhouette(as.numeric(ddSEQ.clusters), ddSEQ.dd))$avg.width


# Dropseq Clustering ####
# loading Annotated seurat objects #
load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/dropseq_seu.obj.RData")
Dropseq.hsap.obj <- dropseq
rm(dropseq)
Dropseq.hsap.obj@ident <- Dropseq.hsap.obj@meta.data$nnet2
names(Dropseq.hsap.obj@ident) <- rownames(Dropseq.hsap.obj@meta.data)
Dropseq.hsap.metadata <- Dropseq.hsap.obj@meta.data

load("/Volumes/Ati-Archive/HCA/SCE_Robjects/Downsampled/DS_secondRound/Separated_DS_datasets/Dropseq_DS_20K.RData")
Dropseq.20K <- mapIDs(Dropseq.20K, "hsap")
colnames(Dropseq.20K) <- gsub("\\.", "_", colnames(Dropseq.20K))
common.cells <- intersect(colnames(Dropseq.20K), colnames(Dropseq.hsap.obj@scale.data))
common.genes <- intersect(rownames(Dropseq.20K), rownames(Dropseq.hsap.obj@scale.data))
Dropseq.20K.clean <- Dropseq.20K[common.genes, common.cells]

Dropseq.DS.obj <- CreateSeuratObject(raw.data = as.matrix(Dropseq.20K.clean), project = "Dropseq")
Dropseq.DS.obj <- AddMetaData(object = Dropseq.DS.obj, metadata = Dropseq.hsap.obj@ident, col.name = "clusters")

Dropseq.DS.obj <- NormalizeData(object = Dropseq.DS.obj, normalization.method = "LogNormalize", scale.factor = 10000)
dim(Dropseq.DS.obj@data)
Dropseq.DS.obj <- FindVariableGenes(object = Dropseq.DS.obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.275, x.high.cutoff = 5, y.cutoff = 0.5)
length(Dropseq.DS.obj@var.genes)
Dropseq.DS.obj <- ScaleData(object = Dropseq.DS.obj)#, vars.to.regress = c("nUMI", "percent.mito"))
Dropseq.DS.obj <- RunPCA(object = Dropseq.DS.obj, pc.genes = Dropseq.DS.obj@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = Dropseq.DS.obj, dim.1 = 1, dim.2 = 2)
#Dropseq.DS.obj <- JackStraw(object = Dropseq.DS.obj, num.replicate = 100, display.progress = FALSE)
#JackStrawPlot(object = Dropseq.DS.obj, PCs = 1:12)
PCElbowPlot(object = Dropseq.DS.obj)
Dropseq.DS.obj <- FindClusters(object = Dropseq.DS.obj, reduction.type = "pca", dims.use = 1:8, resolution = 0.6, print.output = 0, save.SNN = TRUE)
Dropseq.DS.obj <- RunTSNE(object = Dropseq.DS.obj, dims.use = 1:8, do.fast = TRUE)

Dropseq.colors <- c("blueviolet","aquamarine","green4","maroon","orange","black","red","coral2","deepskyblue3", "grey")
png("Dropseq_clustering_DS20K_clusters.png", width = 8, height = 6, units= "in",res = 600)
TSNEPlot(object = Dropseq.DS.obj, group.by = "clusters", colors.use = Dropseq.colors, pt.size = 2.5)
dev.off()
save(Dropseq.DS.obj, file="Dropseq_clustering_DS20K_dim1to8_res06_OBJ.RData")
png(file="/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/elbowPlots/Dropseq_clustering_DS20K_clusters_elbowplot.png", width = 6, height = 6, units= "in",res = 600)
PCElbowPlot(object = Dropseq.DS.obj)
dev.off()


library(scran)
library(igraph)
Dropseq.graph <- graph.adjacency(Dropseq.DS.obj@snn, weighted = T, mode="undirected")
Dropseq.clusters <- Dropseq.DS.obj@ident
Dropseq.graph.mod.score <- modularity(Dropseq.graph, Dropseq.clusters)
Dropseq.cluster.mod.score <- clusterModularity(Dropseq.graph, Dropseq.clusters)#, get.values=FALSE) 
rm(Dropseq.hsap.obj)
rm(Dropseq.20K)
png(file="Dropseq_clustering_DS20K_graph.png", width = 8, height = 6, units= "in",res = 600)
plot.igraph(Dropseq.graph, layout= as.matrix(Dropseq.DS.obj@dr$tsne@cell.embeddings), edge.width = E(graph = Dropseq.graph)$weight, vertex.label = NA,vertex.size = 0 )
dev.off()

#Silhouette
Dropseq.PCA.data <- Dropseq.DS.obj@dr$pca@cell.embeddings[,1:8]
Dropseq.dd <- as.matrix(dist(Dropseq.PCA.data))
Dropseq.sil <- summary(silhouette(as.numeric(Dropseq.clusters), Dropseq.dd))$avg.width


# ICELL8 Clustering ####
# loading Annotated seurat objects #
load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/icell8_seu.obj.RData")
ICELL8.hsap.obj <- icell8
rm(icell8)
ICELL8.hsap.obj@ident <- ICELL8.hsap.obj@meta.data$nnet2
names(ICELL8.hsap.obj@ident) <- rownames(ICELL8.hsap.obj@meta.data)
ICELL8.hsap.metadata <- ICELL8.hsap.obj@meta.data

load("/Volumes/Ati-Archive/HCA/SCE_Robjects/Downsampled/DS_secondRound/Separated_DS_datasets/ICELL8_DS_20K.RData")
ICELL8.20K <- mapIDs(ICELL8.20K, "hsap")
colnames(ICELL8.20K) <- gsub("\\.", "_", colnames(ICELL8.20K))
common.cells <- intersect(colnames(ICELL8.20K), colnames(ICELL8.hsap.obj@scale.data))
common.genes <- intersect(rownames(ICELL8.20K), rownames(ICELL8.hsap.obj@scale.data))
ICELL8.20K.clean <- ICELL8.20K[common.genes, common.cells]

ICELL8.DS.obj <- CreateSeuratObject(raw.data = as.matrix(ICELL8.20K.clean), project = "ICELL8")
ICELL8.DS.obj <- AddMetaData(object = ICELL8.DS.obj, metadata = ICELL8.hsap.obj@ident, col.name = "clusters")

ICELL8.DS.obj <- NormalizeData(object = ICELL8.DS.obj, normalization.method = "LogNormalize", scale.factor = 10000)
dim(ICELL8.DS.obj@data)
ICELL8.DS.obj <- FindVariableGenes(object = ICELL8.DS.obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.3, x.high.cutoff = 5, y.cutoff = 0.5)
length(ICELL8.DS.obj@var.genes)
ICELL8.DS.obj <- ScaleData(object = ICELL8.DS.obj)#, vars.to.regress = c("nUMI", "percent.mito"))
ICELL8.DS.obj <- RunPCA(object = ICELL8.DS.obj, pc.genes = ICELL8.DS.obj@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = ICELL8.DS.obj, dim.1 = 1, dim.2 = 2)
#ICELL8.DS.obj <- JackStraw(object = ICELL8.DS.obj, num.replicate = 100, display.progress = FALSE)
#JackStrawPlot(object = ICELL8.DS.obj, PCs = 1:12)
PCElbowPlot(object = ICELL8.DS.obj)
ICELL8.DS.obj <- FindClusters(object = ICELL8.DS.obj, reduction.type = "pca", dims.use = 1:8, resolution = 0.6, print.output = 0, save.SNN = TRUE)
ICELL8.DS.obj <- RunTSNE(object = ICELL8.DS.obj, dims.use = 1:8, do.fast = TRUE)

ICELL8.colors <- c("blueviolet","aquamarine","green4","maroon","orange","black","red","coral2","deepskyblue3", "grey")
png("ICELL8_clustering_DS20K_clusters.png", width = 8, height = 6, units= "in",res = 600)
TSNEPlot(object = ICELL8.DS.obj, group.by = "clusters", colors.use = ICELL8.colors, pt.size = 2.5)
dev.off()
save(ICELL8.DS.obj, file="ICELL8_clustering_DS20K_dim1to8_res06_OBJ.RData")
png(file="/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/elbowPlots/ICELL8_clustering_DS20K_clusters_elbowplot.png", width = 6, height = 6, units= "in",res = 600)
PCElbowPlot(object = ICELL8.DS.obj)
dev.off()


library(scran)
library(igraph)
ICELL8.graph <- graph.adjacency(ICELL8.DS.obj@snn, weighted = T, mode="undirected")
ICELL8.clusters <- ICELL8.DS.obj@ident
ICELL8.graph.mod.score <- modularity(ICELL8.graph, ICELL8.clusters)
ICELL8.cluster.mod.score <- clusterModularity(ICELL8.graph, ICELL8.clusters)#, get.values=FALSE) 
rm(ICELL8.hsap.obj)
rm(ICELL8.20K)
png(file="ICELL8_clustering_DS20K_graph.png", width = 8, height = 6, units= "in",res = 600)
plot.igraph(ICELL8.graph, layout= as.matrix(ICELL8.DS.obj@dr$tsne@cell.embeddings), edge.width = E(graph = ICELL8.graph)$weight, vertex.label = NA,vertex.size = 0 )
dev.off()
#Silhouette
ICELL8.PCA.data <- ICELL8.DS.obj@dr$pca@cell.embeddings[,1:8]
ICELL8.dd <- as.matrix(dist(ICELL8.PCA.data))
ICELL8.sil <- summary(silhouette(as.numeric(ICELL8.clusters), ICELL8.dd))$avg.width


# inDrop Clustering ####
# loading Annotated seurat objects #
load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/onecb_seu.obj.RData")
inDrop.hsap.obj <- onecb
rm(onecb)
inDrop.hsap.obj@ident <- inDrop.hsap.obj@meta.data$nnet2
names(inDrop.hsap.obj@ident) <- rownames(inDrop.hsap.obj@meta.data)
inDrop.hsap.metadata <- inDrop.hsap.obj@meta.data

load("/Volumes/Ati-Archive/HCA/SCE_Robjects/Downsampled/DS_secondRound/Separated_DS_datasets/inDrop_DS_20K.RData")
inDrop.20K <- mapIDs(inDrop.20K, "hsap")
colnames(inDrop.20K) <- gsub("\\.", "_", colnames(inDrop.20K))
common.cells <- intersect(colnames(inDrop.20K), colnames(inDrop.hsap.obj@scale.data))
common.genes <- intersect(rownames(inDrop.20K), rownames(inDrop.hsap.obj@scale.data))
inDrop.20K.clean <- inDrop.20K[common.genes, common.cells]

inDrop.DS.obj <- CreateSeuratObject(raw.data = as.matrix(inDrop.20K.clean), project = "inDrop")
inDrop.DS.obj <- AddMetaData(object = inDrop.DS.obj, metadata = inDrop.hsap.obj@ident, col.name = "clusters")

inDrop.DS.obj <- NormalizeData(object = inDrop.DS.obj, normalization.method = "LogNormalize", scale.factor = 10000)
dim(inDrop.DS.obj@data)
inDrop.DS.obj <- FindVariableGenes(object = inDrop.DS.obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.24, x.high.cutoff = 5, y.cutoff = 0.5)
length(inDrop.DS.obj@var.genes)
inDrop.DS.obj <- ScaleData(object = inDrop.DS.obj)#, vars.to.regress = c("nUMI", "percent.mito"))
inDrop.DS.obj <- RunPCA(object = inDrop.DS.obj, pc.genes = inDrop.DS.obj@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = inDrop.DS.obj, dim.1 = 1, dim.2 = 2)
#inDrop.DS.obj <- JackStraw(object = inDrop.DS.obj, num.replicate = 100, display.progress = FALSE)
#JackStrawPlot(object = inDrop.DS.obj, PCs = 1:12)
PCElbowPlot(object = inDrop.DS.obj)
inDrop.DS.obj <- FindClusters(object = inDrop.DS.obj, reduction.type = "pca", dims.use = 1:8, resolution = 0.6, print.output = 0, save.SNN = TRUE)
inDrop.DS.obj <- RunTSNE(object = inDrop.DS.obj, dims.use = 1:8, do.fast = TRUE)

inDrop.colors <- c("blueviolet","aquamarine","green4","maroon","orange","black","red","deepskyblue3", "grey")
png("inDrop_clustering_DS20K_clusters.png", width = 8, height = 6, units= "in",res = 600)
TSNEPlot(object = inDrop.DS.obj, group.by = "clusters", colors.use = inDrop.colors, pt.size = 2.5)
dev.off()
save(inDrop.DS.obj, file="inDrop_clustering_DS20K_dim1to8_res06_OBJ.RData")
png(file="/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/elbowPlots/inDrop_clustering_DS20K_clusters_elbowplot.png", width = 6, height = 6, units= "in",res = 600)
PCElbowPlot(object = inDrop.DS.obj)
dev.off()


library(scran)
library(igraph)
inDrop.graph <- graph.adjacency(inDrop.DS.obj@snn, weighted = T, mode="undirected")
inDrop.clusters <- inDrop.DS.obj@ident
inDrop.graph.mod.score <- modularity(inDrop.graph, inDrop.clusters)
inDrop.cluster.mod.score <- clusterModularity(inDrop.graph, inDrop.clusters)#, get.values=FALSE) 
rm(inDrop.hsap.obj)
rm(inDrop.20K)
png(file="inDrop_clustering_DS20K_graph.png", width = 8, height = 6, units= "in",res = 600)
plot.igraph(inDrop.graph, layout= as.matrix(inDrop.DS.obj@dr$tsne@cell.embeddings), edge.width = E(graph = inDrop.graph)$weight, vertex.label = NA,vertex.size = 0 )
dev.off()
#Silhouette
inDrop.PCA.data <- inDrop.DS.obj@dr$pca@cell.embeddings[,1:8]
inDrop.dd <- as.matrix(dist(inDrop.PCA.data))
inDrop.sil <- summary(silhouette(as.numeric(inDrop.clusters), inDrop.dd))$avg.width



modularity.scores <- c(CELseq2.graph.mod.score, MARSseq.graph.mod.score, QUARTZseq.graph.mod.score, SCRBseq.graph.mod.score, SMARTseq2.graph.mod.score,
                       C1HTsmall.graph.mod.score, C1HTmedium.graph.mod.score, Chromium.graph.mod.score, ChromiumNuclei.graph.mod.score, 
                       ddSEQ.graph.mod.score, Dropseq.graph.mod.score, ICELL8.graph.mod.score, inDrop.graph.mod.score)

names(modularity.scores) <- c("CEL-Seq2","MARS-Seq","Quartz-Seq2","mcSCRB-Seq", "Smart-Seq2", "C1HT-Small","C1HT-Medium","Chromium", "Chromium(sn)","ddSEQ","Drop-Seq","ICELL8","inDrop") 

save(modularity.scores, file = "modularity_scores.RData")
#modularity.scores <- modularity.scores -0.4
mod.df <- as.data.frame(modularity.scores)
mod.df <- cbind(mod.df, tech=rownames(mod.df))
#mod.df <- mod.df[order(-mod.df$modularity.scores),]
mod.df$tech <- factor(mod.df$tech, levels = mod.df$tech)

#selected.color = c("tan2", "tomato1", "red3", "palevioletred1", "maroon2","darkmagenta", "#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933")
selected.color <- c("red3","green3","gold","blue1","palevioletred1","darkmagenta","sienna2","darkorange4", "turquoise3","darkgrey","#117733","maroon2","#999933")

pdf("Modularity_scores_ordered.pdf")
ggplot(mod.df, aes(x = tech, y = modularity.scores, fill = tech)) + geom_bar(stat="identity", fill =selected.color) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()



silhouette.scores <- c(CELseq2.sil, MARSseq.sil, QUARTZseq.sil, SCRBseq.sil, SMARTseq2.sil,
                       C1HTsmall.sil, C1HTmedium.sil, Chromium.sil, ChromiumNuclei.sil, 
                       ddSEQ.sil, Dropseq.sil, ICELL8.sil, inDrop.sil)

names(silhouette.scores) <- c("CEL-Seq2","MARS-Seq","Quartz-Seq2","mcSCRB-Seq", "Smart-Seq2", "C1HT-Small","C1HT-Medium","Chromium", "Chromium(sn)","ddSEQ","Drop-Seq","ICELL8","inDrop") 

save(silhouette.scores, file = "Silhouette_scores.RData")

sil.df <- as.data.frame(silhouette.scores)
sil.df <- cbind(sil.df, tech=rownames(sil.df))
#sil.df <- sil.df[order(-sil.df$silhouette.scores),]
sil.df$tech <- factor(sil.df$tech, levels = sil.df$tech)

#selected.color = c("tan2", "tomato1", "red3", "palevioletred1", "maroon2","darkmagenta", "#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933")
selected.color <- c("red3","green3","blue1","darkmagenta","sienna2","darkorange4", "turquoise3","darkgrey","maroon2","#999933")
#colors.for.silhouette <- ("palevioletred1","gold","#117733",)
pdf("Silhouette_scores.pdf")
ggplot(sil.df, aes(x = tech, y = silhouette.scores, fill = tech)) + geom_bar(stat="identity", fill =selected.color) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

#cluster accuracy measures
load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/clustering_accuracy/clus.accuracy.RData")
clus_accuracy
names(clus_accuracy) <- c("C1HT-Medium","C1HT-Small","CEL-Seq2","Chromium", "Chromium(sn)","ddSEQ","Drop-Seq","ICELL8","inDrop","MARS-Seq","Quartz-Seq2","mcSCRB-Seq", "Smart-Seq2")
clust.acc.df <- as.data.frame(clus_accuracy)
clust.acc.df <- cbind(clust.acc.df, tech=rownames(clust.acc.df))
#clust.acc.df <- clust.acc.df[order(-clust.acc.df$silhouette.scores),]
clust.acc.df$tech <- factor(clust.acc.df$tech, levels = clust.acc.df$tech)


#bind.mod.df <- mod.df
#colnames(bind.mod.df) <- c("score", "tech")
#bind.mod.df <- bind.mod.df[order(-bind.mod.df$score),]
bind.clust.acc.df <- clust.acc.df
colnames(bind.clust.acc.df) <- c("score", "tech")
bind.clust.acc.df <- bind.clust.acc.df[order(-bind.clust.acc.df$score),]

bind.sil.df <- sil.df
colnames(bind.sil.df) <- c("score", "tech")
bind.sil.df <- bind.sil.df[order(-bind.sil.df$score),]
#mod.sil.df <- rbind(bind.mod.df, bind.sil.df)
#mod.sil.df <- cbind(mod.sil.df, measure= c(rep("Modularity", nrow(bind.mod.df)), rep("Silhouette", nrow(bind.sil.df))))
#mod.sil.df$tech <- factor(mod.sil.df$tech, levels = rownames(bind.mod.df))#[13:1]
#ggplot(mod.sil.df, aes(x = tech, y = score, fill = measure)) + geom_bar(stat="identity", position = "dodge") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
clust.acc.sil.df <- rbind(bind.clust.acc.df, bind.sil.df[rownames(bind.clust.acc.df),])
clust.acc.sil.df <- cbind(clust.acc.sil.df, measure= c(rep("Cluster Accuracy", nrow(bind.clust.acc.df)), rep("Silhouette", nrow(bind.sil.df))))
clust.acc.sil.df$tech <- factor(clust.acc.sil.df$tech, levels = rownames(bind.clust.acc.df))#[13:1]
ggplot(clust.acc.sil.df, aes(x = tech, y = score, fill = measure)) + geom_bar(stat="identity", position = "dodge") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))


#selected.color.2 <- selected.color[13:1]
#mod.sil.df2 <- mod.sil.df
#mod.sil.df2$score[which(mod.sil.df2$measure == "Modularity")] <- mod.sil.df2$score[which(mod.sil.df2$measure == "Modularity")] * -1
#breaks2 <- seq(-1,1,0.5)
#names(breaks2) <- abs(breaks2)
clust.acc.sil.df2 <- clust.acc.sil.df
clust.acc.sil.df2$score[which(clust.acc.sil.df2$measure == "Cluster Accuracy")] <- clust.acc.sil.df2$score[which(clust.acc.sil.df2$measure == "Cluster Accuracy")] * -1
breaks2 <- seq(-1,1,0.5)
names(breaks2) <- abs(breaks2)
expand_data <- data.frame(measure = c("Cluster Accuracy", "Silhouette"),
                          score = c(-1,1))
clust.acc.sil.df2$tech <- factor(clust.acc.sil.df2$tech, levels = rev(rownames(bind.clust.acc.df)))

png(file="ClusterAccuracy_Silhouette_allTechs_orderedbyCA.png", width = 8, height = 6, units= "in",res = 600)
ggplot(clust.acc.sil.df2) +
  geom_col(aes(x = tech, y = score, fill = measure)) +
  geom_hline(yintercept = 0, colour = "black", size = 1) +
  geom_blank(aes(y = score), expand_data) +
  coord_flip() +
  facet_grid(~ measure, scales = "free_x") +
  scale_y_continuous(breaks = breaks2, labels = names(breaks2), 
                     expand = c(0, 0)) +
  theme(panel.spacing.x = unit(0, "lines"), axis.text.x = element_text(size = 15, face = "bold", hjust = 1),axis.text.y = element_text(size = 15, face = "bold", hjust = 1),
        strip.background = element_rect(colour = "black"), strip.text = element_text(face="bold", size=15)) +scale_fill_manual(values= c("cadetblue3", "lightpink2"))
dev.off() #geom_blank(aes(y = score), expand_data2) +
