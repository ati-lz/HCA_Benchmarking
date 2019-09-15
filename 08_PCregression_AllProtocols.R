library(ggExtra)
library(plyr)
library(gplots)
library(ggplot2)
library(stringr)
library(tidyr)
library(data.table)
library(corrplot)
library(stringi)
library(kBET)
library(plyr)
library(gtools)

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


load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/techniques_colors_2.RData")
selected.color

# Correlation plots
load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/celseq_seu.obj.RData")
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/Final/CELseq2_clustering_DS20K_dim1to8_res06_OBJ.RData")
comm.CELseq2.cells <- intersect(colnames(CELseq2.DS.obj@raw.data), rownames(celseq@meta.data))
CELseq2.annotation.metadata <- celseq@meta.data[comm.CELseq2.cells,]; rm(celseq)
CELseq2.HEK.mat <- as.data.frame(t(CELseq2.DS.obj@raw.data[,rownames(CELseq2.annotation.metadata[which(CELseq2.annotation.metadata$clean.id == "HEK"),])]))
CELseq2.replicates <- CELseq2.annotation.metadata[rownames(CELseq2.HEK.mat), "Library"]

CELseq2.HEK.mat.agg <- cbind(CELseq2.HEK.mat, CELseq2.replicates)
CELseq2.HEK.mat.agg.mean <- aggregate(CELseq2.HEK.mat.agg[,1:ncol(CELseq2.HEK.mat.agg)-1], by=list(CELseq2.replicates), FUN = mean)
CELseq2.HEK.mat.agg.mean.f <- as.matrix(CELseq2.HEK.mat.agg.mean)
rownames(CELseq2.HEK.mat.agg.mean.f) <- CELseq2.HEK.mat.agg.mean.f[,1]
CELseq2.HEK.mat.agg.mean.f <- CELseq2.HEK.mat.agg.mean.f[,-1]
class(CELseq2.HEK.mat.agg.mean.f) <- "numeric"
CELseq2.rep.cor <- cor(t(CELseq2.HEK.mat.agg.mean.f))
CELseq2.pca.data <- prcomp(CELseq2.HEK.mat, center=TRUE) #compute PCA representation of the data
CELseq2.batch.pca <- pcRegression(CELseq2.pca.data, CELseq2.replicates)
PC.regressions.all["CEL-Seq2"] <- CELseq2.batch.pca$pcRegscale


load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/marsseq_seu.obj.RData")
marsseq
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/Final/MARSseq_clustering_DS20K_dim1to8_res06_OBJ.RData")
comm.MARSseq.cells <- intersect(colnames(MARSseq.DS.obj@raw.data), rownames(marsseq@meta.data))
MARSseq.annotation.metadata <- marsseq@meta.data[comm.MARSseq.cells,]; rm(marsseq)
MARSseq.HEK.mat <- as.data.frame(t(MARSseq.DS.obj@raw.data[,rownames(MARSseq.annotation.metadata[which(MARSseq.annotation.metadata$clean.id == "HEK"),])]))
MARSseq.replicates <- MARSseq.annotation.metadata[rownames(MARSseq.HEK.mat), "Library"]

MARSseq.HEK.mat.agg <- cbind(MARSseq.HEK.mat, MARSseq.replicates)
MARSseq.HEK.mat.agg.mean <- aggregate(MARSseq.HEK.mat.agg[,1:ncol(MARSseq.HEK.mat.agg)-1], by=list(MARSseq.replicates), FUN = mean)
MARSseq.HEK.mat.agg.mean.f <- as.matrix(MARSseq.HEK.mat.agg.mean)
rownames(MARSseq.HEK.mat.agg.mean.f) <- MARSseq.HEK.mat.agg.mean.f[,1]
MARSseq.HEK.mat.agg.mean.f <- MARSseq.HEK.mat.agg.mean.f[,-1]
class(MARSseq.HEK.mat.agg.mean.f) <- "numeric"
MARSseq.rep.cor <- cor(t(MARSseq.HEK.mat.agg.mean.f))
MARSseq.pca.data <- prcomp(MARSseq.HEK.mat, center=TRUE) #compute PCA representation of the data
MARSseq.batch.pca <- pcRegression(MARSseq.pca.data, MARSseq.replicates)
PC.regressions.all["MARS-Seq"] <- MARSseq.batch.pca$pcRegscale


load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/quartzseq_seu.obj.RData")
quartzseq
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/Final/QUARTZseq_clustering_DS20K_dim1to8_res06_OBJ.RData")
comm.QUARTZseq.cells <- intersect(colnames(QUARTZseq.DS.obj@raw.data), rownames(quartzseq@meta.data))
QUARTZseq.annotation.metadata <- quartzseq@meta.data[comm.QUARTZseq.cells,]; rm(quartzseq)
QUARTZseq.HEK.mat <- as.data.frame(t(QUARTZseq.DS.obj@raw.data[,rownames(QUARTZseq.annotation.metadata[which(QUARTZseq.annotation.metadata$clean.id == "HEK"),])]))
QUARTZseq.replicates <- QUARTZseq.annotation.metadata[rownames(QUARTZseq.HEK.mat), "Library"]

QUARTZseq.HEK.mat.agg <- cbind(QUARTZseq.HEK.mat, QUARTZseq.replicates)
QUARTZseq.HEK.mat.agg.mean <- aggregate(QUARTZseq.HEK.mat.agg[,1:ncol(QUARTZseq.HEK.mat.agg)-1], by=list(QUARTZseq.replicates), FUN = mean)
QUARTZseq.HEK.mat.agg.mean.f <- as.matrix(QUARTZseq.HEK.mat.agg.mean)
rownames(QUARTZseq.HEK.mat.agg.mean.f) <- QUARTZseq.HEK.mat.agg.mean.f[,1]
QUARTZseq.HEK.mat.agg.mean.f <- QUARTZseq.HEK.mat.agg.mean.f[,-1]
class(QUARTZseq.HEK.mat.agg.mean.f) <- "numeric"
QUARTZseq.rep.cor <- cor(t(QUARTZseq.HEK.mat.agg.mean.f))
QUARTZseq.pca.data <- prcomp(QUARTZseq.HEK.mat, center=TRUE) #compute PCA representation of the data
QUARTZseq.batch.pca <- pcRegression(QUARTZseq.pca.data, QUARTZseq.replicates)
PC.regressions.all["Quartz-Seq2"] <- QUARTZseq.batch.pca$pcRegscale


load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/scrbseq_seu.obj.RData")
scrbseq
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/Final/SCRBseq_clustering_DS20K_dim1to8_res06_OBJ.RData")
comm.SCRBseq.cells <- intersect(colnames(SCRBseq.DS.obj@raw.data), rownames(scrbseq@meta.data))
SCRBseq.annotation.metadata <- scrbseq@meta.data[comm.SCRBseq.cells,]; rm(scrbseq)
SCRBseq.HEK.mat <- as.data.frame(t(SCRBseq.DS.obj@raw.data[,rownames(SCRBseq.annotation.metadata[which(SCRBseq.annotation.metadata$clean.id == "HEK"),])]))
SCRBseq.replicates <- SCRBseq.annotation.metadata[rownames(SCRBseq.HEK.mat), "Library"]
names(SCRBseq.replicates) <- rownames(SCRBseq.HEK.mat)
SCRBseq.rnd.selected.replicates <- c()
for (replicate in levels(SCRBseq.replicates)){
  replicate.cells <- names(SCRBseq.replicates[which(SCRBseq.replicates == replicate)])
  if (length(replicate.cells) >= 5){
    rand.select.rep <- sample(replicate.cells, 5)
    SCRBseq.rnd.selected.replicates <- c(SCRBseq.rnd.selected.replicates, rand.select.rep)
  }}
SCRBseq.rnd.selected.replicates <- SCRBseq.replicates[SCRBseq.rnd.selected.replicates]
SCRBseq.rnd.selected.replicates <- factor(SCRBseq.rnd.selected.replicates)



load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/smartseq_seu.obj.RData")
SS2.annotation <- read.table("/Volumes/Ati-Archive/HCA/SCE_Robjects/SMARTseq2/Smartseq2_HCA_annotation.txt", header = T)
cell.BCs <- paste("SMARTseqFINAL_allLanes_", SS2.annotation$Index1_sequence, SS2.annotation$Index2_sequence_HiSeq2500, sep = "")
BCs.plate <- cbind(cell.BCs, SS2.annotation$HCA_plate)
BCs.plate.batch <- c()
for (i in 1:nrow(BCs.plate)){
  if (BCs.plate[i,2] == "2"){BCs.plate.batch[i] <- "batch1"}
  if (BCs.plate[i,2] == "3" | BCs.plate[i,2] == "4"){BCs.plate.batch[i] <- "batch2"}
  if (BCs.plate[i,2] == "5" | BCs.plate[i,2] == "6"){BCs.plate.batch[i] <- "batch3"}
  if (BCs.plate[i,2] == "7" | BCs.plate[i,2] == "8"){BCs.plate.batch[i] <- "batch4"}
  if (BCs.plate[i,2] == "9"){BCs.plate.batch[i] <- "batch5"}
}
SS2.batches <- cbind(BCs.plate, BCs.plate.batch)
rownames(SS2.batches) <- SS2.batches[, "cell.BCs"]
SS2.batches <- as.factor(SS2.batches[,"BCs.plate.batch"])
smartseq
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/Final/SMARTseq2_clustering_DS20K_dim1to8_res06_OBJ.RData")
comm.SMARTseq2.cells <- intersect(colnames(SMARTseq2.DS.obj@raw.data), rownames(smartseq@meta.data))
SMARTseq2.annotation.metadata <- smartseq@meta.data[comm.SMARTseq2.cells,]; rm(smartseq)
SMARTseq2.HEK.mat <- as.data.frame(t(SMARTseq2.DS.obj@raw.data[,rownames(SMARTseq2.annotation.metadata[which(SMARTseq2.annotation.metadata$clean.id == "HEK"),])]))
SMARTseq2.replicates <- SS2.batches[rownames(SMARTseq2.HEK.mat)]

SMARTseq2.HEK.mat.agg <- cbind(SMARTseq2.HEK.mat, SMARTseq2.replicates)
SMARTseq2.HEK.mat.agg.mean <- aggregate(SMARTseq2.HEK.mat.agg[,1:ncol(SMARTseq2.HEK.mat.agg)-1], by=list(SMARTseq2.replicates), FUN = mean)
SMARTseq2.HEK.mat.agg.mean.f <- as.matrix(SMARTseq2.HEK.mat.agg.mean)
rownames(SMARTseq2.HEK.mat.agg.mean.f) <- SMARTseq2.HEK.mat.agg.mean.f[,1]
SMARTseq2.HEK.mat.agg.mean.f <- SMARTseq2.HEK.mat.agg.mean.f[,-1]
class(SMARTseq2.HEK.mat.agg.mean.f) <- "numeric"
SMARTseq2.rep.cor <- cor(t(SMARTseq2.HEK.mat.agg.mean.f))
SMARTseq2.pca.data <- prcomp(SMARTseq2.HEK.mat, center=TRUE) #compute PCA representation of the data
SMARTseq2.batch.pca <- pcRegression(SMARTseq2.pca.data, SMARTseq2.replicates)
PC.regressions.all["Smart-Seq2"] <- SMARTseq2.batch.pca$pcRegscale


load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/c1ht.s_seu.obj.RData")
c1ht.s
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/Final/C1HTsmall_clustering_DS20K_dim1to8_res06_OBJ.RData")
comm.C1HTsmall.cells <- intersect(colnames(C1HTsmall.DS.obj@raw.data), rownames(c1ht.s@meta.data))
C1HTsmall.annotation.metadata <- c1ht.s@meta.data[comm.C1HTsmall.cells,]; rm(c1ht.s)
C1HTsmall.HEK.mat <- as.data.frame(t(C1HTsmall.DS.obj@raw.data[,rownames(C1HTsmall.annotation.metadata[which(C1HTsmall.annotation.metadata$clean.id == "HEK"),])]))
C1HTsmall.replicates <- as.character(C1HTsmall.annotation.metadata[rownames(C1HTsmall.HEK.mat), "Library"])
for (i in 1:length(C1HTsmall.replicates)){
  id = C1HTsmall.replicates[i]
  id_parts = unlist(strsplit(id, split ="_"))
  if (length(id_parts) == 4 & grepl("Col", id_parts[3])){
    C1HTsmall.replicates[i] <- id_parts[1]
  }
}

C1HTsmall.HEK.mat.agg <- cbind(C1HTsmall.HEK.mat, C1HTsmall.replicates)
C1HTsmall.HEK.mat.agg.mean <- aggregate(C1HTsmall.HEK.mat.agg[,1:ncol(C1HTsmall.HEK.mat.agg)-1], by=list(C1HTsmall.replicates), FUN = mean)
C1HTsmall.HEK.mat.agg.mean.f <- as.matrix(C1HTsmall.HEK.mat.agg.mean)
rownames(C1HTsmall.HEK.mat.agg.mean.f) <- C1HTsmall.HEK.mat.agg.mean.f[,1]
C1HTsmall.HEK.mat.agg.mean.f <- C1HTsmall.HEK.mat.agg.mean.f[,-1]
class(C1HTsmall.HEK.mat.agg.mean.f) <- "numeric"
C1HTsmall.rep.cor <- cor(t(C1HTsmall.HEK.mat.agg.mean.f))
C1HTsmall.pca.data <- prcomp(C1HTsmall.HEK.mat, center=TRUE) #compute PCA representation of the data
C1HTsmall.batch.pca <- pcRegression(C1HTsmall.pca.data, C1HTsmall.replicates)
PC.regressions.all["C1HT-small"] <- C1HTsmall.batch.pca$pcRegscale


load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/c1ht.m_seu.obj.RData")
c1ht.m
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/Final/C1HTmedium_clustering_DS20K_dim1to8_res06_OBJ.RData")
comm.C1HTmedium.cells <- intersect(colnames(C1HTmedium.DS.obj@raw.data), rownames(c1ht.m@meta.data))
C1HTmedium.annotation.metadata <- c1ht.m@meta.data[comm.C1HTmedium.cells,]; rm(c1ht.m)
C1HTmedium.HEK.mat <- as.data.frame(t(C1HTmedium.DS.obj@raw.data[,rownames(C1HTmedium.annotation.metadata[which(C1HTmedium.annotation.metadata$clean.id == "HEK"),])]))
C1HTmedium.replicates <- as.character(C1HTmedium.annotation.metadata[rownames(C1HTmedium.HEK.mat), "Library"])
for (i in 1:length(C1HTmedium.replicates)){
  id = C1HTmedium.replicates[i]
  id_parts = unlist(strsplit(id, split ="_"))
  if (length(id_parts) == 4 & grepl("Col", id_parts[3])){
    C1HTmedium.replicates[i] <- id_parts[1]
  }
}

C1HTmedium.HEK.mat.agg <- cbind(C1HTmedium.HEK.mat, C1HTmedium.replicates)
C1HTmedium.HEK.mat.agg.mean <- aggregate(C1HTmedium.HEK.mat.agg[,1:ncol(C1HTmedium.HEK.mat.agg)-1], by=list(C1HTmedium.replicates), FUN = mean)
C1HTmedium.HEK.mat.agg.mean.f <- as.matrix(C1HTmedium.HEK.mat.agg.mean)
rownames(C1HTmedium.HEK.mat.agg.mean.f) <- C1HTmedium.HEK.mat.agg.mean.f[,1]
C1HTmedium.HEK.mat.agg.mean.f <- C1HTmedium.HEK.mat.agg.mean.f[,-1]
class(C1HTmedium.HEK.mat.agg.mean.f) <- "numeric"
C1HTmedium.rep.cor <- cor(t(C1HTmedium.HEK.mat.agg.mean.f))
C1HTmedium.pca.data <- prcomp(C1HTmedium.HEK.mat, center=TRUE) #compute PCA representation of the data
C1HTmedium.batch.pca <- pcRegression(C1HTmedium.pca.data, C1HTmedium.replicates)
PC.regressions.all["C1HT-medium"] <- C1HTmedium.batch.pca$pcRegscale


library(rowr)

dim(ChromV2.all.DS20K) # From the Chrom_Comparison_rev.R code
dim(ChromV2.all.DS20K.HEKs) # From the Chrom_Comparison_rev.R code
dim(ChromV2.all.metadata) # From the Chrom_Comparison_rev.R code
Chromium.HEK.mat <- ChromV2.all.DS20K.HEKs
Chromium.HEK.mat <- mapIDs(Chromium.HEK.mat, "hsap")
Chromium.HEK.mat <- t(Chromium.HEK.mat)
Chromium.replicates <- ChromV2.all.metadata[rownames(Chromium.HEK.mat), "Library"]

Chromium.HEK.mat.agg <- cbind(Chromium.HEK.mat, Chromium.replicates)
Chromium.HEK.mat.agg.mean <- aggregate(Chromium.HEK.mat.agg[,1:ncol(Chromium.HEK.mat.agg)-1], by=list(Chromium.replicates), FUN = mean)
Chromium.HEK.mat.agg.mean.f <- as.matrix(Chromium.HEK.mat.agg.mean)
rownames(Chromium.HEK.mat.agg.mean.f) <- Chromium.HEK.mat.agg.mean.f[,1]
Chromium.HEK.mat.agg.mean.f <- Chromium.HEK.mat.agg.mean.f[,-1]
class(Chromium.HEK.mat.agg.mean.f) <- "numeric"
Chromium.rep.cor <- cor(t(Chromium.HEK.mat.agg.mean.f))
Chromium.pca.data <- prcomp(Chromium.HEK.mat, center=TRUE) #compute PCA representation of the data
Chromium.batch.pca <- pcRegression(Chromium.pca.data, Chromium.replicates)
PC.regressions.all["Chromium"] <- Chromium.batch.pca$pcRegscale


#Chromium Ref
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/reference_10x/Ref10X_metadata.RData")
dim(Ref10X.metadata)
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/reference_10x/Ref10X_DS20K_expressionMat_HEK.RData")
dim(Ref10X.DS20K.HEK)
Ref10X.DS20K.HEK <- mapIDs(Ref10X.DS20K.HEK, "hsap")
Ref10X.DS20K.HEK <- t(Ref10X.DS20K.HEK)
#comm.ChromRef.cells <- intersect(rownames(DSth.df.HEK.Ref10X), rownames(Ref10X.metadata))
ChromRef.annotation.metadata <- Ref10X.metadata[rownames(Ref10X.DS20K.HEK),]
ChromRef.HEK.mat <- Ref10X.DS20K.HEK
ChromRef.replicates <- ChromRef.annotation.metadata[rownames(ChromRef.HEK.mat), "Library"]

ChromRef.HEK.mat.agg <- cbind(ChromRef.HEK.mat, ChromRef.replicates)
ChromRef.HEK.mat.agg.mean <- aggregate(ChromRef.HEK.mat.agg[,1:ncol(ChromRef.HEK.mat.agg)-1], by=list(ChromRef.replicates), FUN = mean)
ChromRef.HEK.mat.agg.mean.f <- as.matrix(ChromRef.HEK.mat.agg.mean)
rownames(ChromRef.HEK.mat.agg.mean.f) <- ChromRef.HEK.mat.agg.mean.f[,1]
ChromRef.HEK.mat.agg.mean.f <- ChromRef.HEK.mat.agg.mean.f[,-1]
class(ChromRef.HEK.mat.agg.mean.f) <- "numeric"
ChromRef.rep.cor <- cor(t(ChromRef.HEK.mat.agg.mean.f))
ChromRef.pca.data <- prcomp(ChromRef.HEK.mat, center=TRUE) #compute PCA representation of the data
ChromRef.batch.pca <- pcRegression(ChromRef.pca.data, ChromRef.replicates)
PC.regressions.all["Chromium(sn)"] <- ChromRef.batch.pca$pcRegscale


load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/nuclei_seu.obj.RData")
nuclei
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/Final/ChromiumNuclei_clustering_DS20K_dim1to8_res06_OBJ.RData")
comm.ChromNuclei.cells <- intersect(colnames(ChromiumNuclei.DS.obj@raw.data), rownames(nuclei@meta.data))
ChromNuclei.annotation.metadata <- nuclei@meta.data[comm.ChromNuclei.cells,]; rm(nuclei)
ChromNuclei.HEK.mat <- as.data.frame(t(ChromiumNuclei.DS.obj@raw.data[,rownames(ChromNuclei.annotation.metadata[which(ChromNuclei.annotation.metadata$clean.id == "HEK"),])]))
ChromNuclei.replicates <- ChromNuclei.annotation.metadata[rownames(ChromNuclei.HEK.mat), "Library"]

ChromNuclei.HEK.mat.agg <- cbind(ChromNuclei.HEK.mat, ChromNuclei.replicates)
ChromNuclei.HEK.mat.agg.mean <- aggregate(ChromNuclei.HEK.mat.agg[,1:ncol(ChromNuclei.HEK.mat.agg)-1], by=list(ChromNuclei.replicates), FUN = mean)
ChromNuclei.HEK.mat.agg.mean.f <- as.matrix(ChromNuclei.HEK.mat.agg.mean)
rownames(ChromNuclei.HEK.mat.agg.mean.f) <- ChromNuclei.HEK.mat.agg.mean.f[,1]
ChromNuclei.HEK.mat.agg.mean.f <- ChromNuclei.HEK.mat.agg.mean.f[,-1]
class(ChromNuclei.HEK.mat.agg.mean.f) <- "numeric"
ChromNuclei.rep.cor <- cor(t(ChromNuclei.HEK.mat.agg.mean.f))
ChromNuclei.pca.data <- prcomp(ChromNuclei.HEK.mat, center=TRUE) #compute PCA representation of the data
ChromNuclei.batch.pca <- pcRegression(ChromNuclei.pca.data, ChromNuclei.replicates)
PC.regressions.all["Chromium(sn)"] <- ChromNuclei.batch.pca$pcRegscale


load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/ddseq_seu.obj.RData")
ddseq
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/Final/ddSEQ_clustering_DS20K_dim1to8_res06_OBJ.RData")
comm.ddSEQ.cells <- intersect(colnames(ddSEQ.DS.obj@raw.data), rownames(ddseq@meta.data))
ddSEQ.annotation.metadata <- ddseq@meta.data[comm.ddSEQ.cells,]; rm(ddseq)
ddSEQ.HEK.mat <- as.data.frame(t(ddSEQ.DS.obj@raw.data[,rownames(ddSEQ.annotation.metadata[which(ddSEQ.annotation.metadata$clean.id == "HEK"),])]))
ddSEQ.replicates <- ddSEQ.annotation.metadata[rownames(ddSEQ.HEK.mat), "Library"]

ddSEQ.HEK.mat.agg <- cbind(ddSEQ.HEK.mat, ddSEQ.replicates)
ddSEQ.HEK.mat.agg.mean <- aggregate(ddSEQ.HEK.mat.agg[,1:ncol(ddSEQ.HEK.mat.agg)-1], by=list(ddSEQ.replicates), FUN = mean)
ddSEQ.HEK.mat.agg.mean.f <- as.matrix(ddSEQ.HEK.mat.agg.mean)
rownames(ddSEQ.HEK.mat.agg.mean.f) <- ddSEQ.HEK.mat.agg.mean.f[,1]
ddSEQ.HEK.mat.agg.mean.f <- ddSEQ.HEK.mat.agg.mean.f[,-1]
class(ddSEQ.HEK.mat.agg.mean.f) <- "numeric"
ddSEQ.rep.cor <- cor(t(ddSEQ.HEK.mat.agg.mean.f))
ddSEQ.pca.data <- prcomp(ddSEQ.HEK.mat, center=TRUE) #compute PCA representation of the data
ddSEQ.batch.pca <- pcRegression(ddSEQ.pca.data, ddSEQ.replicates)
PC.regressions.all["ddSEQ"] <- ddSEQ.batch.pca$pcRegscale


load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/dropseq_seu.obj.RData")
dropseq
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/Final/Dropseq_clustering_DS20K_dim1to8_res06_OBJ.RData")
comm.Dropseq.cells <- intersect(colnames(Dropseq.DS.obj@raw.data), rownames(dropseq@meta.data))
Dropseq.annotation.metadata <- dropseq@meta.data[comm.Dropseq.cells,]; rm(dropseq)
Dropseq.HEK.mat <- as.data.frame(t(Dropseq.DS.obj@raw.data[,rownames(Dropseq.annotation.metadata[which(Dropseq.annotation.metadata$clean.id == "HEK"),])]))
Dropseq.replicates <- Dropseq.annotation.metadata[rownames(Dropseq.HEK.mat), "Library"]

Dropseq.HEK.mat.agg <- cbind(Dropseq.HEK.mat, Dropseq.replicates)
Dropseq.HEK.mat.agg.mean <- aggregate(Dropseq.HEK.mat.agg[,1:ncol(Dropseq.HEK.mat.agg)-1], by=list(Dropseq.replicates), FUN = mean)
Dropseq.HEK.mat.agg.mean.f <- as.matrix(Dropseq.HEK.mat.agg.mean)
rownames(Dropseq.HEK.mat.agg.mean.f) <- Dropseq.HEK.mat.agg.mean.f[,1]
Dropseq.HEK.mat.agg.mean.f <- Dropseq.HEK.mat.agg.mean.f[,-1]
class(Dropseq.HEK.mat.agg.mean.f) <- "numeric"
Dropseq.rep.cor <- cor(t(Dropseq.HEK.mat.agg.mean.f))
Dropseq.pca.data <- prcomp(Dropseq.HEK.mat, center=TRUE) #compute PCA representation of the data
Dropseq.batch.pca <- pcRegression(Dropseq.pca.data, Dropseq.replicates)
PC.regressions.all["Drop-Seq"] <- Dropseq.batch.pca$pcRegscale


load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/icell8_seu.obj.RData")
icell8
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/Final/ICELL8_clustering_DS20K_dim1to8_res06_OBJ.RData")
comm.ICELL8.cells <- intersect(colnames(ICELL8.DS.obj@raw.data), rownames(icell8@meta.data))
ICELL8.annotation.metadata <- icell8@meta.data[comm.ICELL8.cells,]; rm(icell8)
ICELL8.HEK.mat <- as.data.frame(t(ICELL8.DS.obj@raw.data[,rownames(ICELL8.annotation.metadata[which(ICELL8.annotation.metadata$clean.id == "HEK"),])]))
ICELL8.replicates <- ICELL8.annotation.metadata[rownames(ICELL8.HEK.mat), "Library"]

ICELL8.HEK.mat.agg <- cbind(ICELL8.HEK.mat, ICELL8.replicates) #till here is the prev version
ICELL8.HEK.mat.agg.mean <- aggregate(ICELL8.HEK.mat.agg[,1:ncol(ICELL8.HEK.mat.agg)-1], by=list(ICELL8.replicates), FUN = mean)
ICELL8.HEK.mat.agg.mean.f <- as.matrix(ICELL8.HEK.mat.agg.mean)
rownames(ICELL8.HEK.mat.agg.mean.f) <- ICELL8.HEK.mat.agg.mean.f[,1]
ICELL8.HEK.mat.agg.mean.f <- ICELL8.HEK.mat.agg.mean.f[,-1]
class(ICELL8.HEK.mat.agg.mean.f) <- "numeric"
ICELL8.rep.cor <- cor(t(ICELL8.HEK.mat.agg.mean.f))
ICELL8.pca.data <- prcomp(ICELL8.HEK.mat, center=TRUE) #compute PCA representation of the data
ICELL8.batch.pca <- pcRegression(ICELL8.pca.data, ICELL8.replicates)
PC.regressions.all["ICELL8"] <- ICELL8.batch.pca$pcRegscale


load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/onecb_seu.obj.RData")
inDrop.hsap.obj <- onecb
rm(onecb)
inDrop.hsap.metadata <- inDrop.hsap.obj@meta.data


CELseq2.HEK.mat.agg.mean.df <- as.data.frame(CELseq2.HEK.mat.agg.mean.f) 
MARSseq.HEK.mat.agg.mean.df <- as.data.frame(MARSseq.HEK.mat.agg.mean.f) 
QUARTZseq.HEK.mat.agg.mean.df <- as.data.frame(QUARTZseq.HEK.mat.agg.mean.f) 
#SCRBseq.HEK.mat.agg.mean.df <- as.data.frame(SCRBseq.HEK.mat.agg.mean.f) 
SMARTseq2.HEK.mat.agg.mean.df <- as.data.frame(SMARTseq2.HEK.mat.agg.mean.f) 
C1HTsmall.HEK.mat.agg.mean.df <- as.data.frame(C1HTsmall.HEK.mat.agg.mean.f) 
C1HTmedium.HEK.mat.agg.mean.df <- as.data.frame(C1HTmedium.HEK.mat.agg.mean.f) 
#Chromium.HEK.mat.agg.mean.df <- as.data.frame(Chromium.HEK.mat.agg.mean.f) 
ChromRef.HEK.mat.agg.mean.df <- as.data.frame(ChromRef.HEK.mat.agg.mean.f) 
ChromNuclei.HEK.mat.agg.mean.df <- as.data.frame(ChromNuclei.HEK.mat.agg.mean.f) 
ddSEQ.HEK.mat.agg.mean.df <- as.data.frame(ddSEQ.HEK.mat.agg.mean.f) 
Dropseq.HEK.mat.agg.mean.df <- as.data.frame(Dropseq.HEK.mat.agg.mean.f) 
ICELL8.HEK.mat.agg.mean.df <- as.data.frame(ICELL8.HEK.mat.agg.mean.f) 

library(dplyr)
all.agg.means.df <- bind_rows(CELseq2.HEK.mat.agg.mean.df, MARSseq.HEK.mat.agg.mean.df, QUARTZseq.HEK.mat.agg.mean.df, 
                              SMARTseq2.HEK.mat.agg.mean.df, C1HTsmall.HEK.mat.agg.mean.df,C1HTmedium.HEK.mat.agg.mean.df, 
                              ChromRef.HEK.mat.agg.mean.df, ChromNuclei.HEK.mat.agg.mean.df,
                              ddSEQ.HEK.mat.agg.mean.df, Dropseq.HEK.mat.agg.mean.df, ICELL8.HEK.mat.agg.mean.df)
all.agg.means.df[is.na(all.agg.means.df)] <- 0
tech_names <- c(rownames(CELseq2.HEK.mat.agg.mean.df), rownames(MARSseq.HEK.mat.agg.mean.df), rownames(QUARTZseq.HEK.mat.agg.mean.df),
                rownames(SMARTseq2.HEK.mat.agg.mean.df), rownames(C1HTsmall.HEK.mat.agg.mean.df), 
                rownames(C1HTmedium.HEK.mat.agg.mean.df), rownames(ChromRef.HEK.mat.agg.mean.df), rownames(ChromNuclei.HEK.mat.agg.mean.df), 
                rownames(ddSEQ.HEK.mat.agg.mean.df), rownames(Dropseq.HEK.mat.agg.mean.df), rownames(ICELL8.HEK.mat.agg.mean.df))
techs <- c(rep("CELseq2", nrow(CELseq2.HEK.mat.agg.mean.df)), rep("MARSseq", nrow(MARSseq.HEK.mat.agg.mean.df)),rep("QUARTZseq", nrow(QUARTZseq.HEK.mat.agg.mean.df)),
           rep("SMARTseq2", nrow(SMARTseq2.HEK.mat.agg.mean.df)),rep("C1HTsmall", nrow(C1HTsmall.HEK.mat.agg.mean.df)),
           rep("C1HTmedium", nrow(C1HTmedium.HEK.mat.agg.mean.df)),rep("ChromRef", nrow(ChromRef.HEK.mat.agg.mean.df)),rep("ChromNuclei", nrow(ChromNuclei.HEK.mat.agg.mean.df)),
           rep("ddSEQ", nrow(ddSEQ.HEK.mat.agg.mean.df)),rep("Dropseq", nrow(Dropseq.HEK.mat.agg.mean.df)),rep("ICELL8", nrow(ICELL8.HEK.mat.agg.mean.df)))

rownames.prep <- paste(techs, tech_names, sep = "_")
rownames(all.agg.means.df) <- rownames.prep

all.techs.rep.cor <- cor(t(all.agg.means.df))

png("/Volumes/Ati-Archive/HCA/SCE_Robjects/revision/replicates/AllReaplicate_Correlation_corrplot_DS20K.png", width = 15, height = 15, units = "in", res = 600)
corrplot(all.techs.rep.cor, method = "square", type = "upper", order = "hclust", tl.cex = 1,tl.col = "black")
dev.off()

my_palette <- colorRampPalette(c("yellow","orange", "green", "#6A1B9A"))(n = 100)
col.colors <- c(rep("red3", table(techs)[3]),rep("green3", table(techs)[9]),rep("gold", table(techs)[10]),
                rep("palevioletred1", table(techs)[11]),rep("darkmagenta", table(techs)[2]),rep("sienna2", table(techs)[1]),rep("darkorange4", table(techs)[5]),
                rep("turquoise3", table(techs)[4]),rep("darkgrey", table(techs)[6]),rep("#117733", table(techs)[7]),rep("maroon2", table(techs)[8]))

save(all.techs.rep.cor, file= "/Volumes/Ati-Archive/HCA/SCE_Robjects/revision/replicates/AllReaplicate_Correlation_heatmap_DS20K_HEK_withChromRef.RData")
png("/Volumes/Ati-Archive/HCA/SCE_Robjects/revision/replicates/AllReaplicate_Correlation_heatmap_DS20K_HEK_withChromRef_2.png", width = 15, height = 15, units = "in", res = 600)
heatmap.2(all.techs.rep.cor, Rowv = T, Colv= T, dendrogram = "both", trace = "none", col = my_palette, cexRow = 1, cexCol = 1, margins = c(12, 12), ColSideColors = col.colors, 
          RowSideColors = col.colors,key = T,hclustfun = function(x) hclust(x,method = 'ward.D2'), distfun = function(x) dist(x,method = 'euclidean')) #, colsep = 1:ncol(Dropseq.rep.cor), rowsep = 1:nrow(Dropseq.rep.cor))
dev.off()


########################## Bcells #######################################################################

### Correlation plots
load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/celseq_seu.obj.RData")
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/Final/CELseq2_clustering_DS20K_dim1to8_res06_OBJ.RData")
comm.CELseq2.cells <- intersect(colnames(CELseq2.DS.obj@raw.data), rownames(celseq@meta.data))
CELseq2.annotation.metadata <- celseq@meta.data[comm.CELseq2.cells,]; rm(celseq)
CELseq2.B.mat <- as.data.frame(t(CELseq2.DS.obj@raw.data[,rownames(CELseq2.annotation.metadata[which(CELseq2.annotation.metadata$clean.id == "B"),])]))
CELseq2.replicates <- CELseq2.annotation.metadata[rownames(CELseq2.B.mat), "Library"]
CELseq2.B.mat.agg <- cbind(CELseq2.B.mat, CELseq2.replicates)
CELseq2.B.mat.agg.mean <- aggregate(CELseq2.B.mat.agg[,1:ncol(CELseq2.B.mat.agg)-1], by=list(CELseq2.replicates), FUN = mean)
CELseq2.B.mat.agg.mean.f <- as.matrix(CELseq2.B.mat.agg.mean)
rownames(CELseq2.B.mat.agg.mean.f) <- CELseq2.B.mat.agg.mean.f[,1]
CELseq2.B.mat.agg.mean.f <- CELseq2.B.mat.agg.mean.f[,-1]
class(CELseq2.B.mat.agg.mean.f) <- "numeric"
CELseq2.rep.cor <- cor(t(CELseq2.B.mat.agg.mean.f))

CELseq2.pca.data <- prcomp(CELseq2.B.mat, center=TRUE) #compute PCA representation of the data
CELseq2.batch.pca <- pcRegression(CELseq2.pca.data, CELseq2.replicates)
PC.regressions.all["CEL-Seq2"] <- CELseq2.batch.pca$pcRegscale


load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/marsseq_seu.obj.RData")
marsseq
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/Final/MARSseq_clustering_DS20K_dim1to8_res06_OBJ.RData")
comm.MARSseq.cells <- intersect(colnames(MARSseq.DS.obj@raw.data), rownames(marsseq@meta.data))
MARSseq.annotation.metadata <- marsseq@meta.data[comm.MARSseq.cells,]; rm(marsseq)
MARSseq.B.mat <- as.data.frame(t(MARSseq.DS.obj@raw.data[,rownames(MARSseq.annotation.metadata[which(MARSseq.annotation.metadata$clean.id == "B"),])]))
MARSseq.replicates <- MARSseq.annotation.metadata[rownames(MARSseq.B.mat), "Library"]
MARSseq.B.mat.agg <- cbind(MARSseq.B.mat, MARSseq.replicates)
MARSseq.B.mat.agg.mean <- aggregate(MARSseq.B.mat.agg[,1:ncol(MARSseq.B.mat.agg)-1], by=list(MARSseq.replicates), FUN = mean)
MARSseq.B.mat.agg.mean.f <- as.matrix(MARSseq.B.mat.agg.mean)
rownames(MARSseq.B.mat.agg.mean.f) <- MARSseq.B.mat.agg.mean.f[,1]
MARSseq.B.mat.agg.mean.f <- MARSseq.B.mat.agg.mean.f[,-1]
class(MARSseq.B.mat.agg.mean.f) <- "numeric"
MARSseq.rep.cor <- cor(t(MARSseq.B.mat.agg.mean.f))

MARSseq.pca.data <- prcomp(MARSseq.B.mat, center=TRUE) #compute PCA representation of the data
MARSseq.batch.pca <- pcRegression(MARSseq.pca.data, MARSseq.replicates)
PC.regressions.all["MARS-Seq"] <- MARSseq.batch.pca$pcRegscale


load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/quartzseq_seu.obj.RData")
quartzseq
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/Final/QUARTZseq_clustering_DS20K_dim1to8_res06_OBJ.RData")
comm.QUARTZseq.cells <- intersect(colnames(QUARTZseq.DS.obj@raw.data), rownames(quartzseq@meta.data))
QUARTZseq.annotation.metadata <- quartzseq@meta.data[comm.QUARTZseq.cells,]; rm(quartzseq)
QUARTZseq.B.mat <- as.data.frame(t(QUARTZseq.DS.obj@raw.data[,rownames(QUARTZseq.annotation.metadata[which(QUARTZseq.annotation.metadata$clean.id == "B"),])]))
QUARTZseq.replicates <- QUARTZseq.annotation.metadata[rownames(QUARTZseq.B.mat), "Library"]

QUARTZseq.B.mat.agg <- cbind(QUARTZseq.B.mat, QUARTZseq.replicates)
QUARTZseq.B.mat.agg.mean <- aggregate(QUARTZseq.B.mat.agg[,1:ncol(QUARTZseq.B.mat.agg)-1], by=list(QUARTZseq.replicates), FUN = mean)
QUARTZseq.B.mat.agg.mean.f <- as.matrix(QUARTZseq.B.mat.agg.mean)
rownames(QUARTZseq.B.mat.agg.mean.f) <- QUARTZseq.B.mat.agg.mean.f[,1]
QUARTZseq.B.mat.agg.mean.f <- QUARTZseq.B.mat.agg.mean.f[,-1]
class(QUARTZseq.B.mat.agg.mean.f) <- "numeric"
QUARTZseq.rep.cor <- cor(t(QUARTZseq.B.mat.agg.mean.f))

QUARTZseq.pca.data <- prcomp(QUARTZseq.B.mat, center=TRUE) #compute PCA representation of the data
QUARTZseq.batch.pca <- pcRegression(QUARTZseq.pca.data, QUARTZseq.replicates)
PC.regressions.all["Quartz-Seq2"] <- QUARTZseq.batch.pca$pcRegscale


load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/scrbseq_seu.obj.RData")
scrbseq
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/Final/SCRBseq_clustering_DS20K_dim1to8_res06_OBJ.RData")
comm.SCRBseq.cells <- intersect(colnames(SCRBseq.DS.obj@raw.data), rownames(scrbseq@meta.data))
SCRBseq.annotation.metadata <- scrbseq@meta.data[comm.SCRBseq.cells,]; rm(scrbseq)
SCRBseq.B.mat <- as.data.frame(t(SCRBseq.DS.obj@raw.data[,rownames(SCRBseq.annotation.metadata[which(SCRBseq.annotation.metadata$clean.id == "B"),])]))
SCRBseq.replicates <- SCRBseq.annotation.metadata[rownames(SCRBseq.B.mat), "Library"]

SCRBseq.B.mat.agg <- cbind(SCRBseq.B.mat, SCRBseq.replicates)
SCRBseq.B.mat.agg.mean <- aggregate(SCRBseq.B.mat.agg[,1:ncol(SCRBseq.B.mat.agg)-1], by=list(SCRBseq.replicates), FUN = mean)
SCRBseq.B.mat.agg.mean.f <- as.matrix(SCRBseq.B.mat.agg.mean)
rownames(SCRBseq.B.mat.agg.mean.f) <- SCRBseq.B.mat.agg.mean.f[,1]
SCRBseq.B.mat.agg.mean.f <- SCRBseq.B.mat.agg.mean.f[,-1]
class(SCRBseq.B.mat.agg.mean.f) <- "numeric"
SCRBseq.rep.cor <- cor(t(SCRBseq.B.mat.agg.mean.f))
#png("/Volumes/Ati-Archive/HCA/SCE_Robjects/revision/replicates/SCRBseq_replicates_corrplot.png", width = 6, height = 6, units = "in", res = 600)
#corrplot(SCRBseq.rep.cor, method = "square", type = "upper")
#dev.off()
SCRBseq.pca.data <- prcomp(SCRBseq.B.mat, center=TRUE) #compute PCA representation of the data
SCRBseq.batch.pca <- pcRegression(SCRBseq.pca.data, SCRBseq.replicates)
PC.regressions.all["mcSCRB-Seq"] <- SCRBseq.batch.pca$pcRegscale



load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/smartseq_seu.obj.RData")
SS2.annotation <- read.table("/Volumes/Ati-Archive/HCA/SCE_Robjects/SMARTseq2/Smartseq2_HCA_annotation.txt", header = T)
cell.BCs <- paste("SMARTseqFINAL_allLanes_", SS2.annotation$Index1_sequence, SS2.annotation$Index2_sequence_HiSeq2500, sep = "")
BCs.plate <- cbind(cell.BCs, SS2.annotation$HCA_plate)
BCs.plate.batch <- c()
for (i in 1:nrow(BCs.plate)){
  if (BCs.plate[i,2] == "2"){BCs.plate.batch[i] <- "batch1"}
  if (BCs.plate[i,2] == "3" | BCs.plate[i,2] == "4"){BCs.plate.batch[i] <- "batch2"}
  if (BCs.plate[i,2] == "5" | BCs.plate[i,2] == "6"){BCs.plate.batch[i] <- "batch3"}
  if (BCs.plate[i,2] == "7" | BCs.plate[i,2] == "8"){BCs.plate.batch[i] <- "batch4"}
  if (BCs.plate[i,2] == "9"){BCs.plate.batch[i] <- "batch5"}
}
SS2.batches <- cbind(BCs.plate, BCs.plate.batch)
rownames(SS2.batches) <- SS2.batches[, "cell.BCs"]
SS2.batches <- as.factor(SS2.batches[,"BCs.plate.batch"])
smartseq
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/Final/SMARTseq2_clustering_DS20K_dim1to8_res06_OBJ.RData")
comm.SMARTseq2.cells <- intersect(colnames(SMARTseq2.DS.obj@raw.data), rownames(smartseq@meta.data))
SMARTseq2.annotation.metadata <- smartseq@meta.data[comm.SMARTseq2.cells,]; rm(smartseq)
SMARTseq2.B.mat <- as.data.frame(t(SMARTseq2.DS.obj@raw.data[,rownames(SMARTseq2.annotation.metadata[which(SMARTseq2.annotation.metadata$clean.id == "B"),])]))
SMARTseq2.replicates <- SS2.batches[rownames(SMARTseq2.B.mat)]

SMARTseq2.B.mat.agg <- cbind(SMARTseq2.B.mat, SMARTseq2.replicates)
SMARTseq2.B.mat.agg.mean <- aggregate(SMARTseq2.B.mat.agg[,1:ncol(SMARTseq2.B.mat.agg)-1], by=list(SMARTseq2.replicates), FUN = mean)
SMARTseq2.B.mat.agg.mean.f <- as.matrix(SMARTseq2.B.mat.agg.mean)
rownames(SMARTseq2.B.mat.agg.mean.f) <- SMARTseq2.B.mat.agg.mean.f[,1]
SMARTseq2.B.mat.agg.mean.f <- SMARTseq2.B.mat.agg.mean.f[,-1]
class(SMARTseq2.B.mat.agg.mean.f) <- "numeric"
SMARTseq2.rep.cor <- cor(t(SMARTseq2.B.mat.agg.mean.f))

SMARTseq2.pca.data <- prcomp(SMARTseq2.B.mat, center=TRUE) #compute PCA representation of the data
SMARTseq2.batch.pca <- pcRegression(SMARTseq2.pca.data, SMARTseq2.replicates)
PC.regressions.all["Smart-Seq2"] <- SMARTseq2.batch.pca$pcRegscale


load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/c1ht.s_seu.obj.RData")
c1ht.s
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/Final/C1HTsmall_clustering_DS20K_dim1to8_res06_OBJ.RData")
comm.C1HTsmall.cells <- intersect(colnames(C1HTsmall.DS.obj@raw.data), rownames(c1ht.s@meta.data))
C1HTsmall.annotation.metadata <- c1ht.s@meta.data[comm.C1HTsmall.cells,]; rm(c1ht.s)
C1HTsmall.B.mat <- as.data.frame(t(C1HTsmall.DS.obj@raw.data[,rownames(C1HTsmall.annotation.metadata[which(C1HTsmall.annotation.metadata$clean.id == "B"),])]))
C1HTsmall.replicates <- as.character(C1HTsmall.annotation.metadata[rownames(C1HTsmall.B.mat), "Library"])
for (i in 1:length(C1HTsmall.replicates)){
  id = C1HTsmall.replicates[i]
  id_parts = unlist(strsplit(id, split ="_"))
  if (length(id_parts) == 4 & grepl("Col", id_parts[3])){
    C1HTsmall.replicates[i] <- id_parts[1]
  }
}

C1HTsmall.B.mat.agg <- cbind(C1HTsmall.B.mat, C1HTsmall.replicates)
C1HTsmall.B.mat.agg.mean <- aggregate(C1HTsmall.B.mat.agg[,1:ncol(C1HTsmall.B.mat.agg)-1], by=list(C1HTsmall.replicates), FUN = mean)
C1HTsmall.B.mat.agg.mean.f <- as.matrix(C1HTsmall.B.mat.agg.mean)
rownames(C1HTsmall.B.mat.agg.mean.f) <- C1HTsmall.B.mat.agg.mean.f[,1]
C1HTsmall.B.mat.agg.mean.f <- C1HTsmall.B.mat.agg.mean.f[,-1]
class(C1HTsmall.B.mat.agg.mean.f) <- "numeric"
C1HTsmall.rep.cor <- cor(t(C1HTsmall.B.mat.agg.mean.f))

C1HTsmall.pca.data <- prcomp(C1HTsmall.B.mat, center=TRUE) #compute PCA representation of the data
C1HTsmall.batch.pca <- pcRegression(C1HTsmall.pca.data, C1HTsmall.replicates)
PC.regressions.all["C1HT-small"] <- C1HTsmall.batch.pca$pcRegscale


load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/c1ht.m_seu.obj.RData")
c1ht.m
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/Final/C1HTmedium_clustering_DS20K_dim1to8_res06_OBJ.RData")
comm.C1HTmedium.cells <- intersect(colnames(C1HTmedium.DS.obj@raw.data), rownames(c1ht.m@meta.data))
C1HTmedium.annotation.metadata <- c1ht.m@meta.data[comm.C1HTmedium.cells,]; rm(c1ht.m)
C1HTmedium.B.mat <- as.data.frame(t(C1HTmedium.DS.obj@raw.data[,rownames(C1HTmedium.annotation.metadata[which(C1HTmedium.annotation.metadata$clean.id == "B"),])]))
C1HTmedium.replicates <- as.character(C1HTmedium.annotation.metadata[rownames(C1HTmedium.B.mat), "Library"])
for (i in 1:length(C1HTmedium.replicates)){
  id = C1HTmedium.replicates[i]
  id_parts = unlist(strsplit(id, split ="_"))
  if (length(id_parts) == 4 & grepl("Col", id_parts[3])){
    C1HTmedium.replicates[i] <- id_parts[1]
  }
}

C1HTmedium.B.mat.agg <- cbind(C1HTmedium.B.mat, C1HTmedium.replicates)
C1HTmedium.B.mat.agg.mean <- aggregate(C1HTmedium.B.mat.agg[,1:ncol(C1HTmedium.B.mat.agg)-1], by=list(C1HTmedium.replicates), FUN = mean)
C1HTmedium.B.mat.agg.mean.f <- as.matrix(C1HTmedium.B.mat.agg.mean)
rownames(C1HTmedium.B.mat.agg.mean.f) <- C1HTmedium.B.mat.agg.mean.f[,1]
C1HTmedium.B.mat.agg.mean.f <- C1HTmedium.B.mat.agg.mean.f[,-1]
class(C1HTmedium.B.mat.agg.mean.f) <- "numeric"
C1HTmedium.rep.cor <- cor(t(C1HTmedium.B.mat.agg.mean.f))

C1HTmedium.pca.data <- prcomp(C1HTmedium.B.mat, center=TRUE) #compute PCA representation of the data
C1HTmedium.batch.pca <- pcRegression(C1HTmedium.pca.data, C1HTmedium.replicates)
PC.regressions.all["C1HT-medium"] <- C1HTmedium.batch.pca$pcRegscale


library(rowr)

dim(ChromV2.all.DS20K) # From the Chrom_Comparison_rev.R code
dim(ChromV2.all.DS20K.B) # From the Chrom_Comparison_rev.R code
dim(ChromV2.all.metadata) # From the Chrom_Comparison_rev.R code
Chromium.B.mat <- ChromV2.all.DS20K.B
Chromium.B.mat <- mapIDs(Chromium.B.mat, "hsap")
Chromium.B.mat <- t(Chromium.B.mat)
Chromium.replicates <- ChromV2.all.metadata[rownames(Chromium.B.mat), "Library"]

Chromium.B.mat.agg <- cbind(Chromium.B.mat, Chromium.replicates)
Chromium.B.mat.agg.mean <- aggregate(Chromium.B.mat.agg[,1:ncol(Chromium.B.mat.agg)-1], by=list(Chromium.replicates), FUN = mean)
Chromium.B.mat.agg.mean.f <- as.matrix(Chromium.B.mat.agg.mean)
rownames(Chromium.B.mat.agg.mean.f) <- Chromium.B.mat.agg.mean.f[,1]
Chromium.B.mat.agg.mean.f <- Chromium.B.mat.agg.mean.f[,-1]
class(Chromium.B.mat.agg.mean.f) <- "numeric"
Chromium.rep.cor <- cor(t(Chromium.B.mat.agg.mean.f))

Chromium.pca.data <- prcomp(Chromium.B.mat, center=TRUE) #compute PCA representation of the data
Chromium.batch.pca <- pcRegression(Chromium.pca.data, Chromium.replicates)
PC.regressions.all["Chromium"] <- Chromium.batch.pca$pcRegscale


load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/nuclei_seu.obj.RData")
nuclei
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/Final/ChromiumNuclei_clustering_DS20K_dim1to8_res06_OBJ.RData")
comm.ChromNuclei.cells <- intersect(colnames(ChromiumNuclei.DS.obj@raw.data), rownames(nuclei@meta.data))
ChromNuclei.annotation.metadata <- nuclei@meta.data[comm.ChromNuclei.cells,]; rm(nuclei)
ChromNuclei.B.mat <- as.data.frame(t(ChromiumNuclei.DS.obj@raw.data[,rownames(ChromNuclei.annotation.metadata[which(ChromNuclei.annotation.metadata$clean.id == "B"),])]))
ChromNuclei.replicates <- ChromNuclei.annotation.metadata[rownames(ChromNuclei.B.mat), "Library"]

ChromNuclei.B.mat.agg <- cbind(ChromNuclei.B.mat, ChromNuclei.replicates)
ChromNuclei.B.mat.agg.mean <- aggregate(ChromNuclei.B.mat.agg[,1:ncol(ChromNuclei.B.mat.agg)-1], by=list(ChromNuclei.replicates), FUN = mean)
ChromNuclei.B.mat.agg.mean.f <- as.matrix(ChromNuclei.B.mat.agg.mean)
rownames(ChromNuclei.B.mat.agg.mean.f) <- ChromNuclei.B.mat.agg.mean.f[,1]
ChromNuclei.B.mat.agg.mean.f <- ChromNuclei.B.mat.agg.mean.f[,-1]
class(ChromNuclei.B.mat.agg.mean.f) <- "numeric"
ChromNuclei.rep.cor <- cor(t(ChromNuclei.B.mat.agg.mean.f))

ChromNuclei.pca.data <- prcomp(ChromNuclei.B.mat, center=TRUE) #compute PCA representation of the data
ChromNuclei.batch.pca <- pcRegression(ChromNuclei.pca.data, ChromNuclei.replicates)
PC.regressions.all["Chromium(sn)"] <- ChromNuclei.batch.pca$pcRegscale


load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/ddseq_seu.obj.RData")
ddseq
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/Final/ddSEQ_clustering_DS20K_dim1to8_res06_OBJ.RData")
comm.ddSEQ.cells <- intersect(colnames(ddSEQ.DS.obj@raw.data), rownames(ddseq@meta.data))
ddSEQ.annotation.metadata <- ddseq@meta.data[comm.ddSEQ.cells,]; rm(ddseq)
ddSEQ.B.mat <- as.data.frame(t(ddSEQ.DS.obj@raw.data[,rownames(ddSEQ.annotation.metadata[which(ddSEQ.annotation.metadata$clean.id == "B"),])]))
ddSEQ.replicates <- ddSEQ.annotation.metadata[rownames(ddSEQ.B.mat), "Library"]

ddSEQ.B.mat.agg <- cbind(ddSEQ.B.mat, ddSEQ.replicates)
ddSEQ.B.mat.agg.mean <- aggregate(ddSEQ.B.mat.agg[,1:ncol(ddSEQ.B.mat.agg)-1], by=list(ddSEQ.replicates), FUN = mean)
ddSEQ.B.mat.agg.mean.f <- as.matrix(ddSEQ.B.mat.agg.mean)
rownames(ddSEQ.B.mat.agg.mean.f) <- ddSEQ.B.mat.agg.mean.f[,1]
ddSEQ.B.mat.agg.mean.f <- ddSEQ.B.mat.agg.mean.f[,-1]
class(ddSEQ.B.mat.agg.mean.f) <- "numeric"
ddSEQ.rep.cor <- cor(t(ddSEQ.B.mat.agg.mean.f))

ddSEQ.pca.data <- prcomp(ddSEQ.B.mat, center=TRUE) #compute PCA representation of the data
ddSEQ.batch.pca <- pcRegression(ddSEQ.pca.data, ddSEQ.replicates)
PC.regressions.all["ddSEQ"] <- ddSEQ.batch.pca$pcRegscale


load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/dropseq_seu.obj.RData")
dropseq
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/Final/Dropseq_clustering_DS20K_dim1to8_res06_OBJ.RData")
comm.Dropseq.cells <- intersect(colnames(Dropseq.DS.obj@raw.data), rownames(dropseq@meta.data))
Dropseq.annotation.metadata <- dropseq@meta.data[comm.Dropseq.cells,]; rm(dropseq)
Dropseq.B.mat <- as.data.frame(t(Dropseq.DS.obj@raw.data[,rownames(Dropseq.annotation.metadata[which(Dropseq.annotation.metadata$clean.id == "B"),])]))
Dropseq.replicates <- Dropseq.annotation.metadata[rownames(Dropseq.B.mat), "Library"]

Dropseq.B.mat.agg <- cbind(Dropseq.B.mat, Dropseq.replicates)
Dropseq.B.mat.agg.mean <- aggregate(Dropseq.B.mat.agg[,1:ncol(Dropseq.B.mat.agg)-1], by=list(Dropseq.replicates), FUN = mean)
Dropseq.B.mat.agg.mean.f <- as.matrix(Dropseq.B.mat.agg.mean)
rownames(Dropseq.B.mat.agg.mean.f) <- Dropseq.B.mat.agg.mean.f[,1]
Dropseq.B.mat.agg.mean.f <- Dropseq.B.mat.agg.mean.f[,-1]
class(Dropseq.B.mat.agg.mean.f) <- "numeric"
Dropseq.rep.cor <- cor(t(Dropseq.B.mat.agg.mean.f))

Dropseq.pca.data <- prcomp(Dropseq.B.mat, center=TRUE) #compute PCA representation of the data
Dropseq.batch.pca <- pcRegression(Dropseq.pca.data, Dropseq.replicates)
PC.regressions.all["Drop-Seq"] <- Dropseq.batch.pca$pcRegscale


load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/icell8_seu.obj.RData")
icell8
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/modularity/Final/ICELL8_clustering_DS20K_dim1to8_res06_OBJ.RData")
comm.ICELL8.cells <- intersect(colnames(ICELL8.DS.obj@raw.data), rownames(icell8@meta.data))
ICELL8.annotation.metadata <- icell8@meta.data[comm.ICELL8.cells,]; rm(icell8)
ICELL8.B.mat <- as.data.frame(t(ICELL8.DS.obj@raw.data[,rownames(ICELL8.annotation.metadata[which(ICELL8.annotation.metadata$clean.id == "B"),])]))
ICELL8.replicates <- ICELL8.annotation.metadata[rownames(ICELL8.B.mat), "Library"]

ICELL8.B.mat.agg <- cbind(ICELL8.B.mat, ICELL8.replicates) #till here is the prev version
ICELL8.B.mat.agg.mean <- aggregate(ICELL8.B.mat.agg[,1:ncol(ICELL8.B.mat.agg)-1], by=list(ICELL8.replicates), FUN = mean)
ICELL8.B.mat.agg.mean.f <- as.matrix(ICELL8.B.mat.agg.mean)
rownames(ICELL8.B.mat.agg.mean.f) <- ICELL8.B.mat.agg.mean.f[,1]
ICELL8.B.mat.agg.mean.f <- ICELL8.B.mat.agg.mean.f[,-1]
class(ICELL8.B.mat.agg.mean.f) <- "numeric"
ICELL8.rep.cor <- cor(t(ICELL8.B.mat.agg.mean.f))

ICELL8.pca.data <- prcomp(ICELL8.B.mat, center=TRUE) #compute PCA representation of the data
ICELL8.batch.pca <- pcRegression(ICELL8.pca.data, ICELL8.replicates)
PC.regressions.all["ICELL8"] <- ICELL8.batch.pca$pcRegscale


load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/onecb_seu.obj.RData")
inDrop.hsap.obj <- onecb
rm(onecb)
inDrop.hsap.metadata <- inDrop.hsap.obj@meta.data



CELseq2.B.mat.agg.mean.df <- as.data.frame(CELseq2.B.mat.agg.mean.f) 
MARSseq.B.mat.agg.mean.df <- as.data.frame(MARSseq.B.mat.agg.mean.f) 
QUARTZseq.B.mat.agg.mean.df <- as.data.frame(QUARTZseq.B.mat.agg.mean.f) 
SCRBseq.B.mat.agg.mean.df <- as.data.frame(SCRBseq.B.mat.agg.mean.f) 
SMARTseq2.B.mat.agg.mean.df <- as.data.frame(SMARTseq2.B.mat.agg.mean.f) 
C1HTsmall.B.mat.agg.mean.df <- as.data.frame(C1HTsmall.B.mat.agg.mean.f) 
C1HTmedium.B.mat.agg.mean.df <- as.data.frame(C1HTmedium.B.mat.agg.mean.f) 
#Chromium.B.mat.agg.mean.df <- as.data.frame(Chromium.B.mat.agg.mean.f) 
ChromNuclei.B.mat.agg.mean.df <- as.data.frame(ChromNuclei.B.mat.agg.mean.f) 
ddSEQ.B.mat.agg.mean.df <- as.data.frame(ddSEQ.B.mat.agg.mean.f) 
Dropseq.B.mat.agg.mean.df <- as.data.frame(Dropseq.B.mat.agg.mean.f) 
ICELL8.B.mat.agg.mean.df <- as.data.frame(ICELL8.B.mat.agg.mean.f) 

library(dplyr)
all.agg.means.df <- bind_rows(CELseq2.B.mat.agg.mean.df, MARSseq.B.mat.agg.mean.df, QUARTZseq.B.mat.agg.mean.df, 
                              SCRBseq.B.mat.agg.mean.df, SMARTseq2.B.mat.agg.mean.df, C1HTsmall.B.mat.agg.mean.df,
                              C1HTmedium.B.mat.agg.mean.df, ChromNuclei.B.mat.agg.mean.df,
                              ddSEQ.B.mat.agg.mean.df, Dropseq.B.mat.agg.mean.df, ICELL8.B.mat.agg.mean.df)
all.agg.means.df[is.na(all.agg.means.df)] <- 0
tech_names <- c(rownames(CELseq2.B.mat.agg.mean.df), rownames(MARSseq.B.mat.agg.mean.df), rownames(QUARTZseq.B.mat.agg.mean.df),
                rownames(SCRBseq.B.mat.agg.mean.df), rownames(SMARTseq2.B.mat.agg.mean.df), rownames(C1HTsmall.B.mat.agg.mean.df), 
                rownames(C1HTmedium.B.mat.agg.mean.df), rownames(ChromNuclei.B.mat.agg.mean.df), 
                rownames(ddSEQ.B.mat.agg.mean.df), rownames(Dropseq.B.mat.agg.mean.df), rownames(ICELL8.B.mat.agg.mean.df))
techs <- c(rep("CELseq2", nrow(CELseq2.B.mat.agg.mean.df)), rep("MARSseq", nrow(MARSseq.B.mat.agg.mean.df)),rep("QUARTZseq", nrow(QUARTZseq.B.mat.agg.mean.df)),
           rep("SCRBseq", nrow(SCRBseq.B.mat.agg.mean.df)),rep("SMARTseq2", nrow(SMARTseq2.B.mat.agg.mean.df)),rep("C1HTsmall", nrow(C1HTsmall.B.mat.agg.mean.df)),
           rep("C1HTmedium", nrow(C1HTmedium.B.mat.agg.mean.df)),rep("ChromNuclei", nrow(ChromNuclei.B.mat.agg.mean.df)),
           rep("ddSEQ", nrow(ddSEQ.B.mat.agg.mean.df)),rep("Dropseq", nrow(Dropseq.B.mat.agg.mean.df)),rep("ICELL8", nrow(ICELL8.B.mat.agg.mean.df)))

rownames.prep <- paste(techs, tech_names, sep = "_")
rownames(all.agg.means.df) <- rownames.prep

all.techs.rep.cor <- cor(t(all.agg.means.df))

png("/Volumes/Ati-Archive/HCA/SCE_Robjects/revision/replicates/AllReaplicate_Correlation_corrplot_DS20K.png", width = 15, height = 15, units = "in", res = 600)
corrplot(all.techs.rep.cor, method = "square", type = "upper", order = "hclust", tl.cex = 1,tl.col = "black")
dev.off()

my_palette <- colorRampPalette(c("yellow","orange", "green", "#6A1B9A"))(n = 100)
col.colors <- c(rep("red3", table(techs)[3]),rep("green3", table(techs)[8]),rep("gold", table(techs)[9]),rep("blue1", table(techs)[10]),
                rep("palevioletred1", table(techs)[11]),rep("darkmagenta", table(techs)[2]),rep("sienna2", table(techs)[1]),
                rep("turquoise3", table(techs)[4]),rep("darkgrey", table(techs)[5]),rep("#117733", table(techs)[6]),rep("maroon2", table(techs)[7]))

save(all.techs.rep.cor, file = "/Volumes/Ati-Archive/HCA/SCE_Robjects/revision/replicates/AllReaplicate_Correlation_heatmap_DS20K_Bcells.RData")
png("/Volumes/Ati-Archive/HCA/SCE_Robjects/revision/replicates/AllReaplicate_Correlation_heatmap_DS20K_Bcells_withoutkey_2.png", width = 15, height = 15, units = "in", res = 600)
heatmap.2(all.techs.rep.cor, Rowv = T, Colv= T, dendrogram = "both", trace = "none", col = my_palette, cexRow = 1, cexCol = 1, margins = c(12, 12), ColSideColors = col.colors, 
          RowSideColors = col.colors,key = F ,hclustfun = function(x) hclust(x,method = 'ward.D2'), distfun = function(x) dist(x,method = 'euclidean')) #, colsep = 1:ncol(Dropseq.rep.cor), rowsep = 1:nrow(Dropseq.rep.cor))
dev.off()
