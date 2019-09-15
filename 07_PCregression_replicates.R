library(kBET)
library(data.table)
library(Seurat)

load('Biomart_mmus_mapping_table.RData')
load('Biomart_hsap_mapping_table.RData')
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


### PC Regressions
set.seed(123)
PC.regressions.all <- c(rep(NA,11))
names(PC.regressions.all) <- c("CEL-Seq2", "MARS-Seq", "Quartz-Seq2", "Smart-Seq2", "C1HT-small", "C1HT-medium", "Chromium-Ref","Chromium(sn)", "ddSEQ", "Drop-Seq", "ICELL8")
PC.R2.df <- data.frame(tech = NA, PCs=NA, R2vals=NA)

load("Human/celseq_seu.obj.RData")
celseq <- UpdateSeuratObject(celseq)
CELseq2.annotation.metadata <- celseq@meta.data
CELseq2.all.mat <- t(as.matrix(celseq@assays$RNA@counts))
CELseq2.cell.embeddings <- celseq@reductions$pca@cell.embeddings
CELseq2.sdev <- celseq@reductions$pca@stdev
CELseq2.replicates <- CELseq2.annotation.metadata[rownames(CELseq2.cell.embeddings), "Library"]
names(CELseq2.replicates) <- rownames(CELseq2.cell.embeddings)
CELseq2.cells.per.replicate <- table(CELseq2.replicates)
CELseq2.selected.replicates <- names(CELseq2.cells.per.replicate[which(CELseq2.cells.per.replicate >= 5)])
CELseq2.all.mat.selectedReps <- CELseq2.all.mat[names(CELseq2.replicates[which(CELseq2.replicates %in% CELseq2.selected.replicates)]),]
#CELseq2.pca.data <- prcomp(CELseq2.all.mat.selectedReps, center=TRUE) #compute PCA representation of the data
#CELseq2.batch.pca <- pcRegression(CELseq2.pca.data, CELseq2.replicates[rownames(CELseq2.all.mat.selectedReps)])
#PC.regressions.all["CEL-Seq2"] <- CELseq2.batch.pca$pcRegscale
CELseq2.pcRegression.input <- list(x=CELseq2.cell.embeddings, sdev=CELseq2.sdev)
CELseq2.batch.pca <- pcRegression(CELseq2.pcRegression.input, CELseq2.replicates[rownames(CELseq2.all.mat.selectedReps)])
CELseq2.PC.R2.df <- data.frame(tech = rep("CEL-Seq2", 20), PCs = rownames(CELseq2.batch.pca$r2)[1:20], R2vals= CELseq2.batch.pca$r2[1:20,"R.squared"])
PC.R2.df <- rbind(PC.R2.df,CELseq2.PC.R2.df)
#save(CELseq2.batch.pca, file= "PC_regression/CELseq2.PC.regression.RData")


load("Human/marsseq_seu.obj.RData")
marsseq <- UpdateSeuratObject(marsseq)
MARSseq.annotation.metadata <- marsseq@meta.data
MARSseq.all.mat <- t(as.matrix(marsseq@assays$RNA@counts))
MARSseq.cell.embeddings <- marsseq@reductions$pca@cell.embeddings
MARSseq.sdev <- marsseq@reductions$pca@stdev
MARSseq.replicates <- MARSseq.annotation.metadata[rownames(MARSseq.cell.embeddings), "Library"]
names(MARSseq.replicates) <- rownames(MARSseq.cell.embeddings)
MARSseq.cells.per.replicate <- table(MARSseq.replicates)
MARSseq.selected.replicates <- names(MARSseq.cells.per.replicate[which(MARSseq.cells.per.replicate >= 5)])
MARSseq.all.mat.selectedReps <- MARSseq.all.mat[names(MARSseq.replicates[which(MARSseq.replicates %in% MARSseq.selected.replicates)]),]
#MARSseq.pca.data <- prcomp(MARSseq.all.mat.selectedReps, center=TRUE) #compute PCA representation of the data
#MARSseq.batch.pca <- pcRegression(MARSseq.pca.data, MARSseq.replicates[rownames(MARSseq.all.mat.selectedReps)])
#PC.regressions.all["MARS-Seq"] <- MARSseq.batch.pca$pcRegscale
MARSseq.pcRegression.input <- list(x=MARSseq.cell.embeddings, sdev=MARSseq.sdev)
MARSseq.batch.pca <- pcRegression(MARSseq.pcRegression.input, MARSseq.replicates[rownames(MARSseq.all.mat.selectedReps)])
MARSseq.PC.R2.df <- data.frame(tech = rep("MARS-Seq", 20), PCs = rownames(MARSseq.batch.pca$r2)[1:20], R2vals= MARSseq.batch.pca$r2[1:20,"R.squared"])
PC.R2.df <- rbind(PC.R2.df,MARSseq.PC.R2.df)
#save(MARSseq.batch.pca, file= "PC_regression/MARSseq.PC.regression.RData")



load("Human/QUARTZseq_seu.obj.RData")
quartzseq <- UpdateSeuratObject(quartzseq)
QUARTZseq.annotation.metadata <- quartzseq@meta.data
QUARTZseq.all.mat <- t(as.matrix(quartzseq@assays$RNA@counts))
QUARTZseq.cell.embeddings <- quartzseq@reductions$pca@cell.embeddings
QUARTZseq.sdev <- quartzseq@reductions$pca@stdev
QUARTZseq.replicates <- QUARTZseq.annotation.metadata[rownames(QUARTZseq.cell.embeddings), "Library"]
names(QUARTZseq.replicates) <- rownames(QUARTZseq.cell.embeddings)
QUARTZseq.cells.per.replicate <- table(QUARTZseq.replicates)
QUARTZseq.selected.replicates <- names(QUARTZseq.cells.per.replicate[which(QUARTZseq.cells.per.replicate >= 5)])
QUARTZseq.all.mat.selectedReps <- QUARTZseq.all.mat[names(QUARTZseq.replicates[which(QUARTZseq.replicates %in% QUARTZseq.selected.replicates)]),]
#QUARTZseq.pca.data <- prcomp(QUARTZseq.all.mat.selectedReps, center=TRUE) #compute PCA representation of the data
#QUARTZseq.batch.pca <- pcRegression(QUARTZseq.pca.data, QUARTZseq.replicates[rownames(QUARTZseq.all.mat.selectedReps)])
QUARTZseq.pcRegression.input <- list(x=QUARTZseq.cell.embeddings, sdev=QUARTZseq.sdev)
QUARTZseq.batch.pca <- pcRegression(QUARTZseq.pcRegression.input, QUARTZseq.replicates[rownames(QUARTZseq.all.mat.selectedReps)])
PC.regressions.all["Quartz-Seq2"] <- QUARTZseq.batch.pca$pcRegscale
QUARTZseq.PC.R2.df <- data.frame(tech = rep("Quartz-Seq2", 20), PCs = rownames(QUARTZseq.batch.pca$r2)[1:20], R2vals= QUARTZseq.batch.pca$r2[1:20,"R.squared"])
PC.R2.df <- rbind(PC.R2.df,QUARTZseq.PC.R2.df)
#save(QUARTZseq.batch.pca, file= "PC_regression/QUARTZseq.PC.regression.RData")

load("Human/scrbseq_seu.obj.RData")
scrbseq <- UpdateSeuratObject(scrbseq)
SCRBseq.annotation.metadata <- scrbseq@meta.data
SCRBseq.all.mat <- t(as.matrix(scrbseq@assays$RNA@counts))
SCRBseq.cell.embeddings <- scrbseq@reductions$pca@cell.embeddings
SCRBseq.sdev <- scrbseq@reductions$pca@stdev
SCRBseq.replicates <- SCRBseq.annotation.metadata[rownames(SCRBseq.cell.embeddings), "Library"]
names(SCRBseq.replicates) <- rownames(SCRBseq.cell.embeddings)
SCRBseq.cells.per.replicate <- table(SCRBseq.replicates)
SCRBseq.selected.replicates <- names(SCRBseq.cells.per.replicate[which(SCRBseq.cells.per.replicate >= 5)])
SCRBseq.all.mat.selectedReps <- SCRBseq.all.mat[names(SCRBseq.replicates[which(SCRBseq.replicates %in% SCRBseq.selected.replicates)]),]
#SCRBseq.pca.data <- prcomp(SCRBseq.all.mat.selectedReps, center=TRUE) #compute PCA representation of the data
#SCRBseq.batch.pca <- pcRegression(SCRBseq.pca.data, SCRBseq.replicates[rownames(SCRBseq.all.mat.selectedReps)])
SCRBseq.pcRegression.input <- list(x=SCRBseq.cell.embeddings, sdev=SCRBseq.sdev)
SCRBseq.batch.pca <- pcRegression(SCRBseq.pcRegression.input, SCRBseq.replicates[rownames(SCRBseq.all.mat.selectedReps)])
PC.regressions.all["mcSCRB-Seq"] <- SCRBseq.batch.pca$pcRegscale
SCRBseq.PC.R2.df <- data.frame(tech = rep("mcSCRB-Seq", 20), PCs = rownames(SCRBseq.batch.pca$r2)[1:20], R2vals= SCRBseq.batch.pca$r2[1:20,"R.squared"])
PC.R2.df <- rbind(PC.R2.df,SCRBseq.PC.R2.df)
#save(SCRBseq.batch.pca, file= "PC_regression/SCRBseq.PC.regression.RData")


load("Human/smartseq_seu.obj.RData")
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
smartseq <- UpdateSeuratObject(smartseq)
SMARTseq2.annotation.metadata <- smartseq@meta.data
SMARTseq2.all.mat <- t(as.matrix(smartseq@assays$RNA@counts))
SMARTseq2.cell.embeddings <- smartseq@reductions$pca@cell.embeddings
SMARTseq2.sdev <- smartseq@reductions$pca@stdev
SMARTseq2.replicates <- SS2.batches[rownames(SMARTseq2.cell.embeddings)]
#names(SMARTseq2.replicates) <- rownames(SMARTseq2.cell.embeddings)
SMARTseq2.cells.per.replicate <- table(SMARTseq2.replicates)
SMARTseq2.selected.replicates <- names(SMARTseq2.cells.per.replicate[which(SMARTseq2.cells.per.replicate >= 5)])
SMARTseq2.all.mat.selectedReps <- SMARTseq2.all.mat[names(SMARTseq2.replicates[which(SMARTseq2.replicates %in% SMARTseq2.selected.replicates)]),]
#SMARTseq2.pca.data <- prcomp(SMARTseq2.all.mat.selectedReps, center=TRUE) #compute PCA representation of the data
#SMARTseq2.batch.pca <- pcRegression(SMARTseq2.pca.data, SMARTseq2.replicates[rownames(SMARTseq2.all.mat.selectedReps)])
SMARTseq2.pcRegression.input <- list(x=SMARTseq2.cell.embeddings[rownames(SMARTseq2.all.mat.selectedReps),], sdev=SMARTseq2.sdev)
SMARTseq2.batch.pca <- pcRegression(SMARTseq2.pcRegression.input, SMARTseq2.replicates[rownames(SMARTseq2.all.mat.selectedReps)])
PC.regressions.all["Smart-Seq2"] <- SMARTseq2.batch.pca$pcRegscale
SMARTseq2.PC.R2.df <- data.frame(tech = rep("Smart-Seq2", 20), PCs = rownames(SMARTseq2.batch.pca$r2)[1:20], R2vals= SMARTseq2.batch.pca$r2[1:20,"R.squared"])
PC.R2.df <- rbind(PC.R2.df,SMARTseq2.PC.R2.df)
#save(SMARTseq2.batch.pca, file= "PC_regression/SMARTseq2.PC.regression.RData")

load("Human/c1ht.s_seu.obj.RData")
c1ht.s <- UpdateSeuratObject(c1ht.s)
C1HTsmall.annotation.metadata <- c1ht.s@meta.data
C1HTsmall.all.mat <- t(as.matrix(c1ht.s@assays$RNA@counts))
C1HTsmall.cell.embeddings <- c1ht.s@reductions$pca@cell.embeddings
C1HTsmall.sdev <- c1ht.s@reductions$pca@stdev
C1HTsmall.replicates <- as.character(C1HTsmall.annotation.metadata[rownames(C1HTsmall.cell.embeddings), "Library"])
names(C1HTsmall.replicates) <- rownames(C1HTsmall.cell.embeddings)
for (i in 1:length(C1HTsmall.replicates)){
  id = C1HTsmall.replicates[i]
  id_parts = unlist(strsplit(id, split ="_"))
  if (length(id_parts) == 4 & grepl("Col", id_parts[3])){
    C1HTsmall.replicates[i] <- id_parts[1]
  }}
C1HTsmall.cells.per.replicate <- table(C1HTsmall.replicates)
C1HTsmall.selected.replicates <- names(C1HTsmall.cells.per.replicate[which(C1HTsmall.cells.per.replicate >= 5)])
C1HTsmall.all.mat.selectedReps <- C1HTsmall.all.mat[names(C1HTsmall.replicates[which(C1HTsmall.replicates %in% C1HTsmall.selected.replicates)]),]
#C1HTsmall.pca.data <- prcomp(C1HTsmall.all.mat.selectedReps, center=TRUE) #compute PCA representation of the data
#C1HTsmall.batch.pca <- pcRegression(C1HTsmall.pca.data, C1HTsmall.replicates[rownames(C1HTsmall.all.mat.selectedReps)])
C1HTsmall.pcRegression.input <- list(x=C1HTsmall.cell.embeddings, sdev=C1HTsmall.sdev)
C1HTsmall.batch.pca <- pcRegression(C1HTsmall.pcRegression.input, C1HTsmall.replicates[rownames(C1HTsmall.all.mat.selectedReps)])
PC.regressions.all["C1HT-small"] <- C1HTsmall.batch.pca$pcRegscale
C1HTsmall.PC.R2.df <- data.frame(tech = rep("C1HT-small", 20), PCs = rownames(C1HTsmall.batch.pca$r2)[1:20], R2vals= C1HTsmall.batch.pca$r2[1:20,"R.squared"])
PC.R2.df <- rbind(PC.R2.df,C1HTsmall.PC.R2.df)
#save(C1HTsmall.batch.pca, file= "PC_regression/C1HTsmall.PC.regression.RData")


load("Human/c1ht.m_seu.obj.RData")
c1ht.m <- UpdateSeuratObject(c1ht.m)
C1HTmedium.annotation.metadata <- c1ht.m@meta.data
C1HTmedium.all.mat <- t(as.matrix(c1ht.m@assays$RNA@counts))
C1HTmedium.cell.embeddings <- c1ht.m@reductions$pca@cell.embeddings
C1HTmedium.sdev <- c1ht.m@reductions$pca@stdev
C1HTmedium.replicates <- as.character(C1HTmedium.annotation.metadata[rownames(C1HTmedium.cell.embeddings), "Library"])
names(C1HTmedium.replicates) <- rownames(C1HTmedium.cell.embeddings)
for (i in 1:length(C1HTmedium.replicates)){
  id = C1HTmedium.replicates[i]
  id_parts = unlist(strsplit(id, split ="_"))
  if (length(id_parts) == 4 & grepl("Col", id_parts[3])){
    C1HTmedium.replicates[i] <- id_parts[1]
  }}
C1HTmedium.cells.per.replicate <- table(C1HTmedium.replicates)
C1HTmedium.selected.replicates <- names(C1HTmedium.cells.per.replicate[which(C1HTmedium.cells.per.replicate >= 5)])
C1HTmedium.all.mat.selectedReps <- C1HTmedium.all.mat[names(C1HTmedium.replicates[which(C1HTmedium.replicates %in% C1HTmedium.selected.replicates)]),]
#C1HTmedium.pca.data <- prcomp(C1HTmedium.all.mat.selectedReps, center=TRUE) #compute PCA representation of the data
#C1HTmedium.batch.pca <- pcRegression(C1HTmedium.pca.data, C1HTmedium.replicates[rownames(C1HTmedium.all.mat.selectedReps)])
C1HTmedium.pcRegression.input <- list(x=C1HTmedium.cell.embeddings, sdev=C1HTmedium.sdev)
C1HTmedium.batch.pca <- pcRegression(C1HTmedium.pcRegression.input, C1HTmedium.replicates[rownames(C1HTmedium.all.mat.selectedReps)])
PC.regressions.all["C1HT-medium"] <- C1HTmedium.batch.pca$pcRegscale
C1HTmedium.PC.R2.df <- data.frame(tech = rep("C1HT-medium", 20), PCs = rownames(C1HTmedium.batch.pca$r2)[1:20], R2vals= C1HTmedium.batch.pca$r2[1:20,"R.squared"])
PC.R2.df <- rbind(PC.R2.df,C1HTmedium.PC.R2.df)
#save(C1HTmedium.batch.pca, file= "PC_regression/C1HTmedium.PC.regression.RData")

#Chromium Ref
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/reference_10x/Ref10X_metadata.RData")
dim(Ref10X.metadata)
Ref10X.replicates <- Ref10X.metadata[, "Library"]
names(Ref10X.replicates) <- rownames(Ref10X.metadata)

load("PC_regression/10XRef_seuratObj_PCA.RData")
Ref10X.cell.embedings <- PCA.10XRef[[1]]
colnames(Ref10X.cell.embedings) <- colnames(C1HTmedium.cell.embeddings) # because for 10Xref colnames didnt have _ in the PC names and sinc eits just 20 for all, I got it from previous one
Ref10X.sdev <- PCA.10XRef[[2]]
Ref10X.pcRegression.input <- list(x= Ref10X.cell.embedings, sdev=Ref10X.sdev)
Ref10X.batch.pca <- pcRegression(Ref10X.pcRegression.input, Ref10X.replicates[rownames(Ref10X.cell.embedings)])
Ref10X.PC.R2.df <- data.frame(tech = rep("Chromium-Ref", 20), PCs = rownames(Ref10X.batch.pca$r2)[1:20], R2vals= Ref10X.batch.pca$r2[1:20,"R.squared"])
PC.R2.df <- rbind(PC.R2.df,Ref10X.PC.R2.df)



load("Human/nuclei_seu.obj.RData")
nuclei <- UpdateSeuratObject(nuclei)
ChromNuclei.annotation.metadata <- nuclei@meta.data
ChromNuclei.all.mat <- t(as.matrix(nuclei@assays$RNA@counts))
ChromNuclei.cell.embeddings <- nuclei@reductions$pca@cell.embeddings
ChromNuclei.sdev <- nuclei@reductions$pca@stdev
ChromNuclei.replicates <- ChromNuclei.annotation.metadata[rownames(ChromNuclei.cell.embeddings), "Library"]
names(ChromNuclei.replicates) <- rownames(ChromNuclei.cell.embeddings)
ChromNuclei.cells.per.replicate <- table(ChromNuclei.replicates)
ChromNuclei.selected.replicates <- names(ChromNuclei.cells.per.replicate[which(ChromNuclei.cells.per.replicate >= 5)])
ChromNuclei.all.mat.selectedReps <- ChromNuclei.all.mat[names(ChromNuclei.replicates[which(ChromNuclei.replicates %in% ChromNuclei.selected.replicates)]),]
#ChromNuclei.pca.data <- prcomp(ChromNuclei.all.mat.selectedReps, center=TRUE) #compute PCA representation of the data
#ChromNuclei.batch.pca <- pcRegression(ChromNuclei.pca.data, ChromNuclei.replicates[rownames(ChromNuclei.all.mat.selectedReps)])
ChromNuclei.pcRegression.input <- list(x=ChromNuclei.cell.embeddings, sdev=ChromNuclei.sdev)
ChromNuclei.batch.pca <- pcRegression(ChromNuclei.pcRegression.input, ChromNuclei.replicates[rownames(ChromNuclei.all.mat.selectedReps)])
PC.regressions.all["Chromium(sn)"] <- ChromNuclei.batch.pca$pcRegscale
ChromNuclei.PC.R2.df <- data.frame(tech = rep("Chromium(sn)", 20), PCs = rownames(ChromNuclei.batch.pca$r2)[1:20], R2vals= ChromNuclei.batch.pca$r2[1:20,"R.squared"])
PC.R2.df <- rbind(PC.R2.df,ChromNuclei.PC.R2.df)
#save(ChromNuclei.batch.pca, file= "PC_regression/ChromNuclei.PC.regression.RData")

load("Human/ddseq_seu.obj.RData")
ddseq <- UpdateSeuratObject(ddseq)
ddSEQ.annotation.metadata <- ddseq@meta.data
ddSEQ.all.mat <- t(as.matrix(ddseq@assays$RNA@counts))
ddSEQ.cell.embeddings <- ddseq@reductions$pca@cell.embeddings
ddSEQ.sdev <- ddseq@reductions$pca@stdev
ddSEQ.replicates <- ddSEQ.annotation.metadata[rownames(ddSEQ.cell.embeddings), "Library"]
names(ddSEQ.replicates) <- rownames(ddSEQ.cell.embeddings)
ddSEQ.cells.per.replicate <- table(ddSEQ.replicates)
ddSEQ.selected.replicates <- names(ddSEQ.cells.per.replicate[which(ddSEQ.cells.per.replicate >= 5)])
ddSEQ.all.mat.selectedReps <- ddSEQ.all.mat[names(ddSEQ.replicates[which(ddSEQ.replicates %in% ddSEQ.selected.replicates)]),]
#ddSEQ.pca.data <- prcomp(ddSEQ.all.mat.selectedReps, center=TRUE) #compute PCA representation of the data
#ddSEQ.batch.pca <- pcRegression(ddSEQ.pca.data, ddSEQ.replicates[rownames(ddSEQ.all.mat.selectedReps)])
ddSEQ.pcRegression.input <- list(x=ddSEQ.cell.embeddings, sdev=ddSEQ.sdev)
ddSEQ.batch.pca <- pcRegression(ddSEQ.pcRegression.input, ddSEQ.replicates[rownames(ddSEQ.all.mat.selectedReps)])
PC.regressions.all["ddSEQ"] <- ddSEQ.batch.pca$pcRegscale
ddSEQ.PC.R2.df <- data.frame(tech = rep("ddSEQ", 20), PCs = rownames(ddSEQ.batch.pca$r2)[1:20], R2vals= ddSEQ.batch.pca$r2[1:20,"R.squared"])
PC.R2.df <- rbind(PC.R2.df,ddSEQ.PC.R2.df)
#save(ddSEQ.batch.pca, file= "PC_regression/ddSEQ.PC.regression.RData")


load("Human/dropseq_seu.obj.RData")
dropseq <- UpdateSeuratObject(dropseq)
Dropseq.annotation.metadata <- dropseq@meta.data
Dropseq.all.mat <- t(as.matrix(dropseq@assays$RNA@counts))
Dropseq.cell.embeddings <- dropseq@reductions$pca@cell.embeddings
Dropseq.sdev <- dropseq@reductions$pca@stdev
Dropseq.replicates <- Dropseq.annotation.metadata[rownames(Dropseq.cell.embeddings), "Library"]
names(Dropseq.replicates) <- rownames(Dropseq.cell.embeddings)
Dropseq.cells.per.replicate <- table(Dropseq.replicates)
Dropseq.selected.replicates <- names(Dropseq.cells.per.replicate[which(Dropseq.cells.per.replicate >= 5)])
Dropseq.all.mat.selectedReps <- Dropseq.all.mat[names(Dropseq.replicates[which(Dropseq.replicates %in% Dropseq.selected.replicates)]),]
#Dropseq.pca.data <- prcomp(Dropseq.all.mat.selectedReps, center=TRUE) #compute PCA representation of the data
#Dropseq.batch.pca <- pcRegression(Dropseq.pca.data, Dropseq.replicates[rownames(Dropseq.all.mat.selectedReps)])
Dropseq.pcRegression.input <- list(x=Dropseq.cell.embeddings, sdev=Dropseq.sdev)
Dropseq.batch.pca <- pcRegression(Dropseq.pcRegression.input, Dropseq.replicates[rownames(Dropseq.all.mat.selectedReps)])
PC.regressions.all["Drop-Seq"] <- Dropseq.batch.pca$pcRegscale
Dropseq.PC.R2.df <- data.frame(tech = rep("Drop-Seq", 20), PCs = rownames(Dropseq.batch.pca$r2)[1:20], R2vals= Dropseq.batch.pca$r2[1:20,"R.squared"])
PC.R2.df <- rbind(PC.R2.df,Dropseq.PC.R2.df)
#save(Dropseq.batch.pca, file= "PC_regression/Dropseq.PC.regression.RData")

load("Human/icell8_seu.obj.RData")
icell8 <- UpdateSeuratObject(icell8)
ICELL8.annotation.metadata <- icell8@meta.data
ICELL8.all.mat <- t(as.matrix(icell8@assays$RNA@counts))
ICELL8.cell.embeddings <- icell8@reductions$pca@cell.embeddings
ICELL8.sdev <- icell8@reductions$pca@stdev
ICELL8.replicates <- ICELL8.annotation.metadata[rownames(ICELL8.cell.embeddings), "Library"]
names(ICELL8.replicates) <- rownames(ICELL8.cell.embeddings)
ICELL8.cells.per.replicate <- table(ICELL8.replicates)
ICELL8.selected.replicates <- names(ICELL8.cells.per.replicate[which(ICELL8.cells.per.replicate >= 5)])
ICELL8.all.mat.selectedReps <- ICELL8.all.mat[names(ICELL8.replicates[which(ICELL8.replicates %in% ICELL8.selected.replicates)]),]
#ICELL8.pca.data <- prcomp(ICELL8.all.mat.selectedReps, center=TRUE) #compute PCA representation of the data
#ICELL8.batch.pca <- pcRegression(ICELL8.pca.data, ICELL8.replicates[rownames(ICELL8.all.mat.selectedReps)])
ICELL8.pcRegression.input <- list(x=ICELL8.cell.embeddings, sdev=ICELL8.sdev)
ICELL8.batch.pca <- pcRegression(ICELL8.pcRegression.input, ICELL8.replicates[rownames(ICELL8.all.mat.selectedReps)])
PC.regressions.all["ICELL8"] <- ICELL8.batch.pca$pcRegscale
ICELL8.PC.R2.df <- data.frame(tech = rep("ICELL8", 20), PCs = rownames(ICELL8.batch.pca$r2)[1:20], R2vals= ICELL8.batch.pca$r2[1:20,"R.squared"])
PC.R2.df <- rbind(PC.R2.df,ICELL8.PC.R2.df)
#save(ICELL8.batch.pca, file= "PC_regression/ICELL8.PC.regression.RData")


PC.regressions.all.df <- melt(PC.regressions.all[-c(13)])
PC.regressions.all.df <- cbind(PC.regressions.all.df, tech = rownames(PC.regressions.all.df))
PC.regressions.all.df$tech <- factor(PC.regressions.all.df$tech, levels = rownames(PC.regressions.all.df))
PC.regressions.all.df.log <- PC.regressions.all.df
#PC.regressions.all.df.log$value <- log(PC.regressions.all.df.log$value + 1)
library(ppls)
ddd <- cbind(PC.regressions.all.df$value,PC.regressions.all.df$value)
rownames(ddd) <- rownames(PC.regressions.all.df)
my_palette.2 <- colorRampPalette(c("turquoise", "royalblue1","yellow2","springgreen","sienna2", "hotpink"))(n = 100)
#rc <- rainbow(nrow(ddd), start=0, end=.8)
color.vec = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
png("AllReaplicate_PCRegression_heatmap_2.png", width = 15, height = 15, units = "in", res = 600)
heatmap.2(ddd, Rowv = T, Colv = F, trace = "none", col = my_palette.2)
dev.off()

PC.regressions.all.df.log$value <- log(PC.regressions.all.df.log$value + 1)
png("AllReaplicate_PCRegression_barplot_2.png", width = 15, height = 15, units = "in", res = 600)
ggplot(data=PC.regressions.all.df.log, aes(x = tech, y = value, fill = tech)) + geom_bar(stat = "identity") +
  theme (axis.text.x = element_text(angle = 50, hjust = 1, size = 12, face = "bold"),axis.text.y = element_text(size = 12, face = "bold"), legend.position = "none",
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
         axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=selected.color.test)# + scale_y_continuous(trans = "log") #c("red3","green3","gold","blue1","palevioletred1","darkmagenta","sienna2","turquoise3","darkgrey","#117733","maroon2")
dev.off()

PC.R2.df <- PC.R2.df[-1,]
load("techniques_colors_2.RData")
names(selected.color)[8] <- "Chromium-Ref"
selected.color <- selected.color[-c(13)]
names(selected.color)[7] <- "Chromium-Ref"
library(ggplot2)
PC.R2.df$PCs <- factor(PC.R2.df$PCs, levels = PC.R2.df$PCs[1:20])
PC.R2.df$tech <- factor(PC.R2.df$tech, levels = c("CEL-Seq2","MARS-Seq","Quartz-Seq2","mcSCRB-Seq","Smart-Seq2","C1HT-small","C1HT-medium","Chromium-Ref","Chromium(sn)","ddSEQ","Drop-Seq","ICELL8"))
png("AllReaplicate_PCRegression_PCvsR2_allcells_HVGs.png", width = 7, height = 8, units = "in", res = 600)
ggplot(PC.R2.df, aes(x=PCs, y=R2vals, color = tech, group=tech)) + geom_point()+geom_line(size=1)  +facet_grid(tech ~ .) +  scale_color_manual(values=selected.color) +
  theme (axis.text.x = element_text(angle = 90, hjust = 0.5,vjust = 0.3, size = 20, face = "bold", colour = "black"),axis.text.y = element_text(size = 7, face = "bold", colour = "black"), legend.position = "none", strip.background = element_blank(), strip.text = element_blank()) +
  ylab("R.squared")
dev.off()


PC.R2.df.test <- PC.R2.df
PC.R2.df.test$PCs <- factor(PC.R2.df.test$PCs, levels = PC.R2.df.test$PCs[1:20])
selected.color.4 <- c(selected.color, "darkorange4")
names(selected.color.4)[11] <- "Chromium-Ref"
png("AllReaplicate_PCRegression_PCvsR2_allcells.png", width = 7, height = 8, units = "in", res = 600)
ggplot(PC.R2.df.test, aes(x=PCs, y=R2vals, color = tech, group=tech)) + geom_point()+geom_line()  +facet_grid(tech ~ .) +  scale_color_manual(values=selected.color.4) +
  theme (axis.text.x = element_text(angle = 90, hjust = 0.5,vjust = 0.3, size = 10, face = "bold", colour = "black"),axis.text.y = element_text(size = 6, face = "bold", colour = "black"), legend.position = "none", strip.background = element_blank(), strip.text = element_blank()) +
  ylab("R.squared")
dev.off()
