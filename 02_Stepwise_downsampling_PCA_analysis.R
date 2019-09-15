Robjects.DS <- c("CELseq2.hsap.full.SCE.jointDSmat.Robj","MARSseq.hsap.full.SCE.jointDSmat.Robj","QUARTZseq.hsap.full.SCE.jointDSmat.Robj","SCRBseq.hsap.full.SCE.jointDSmat.Robj",
              "SMARTseqFINAL.hsap.full.SCE.jointDSmat.Robj", "C1HTsmall.hsap.full.SCE.jointDSmat.Robj","C1HTmedium.hsap.full.SCE.jointDSmat.Robj","10X2x5K.hsap.full.SCE.jointDSmat.Robj", "Nuclei10X.hsap.full.SCE.jointDSmat.Robj",
              "ddSEQ.hsap.full.SCE.jointDSmat.Robj", "Dropseq.hsap.full.SCE.1000cellsPerPool.Robj", "ICELL8.hsap.full.SCE.jointDSmat.Robj", "1CB.hsap.full.SCE.AllReads.Robj")
Seurat_objs.hsap <- c("celseq_seu.obj.RData","marsseq_seu.obj.RData", "quartzseq_seu.obj.RData", "scrbseq_seu.obj.RData",
                      "smartseq_seu.obj.RData", "c1ht.s_seu.obj.RData", "c1ht.m_seu.obj.RData", "chromium_seu.obj.RData", "nuclei_seu.obj.RData",
                      "ddseq_seu.obj.RData", "dropseq_seu.obj.RData", "icell8_seu.obj.RData", "onecb_seu.obj.RData")

Robjects.DS = c("/Volumes/Ati-Archive/HCA/SCE_Robjects/Downsampled/DS_secondRound/Human/CELseq2.hsap.full.SCE.jointDSmat.Robj", "/Volumes/Ati-Archive/HCA/SCE_Robjects/Downsampled/DS_secondRound/Human/MARSseq.hsap.full.SCE.jointDSmat.Robj", "/Volumes/Ati-Archive/HCA/SCE_Robjects/Downsampled/DS_secondRound/Human/QUARTZseq.hsap.full.SCE.jointDSmat.Robj")
Seurat_objs.hsap = c("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/celseq_seu.obj.RData", "/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/marsseq_seu.obj.RData", "/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/quartzseq_seu.obj.RData")
#Seurat_objs.mmus = c("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Mouse/celseq_seu.obj.RData", "/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Mouse/marsseq_seu.obj.RData", "/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Mouse/quartzseq_seu.obj.RData")


for(i in 1:length(Robjects.DS)){
  load(Robjects.DS[i])
  protocol.name= unlist(strsplit(Robjects.DS[i], split = "/"))[9]
  protocol.name= unlist(strsplit(protocol.name, split = "\\."))[1]
  #protocol.name= unlist(strsplit(Robjects.DS[i], split = "."))[1]
  print(protocol.name)
  protocol.DS <- output.readcount.umicount.joint.mats
  rm(output.readcount.umicount.joint.mats)
  protocol.DS.UMI <- protocol.DS$UMI
  if (protocol.name == "SMARTseqFINAL"){
    protocol.name <- "SMARTseq2"
    protocol.DS.UMI <- protocol.DS$Reads
  }
  DS.data.name <- paste(protocol.name,"DS.UMI", sep = ".")
  assign(DS.data.name, protocol.DS.UMI)
  rm(protocol.DS)
  rm(protocol.DS.UMI)
  
  load(Seurat_objs.hsap[i])
  seurat.obj.name.hsap <- unlist(strsplit(Seurat_objs.hsap[i],split = "/"))[8]
  seurat.obj.name.hsap <- unlist(strsplit(seurat.obj.name.hsap,split = "_"))[1]
  #seurat.obj.name.hsap <- unlist(strsplit(Seurat_objs.hsap[i],split = "_"))[1]
  seurat.obj.hsap <- get(seurat.obj.name.hsap) 
  protocol.hsap.metadata <- seurat.obj.hsap@meta.data
  metadata.obj.name <- paste(protocol.name,"hsap.metadata", sep = ".")
  assign(metadata.obj.name, protocol.hsap.metadata)
  rm(seurat.obj.hsap)
  
  # HEK cells
  protocol.HEK <- rownames(protocol.hsap.metadata)[which(protocol.hsap.metadata$clean.id == "HEK")]
  assign(paste(protocol.name,"HEK", sep = "."), protocol.HEK)
  
  # Monocytes 
  protocol.Monocytes <- rownames(protocol.hsap.metadata)[which(protocol.hsap.metadata$clean.id == "Monocytes")]
  assign(paste(protocol.name,"Monocytes", sep = "."), protocol.Monocytes)
  
  # B cells 
  protocol.Bcells <- rownames(protocol.hsap.metadata)[which(protocol.hsap.metadata$clean.id == "B")]
  assign(paste(protocol.name,"Bcells", sep = "."), protocol.Bcells)
  
}




techniques <- c("CELseq2","MARSseq","QUARTZseq","SCRBseq","SMARTseq2", "C1HTsmall","C1HTmedium","Chromium", "ChromiumNuclei","ddSEQ","Dropseq","ICELL8","inDrop")

# HEK
DSth.df.HEK <- data.frame()
techs.HEK.20K.list <- list()
for (tech in techniques){
  print(tech)
  tech.DS.UMI <- get(paste(tech,".DS.UMI", sep = ""))
  tech.HEK <- get(paste(tech,".HEK", sep = ""))
  for (DSth in names(tech.DS.UMI)){
    print(paste(tech, DSth, sep = "_"))
    colnames(tech.DS.UMI[[DSth]]) <- gsub(x = colnames(tech.DS.UMI[[DSth]]), pattern = "\\.", replacement = "_")
    comm.cells <- intersect(tech.HEK, colnames(tech.DS.UMI[[DSth]]))
    DS.mat.HEKS <- tech.DS.UMI[[DSth]][, comm.cells]
    if (DSth == "downsampled_20000"){
      DS.mat.HEKS.dup <- DS.mat.HEKS
      DS.mat.HEKS.dup$gene_id <- rownames(DS.mat.HEKS.dup)
      techs.HEK.20K.list[[tech]] <- as.data.frame((DS.mat.HEKS.dup))}
    DS.gene.distribution.HEK <- colSums(DS.mat.HEKS[,]>0)
    DS.UMI.distribution.HEK <- colSums(DS.mat.HEKS)
    DS.labels <- rep(DSth, length(DS.gene.distribution.HEK))
    DSth_number = as.numeric(unlist(strsplit(DSth, "_"))[[2]])
    DSth.number.vec = rep(DSth_number, length(DS.gene.distribution.HEK))
    DS.labels <- rep(DSth, length(DS.gene.distribution.HEK))
    DS.techs <- rep(tech, length(DS.gene.distribution.HEK))
    DS.df <- data.frame(nGenes = DS.gene.distribution.HEK, nUMIs = DS.UMI.distribution.HEK, DSthNum = DSth.number.vec, DStech = DS.techs)
    DSth.df.HEK <- rbind(DSth.df.HEK, DS.df)
    print(dim(DSth.df.HEK))
  }
}

save(DSth.df.HEK, file="/Stepwise_DS_plotdf_HEK.RData")


# Monocytes
DSth.df.Monocytes <- data.frame()
techs.Monocytes.20K.list <- list()
for (tech in techniques){
  print(tech)
  tech.DS.UMI <- get(paste(tech,".DS.UMI", sep = ""))
  tech.Monocytes <- get(paste(tech,".Monocytes", sep = ""))
  for (DSth in names(tech.DS.UMI)){
    print(paste(tech, DSth, sep = "_"))
    colnames(tech.DS.UMI[[DSth]]) <- gsub(x = colnames(tech.DS.UMI[[DSth]]), pattern = "\\.", replacement = "_")
    comm.cells <- intersect(tech.Monocytes, colnames(tech.DS.UMI[[DSth]]))
    DS.mat.MonocytesS <- tech.DS.UMI[[DSth]][, comm.cells]
    if (DSth == "downsampled_20000"){
      DS.mat.MonocytesS.dup <- DS.mat.MonocytesS
      DS.mat.MonocytesS.dup$gene_id <- rownames(DS.mat.MonocytesS.dup)
      techs.Monocytes.20K.list[[tech]] <- as.data.frame((DS.mat.MonocytesS.dup))}
    DS.gene.distribution.Monocytes <- colSums(DS.mat.MonocytesS[,]>0)
    DS.UMI.distribution.Monocytes <- colSums(DS.mat.MonocytesS)
    DS.labels <- rep(DSth, length(DS.gene.distribution.Monocytes))
    DSth_number = as.numeric(unlist(strsplit(DSth, "_"))[[2]])
    DSth.number.vec = rep(DSth_number, length(DS.gene.distribution.Monocytes))
    DS.labels <- rep(DSth, length(DS.gene.distribution.Monocytes))
    DS.techs <- rep(tech, length(DS.gene.distribution.Monocytes))
    DS.df <- data.frame(nGenes = DS.gene.distribution.Monocytes, nUMIs = DS.UMI.distribution.Monocytes, DSthNum = DSth.number.vec, DStech = DS.techs)
    DSth.df.Monocytes <- rbind(DSth.df.Monocytes, DS.df)
    print(dim(DSth.df.Monocytes))
  }
}

save(DSth.df.Monocytes, file="Stepwise_DS_plotdf_Monocytes.RData")


#B cells
DSth.df.Bcells <- data.frame()
techs.Bcells.20K.list <- list()
for (tech in techniques){
  print(tech)
  tech.DS.UMI <- get(paste(tech,".DS.UMI", sep = ""))
  tech.Bcells <- get(paste(tech,".Bcells", sep = ""))
  for (DSth in names(tech.DS.UMI)){
    print(paste(tech, DSth, sep = "_"))
    colnames(tech.DS.UMI[[DSth]]) <- gsub(x = colnames(tech.DS.UMI[[DSth]]), pattern = "\\.", replacement = "_")
    comm.cells <- intersect(tech.Bcells, colnames(tech.DS.UMI[[DSth]]))
    DS.mat.BcellsS <- tech.DS.UMI[[DSth]][, comm.cells]
    if (DSth == "downsampled_20000"){
      DS.mat.BcellsS.dup <- DS.mat.BcellsS
      DS.mat.BcellsS.dup$gene_id <- rownames(DS.mat.BcellsS.dup)
      techs.Bcells.20K.list[[tech]] <- as.data.frame((DS.mat.BcellsS.dup))}
    DS.gene.distribution.Bcells <- colSums(DS.mat.BcellsS[,]>0)
    DS.UMI.distribution.Bcells <- colSums(DS.mat.BcellsS)
    DS.labels <- rep(DSth, length(DS.gene.distribution.Bcells))
    DSth_number = as.numeric(unlist(strsplit(DSth, "_"))[[2]])
    DSth.number.vec = rep(DSth_number, length(DS.gene.distribution.Bcells))
    DS.labels <- rep(DSth, length(DS.gene.distribution.Bcells))
    DS.techs <- rep(tech, length(DS.gene.distribution.Bcells))
    DS.df <- data.frame(nGenes = DS.gene.distribution.Bcells, nUMIs = DS.UMI.distribution.Bcells, DSthNum = DSth.number.vec, DStech = DS.techs)
    DSth.df.Bcells <- rbind(DSth.df.Bcells, DS.df)
    print(dim(DSth.df.Bcells))
  }
}

save(DSth.df.Bcells, file="Stepwise_DS_plotdf_Bcells.RData")





# The PCA for celltypes HEK
library(plyr)
techs.HEK.20K.df=join_all(techs.HEK.20K.list, by = "gene_id", type = 'full')
rownames(techs.HEK.20K.df) <- techs.HEK.20K.df[,"gene_id"]
gene.id.col <- which(colnames(techs.HEK.20K.df) == "gene_id")
techs.HEK.20K.df <- techs.HEK.20K.df[,-gene.id.col]
dim(techs.HEK.20K.df)
techs.HEK.20K.df <- na.omit(techs.HEK.20K.df)
HEK.seurat <- CreateSeuratObject(techs.HEK.20K.df, min.cells = 5, min.genes = 1, normalization.method = "LogNormalize", scale.factor = 10000)#, total.expr = 1e4, project = "Dicroce")
HEK.seurat <- ScaleData(object= HEK.seurat)
HEK.seurat <- FindVariableGenes(object = HEK.seurat, mean.function = ExpMean, dispersion.function = LogVMR, 
                                x.low.cutoff = 0.2, x.high.cutoff = 6, y.cutoff = 0.5)
HEK.seurat <- RunPCA(HEK.seurat, pc.genes = HEK.seurat@var.genes, do.print = F)
HEK.colors.nUMI = HEK.seurat@meta.data[names(HEK.seurat@ident), "nUMI"]
names(HEK.colors.nUMI) <- names(HEK.seurat@ident)
HEK.data.plot <- as.data.frame(HEK.seurat@dr$pca@cell.embeddings[,1:8])
HEK.data.plot <- cbind(HEK.data.plot, nUMI=HEK.colors.nUMI, techs=HEK.seurat@meta.data$orig.ident)

save(HEK.data.plot, file="all_techs_stepwise_DS20K_PCAseurat_HEK_dataPlot.RData")
save(HEK.seurat, file="all_techs_stepwise_DS20K_PCAseurat_HEK_SeuratObj.RData")



# The PCA for celltypes Monocytes
techs.Monocytes.20K.df=join_all(techs.Monocytes.20K.list, by = "gene_id", type = 'full')
rownames(techs.Monocytes.20K.df) <- techs.Monocytes.20K.df[,"gene_id"]
gene.id.col <- which(colnames(techs.Monocytes.20K.df) == "gene_id")
techs.Monocytes.20K.df <- techs.Monocytes.20K.df[,-gene.id.col]
dim(techs.Monocytes.20K.df)
techs.Monocytes.20K.df <- na.omit(techs.Monocytes.20K.df)
Monocytes.seurat <- CreateSeuratObject(techs.Monocytes.20K.df, min.cells = 5, min.genes = 1, normalization.method = "LogNormalize", scale.factor = 10000)#, total.expr = 1e4, project = "Dicroce")
Monocytes.seurat <- ScaleData(object= Monocytes.seurat)
Monocytes.seurat <- FindVariableGenes(object = Monocytes.seurat, mean.function = ExpMean, dispersion.function = LogVMR, 
                                      x.low.cutoff = 0.2, x.high.cutoff = 6, y.cutoff = 0.5)
Monocytes.seurat <- RunPCA(Monocytes.seurat, pc.genes = Monocytes.seurat@var.genes, do.print = F)
Monocytes.colors.nUMI = Monocytes.seurat@meta.data[names(Monocytes.seurat@ident), "nUMI"]
names(Monocytes.colors.nUMI) <- names(Monocytes.seurat@ident)
Monocytes.data.plot <- as.data.frame(Monocytes.seurat@dr$pca@cell.embeddings[,1:8])
Monocytes.data.plot <- cbind(Monocytes.data.plot, nUMI=Monocytes.colors.nUMI, techs=Monocytes.seurat@meta.data$orig.ident)

save(Monocytes.data.plot, file="all_techs_stepwise_DS20K_PCAseurat_Monocytes_dataPlot.RData")
save(Monocytes.seurat, file="all_techs_stepwise_DS20K_PCAseurat_Monocytes_SeuratObj.RData")


# The PCA for celltypes Bcells
techs.Bcells.20K.df=join_all(techs.Bcells.20K.list, by = "gene_id", type = 'full')
rownames(techs.Bcells.20K.df) <- techs.Bcells.20K.df[,"gene_id"]
gene.id.col <- which(colnames(techs.Bcells.20K.df) == "gene_id")
techs.Bcells.20K.df <- techs.Bcells.20K.df[,-gene.id.col]
dim(techs.Bcells.20K.df)
techs.Bcells.20K.df <- na.omit(techs.Bcells.20K.df)
Bcells.seurat <- CreateSeuratObject(techs.Bcells.20K.df, min.cells = 5, min.genes = 1, normalization.method = "LogNormalize", scale.factor = 10000)#, total.expr = 1e4, project = "Dicroce")
Bcells.seurat <- ScaleData(object= Bcells.seurat)
Bcells.seurat <- FindVariableGenes(object = Bcells.seurat, mean.function = ExpMean, dispersion.function = LogVMR, 
                                   x.low.cutoff = 0.2, x.high.cutoff = 6, y.cutoff = 0.5)
Bcells.seurat <- RunPCA(Bcells.seurat, pc.genes = Bcells.seurat@var.genes, do.print = F)
Bcells.colors.nUMI = Bcells.seurat@meta.data[names(Bcells.seurat@ident), "nUMI"]
names(Bcells.colors.nUMI) <- names(Bcells.seurat@ident)
Bcells.data.plot <- as.data.frame(Bcells.seurat@dr$pca@cell.embeddings[,1:8])
Bcells.data.plot <- cbind(Bcells.data.plot, nUMI=Bcells.colors.nUMI, techs=Bcells.seurat@meta.data$orig.ident)

save(Bcells.data.plot, file="all_techs_stepwise_DS20K_PCAseurat_Bcells_dataPlot.RData")
save(Bcells.seurat, file="all_techs_stepwise_DS20K_PCAseurat_Bcells_SeuratObj.RData")


  
  