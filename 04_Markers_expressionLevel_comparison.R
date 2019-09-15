
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


Robjects <- c("CELseq2.hsap.full.SCE.Robj","MARSseq.hsap.full.SCE.Robj","QUARTZseq.hsap.full.SCE.Robj","SCRBseq.hsap.full.SCE.Robj",
              "SMARTseqFINAL.hsap.full.SCE.Robj", "C1HTsmall.hsap.full.SCE.Robj","C1HTmedium.hsap.full.SCE.Robj","10X2x5K.hsap.full.SCE.Robj", "Nuclei10X.hsap.full.SCE.Robj",
              "ddSEQ.hsap.full.SCE.Robj", "Dropseq.hsap.full.SCE.1000cellsPerPool.Robj", "ICELL8.hsap.full.SCE.Robj", "1CB.hsap.full.SCE.AllReads.Robj")
Seurat_objs.hsap <- c("celseq_seu.obj.RData","marsseq_seu.obj.RData", "quartzseq_seu.obj.RData", "scrbseq_seu.obj.RData",
                      "smartseq_seu.obj.RData", "c1ht.s_seu.obj.RData", "c1ht.m_seu.obj.RData", "chromium_seu.obj.RData", "nuclei_seu.obj.RData",
                      "ddseq_seu.obj.RData", "dropseq_seu.obj.RData", "icell8_seu.obj.RData", "onecb_seu.obj.RData")
Seurat_objs.mmus <- c("celseq_seu.obj.mmus.RData","marsseq_seu.obj.mmus.RData", "quartzseq_seu.obj.mmus.RData", "scrbseq_seu.obj.mmus.RData",
                      "smartseq_seu.obj.mmus.RData", "c1ht.s_seu.obj.mmus.RData", "c1ht.m_seu.obj.mmus.RData", "chromium_seu.obj.mmus.RData", "nuclei_seu.obj.mmus.RData",
                      "ddseq_seu.obj.mmus.RData", "dropseq_seu.obj.mmus.RData", "icell8_seu.obj.mmus.RData", "onecb_seu.obj.mmus.RData")

Robjects.DS = c("/Volumes/Ati-Archive/HCA/SCE_Robjects/Downsampled/DS_secondRound/Human/CELseq2.hsap.full.SCE.jointDSmat.Robj", "/Volumes/Ati-Archive/HCA/SCE_Robjects/Downsampled/DS_secondRound/Human/MARSseq.hsap.full.SCE.jointDSmat.Robj", "/Volumes/Ati-Archive/HCA/SCE_Robjects/Downsampled/DS_secondRound/Human/QUARTZseq.hsap.full.SCE.jointDSmat.Robj")
Seurat_objs.hsap = c("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/celseq_seu.obj.RData", "/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/marsseq_seu.obj.RData", "/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/quartzseq_seu.obj.RData")

set.seed(123)
for(i in 1:length(Robjects.DS)){
  load(Robjects.DS[i])
  protocol.name= unlist(strsplit(Robjects.DS[i], split = "/"))[9]
  protocol.name= unlist(strsplit(protocol.name, split = "\\."))[1]
  #protocol.name= unlist(strsplit(Robjects.DS[i], split = "."))[1]
  print(protocol.name)
  protocol.DS.20K <- output.readcount.umicount.joint.mats$UMI$downsampled_20000
  if (protocol.name == "SMARTseqFINAL"){
    protocol.name <- "SMARTseq2"
    protocol.DS.20K <- output.readcount.umicount.joint.mats$Reads$downsampled_20000
  }
  rm(output.readcount.umicount.joint.mats)
  DS.data.name <- paste(protocol.name,"DS.20K", sep = ".")
  assign(DS.data.name, protocol.DS.20K)
  rm(protocol.DS.20K)
  
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
  if(length(protocol.HEK) > 50){ protocol.HEK <- sample(protocol.HEK,50)}
  assign(paste(protocol.name,"HEK", sep = "."), protocol.HEK)
  
  # Monocytes 
  protocol.Monocytes <- rownames(protocol.hsap.metadata)[which(protocol.hsap.metadata$clean.id == "Monocytes")]
  if(length(protocol.Monocytes) > 50){ protocol.Monocytes <- sample(protocol.Monocytes,50)}
  assign(paste(protocol.name,"Monocytes", sep = "."), protocol.Monocytes)
  
  # B cells 
  protocol.Bcells <- rownames(protocol.hsap.metadata)[which(protocol.hsap.metadata$clean.id == "B")]
  if(length(protocol.Bcells) > 50){ protocol.Bcells <- sample(protocol.Bcells,50)}
  assign(paste(protocol.name,"Bcells", sep = "."), protocol.Bcells)
}


load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/gene_cl.ref3.RData")
ref.markers <- gene_cl.ref

# preparing  cells downsampleded matrix and plot
techniques <- c("CELseq2","MARSseq","QUARTZseq","SCRBseq","SMARTseq2", "C1HTsmall","C1HTmedium","Chromium", "ChromiumNuclei","ddSEQ","Dropseq","ICELL8","inDrop")
for (protocol in techniques){
  protocol.20K <- get(paste(protocol,"DS.20K", sep = "."))
  protocol.20K <- mapIDs(protocol.20K, "hsap")
  colnames(protocol.20K) <- gsub(x = colnames(protocol.20K), pattern = "\\.", replacement = "_")
  protocol.Bcells.common <- intersect(colnames(protocol.20K),protocol.Bcells)
  protocol.Monocytes.common <- intersect(colnames(protocol.20K),protocol.Monocytes)
  protocol.HEK.common <- intersect(colnames(protocol.20K),protocol.HEK)
  assign(paste(protocol,"Bcells.common", sep = "."), protocol.Bcells.common)
  assign(paste(protocol,"Monocytes.common", sep = "."), protocol.Monocytes.common)
  assign(paste(protocol,"HEK.common", sep = "."), protocol.HEK.common)
}


# Bcells heatmap
Bcells.common.markers <- Reduce(intersect, list(ref.markers[["B cells"]], rownames(CELseq2.20K), rownames(MARSseq.20K), rownames(QUARTZseq.20K), rownames(SCRBseq.20K), 
                                                rownames(SMARTseq2.20K) , rownames(C1HTsmall.20K), rownames(C1HTmedium.20K), rownames(Chromium.20K), rownames(ChromiumNuclei.20K),
                                                rownames(ddSEQ.20K) , rownames(Dropseq.20K), rownames(ICELL8.20K), rownames(inDrop.20K)))
print(paste("length common Bcell markers = ", length(Bcells.common), sep=""))

CELseq2.Bcell.mat <- CELseq2.20K[Bcells.common.markers, CELseq2.Bcells.common]
MARSseq.Bcell.mat <- MARSseq.20K[Bcells.common.markers, MARSseq.Bcells.common]
QUARTZseq.Bcell.mat <- QUARTZseq.20K[Bcells.common.markers, QUARTZseq.Bcells.common]
SCRBseq.Bcell.mat <- SCRBseq.20K[Bcells.common.markers, SCRBseq.Bcells.common]
SMARTseq2.Bcell.mat <- SMARTseq2.20K[Bcells.common.markers, SMARTseq2.Bcells.common]
C1HTsmall.Bcell.mat <- C1HTsmall.20K[Bcells.common.markers, C1HTsmall.Bcells.common]
C1HTmedium.Bcell.mat <- C1HTmedium.20K[Bcells.common.markers, C1HTmedium.Bcells.common]
Chromium.Bcell.mat <- Chromium.20K[Bcells.common.markers, Chromium.Bcells.common]
ChromiumNuclei.Bcell.mat <- ChromiumNuclei.20K[Bcells.common.markers, ChromiumNuclei.Bcells.common]
ddSEQ.Bcell.mat <- ddSEQ.20K[Bcells.common.markers, ddSEQ.Bcells.common]
Dropseq.Bcell.mat <- Dropseq.20K[Bcells.common.markers, Dropseq.Bcells.common]
ICELL8.Bcell.mat <- ICELL8.20K[Bcells.common.markers, ICELL8.Bcells.common]
inDrop.Bcell.mat <- inDrop.20K[Bcells.common.markers, inDrop.Bcells.common]


library(dplyr)
library(ggplot2)
library(gplots)

Bcell.heatmap.df <- log(as.matrix(bind_cols(CELseq2.Bcell.mat, MARSseq.Bcell.mat, QUARTZseq.Bcell.mat, SCRBseq.Bcell.mat,SMARTseq2.Bcell.mat,
                                            C1HTsmall.Bcell.mat, C1HTmedium.Bcell.mat, Chromium.Bcell.mat, ChromiumNuclei.Bcell.mat,
                                            ddSEQ.Bcell.mat,Dropseq.Bcell.mat,ICELL8.Bcell.mat, inDrop.Bcell.mat)) + 1)
rownames(Bcell.heatmap.df) <- Bcells.common.markers

techniques <- c("CELseq2","MARSseq","QUARTZseq","SCRBseq","SMARTseq2", "C1HTsmall","C1HTmedium","Chromium", "ChromiumNuclei","ddSEQ","Dropseq","ICELL8","inDrop") 
colseps <- c()
tech.ncols <- c()
vec.Bcell.matrices <- vector("list", length = length(techniques))
for (i in 1:length(techniques)){
  tech = techniques[i]
  tech.mat <- get(paste(tech,".Bcell.mat", sep = ""))
  vec.Bcell.matrices[[i]] <- log(as.matrix(tech.mat)+1)
  num.col <- ncol(tech.mat)
  tech.ncols <- c(tech.ncols,num.col)
  if (length(colseps) < 1){colseps <- c(num.col)}
  else{colseps <- c(colseps, colseps[i-1] + num.col)}
}

#selected.color = c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499")
selected.color = c("tan2", "tomato1", "red3", "palevioletred1","purple", "maroon2","darkmagenta", "#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933")

col.sep.color <- rep(selected.color, tech.ncols)

my_palette <- colorRampPalette(c("#FFECB3", "#E85285", "#6A1B9A"))(n = 299)
png("Markers_heatmap_DS20K_50cells_Bcells_normal_V2.png", width = 8, height = 8, units = "in", res = 600)
heatmap.2(Bcell.heatmap.df, Rowv = T, Colv = F, trace = "none", dendrogram = c("none"),
          colsep = colseps, sepcolor = "black",sepwidth = c(0.8,0.8), key= F,
          ColSideColors= col.sep.color, labCol = F, col = my_palette, cexRow = 0.4, main = "B Cells markers downsampled to 20K")
#legend(-0.1,0.8,legend=c("CEL-Seq2", "MARS-Seq", "Quartz-Seq2", "SCRB-Seq","Smart-Seq2", "C1HT-small", "C1HT-medium", "Chromium", "Chromium(sn)", "ddSEQ", "Drop-Seq", "ICELL8","inDrop"),
#       fill=selected.color,border=FALSE, bty="n", y.intersp = 0.8, cex=0.7, xpd = T)

dev.off()

mean.tech.matrices.Bcell <- unlist(lapply(vec.Bcell.matrices, function(x) mean(x)))
names(mean.tech.matrices.Bcell) <- techniques



# Monocytes heatmap
Monocytes.common.markers <- Reduce(intersect, list(ref.markers[["CD14+ Monocytes"]], rownames(CELseq2.20K), rownames(MARSseq.20K), rownames(QUARTZseq.20K), rownames(SCRBseq.20K), 
                                                   rownames(SMARTseq2.20K) , rownames(C1HTsmall.20K), rownames(C1HTmedium.20K), rownames(Chromium.20K), rownames(ChromiumNuclei.20K),
                                                   rownames(ddSEQ.20K) , rownames(Dropseq.20K), rownames(ICELL8.20K), rownames(inDrop.20K)))
print(paste("length common Monocytes markers = ", length(Monocytes.common), sep=""))

CELseq2.Monocytes.mat <- CELseq2.20K[Monocytes.common.markers, CELseq2.Monocytes.common]
MARSseq.Monocytes.mat <- MARSseq.20K[Monocytes.common.markers, MARSseq.Monocytes.common]
QUARTZseq.Monocytes.mat <- QUARTZseq.20K[Monocytes.common.markers, QUARTZseq.Monocytes.common]
SCRBseq.Monocytes.mat <- SCRBseq.20K[Monocytes.common.markers, SCRBseq.Monocytes.common]
SMARTseq2.Monocytes.mat <- SMARTseq2.20K[Monocytes.common.markers, SMARTseq2.Monocytes.common]
C1HTsmall.Monocytes.mat <- C1HTsmall.20K[Monocytes.common.markers, C1HTsmall.Monocytes.common]
C1HTmedium.Monocytes.mat <- C1HTmedium.20K[Monocytes.common.markers, C1HTmedium.Monocytes.common]
Chromium.Monocytes.mat <- Chromium.20K[Monocytes.common.markers, Chromium.Monocytes.common]
ChromiumNuclei.Monocytes.mat <- ChromiumNuclei.20K[Monocytes.common.markers, ChromiumNuclei.Monocytes.common]
ddSEQ.Monocytes.mat <- ddSEQ.20K[Monocytes.common.markers, ddSEQ.Monocytes.common]
Dropseq.Monocytes.mat <- Dropseq.20K[Monocytes.common.markers, Dropseq.Monocytes.common]
ICELL8.Monocytes.mat <- ICELL8.20K[Monocytes.common.markers, ICELL8.Monocytes.common]
inDrop.Monocytes.mat <- inDrop.20K[Monocytes.common.markers, inDrop.Monocytes.common]

Monocytes.heatmap.df <- log(as.matrix(bind_cols(CELseq2.Monocytes.mat, MARSseq.Monocytes.mat, QUARTZseq.Monocytes.mat, SCRBseq.Monocytes.mat,SMARTseq2.Monocytes.mat,
                                                C1HTsmall.Monocytes.mat, C1HTmedium.Monocytes.mat, Chromium.Monocytes.mat, ChromiumNuclei.Monocytes.mat,
                                                ddSEQ.Monocytes.mat,Dropseq.Monocytes.mat,ICELL8.Monocytes.mat, inDrop.Monocytes.mat)) + 1)
rownames(Monocytes.heatmap.df) <- Monocytes.common.markers

techniques <- c("CELseq2","MARSseq","QUARTZseq","SCRBseq","SMARTseq2", "C1HTsmall","C1HTmedium","Chromium", "ChromiumNuclei","ddSEQ","Dropseq","ICELL8","inDrop") 
colseps <- c()
tech.ncols <- c()
vec.Monocytes.matrices <- vector("list", length = length(techniques))
for (i in 1:length(techniques)){
  tech = techniques[i]
  tech.mat <- get(paste(tech,".Monocytes.mat", sep = ""))
  vec.Monocytes.matrices[[i]] <- log(as.matrix(tech.mat)+1)
  num.col <- ncol(tech.mat)
  tech.ncols <- c(tech.ncols,num.col)
  if (length(colseps) < 1){colseps <- c(num.col)}
  else{colseps <- c(colseps, colseps[i-1] + num.col)}
}

#selected.color = c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499")
#selected.color = c("tan2", "tomato1", "red3", "palevioletred1", "maroon2","darkmagenta", "#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933")
load("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/techniques_colors_2.RData")

col.sep.color <- rep(selected.color, tech.ncols)

my_palette <- colorRampPalette(c("#FFECB3", "#E85285", "#6A1B9A"))(n = 299)
png("Markers_heatmap_DS20K_50cells_Monocytes_normal_V2.png", width = 8, height = 8, units = "in", res = 600)
heatmap.2(Monocytes.heatmap.df, Rowv = T, Colv = F, trace = "none", dendrogram = c("none"),
          colsep = colseps, sepcolor = "black",sepwidth = c(0.8,0.8), key = T,cexCol = 0.7,
          ColSideColors= col.sep.color, labCol = F, col = my_palette, cexRow = 0.4, main = "Monocytes markers downsampled to 20K")
#legend(-0.1,0.8,legend=c("CEL-Seq2", "MARS-Seq", "Quartz-Seq2", "SCRB-Seq","Smart-Seq2", "C1HT-small", "C1HT-medium", "Chromium", "Chromium(sn)", "ddSEQ", "Drop-Seq", "ICELL8","inDrop"),
#       fill=selected.color,border=FALSE, bty="n", y.intersp = 0.8, cex=0.7, xpd = T)

dev.off()

mean.tech.matrices.Monocytes <- unlist(lapply(vec.Monocytes.matrices, function(x) mean(x)))
names(mean.tech.matrices.Monocytes) <- techniques



# HEK heatmap
HEK.common.markers <- Reduce(intersect, list(ref.markers[["HEK cells"]], rownames(CELseq2.20K), rownames(MARSseq.20K), rownames(QUARTZseq.20K), rownames(SCRBseq.20K), 
                                             rownames(SMARTseq2.20K) , rownames(C1HTsmall.20K), rownames(C1HTmedium.20K), rownames(Chromium.20K), rownames(ChromiumNuclei.20K),
                                             rownames(ddSEQ.20K) , rownames(Dropseq.20K), rownames(ICELL8.20K), rownames(inDrop.20K)))
#print(paste("length common HEK markers = ", length(HEK.common), sep=""))

CELseq2.HEK.mat <- CELseq2.20K[HEK.common.markers, CELseq2.HEK.common]
MARSseq.HEK.mat <- MARSseq.20K[HEK.common.markers, MARSseq.HEK.common]
QUARTZseq.HEK.mat <- QUARTZseq.20K[HEK.common.markers, QUARTZseq.HEK.common]
SCRBseq.HEK.mat <- SCRBseq.20K[HEK.common.markers, SCRBseq.HEK.common]
SMARTseq2.HEK.mat <- SMARTseq2.20K[HEK.common.markers, SMARTseq2.HEK.common]
C1HTsmall.HEK.mat <- C1HTsmall.20K[HEK.common.markers, C1HTsmall.HEK.common]
C1HTmedium.HEK.mat <- C1HTmedium.20K[HEK.common.markers, C1HTmedium.HEK.common]
Chromium.HEK.mat <- Chromium.20K[HEK.common.markers, Chromium.HEK.common]
ChromiumNuclei.HEK.mat <- ChromiumNuclei.20K[HEK.common.markers, ChromiumNuclei.HEK.common]
ddSEQ.HEK.mat <- ddSEQ.20K[HEK.common.markers, ddSEQ.HEK.common]
Dropseq.HEK.mat <- Dropseq.20K[HEK.common.markers, Dropseq.HEK.common]
ICELL8.HEK.mat <- ICELL8.20K[HEK.common.markers, ICELL8.HEK.common]
inDrop.HEK.mat <- inDrop.20K[HEK.common.markers, inDrop.HEK.common]


HEK.heatmap.df <- log(as.matrix(bind_cols(CELseq2.HEK.mat, MARSseq.HEK.mat, QUARTZseq.HEK.mat, SCRBseq.HEK.mat,SMARTseq2.HEK.mat,
                                          C1HTsmall.HEK.mat, C1HTmedium.HEK.mat, Chromium.HEK.mat, ChromiumNuclei.HEK.mat,
                                          ddSEQ.HEK.mat,Dropseq.HEK.mat,ICELL8.HEK.mat, inDrop.HEK.mat)) + 1)
rownames(HEK.heatmap.df) <- HEK.common.markers

techniques <- c("CELseq2","MARSseq","QUARTZseq","SCRBseq","SMARTseq2", "C1HTsmall","C1HTmedium","Chromium", "ChromiumNuclei","ddSEQ","Dropseq","ICELL8","inDrop") 
colseps <- c()
tech.ncols <- c()
vec.HEK.matrices <- vector("list", length = length(techniques))
for (i in 1:length(techniques)){
  tech = techniques[i]
  tech.mat <- get(paste(tech,".HEK.mat", sep = ""))
  vec.HEK.matrices[[i]] <- log(as.matrix(tech.mat)+1)
  num.col <- ncol(tech.mat)
  tech.ncols <- c(tech.ncols,num.col)
  if (length(colseps) < 1){colseps <- c(num.col)}
  else{colseps <- c(colseps, colseps[i-1] + num.col)}
}

#selected.color = c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499")
#selected.color = c("tan2", "tomato1", "red3", "palevioletred1", "maroon2","darkmagenta", "#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933")

col.sep.color <- rep(selected.color, tech.ncols)

my_palette <- colorRampPalette(c("#FFECB3", "#E85285", "#6A1B9A"))(n = 299)
png("Markers_heatmap_DS20K_50cells_HEK_normal_V2_colorkey.png", width = 8, height = 8, units = "in", res = 600)
heatmap.2(HEK.heatmap.df, Rowv = T, Colv = F, trace = "none", dendrogram = c("none"),
          colsep = colseps, sepcolor = "black",sepwidth = c(0.8,0.8), key = T,density.info="none",
          ColSideColors= col.sep.color, labCol = F, col = my_palette, cexRow = 0.4, main = "HEK cells markers downsampled to 20K")
#legend(-0.1,0.8,legend=c("CEL-Seq2", "MARS-Seq", "Quartz-Seq2", "SCRB-Seq","Smart-Seq2", "C1HT-small", "C1HT-medium", "Chromium", "Chromium(sn)", "ddSEQ", "Drop-Seq", "ICELL8","inDrop"),
#       fill=selected.color,border=FALSE, bty="n", y.intersp = 0.8, cex=0.7, xpd = T)

dev.off()

mean.tech.matrices.HEK <- unlist(lapply(vec.HEK.matrices, function(x) mean(x)))
names(mean.tech.matrices.HEK) <- techniques

library(devEMF)


#The overall mean hetamap per tech per celltype
whole.mean.df <- cbind(Monocytes=mean.tech.matrices.Monocytes, Bcells=mean.tech.matrices.Bcell, HEK=mean.tech.matrices.HEK)
rownames(whole.mean.df) <- c("CEL-Seq2", "MARS-Seq", "Quartz-Seq2", "gmcSCRB-seq","Smart-seq2", "C1HT-Small", "C1HT-Medium", "Chromium", "Chromium(sn)", "ddSEQ", "Drop-seq", "ICELL8","inDrop")
#col.sep.color = c("slateblue1", "orchid1")#, "red", "olivedrab3")
my_palette <- colorRampPalette(c("#FFECB3", "#E85285", "#6A1B9A"))(n = 20)
png("Markers_comparison_MeanPerTechPerCelltype_50cells_V3.png", width = 10, height = 12, units = "in", res=600)
#pdf("Markers_comparison_MeanPerTechPerCelltype_50cells_V2_colorkey.pdf")
heatmap.2(whole.mean.df, Rowv = T, Colv = F, trace = "none", dendrogram = c("none"), col = my_palette, colsep = c(1,2),sepcolor = "black",sepwidth = c(0.02,0.02), cexRow = 3.5, cexCol = 3, srtCol= 50, offsetCol = 0.5,adjCol= 1, margins = c(10,21), key = F,density.info="none") #ColSideColors= col.sep.color,
dev.off()

save(whole.mean.df, file= "HEATMAP_Markers_comparison_MeanPerTechPerCelltype_50cells.RData")
load("HEATMAP_Markers_comparison_MeanPerTechPerCelltype_50cells.RData")


