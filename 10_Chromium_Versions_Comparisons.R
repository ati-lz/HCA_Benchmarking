### Here we saw that for Chrom V3 non-viability (Lib-92) after DS to 20K we dont have Monocytes and Bcells so we decided to use 10K DS

##### Number of genes and Number of UMIs boxplots ########################################

library(dplyr)
library(Seurat)
#Chromium V2 small scale:
load("chromium_seu.obj.RData")
Chromium.hsap.obj <- chromium
rm(chromium)


head(Chromium.hsap.obj@meta.data)
unique(Chromium.hsap.obj@meta.data$Library)

ChromV2.64221.metadata <- Chromium.hsap.obj@meta.data[,c("nGene","nUMI","orig.ident","nTReads","nExonReads", "Library", "clean.id")]
ChromV2.all.metadata <- ChromV2.64221.metadata
ChromV2.all.expressionMat <- Chromium.hsap.obj@raw.data

load("10X2x5K.hsap.full.SCE.jointDSmat.Robj")
ChromV2.64221.DS10K <- output.readcount.umicount.joint.mats$UMI$downsampled_10000
rm(output.readcount.umicount.joint.mats)
ChromV2.64221.common.cells <- intersect(colnames(ChromV2.64221.DS10K), rownames(ChromV2.64221.metadata))
ChromV2.64221.DS10K <- ChromV2.64221.DS10K[, ChromV2.64221.common.cells]

#exon reads downsampleds
ChromV2.64221.dgecounts <- readRDS("10X2x5K.64221.dgecounts.rds")
ChromV2.64221.DS10K.exon <- ChromV2.64221.dgecounts$umicount$exon$downsampling$downsampled_10000
rm(ChromV2.64221.dgecounts)
colnames(ChromV2.64221.DS10K.exon) <- paste("10X2x5K_64221", colnames(ChromV2.64221.DS10K.exon), sep = "_")
ChromV2.64221.DS10K.exon <- ChromV2.64221.DS10K.exon[,colnames(ChromV2.64221.DS10K)]


library(rowr)
ChromV2.all.DS10K <-  ChromV2.64221.DS10K
#rownames(ChromV2.all.DS10K) <- rownames(ChromV2.64220.DS10K)
#colnames(ChromV2.all.DS10K) <- c(colnames(ChromV2.64220.DS10K), colnames(ChromV2.64221.DS10K))
ChromV2.all.metadata
common.cells.all <- intersect(colnames(ChromV2.all.DS10K), rownames(ChromV2.all.metadata))
ChromV2.all.metadata.inex <- ChromV2.all.metadata[common.cells.all,]
ChromV2.all.DS10K.HEKs <- ChromV2.all.DS10K[,rownames(ChromV2.all.metadata.inex[which(ChromV2.all.metadata.inex$clean.id == "HEK"),])]
ChromV2.DS10K.HEK.obj <- CreateSeuratObject(as.matrix(ChromV2.all.DS10K.HEKs), project = "Chromium")
head(ChromV2.DS10K.HEK.obj@meta.data)
ChromV2.all.DS10K.Monocytes <- ChromV2.all.DS10K[,rownames(ChromV2.all.metadata.inex[which(ChromV2.all.metadata.inex$clean.id == "Monocytes"),])]
ChromV2.DS10K.Monocytes.obj <- CreateSeuratObject(as.matrix(ChromV2.all.DS10K.Monocytes), project = "Chromium")
head(ChromV2.DS10K.Monocytes.obj@meta.data)
ChromV2.all.DS10K.B <- ChromV2.all.DS10K[,rownames(ChromV2.all.metadata.inex[which(ChromV2.all.metadata.inex$clean.id == "B"),])]
ChromV2.DS10K.B.obj <- CreateSeuratObject(as.matrix(ChromV2.all.DS10K.B), project = "Chromium")
head(ChromV2.DS10K.B.obj@meta.data)


ChromV2.all.DS10K.exon <- ChromV2.64221.DS10K.exon
#rownames(ChromV2.all.DS10K.exon) <- rownames(ChromV2.64220.DS10K.exon)
#colnames(ChromV2.all.DS10K.exon) <- c(colnames(ChromV2.64220.DS10K.exon), colnames(ChromV2.64221.DS10K.exon))
ChromV2.all.metadata
common.cells.all.exon <- intersect(colnames(ChromV2.all.DS10K.exon), rownames(ChromV2.all.metadata))
ChromV2.all.metadata.exon <- ChromV2.all.metadata[common.cells.all.exon,]
ChromV2.all.DS10K.HEKs.exon <- ChromV2.all.DS10K.exon[,rownames(ChromV2.all.metadata.exon[which(ChromV2.all.metadata.exon$clean.id == "HEK"),])]
ChromV2.DS10K.HEK.exon.obj <- CreateSeuratObject(as.matrix(ChromV2.all.DS10K.HEKs.exon), project = "Chromium")
head(ChromV2.DS10K.HEK.exon.obj@meta.data)
colnames(ChromV2.DS10K.HEK.exon.obj@meta.data)[c(2,3)] <- c("EXON_nCount_RNA", "EXON_nFeature_RNA")
ChromV2.final.table.HEK <- cbind(ChromV2.DS10K.HEK.obj@meta.data, ChromV2.DS10K.HEK.exon.obj@meta.data[,2:3])

ChromV2.all.DS10K.Monocytes.exon <- ChromV2.all.DS10K.exon[,rownames(ChromV2.all.metadata.exon[which(ChromV2.all.metadata.exon$clean.id == "Monocytes"),])]
ChromV2.DS10K.Monocytes.exon.obj <- CreateSeuratObject(as.matrix(ChromV2.all.DS10K.Monocytes.exon), project = "Chromium")
head(ChromV2.DS10K.Monocytes.exon.obj@meta.data)
colnames(ChromV2.DS10K.Monocytes.exon.obj@meta.data)[c(2,3)] <- c("EXON_nCount_RNA", "EXON_nFeature_RNA")
ChromV2.final.table.Monocytes <- cbind(ChromV2.DS10K.Monocytes.obj@meta.data, ChromV2.DS10K.Monocytes.exon.obj@meta.data[,2:3])

ChromV2.all.DS10K.B.exon <- ChromV2.all.DS10K.exon[,rownames(ChromV2.all.metadata.exon[which(ChromV2.all.metadata.exon$clean.id == "B"),])]
ChromV2.DS10K.B.exon.obj <- CreateSeuratObject(as.matrix(ChromV2.all.DS10K.B.exon), project = "Chromium")
head(ChromV2.DS10K.B.exon.obj@meta.data)
colnames(ChromV2.DS10K.B.exon.obj@meta.data)[c(2,3)] <- c("EXON_nCount_RNA", "EXON_nFeature_RNA")
ChromV2.final.table.B <- cbind(ChromV2.DS10K.B.obj@meta.data, ChromV2.DS10K.B.exon.obj@meta.data[,2:3])


#Chrom Nuclei
load("nuclei_seu.obj.RData")
ChromNuclei.hsap.obj <- nuclei
rm(nuclei)
ChromNuclei.metadata <- ChromNuclei.hsap.obj@meta.data[,c("nGene","nUMI","orig.ident","nTReads","nExonReads", "Library", "clean.id")]

load("Nuclei10X.hsap.full.SCE.jointDSmat.Robj")
ChromNuclei.DS10K <- output.readcount.umicount.joint.mats$UMI$downsampled_10000
rm(output.readcount.umicount.joint.mats)
ChromNuclei.common.cells <- intersect(colnames(ChromNuclei.DS10K), rownames(ChromNuclei.metadata))
ChromNuclei.DS10K <- ChromNuclei.DS10K[, ChromNuclei.common.cells]
ChromNuclei.metadata.inex <- ChromNuclei.metadata[ChromNuclei.common.cells,]
ChromNuclei.DS10K.HEKs <- ChromNuclei.DS10K[,rownames(ChromNuclei.metadata.inex[which(ChromNuclei.metadata.inex$clean.id == "HEK"),])]
ChromNuclei.DS10K.HEK.obj <- CreateSeuratObject(as.matrix(ChromNuclei.DS10K.HEKs), project = "Chromium")
head(ChromNuclei.DS10K.HEK.obj@meta.data)
ChromNuclei.DS10K.Monocytes <- ChromNuclei.DS10K[,rownames(ChromNuclei.metadata.inex[which(ChromNuclei.metadata.inex$clean.id == "Monocytes"),])]
ChromNuclei.DS10K.Monocytes.obj <- CreateSeuratObject(as.matrix(ChromNuclei.DS10K.Monocytes), project = "Chromium")
head(ChromNuclei.DS10K.Monocytes.obj@meta.data)
ChromNuclei.DS10K.B <- ChromNuclei.DS10K[,rownames(ChromNuclei.metadata.inex[which(ChromNuclei.metadata.inex$clean.id == "B"),])]
ChromNuclei.DS10K.B.obj <- CreateSeuratObject(as.matrix(ChromNuclei.DS10K.B), project = "Chromium")
head(ChromNuclei.DS10K.B.obj@meta.data)

load("Nuclei10X.hsap.full.SCE.jointDSmat.exonReads.Robj")
ChromNuclei.DS10K.exon <- final.sample.merged.mat$downsampled_10000
rm(final.sample.merged.mat)
ChromNuclei.common.cells.exon <- intersect(colnames(ChromNuclei.DS10K.exon), rownames(ChromNuclei.metadata))
ChromNuclei.DS10K.exon <- ChromNuclei.DS10K.exon[, ChromNuclei.common.cells.exon]
ChromNuclei.metadata.exon <- ChromNuclei.metadata[ChromNuclei.common.cells.exon,]
ChromNuclei.DS10K.exon.HEKs <- ChromNuclei.DS10K.exon[,rownames(ChromNuclei.metadata.exon[which(ChromNuclei.metadata.exon$clean.id == "HEK"),])]
ChromNuclei.DS10K.exon.HEK.obj <- CreateSeuratObject(as.matrix(ChromNuclei.DS10K.exon.HEKs), project = "Chromium")
head(ChromNuclei.DS10K.exon.HEK.obj@meta.data)
colnames(ChromNuclei.DS10K.exon.HEK.obj@meta.data)[c(2,3)] <- c("EXON_nCount_RNA", "EXON_nFeature_RNA")
ChromNuclei.final.table.HEK <- cbind(ChromNuclei.DS10K.HEK.obj@meta.data, ChromNuclei.DS10K.exon.HEK.obj@meta.data[,2:3])

ChromNuclei.DS10K.exon.Monocytes <- ChromNuclei.DS10K.exon[,rownames(ChromNuclei.metadata.exon[which(ChromNuclei.metadata.exon$clean.id == "Monocytes"),])]
ChromNuclei.DS10K.exon.Monocytes.obj <- CreateSeuratObject(as.matrix(ChromNuclei.DS10K.exon.Monocytes), project = "Chromium")
head(ChromNuclei.DS10K.exon.Monocytes.obj@meta.data)
colnames(ChromNuclei.DS10K.exon.Monocytes.obj@meta.data)[c(2,3)] <- c("EXON_nCount_RNA", "EXON_nFeature_RNA")
ChromNuclei.final.table.Monocytes <- cbind(ChromNuclei.DS10K.Monocytes.obj@meta.data, ChromNuclei.DS10K.exon.Monocytes.obj@meta.data[,2:3])

ChromNuclei.DS10K.exon.B <- ChromNuclei.DS10K.exon[,rownames(ChromNuclei.metadata.exon[which(ChromNuclei.metadata.exon$clean.id == "B"),])]
ChromNuclei.DS10K.exon.B.obj <- CreateSeuratObject(as.matrix(ChromNuclei.DS10K.exon.B), project = "Chromium")
head(ChromNuclei.DS10K.exon.B.obj@meta.data)
colnames(ChromNuclei.DS10K.exon.B.obj@meta.data)[c(2,3)] <- c("EXON_nCount_RNA", "EXON_nFeature_RNA")
ChromNuclei.final.table.B <- cbind(ChromNuclei.DS10K.B.obj@meta.data, ChromNuclei.DS10K.exon.B.obj@meta.data[,2:3])


#Chrom V3
ChromV3.final.table <- data.frame()
load("chromiumV3.Without_Viability_seu.obj.RData")
ChromV3.hsap.obj <- chromium
rm(chromium)
#ChromV3.metadata <- ChromV3.hsap.obj@meta.data[,c("nGene","nUMI","orig.ident","nTReads","nExonReads", "Library", "clean.id")]
ChromV3.metadata <- ChromV3.hsap.obj@meta.data[,c("nFeature_RNA","nCount_RNA","orig.ident","clean.id")]#,"nTReads","nExonReads", "Library", "clean.id")]
colnames(ChromV3.metadata) <- c("nGene","nUMI","orig.ident", "clean.id")#,"nTReads","nExonReads", "Library", "clean.id")

#load("10XV3Rev.hsap.full.SCE.jointDSmat_5Kcells.Robj")
#ChromV3.DS10K <- output.readcount.umicount.joint.mats$UMI$downsampled_50000
#rm(output.readcount.umicount.joint.mats)
Lib92.ds.data <- readRDS("10XV3Rev_lib92_only/10XV3Rev.dgecounts.rds")
ChromV3.DS10K <- as.matrix(Lib92.ds.data$umicount$inex$downsampling$downsampled_10000)
colnames(ChromV3.DS10K) <- paste("X10XV3Rev_AN4492_", colnames(ChromV3.DS10K), sep = "")
ChromV3.common.cells <- intersect(colnames(ChromV3.DS10K), rownames(ChromV3.metadata))
ChromV3.DS10K <- ChromV3.DS10K[, ChromV3.common.cells]
ChromV3.metadata.inex <- ChromV3.metadata[ChromV3.common.cells,]
ChromV3.DS10K.HEKs <- ChromV3.DS10K[,rownames(ChromV3.metadata.inex[which(ChromV3.metadata.inex$clean.id == "HEK"),])]
ChromV3.DS10K.HEK.obj <- CreateSeuratObject(as.matrix(ChromV3.DS10K.HEKs), project = "Chromium")
head(ChromV3.DS10K.HEK.obj@meta.data)
ChromV3.DS10K.Monocytes <- ChromV3.DS10K[,rownames(ChromV3.metadata.inex[which(ChromV3.metadata.inex$clean.id == "Monocytes"),])]
ChromV3.DS10K.Monocytes.obj <- CreateSeuratObject(as.matrix(ChromV3.DS10K.Monocytes), project = "Chromium")
head(ChromV3.DS10K.Monocytes.obj@meta.data)
ChromV3.DS10K.B <- ChromV3.DS10K[,rownames(ChromV3.metadata.inex[which(ChromV3.metadata.inex$clean.id == "B"),])]
ChromV3.DS10K.B.obj <- CreateSeuratObject(as.matrix(ChromV3.DS10K.B), project = "Chromium")
head(ChromV3.DS10K.B.obj@meta.data)

#load("10XV3Rev.hsap.full.SCE.jointDSmat.exonReads_10Kcells.Robj")
#ChromV3.DS10K.exon <- final.sample.merged.mat$downsampled_10000
#rm(final.sample.merged.mat)
ChromV3.DS10K.exon <- as.matrix(Lib92.ds.data$umicount$exon$downsampling$downsampled_10000)
colnames(ChromV3.DS10K.exon) <- paste("X10XV3Rev_AN4492_", colnames(ChromV3.DS10K.exon), sep = "")
ChromV3.common.cells.exon <- intersect(colnames(ChromV3.DS10K.exon), rownames(ChromV3.metadata))
ChromV3.DS10K.exon <- ChromV3.DS10K.exon[, ChromV3.common.cells.exon]
ChromV3.metadata.exon <- ChromV3.metadata[ChromV3.common.cells.exon,]
ChromV3.DS10K.exon.HEKs <- ChromV3.DS10K.exon[,rownames(ChromV3.metadata.exon[which(ChromV3.metadata.exon$clean.id == "HEK"),])]
ChromV3.DS10K.exon.HEK.obj <- CreateSeuratObject(as.matrix(ChromV3.DS10K.exon.HEKs), project = "Chromium")
head(ChromV3.DS10K.exon.HEK.obj@meta.data)
colnames(ChromV3.DS10K.exon.HEK.obj@meta.data)[c(2,3)] <- c("EXON_nCount_RNA", "EXON_nFeature_RNA")
ChromV3.final.table.HEK <- cbind(ChromV3.DS10K.HEK.obj@meta.data, ChromV3.DS10K.exon.HEK.obj@meta.data[,2:3])

ChromV3.DS10K.exon.Monocytes <- ChromV3.DS10K.exon[,rownames(ChromV3.metadata.exon[which(ChromV3.metadata.exon$clean.id == "Monocytes"),])]
ChromV3.DS10K.exon.Monocytes.obj <- CreateSeuratObject(as.matrix(ChromV3.DS10K.exon.Monocytes), project = "Chromium")
head(ChromV3.DS10K.exon.Monocytes.obj@meta.data)
colnames(ChromV3.DS10K.exon.Monocytes.obj@meta.data)[c(2,3)] <- c("EXON_nCount_RNA", "EXON_nFeature_RNA")
ChromV3.final.table.Monocytes <- cbind(ChromV3.DS10K.Monocytes.obj@meta.data, ChromV3.DS10K.exon.Monocytes.obj@meta.data[,2:3])

ChromV3.DS10K.exon.B <- ChromV3.DS10K.exon[,rownames(ChromV3.metadata.exon[which(ChromV3.metadata.exon$clean.id == "B"),])]
ChromV3.DS10K.exon.B.obj <- CreateSeuratObject(as.matrix(ChromV3.DS10K.exon.B), project = "Chromium")
head(ChromV3.DS10K.exon.B.obj@meta.data)
colnames(ChromV3.DS10K.exon.B.obj@meta.data)[c(2,3)] <- c("EXON_nCount_RNA", "EXON_nFeature_RNA")
ChromV3.final.table.B <- cbind(ChromV3.DS10K.B.obj@meta.data, ChromV3.DS10K.exon.B.obj@meta.data[,2:3])

#Plotting
library(data.table)
library(plyr)
Chroms.all.df.HEK <- rbind(ChromV2.final.table.HEK, ChromNuclei.final.table.HEK, ChromV3.final.table.HEK)
Chroms.all.df.nGenes.HEK <- Chroms.all.df.HEK[, c("orig.ident", "nFeature_RNA", "EXON_nFeature_RNA")]
Chroms.all.df.nGenes.HEK.melted <- melt(Chroms.all.df.nGenes.HEK)
Chroms.all.df.nGenes.HEK.melted$orig.ident <- mapvalues(Chroms.all.df.nGenes.HEK.melted$orig.ident,
                                                        from = c("10X2x5K", "Nuclei10X", "X10XV3Rev"),
                                                        to= c("Chromium (V2)", "Chromium (sn)", "Chromium (V3)"))
Chroms.all.df.nGenes.HEK.melted$orig.ident <- factor(Chroms.all.df.nGenes.HEK.melted$orig.ident, levels = c("Chromium (sn)", "Chromium (V2)","Chromium (V3)"))
#UMI
Chroms.all.df.nUMIs.HEK <- Chroms.all.df.HEK[, c("orig.ident", "nCount_RNA", "EXON_nCount_RNA")]
Chroms.all.df.nUMIs.HEK.melted <- melt(Chroms.all.df.nUMIs.HEK)
Chroms.all.df.nUMIs.HEK.melted$orig.ident <- mapvalues(Chroms.all.df.nUMIs.HEK.melted$orig.ident,
                                                       from = c("10X2x5K", "Nuclei10X", "X10XV3Rev"),
                                                       to= c("Chromium (V2)", "Chromium (sn)", "Chromium (V3)"))
Chroms.all.df.nUMIs.HEK.melted$orig.ident <- factor(Chroms.all.df.nUMIs.HEK.melted$orig.ident, levels = c("Chromium (sn)", "Chromium (V2)","Chromium (V3)"))

#ggplot(Chroms.all.df.nGenes.HEK.melted, aes(x= orig.ident, y=value)) + geom_boxplot(aes(fill = variable),position = position_dodge(0.7)) +
#  scale_fill_brewer(palette="OrRd")
library(patternplot)
library(ggplot2)
pattern.type<-c('hdashes', 'crosshatch')
pattern.type=c('blank', 'nwlines')
pattern.color=c('grey20', 'grey20')
background.color=c('lightblue','lightblue')
density<-c(5,5)
png("ChromCompare_nGenes_detected_HEK_ChromV2-21_ChromV3-Lib92_DS10K.png", width = 12, height = 12, units = "in", res = 600)
patternboxplot(Chroms.all.df.nGenes.HEK.melted,Chroms.all.df.nGenes.HEK.melted$orig.ident, Chroms.all.df.nGenes.HEK.melted$value, group=Chroms.all.df.nGenes.HEK.melted$variable,pattern.type=pattern.type, 
               pattern.color=pattern.color, background.color=background.color, 
               density=density, pixel=1.2, xlab='', ylab='nGenes', 
               legend.h=5, legend.y.pos=0.5, legend.ratio1=0.1,
               legend.x.pos=0.1) +
  theme (axis.text.x = element_text(angle = 0, hjust = 0.5, size = 30, face = "bold", colour = "black"),axis.text.y = element_text(angle = 45,size = 40, face = "bold", colour = "black"), legend.position = "none",
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
         axis.line = element_line(colour = "black"))
dev.off()
png("ChromCompare_nUMIs_detected_HEK_ChromV2-21_ChromV3-Lib92_DS10K.png", width = 12, height = 12, units = "in", res = 600)
patternboxplot(Chroms.all.df.nUMIs.HEK.melted,Chroms.all.df.nUMIs.HEK.melted$orig.ident, Chroms.all.df.nUMIs.HEK.melted$value, group=Chroms.all.df.nUMIs.HEK.melted$variable,pattern.type=pattern.type, 
               pattern.color=pattern.color, background.color=background.color, 
               density=density, pixel=1.2, xlab='', ylab='nUMIs', 
               legend.h=5, legend.y.pos=0.5, legend.ratio1=0.1,
               legend.x.pos=0.1) +
  theme (axis.text.x = element_text(angle = 0, hjust = 0.5, size = 30, face = "bold", colour = "black"),axis.text.y = element_text(angle = 45,size = 40, face = "bold", colour = "black"), legend.position = "none",
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
         axis.line = element_line(colour = "black"))
dev.off()


# Monocytes
Chroms.all.df.Monocytes <- rbind(ChromV2.final.table.Monocytes, ChromNuclei.final.table.Monocytes, ChromV3.final.table.Monocytes)
Chroms.all.df.nGenes.Monocytes <- Chroms.all.df.Monocytes[, c("orig.ident", "nFeature_RNA", "EXON_nFeature_RNA")]
Chroms.all.df.nGenes.Monocytes.melted <- melt(Chroms.all.df.nGenes.Monocytes)
Chroms.all.df.nGenes.Monocytes.melted$orig.ident <- mapvalues(Chroms.all.df.nGenes.Monocytes.melted$orig.ident,
                                                              from = c("10X2x5K", "Nuclei10X", "X10XV3Rev"),
                                                              to= c("Chromium (V2)", "Chromium (sn)", "Chromium (V3)"))
Chroms.all.df.nGenes.Monocytes.melted$orig.ident <- factor(Chroms.all.df.nGenes.Monocytes.melted$orig.ident, levels = c("Chromium (sn)", "Chromium (V2)","Chromium (V3)"))
#UMI
Chroms.all.df.nUMIs.Monocytes <- Chroms.all.df.Monocytes[, c("orig.ident", "nCount_RNA", "EXON_nCount_RNA")]
Chroms.all.df.nUMIs.Monocytes.melted <- melt(Chroms.all.df.nUMIs.Monocytes)
Chroms.all.df.nUMIs.Monocytes.melted$orig.ident <- mapvalues(Chroms.all.df.nUMIs.Monocytes.melted$orig.ident,
                                                             from = c("10X2x5K", "Nuclei10X", "X10XV3Rev"),
                                                             to= c("Chromium (V2)", "Chromium (sn)", "Chromium (V3)"))
Chroms.all.df.nUMIs.Monocytes.melted$orig.ident <- factor(Chroms.all.df.nUMIs.Monocytes.melted$orig.ident, levels = c("Chromium (sn)", "Chromium (V2)","Chromium (V3)"))

library(patternplot)
pattern.type<-c('hdashes', 'crosshatch')
pattern.type=c('blank', 'nwlines')
pattern.color=c('grey20', 'grey20')
background.color=c('sandybrown','sandybrown')
density<-c(5,5)
png("ChromCompare_nGenes_detected_Monocytes_ChromV2-21_ChromV3-Lib92_DS10K.png", width = 12, height = 12, units = "in", res = 600)
patternboxplot(Chroms.all.df.nGenes.Monocytes.melted,Chroms.all.df.nGenes.Monocytes.melted$orig.ident, Chroms.all.df.nGenes.Monocytes.melted$value, group=Chroms.all.df.nGenes.Monocytes.melted$variable,pattern.type=pattern.type, 
               pattern.color=pattern.color, background.color=background.color, 
               density=density, pixel=1.2, xlab='', ylab='nGenes', 
               legend.h=5, legend.y.pos=0.5, legend.ratio1=0.1,
               legend.x.pos=0.1) +
  theme (axis.text.x = element_text(angle = 0, hjust = 0.5, size = 30, face = "bold", colour = "black"),axis.text.y = element_text(angle = 45,size = 40, face = "bold", colour = "black"), legend.position = "none",
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
         axis.line = element_line(colour = "black"))
dev.off()
png("ChromCompare_nUMIs_detected_Monocytes_ChromV2-21_ChromV3-Lib92_DS10K.png", width = 12, height = 12, units = "in", res = 600)
patternboxplot(Chroms.all.df.nUMIs.Monocytes.melted,Chroms.all.df.nUMIs.Monocytes.melted$orig.ident, Chroms.all.df.nUMIs.Monocytes.melted$value, group=Chroms.all.df.nUMIs.Monocytes.melted$variable,pattern.type=pattern.type, 
               pattern.color=pattern.color, background.color=background.color, 
               density=density, pixel=1.2, xlab='', ylab='nUMIs', 
               legend.h=5, legend.y.pos=0.5, legend.ratio1=0.1,
               legend.x.pos=0.1) +
  theme (axis.text.x = element_text(angle = 0, hjust = 0.5, size = 30, face = "bold", colour = "black"),axis.text.y = element_text(angle = 45,size = 40, face = "bold", colour = "black"), legend.position = "none",
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
         axis.line = element_line(colour = "black"))
dev.off()


# B cells
Chroms.all.df.B <- rbind(ChromV2.final.table.B, ChromNuclei.final.table.B, ChromV3.final.table.B)
Chroms.all.df.nGenes.B <- Chroms.all.df.B[, c("orig.ident", "nFeature_RNA", "EXON_nFeature_RNA")]
Chroms.all.df.nGenes.B.melted <- melt(Chroms.all.df.nGenes.B)
Chroms.all.df.nGenes.B.melted$orig.ident <- mapvalues(Chroms.all.df.nGenes.B.melted$orig.ident,
                                                      from = c("10X2x5K", "Nuclei10X", "X10XV3Rev"),
                                                      to= c("Chromium (V2)", "Chromium (sn)", "Chromium (V3)"))
Chroms.all.df.nGenes.B.melted$orig.ident <- factor(Chroms.all.df.nGenes.B.melted$orig.ident, levels = c("Chromium (sn)", "Chromium (V2)","Chromium (V3)"))
#UMI
Chroms.all.df.nUMIs.B <- Chroms.all.df.B[, c("orig.ident", "nCount_RNA", "EXON_nCount_RNA")]
Chroms.all.df.nUMIs.B.melted <- melt(Chroms.all.df.nUMIs.B)
Chroms.all.df.nUMIs.B.melted$orig.ident <- mapvalues(Chroms.all.df.nUMIs.B.melted$orig.ident,
                                                     from = c("10X2x5K", "Nuclei10X", "X10XV3Rev"),
                                                     to= c("Chromium (V2)", "Chromium (sn)", "Chromium (V3)"))
Chroms.all.df.nUMIs.B.melted$orig.ident <- factor(Chroms.all.df.nUMIs.B.melted$orig.ident, levels = c("Chromium (sn)", "Chromium (V2)","Chromium (V3)"))

library(patternplot)
pattern.type<-c('hdashes', 'crosshatch')
pattern.type=c('blank', 'nwlines')
pattern.color=c('grey20', 'grey20')
background.color=c('palegreen3','palegreen3')
density<-c(5,5)
png("ChromCompare_nGenes_detected_Bcells_ChromV2-21_ChromV3-Lib92_DS10K.png", width = 12, height = 12, units = "in", res = 600)
patternboxplot(Chroms.all.df.nGenes.B.melted,Chroms.all.df.nGenes.B.melted$orig.ident, Chroms.all.df.nGenes.B.melted$value, group=Chroms.all.df.nGenes.B.melted$variable,pattern.type=pattern.type, 
               pattern.color=pattern.color, background.color=background.color, 
               density=density, pixel=1.2, xlab='', ylab='nGenes', 
               legend.h=5, legend.y.pos=0.5, legend.ratio1=0.1,
               legend.x.pos=0.1) +
  theme (axis.text.x = element_text(angle = 0, hjust = 0.5, size = 30, face = "bold", colour = "black"),axis.text.y = element_text(angle = 45,size = 40, face = "bold", colour = "black"), legend.position = "none",
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
         axis.line = element_line(colour = "black"))
dev.off()
png("ChromCompare_nUMIs_detected_Bcells_ChromV2-21_ChromV3-Lib92_DS10K.png", width = 12, height = 12, units = "in", res = 600)
patternboxplot(Chroms.all.df.nUMIs.B.melted,Chroms.all.df.nUMIs.B.melted$orig.ident, Chroms.all.df.nUMIs.B.melted$value, group=Chroms.all.df.nUMIs.B.melted$variable,pattern.type=pattern.type, 
               pattern.color=pattern.color, background.color=background.color, 
               density=density, pixel=1.2, xlab='', ylab='nUMIs', 
               legend.h=5, legend.y.pos=0.5, legend.ratio1=0.1,
               legend.x.pos=0.1) +
  theme (axis.text.x = element_text(angle = 0, hjust = 0.5, size = 30, face = "bold", colour = "black"),axis.text.y = element_text(angle = 45,size = 40, face = "bold", colour = "black"), legend.position = "none",
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
         axis.line = element_line(colour = "black"))
dev.off()

###End###




##### CUULATIVE COMPARISONS #############

#HEKS
set.seed(123)
dim(ChromV2.all.DS10K.HEKs)
ChromV2.sampled.HEKs <- sample(colnames(ChromV2.all.DS10K.HEKs),50)
ChromV2.all.DS10K.HEKs.sampled <- ChromV2.all.DS10K.HEKs[,ChromV2.sampled.HEKs]
ChromV2.HEKs.expressed.genes <- which(rowSums(ChromV2.all.DS10K.HEKs.sampled) > 0)
ChromV2.all.DS10K.HEKs.sampled <- ChromV2.all.DS10K.HEKs.sampled[ChromV2.HEKs.expressed.genes,]

HEK.plot.df <- data.frame()
tech.cumul.gene.numbers <- rep(NA, ncol(ChromV2.all.DS10K.HEKs.sampled))
for (cell in 1:ncol(ChromV2.all.DS10K.HEKs.sampled)){
  print(cell)
  sample.size = cell
  sample.gene.numbers <- rep(NA,50)
  for (i in 1:50){
    selected.cells <- sample(colnames(ChromV2.all.DS10K.HEKs.sampled), sample.size)
    sample.gene.names <- rep(NA, 70000)
    for (cell2 in selected.cells){
      vec.start.point <- (length(sample.gene.names[which(is.na(sample.gene.names) == FALSE)])) +1
      detected.genes <- rownames(ChromV2.all.DS10K.HEKs.sampled[which(ChromV2.all.DS10K.HEKs.sampled[,cell2] > 0),])
      new.genes <- setdiff(detected.genes,sample.gene.names)
      vec.end.point <- vec.start.point + length(new.genes) -1
      sample.gene.names[vec.start.point:vec.end.point] <- new.genes
      #detected.genes <- rownames(DS.mat.HEKS[which(DS.mat.HEKS[,cell2] > 0),])
      #sample.gene.names <- c(sample.gene.names,detected.genes)
    }
    sample.gene.names <- na.omit(sample.gene.names)
    sample.gene.names.uniq <- unique(sample.gene.names)
    num.uniq.genes <- length(sample.gene.names.uniq)
    if (cell == 46){ ChromV2.all.cumul.genes.46 <-  sample.gene.names.uniq}
    sample.gene.numbers[i] <- num.uniq.genes
  }
  sample.gene.numbers <- na.omit(sample.gene.numbers)
  print(paste(cell," = ", mean(sample.gene.numbers), sep = ""))
  tech.cumul.gene.numbers[cell] <- mean(sample.gene.numbers)
}
HEK.plot.df <- rbind(HEK.plot.df, data.frame(Cumul= tech.cumul.gene.numbers, tech= rep("ChromV2", length(tech.cumul.gene.numbers)), cell.num= seq(1:length(tech.cumul.gene.numbers))))
#save(HEK.plot.df, file ="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/Cumulative_gene_final/Cumulative_gene_dist_DS10K_HEK_dataPlot.RData")

dim(ChromNuclei.DS10K.HEKs)
ChromNuclei.sampled.HEKs <- sample(colnames(ChromNuclei.DS10K.HEKs),46)
ChromNuclei.DS10K.HEKs.sampled <- ChromNuclei.DS10K.HEKs[,ChromNuclei.sampled.HEKs]
ChromNuclei.HEKs.expressed.genes <- which(rowSums(ChromNuclei.DS10K.HEKs.sampled) > 0)
ChromNuclei.DS10K.HEKs.sampled <- ChromNuclei.DS10K.HEKs.sampled[ChromNuclei.HEKs.expressed.genes,]

tech.cumul.gene.numbers <- rep(NA, ncol(ChromNuclei.DS10K.HEKs.sampled))
for (cell in 1:ncol(ChromNuclei.DS10K.HEKs.sampled)){
  print(cell)
  sample.size = cell
  sample.gene.numbers <- rep(NA,50)
  for (i in 1:50){
    selected.cells <- sample(colnames(ChromNuclei.DS10K.HEKs.sampled), sample.size)
    sample.gene.names <- rep(NA, 70000)
    for (cell2 in selected.cells){
      vec.start.point <- (length(sample.gene.names[which(is.na(sample.gene.names) == FALSE)])) +1
      detected.genes <- rownames(ChromNuclei.DS10K.HEKs.sampled[which(ChromNuclei.DS10K.HEKs.sampled[,cell2] > 0),])
      new.genes <- setdiff(detected.genes,sample.gene.names)
      vec.end.point <- vec.start.point + length(new.genes) -1
      sample.gene.names[vec.start.point:vec.end.point] <- new.genes
      #detected.genes <- rownames(DS.mat.HEKS[which(DS.mat.HEKS[,cell2] > 0),])
      #sample.gene.names <- c(sample.gene.names,detected.genes)
    }
    sample.gene.names <- na.omit(sample.gene.names)
    sample.gene.names.uniq <- unique(sample.gene.names)
    num.uniq.genes <- length(sample.gene.names.uniq)
    if (cell == 46){ ChromNuclei.all.cumul.genes.46 <-  sample.gene.names.uniq}
    sample.gene.numbers[i] <- num.uniq.genes
  }
  sample.gene.numbers <- na.omit(sample.gene.numbers)
  print(paste(cell," = ", mean(sample.gene.numbers), sep = ""))
  tech.cumul.gene.numbers[cell] <- mean(sample.gene.numbers)
}
HEK.plot.df <- rbind(HEK.plot.df, data.frame(Cumul= tech.cumul.gene.numbers, tech= rep("ChromNuclei", length(tech.cumul.gene.numbers)), cell.num= seq(1:length(tech.cumul.gene.numbers))))
#save(HEK.plot.df, file ="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/Cumulative_gene_final/Cumulative_gene_dist_DS10K_HEK_dataPlot.RData")


dim(ChromV3.DS10K.HEKs)
ChromV3.sampled.HEKs <- sample(colnames(ChromV3.DS10K.HEKs),50)
ChromV3.DS10K.HEKs.sampled <- ChromV3.DS10K.HEKs[,ChromV3.sampled.HEKs]
ChromV3.HEKs.expressed.genes <- which(rowSums(ChromV3.DS10K.HEKs.sampled) > 0)
ChromV3.DS10K.HEKs.sampled <- ChromV3.DS10K.HEKs.sampled[ChromV3.HEKs.expressed.genes,]

tech.cumul.gene.numbers <- rep(NA, ncol(ChromV3.DS10K.HEKs.sampled))
for (cell in 1:ncol(ChromV3.DS10K.HEKs.sampled)){
  print(cell)
  sample.size = cell
  sample.gene.numbers <- rep(NA,50)
  for (i in 1:50){
    selected.cells <- sample(colnames(ChromV3.DS10K.HEKs.sampled), sample.size)
    sample.gene.names <- rep(NA, 70000)
    for (cell2 in selected.cells){
      vec.start.point <- (length(sample.gene.names[which(is.na(sample.gene.names) == FALSE)])) +1
      detected.genes <- rownames(ChromV3.DS10K.HEKs.sampled[which(ChromV3.DS10K.HEKs.sampled[,cell2] > 0),])
      new.genes <- setdiff(detected.genes,sample.gene.names)
      vec.end.point <- vec.start.point + length(new.genes) -1
      sample.gene.names[vec.start.point:vec.end.point] <- new.genes
      #detected.genes <- rownames(DS.mat.HEKS[which(DS.mat.HEKS[,cell2] > 0),])
      #sample.gene.names <- c(sample.gene.names,detected.genes)
    }
    sample.gene.names <- na.omit(sample.gene.names)
    sample.gene.names.uniq <- unique(sample.gene.names)
    num.uniq.genes <- length(sample.gene.names.uniq)
    if (cell == 46){ ChromV3.all.cumul.genes.46 <-  sample.gene.names.uniq}
    sample.gene.numbers[i] <- num.uniq.genes
  }
  sample.gene.numbers <- na.omit(sample.gene.numbers)
  print(paste(cell," = ", mean(sample.gene.numbers), sep = ""))
  tech.cumul.gene.numbers[cell] <- mean(sample.gene.numbers)
}
HEK.plot.df <- rbind(HEK.plot.df, data.frame(Cumul= tech.cumul.gene.numbers, tech= rep("ChromV3", length(tech.cumul.gene.numbers)), cell.num= seq(1:length(tech.cumul.gene.numbers))))
#save(HEK.plot.df, file ="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/Cumulative_gene_final/Cumulative_gene_dist_DS10K_HEK_dataPlot.RData")

png("ChromCompare_CumulativenGenes_HEK_ChromV2-21_ChromV3-Lib92_DS10K.png", width = 12, height = 12, units = "in", res = 600)
ggplot(HEK.plot.df, aes(x=cell.num, y=Cumul, group=tech)) + geom_point(aes(color=tech, size=3)) + geom_smooth(se = F, aes(color=tech), size=3) +
  theme_bw() + theme (axis.title.y = element_blank(), axis.text.y = element_text(angle = 45,size = 40, face = "bold", colour = "black"), axis.text.x = element_text(size = 40, face = "bold", angle= 50,hjust=1,vjust=1, colour = "black"), plot.title = element_text(size = 20, face = "bold"), panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
dev.off()


set.seed(123)
dim(ChromV2.all.DS10K.Monocytes)
ChromV2.sampled.Monocytes <- sample(colnames(ChromV2.all.DS10K.Monocytes),50)
ChromV2.all.DS10K.Monocytes.sampled <- ChromV2.all.DS10K.Monocytes[,ChromV2.sampled.Monocytes]
ChromV2.Monocytes.expressed.genes <- which(rowSums(ChromV2.all.DS10K.Monocytes.sampled) > 0)
ChromV2.all.DS10K.Monocytes.sampled <- ChromV2.all.DS10K.Monocytes.sampled[ChromV2.Monocytes.expressed.genes,]

Monocytes.plot.df <- data.frame()
tech.cumul.gene.numbers <- rep(NA, ncol(ChromV2.all.DS10K.Monocytes.sampled))
for (cell in 1:ncol(ChromV2.all.DS10K.Monocytes.sampled)){
  print(cell)
  sample.size = cell
  sample.gene.numbers <- rep(NA,50)
  for (i in 1:50){
    selected.cells <- sample(colnames(ChromV2.all.DS10K.Monocytes.sampled), sample.size)
    sample.gene.names <- rep(NA, 70000)
    for (cell2 in selected.cells){
      vec.start.point <- (length(sample.gene.names[which(is.na(sample.gene.names) == FALSE)])) +1
      detected.genes <- rownames(ChromV2.all.DS10K.Monocytes.sampled[which(ChromV2.all.DS10K.Monocytes.sampled[,cell2] > 0),])
      new.genes <- setdiff(detected.genes,sample.gene.names)
      vec.end.point <- vec.start.point + length(new.genes) -1
      sample.gene.names[vec.start.point:vec.end.point] <- new.genes
      #detected.genes <- rownames(DS.mat.Monocytes[which(DS.mat.Monocytes[,cell2] > 0),])
      #sample.gene.names <- c(sample.gene.names,detected.genes)
    }
    sample.gene.names <- na.omit(sample.gene.names)
    sample.gene.names.uniq <- unique(sample.gene.names)
    num.uniq.genes <- length(sample.gene.names.uniq)
    if (cell == 50){ ChromV2.all.cumul.genes.Monocytes.50 <-  sample.gene.names.uniq}
    sample.gene.numbers[i] <- num.uniq.genes
  }
  sample.gene.numbers <- na.omit(sample.gene.numbers)
  print(paste(cell," = ", mean(sample.gene.numbers), sep = ""))
  tech.cumul.gene.numbers[cell] <- mean(sample.gene.numbers)
}
Monocytes.plot.df <- rbind(Monocytes.plot.df, data.frame(Cumul= tech.cumul.gene.numbers, tech= rep("ChromV2", length(tech.cumul.gene.numbers)), cell.num= seq(1:length(tech.cumul.gene.numbers))))
#save(Monocytes.plot.df, file ="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/Cumulative_gene_final/Cumulative_gene_dist_DS10K_Monocytes_dataPlot.RData")

dim(ChromNuclei.DS10K.Monocytes)
ChromNuclei.sampled.Monocytes <- sample(colnames(ChromNuclei.DS10K.Monocytes),50)
ChromNuclei.DS10K.Monocytes.sampled <- ChromNuclei.DS10K.Monocytes[,ChromNuclei.sampled.Monocytes]
ChromNuclei.Monocytes.expressed.genes <- which(rowSums(ChromNuclei.DS10K.Monocytes.sampled) > 0)
ChromNuclei.DS10K.Monocytes.sampled <- ChromNuclei.DS10K.Monocytes.sampled[ChromNuclei.Monocytes.expressed.genes,]

tech.cumul.gene.numbers <- rep(NA, ncol(ChromNuclei.DS10K.Monocytes.sampled))
for (cell in 1:ncol(ChromNuclei.DS10K.Monocytes.sampled)){
  print(cell)
  sample.size = cell
  sample.gene.numbers <- rep(NA,50)
  for (i in 1:50){
    selected.cells <- sample(colnames(ChromNuclei.DS10K.Monocytes.sampled), sample.size)
    sample.gene.names <- rep(NA, 70000)
    for (cell2 in selected.cells){
      vec.start.point <- (length(sample.gene.names[which(is.na(sample.gene.names) == FALSE)])) +1
      detected.genes <- rownames(ChromNuclei.DS10K.Monocytes.sampled[which(ChromNuclei.DS10K.Monocytes.sampled[,cell2] > 0),])
      new.genes <- setdiff(detected.genes,sample.gene.names)
      vec.end.point <- vec.start.point + length(new.genes) -1
      sample.gene.names[vec.start.point:vec.end.point] <- new.genes
      #detected.genes <- rownames(DS.mat.Monocytes[which(DS.mat.Monocytes[,cell2] > 0),])
      #sample.gene.names <- c(sample.gene.names,detected.genes)
    }
    sample.gene.names <- na.omit(sample.gene.names)
    sample.gene.names.uniq <- unique(sample.gene.names)
    num.uniq.genes <- length(sample.gene.names.uniq)
    if (cell == 50){ ChromNuclei.all.cumul.genes.Monocytes.50 <-  sample.gene.names.uniq}
    sample.gene.numbers[i] <- num.uniq.genes
  }
  sample.gene.numbers <- na.omit(sample.gene.numbers)
  print(paste(cell," = ", mean(sample.gene.numbers), sep = ""))
  tech.cumul.gene.numbers[cell] <- mean(sample.gene.numbers)
}
Monocytes.plot.df <- rbind(Monocytes.plot.df, data.frame(Cumul= tech.cumul.gene.numbers, tech= rep("ChromNuclei", length(tech.cumul.gene.numbers)), cell.num= seq(1:length(tech.cumul.gene.numbers))))
#save(Monocytes.plot.df, file ="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/Cumulative_gene_final/Cumulative_gene_dist_DS10K_Monocytes_dataPlot.RData")


dim(ChromV3.DS10K.Monocytes)
ChromV3.sampled.Monocytes <- sample(colnames(ChromV3.DS10K.Monocytes),50)
ChromV3.DS10K.Monocytes.sampled <- ChromV3.DS10K.Monocytes[,ChromV3.sampled.Monocytes]
ChromV3.Monocytes.expressed.genes <- which(rowSums(ChromV3.DS10K.Monocytes.sampled) > 0)
ChromV3.DS10K.Monocytes.sampled <- ChromV3.DS10K.Monocytes.sampled[ChromV3.Monocytes.expressed.genes,]

tech.cumul.gene.numbers <- rep(NA, ncol(ChromV3.DS10K.Monocytes.sampled))
for (cell in 1:ncol(ChromV3.DS10K.Monocytes.sampled)){
  print(cell)
  sample.size = cell
  sample.gene.numbers <- rep(NA,50)
  for (i in 1:50){
    selected.cells <- sample(colnames(ChromV3.DS10K.Monocytes.sampled), sample.size)
    sample.gene.names <- rep(NA, 70000)
    for (cell2 in selected.cells){
      vec.start.point <- (length(sample.gene.names[which(is.na(sample.gene.names) == FALSE)])) +1
      detected.genes <- rownames(ChromV3.DS10K.Monocytes.sampled[which(ChromV3.DS10K.Monocytes.sampled[,cell2] > 0),])
      new.genes <- setdiff(detected.genes,sample.gene.names)
      vec.end.point <- vec.start.point + length(new.genes) -1
      sample.gene.names[vec.start.point:vec.end.point] <- new.genes
      #detected.genes <- rownames(DS.mat.Monocytes[which(DS.mat.Monocytes[,cell2] > 0),])
      #sample.gene.names <- c(sample.gene.names,detected.genes)
    }
    sample.gene.names <- na.omit(sample.gene.names)
    sample.gene.names.uniq <- unique(sample.gene.names)
    num.uniq.genes <- length(sample.gene.names.uniq)
    if (cell == 50){ ChromV3.all.cumul.genes.Monocytes.50 <-  sample.gene.names.uniq}
    sample.gene.numbers[i] <- num.uniq.genes
  }
  sample.gene.numbers <- na.omit(sample.gene.numbers)
  print(paste(cell," = ", mean(sample.gene.numbers), sep = ""))
  tech.cumul.gene.numbers[cell] <- mean(sample.gene.numbers)
}
Monocytes.plot.df <- rbind(Monocytes.plot.df, data.frame(Cumul= tech.cumul.gene.numbers, tech= rep("ChromV3", length(tech.cumul.gene.numbers)), cell.num= seq(1:length(tech.cumul.gene.numbers))))
#save(Monocytes.plot.df, file ="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/Cumulative_gene_final/Cumulative_gene_dist_DS10K_Monocytes_dataPlot.RData")

png("ChromCompare_CumulativenGenes_Monocytes_ChromV2-21_ChromV3-Lib92_DS10K.png", width = 12, height = 12, units = "in", res = 600)
ggplot(Monocytes.plot.df, aes(x=cell.num, y=Cumul, group=tech)) + geom_point(aes(color=tech, size=3)) + geom_smooth(se = F, aes(color=tech), size = 3) +
  theme_bw() + theme (axis.title.y = element_blank(), axis.text.y = element_text(angle = 45,size = 40, face = "bold", colour = "black"), axis.text.x = element_text(size = 40, face = "bold", angle= 50,hjust=1,vjust=1, colour = "black"), plot.title = element_text(size = 20, face = "bold"), panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
dev.off()



set.seed(123)
dim(ChromV2.all.DS10K.B)
ChromV2.sampled.B <- sample(colnames(ChromV2.all.DS10K.B),50)
ChromV2.all.DS10K.B.sampled <- ChromV2.all.DS10K.B[,ChromV2.sampled.B]
ChromV2.B.expressed.genes <- which(rowSums(ChromV2.all.DS10K.B.sampled) > 0)
ChromV2.all.DS10K.B.sampled <- ChromV2.all.DS10K.B.sampled[ChromV2.B.expressed.genes,]

B.plot.df <- data.frame()
tech.cumul.gene.numbers <- rep(NA, ncol(ChromV2.all.DS10K.B.sampled))
for (cell in 1:ncol(ChromV2.all.DS10K.B.sampled)){
  print(cell)
  sample.size = cell
  sample.gene.numbers <- rep(NA,50)
  for (i in 1:50){
    selected.cells <- sample(colnames(ChromV2.all.DS10K.B.sampled), sample.size)
    sample.gene.names <- rep(NA, 70000)
    for (cell2 in selected.cells){
      vec.start.point <- (length(sample.gene.names[which(is.na(sample.gene.names) == FALSE)])) +1
      detected.genes <- rownames(ChromV2.all.DS10K.B.sampled[which(ChromV2.all.DS10K.B.sampled[,cell2] > 0),])
      new.genes <- setdiff(detected.genes,sample.gene.names)
      vec.end.point <- vec.start.point + length(new.genes) -1
      sample.gene.names[vec.start.point:vec.end.point] <- new.genes
      #detected.genes <- rownames(DS.mat.B[which(DS.mat.B[,cell2] > 0),])
      #sample.gene.names <- c(sample.gene.names,detected.genes)
    }
    sample.gene.names <- na.omit(sample.gene.names)
    sample.gene.names.uniq <- unique(sample.gene.names)
    num.uniq.genes <- length(sample.gene.names.uniq)
    if (cell == 13){ ChromV2.all.cumul.genes.B.13 <-  sample.gene.names.uniq}
    sample.gene.numbers[i] <- num.uniq.genes
  }
  sample.gene.numbers <- na.omit(sample.gene.numbers)
  print(paste(cell," = ", mean(sample.gene.numbers), sep = ""))
  tech.cumul.gene.numbers[cell] <- mean(sample.gene.numbers)
}
B.plot.df <- rbind(B.plot.df, data.frame(Cumul= tech.cumul.gene.numbers, tech= rep("ChromV2", length(tech.cumul.gene.numbers)), cell.num= seq(1:length(tech.cumul.gene.numbers))))
#save(B.plot.df, file ="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/Cumulative_gene_final/Cumulative_gene_dist_DS10K_B_dataPlot.RData")

dim(ChromNuclei.DS10K.B)
ChromNuclei.sampled.B <- sample(colnames(ChromNuclei.DS10K.B),50)
ChromNuclei.DS10K.B.sampled <- ChromNuclei.DS10K.B[,ChromNuclei.sampled.B]
ChromNuclei.B.expressed.genes <- which(rowSums(ChromNuclei.DS10K.B.sampled) > 0)
ChromNuclei.DS10K.B.sampled <- ChromNuclei.DS10K.B.sampled[ChromNuclei.B.expressed.genes,]

tech.cumul.gene.numbers <- rep(NA, ncol(ChromNuclei.DS10K.B.sampled))
for (cell in 1:ncol(ChromNuclei.DS10K.B.sampled)){
  print(cell)
  sample.size = cell
  sample.gene.numbers <- rep(NA,50)
  for (i in 1:50){
    selected.cells <- sample(colnames(ChromNuclei.DS10K.B.sampled), sample.size)
    sample.gene.names <- rep(NA, 70000)
    for (cell2 in selected.cells){
      vec.start.point <- (length(sample.gene.names[which(is.na(sample.gene.names) == FALSE)])) +1
      detected.genes <- rownames(ChromNuclei.DS10K.B.sampled[which(ChromNuclei.DS10K.B.sampled[,cell2] > 0),])
      new.genes <- setdiff(detected.genes,sample.gene.names)
      vec.end.point <- vec.start.point + length(new.genes) -1
      sample.gene.names[vec.start.point:vec.end.point] <- new.genes
      #detected.genes <- rownames(DS.mat.B[which(DS.mat.B[,cell2] > 0),])
      #sample.gene.names <- c(sample.gene.names,detected.genes)
    }
    sample.gene.names <- na.omit(sample.gene.names)
    sample.gene.names.uniq <- unique(sample.gene.names)
    num.uniq.genes <- length(sample.gene.names.uniq)
    if (cell == 13){ ChromNuclei.all.cumul.genes.B.13 <-  sample.gene.names.uniq}
    sample.gene.numbers[i] <- num.uniq.genes
  }
  sample.gene.numbers <- na.omit(sample.gene.numbers)
  print(paste(cell," = ", mean(sample.gene.numbers), sep = ""))
  tech.cumul.gene.numbers[cell] <- mean(sample.gene.numbers)
}
B.plot.df <- rbind(B.plot.df, data.frame(Cumul= tech.cumul.gene.numbers, tech= rep("ChromNuclei", length(tech.cumul.gene.numbers)), cell.num= seq(1:length(tech.cumul.gene.numbers))))
#save(B.plot.df, file ="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/Cumulative_gene_final/Cumulative_gene_dist_DS10K_B_dataPlot.RData")


dim(ChromV3.DS10K.B)
ChromV3.sampled.B <- sample(colnames(ChromV3.DS10K.B),13)
ChromV3.DS10K.B.sampled <- ChromV3.DS10K.B[,ChromV3.sampled.B]
ChromV3.B.expressed.genes <- which(rowSums(ChromV3.DS10K.B.sampled) > 0)
ChromV3.DS10K.B.sampled <- ChromV3.DS10K.B.sampled[ChromV3.B.expressed.genes,]

tech.cumul.gene.numbers <- rep(NA, ncol(ChromV3.DS10K.B.sampled))
for (cell in 1:ncol(ChromV3.DS10K.B.sampled)){
  print(cell)
  sample.size = cell
  sample.gene.numbers <- rep(NA,50)
  for (i in 1:50){
    selected.cells <- sample(colnames(ChromV3.DS10K.B.sampled), sample.size)
    sample.gene.names <- rep(NA, 70000)
    for (cell2 in selected.cells){
      vec.start.point <- (length(sample.gene.names[which(is.na(sample.gene.names) == FALSE)])) +1
      detected.genes <- rownames(ChromV3.DS10K.B.sampled[which(ChromV3.DS10K.B.sampled[,cell2] > 0),])
      new.genes <- setdiff(detected.genes,sample.gene.names)
      vec.end.point <- vec.start.point + length(new.genes) -1
      sample.gene.names[vec.start.point:vec.end.point] <- new.genes
      #detected.genes <- rownames(DS.mat.B[which(DS.mat.B[,cell2] > 0),])
      #sample.gene.names <- c(sample.gene.names,detected.genes)
    }
    sample.gene.names <- na.omit(sample.gene.names)
    sample.gene.names.uniq <- unique(sample.gene.names)
    num.uniq.genes <- length(sample.gene.names.uniq)
    if (cell == 13){ ChromV3.all.cumul.genes.B.13 <-  sample.gene.names.uniq}
    sample.gene.numbers[i] <- num.uniq.genes
  }
  sample.gene.numbers <- na.omit(sample.gene.numbers)
  print(paste(cell," = ", mean(sample.gene.numbers), sep = ""))
  tech.cumul.gene.numbers[cell] <- mean(sample.gene.numbers)
}
B.plot.df <- rbind(B.plot.df, data.frame(Cumul= tech.cumul.gene.numbers, tech= rep("ChromV3", length(tech.cumul.gene.numbers)), cell.num= seq(1:length(tech.cumul.gene.numbers))))
#save(B.plot.df, file ="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/Cumulative_gene_final/Cumulative_gene_dist_DS10K_B_dataPlot.RData")

png("ChromCompare_CumulativenGenes_B_ChromV2-21_ChromV3-Lib92_DS10K.png", width = 12, height = 12, units = "in", res = 600)
ggplot(B.plot.df, aes(x=cell.num, y=Cumul, group=tech)) + geom_point(aes(color=tech, size =3)) + geom_smooth(se = F, aes(color=tech), size=3) +
  theme_bw() + theme (axis.title.y = element_blank(), axis.text.y = element_text(angle = 45,size = 40, face = "bold", colour = "black"), axis.text.x = element_text(size = 40, face = "bold", angle= 50,hjust=1,vjust=1, colour = "black"), plot.title = element_text(size = 20, face = "bold"), panel.border = element_blank(), panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
dev.off()

###End###


#### Venn diagrams ####
areaV2.HEK <- length(ChromV2.all.cumul.genes.46)
areaV3.HEK <- length(ChromV3.all.cumul.genes.46)
areaSN.HEK <- length(ChromNuclei.all.cumul.genes.46)
nV2V3.HEK <- length(intersect(ChromV2.all.cumul.genes.46, ChromV3.all.cumul.genes.46))
nV3SN.HEK <- length(intersect(ChromV3.all.cumul.genes.46, ChromNuclei.all.cumul.genes.46))
nV2SN.HEK <- length(intersect(ChromV2.all.cumul.genes.46, ChromNuclei.all.cumul.genes.46))
nV2V3SN.HEK <- length(Reduce(intersect,list(ChromV3.all.cumul.genes.46, ChromV2.all.cumul.genes.46, ChromNuclei.all.cumul.genes.46)))

library(VennDiagram)
png("ChromCompare_VennDiagram_HEK.png", width = 12, height = 12, units = "in", res = 600)
draw.triple.venn(areaV2.HEK, areaV3.HEK, areaSN.HEK, nV2V3.HEK, nV3SN.HEK, nV2SN.HEK, nV2V3SN.HEK, category =c("ChromV2", "ChromV3", "ChromNuclei"), fill = c("red", "blue","green"), cex = 4 , cat.cex=0.4, cat.col = c("white", "white", "white"))  
dev.off()

areaV2.Monocytes <- length(ChromV2.all.cumul.genes.Monocytes.50)
areaV3.Monocytes <- length(ChromV3.all.cumul.genes.Monocytes.50)
areaSN.Monocytes <- length(ChromNuclei.all.cumul.genes.Monocytes.50)
nV2V3.Monocytes <- length(intersect(ChromV2.all.cumul.genes.Monocytes.50, ChromV3.all.cumul.genes.Monocytes.50))
nV3SN.Monocytes <- length(intersect(ChromV3.all.cumul.genes.Monocytes.50, ChromNuclei.all.cumul.genes.Monocytes.50))
nV2SN.Monocytes <- length(intersect(ChromV2.all.cumul.genes.Monocytes.50, ChromNuclei.all.cumul.genes.Monocytes.50))
nV2V3SN.Monocytes <- length(Reduce(intersect,list(ChromV3.all.cumul.genes.Monocytes.50, ChromV2.all.cumul.genes.Monocytes.50, ChromNuclei.all.cumul.genes.Monocytes.50)))

library(VennDiagram)
png("ChromCompare_VennDiagram_Monocytes.png", width = 12, height = 12, units = "in", res = 600)
draw.triple.venn(areaV2.Monocytes, areaV3.Monocytes, areaSN.Monocytes, nV2V3.Monocytes, nV3SN.Monocytes, nV2SN.Monocytes, nV2V3SN.Monocytes, category =c("ChromV2", "ChromV3", "ChromNuclei"), fill = c("red", "blue","green"), cex = 4 , cat.cex=0.4, cat.col = c("white", "white", "white"))  
dev.off()


areaV2.B <- length(ChromV2.all.cumul.genes.B.13)
areaV3.B <- length(ChromV3.all.cumul.genes.B.13)
areaSN.B <- length(ChromNuclei.all.cumul.genes.B.13)
nV2V3.B <- length(intersect(ChromV2.all.cumul.genes.B.13, ChromV3.all.cumul.genes.B.13))
nV3SN.B <- length(intersect(ChromV3.all.cumul.genes.B.13, ChromNuclei.all.cumul.genes.B.13))
nV2SN.B <- length(intersect(ChromV2.all.cumul.genes.B.13, ChromNuclei.all.cumul.genes.B.13))
nV2V3SN.B <- length(Reduce(intersect,list(ChromV3.all.cumul.genes.B.13, ChromV2.all.cumul.genes.B.13, ChromNuclei.all.cumul.genes.B.13)))

library(VennDiagram)
png("ChromCompare_VennDiagram_B.png", width = 12, height = 12, units = "in", res = 600)
draw.triple.venn(areaV2.B, areaV3.B, areaSN.B, nV2V3.B, nV3SN.B, nV2SN.B, nV2V3SN.B, category =c("ChromV2", "ChromV3", "ChromNuclei"), fill = c("red", "blue","green"), cex = 4 , cat.cex=0.4, cat.col = c("white", "white", "white"))  
dev.off()



#GO analysis
library(simpleGO)
library(Hmisc)

# V2 and V3 only
HEK.V2.V3 <- intersect(ChromV2.all.cumul.genes.46, ChromV3.all.cumul.genes.46)
HEK.V2.V3.minusSN <- setdiff(HEK.V2.V3, ChromNuclei.all.cumul.genes.46)
HEK.V2.V3.minusSN <- as.character(lapply(HEK.V2.V3.minusSN, function (x) unlist(strsplit(x, split = ".", fixed = T))[1]))
#HEK.V2.V3.minusSN <- mapIDs(HEK.V2.V3.minusSN)

Monocytes.V2.V3 <- intersect(ChromV2.all.cumul.genes.Monocytes.50, ChromV3.all.cumul.genes.Monocytes.50)
Monocytes.V2.V3.minusSN <- setdiff(Monocytes.V2.V3, ChromNuclei.all.cumul.genes.Monocytes.50)
Monocytes.V2.V3.minusSN <- as.character(lapply(Monocytes.V2.V3.minusSN, function (x) unlist(strsplit(x, split = ".", fixed = T))[1]))
#Monocytes.V2.V3.minusSN <- mapIDs(Monocytes.V2.V3.minusSN)

B.V2.V3 <- intersect(ChromV2.all.cumul.genes.B.13, ChromV3.all.cumul.genes.B.13)
B.V2.V3.minusSN <- setdiff(B.V2.V3, ChromNuclei.all.cumul.genes.B.13)
B.V2.V3.minusSN <- as.character(lapply(B.V2.V3.minusSN, function (x) unlist(strsplit(x, split = ".", fixed = T))[1]))
#B.V2.V3.minusSN <- mapIDs(B.V2.V3.minusSN)

V2.V3.minusSN.all.intersect <- Reduce(intersect, list(HEK.V2.V3.minusSN,Monocytes.V2.V3.minusSN, B.V2.V3.minusSN))
results.V2.V3.minusSN.all.intersect=simpleGO(gene.list=V2.V3.minusSN.all.intersect)
results.V2.V3.minusSN.all.intersect[[1]][1:3]

V2.V3.minusSN.all.union <- Reduce(union, list(HEK.V2.V3.minusSN,Monocytes.V2.V3.minusSN, B.V2.V3.minusSN))
results.V2.V3.minusSN.all.union=simpleGO_ati(gene.list=V2.V3.minusSN.all.union, excel.export = "GO_V2V3minusSN_all_union.xlsx")
results.V2.V3.minusSN.all.union[[1]][,1:3]

GO.V2.V3.minusSN.df <- data.frame(matrix(ncol = 4))
colnames(GO.V2.V3.minusSN.df) <- c("BH p-value", "Overlap", "GO Signature", "Genes")
for (i in 1:length(results.V2.V3.minusSN.all.union)){
  df <- results.V2.V3.minusSN.all.union[[i]]
  GO.V2.V3.minusSN.df <- rbind(GO.V2.V3.minusSN.df, df)
}
GO.V2.V3.minusSN.df <- GO.V2.V3.minusSN.df[-1,]
GO.V2.V3.minusSN.df.orderd <- GO.V2.V3.minusSN.df[order(GO.V2.V3.minusSN.df$`BH p-value`),]
write.xlsx(GO.V2.V3.minusSN.df.orderd, file= "GO_V2V3minusSN_union_SC.xlsx", colnames= T, row.names=T)


#SN only
HEK.both.V2v3 <- union(ChromV2.all.cumul.genes.46, ChromV3.all.cumul.genes.46)
HEK.SN.only <- setdiff(ChromNuclei.all.cumul.genes.46, HEK.both.V2v3)
HEK.SN.only <- as.character(lapply(HEK.SN.only, function (x) unlist(strsplit(x, split = ".", fixed = T))[1]))
#HEK.SN.only <- mapIDs(HEK.SN.only)

Monocytes.both.V2v3 <- union(ChromV2.all.cumul.genes.Monocytes.50, ChromV3.all.cumul.genes.Monocytes.50)
Monocytes.SN.only <- setdiff(ChromNuclei.all.cumul.genes.Monocytes.50, Monocytes.both.V2v3)
Monocytes.SN.only <- as.character(lapply(Monocytes.SN.only, function (x) unlist(strsplit(x, split = ".", fixed = T))[1]))
#Monocytes.SN.only <- mapIDs(Monocytes.SN.only)

B.both.V2v3 <- union(ChromV2.all.cumul.genes.B.13, ChromV3.all.cumul.genes.B.13)
B.SN.only <- setdiff(ChromNuclei.all.cumul.genes.B.13, B.both.V2v3)
B.SN.only <- as.character(lapply(B.SN.only, function (x) unlist(strsplit(x, split = ".", fixed = T))[1]))
#B.SN.only <- mapIDs(B.SN.only)


SN.only.all.intersect <- Reduce(intersect, list(HEK.SN.only,Monocytes.SN.only, B.SN.only))
results.SN.only.all.intersect=simpleGO(gene.list=SN.only.all.intersect)
results.SN.only.all.intersect[[1]][1:3]

SN.only.all.union <- Reduce(union, list(HEK.SN.only,Monocytes.SN.only, B.SN.only))
results.SN.only.all.union=simpleGO(gene.list=SN.only.all.union)
results.SN.only.all.union[[1]][1:3]

GO.SN.only.df <- data.frame(matrix(ncol = 4))
colnames(GO.SN.only.df) <- c("BH p-value", "Overlap", "GO Signature", "Genes")
for (i in 1:length(results.SN.only.all.union)){
  df <- results.SN.only.all.union[[i]]
  GO.SN.only.df <- rbind(GO.SN.only.df, df)
}
GO.SN.only.df <- GO.SN.only.df[-1,]
GO.SN.only.df.orderd <- GO.SN.only.df[order(GO.SN.only.df$`BH p-value`),]
write.xlsx(GO.SN.only.df.orderd, file= "GO_SNonly_union_SN.xlsx", colnames= T, row.names=T)


#All intersections:
HEK.all.common <- Reduce(intersect, list(ChromV2.all.cumul.genes.46, ChromV3.all.cumul.genes.46,ChromNuclei.all.cumul.genes.46))
HEK.all.common <- as.character(lapply(HEK.all.common, function (x) unlist(strsplit(x, split = ".", fixed = T))[1]))

Monocytes.all.common <- Reduce(intersect, list(ChromV2.all.cumul.genes.Monocytes.50, ChromV3.all.cumul.genes.Monocytes.50, ChromNuclei.all.cumul.genes.Monocytes.50))
Monocytes.all.common <- as.character(lapply(Monocytes.all.common, function (x) unlist(strsplit(x, split = ".", fixed = T))[1]))

B.all.common <- Reduce(intersect, list(ChromV2.all.cumul.genes.B.13, ChromV3.all.cumul.genes.B.13,ChromNuclei.all.cumul.genes.B.13))
B.all.common <- as.character(lapply(B.all.common, function (x) unlist(strsplit(x, split = ".", fixed = T))[1]))


all.common.intersect <- Reduce(intersect, list(HEK.all.common,Monocytes.all.common, B.all.common))
results.all.common.intersect=simpleGO(gene.list=all.common.intersect)
results.all.common.intersect[[1]][1:3]

#SN.only.all.union <- Reduce(union, list(HEK.all.common,Monocytes.all.common, B.all.common))
#results.SN.only.all.union=simpleGO(gene.list=SN.only.all.union)
#results.SN.only.all.union[[1]][1:3]

GO.all.common.df <- data.frame(matrix(ncol = 4))
colnames(GO.all.common.df) <- c("BH p-value", "Overlap", "GO Signature", "Genes")
for (i in 1:length(results.all.common.intersect)){
  df <- results.all.common.intersect[[i]]
  GO.all.common.df <- rbind(GO.all.common.df, df)
}
GO.all.common.df <- GO.all.common.df[-1,]
GO.all.common.df.orderd <- GO.all.common.df[order(GO.all.common.df$`BH p-value`),]
write.xlsx(GO.all.common.df.orderd, file= "GO_allChromSN_intersect.xlsx", colnames= T, row.names=T)

load('/Volumes/Ati-Archive/HCA/SC_Protocols/Version3/R_analysis/Biomart_hsap_mapping_table.RData')
mapIDs <- function (gene_set){
  mat.rownames <- as.character(lapply(gene_set, function (x) unlist(strsplit(x, split = ".", fixed = T))[1]))
  hsap.genes.with.ids <- intersect(mat.rownames, rownames(hsap.GID.mapping.final))
  hsap.rownames.mapped <- mat.rownames
  names(hsap.rownames.mapped) <- hsap.rownames.mapped
  hsap.rownames.mapped[hsap.genes.with.ids] <- hsap.GID.mapping.final[hsap.genes.with.ids, "hgnc_symbol"]
  hsap.rownames.mapped <- as.character(hsap.rownames.mapped)
  #hsap.rownames.mapped <- tolower(hsap.rownames.mapped)
  #test.set.mapped <- capitalize(hsap.rownames.mapped)
  return(hsap.rownames.mapped)
}

###END###


