#Quartz-Seq2

load("/QUARTZseq.hsap.full.SCE.jointDSmat_1.Robj")
QUARTZseq.1 <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
QUARTZseq.DS.UMI.1 <- QUARTZseq.1$UMI
rm(QUARTZseq.1)

load("/QUARTZseq.hsap.full.SCE.jointDSmat_2.Robj")
QUARTZseq.2 <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
QUARTZseq.DS.UMI.2 <- QUARTZseq.2$UMI
rm(QUARTZseq.2)

load("/QUARTZseq.hsap.full.SCE.jointDSmat_3.Robj")
QUARTZseq.3 <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
QUARTZseq.DS.UMI.3 <- QUARTZseq.3$UMI
rm(QUARTZseq.3)

load("/QUARTZseq.hsap.full.SCE.jointDSmat_4.Robj")
QUARTZseq.4 <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
QUARTZseq.DS.UMI.4 <- QUARTZseq.4$UMI
rm(QUARTZseq.4)

load("/QUARTZseq.hsap.full.SCE.jointDSmat_5.Robj")
QUARTZseq.5 <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
QUARTZseq.DS.UMI.5 <- QUARTZseq.5$UMI
rm(QUARTZseq.5)

load("quartzseq_seu.obj.Rdata")
QUARTZseq.hsap.obj <- quartzseq
rm(quartzseq)
QUARTZseq.hsap.metadata <- QUARTZseq.hsap.obj@meta.data
rm(QUARTZseq.hsap.obj)

QUARTZseq.HEK <- rownames(QUARTZseq.hsap.metadata)[which(QUARTZseq.hsap.metadata$clean.id == "HEK")]
QUARTZseq.Monocytes <- rownames(QUARTZseq.hsap.metadata)[which(QUARTZseq.hsap.metadata$clean.id == "Monocytes")]
QUARTZseq.Bcells <- rownames(QUARTZseq.hsap.metadata)[which(QUARTZseq.hsap.metadata$clean.id == "B")]

repeats <- c("1","2", "3","4","5")
DSth.df.HEK <- data.frame()
for (rep in repeats){
  print(rep)
  rep.DS.UMI <- get(paste("QUARTZseq.DS.UMI.",rep, sep = ""))
  rep.HEK <- QUARTZseq.HEK
  for (DSth in names(rep.DS.UMI)){
    print(paste(rep, DSth, sep = "_"))
    colnames(rep.DS.UMI[[DSth]]) <- gsub(x = colnames(rep.DS.UMI[[DSth]]), pattern = "\\.", replacement = "_")
    colnames(rep.DS.UMI[[DSth]]) <- sub("human_","",colnames(rep.DS.UMI[[DSth]]))
    comm.cells <- intersect(rep.HEK, colnames(rep.DS.UMI[[DSth]]))
    DS.mat.HEKS <- rep.DS.UMI[[DSth]][, comm.cells]
    DS.gene.distribution.HEK <- colSums(DS.mat.HEKS[,]>0)
    DS.UMI.distribution.HEK <- colSums(DS.mat.HEKS)
    DS.labels <- rep(DSth, length(DS.gene.distribution.HEK))
    DSth_number = as.numeric(unlist(strsplit(DSth, "_"))[[2]])
    DSth.number.vec = rep(DSth_number, length(DS.gene.distribution.HEK))
    DS.labels <- rep(DSth, length(DS.gene.distribution.HEK))
    DS.reps <- rep(rep, length(DS.gene.distribution.HEK))
    DS.df <- data.frame(nGenes = DS.gene.distribution.HEK, nUMIs = DS.UMI.distribution.HEK, DSthNum = DSth.number.vec, DSrep = DS.reps)
    DSth.df.HEK <- rbind(DSth.df.HEK, DS.df)
    print(dim(DSth.df.HEK))
  }
}

DSth.df.Monocytes <- data.frame()
for (rep in repeats){
  print(rep)
  rep.DS.UMI <- get(paste("QUARTZseq.DS.UMI.",rep, sep = ""))
  rep.Monocytes <- QUARTZseq.Monocytes
  for (DSth in names(rep.DS.UMI)){
    print(paste(rep, DSth, sep = "_"))
    colnames(rep.DS.UMI[[DSth]]) <- gsub(x = colnames(rep.DS.UMI[[DSth]]), pattern = "\\.", replacement = "_")
    colnames(rep.DS.UMI[[DSth]]) <- sub("human_","",colnames(rep.DS.UMI[[DSth]]))
    comm.cells <- intersect(rep.Monocytes, colnames(rep.DS.UMI[[DSth]]))
    DS.mat.MonocytesS <- rep.DS.UMI[[DSth]][, comm.cells]
    DS.gene.distribution.Monocytes <- colSums(DS.mat.MonocytesS[,]>0)
    DS.UMI.distribution.Monocytes <- colSums(DS.mat.MonocytesS)
    DS.labels <- rep(DSth, length(DS.gene.distribution.Monocytes))
    DSth_number = as.numeric(unlist(strsplit(DSth, "_"))[[2]])
    DSth.number.vec = rep(DSth_number, length(DS.gene.distribution.Monocytes))
    DS.labels <- rep(DSth, length(DS.gene.distribution.Monocytes))
    DS.reps <- rep(rep, length(DS.gene.distribution.Monocytes))
    DS.df <- data.frame(nGenes = DS.gene.distribution.Monocytes, nUMIs = DS.UMI.distribution.Monocytes, DSthNum = DSth.number.vec, DSrep = DS.reps)
    DSth.df.Monocytes <- rbind(DSth.df.Monocytes, DS.df)
    print(dim(DSth.df.Monocytes))
  }
}


DSth.df.Bcells <- data.frame()
for (rep in repeats){
  print(rep)
  rep.DS.UMI <- get(paste("QUARTZseq.DS.UMI.",rep, sep = ""))
  rep.Bcells <- QUARTZseq.Bcells
  for (DSth in names(rep.DS.UMI)){
    print(paste(rep, DSth, sep = "_"))
    colnames(rep.DS.UMI[[DSth]]) <- gsub(x = colnames(rep.DS.UMI[[DSth]]), pattern = "\\.", replacement = "_")
    colnames(rep.DS.UMI[[DSth]]) <- sub("human_","",colnames(rep.DS.UMI[[DSth]]))
    comm.cells <- intersect(rep.Bcells, colnames(rep.DS.UMI[[DSth]]))
    DS.mat.BcellsS <- rep.DS.UMI[[DSth]][, comm.cells]
    DS.gene.distribution.Bcells <- colSums(DS.mat.BcellsS[,]>0)
    DS.UMI.distribution.Bcells <- colSums(DS.mat.BcellsS)
    DS.labels <- rep(DSth, length(DS.gene.distribution.Bcells))
    DSth_number = as.numeric(unlist(strsplit(DSth, "_"))[[2]])
    DSth.number.vec = rep(DSth_number, length(DS.gene.distribution.Bcells))
    DS.labels <- rep(DSth, length(DS.gene.distribution.Bcells))
    DS.reps <- rep(rep, length(DS.gene.distribution.Bcells))
    DS.df <- data.frame(nGenes = DS.gene.distribution.Bcells, nUMIs = DS.UMI.distribution.Bcells, DSthNum = DSth.number.vec, DSrep = DS.reps)
    DSth.df.Bcells <- rbind(DSth.df.Bcells, DS.df)
    print(dim(DSth.df.Bcells))
  }
}



head(DSth.df.HEK)
DSth.df.HEK$DSthNum <- factor(DSth.df.HEK$DSthNum)


png("/HEK_Repetative_DS_Boxplot_QUARTZseq_UMI.png", width = 6, height = 4.5, units = "in", res = 600)
ggplot(data=DSth.df.HEK, aes(x=DSthNum, y=nUMIs, fill=DSrep)) + geom_boxplot() +
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(500, 1000, 2000, 4000, 8000, 14000, 25000)) +theme (axis.text.x = element_text(angle = 50, hjust = 1, size = 20, face = "bold", colour = "black"),axis.text.y = element_text(size = 20, face = "bold", colour = "black"),
                                                                                                                          panel.spacing.x = unit(0, "null"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.title = element_blank(),
                                                                                                                          axis.line = element_line(colour = "black"),legend.position = "none") #+ ggtitle("QUARTZseq HEK")
dev.off()

png("/HEK_Repetative_DS_Boxplot_QUARTZseq_Gene.png", width = 6, height = 4.5, units = "in", res = 600)
ggplot(data=DSth.df.HEK, aes(x=DSthNum, y=nGenes, fill=DSrep)) + geom_boxplot() +
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(500, 1000, 2000, 4000,6000, 8000)) +theme (axis.text.x = element_text(angle = 50, hjust = 1, size = 20, face = "bold", colour = "black"),axis.text.y = element_text(size = 20, face = "bold", colour = "black"),
                                                                                                                 panel.spacing.x = unit(0, "null"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.title = element_blank(),
                                                                                                                 axis.line = element_line(colour = "black"),legend.position = "none") #+ ggtitle("QUARTZseq HEK")
dev.off()

head(DSth.df.Monocytes)
DSth.df.Monocytes$DSthNum <- factor(DSth.df.Monocytes$DSthNum)
png("/Monocytes_Repetative_DS_Boxplot_QUARTZseq_UMI.png", width = 6, height = 4.5, units = "in", res = 600)
ggplot(data=DSth.df.Monocytes, aes(x=DSthNum, y=nUMIs, fill=DSrep)) + geom_boxplot() +
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(500, 1000, 2000, 4000, 8000, 14000, 25000)) +theme (axis.text.x = element_text(angle = 50, hjust = 1, size = 20, face = "bold", colour = "black"),axis.text.y = element_text(size = 20, face = "bold", colour = "black"),
                                                                                                                          panel.spacing.x = unit(0, "null"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.title = element_blank(),
                                                                                                                          axis.line = element_line(colour = "black"),legend.position = "none")# + ggtitle("QUARTZseq Monocytes")
dev.off()
head(DSth.df.Monocytes)
DSth.df.Monocytes$DSthNum <- factor(DSth.df.Monocytes$DSthNum)
png("/Monocytes_Repetative_DS_Boxplot_QUARTZseq_Gene.png", width = 6, height = 4.5, units = "in", res = 600)
ggplot(data=DSth.df.Monocytes, aes(x=DSthNum, y=nGenes, fill=DSrep)) + geom_boxplot() +
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(500, 1000, 2000, 4000,6000, 8000)) +theme (axis.text.x = element_text(angle = 50, hjust = 1, size = 20, face = "bold", colour = "black"),axis.text.y = element_text(size = 20, face = "bold", colour = "black"),
                                                                                                                 panel.spacing.x = unit(0, "null"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.title = element_blank(),
                                                                                                                 axis.line = element_line(colour = "black"),legend.position = "none")# + ggtitle("QUARTZseq Monocytes")
dev.off()

head(DSth.df.Bcells)
DSth.df.Bcells$DSthNum <- factor(DSth.df.Bcells$DSthNum)
png("/Bcells_Repetative_DS_Boxplot_QUARTZseq_UMI.png", , width = 6, height = 4.5, units = "in", res = 600)
ggplot(data=DSth.df.Bcells, aes(x=DSthNum, y=nUMIs, fill=DSrep)) + geom_boxplot() +
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(500, 1000, 2000, 4000, 8000, 14000, 25000)) +theme (axis.text.x = element_text(angle = 50, hjust = 1, size = 20, face = "bold", colour = "black"),axis.text.y = element_text(size = 20, face = "bold", colour = "black"),
                                                                                                                          panel.spacing.x = unit(0, "null"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.title = element_blank(),
                                                                                                                          axis.line = element_line(colour = "black"),legend.position = "none")# + ggtitle("QUARTZseq Bcells")
dev.off()
head(DSth.df.Bcells)
DSth.df.Bcells$DSthNum <- factor(DSth.df.Bcells$DSthNum)
png("/Bcells_Repetative_DS_Boxplot_QUARTZseq_Gene.png", , width = 6, height = 4.5, units = "in", res = 600)
ggplot(data=DSth.df.Bcells, aes(x=DSthNum, y=nGenes, fill=DSrep)) + geom_boxplot() +
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(500, 1000, 2000, 4000,6000, 8000)) +theme (axis.text.x = element_text(angle = 50, hjust = 1, size = 20, face = "bold", colour = "black"),axis.text.y = element_text(size = 20, face = "bold", colour = "black"),
                                                                                                                 panel.spacing.x = unit(0, "null"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.title = element_blank(),
                                                                                                                 axis.line = element_line(colour = "black"),legend.position = "none")# + ggtitle("QUARTZseq Bcells")
dev.off()


#######################################################################

# MARS-seq

load("/MARSseq.hsap.full.SCE.jointDSmat_1.Robj")
MARSseq.1 <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
MARSseq.DS.UMI.1 <- MARSseq.1$UMI
rm(MARSseq.1)

load("/MARSseq.hsap.full.SCE.jointDSmat_2.Robj")
MARSseq.2 <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
MARSseq.DS.UMI.2 <- MARSseq.2$UMI
rm(MARSseq.2)

load("/MARSseq.hsap.full.SCE.jointDSmat_3.Robj")
MARSseq.3 <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
MARSseq.DS.UMI.3 <- MARSseq.3$UMI
rm(MARSseq.3)

load("/MARSseq.hsap.full.SCE.jointDSmat_4.Robj")
MARSseq.4 <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
MARSseq.DS.UMI.4 <- MARSseq.4$UMI
rm(MARSseq.4)

load("/MARSseq.hsap.full.SCE.jointDSmat_5.Robj")
MARSseq.5 <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
MARSseq.DS.UMI.5 <- MARSseq.5$UMI
rm(MARSseq.5)

load("marsseq_seu.obj.Rdata")
MARSseq.hsap.obj <- marsseq
rm(marsseq)
MARSseq.hsap.metadata <- MARSseq.hsap.obj@meta.data
rm(MARSseq.hsap.obj)

MARSseq.HEK <- rownames(MARSseq.hsap.metadata)[which(MARSseq.hsap.metadata$clean.id == "HEK")]
MARSseq.Monocytes <- rownames(MARSseq.hsap.metadata)[which(MARSseq.hsap.metadata$clean.id == "Monocytes")]
MARSseq.Bcells <- rownames(MARSseq.hsap.metadata)[which(MARSseq.hsap.metadata$clean.id == "B")]

repeats <- c("1","2", "3","4","5")
DSth.df.HEK <- data.frame()
for (rep in repeats){
  print(rep)
  rep.DS.UMI <- get(paste("MARSseq.DS.UMI.",rep, sep = ""))
  rep.HEK <- MARSseq.HEK
  for (DSth in names(rep.DS.UMI)){
    print(paste(rep, DSth, sep = "_"))
    colnames(rep.DS.UMI[[DSth]]) <- gsub(x = colnames(rep.DS.UMI[[DSth]]), pattern = "\\.", replacement = "_")
    colnames(rep.DS.UMI[[DSth]]) <- sub("human_","",colnames(rep.DS.UMI[[DSth]]))
    comm.cells <- intersect(rep.HEK, colnames(rep.DS.UMI[[DSth]]))
    DS.mat.HEKS <- rep.DS.UMI[[DSth]][, comm.cells]
    DS.gene.distribution.HEK <- colSums(DS.mat.HEKS[,]>0)
    DS.UMI.distribution.HEK <- colSums(DS.mat.HEKS)
    DS.labels <- rep(DSth, length(DS.gene.distribution.HEK))
    DSth_number = as.numeric(unlist(strsplit(DSth, "_"))[[2]])
    DSth.number.vec = rep(DSth_number, length(DS.gene.distribution.HEK))
    DS.labels <- rep(DSth, length(DS.gene.distribution.HEK))
    DS.reps <- rep(rep, length(DS.gene.distribution.HEK))
    DS.df <- data.frame(nGenes = DS.gene.distribution.HEK, nUMIs = DS.UMI.distribution.HEK, DSthNum = DSth.number.vec, DSrep = DS.reps)
    DSth.df.HEK <- rbind(DSth.df.HEK, DS.df)
    print(dim(DSth.df.HEK))
  }
}

DSth.df.Monocytes <- data.frame()
for (rep in repeats){
  print(rep)
  rep.DS.UMI <- get(paste("MARSseq.DS.UMI.",rep, sep = ""))
  rep.Monocytes <- MARSseq.Monocytes
  for (DSth in names(rep.DS.UMI)){
    print(paste(rep, DSth, sep = "_"))
    colnames(rep.DS.UMI[[DSth]]) <- gsub(x = colnames(rep.DS.UMI[[DSth]]), pattern = "\\.", replacement = "_")
    colnames(rep.DS.UMI[[DSth]]) <- sub("human_","",colnames(rep.DS.UMI[[DSth]]))
    comm.cells <- intersect(rep.Monocytes, colnames(rep.DS.UMI[[DSth]]))
    DS.mat.MonocytesS <- rep.DS.UMI[[DSth]][, comm.cells]
    DS.gene.distribution.Monocytes <- colSums(DS.mat.MonocytesS[,]>0)
    DS.UMI.distribution.Monocytes <- colSums(DS.mat.MonocytesS)
    DS.labels <- rep(DSth, length(DS.gene.distribution.Monocytes))
    DSth_number = as.numeric(unlist(strsplit(DSth, "_"))[[2]])
    DSth.number.vec = rep(DSth_number, length(DS.gene.distribution.Monocytes))
    DS.labels <- rep(DSth, length(DS.gene.distribution.Monocytes))
    DS.reps <- rep(rep, length(DS.gene.distribution.Monocytes))
    DS.df <- data.frame(nGenes = DS.gene.distribution.Monocytes, nUMIs = DS.UMI.distribution.Monocytes, DSthNum = DSth.number.vec, DSrep = DS.reps)
    DSth.df.Monocytes <- rbind(DSth.df.Monocytes, DS.df)
    print(dim(DSth.df.Monocytes))
  }
}


DSth.df.Bcells <- data.frame()
for (rep in repeats){
  print(rep)
  rep.DS.UMI <- get(paste("MARSseq.DS.UMI.",rep, sep = ""))
  rep.Bcells <- MARSseq.Bcells
  for (DSth in names(rep.DS.UMI)){
    print(paste(rep, DSth, sep = "_"))
    colnames(rep.DS.UMI[[DSth]]) <- gsub(x = colnames(rep.DS.UMI[[DSth]]), pattern = "\\.", replacement = "_")
    colnames(rep.DS.UMI[[DSth]]) <- sub("human_","",colnames(rep.DS.UMI[[DSth]]))
    comm.cells <- intersect(rep.Bcells, colnames(rep.DS.UMI[[DSth]]))
    DS.mat.BcellsS <- rep.DS.UMI[[DSth]][, comm.cells]
    DS.gene.distribution.Bcells <- colSums(DS.mat.BcellsS[,]>0)
    DS.UMI.distribution.Bcells <- colSums(DS.mat.BcellsS)
    DS.labels <- rep(DSth, length(DS.gene.distribution.Bcells))
    DSth_number = as.numeric(unlist(strsplit(DSth, "_"))[[2]])
    DSth.number.vec = rep(DSth_number, length(DS.gene.distribution.Bcells))
    DS.labels <- rep(DSth, length(DS.gene.distribution.Bcells))
    DS.reps <- rep(rep, length(DS.gene.distribution.Bcells))
    DS.df <- data.frame(nGenes = DS.gene.distribution.Bcells, nUMIs = DS.UMI.distribution.Bcells, DSthNum = DSth.number.vec, DSrep = DS.reps)
    DSth.df.Bcells <- rbind(DSth.df.Bcells, DS.df)
    print(dim(DSth.df.Bcells))
  }
}


head(DSth.df.HEK)
DSth.df.HEK$DSthNum <- factor(DSth.df.HEK$DSthNum)

png("/HEK_Repetative_DS_Boxplot_MARSseq_UMI.png", width = 6, height = 4.5, units = "in", res = 600)
ggplot(data=DSth.df.HEK, aes(x=DSthNum, y=nUMIs, fill=DSrep)) + geom_boxplot() +
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(500, 1000, 2000, 4000, 8000, 14000, 25000)) +theme (axis.text.x = element_text(angle = 50, hjust = 1, size = 20, face = "bold", colour = "black"),axis.text.y = element_text(size = 20, face = "bold", colour = "black"),
                                                                                                                          panel.spacing.x = unit(0, "null"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.title = element_blank(),
                                                                                                                          axis.line = element_line(colour = "black"),legend.position = "none") #+ ggtitle("MARSseq HEK")
dev.off()
png("/HEK_Repetative_DS_Boxplot_MARSseq_Gene.png", width = 6, height = 4.5, units = "in", res = 600)
ggplot(data=DSth.df.HEK, aes(x=DSthNum, y=nGenes, fill=DSrep)) + geom_boxplot() +
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(500, 1000, 2000, 4000,6000, 8000)) +theme (axis.text.x = element_text(angle = 50, hjust = 1, size = 20, face = "bold", colour = "black"),axis.text.y = element_text(size = 20, face = "bold", colour = "black"),
                                                                                                                 panel.spacing.x = unit(0, "null"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.title = element_blank(),
                                                                                                                 axis.line = element_line(colour = "black"),legend.position = "none") #+ ggtitle("MARSseq HEK")
dev.off()


head(DSth.df.Monocytes)
DSth.df.Monocytes$DSthNum <- factor(DSth.df.Monocytes$DSthNum)
png("/Monocytes_Repetative_DS_Boxplot_MARSseq_UMI.png", width = 6, height = 4.5, units = "in", res = 600)
ggplot(data=DSth.df.Monocytes, aes(x=DSthNum, y=nUMIs, fill=DSrep)) + geom_boxplot() +
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(500, 1000, 2000, 4000, 8000, 14000, 25000)) +theme (axis.text.x = element_text(angle = 50, hjust = 1, size = 20, face = "bold", colour = "black"),axis.text.y = element_text(size = 20, face = "bold", colour = "black"),
                                                                                                                          panel.spacing.x = unit(0, "null"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.title = element_blank(),
                                                                                                                          axis.line = element_line(colour = "black"),legend.position = "none")# + ggtitle("MARSseq Monocytes")
dev.off()
png("/Monocytes_Repetative_DS_Boxplot_MARSseq_Gene.png", width = 6, height = 4.5, units = "in", res = 600)
ggplot(data=DSth.df.Monocytes, aes(x=DSthNum, y=nGenes, fill=DSrep)) + geom_boxplot() +
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(500, 1000, 2000, 4000,6000, 8000)) +theme (axis.text.x = element_text(angle = 50, hjust = 1, size = 20, face = "bold", colour = "black"),axis.text.y = element_text(size = 20, face = "bold", colour = "black"),
                                                                                                                 panel.spacing.x = unit(0, "null"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.title = element_blank(),
                                                                                                                 axis.line = element_line(colour = "black"),legend.position = "none")# + ggtitle("MARSseq Monocytes")
dev.off()


head(DSth.df.Bcells)
DSth.df.Bcells$DSthNum <- factor(DSth.df.Bcells$DSthNum)
png("/Bcells_Repetative_DS_Boxplot_MARSseq_UMI.png", , width = 6, height = 4.5, units = "in", res = 600)
ggplot(data=DSth.df.Bcells, aes(x=DSthNum, y=nUMIs, fill=DSrep)) + geom_boxplot() +
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(500, 1000, 2000, 4000, 8000, 14000, 25000)) +theme (axis.text.x = element_text(angle = 50, hjust = 1, size = 20, face = "bold", colour = "black"),axis.text.y = element_text(size = 20, face = "bold", colour = "black"),
                                                                                                                          panel.spacing.x = unit(0, "null"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.title = element_blank(),
                                                                                                                          axis.line = element_line(colour = "black"),legend.position = "none")# + ggtitle("MARSseq Bcells")
dev.off()
png("/Bcells_Repetative_DS_Boxplot_MARSseq_Gene.png", , width = 6, height = 4.5, units = "in", res = 600)
ggplot(data=DSth.df.Bcells, aes(x=DSthNum, y=nGenes, fill=DSrep)) + geom_boxplot() +
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(500, 1000, 2000, 4000,6000, 8000)) +theme (axis.text.x = element_text(angle = 50, hjust = 1, size = 20, face = "bold", colour = "black"),axis.text.y = element_text(size = 20, face = "bold", colour = "black"),
                                                                                                                 panel.spacing.x = unit(0, "null"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.title = element_blank(),
                                                                                                                 axis.line = element_line(colour = "black"),legend.position = "none")# + ggtitle("MARSseq Bcells")
dev.off()



#######################################################################

# ddSEQ

load("/ddSEQ.hsap.full.SCE.jointDSmat_1.Robj")
ddSEQ.1 <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
ddSEQ.DS.UMI.1 <- ddSEQ.1$UMI
rm(ddSEQ.1)

load("/ddSEQ.hsap.full.SCE.jointDSmat_2.Robj")
ddSEQ.2 <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
ddSEQ.DS.UMI.2 <- ddSEQ.2$UMI
rm(ddSEQ.2)

load("/ddSEQ.hsap.full.SCE.jointDSmat_3.Robj")
ddSEQ.3 <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
ddSEQ.DS.UMI.3 <- ddSEQ.3$UMI
rm(ddSEQ.3)

load("/ddSEQ.hsap.full.SCE.jointDSmat_4.Robj")
ddSEQ.4 <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
ddSEQ.DS.UMI.4 <- ddSEQ.4$UMI
rm(ddSEQ.4)

load("/ddSEQ.hsap.full.SCE.jointDSmat_5.Robj")
ddSEQ.5 <- output.readcount.umicount.joint.mats
rm(output.readcount.umicount.joint.mats)
ddSEQ.DS.UMI.5 <- ddSEQ.5$UMI
rm(ddSEQ.5)

load("ddSEQ_seu.obj.Rdata")
ddSEQ.hsap.obj <- ddseq
rm(ddseq)
ddSEQ.hsap.metadata <- ddSEQ.hsap.obj@meta.data
rm(ddSEQ.hsap.obj)

ddSEQ.HEK <- rownames(ddSEQ.hsap.metadata)[which(ddSEQ.hsap.metadata$clean.id == "HEK")]
ddSEQ.Monocytes <- rownames(ddSEQ.hsap.metadata)[which(ddSEQ.hsap.metadata$clean.id == "Monocytes")]
ddSEQ.Bcells <- rownames(ddSEQ.hsap.metadata)[which(ddSEQ.hsap.metadata$clean.id == "B")]

repeats <- c("1","2", "3","4","5")
DSth.df.HEK <- data.frame()
for (rep in repeats){
  print(rep)
  rep.DS.UMI <- get(paste("ddSEQ.DS.UMI.",rep, sep = ""))
  rep.HEK <- ddSEQ.HEK
  for (DSth in names(rep.DS.UMI)){
    print(paste(rep, DSth, sep = "_"))
    colnames(rep.DS.UMI[[DSth]]) <- gsub(x = colnames(rep.DS.UMI[[DSth]]), pattern = "\\.", replacement = "_")
    colnames(rep.DS.UMI[[DSth]]) <- sub("human_","",colnames(rep.DS.UMI[[DSth]]))
    comm.cells <- intersect(rep.HEK, colnames(rep.DS.UMI[[DSth]]))
    DS.mat.HEKS <- rep.DS.UMI[[DSth]][, comm.cells]
    DS.gene.distribution.HEK <- colSums(DS.mat.HEKS[,]>0)
    DS.UMI.distribution.HEK <- colSums(DS.mat.HEKS)
    DS.labels <- rep(DSth, length(DS.gene.distribution.HEK))
    DSth_number = as.numeric(unlist(strsplit(DSth, "_"))[[2]])
    DSth.number.vec = rep(DSth_number, length(DS.gene.distribution.HEK))
    DS.labels <- rep(DSth, length(DS.gene.distribution.HEK))
    DS.reps <- rep(rep, length(DS.gene.distribution.HEK))
    DS.df <- data.frame(nGenes = DS.gene.distribution.HEK, nUMIs = DS.UMI.distribution.HEK, DSthNum = DSth.number.vec, DSrep = DS.reps)
    DSth.df.HEK <- rbind(DSth.df.HEK, DS.df)
    print(dim(DSth.df.HEK))
  }
}

DSth.df.Monocytes <- data.frame()
for (rep in repeats){
  print(rep)
  rep.DS.UMI <- get(paste("ddSEQ.DS.UMI.",rep, sep = ""))
  rep.Monocytes <- ddSEQ.Monocytes
  for (DSth in names(rep.DS.UMI)){
    print(paste(rep, DSth, sep = "_"))
    colnames(rep.DS.UMI[[DSth]]) <- gsub(x = colnames(rep.DS.UMI[[DSth]]), pattern = "\\.", replacement = "_")
    colnames(rep.DS.UMI[[DSth]]) <- sub("human_","",colnames(rep.DS.UMI[[DSth]]))
    comm.cells <- intersect(rep.Monocytes, colnames(rep.DS.UMI[[DSth]]))
    DS.mat.MonocytesS <- rep.DS.UMI[[DSth]][, comm.cells]
    DS.gene.distribution.Monocytes <- colSums(DS.mat.MonocytesS[,]>0)
    DS.UMI.distribution.Monocytes <- colSums(DS.mat.MonocytesS)
    DS.labels <- rep(DSth, length(DS.gene.distribution.Monocytes))
    DSth_number = as.numeric(unlist(strsplit(DSth, "_"))[[2]])
    DSth.number.vec = rep(DSth_number, length(DS.gene.distribution.Monocytes))
    DS.labels <- rep(DSth, length(DS.gene.distribution.Monocytes))
    DS.reps <- rep(rep, length(DS.gene.distribution.Monocytes))
    DS.df <- data.frame(nGenes = DS.gene.distribution.Monocytes, nUMIs = DS.UMI.distribution.Monocytes, DSthNum = DSth.number.vec, DSrep = DS.reps)
    DSth.df.Monocytes <- rbind(DSth.df.Monocytes, DS.df)
    print(dim(DSth.df.Monocytes))
  }
}


DSth.df.Bcells <- data.frame()
for (rep in repeats){
  print(rep)
  rep.DS.UMI <- get(paste("ddSEQ.DS.UMI.",rep, sep = ""))
  rep.Bcells <- ddSEQ.Bcells
  for (DSth in names(rep.DS.UMI)){
    print(paste(rep, DSth, sep = "_"))
    colnames(rep.DS.UMI[[DSth]]) <- gsub(x = colnames(rep.DS.UMI[[DSth]]), pattern = "\\.", replacement = "_")
    colnames(rep.DS.UMI[[DSth]]) <- sub("human_","",colnames(rep.DS.UMI[[DSth]]))
    comm.cells <- intersect(rep.Bcells, colnames(rep.DS.UMI[[DSth]]))
    DS.mat.BcellsS <- rep.DS.UMI[[DSth]][, comm.cells]
    DS.gene.distribution.Bcells <- colSums(DS.mat.BcellsS[,]>0)
    DS.UMI.distribution.Bcells <- colSums(DS.mat.BcellsS)
    DS.labels <- rep(DSth, length(DS.gene.distribution.Bcells))
    DSth_number = as.numeric(unlist(strsplit(DSth, "_"))[[2]])
    DSth.number.vec = rep(DSth_number, length(DS.gene.distribution.Bcells))
    DS.labels <- rep(DSth, length(DS.gene.distribution.Bcells))
    DS.reps <- rep(rep, length(DS.gene.distribution.Bcells))
    DS.df <- data.frame(nGenes = DS.gene.distribution.Bcells, nUMIs = DS.UMI.distribution.Bcells, DSthNum = DSth.number.vec, DSrep = DS.reps)
    DSth.df.Bcells <- rbind(DSth.df.Bcells, DS.df)
    print(dim(DSth.df.Bcells))
  }
}


head(DSth.df.HEK)
DSth.df.HEK$DSthNum <- factor(DSth.df.HEK$DSthNum)
#HEK_aggregate <- aggregate(DSth.df.HEK$nGenes, list(DSth.df.HEK$DSrep, DSth.df.HEK$DSthNum), mean)
#colnames(HEK_aggregate) <- c("DSrep", "DSthNum", "MeanNGenes")
#png("/Volumes/Ati-Archive/HCA/SCE_Robjects/stepwise_analysis_final/Final_plots/HEK_stepwise_DS_READSvsGENES_V3.png")
#ggplot(HEK_aggregate, aes(x=DSthNum, y=MeanNGenes, group =DSrep))  + geom_point( aes(color=DSrep),size =4) + geom_smooth(method = "lm", formula = y ~ log(x), se = F, aes(color=DSrep), size=1) +
#  scale_x_continuous(breaks = c(5000,10000,15000,20000,50000))# + scale_color_manual(values=selected.color)# geom_smooth(method = "lm", formula = y ~ log(x), se = F, aes(color=DStech), size=0.5)
#dev.off()

png("/HEK_Repetative_DS_Boxplot_ddSEQ_UMI.png", width = 6, height = 4.5, units = "in", res = 600)
ggplot(data=DSth.df.HEK, aes(x=DSthNum, y=nUMIs, fill=DSrep)) + geom_boxplot() +
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(500, 1000, 2000, 4000, 8000, 14000, 25000)) +theme (axis.text.x = element_text(angle = 50, hjust = 1, size = 20, face = "bold", colour = "black"),axis.text.y = element_text(size = 20, face = "bold", colour = "black"),
                                                                                                                          panel.spacing.x = unit(0, "null"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.title = element_blank(),
                                                                                                                          axis.line = element_line(colour = "black"),legend.position = "none") #+ ggtitle("ddSEQ HEK")
dev.off()
png("/HEK_Repetative_DS_Boxplot_ddSEQ_Gene.png", width = 6, height = 4.5, units = "in", res = 600)
ggplot(data=DSth.df.HEK, aes(x=DSthNum, y=nGenes, fill=DSrep)) + geom_boxplot() +
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(500, 1000, 2000, 4000,6000, 8000)) +theme (axis.text.x = element_text(angle = 50, hjust = 1, size = 20, face = "bold", colour = "black"),axis.text.y = element_text(size = 20, face = "bold", colour = "black"),
                                                                                                                 panel.spacing.x = unit(0, "null"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.title = element_blank(),
                                                                                                                 axis.line = element_line(colour = "black"),legend.position = "none") #+ ggtitle("ddSEQ HEK")
dev.off()



head(DSth.df.Monocytes)
DSth.df.Monocytes$DSthNum <- factor(DSth.df.Monocytes$DSthNum)
png("/Monocytes_Repetative_DS_Boxplot_ddSEQ_UMI.png", width = 6, height = 4.5, units = "in", res = 600)
ggplot(data=DSth.df.Monocytes, aes(x=DSthNum, y=nUMIs, fill=DSrep)) + geom_boxplot() +
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(500, 1000, 2000, 4000, 8000, 14000, 25000)) +theme (axis.text.x = element_text(angle = 50, hjust = 1, size = 20, face = "bold", colour = "black"),axis.text.y = element_text(size = 20, face = "bold", colour = "black"),
                                                                                                                          panel.spacing.x = unit(0, "null"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.title = element_blank(),
                                                                                                                          axis.line = element_line(colour = "black"),legend.position = "none")# + ggtitle("ddSEQ Monocytes")
dev.off()
png("/Monocytes_Repetative_DS_Boxplot_ddSEQ_Gene.png", width = 6, height = 4.5, units = "in", res = 600)
ggplot(data=DSth.df.Monocytes, aes(x=DSthNum, y=nGenes, fill=DSrep)) + geom_boxplot() +
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(500, 1000, 2000, 4000,6000, 8000)) +theme (axis.text.x = element_text(angle = 50, hjust = 1, size = 20, face = "bold", colour = "black"),axis.text.y = element_text(size = 20, face = "bold", colour = "black"),
                                                                                                                 panel.spacing.x = unit(0, "null"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.title = element_blank(),
                                                                                                                 axis.line = element_line(colour = "black"),legend.position = "none")# + ggtitle("ddSEQ Monocytes")
dev.off()

head(DSth.df.Bcells)
DSth.df.Bcells$DSthNum <- factor(DSth.df.Bcells$DSthNum)
png("/Bcells_Repetative_DS_Boxplot_ddSEQ_UMI.png", , width = 6, height = 4.5, units = "in", res = 600)
ggplot(data=DSth.df.Bcells, aes(x=DSthNum, y=nUMIs, fill=DSrep)) + geom_boxplot() +
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(500, 1000, 2000, 4000, 8000, 14000, 25000)) +theme (axis.text.x = element_text(angle = 50, hjust = 1, size = 20, face = "bold", colour = "black"),axis.text.y = element_text(size = 20, face = "bold", colour = "black"),
                                                                                                                          panel.spacing.x = unit(0, "null"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.title = element_blank(),
                                                                                                                          axis.line = element_line(colour = "black"),legend.position = "none")# + ggtitle("ddSEQ Bcells")
dev.off()
png("/Bcells_Repetative_DS_Boxplot_ddSEQ_Gene.png", , width = 6, height = 4.5, units = "in", res = 600)
ggplot(data=DSth.df.Bcells, aes(x=DSthNum, y=nGenes, fill=DSrep)) + geom_boxplot() +
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(500, 1000, 2000, 4000,6000, 8000)) +theme (axis.text.x = element_text(angle = 50, hjust = 1, size = 20, face = "bold", colour = "black"),axis.text.y = element_text(size = 20, face = "bold", colour = "black"),
                                                                                                                 panel.spacing.x = unit(0, "null"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.title = element_blank(),
                                                                                                                 axis.line = element_line(colour = "black"),legend.position = "none")# + ggtitle("ddSEQ Bcells")
dev.off()



