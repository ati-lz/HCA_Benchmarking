#Revision

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

# Replicate as comparing processing units 
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/stepwise_analysis_final/Stepwise_DS_plotdf_HEK.RData")
head(DSth.df.HEK)
DSth.df.HEK <- DSth.df.HEK[which(DSth.df.HEK$DSthNum == 20000),]
DSth.df.HEK$DStech <- mapvalues(DSth.df.HEK$DStech, 
                                from = c("CELseq2","MARSseq","QUARTZseq","SCRBseq", "SMARTseq2", "C1HTsmall","C1HTmedium","Chromium", "ChromiumNuclei","ddSEQ","Dropseq","ICELL8","inDrop"), 
                                to=c("CEL-Seq2", "MARS-Seq", "Quartz-Seq2", "mcSCRB-Seq","Smart-Seq2", "C1HT-small", "C1HT-medium", "Chromium", "Chromium(sn)", "ddSEQ", "Drop-Seq", "ICELL8","inDrop")) #"SMARTseqFINAL",  "Smart-Seq2",

libraries <- lapply(rownames(DSth.df.HEK), function(x) unlist(strsplit(x, "_")))
libraries_f <- sapply(libraries, function(x) paste(x[3:length(x)-1], sep="_"))
replicate <- stri_join_list(libraries_f, sep = "_", collapse = NULL)
for (i in 1:length(replicate)){
  id = replicate[i]
  id_parts = unlist(strsplit(id, split ="_"))
  if (length(id_parts) == 4 & grepl("Col", id_parts[3])){
    replicate[i] <- id_parts[1]
  }
}

DSth.df.HEK <- cbind(DSth.df.HEK, replicate)
DSth.df.HEK.bk <- DSth.df.HEK

#Add Chromium Reference DSth.df.HEK as well
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/reference_10x/Ref10X_DSth_df_HEK.RData")
DSth.df.HEK <- rbind(DSth.df.HEK, DSth.df.HEK.Ref10X)

#fixing the confusion between the SCRB-seq replicate labels and Dropseq replicate leables
new.SCRBseq.levels <- unique(paste("SCRB-", DSth.df.HEK[which(DSth.df.HEK$DStech == "mcSCRB-Seq"), "replicate"], sep = ""))
levels(DSth.df.HEK$replicate) <- c(levels(DSth.df.HEK$replicate), new.SCRBseq.levels)
DSth.df.HEK[which(DSth.df.HEK$DStech == "mcSCRB-Seq"), "replicate"] <- paste("SCRB-", DSth.df.HEK[which(DSth.df.HEK$DStech == "mcSCRB-Seq"), "replicate"], sep = "")

#adding the SMARTseq replicates
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
rownames(DSth.df.HEK)[which(DSth.df.HEK$DStech == "Smart-Seq2")] <- strtrim(rownames(DSth.df.HEK[which(DSth.df.HEK$DStech == "Smart-Seq2"),]),39)
levels(DSth.df.HEK$replicate) <- c(levels(DSth.df.HEK$replicate), c("batch1", "batch2", "batch3", "batch4", "batch5"))
DSth.df.HEK[which(DSth.df.HEK$replicate == "allLanes" & DSth.df.HEK$DSthNum == 20000),"replicate"] <- SS2.batches[rownames(DSth.df.HEK[which(DSth.df.HEK$replicate == "allLanes" & DSth.df.HEK$DSthNum == 20000),])]


DSth.df.HEK <- DSth.df.HEK[-which(DSth.df.HEK$DStech =="Chromium"),]
DSth.df.HEK <- DSth.df.HEK[-which(DSth.df.HEK$DStech =="inDrop"),]
DSth.df.HEK$DStech <- factor(DSth.df.HEK$DStech)
DSth.df.HEK$replicate <- factor(DSth.df.HEK$replicate)

chromRef.color <- c("darkorange4")
names(chromRef.color) <- c("Reference-Chromium")
selected.color.3 <- c(selected.color, chromRef.color)

DSth.df.HEK.bk <- DSth.df.HEK
#DSth.df.HEK.bk.20K <- DSth.df.HEK.bk[which(DSth.df.HEK.bk$DSthNum == 20000),] # Cause I did this step at the beginning to filter only the 20000 cells
Num.DS20.HEK.each.Replicate <- colSums(table(DSth.df.HEK.bk$DStech, DSth.df.HEK.bk$replicate))
selected.batches <- names(Num.DS20.HEK.each.Replicate[which(Num.DS20.HEK.each.Replicate >= 5)])

DSth.df.HEK.bk.20K.final <- DSth.df.HEK.bk[which(DSth.df.HEK.bk$replicate %in% selected.batches),]
#Excluding inDrop and SCRBseq since after all, they left only with one replicate
DSth.df.HEK.bk.20K.final <- DSth.df.HEK.bk.20K.final[-which(DSth.df.HEK.bk.20K.final$DStech %in% c("inDrop", "mcSCRB-Seq")),]
DSth.df.HEK.bk.20K.final$replicate <- factor(DSth.df.HEK.bk.20K.final$replicate)
colSums(table(DSth.df.HEK.bk.20K.final$DStech, DSth.df.HEK.bk.20K.final$replicate))

pdf("replicates_nGene_replicatesHaving5cellsAtLeast_HEK_3.pdf", width = 30)
ggplot(DSth.df.HEK.bk.20K.final, aes(x=replicate, y=nGenes)) + geom_boxplot(aes(fill=DStech)) + facet_grid(. ~ DStech, scales = "free", space = "free") + theme_bw() + theme (axis.text.x = element_text(angle = 50, hjust = 1, size = 15, face = "bold"),axis.text.y = element_text(size = 15, face = "bold"), legend.position = "none",
                                                                                                                                                                              panel.spacing = unit(0.2, "lines"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = NA, color = "black"),
                                                                                                                                                                              axis.line = element_line(colour = "black"))+
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(500, 1000, 2000, 4000, 8000, 15000, 25000, 40000))  +
  scale_fill_manual(values=selected.color.3)+ ggtitle("HEK Downsampled to 20K")
dev.off()







####### BCELLs #######
# Replicate as comparing processing units 
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/stepwise_analysis_final/Stepwise_DS_plotdf_Bcells.RData")
head(DSth.df.Bcells)
DSth.df.B <- DSth.df.Bcells[which(DSth.df.Bcells$DSthNum == 20000),]
DSth.df.B$DStech <- mapvalues(DSth.df.B$DStech, 
                              from = c("CELseq2","MARSseq","QUARTZseq","SCRBseq", "SMARTseq2", "C1HTsmall","C1HTmedium","Chromium", "ChromiumNuclei","ddSEQ","Dropseq","ICELL8","inDrop"), 
                              to=c("CEL-Seq2", "MARS-Seq", "Quartz-Seq2", "mcSCRB-Seq","Smart-Seq2", "C1HT-small", "C1HT-medium", "Chromium", "Chromium(sn)", "ddSEQ", "Drop-Seq", "ICELL8","inDrop")) #"SMARTseqFINAL",  "Smart-Seq2",

libraries <- lapply(rownames(DSth.df.B), function(x) unlist(strsplit(x, "_")))
libraries_f <- sapply(libraries, function(x) paste(x[3:length(x)-1], sep="_"))
replicate <- stri_join_list(libraries_f, sep = "_", collapse = NULL)
for (i in 1:length(replicate)){
  id = replicate[i]
  id_parts = unlist(strsplit(id, split ="_"))
  if (length(id_parts) == 4 & grepl("Col", id_parts[3])){
    replicate[i] <- id_parts[1]
  }
}

DSth.df.B <- cbind(DSth.df.B, replicate)
DSth.df.B.bk <- DSth.df.B

#Add Chromium Reference DSth.df.B as well (*** Here for Bcells after 20K downsampling we have only 6 cells...)
load("/Volumes/Ati-Archive/HCA/SCE_Robjects/reference_10x/Ref10X_DSth_df_Bcells.RData")
#DSth.df.B.Ref10X$replicate <- as.character(DSth.df.B.Ref10X$replicate)
DSth.df.B <- rbind(DSth.df.B, DSth.df.B.Ref10X)

#fixing the confusion between the SCRB-seq replicate labels and Dropseq replicate leables
new.SCRBseq.levels <- unique(paste("SCRB-", DSth.df.B[which(DSth.df.B$DStech == "mcSCRB-Seq"), "replicate"], sep = ""))
levels(DSth.df.B$replicate) <- c(levels(DSth.df.B$replicate), new.SCRBseq.levels)
DSth.df.B[which(DSth.df.B$DStech == "mcSCRB-Seq"), "replicate"] <- paste("SCRB-", DSth.df.B[which(DSth.df.B$DStech == "mcSCRB-Seq"), "replicate"], sep = "")

#adding the SMARTseq replicates
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

rownames(DSth.df.B)[which(DSth.df.B$DStech == "Smart-Seq2")] <- strtrim(rownames(DSth.df.B[which(DSth.df.B$DStech == "Smart-Seq2"),]),39)
levels(DSth.df.B$replicate) <- c(levels(DSth.df.B$replicate), c("batch1", "batch2", "batch3", "batch4", "batch5"))
DSth.df.B[which(DSth.df.B$replicate == "allLanes" & DSth.df.B$DSthNum == 20000),"replicate"] <- SS2.batches[rownames(DSth.df.B[which(DSth.df.B$replicate == "allLanes" & DSth.df.B$DSthNum == 20000),])]

DSth.df.B <- DSth.df.B[-which(DSth.df.B$DStech =="Chromium"),]
DSth.df.B <- DSth.df.B[-which(DSth.df.B$DStech =="inDrop"),]
DSth.df.B$DStech <- factor(DSth.df.B$DStech)
DSth.df.B$replicate <- factor(DSth.df.B$replicate)

chromRef.color <- c("darkorange4")
names(chromRef.color) <- c("Reference-Chromium")
selected.color.3 <- c(selected.color, chromRef.color)
DSth.df.B.bk <- DSth.df.B
Num.DS20.B.each.Replicate <- colSums(table(DSth.df.B.bk$DStech, DSth.df.B.bk$replicate))
selected.batches <- names(Num.DS20.B.each.Replicate[which(Num.DS20.B.each.Replicate >= 5)])

DSth.df.B.bk.20K.final <- DSth.df.B.bk[which(DSth.df.B.bk$replicate %in% selected.batches),]
#Excluding inDrop since after all, they left only with one replicate
#DSth.df.B.bk.20K.final <- DSth.df.B.bk.20K.final[-which(DSth.df.B.bk.20K.final$DStech %in% c("inDrop")),]
DSth.df.B.bk.20K.final$replicate <- factor(DSth.df.B.bk.20K.final$replicate)
colSums(table(DSth.df.B.bk.20K.final$DStech, DSth.df.B.bk.20K.final$replicate))

pdf("replicates_nGene_replicatesHaving5cellsAtLeast_Bcells_3.pdf", width = 30)
ggplot(DSth.df.B.bk.20K.final, aes(x=replicate, y=nGenes)) + geom_boxplot(aes(fill=DStech)) + facet_grid(. ~ DStech, scales = "free", space = "free") + theme_bw() + theme (axis.text.x = element_text(angle = 50, hjust = 1, size = 15, face = "bold"),axis.text.y = element_text(size = 15, face = "bold"), legend.position = "none",
                                                                                                                                                                            panel.spacing = unit(0.2, "lines"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = NA, color = "black"),
                                                                                                                                                                            axis.line = element_line(colour = "black"))+
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(500, 1000, 2000, 4000, 8000, 15000, 25000, 40000))  +
  scale_fill_manual(values=selected.color.3)+ ggtitle("HEK Downsampled to 20K")
dev.off()


pdf("replicates_nUMI_2.pdf", width = 30)
ggplot(DSth.df.B.bk, aes(x=replicate, y=nUMIs)) + geom_boxplot(aes(fill=DStech)) + facet_grid(. ~ DStech, scales = "free", space = "free")+theme (axis.text.x = element_text(angle = 50, hjust = 1, size = 15, face = "bold"),axis.text.y = element_text(size = 15, face = "bold"), legend.position = "none",
                                                                                                                                                  panel.spacing = unit(0.2, "lines"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = NA, color = "black"),
                                                                                                                                                  axis.line = element_line(colour = "black"))+
  scale_y_continuous(trans = "log", labels= scales::comma, breaks = c(500, 1000, 2000, 4000, 8000, 15000, 25000, 40000))  +
  scale_fill_manual(values=selected.color)+ ggtitle("B Downsampled to 20K")
dev.off()

