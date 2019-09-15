Species.pct.df <- data.frame()
Species.pct.df.num <- data.frame()


Robjects <- c("CELseq2.hsap.full.SCE.Robj","MARSseq.hsap.full.SCE.Robj","QUARTZseq.hsap.full.SCE.Robj","SCRBseq.hsap.full.SCE.Robj",
"SMARTseqFINAL.hsap.full.SCE.Robj", "C1HTsmall.hsap.full.SCE.Robj","C1HTmedium.hsap.full.SCE.Robj","10X2x5K.hsap.full.SCE.Robj", "Nuclei10X.hsap.full.SCE.Robj",
"ddSEQ.hsap.full.SCE.Robj", "Dropseq.hsap.full.SCE.1000cellsPerPool.Robj", "ICELL8.hsap.full.SCE.Robj", "1CB.hsap.full.SCE.AllReads.Robj")
Seurat_objs.hsap <- c("celseq_seu.obj.RData","marsseq_seu.obj.RData", "quartzseq_seu.obj.RData", "scrbseq_seu.obj.RData",
"smartseq_seu.obj.RData", "c1ht.s_seu.obj.RData", "c1ht.m_seu.obj.RData", "chromium_seu.obj.RData", "nuclei_seu.obj.RData",
"ddseq_seu.obj.RData", "dropseq_seu.obj.RData", "icell8_seu.obj.RData", "onecb_seu.obj.RData")
Seurat_objs.mmus <- c("celseq_seu.obj.mmus.RData","marsseq_seu.obj.mmus.RData", "quartzseq_seu.obj.mmus.RData", "scrbseq_seu.obj.mmus.RData",
"smartseq_seu.obj.mmus.RData", "c1ht.s_seu.obj.mmus.RData", "c1ht.m_seu.obj.mmus.RData", "chromium_seu.obj.mmus.RData", "nuclei_seu.obj.mmus.RData",
"ddseq_seu.obj.mmus.RData", "dropseq_seu.obj.mmus.RData", "icell8_seu.obj.mmus.RData", "onecb_seu.obj.mmus.RData")


#Robjects = c("/Users/alafzi/Dropbox (Personal)/HCA_ELI/New_datasets/CELseq2.hsap.full.SCE.Robj", "/Users/alafzi/Dropbox (Personal)/HCA_ELI/New_datasets/MARSseq.hsap.full.SCE.Robj", "/Users/alafzi/Dropbox (Personal)/HCA_ELI/New_datasets/QUARTZseq.hsap.full.SCE.Robj")
#Seurat_objs.hsap = c("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/celseq_seu.obj.RData", "/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/marsseq_seu.obj.RData", "/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Human/quartzseq_seu.obj.RData")
#Seurat_objs.mmus = c("/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Mouse/celseq_seu.obj.RData", "/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Mouse/marsseq_seu.obj.RData", "/Users/alafzi/Dropbox (Personal)/HCA_ELI/Eli_analysis/Mouse/quartzseq_seu.obj.RData")

species_composition <- function (Protocol_name,protocol.metadata, protocol.cells.hsap, protocol.cells.mmus){
    protocol.filter.cells <- which(protocol.metadata$nTReads < 10000)
    protocol.metadata <- protocol.metadata[-protocol.filter.cells,]
    if (Protocol_name == "SCRBseq" | Protocol_name == "ICELL8"){
        protocol.cells.mmus <- which(protocol.metadata$Species == "Mouse")
    }
    protocol.cells.dog <- which(protocol.metadata$Species == "Dog")
    protocol.cells.doublet <- which(protocol.metadata$Species == "Doublet")
    #protocol.all.cells.filtered <- sum(length(protocol.cells.hsap), length(protocol.cells.mmus), length(protocol.cells.dog), length(protocol.cells.doublet))
    protocol.all.cells.filtered <- sum(length(protocol.cells.hsap), length(protocol.cells.mmus), length(protocol.cells.dog))

    protocol.pct.df <- data.frame(tech= rep(Protocol_name, 3), Species = c("Human", "Mouse", "Dog"),
                             pct= c(length(protocol.cells.hsap)/protocol.all.cells.filtered, length(protocol.cells.mmus)/protocol.all.cells.filtered,
                                    length(protocol.cells.dog)/protocol.all.cells.filtered))

    protocol.pct.df.num <- data.frame(tech= rep(Protocol_name, 3), Species = c("Human", "Mouse", "Dog"),
                             pct= c(length(protocol.cells.hsap), length(protocol.cells.mmus),
                                    length(protocol.cells.dog)))

    outlist <- list(protocol.pct.df,protocol.pct.df.num)

    return(outlist)
  
}



Species.pct.df <- data.frame()
Species.pct.df.num <- data.frame()

for(i in 1:length(Robjects)){
  load(Robjects[i])
  #protocol.name= unlist(strsplit(Robjects[i], split = "/"))[7]
  #protocol.name= unlist(strsplit(protocol.name, split = "\\."))[1]
  protocol.name= unlist(strsplit(Robjects[i], split = "."))[1]
  print(protocol.name)
  protocol.hsap <- full.SCE.hsap
  rm(full.SCE.hsap)
  protocol.metadata <- colData(protocol.hsap)
  rm(protocol.hsap)
  load(Seurat_objs.hsap[i])
  #seurat.obj.name.hsap <- unlist(strsplit(Seurat_objs.hsap[i],split = "/"))[8]
  #seurat.obj.name.hsap <- unlist(strsplit(seurat.obj.name.hsap,split = "_"))[1]
  seurat.obj.name.hsap <- unlist(strsplit(Seurat_objs.hsap[i],split = "_"))[1]
  seurat.obj.hsap <- get(seurat.obj.name.hsap) 
  protocol.cells.hsap <- seurat.obj.hsap@cell.names
  rm(seurat.obj.hsap)
  load(Seurat_objs.mmus[i])
  #seurat.obj.name.mmus <- unlist(strsplit(Seurat_objs.mmus[i],split = "/"))[8]
  #seurat.obj.name.mmus <- unlist(strsplit(seurat.obj.name.mmus,split = "_"))[1]
  seurat.obj.name.mmus <- unlist(strsplit(Seurat_objs.mmus[i],split = "_"))[1]
  seurat.obj.mmus <- get(seurat.obj.name.mmus)
  protocol.cells.mmus <- seurat.obj.mmus@cell.names
  rm(seurat.obj.mmus)
  protocol.out <- species_composition(protocol.name,protocol.metadata, protocol.cells.hsap, protocol.cells.mmus)
  Species.pct.df <- rbind(Species.pct.df, protocol.out[[1]])
  Species.pct.df.num <- rbind(Species.pct.df.num, protocol.out[[2]])
  
}


##############################################################################################################

library(data.table)
library(reshape2)
library(scales)
library(ggplot2)

Species.pct.df$Species <- factor(Species.pct.df$Species, levels = c("Dog", "Mouse", "Human"))
pdf("Species_composition_withoutDoublet.pdf")
ggplot(data=Species.pct.df, aes(tech, pct, fill= Species)) + geom_bar(stat = "identity", width = 0.5, position = "stack") + theme_bw() +
  theme (axis.text.x = element_text(angle = 50, hjust = 1, size = 15, face = "bold", colour = "black"),axis.text.y = element_text( hjust = 1, size = 20, face = "bold", colour = "black"),
           panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(labels = percent_format()) +scale_fill_manual(values= c("brown1", "skyblue2", "chartreuse4"))
dev.off()
