
"
Runs scde for each technique

Usage:
dropout_wrapper_snakemake.R --technology <tech_id> --seuratObj_path <path_IN> --DS_path <DATA_IN> --output_path <FILE_OUT>

Options:
--technology <tech_id>                      The scRNASeq technology
--seuratObj_path <path_IN>                  path of the techs seurat object
--Monocyte_annot <DATA_IN>                  the way monocytes are annotated
--DS_path <DATA_IN>                         path to the techs downsampled data
--output_path <FILE_OUT>                    The output location for the results of the scde 

-h --help                                   show this
-v --version                                print version and stop

" -> doc

suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggExtra))
suppressPackageStartupMessages(library(cardelino))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(data.table))

main <- function(technology, seuratObj_path, DS_path, output_path) {

    load(DS_path)
    techs.DS <- output.readcount.umicount.joint.mats
    rm(output.readcount.umicount.joint.mats)
    if (technology == "SMARTseqFINAL"){
      techs.DS.UMI <- techs.DS$Reads
    } else {techs.DS.UMI <- techs.DS$UMI}
    rm(techs.DS)

    tech.DS.UMI <- techs.DS.UMI
    DSth = "downsampled_20000"
    colnames(tech.DS.UMI[[DSth]]) <- gsub(x = colnames(tech.DS.UMI[[DSth]]), pattern = "\\.", replacement = "_")
    DS.mat.20K <- tech.DS.UMI[[DSth]]
    DS.mat.20K <- mapIDs(DS.mat.20K, "hsap")


    seurat.obj = load(seuratObj_path)
    techs.obj = get(seurat.obj)
    rm(seurat.obj)
    techs.metadata <- techs.obj@meta.data
    common.seurat.DS20K.cells <- intersect(rownames(techs.metadata), colnames(DS.mat.20K))
    techs.metadata <- techs.metadata[common.seurat.DS20K.cells,]
  
    techs.HEK <- rownames(techs.metadata)[which(techs.metadata$clean.id == "HEK")]
    if(length(techs.HEK) > 50){ techs.HEK <- sample(techs.HEK,50)}
  
    techs.Monocytes <- rownames(techs.metadata)[which(techs.metadata$clean.id == "Monocytes")]
    if(length(techs.Monocytes) > 50){ techs.Monocytes <- sample(techs.Monocytes,50)}
  
    techs.Bcells <- rownames(techs.metadata)[which(techs.metadata$clean.id == "B")]
    if(length(techs.Bcells) > 50){ techs.Bcells <- sample(techs.Bcells,50)}
  
  
    DS.mat.HEKS <- DS.mat.20K[, techs.HEK]
    DS.mat.Monocytes <- DS.mat.20K[, techs.Monocytes]
    DS.mat.Bcells <- DS.mat.20K[, techs.Bcells]



    print("calculating the cumulatives For HEK")
    HEK.plot.df <- data.frame()
    tech.cumul.gene.numbers <- rep(NA, ncol(DS.mat.HEKS))
    for (cell in 1:ncol(DS.mat.HEKS)){
      print(cell)
      sample.size = cell
      sample.gene.numbers <- rep(NA,50)
      for (i in 1:50){
        selected.cells <- sample(colnames(DS.mat.HEKS), sample.size)
        sample.gene.names <- rep(NA, 70000)
        for (cell2 in selected.cells){
          vec.start.point <- (length(sample.gene.names[which(is.na(sample.gene.names) == FALSE)])) +1
          detected.genes <- rownames(DS.mat.HEKS[which(DS.mat.HEKS[,cell2] > 0),])
          new.genes <- setdiff(detected.genes,sample.gene.names)
          vec.end.point <- vec.start.point + length(new.genes) -1
          sample.gene.names[vec.start.point:vec.end.point] <- new.genes
          #detected.genes <- rownames(DS.mat.HEKS[which(DS.mat.HEKS[,cell2] > 0),])
          #sample.gene.names <- c(sample.gene.names,detected.genes)
          }
        sample.gene.names <- na.omit(sample.gene.names)
        num.uniq.genes <- length(unique(sample.gene.names))
        sample.gene.numbers[i] <- num.uniq.genes
      }
      sample.gene.numbers <- na.omit(sample.gene.numbers)
      print(paste(cell," = ", mean(sample.gene.numbers), sep = ""))
      tech.cumul.gene.numbers[cell] <- mean(sample.gene.numbers)
    }
  HEK.plot.df <- rbind(HEK.plot.df, data.frame(Cumul= tech.cumul.gene.numbers, tech= rep(technology, length(tech.cumul.gene.numbers)), cell.num= seq(1:length(tech.cumul.gene.numbers))))
  #save(HEK.plot.df, file ="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/Cumulative_gene_final/Cumulative_gene_dist_DS20K_HEK_dataPlot.RData")
  
  
  print("calculating the cumulatives For Monocytes")
  Monocytes.plot.df <- data.frame()
  tech.cumul.gene.numbers <- rep(NA, ncol(DS.mat.Monocytes))
  for (cell in 1:ncol(DS.mat.Monocytes)){
    print(cell)
    sample.size = cell
    sample.gene.numbers <- rep(NA,50)
    for (i in 1:50){
      selected.cells <- sample(colnames(DS.mat.Monocytes), sample.size)
      sample.gene.names <- rep(NA, 70000)
      for (cell2 in selected.cells){
        vec.start.point <- (length(sample.gene.names[which(is.na(sample.gene.names) == FALSE)])) +1
        detected.genes <- rownames(DS.mat.Monocytes[which(DS.mat.Monocytes[,cell2] > 0),])
        new.genes <- setdiff(detected.genes,sample.gene.names)
        vec.end.point <- vec.start.point + length(new.genes) -1
        sample.gene.names[vec.start.point:vec.end.point] <- new.genes
        #detected.genes <- rownames(DS.mat.Monocytes[which(DS.mat.Monocytes[,cell2] > 0),])
        #sample.gene.names <- c(sample.gene.names,detected.genes)
        }
      sample.gene.names <- na.omit(sample.gene.names)
      if (cell== ncol(DS.mat.Monocytes) & technology=="CELseq2"){CELseq2.monocytes.total.genenames <- sample.gene.names}
      if (cell== ncol(DS.mat.Monocytes) & technology=="QUARTZseq"){QUARTZseq.monocytes.total.genenames <- sample.gene.names}
      num.uniq.genes <- length(unique(sample.gene.names))
      sample.gene.numbers[i] <- num.uniq.genes
    }
    sample.gene.numbers <- na.omit(sample.gene.numbers)
    print(paste(cell," = ", mean(sample.gene.numbers), sep = ""))
    tech.cumul.gene.numbers[cell] <- mean(sample.gene.numbers)
  }
  Monocytes.plot.df <- rbind(Monocytes.plot.df, data.frame(Cumul= tech.cumul.gene.numbers, tech= rep(technology, length(tech.cumul.gene.numbers)), cell.num= seq(1:length(tech.cumul.gene.numbers))))
  #save(Monocytes.plot.df, file ="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/Cumulative_gene_final/Cumulative_gene_dist_DS20K_Monocytes_dataPlot.RData")

  print("calculating the cumulatives For Bcells")
  Bcells.plot.df <- data.frame()
  tech.cumul.gene.numbers <- rep(NA, ncol(DS.mat.Bcells))
  for (cell in 1:ncol(DS.mat.Bcells)){
    print(cell)
    sample.size = cell
    sample.gene.numbers <- rep(NA,50)
    for (i in 1:50){
      selected.cells <- sample(colnames(DS.mat.Bcells), sample.size)
      sample.gene.names <- rep(NA, 70000)
      for (cell2 in selected.cells){
        vec.start.point <- (length(sample.gene.names[which(is.na(sample.gene.names) == FALSE)])) +1
        detected.genes <- rownames(DS.mat.Bcells[which(DS.mat.Bcells[,cell2] > 0),])
        new.genes <- setdiff(detected.genes,sample.gene.names)
        vec.end.point <- vec.start.point + length(new.genes) -1
        sample.gene.names[vec.start.point:vec.end.point] <- new.genes
        #detected.genes <- rownames(DS.mat.Bcells[which(DS.mat.Bcells[,cell2] > 0),])
        #sample.gene.names <- c(sample.gene.names,detected.genes)
        }
      sample.gene.names <- na.omit(sample.gene.names)
      if (cell== ncol(DS.mat.Bcells) & technology=="CELseq2"){CELseq2.Bcells.total.genenames <- sample.gene.names}
      if (cell== ncol(DS.mat.Bcells) & technology=="QUARTZseq"){QUARTZseq.Bcells.total.genenames <- sample.gene.names}
      num.uniq.genes <- length(unique(sample.gene.names))
      sample.gene.numbers[i] <- num.uniq.genes
    }
    sample.gene.numbers <- na.omit(sample.gene.numbers)
    print(paste(cell," = ", mean(sample.gene.numbers), sep = ""))
    tech.cumul.gene.numbers[cell] <- mean(sample.gene.numbers)
  }
  Bcells.plot.df <- rbind(Bcells.plot.df, data.frame(Cumul= tech.cumul.gene.numbers, tech= rep(technology, length(tech.cumul.gene.numbers)), cell.num= seq(1:length(tech.cumul.gene.numbers))))
  #save(Bcells.plot.df, file ="/project/devel/alafzi/SC_Protocols/Version3/R_analysis/Cumulative_gene_final/Cumulative_gene_dist_DS20K_Bcells_dataPlot.RData")

  if (technology=="CELseq2"){output.list <- list(HEK.plot.df, Monocytes.plot.df, Bcells.plot.df, CELseq2.monocytes.total.genenames,CELseq2.Bcells.total.genenames)}
  else if (technology=="QUARTZseq"){output.list <- list(HEK.plot.df, Monocytes.plot.df, Bcells.plot.df,QUARTZseq.monocytes.total.genenames,QUARTZseq.Bcells.total.genenames)}
  else {output.list <- list(HEK.plot.df, Monocytes.plot.df, Bcells.plot.df)}
  
  save(output.list, file = paste(output_path,"/", technology,"_Cumulatives_V2.Robj", sep = ""))

}


load("Biomart_hsap_mapping_table.RData")
load("Biomart_mmus_mapping_table.RData")
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
#### end ####





## Get command line options
opt <- docopt::docopt(doc, version = "version 0.0.1\n")

#message("working directory: ", getwd(), "\n")
#message("input file: ", opt$input_file, "\n")
#message("output file: ", opt$output_file, "\n")


## Run main function
main(opt$technology, opt$seuratObj_path, opt$DS_path, opt$output_path)
