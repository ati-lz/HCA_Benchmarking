
"
Runs scde for each technique

Usage:
dropout_wrapper_snakemake.R --technology <tech_id> --seuratObj_path <path_IN>  --Monocyte_annot <DATA_IN> --DS_path <DATA_IN> --output_path <FILE_OUT>

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
suppressPackageStartupMessages(library(scde))
suppressPackageStartupMessages(library(data.table))

set.seed(123)
main <- function(technology, seuratObj_path, Monocyte_annot, DS_path, output_path) {
  
  load(seuratObj_path)
  techs.obj <- data
  rm(data)
  techs.metadata <- techs.obj@meta.data
  
  techs.HEK <- names(techs.obj@ident)[which(techs.obj@ident == "HEK cells")]
  if(length(techs.HEK) > 50){ techs.HEK <- sample(techs.HEK,50)}
  
  techs.Monocytes <- names(techs.obj@ident)[which(techs.obj@ident == Monocyte_annot)]
  if(length(techs.Monocytes) > 50){ techs.Monocytes <- sample(techs.Monocytes,50)}
  
  techs.Bcells <- names(techs.obj@ident)[which(techs.obj@ident == "B cells")]
  if(length(techs.Bcells) > 50){ techs.Bcells <- sample(techs.Bcells,50)}
  
  load(DS_path)
  techs.DS <- output.readcount.umicount.joint.mats
  rm(output.readcount.umicount.joint.mats)
  techs.DS.UMI <- techs.DS$UMI
  rm(techs.DS)
  
  tech.DS.UMI <- techs.DS.UMI
  tech.HEK <- techs.HEK
  DSth = "downsampled_20000"
  colnames(tech.DS.UMI[[DSth]]) <- gsub(x = colnames(tech.DS.UMI[[DSth]]), pattern = "\\.", replacement = "_")
  comm.cells <- intersect(tech.HEK, colnames(tech.DS.UMI[[DSth]]))
  DS.mat.HEKS <- tech.DS.UMI[[DSth]][, comm.cells]
  
  techs.raw.HEK <- DS.mat.HEKS
  techs.raw.HEK <- mapIDs(techs.raw.HEK, "hsap")
  techs.common.cells.HEK <- intersect(colnames(techs.raw.HEK),colnames(techs.obj@scale.data))
  techs.common.genes.HEK <- intersect(rownames(techs.raw.HEK),rownames(techs.obj@scale.data))
  techs.raw.HEK <- techs.raw.HEK[techs.common.genes.HEK,techs.common.cells.HEK]
  techs.raw.HEK <- as.data.frame(lapply(techs.raw.HEK, as.integer))
  
  techs.o.ifm <- scde.error.models(counts = techs.raw.HEK, n.cores = 10, threshold.segmentation = FALSE, save.crossfit.plots = F, save.model.plots = F, verbose = 1)
  #save(techs.o.ifm)
  print(paste(technology,"error model done", sep= "")
        
        #techs priors
        techs.o.prior <- scde.expression.prior(models = techs.o.ifm, counts = techs.raw.HEK, length.out = 400, show.plot = FALSE)
        #save(techs.o.prior)
        print(paste(technology,"Prior done", sep= "")
              
              #techs failure prob
              techs.self.fail <- scde.failure.probability(models = techs.o.ifm, counts = techs.raw.HEK)
              #save(techs.self.fail)
              print(paste(technology,"failure prob done", sep= "")
                    
                    techs.fail.curves <- scde.failure.probability(models = techs.o.ifm, magnitudes = log((10^techs.o.prior$x)-1))
                    #save(techs.fail.curves)
                    print(paste(technology,"fail curves done", sep= "")
                          
                          output.list <- list(techs.o.ifm, techs.o.prior, techs.self.fail, techs.fail.curves)
                          
                          save(output.list, file = paste(output_path,"/", technology,".SCDEoutputs.Robj", sep = ""))
                          
                          return("Done")
                          
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
main(opt$technology, opt$seuratObj_path, opt$Monocyte_annot, opt$DS_path, opt$output_path)
