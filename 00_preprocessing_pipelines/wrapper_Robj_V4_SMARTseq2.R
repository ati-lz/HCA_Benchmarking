
"
Wrap all the output from previous steps of the pipline and produce a single SCE object for each sc technology

Usage:
Robject_wrapper.R --hsapExp <DATA_IN> --mmusExp <DATA_IN> --species <DATA_IN> --output_SCEobj <FILE_OUT> [--technology <sc_technology> ]

Options:
-hsapExp --hsapExp <DATA_IN>                list of expression object from Human mapping for all pools
-mmusExp --mmusExp <DATA_IN>                list of expression object from Mouse mapping for all pools
-species --species <DATA_IN>                list of species infor from species deconvolution step from mixed mapping for all pools
-Output --output_SCEobj <FILE_OUT>          The output location for two R SCE object that has expression data mapped to human ref and mouse Ref reparately
--technology <tech_id>                      The scRNASeq technology
-h --help                                   show this
-v --version                                print version and stop

" -> doc

suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(vcfR))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(cardelino))
#suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(data.table))

main <- function(hsapExp,mmusExp, species,output_SCEobj, technology) {
  
  #ng_list <- scan("/Volumes/Ati-Archive/HCA/donor_id/GeneNumber_list.txt", what="", sep=" ")
  #nUMI_list <- scan("/Volumes/Ati-Archive/HCA/donor_id/UMINumber_list.txt", what="", sep=" ")
  #nread_list <- scan("/Volumes/Ati-Archive/HCA/donor_id/ReadNumber_list.txt", what="", sep=" ")
  #species_list <- scan("/Volumes/Ati-Archive/HCA/donor_id/species_list.txt", what="", sep=" ")
  #hsapExp_list <- scan("/Volumes/Ati-Archive/HCA/donor_id/Hsap_expression_list.txt", what="", sep=" ")
  #mmusExp_list <- scan("/Volumes/Ati-Archive/HCA/donor_id/Mmus_expression_list.txt", what="", sep=" ")
  
  #hsapExp = "/Volumes/Ati-Archive/HCA/donor_id/Hsap_expression_list.txt"
  #mmusExp = "/Volumes/Ati-Archive/HCA/donor_id/Mmus_expression_list.txt"
  #nReads = "/Volumes/Ati-Archive/HCA/donor_id/ReadNumber_list.txt"
  #nUMI = "/Volumes/Ati-Archive/HCA/donor_id/UMINumber_list.txt"
  #nGene = "/Volumes/Ati-Archive/HCA/donor_id/GeneNumber_list.txt"
  #nFeatures = "/Volumes/Ati-Archive/HCA/donor_id/Features_list.txt"
  #species = "/Volumes/Ati-Archive/HCA/donor_id/species_list.txt"

  
  print("Wrapping starts")
  hsapExp_list <- scan(hsapExp, what="", sep=" ")
  mmusExp_list <- scan(mmusExp, what="", sep=" ")
  species_list <- scan(species, what="", sep=" ")
  number_of_samples <- length(hsapExp_list)
  
  #SCEobj.list <- vector("list", number_of_samples)
  hsap.ExpsMat.list <- list()
  hsap.metadata.list <- list()
  mmus.ExpsMat.list <- list()
  mmus.metadata.list <- list()
  for (sample in 1:number_of_samples){    
    sample.species.file <- read.table(species_list[sample])
    sample.hsapExp.obj <- readRDS(hsapExp_list[sample])
    sample.mmusExp.obj <- readRDS(mmusExp_list[sample])
    sample.species.file[,1] <- gsub(x=sample.species.file[,1], pattern="SMARTseqZUMI", replacement = "SMARTseq2")
    sample.ID.pre <- unlist(strsplit(unlist(strsplit(hsapExp_list[sample], "/"))[2],"_")) #for local 6
    print(sample.ID.pre)
    sample.ID <- paste(sample.ID.pre[-length(sample.ID.pre)], collapse = "_")
    print(sample.ID)
    
    #sample.exon.reads <- as.data.frame(sample.features.file[which(sample.features.file$type == "exon"), "N"], row.names = create_cell_IDs(sample.features.file[which(sample.features.file$type == "exon"), "RG"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)); colnames(sample.exon.reads) <- "nExonReads" 
    #sample.intron.reads <- as.data.frame(sample.features.file[which(sample.features.file$type == "intron"), "N"], row.names = create_cell_IDs(sample.features.file[which(sample.features.file$type == "intron"), "RG"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)); colnames(sample.exon.reads) <- "nIntronReads" 
    #sample.intergenic.reads <- as.data.frame(sample.features.file[which(sample.features.file$type == "Intergenic"), "N"], row.names = create_cell_IDs(sample.features.file[which(sample.features.file$type == "Intergenic"), "RG"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)); colnames(sample.exon.reads) <- "nIntergenicReads" 
    #sample.Unmapped.reads <- as.data.frame(sample.features.file[which(sample.features.file$type == "Unmapped"), "N"], row.names = create_cell_IDs(sample.features.file[which(sample.features.file$type == "Unmapped"), "RG"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)); colnames(sample.exon.reads) <- "nUnmappedReads" 
    #sample.Ambiguity.reads <- as.data.frame(sample.features.file[which(sample.features.file$type == "Ambiguity"), "N"], row.names = create_cell_IDs(sample.features.file[which(sample.features.file$type == "Ambiguity"), "RG"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)); colnames(sample.exon.reads) <- "nAmbiguityReads" 
    #sample.ng <- as.data.frame(sample.ng.file[which(sample.ng.file$featureType == "Intron+Exon"), "Count"], row.names = create_cell_IDs(sample.ng.file[which(sample.ng.file$featureType == "Intron+Exon"), "SampleID"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)); colnames(sample.ng) <- "nGenes"
    #sample.nUMI <- as.data.frame(sample.nUMI.file[which(sample.nUMI.file$featureType == "Intron+Exon"), "Count"], row.names = create_cell_IDs(sample.nUMI.file[which(sample.nUMI.file$featureType == "Intron+Exon"), "SampleID"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)); colnames(sample.nUMI) <- "nUMIs"
    #sample.nread <- as.data.frame(sample.nread.file$Total, row.names = create_cell_IDs(sample.nread.file$RG, id.type = "cell_Barcode",tech = technology, lib = sample.ID)); colnames(sample.nread) <- "nReads"
    #sample.species <- as.data.frame(sample.species.file[,2], row.names = create_cell_IDs(sample.species.file, id.type = "standard"), col.names = c("species")); colnames(sample.species) <- "species"
    #sample.species <- as.data.frame(sample.species[rownames(sample.nread),],row.names = rownames(sample.nread)); colnames(sample.species) <- "Species"
    
    sample.hsap.readMat.exon <- as.matrix(sample.hsapExp.obj$readcount$exon$all)
    sample.hsap.readMat.intron <- as.matrix(sample.hsapExp.obj$readcount$intron$all)
    sample.hsap.readMat.inex <- as.matrix(sample.hsapExp.obj$readcount$inex$all)
    sample.hsap.exon.reads <- data.frame(nExonReads=colSums(sample.hsap.readMat.exon), cellID=create_cell_IDs(colnames(sample.hsap.readMat.exon), id.type = "cell_Barcode",tech = technology, lib = sample.ID))#; colnames(sample.exon.reads) <- "nExonReads" 
    sample.hsap.intron.reads <- data.frame(nIntronReads=colSums(sample.hsap.readMat.intron), cellID=create_cell_IDs(colnames(sample.hsap.readMat.intron), id.type = "cell_Barcode",tech = technology, lib = sample.ID))#; colnames(sample.exon.reads) <- "nExonReads"  
    #sample.hsap.intergenic.reads <- data.frame(nIntergenicReads=sample.hfeatures.file[which(sample.hfeatures.file$type == "Intergenic"), "N"], cellID = create_cell_IDs(sample.hfeatures.file[which(sample.hfeatures.file$type == "Intergenic"), "RG"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)) 
    #sample.hsap.Unmapped.reads <- data.frame(nUnmappedReads=sample.hfeatures.file[which(sample.hfeatures.file$type == "Unmapped"), "N"], cellID = create_cell_IDs(sample.hfeatures.file[which(sample.hfeatures.file$type == "Unmapped"), "RG"], id.type = "cell_Barcode",tech = technology, lib = sample.ID))
    #sample.hsap.Ambiguity.reads <- data.frame(nAmbiguityReads=sample.hfeatures.file[which(sample.hfeatures.file$type == "Ambiguity"), "N"], cellID = create_cell_IDs(sample.hfeatures.file[which(sample.hfeatures.file$type == "Ambiguity"), "RG"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)) 
    #sample.hsap.Multimap.reads <- data.frame(nMultimapReads=sample.hfeatures.file[which(sample.hfeatures.file$type == "Multimapping"), "N"], cellID = create_cell_IDs(sample.hfeatures.file[which(sample.hfeatures.file$type == "Multimapping"), "RG"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)) 
    sample.hsap.nread <- data.frame(nTReads=colSums(sample.hsap.readMat.inex), cellID=create_cell_IDs(colnames(sample.hsap.readMat.inex), id.type = "cell_Barcode",tech = technology, lib = sample.ID))#; colnames(sample.exon.reads) 
    sample.hsap.ng <- data.frame(nGenes=colSums(sample.hsap.readMat.inex[,] > 0), cellID=create_cell_IDs(colnames(sample.hsap.readMat.inex), id.type = "cell_Barcode",tech = technology, lib = sample.ID))#; colnames(sample.exon.reads) 
    #sample.hsap.nUMI <- data.frame(nUMIs=sample.hnUMI.file[which(sample.hnUMI.file$type == "Intron+Exon"), "Count"], cellID = create_cell_IDs(sample.hnUMI.file[which(sample.hnUMI.file$type == "Intron+Exon"), "SampleID"], id.type = "cell_Barcode",tech = technology, lib = sample.ID))
    
    sample.mmus.readMat.exon <- as.matrix(sample.mmusExp.obj$readcount$exon$all)
    sample.mmus.readMat.intron <- as.matrix(sample.mmusExp.obj$readcount$intron$all)
    sample.mmus.readMat.inex <- as.matrix(sample.mmusExp.obj$readcount$inex$all)
    sample.mmus.exon.reads <- data.frame(nExonReads=colSums(sample.mmus.readMat.exon), cellID=create_cell_IDs(colnames(sample.mmus.readMat.exon), id.type = "cell_Barcode",tech = technology, lib = sample.ID))#; colnames(sample.exon.reads) <- "nExonReads" 
    sample.mmus.intron.reads <- data.frame(nIntronReads=colSums(sample.mmus.readMat.intron), cellID=create_cell_IDs(colnames(sample.mmus.readMat.intron), id.type = "cell_Barcode",tech = technology, lib = sample.ID))#; colnames(sample.exon.reads) <- "nExonReads"  
    #sample.mmus.intergenic.reads <- data.frame(nIntergenicReads=sample.hfeatures.file[which(sample.hfeatures.file$type == "Intergenic"), "N"], cellID = create_cell_IDs(sample.hfeatures.file[which(sample.hfeatures.file$type == "Intergenic"), "RG"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)) 
    #sample.mmus.Unmapped.reads <- data.frame(nUnmappedReads=sample.hfeatures.file[which(sample.hfeatures.file$type == "Unmapped"), "N"], cellID = create_cell_IDs(sample.hfeatures.file[which(sample.hfeatures.file$type == "Unmapped"), "RG"], id.type = "cell_Barcode",tech = technology, lib = sample.ID))
    #sample.mmus.Ambiguity.reads <- data.frame(nAmbiguityReads=sample.hfeatures.file[which(sample.hfeatures.file$type == "Ambiguity"), "N"], cellID = create_cell_IDs(sample.hfeatures.file[which(sample.hfeatures.file$type == "Ambiguity"), "RG"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)) 
    #sample.mmus.Multimap.reads <- data.frame(nMultimapReads=sample.hfeatures.file[which(sample.hfeatures.file$type == "Multimapping"), "N"], cellID = create_cell_IDs(sample.hfeatures.file[which(sample.hfeatures.file$type == "Multimapping"), "RG"], id.type = "cell_Barcode",tech = technology, lib = sample.ID)) 
    sample.mmus.nread <- data.frame(nTReads=colSums(sample.mmus.readMat.inex), cellID=create_cell_IDs(colnames(sample.mmus.readMat.inex), id.type = "cell_Barcode",tech = technology, lib = sample.ID))#; colnames(sample.exon.reads) 
    sample.mmus.ng <- data.frame(nGenes=colSums(sample.mmus.readMat.inex[,] > 0), cellID=create_cell_IDs(colnames(sample.mmus.readMat.inex), id.type = "cell_Barcode",tech = technology, lib = sample.ID))#; colnames(sample.exon.reads) 
    #sample.mmus.nUMI <- data.frame(nUMIs=sample.hnUMI.file[which(sample.hnUMI.file$type == "Intron+Exon"), "Count"], cellID = create_cell_IDs(sample.hnUMI.file[which(sample.hnUMI.file$type == "Intron+Exon"), "SampleID"], id.type = "cell_Barcode",tech = technology, lib = sample.ID))
    
    sample.species <- data.frame(Species=sample.species.file[,2], cellID = create_cell_IDs(sample.species.file, id.type = "standard",tech = technology, lib = sample.ID))
    sample.species <- sample.species[which(sample.species$cellID %in% sample.hsap.nread$cellID),]
    sample.species$cellID <- droplevels(sample.species$cellID)
    
    sample.library=data.frame(Library=rep(sample.ID, nrow(sample.hsap.nread)), cellID=sample.hsap.nread$cellID)
    
    #sample.hsap.metadata <- join_all(list(sample.hsap.nread, sample.hsap.exon.reads,sample.hsap.intron.reads, 
    #                                      sample.hsap.intergenic.reads, sample.hsap.Unmapped.reads, sample.hsap.Ambiguity.reads,
    #                                     sample.hsap.nUMI, sample.hsap.ng, sample.species, sample.library), by = "cellID", type = 'full')
    sample.hsap.metadata <- join_all(list(sample.hsap.nread, sample.hsap.exon.reads,sample.hsap.intron.reads, 
                                          sample.hsap.ng, sample.species, sample.library), by = "cellID", type = 'full')
    
    
    #print("is the problem here?")
    #print(class(sample.hfeatures.file))
    #print(head(sample.hfeatures.file))
    #print(dim(sample.hfeatures.file))
    #print(class(as.character(sample.hfeatures.file[which(sample.hfeatures.file$type == "Exon"), "RG"])))
    #cellIDtest=create_cell_IDs(as.character(sample.hfeatures.file[which(sample.hfeatures.file$type == "Exon"), "RG"]), id.type = "cell_Barcode",tech = technology, lib = sample.ID)
    #print(cellIDtest[1])
    #print(as.character(sample.hfeatures.file[which(sample.hfeatures.file$type == "Exon"), "RG"])[1:5])
    print(technology)
    print(sample.ID)
    
    rownames(sample.hsap.metadata) <- sample.hsap.metadata$cellID
    sample.hsap.metadata <- sample.hsap.metadata[,-2]
    sample.hsap.metadata <- na.omit(sample.hsap.metadata)
    hsap.metadata.list[[sample.ID]] <- sample.hsap.metadata
    
    print("or here?")
    sample.mmus.metadata <- join_all(list(sample.mmus.nread, sample.mmus.exon.reads,sample.mmus.intron.reads, 
                                          sample.mmus.ng, sample.species, sample.library), by = "cellID", type = 'full')
    rownames(sample.mmus.metadata) <- sample.mmus.metadata$cellID
    sample.mmus.metadata <- sample.mmus.metadata[,-2]
    sample.mmus.metadata <- na.omit(sample.mmus.metadata)
    mmus.metadata.list[[sample.ID]] <- sample.mmus.metadata
    
    
    #sample.hsapExp <- as.data.frame(as.matrix(sample.hsapExp.obj$Intron+Exon$umicounts))
    sample.hsapExp <- as.data.frame(sample.hsap.readMat.inex)
    colnames(sample.hsapExp) <- create_cell_IDs(colnames(sample.hsapExp), id.type = "cell_Barcode",tech = technology, lib = sample.ID)
    sample.hsapExp <- sample.hsapExp[, rownames(sample.hsap.metadata)]
    sample.hsapExp$rn <- rownames(sample.hsapExp)
    hsap.ExpsMat.list[[sample.ID]] <- sample.hsapExp
    #sample.hsapMetadata <- as.data.frame(cbind(nGenes=sample.nread, nUMIs=sample.nUMI, nReads=sample.nread, species=sample.species, library=rep(sample.ID, nrow(sample.nread))))
    
    sample.mmusExp <- as.data.frame(sample.mmus.readMat.inex)
    colnames(sample.mmusExp) <- create_cell_IDs(colnames(sample.mmusExp), id.type = "cell_Barcode",tech = technology, lib = sample.ID)
    sample.mmusExp <- sample.mmusExp[, rownames(sample.mmus.metadata)]
    sample.mmusExp$rn <- rownames(sample.mmusExp)
    mmus.ExpsMat.list[[sample.ID]] <- sample.mmusExp
    #sample.hsapMetadata <- as.data.frame(cbind(nGenes=sample.nread, nUMIs=sample.nUMI, nReads=sample.nread, species=sample.species, library=rep(sample.ID, nrow(sample.nread))))
    
    #sample.sce <- SingleCellExperiment(assays = list(counts = as.matrix(sample.mixedExp)), colData = sample.metadata)
    #SCEobj.list[[sample]] <- sample.sce
    #names(SCEobj.list)[sample] <- sample.ID
    
  }
  print("loop is done")
  all.hsapExp <- join_all(hsap.ExpsMat.list, by = "rn", type = 'full')
  print("hsap join is done")
  rownames(all.hsapExp) <- all.hsapExp$rn
  all.hsapExp <- all.hsapExp[,!(names(all.hsapExp) %in% c("rn"))]
  all.hsapExp[is.na(all.hsapExp)] <- 0
  print("hsap exp is done")
  
  all.mmusExp <- join_all(mmus.ExpsMat.list, by = "rn", type = 'full')
  print("mmus join is done")
  rownames(all.mmusExp) <- all.mmusExp$rn
  all.mmusExp <- all.mmusExp[,!(names(all.mmusExp) %in% c("rn"))]
  print("lets see if I can solve it by changing to matrix") ###### these two lines added for the caught segmentation error
  all.mmusExp <- as.matrix(all.mmusExp)
  all.mmusExp[is.na(all.mmusExp)] <- 0
  print("mmus exp is done")
  
  all.hsap.metadata <- do.call("rbind", hsap.metadata.list)
  print("hsap rbind metadata is done")
  print(rownames(all.hsap.metadata)[1:5])
  rownames(all.hsap.metadata) <- lapply(rownames(all.hsap.metadata), function (x) unlist(strsplit(x, "[.]"))[2])
  
  all.mmus.metadata <- do.call("rbind", mmus.metadata.list)
  print("mmus rbind metadata is done")
  rownames(all.mmus.metadata) <- lapply(rownames(all.mmus.metadata), function (x) unlist(strsplit(x, "[.]"))[2])
  
  print("Making the SingleCellExperiment Objects")
  print(length(colnames(all.hsapExp)))
  print(colnames(all.hsapExp)[1:10])
  print(length(rownames(all.hsap.metadata)))
  print(rownames(all.hsap.metadata)[1:10])
  full.SCE.hsap <- SingleCellExperiment(assays = list(counts = as.matrix(all.hsapExp)), colData = all.hsap.metadata[colnames(all.hsapExp),])
  full.SCE.hsap <- calculateQCMetrics(full.SCE.hsap)
  #libsize.drop.hsap <- isOutlier(full.SCE.hsap$total_counts, nmads=3, type="lower", log=TRUE)
  #feature.drop.hsap <- isOutlier(full.SCE.hsap$total_features, nmads=3, type="lower", log=TRUE)
  #full.SCE.filterred.hsap <- full.SCE.hsap[,!(libsize.drop.hsap | feature.drop.hsap)]
  #data.frame(ByLibSize=sum(libsize.drop.hsap), ByFeature=sum(feature.drop.hsap), Remaining=ncol(full.SCE.filterred.hsap))
  
  print(length(colnames(all.mmusExp)))
  print(colnames(all.mmusExp)[1:10])
  print(length(rownames(all.mmus.metadata)))
  print(rownames(all.mmus.metadata)[1:10])
  full.SCE.mmus <- SingleCellExperiment(assays = list(counts = as.matrix(all.mmusExp)), colData = all.mmus.metadata[colnames(all.mmusExp),])
  full.SCE.mmus <- calculateQCMetrics(full.SCE.mmus)
  #libsize.drop.mmus <- isOutlier(full.SCE.mmus$total_counts, nmads=3, type="lower", log=TRUE)
  #feature.drop.mmus <- isOutlier(full.SCE.mmus$total_features, nmads=3, type="lower", log=TRUE)
  #full.SCE.filterred.mmus <- full.SCE.mmus[,!(libsize.drop.mmus | feature.drop.mmus)]
  #data.frame(ByLibSize=sum(libsize.drop.mmus), ByFeature=sum(feature.drop.mmus), Remaining=ncol(full.SCE.filterred.mmus))
  
  #save(full.SCE.filterred.hsap, file = paste(output_SCEobj,"/hsap.full.SCE.Robj", sep = ""))
  #save(full.SCE.filterred.mmus, file = paste(output_SCEobj,"/mmus.full.SCE.Robj", sep = ""))
  save(full.SCE.hsap, file = paste(output_SCEobj,"/", technology,".hsap.full.SCE.Robj", sep = ""))
  save(full.SCE.mmus, file = paste(output_SCEobj,"/", technology,".mmus.full.SCE.Robj", sep = ""))
  
  return("Done")  
}

#counts(all.sce)[1:5,1:5]
#colData(all.sce)
#rowData(all.sce)

create_cell_IDs <- function(cell.IDs, id.type = "standard", tech = technology, lib = sample.ID){
  if (id.type == "standard"){
    pre.ids <- as.character(cell.IDs[,1])
    IDs.parts <- lapply(pre.ids, function(x) unlist(strsplit(x, split = ".", fixed = T)))
    new.IDs <- lapply(IDs.parts, function(x) paste(x[1], x[3], x[4], sep = "_"))
    return(as.character(new.IDs))}
  if (id.type == "cell_Barcode"){
    pre.ids <- cell.IDs
    new.IDs <- lapply(pre.ids, function(x) paste(tech, lib, x, sep = "_"))
    return(as.character(new.IDs))}
  if (id.type == "vcf_id"){
    pre.ids <- cell.IDs
    IDs.parts1 <- lapply(pre.ids, function(x) unlist(strsplit(x, split = "/", fixed = T))[[4]])
    IDs.parts2 <- lapply(IDs.parts1, function(x) unlist(strsplit(x, split = ".", fixed = T)))
    new.IDs <- lapply(IDs.parts2, function(x) paste(x[1], x[3], x[4], sep = "_"))
  }
}


## Get command line options
opt <- docopt::docopt(doc, version = "version 0.0.1\n")

#message("working directory: ", getwd(), "\n")
#message("input file: ", opt$input_file, "\n")
#message("output file: ", opt$output_file, "\n")


## Run main function
main(opt$hsapExp, opt$mmusExp, 
     opt$species, opt$output_SCEobj, opt$technology)

