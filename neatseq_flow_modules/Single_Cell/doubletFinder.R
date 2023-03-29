.libPaths()
library(optparse)
library(Seurat)
library(DoubletFinder)
library(dplyr)
library(ggplot2)
#library(plotly)
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--inputRDS"), type="character", default = NA,
              help="Path to Seurat object RDS file", metavar = "character"),
  make_option(c("-s", "--Sample"), type="character", default = NA,
              help="Sample's name", metavar = "character"),
  make_option(c("-d", "--dims"), type="numeric", default = 50,
              help="Dimensions to use as input for DoubletFinder (Default is 50). Must include --overwrite_dims", metavar = "character"),
  make_option(c("-o", "--outDir"), type="character", default = NA,
              help="Path to the output directory", metavar = "character"),
  make_option(c("--SCTransform"), action="store_true", default = FALSE,
              help="Will use SCTransform in DoubletFinder instead of standard log-normalization", metavar = "character"),
  make_option(c("--overwrite_dims"), action="store_true", default = FALSE,
              help="Run with manually number of dimensions (Set dims with --dims)", metavar = "character")
);

#Get user information
opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", 
                                    option_list=option_list,
                                    epilogue="\n\nAuthor:Gil Sorek");
opt = optparse::parse_args(opt_parser);

# Set log framework
logfile = paste(opt$outDir,opt$Sample,'_log.txt',sep='')
cat(paste('[',Sys.time(),']: DoubletFinder: ',opt$Sample,sep=''), file=logfile, sep='\n')
writeLog <- function(logfile, msg) { cat(paste('(',strftime(Sys.time(), format = "%H:%M:%S"),'): ',msg,sep=''),file=logfile,append=TRUE,sep='\n') }
saveRDS(opt, paste(opt$outDir,'opt.rds',sep=''))

print("Input var:",quote = F)
print.data.frame(as.data.frame(x = unlist(opt),row.names = names(opt)),right = F,quote = F)

# Functions
find_pk <- function(obj_seurat, num_dims, sct) {
  sweep.res.list <- paramSweep_v3(obj_seurat, PCs = 1:num_dims, sct = sct)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  plt = ggplot(bcmvn, aes(x=pK, y=BCmetric, group=1)) +
			geom_point() +
			geom_line()
  jpeg(paste(opt$outDir,opt$Sample, "_DoubletFinder_find_pK.jpeg", sep = ""), width=1000,height=1000)
  print(plt)
  dev.off()
  
  pK <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric),]$pK))
  nExp_poi <- round(0.075*nrow(obj_seurat@meta.data))
  annotations <- obj_seurat@active.ident ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) ## Assuming 7.5% doublet formation rate:
  df2return <- t(data.frame("pK" = pK, "nExp_poi" = nExp_poi, "nExp_poi.adj" = nExp_poi.adj))
  return(df2return)
}

# Import data
if (!is.na(opt$inputRDS)) {
  writeLog(logfile, paste("Importing Seurat RDS object...",sep=''))
  obj_seurat <- readRDS(opt$inputRDS)
} else {
  writeLog(logfile, paste("ERROR: Seurat object RDS file must be specified [--inputRDS]"))
  stop()
}
if (is.na(opt$Sample)) {
  writeLog(logfile, paste("ERROR: Samples name must be specified [--Sample]"))
  stop()
}
if (is.na(opt$outDir)) {
  writeLog(logfile, paste("ERROR: Output directory must be specified [--outDir]"))
  stop()
}

# Set number of dimensions for DoubletFinder
if (opt$SCTransform) {
  if (obj_seurat@active.assay == "SCT") {
    writeLog(logfile, paste("Active assay detected: ",obj_seurat@active.assay,sep=''))
    SigDims = opt$dims
  } else {
    writeLog(logfile, paste("ERROR: Active assay is not SCT"))
    stop()
  }
} else {
  writeLog(logfile, paste("Active assay detected: ",obj_seurat@active.assay,sep=''))
  if (opt$overwrite_dims) {
    SigDims = opt$dims
  } else {
    pAll <- JS(object = obj_seurat[["pca"]])
    nonSigDims = as.numeric(which(pAll$overall.p.values[,'Score'] > 0.05))
    if (length(nonSigDims) != 0) {
      SigDims = min(nonSigDims)-1
    } else {
      SigDims = dim(pAll$overall.p.values)[1]
    }
  }
}
writeLog(logfile, paste("Using ",SigDims," dimensions as input for DoubletFinder",sep=''))

# Calculate sample's DB parameters
writeLog(logfile, paste("Calculating DB parameters..."))
DB_params <- find_pk(obj_seurat, num_dims = SigDims, opt$SCTransform)
colnames(DB_params) <- opt$Sample
write.csv(DB_params, file = paste(opt$outDir,opt$Sample,"_DB_params.csv",sep=''), row.names = TRUE)

# Run DoubletFinder
writeLog(logfile, paste("Running DoubletFinder..."))
obj_seurat <- doubletFinder_v3(obj_seurat, PCs = 1:SigDims, pN = 0.25, pK = DB_params[1], 
                               nExp = DB_params[3], reuse.pANN = FALSE, sct = opt$SCTransform)

# Save table with cell names, classifications and scores
pANN <- paste0("pANN_0.25_", as.character(DB_params[1]),"_", as.character(DB_params[3]))
DF.classification <- paste0("DF.classifications_0.25_", as.character(DB_params[1]),"_", as.character(DB_params[3]))
DF_classification <- data.frame(row.names = row.names(obj_seurat@meta.data), 
                                'cluster'=obj_seurat@active.ident,
                                'pANN'=obj_seurat@meta.data[,pANN],
                              #  'DF.classification'=obj_seurat@meta.data[,DF.classification1],
                                'DF.classification_homotypic'=obj_seurat@meta.data[,DF.classification])
write.csv(DF_classification, 
          file = paste(opt$outDir,opt$Sample,"_DF_classifications.csv", sep = ""), row.names = TRUE)

# Save Summary table
DF_summary <- table(DF_classification$cluster, DF_classification$DF.classification_homotypic)
DF_summary <- prop.table(DF_summary, 1)*100
write.csv(DF_summary, 
          file = paste(opt$outDir,opt$Sample,"_DF_clusters_percentage.csv", sep = ""), row.names = TRUE)

# Plots
plt1 <- FeaturePlot(obj_seurat, reduction = "umap", pt.size = 0.8, features = pANN)
jpeg(paste(opt$outDir,opt$Sample, "_DF_pANN_FeaturePlot.jpeg", sep = ""),width = 1250, height = 1250)
print(plt1)
dev.off()
plt3 <- DimPlot(obj_seurat, reduction = "umap", pt.size = 0.8, group.by = DF.classification, cols = c("purple", "grey"))
jpeg(paste(opt$outDir,opt$Sample, "_DF.classification_homotypic.jpeg", sep = ""),width = 1250, height = 1250)
print(plt3)
dev.off()
plt4 <- VlnPlot(obj_seurat, features = pANN)
jpeg(paste(opt$outDir,opt$Sample, "_DF_pANN_ViolinPlot.jpeg", sep = ""),width = 1250, height = 1250)
print(plt4)
dev.off()
