.libPaths()
library(optparse)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--inputRDS"), type="character", default = NA,
              help="Path to Seurat object RDS file", metavar = "character"),
  make_option(c("--CPUs"), type="numeric", default = 1,
              help="Number of threads to use (Default is 1)", metavar = "character"),
  make_option(c("--Memory"), type="numeric", default = NA,
              help="Allocate memory for Seurat parallelization in GB (Default is NA)", metavar = "character"),
  make_option(c("--fraction"), type="numeric", default = 0.05,
              help="fraction of cells that a gene needs to be expressed in order to be included (Default is 0.05)", metavar = "character"),
  make_option(c("--KNN"), type="numeric", default = 25,
              help="Nearest-neighbors parameter (Default is 25)", metavar = "character"),
  make_option(c("--GroupBy"), type="character", default = "seurat_clusters,orig.ident",
              help="Specify the columns in meta.data to group by (comma-separated)", metavar = "character"),   
  make_option(c("--reduction"), type="character", default = "umap",
              help="Select the dimensionality reduction to perform KNN on (Default is umap)", metavar = "character"),
  make_option(c("--max_shared"), type="numeric", default = 10,
              help="Maximum number of shared cells between two metacells (Default is 10)", metavar = "character"),
  make_option(c("--IdentGroup"), type="character", default = "seurat_clusters",
              help="Set the Idents of the metacell seurat object (Default is umap)", metavar = "character"),
  make_option(c("--Sample"), type="character", default = NA,
              help="Sample's name", metavar = "character"),
  make_option(c("--outDir"), type="character", default = NA,
              help="Path to the output directory", metavar = "character")
);

#Get user information
opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", 
                                    option_list=option_list,
                                    epilogue="\n\nAuthor:Gil Sorek");
opt = optparse::parse_args(opt_parser);

# Set log framework
logfile = paste(opt$outDir,opt$Sample,'_log.txt',sep='')
cat(paste('[',Sys.time(),']: basic_analysis: ',opt$Sample,sep=''), file=logfile, sep='\n')
writeLog <- function(logfile, msg) { cat(paste('(',strftime(Sys.time(), format = "%H:%M:%S"),'): ',msg,sep=''),file=logfile,append=TRUE,sep='\n') }
saveRDS(opt, paste(opt$outDir,'opt.rds',sep=''))

# Parallelization
if (opt$CPUs>1) {
  writeLog(logfile, paste("Using multicore parallelization with ",opt$CPUs,' threads',sep=''))
  enableWGCNAThreads(nThreads = opt$CPUs)
  if (!is.na(opt$Memory)) {
    writeLog(logfile, paste("Allocating ",as.numeric(opt$Memory),'GB of memory',sep=''))
    options(future.globals.maxSize = as.numeric(opt$Memory) * 1024 ^ 3)
  }
}

print("Input var:",quote = F)
print.data.frame(as.data.frame(x = unlist(opt),row.names = names(opt)),right = F,quote = F)

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


obj_seurat <- SetupForWGCNA(
  obj_seurat,
  gene_select = "fraction", # the gene selection approach
  fraction = opt$fraction, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = opt$Sample # the name of the hdWGCNA experiment
)

obj_seurat <- MetacellsByGroups(
  seurat_obj = obj_seurat,
  group.by = unlist(stringi::stri_split(str = opt$GroupBy,fixed = ',')), # specify the columns in seurat_obj@meta.data to group by
  reduction = opt$reduction, # select the dimensionality reduction to perform KNN on
  k = opt$KNN, # nearest-neighbors parameter
  max_shared = opt$max_shared, # maximum number of shared cells between two metacells
  ident.group = opt$IdentGroup # set the Idents of the metacell seurat object
)


# # normalize metacell expression matrix:
obj_seurat   <- NormalizeMetacells(obj_seurat)
obj_seurat   <- ScaleMetacells(obj_seurat, features=VariableFeatures(obj_seurat))
metacell_obj <- GetMetacellObject(obj_seurat)

# Save results
saveRDS(metacell_obj, file = paste(opt$outDir, opt$Sample, '.rds', sep=''))
writeLog(logfile, paste("Finished"))