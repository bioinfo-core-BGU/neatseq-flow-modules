.libPaths()
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(optparse)
library(Seurat)
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--inputRDS"), type="character", default = NA,
              help="Path to Seurat object RDS file", metavar = "character"),
  make_option(c("-s", "--Sample"), type="character", default = NA,
              help="Sample's name", metavar = "character"),
  make_option(c("-o", "--outDir"), type="character", default = NA,
              help="Path to the output directory", metavar = "character")
);

#Get user information
opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", 
                                    option_list=option_list,
                                    epilogue="\n\nAuthor:Gil Sorek");
opt = optparse::parse_args(opt_parser);

# Set log framework
logfile = paste(opt$outDir,opt$Sample,'_log.txt',sep='')
cat(paste('[',Sys.time(),']: RNA Velocyto: ',opt$Sample,sep=''), file=logfile, sep='\n')
writeLog <- function(logfile, msg) { cat(paste('(',strftime(Sys.time(), format = "%H:%M:%S"),'): ',msg,sep=''),file=logfile,append=TRUE,sep='\n') }
saveRDS(opt, paste(opt$outDir,'opt.rds',sep=''))

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

h5Seurat_file = paste(opt$outDir,opt$Sample,'.h5Seurat',sep='')
obj_seurat[["RNA"]] <- obj_seurat[["spliced"]]
DefaultAssay(obj_seurat) <- "RNA"
SaveH5Seurat(obj_seurat, filename = h5Seurat_file,overwrite = TRUE )
Convert(h5Seurat_file, dest = "h5ad",overwrite = TRUE)