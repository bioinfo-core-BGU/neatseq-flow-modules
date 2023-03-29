.libPaths()
library(optparse)
library(Seurat)
library(phateR)
library(dplyr)
library(ggplot2)
library(batchelor)
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--inputRDS"), type="character", default = NA,
              help="Path to Seurat object RDS file", metavar = "character"),
  make_option(c("-s", "--Sample"), type="character", default = NA,
              help="Sample's name", metavar = "character"),
  make_option(c("--min_cells"), type="numeric", default = 10,
              help="Keep genes expressed in at least <min_cells> cells (Default is 10)", metavar = "character"),
  make_option(c("--min_umi"), type="numeric", default = 1000,
              help="Keep cells with at least <min_umi> UMIs (Default is 1000)", metavar = "character"),
  make_option(c("--knn"), type="numeric", default = 5,
              help="Number of nearest neighbors on which to build kernel (Default is 5)", metavar = "character"),
  make_option(c("--decay"), type="numeric", default = 40,
              help="Sets decay rate of kernel tails (Default is 40)", metavar = "character"),
  make_option(c("--gamma"), type="numeric", default = 1,
              help="Informational distance constant between -1 and 1. ‘gamma=1’ gives the PHATE log potential, ‘gamma=0’ gives a square root potential (Default is 1)", metavar = "character"),
  make_option(c("--t"), type="character", default = 'auto',
              help="Power to which the diffusion operator is powered sets the level of diffusion (Default is auto)", metavar = "character"),
  make_option(c("--Markers"), type="character", default = NA,
              help="Vector of markers to plot with PHATE's embeddings (comma-separated)", metavar = "character"),
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
cat(paste('[',Sys.time(),']: PHATE: ',opt$Sample,sep=''), file=logfile, sep='\n')
writeLog <- function(logfile, msg) { cat(paste('(',strftime(Sys.time(), format = "%H:%M:%S"),'): ',msg,sep=''),file=logfile,append=TRUE,sep='\n') }

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
if (!is.na(opt$Markers)) {
  markers = unlist(stringi::stri_split(str = opt$Markers,fixed = ','))
}
if (opt$gamma >1 | opt$gamma < (-1)) {
  writeLog(logfile, paste("ERROR: PHATE's gamma parameter cannot be greater than 1 or lower than -1 [--gamma]"))
  stop()
}
if (opt$t != 'auto') {
  opt$t = as.numeric(opt$t)
}

# Load data
writeLog(logfile, paste("Loading data..."))
counts = as.matrix(GetAssayData(object = obj_seurat, slot = "counts"))
counts = t(counts)

# Filter Data
writeLog(logfile, paste("Filtering data..."))
writeLog(logfile, paste("Keeping genes expressed in at least ",opt$min_cells," cells",sep=''))
keep_cols <- colSums(counts > 0) > opt$min_cells
writeLog(logfile, paste("Keeping cells with at least ",opt$min_umi," UMIs",sep=''))
keep_rows <- rowSums(counts) > opt$min_umi
counts = counts[keep_rows,keep_cols]

# Normalize data
writeLog(logfile, paste("Normalizing data..."))
counts = library.size.normalize(counts)
counts = sqrt(counts)

# Meta-data
writeLog(logfile, paste("Writing cells meta-data..."))
meta = data.frame(row.names = rownames(counts))
meta$Cluster = obj_seurat$seurat_clusters[rownames(meta)]
if (length(unique(obj_seurat$orig.ident))>1) {
  meta$Sample = obj_seurat$orig.ident[rownames(meta)]
}

# Batch Correction
if ('Sample' %in% colnames(meta)) {
  writeLog(logfile, paste("Running batch correction..."))
  counts_corrected = fastMNN(t(counts), batch=meta$Sample)
  counts = t(assay(counts_corrected))
}

# PHATE
writeLog(logfile, paste("Running PHATE..."))
writeLog(logfile, paste("PHATE: Using knn=",opt$knn))
writeLog(logfile, paste("PHATE: Using decay=",opt$decay))
writeLog(logfile, paste("PHATE: Using gamma=",opt$gamma))
writeLog(logfile, paste("PHATE: Using t=",opt$t))
obj_phate = phate(counts, knn=opt$knn, decay=opt$decay, gamma=opt$gamma, t=opt$t)

# Plots
writeLog(logfile, paste("Plotting PHATE..."))
if ('Sample' %in% colnames(meta)) {
  plt = ggplot(obj_phate) + 
    geom_point(aes(PHATE1, PHATE2, color = meta$Cluster, shape=meta$Sample), size=2) +
	scale_color_brewer(palette="Paired") +
    labs(color='Cluster', shape='Sample')
} else {
  plt = ggplot(obj_phate) + 
    geom_point(aes(PHATE1, PHATE2, color = meta$Cluster), size=2) +
	scale_color_brewer(palette="Paired") +
    labs(color='Cluster')
}
ggsave(paste(opt$outDir,opt$Sample,"_PHATE.jpeg", sep = ""), plt, dpi=600, height=10, width=10)
if (!is.na(opt$Markers)) {
  for (marker in markers) {
    if (marker %in% colnames(counts)) {
	  plt = ggplot(obj_phate) + 
        geom_point(aes(PHATE1, PHATE2, color = counts[,marker]), size=2) +
	    viridis::scale_color_viridis(option="B") +
        labs(color=marker)
	  ggsave(paste(opt$outDir,opt$Sample,"_PHATE_",marker,".jpeg", sep = ""), plt, dpi=600, height=10, width=10)
	}
  }
}

writeLog(logfile, paste("Finished",sep=''))