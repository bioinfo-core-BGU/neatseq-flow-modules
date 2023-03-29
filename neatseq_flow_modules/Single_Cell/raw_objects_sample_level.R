.libPaths()
library(optparse)
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-r", "--CellRanger_h5"), type="character", default = NA,
              help="Path to CellRanger filtered feature bc matrix h5 file", metavar = "character"),
  make_option(c("-b", "--CellBender_h5"), type="character", default = NA,
              help="Path to CellBender filtered feature bc matrix h5 file", metavar = "character"),
  make_option(c("-s", "--Sample"), type="character", default = NA,
              help="Sample's name", metavar = "character"),
  make_option(c("-o", "--outDir"), type="character", default = NA,
              help="Path to the output directory", metavar = "character"),
  make_option(c("--min_cells"), type="numeric", default = 3,
              help="Include features detected in at least this many cells (Default is 3)", metavar = "character"),
  make_option(c("--UMI_low_cutoff"), type="numeric", default = 500,
              help="Exclude cells with nCount_RNA lower than cutoff (Default is 500)", metavar = "character"),
  make_option(c("--Genes_low_cutoff"), type="numeric", default = 300,
              help="Exclude cells with nFeature_RNA lower than cutoff (Default is 300)", metavar = "character"),
  make_option(c("--Genes_high_cutoff"), type="numeric", default = 5000,
              help="Exclude cells with nFeature_RNA higher than cutoff (Default is 5000)", metavar = "character"),
  make_option(c("--log10GenesPerUMI_cutoff"), type="numeric", default = 0.8,
              help="Exclude cells with log10(GenesPerUMI) lower than cutoff (Default is 0.8)", metavar = "character"),
  make_option(c("--MT_pattern"), type="character", default = NA,
              help="Mitochondrial genes pattern used to calculate percent.mt (usually 'MT' for human and 'mt' for mouse)", metavar = "character"),
  make_option(c("--MT_cutoff"), type="numeric", default = 0.05,
              help="Mitochondrial genes cutoff (Default is 0.05)", metavar = "character")
);

#Get user information
opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", 
                                    option_list=option_list,
                                    epilogue="\n\nAuthor:Gil Sorek");
opt = optparse::parse_args(opt_parser);

# Set log framework
logfile = paste(opt$outDir,opt$Sample,'_log.txt',sep='')
cat(paste('[',Sys.time(),']: raw_objects: ',opt$Sample,sep=''), file=logfile, sep='\n')
writeLog <- function(logfile, msg) { cat(paste('(',strftime(Sys.time(), format = "%H:%M:%S"),'): ',msg,sep=''),file=logfile,append=TRUE,sep='\n') }
saveRDS(opt, paste(opt$outDir,'opt.rds',sep=''))

print("Input var:",quote = F)
print.data.frame(as.data.frame(x = unlist(opt),row.names = names(opt)),right = F,quote = F)

# Functions
QC_Metrics <- function(metadata, dir) {
    plt_NCells = metadata %>% ggplot(aes(x=orig.ident, fill=orig.ident)) +
	geom_bar() +
	theme_classic() +
	theme(plot.title = element_text(hjust=0.5, face="bold")) +
	ggtitle("NCells") +
	xlab("")
  ggsave(paste(opt$outDir,dir,'/',opt$Sample,'_NCells.jpeg',sep=''), plt_NCells,dpi=600,width=6,height=8)
  
  plt_UMIperCell = metadata %>% ggplot(aes(x=nCount_RNA, fill=orig.ident)) +
	geom_density(alpha=0.6) +
	scale_x_log10() +
	theme_classic() +
	ylab("Cell density") +
	geom_vline(xintercept = opt$UMI_low_cutoff) +
	geom_label(aes(opt$UMI_low_cutoff,+Inf,vjust=2), label=paste(opt$UMI_low_cutoff), show.legend = FALSE)
  ggsave(paste(opt$outDir,dir,'/',opt$Sample,'_UMI_Counts_per_Cell.jpeg',sep=''), plt_UMIperCell,dpi=600,width=6,height=8)
  
  plt_GenesperCell = metadata %>% ggplot(aes(x=nFeature_RNA, fill=orig.ident)) +
	geom_density(alpha=0.6) +
	scale_x_log10() +
	theme_classic() +
	ylab("Cell density") +
	geom_vline(xintercept = c(opt$Genes_low_cutoff,opt$Genes_high_cutoff)) +
	geom_label(aes(opt$Genes_low_cutoff,+Inf,vjust=2), label=paste(opt$Genes_low_cutoff), show.legend = FALSE) +
	geom_label(aes(opt$Genes_high_cutoff,+Inf,vjust=2), label=paste(opt$Genes_high_cutoff), show.legend = FALSE)
  ggsave(paste(opt$outDir,dir,'/',opt$Sample,'_Genes_per_Cell.jpeg',sep=''), plt_GenesperCell,dpi=600,width=6,height=8)
  
  plt_UMIvsGenes = metadata %>% ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
	geom_point() +
	scale_colour_gradient(low = "gray90", high = "black") +
	stat_smooth(method=lm) +
	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
	geom_vline(xintercept = opt$UMI_low_cutoff) +
	geom_hline(yintercept = c(opt$Genes_low_cutoff,opt$Genes_high_cutoff)) +
	geom_label(aes(opt$UMI_low_cutoff,+Inf,vjust=2), label=paste(opt$UMI_low_cutoff), show.legend = FALSE) +
	geom_label(aes(+Inf,opt$Genes_low_cutoff,hjust=1,vjust=2), label=paste(opt$Genes_low_cutoff), show.legend = FALSE) +
	geom_label(aes(+Inf,opt$Genes_high_cutoff,hjust=1,vjust=2), label=paste(opt$Genes_high_cutoff), show.legend = FALSE)
  ggsave(paste(opt$outDir,dir,'/',opt$Sample,'_UMI_vs_Genes.jpeg',sep=''), plt_UMIvsGenes,dpi=600,width=6,height=8)
  
  plt_MitoCounts = metadata %>% ggplot(aes(x=percent.mt, fill=orig.ident)) +
	geom_density(alpha=0.6) +
	scale_x_log10() + 
  	theme_classic() +
	ylab("Cell density") +
	geom_vline(xintercept = opt$MT_cutoff) +
	geom_label(aes(opt$MT_cutoff,+Inf,vjust=2), label=paste(opt$MT_cutoff), show.legend = FALSE)
  ggsave(paste(opt$outDir,dir,'/',opt$Sample,'_Mito_Ratio.jpeg',sep=''), plt_MitoCounts,dpi=600,width=6,height=8)
  
  plt_Complexity = metadata %>% ggplot(aes(x=log10GenesPerUMI, fill=orig.ident)) +
	geom_density(alpha=0.6) +
	theme_classic() +
	ylab("Cell density") +
	geom_vline(xintercept = opt$log10GenesPerUMI_cutoff) +
	geom_label(aes(opt$log10GenesPerUMI_cutoff,+Inf,vjust=2), label=paste(opt$log10GenesPerUMI_cutoff), show.legend = FALSE)
  ggsave(paste(opt$outDir,dir,'/',opt$Sample,'_log10_Genes_per_UMI.jpeg',sep=''), plt_Complexity,dpi=600,width=6,height=8)
}

# Import data
if (!is.na(opt$CellRanger_h5)) {
  ranger_data <- Read10X_h5(filename = opt$CellRanger_h5, 
                            use.names = TRUE, unique.features = TRUE)
  writeLog(logfile, paste("Import CellRanger: ",dim(ranger_data)[2],' cells',sep=''))
} else {
  writeLog(logfile, paste("ERROR: CellRanger filtered h5 file must be specified [--CellRanger_h5]"))
  stop()
}
if (!is.na(opt$CellBender_h5)) {
  bender_data <- Read10X_h5(filename = opt$CellBender_h5, 
                            use.names = TRUE, unique.features = TRUE)
  writeLog(logfile, paste("Import CellBender: ",dim(bender_data)[2],' cells',sep=''))
} else {
  writeLog(logfile, paste("ERROR: CellBender filtered h5 file must be specified [--CellBender_h5]"))
  stop()
}
if (!is.na(opt$MT_pattern)) {
  MT_pattern = paste("^",opt$MT_pattern,"-",sep='')
  writeLog(logfile, paste("Using pattern '",MT_pattern,"' to detect mitochondrial genes",sep=''))
} else {
  writeLog(logfile, paste("ERROR: Mitochondrial genes pattern must be specified [--MT_pattern]"))
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

# Intersect CellRanger and CellBender outputs
cells_inBoth <- intersect(colnames(ranger_data), colnames(bender_data))
writeLog(logfile, paste("Number of cells in both datasets: ",length(cells_inBoth),sep=''))
writeLog(logfile, paste('Unique cells in CellRanger: ',sum(!(colnames(ranger_data) %in% cells_inBoth)),sep = ''))
writeLog(logfile, paste('Unique cells in CellBender: ',sum(!(colnames(bender_data) %in% cells_inBoth)),sep = ''))
writeLog(logfile, paste("Creating merged Seurat object..."))
writeLog(logfile, paste("Removing features detected only in ",opt$min_cells," or less cells",sep=''))
merged_data = bender_data[,cells_inBoth]
obj_seurat = CreateSeuratObject(counts = merged_data, project = opt$Sample, min.cells = opt$min_cells)
obj_seurat$log10GenesPerUMI = log10(obj_seurat$nFeature_RNA) / log10(obj_seurat$nCount_RNA)
obj_seurat[["percent.mt"]] <- PercentageFeatureSet(obj_seurat, pattern = MT_pattern)
obj_seurat@meta.data[which(is.na(obj_seurat$percent.mt)),'percent.mt'] = 0
obj_seurat[["percent.mt"]] <- obj_seurat[["percent.mt"]] / 100

# Pre-filtering QC Plots
metadata <- obj_seurat@meta.data
if (!dir.exists(paste(opt$outDir,'Raw_Metrics',sep='')))
  dir.create(paste(opt$outDir,'Raw_Metrics',sep=''))
writeLog(logfile, paste("Plotting pre-filtering sample metrics..."))
QC_Metrics(metadata,'Raw_Metrics')

# Filter Object
obj_filtered = subset(obj_seurat, subset = (nCount_RNA >= opt$UMI_low_cutoff) & (nFeature_RNA >= opt$Genes_low_cutoff) & (nFeature_RNA <= opt$Genes_high_cutoff) & (log10GenesPerUMI > opt$log10GenesPerUMI_cutoff) & (percent.mt < opt$MT_cutoff))
writeLog(logfile, paste("Filtered data has ",dim(obj_seurat)[2]," barcodes and ",dim(obj_seurat)[1]," features",sep=''))
writeLog(logfile, paste("Median number of genes detected per cell: ",median(obj_seurat$nFeature_RNA),sep=''))
writeLog(logfile, paste("Median number of UMIs detected per cell: ",median(obj_seurat$nCount_RNA),sep=''))

# Filtered Object QC Plots
metadata <- obj_filtered@meta.data
if (!dir.exists(paste(opt$outDir,'Filtered',sep='')))
  dir.create(paste(opt$outDir,'Filtered',sep=''))
writeLog(logfile, paste("Plotting filtered sample metrics..."))
QC_Metrics(metadata,'Filtered')

# Save results
saveRDS(obj_filtered, file = paste(opt$outDir, opt$Sample, '.rds', sep=''))
writeLog(logfile, paste("Finished"))