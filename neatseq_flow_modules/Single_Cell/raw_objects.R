.libPaths()
library(optparse)
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-r", "--CellRanger_h5"), type="character", default = NA,
              help="Vector of paths to CellRanger filtered feature bc matrix h5 files (comma-separated)", metavar = "character"),
  make_option(c("-b", "--CellBender_h5"), type="character", default = NA,
              help="Vector of paths to CellBender filtered feature bc matrix h5 files (comma-separated)", metavar = "character"),
  make_option(c("-l", "--Loom_Files"), type="character", default = NA,
              help="Vector of paths to Loom files (comma-separated) with un-spliced information", metavar = "character"),
  make_option(c("--Samples"), type="character", default = NA,
              help="Vector of Samples names (comma-separated)", metavar = "character"),
  make_option(c("--outDir"), type="character", default = NA,
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
  make_option(c("--Filtered_plots_only"), action="store_true", default = FALSE,
              help="Plot filtered QC plots only (Default is False)", metavar = "character"),
  make_option(c("--Use_Spiced_as_RNA"), action="store_true", default = FALSE,
              help="Use the Sliced assay (from loom file) as RNA (Default is to use the CellBender counts)", metavar = "character"),
  make_option(c("--MT_pattern"), type="character", default = NA,
              help="Mitochondrial genes pattern used to calculate percent.mt (usually 'MT' for human and 'mt' for mouse)", metavar = "character"),
  make_option(c("--MT_cutoff"), type="numeric", default = 0.05,
              help="Exclude cells with Mitochondrial genes ratio higher than cutoff (Default is 0.05)", metavar = "character")
);

#Get user information
opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", 
                                    option_list=option_list,
                                    epilogue="\n\nAuthor:Gil Sorek");
opt = optparse::parse_args(opt_parser);

# Set log framework
logfile = paste(opt$outDir,'Raw_objects_log.txt',sep='')
cat(paste('[',Sys.time(),']: raw_objects: ',opt$Samples,sep=''), file=logfile, sep='\n')
writeLog <- function(logfile, msg) { cat(paste('(',strftime(Sys.time(), format = "%H:%M:%S"),'): ',msg,sep=''),file=logfile,append=TRUE,sep='\n') }
saveRDS(opt, paste(opt$outDir,'opt.rds',sep=''))

print("Input var:",quote = F)
print.data.frame(as.data.frame(x = unlist(opt),row.names = names(opt)),right = F,quote = F)

# Functions
QC_Metrics <- function(metadata, dir) {
  plt_NCells = metadata %>% ggplot(aes(x=orig.ident, fill=orig.ident)) +
	geom_bar(stat="count") +
	geom_text(aes(label = ..count..), stat = 'count', size=6, vjust=1.3) +
	theme_classic() +
	theme(plot.title = element_text(hjust=0.5, face="bold"), panel.spacing = unit(4, "lines")) +
	ggtitle("Number of Cells") +
	xlab("") +
	facet_wrap(vars(Source))
  ggsave(paste(opt$outDir,dir,'/','NCells.jpeg',sep=''), plt_NCells,dpi=600,width=10,height=8)
  
  plt_UMIperCell = metadata %>% ggplot(aes(x=nCount_RNA, fill=orig.ident)) +
	geom_density(alpha=0.2) +
	scale_x_log10() +
	theme_classic() +
	theme(plot.title = element_text(hjust=0.5, face="bold"), panel.spacing = unit(4, "lines")) +
	ggtitle("Number of UMIs per Cell") +
	xlab("Number of UMIs") +
	ylab("Cell density") +
	geom_vline(xintercept = opt$UMI_low_cutoff) +
	geom_label(aes(opt$UMI_low_cutoff,+Inf,vjust=2), label=paste(opt$UMI_low_cutoff), show.legend = FALSE, color="black", fill="red") +
	facet_wrap(vars(Source), scales="free")
  ggsave(paste(opt$outDir,dir,'/','UMI_Counts_per_Cell.jpeg',sep=''), plt_UMIperCell,dpi=600,width=12,height=8)
  
  plt_GenesperCell = metadata %>% ggplot(aes(x=nFeature_RNA, fill=orig.ident)) +
	geom_density(alpha=0.2) +
	scale_x_log10() +
	theme_classic() +
	theme(plot.title = element_text(hjust=0.5, face="bold"), panel.spacing = unit(4, "lines")) +
	ggtitle("Number of Genes per Cell") +
	xlab("Number of Genes") +
	ylab("Cell density") +
	geom_vline(xintercept = c(opt$Genes_low_cutoff,opt$Genes_high_cutoff)) +
	geom_label(aes(opt$Genes_low_cutoff,+Inf,vjust=2), label=paste(opt$Genes_low_cutoff), show.legend = FALSE, color="black", fill="red") +
	geom_label(aes(opt$Genes_high_cutoff,+Inf,vjust=2), label=paste(opt$Genes_high_cutoff), show.legend = FALSE, color="black", fill="red") +
	facet_wrap(vars(Source), scales="free")
  ggsave(paste(opt$outDir,dir,'/','Genes_per_Cell.jpeg',sep=''), plt_GenesperCell,dpi=600,width=12,height=8)
  
  plt_UMIvsGenes = metadata %>% filter(Source == 'Raw') %>% ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
	geom_point() +
	scale_colour_gradient(low = "gray90", high = "black") +
	stat_smooth(method=lm) +
	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
	theme(plot.title = element_text(hjust=0.5, face="bold")) +
	xlab("Number of UMIs") +
	ylab("Number of Genes") +
	ggtitle("UMIs vs. Genes: Raw data") +
	geom_vline(xintercept = opt$UMI_low_cutoff) +
	geom_hline(yintercept = c(opt$Genes_low_cutoff,opt$Genes_high_cutoff)) +
	geom_label(aes(opt$UMI_low_cutoff,+Inf,vjust=2), label=paste(opt$UMI_low_cutoff), show.legend = FALSE, color="black", fill="red") +
	geom_label(aes(0,opt$Genes_low_cutoff,hjust=-2), label=paste(opt$Genes_low_cutoff), show.legend = FALSE, color="black", fill="red") +
	geom_label(aes(0,opt$Genes_high_cutoff,hjust=-2), label=paste(opt$Genes_high_cutoff), show.legend = FALSE, color="black", fill="red") +
	facet_wrap(~orig.ident, scales="free")
  ggsave(paste(opt$outDir,dir,'/','UMI_vs_Genes_Raw.jpeg',sep=''), plt_UMIvsGenes,dpi=600,width=10,height=8)
  
  plt_UMIvsGenes = metadata %>% filter(Source == 'Filtered') %>% ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
	geom_point() +
	scale_colour_gradient(low = "gray90", high = "black") +
	stat_smooth(method=lm) +
	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
	theme(plot.title = element_text(hjust=0.5, face="bold")) +
	xlab("Number of UMIs") +
	ylab("Number of Genes") +
	ggtitle("UMIs vs. Genes: Filtered data") +
	geom_vline(xintercept = opt$UMI_low_cutoff) +
	geom_hline(yintercept = c(opt$Genes_low_cutoff,opt$Genes_high_cutoff)) +
	geom_label(aes(opt$UMI_low_cutoff,+Inf,vjust=2), label=paste(opt$UMI_low_cutoff), show.legend = FALSE, color="black", fill="red") +
	geom_label(aes(0,opt$Genes_low_cutoff,hjust=-2), label=paste(opt$Genes_low_cutoff), show.legend = FALSE, color="black", fill="red") +
	geom_label(aes(0,opt$Genes_high_cutoff,hjust=-2), label=paste(opt$Genes_high_cutoff), show.legend = FALSE, color="black", fill="red") +
	facet_wrap(~orig.ident, scales="free")
  ggsave(paste(opt$outDir,dir,'/','UMI_vs_Genes_Filtered.jpeg',sep=''), plt_UMIvsGenes,dpi=600,width=10,height=8)
  
  if (var(metadata$percent.mt)>0){
    plt_MitoCounts = metadata %>% ggplot(aes(x=percent.mt, fill=orig.ident)) +
	geom_density(alpha=0.2) +
	scale_x_log10() + 
  	theme_classic() +
	theme(plot.title = element_text(hjust=0.5, face="bold"), panel.spacing = unit(4, "lines")) +
	xlab("Mitochondrial Genes Ratio") +
	ggtitle("Mitochondrial Genes Ratio") +
	ylab("Cell density") +
	geom_vline(xintercept = opt$MT_cutoff) +
	geom_label(aes(opt$MT_cutoff,+Inf,vjust=2), label=paste(opt$MT_cutoff), show.legend = FALSE, color="black", fill="red") +
	facet_wrap(vars(Source), scales="free")
    ggsave(paste(opt$outDir,dir,'/','Mito_Ratio.jpeg',sep=''), plt_MitoCounts,dpi=600,width=12,height=8)
  }
  plt_Complexity = metadata %>% ggplot(aes(x=log10GenesPerUMI, fill=orig.ident)) +
	geom_density(alpha=0.2) +
	theme_classic() +
	theme(plot.title = element_text(hjust=0.5, face="bold"), panel.spacing = unit(4, "lines")) +
	ggtitle("Genes per UMI Ratio") +
	ylab("Cell density") +
	geom_vline(xintercept = opt$log10GenesPerUMI_cutoff) +
	geom_label(aes(opt$log10GenesPerUMI_cutoff,+Inf,vjust=2), label=paste(opt$log10GenesPerUMI_cutoff), show.legend = FALSE, color="black", fill="red") +
	facet_wrap(vars(Source), scales="free")
  ggsave(paste(opt$outDir,dir,'/','log10_Genes_per_UMI.jpeg',sep=''), plt_Complexity,dpi=600,width=12,height=8)
}

# Import data
if (!is.na(opt$CellRanger_h5)) {
  writeLog(logfile, paste("Importing CellRanger objects file names...",sep=''))
  cr_files = unlist(stringi::stri_split(str = opt$CellRanger_h5,fixed = ','))
} else {
  writeLog(logfile, paste("ERROR: CellRanger filtered h5 files must be specified [--CellRanger_h5]"))
  stop()
}
if (!is.na(opt$CellBender_h5)) {
  writeLog(logfile, paste("Importing CellBender objects file names...",sep=''))
  cb_files = unlist(stringi::stri_split(str = opt$CellBender_h5,fixed = ','))
  if (length(cb_files)!=length(cr_files)) {
    writeLog(logfile, paste("ERROR: The number of CellRanger files is not equal to the number of CellBender files"))
	stop()
  }
} else {
  writeLog(logfile, paste("ERROR: CellBender filtered h5 files must be specified [--CellBender_h5]"))
  stop()
}


if (!is.na(opt$Loom_Files)) {
  writeLog(logfile, paste("Importing loom objects file names...",sep=''))
  loom_files = unlist(stringi::stri_split(str = opt$Loom_Files,fixed = ','))
} else {
  loom_files=c()
}

if (!is.na(opt$MT_pattern)) {
  MT_pattern = paste("^",opt$MT_pattern,"-",sep='')
  writeLog(logfile, paste("Using pattern '",MT_pattern,"' to detect mitochondrial genes",sep=''))
} else {
  writeLog(logfile, paste("ERROR: Mitochondrial genes pattern must be specified [--MT_pattern]"))
  stop()
}
if (is.na(opt$Samples)) {
  writeLog(logfile, paste("ERROR: Samples name must be specified [--Sample]"))
  stop()
} else {
  samples = unlist(stringi::stri_split(str = opt$Samples,fixed = ','))
  if (length(cr_files)!=length(samples)){
    writeLog(logfile, paste("ERROR: The number of sample's names is not equal to the number of CellRanger/CellBender files"))
    stop()
  }
}
if (is.na(opt$outDir)) {
  writeLog(logfile, paste("ERROR: Output directory must be specified [--outDir]"))
  stop()
}

obj_list = list()
raw_metadata = data.frame()
writeLog(logfile, paste("Loading Datasets..."))
# Load datasets and intersect
for (sample in samples) {
  writeLog(logfile, paste("Loading Sample: ",sample,sep=''))
  cr_file = cr_files[which(sample==samples)]
  cb_file = cb_files[which(sample==samples)]
  ranger_data <- Read10X_h5(filename = cr_file, 
                            use.names = TRUE, unique.features = TRUE)
  writeLog(logfile, paste("Import CellRanger: ",dim(ranger_data)[2],' cells',sep=''))
  bender_data <- Read10X_h5(filename = cb_file, 
                            use.names = TRUE, unique.features = TRUE)
  writeLog(logfile, paste("Import CellBender: ",dim(bender_data)[2],' cells',sep=''))
  cells_inBoth <- intersect(colnames(ranger_data), colnames(bender_data))
  writeLog(logfile, paste("Number of cells in both datasets: ",length(cells_inBoth),sep=''))
  writeLog(logfile, paste('Unique cells in CellRanger: ',sum(!(colnames(ranger_data) %in% cells_inBoth)),sep = ''))
  writeLog(logfile, paste('Unique cells in CellBender: ',sum(!(colnames(bender_data) %in% cells_inBoth)),sep = ''))
  writeLog(logfile, paste("Creating merged Seurat object..."))
  writeLog(logfile, paste("Removing features detected only in ",opt$min_cells," or less cells",sep=''))
  
  if (length(loom_files)>0){
      if (opt$Use_Spiced_as_RNA){
          cells_inBoth = unlist(lapply(X= cells_inBoth, FUN= function(x) stringi::stri_split_fixed(str=x,pattern="-")[[1]][1]))
      }else{
          merged_data     = bender_data[,cells_inBoth]
          colnames(merged_data) = unlist(lapply(X= colnames(merged_data), FUN= function(x) stringi::stri_split_fixed(str=x,pattern="-")[[1]][1]))
          obj_seurat      = CreateSeuratObject(counts = merged_data, project = sample, min.cells = opt$min_cells)
          cells_inBoth    = names(obj_seurat$orig.ident)
      }
      library(Seurat)
      library(SeuratDisk)
      library(SeuratWrappers)
      
      lo_file = loom_files[which(sample==samples)]
      ldat <- ReadVelocity(file = lo_file)
      obj_seurat_loom <- as.Seurat(x = ldat, project = sample, min.cells = 0)
      obj_seurat_loom$orig.ident = sample
      
      New_cell_names= stringi::stri_replace_first_fixed(str=names(obj_seurat_loom$orig.ident),pattern=paste(sample,":",sep=""),replacement="")
      New_cell_names= stringi::stri_replace_last_fixed(str=New_cell_names,pattern="x",replacement="")
      
      obj_seurat_loom = RenameCells(obj_seurat_loom, new.names =New_cell_names)
      obj_seurat_loom = subset(obj_seurat_loom,cells = cells_inBoth)
      
      if (opt$Use_Spiced_as_RNA){
          obj_seurat <- obj_seurat_loom
          obj_seurat[["RNA"]] <- obj_seurat[["spliced"]]
      }else{
          obj_seurat[["spliced"]]   <- obj_seurat_loom[["spliced"]]
          obj_seurat[["unspliced"]] <- obj_seurat_loom[["unspliced"]]
          obj_seurat[["ambiguous"]] <- obj_seurat_loom[["ambiguous"]]
      }
  }else{
      merged_data = bender_data[,cells_inBoth]
      obj_seurat = CreateSeuratObject(counts = merged_data, project = sample, min.cells = opt$min_cells)
  }
  obj_seurat$log10GenesPerUMI = log10(obj_seurat$nFeature_RNA) / log10(obj_seurat$nCount_RNA)
  obj_seurat[["percent.mt"]] <- PercentageFeatureSet(obj_seurat, pattern = MT_pattern)
  obj_seurat@meta.data[which(is.na(obj_seurat$percent.mt)),'percent.mt'] = 0
  obj_seurat[["percent.mt"]] <- obj_seurat[["percent.mt"]] / 100
  raw_metadata = rbind(raw_metadata, obj_seurat@meta.data)
  obj_list[[sample]] = obj_seurat
  rm(cr_file,cb_file,ranger_data,bender_data,cells_inBoth,merged_data,obj_seurat)
}
raw_metadata$Source = 'Raw'

# Filter Object
filtered_list = list()
filtered_metadata = data.frame()
writeLog(logfile, paste("Filtering Objects..."))
for (sample in samples) {
  writeLog(logfile, paste("Filtering Sample: ",sample,sep=''))
  sample_data = FetchData(obj_list[[sample]], vars = c('nCount_RNA','nFeature_RNA','log10GenesPerUMI','percent.mt'))
  writeLog(logfile, paste("There are ",sum(sample_data$nCount_RNA==0)," cells with zero counts",sep=''))
  sample_cutoffs = data.frame(row.names = rownames(sample_data))
  sample_cutoffs$Sample = sample
  sample_cutoffs[[paste('nCount_RNA>=',opt$UMI_low_cutoff,sep='')]] = sample_data$nCount_RNA >= opt$UMI_low_cutoff
  sample_cutoffs[[paste('nFeature_RNA>=',opt$Genes_low_cutoff,sep='')]] = sample_data$nFeature_RNA >= opt$Genes_low_cutoff
  sample_cutoffs[[paste('nFeature_RNA<=',opt$Genes_high_cutoff,sep='')]] = sample_data$nFeature_RNA <= opt$Genes_high_cutoff
  sample_cutoffs[[paste('log10GenesPerUMI>',opt$log10GenesPerUMI_cutoff,sep='')]] = sample_data$log10GenesPerUMI > opt$log10GenesPerUMI_cutoff
  sample_cutoffs[[paste('percent.mt<',opt$MT_cutoff,sep='')]] = sample_data$percent.mt < opt$MT_cutoff
  sample_cutoffs[which(is.na(sample_cutoffs[[paste('log10GenesPerUMI>',opt$log10GenesPerUMI_cutoff,sep='')]])),paste('log10GenesPerUMI>',opt$log10GenesPerUMI_cutoff,sep='')] = FALSE
  writeLog(logfile, paste("There are ",sum(!(sample_cutoffs[[paste('nCount_RNA>=',opt$UMI_low_cutoff,sep='')]]))," cells with ",paste('nCount_RNA<',opt$UMI_low_cutoff,sep=''),sep=''))
  writeLog(logfile, paste("There are ",sum(!(sample_cutoffs[[paste('nFeature_RNA>=',opt$Genes_low_cutoff,sep='')]]))," cells with ",paste('nFeature_RNA<',opt$Genes_low_cutoff,sep=''),sep=''))
  writeLog(logfile, paste("There are ",sum(!(sample_cutoffs[[paste('nFeature_RNA<=',opt$Genes_high_cutoff,sep='')]]))," cells with ",paste('nFeature_RNA>',opt$Genes_high_cutoff,sep=''),sep=''))
  writeLog(logfile, paste("There are ",sum(!(sample_cutoffs[[paste('log10GenesPerUMI>',opt$log10GenesPerUMI_cutoff,sep='')]]))," cells with ",paste('log10GenesPerUMI<=',opt$log10GenesPerUMI_cutoff,sep=''),sep=''))
  writeLog(logfile, paste("There are ",sum(!(sample_cutoffs[[paste('percent.mt<',opt$MT_cutoff,sep='')]]))," cells with ",paste('percent.mt>',opt$MT_cutoff,sep=''),sep=''))
  write.csv(sample_cutoffs, paste(opt$outDir,opt$Sample,"_cutoffs.csv",sep=''))
  obj_filtered = subset(obj_list[[sample]], subset = (nCount_RNA >= opt$UMI_low_cutoff) & (nFeature_RNA >= opt$Genes_low_cutoff) & (nFeature_RNA <= opt$Genes_high_cutoff) & (log10GenesPerUMI > opt$log10GenesPerUMI_cutoff) & (percent.mt < opt$MT_cutoff))
  writeLog(logfile, paste("Filtered data has ",dim(obj_filtered)[2]," barcodes and ",dim(obj_filtered)[1]," features",sep=''))
  writeLog(logfile, paste("Median number of genes detected per cell: ",median(obj_filtered$nFeature_RNA),sep=''))
  writeLog(logfile, paste("Median number of UMIs detected per cell: ",median(obj_filtered$nCount_RNA),sep=''))
  filtered_metadata = rbind(filtered_metadata, obj_filtered@meta.data)
  filtered_list[[sample]] = obj_filtered
  rm(obj_filtered,sample_cutoffs,sample_data)
}
filtered_metadata$Source = 'Filtered'

# QC Plots
metadata = rbind(raw_metadata, filtered_metadata)
metadata$Source = factor(metadata$Source,levels=c('Raw','Filtered'))
if (!dir.exists(paste(opt$outDir,'QC_Plots',sep='')))
  dir.create(paste(opt$outDir,'QC_Plots',sep=''))
writeLog(logfile, paste("Plotting filtered sample metrics..."))
QC_Metrics(metadata,'QC_Plots')

# Save results
writeLog(logfile, paste("Saving Results..."))
for (sample in samples) {
  saveRDS(filtered_list[[sample]], file = paste(opt$outDir, sample, '.rds', sep=''))
  # saveRDS(obj_list[[sample]],      file = paste(opt$outDir, sample, '.original.rds', sep=''))
}
  saveRDS(metadata,      file = paste(opt$outDir, 'MetaData.rds', sep=''))
writeLog(logfile, paste("Finished"))