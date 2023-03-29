.libPaths()
library(optparse)
library(Seurat)
library(dplyr)
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--inputRDS"), type="character", default = NA,
              help="Vector of paths to Seurat object RDS files (comma-separated)", metavar = "character"),
  make_option(c("--Remove_lowQuality"), action="store_true", default = FALSE,
              help="Remove cells from low quality clusters (Default is False)", metavar = "character"),
  make_option(c("--Samples"), type="character", default = NA,
              help="Vector of Samples names (comma-separated)", metavar = "character"),
  make_option(c("--outDir"), type="character", default = NA,
              help="Path to the output directory", metavar = "character")
);

#Get user information
opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", 
                                    option_list=option_list,
                                    epilogue="\n\nAuthor:Gil Sorek");
opt = optparse::parse_args(opt_parser);

# Set log framework
logfile = paste(opt$outDir,'logs.txt',sep='')
cat(paste('[',Sys.time(),']: remove_lowQuality_clusters: ',opt$Samples,sep=''), file=logfile, sep='\n')
writeLog <- function(logfile, msg) { cat(paste('(',strftime(Sys.time(), format = "%H:%M:%S"),'): ',msg,sep=''),file=logfile,append=TRUE,sep='\n') }
saveRDS(opt, paste(opt$outDir,'opt.rds',sep=''))


print("Input var:",quote = F)
print.data.frame(as.data.frame(x = unlist(opt),row.names = names(opt)),right = F,quote = F)

# Import data
if (!is.na(opt$inputRDS)) {
  writeLog(logfile, paste("Importing Seurat objects file names...",sep=''))
  files = unlist(stringi::stri_split(str = opt$inputRDS,fixed = ','))
} else {
  writeLog(logfile, paste("ERROR: Seurat object RDS files must be specified [--inputRDS]"))
  stop()
}
if (is.na(opt$Samples)) {
  writeLog(logfile, paste("Samples names must be specified [--Samples]"))
  stop()
} else {
  samples = unlist(stringi::stri_split(str = opt$Samples,fixed = ','))
  if (length(files)!=length(samples)){
    writeLog(logfile, paste("ERROR: The number of sample's names is not equal to the number of sample's files"))
    stop()
  }
}
if (is.na(opt$outDir)) {
  writeLog(logfile, paste("Output directory must be specified [--outDir]"))
  stop()
}

summary_DF <- data.frame(row.names = c("num cells before",
                                       "max nFeature", "mean nFeature", "min nFeature", 
                                       "max nCounts","mean nCounts", "min nCounts", 
                                       "max percent.mt","mean percent.mt", "min percent.mt", 
                                       "lowExp clusters", "highExp clusters"))

for (sample in samples) {
  # Load Seurat object
  file = files[which(sample==samples)]
  writeLog(logfile, paste("New sample: ",sample,sep=''))
  writeLog(logfile, paste("inputRDS: ",file,sep=''))
  writeLog(logfile, paste("Importing Seurat RDS object...",sep=''))
  seurat_obj = readRDS(file)
  
  # Store statistics
  col2add <- list(dim(seurat_obj)[2],
                  max(seurat_obj@meta.data$nFeature_RNA), mean(seurat_obj@meta.data$nFeature_RNA), min(seurat_obj@meta.data$nFeature_RNA),
                  max(seurat_obj@meta.data$nCount_RNA), mean(seurat_obj@meta.data$nCount_RNA), min(seurat_obj@meta.data$nCount_RNA),
                  max(seurat_obj@meta.data$percent.mt), mean(seurat_obj@meta.data$percent.mt), min(seurat_obj@meta.data$percent.mt))
  temp_df <- data.frame(row.names = levels(seurat_obj@active.ident))
  
  
  # For loop to find the mean features and counts per cluster
  for(j in row.names(temp_df)) { 
    temp_df[j,"mean_counts"] <-  mean(seurat_obj@meta.data$nCount_RNA[seurat_obj@meta.data$seurat_clusters==j])
    temp_df[j,"mean_Features"] <-  mean(seurat_obj@meta.data$nFeature_RNA[seurat_obj@meta.data$seurat_clusters==j])
  }
  
  # Calculate threshold for nFeatures and nCounts
  features_thresh_high <- (mean(temp_df$mean_Features) + 2*sd(temp_df$mean_Features))
  counts_thresh_high <- (mean(temp_df$mean_counts) + 2*sd(temp_df$mean_counts))
  features_thresh_low <- (mean(temp_df$mean_Features) - 2*sd(temp_df$mean_Features))
  counts_thresh_low <- (mean(temp_df$mean_counts) - 2*sd(temp_df$mean_counts))
  writeLog(logfile, paste("Features threshold: ",features_thresh_high,' / ',features_thresh_low,sep=''))
  writeLog(logfile, paste("Counts threshold: ",counts_thresh_high,' / ',counts_thresh_low,sep=''))
  
  # Decide which clusters are unfit and remove them
  highExp <- list() 
  lowExp <- list()
  for(k in row.names(temp_df)){ 
    if(temp_df[k, "mean_Features"] < features_thresh_low & temp_df[k, "mean_counts"] < counts_thresh_low){
      lowExp <- c(lowExp, k) # list of clusters with nFeatures or nCounts lower than the threshold
    }
    if(temp_df[k, "mean_Features"] > features_thresh_high & temp_df[k, "mean_counts"] > counts_thresh_high){
      highExp <- c(highExp, k) # list of clusters with nFeatures or nCounts lower than the threshold
    }
  }
  
  # Add low/high expression meta-data feature
  seurat_obj@meta.data$lowExp = FALSE
  if (length(lowExp)>0)
    seurat_obj@meta.data$lowExp[seurat_obj@meta.data$seurat_clusters %in% unlist(lowExp)] = TRUE
  seurat_obj@meta.data$highExp = FALSE
  if (length(highExp)>0)
    seurat_obj@meta.data$highExp[seurat_obj@meta.data$seurat_clusters %in% unlist(highExp)] = TRUE
  
  # Plot low quality clusters
  writeLog(logfile, paste("Plotting figures..."))
  cells <- colnames(seurat_obj)
  if (any(seurat_obj$lowExp) & any(seurat_obj$highExp)) {
    highlights = list(cells[seurat_obj$lowExp==TRUE], cells[seurat_obj$highExp==TRUE])
    plt <- DimPlot(object = seurat_obj, label = T, pt.size = 1.5, label.size = 9,
                   cells.highlight = highlights, cols.highlight = c("lightblue","red"), cols = "grey") + 
      theme(legend.position = "none")
  } else if (any(seurat_obj$lowExp)) {
    highlights = list(cells[seurat_obj$lowExp==TRUE])
    plt <- DimPlot(object = seurat_obj, label = T, pt.size = 1.5, label.size = 9,
                   cells.highlight = highlights, cols.highlight = c("lightblue"), cols = "grey") + 
      theme(legend.position = "none")
  } else if (any(seurat_obj$highExp)) {
    highlights = list(cells[seurat_obj$highExp==TRUE])
    plt <- DimPlot(object = seurat_obj, label = T, pt.size = 1.5, label.size = 9,
                   cells.highlight = highlights, cols.highlight = c("red"), cols = "grey") + 
      theme(legend.position = "none")
  }
  if (exists('plt')) {
    jpeg(paste(opt$outDir,sample,'_lowQuality.jpeg',sep=''), width = 1600, height = 1200)
    print(plt)
    dev.off()
  }
  
  # Remove clusters with low expression
  if (isTRUE(opt$Remove_lowQuality)) {
    writeLog(logfile, paste("Removing lowQuality clusters..."))
    if (length(lowExp)>0) {
      writeLog(logfile, paste("Removed low expression clusters: ", paste(unlist(lowExp),collapse = ','), sep=''))
      seurat_obj <- subset(seurat_obj, idents = lowExp, invert = T)
    } else {
      writeLog(logfile, paste("None low expression clusters were found"))
    }
    if (length(highExp)>0) {
      writeLog(logfile, paste("Removed high expression clusters: ", paste(unlist(highExp),collapse = ','), sep=''))
      seurat_obj <- subset(seurat_obj, idents = highExp, invert = T)
    } else {
      writeLog(logfile, paste("None high expression clusters were found"))
    }
  }
  lowExp <- paste(lowExp, collapse = ', ')
  highExp <- paste(highExp, collapse=', ' ) # colapse to create one 'string' of cluster numbers
  
  # Check for #cells after removal and #cells with nFeatures>5K
  col2add <- c(col2add, lowExp, highExp)
  print.data.frame(temp_df)
  writeLog(logfile, paste("-----------------------------------------------------------------"))
  
  # Insert relevant data to summary table
  summary_DF[[sample]] <- col2add
  saveRDS(seurat_obj, file = paste(opt$outDir,sample,'.rds',sep=''))
  rm(j,k, temp_df, col2add, seurat_obj, features_thresh_low, counts_thresh_low, features_thresh_high, counts_thresh_high, plt, cells,highExp,lowExp)
}

summary_DF <- apply(summary_DF, 2, as.character)
row.names(summary_DF) <- c("num cells before", 
                           "max nFeature", "mean nFeature", "min nFeature", 
                           "max nCounts","mean nCounts", "min nCounts", 
                           "max percent.mt","mean percent.mt", "min percent.mt",
                           "lowExp clusters","highExp clusters")
write.csv(summary_DF, file = paste(opt$outDir, "summary_DF.csv", sep = ""), row.names = TRUE)
