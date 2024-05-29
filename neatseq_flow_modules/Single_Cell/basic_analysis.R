.libPaths()
library(optparse)
library(Seurat)
library(dplyr)
library(ggplot2)
library(future)
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--inputRDS"), type="character", default = NA,
              help="Path to Seurat object RDS file", metavar = "character"),
  make_option(c("--VariableFeatures"), type="numeric", default = 2000,
              help="Number of features to select as top variable features (Default is 2000)", metavar = "character"),
  make_option(c("--JackStraw_dims"), type="numeric", default = 50,
              help="Number of PCs to compute significance for in JackStraw (Default is 50)", metavar = "character"),
  make_option(c("--Remove_Clusters"), type="character", default = NA,
              help="Vector of clusters to remove (comma-separated)", metavar = "character"),
  make_option(c("--Exclude_Outliers"), action="store_true", default = FALSE,
              help="Remove extreme low/high expression clusters", metavar = "character"),
  make_option(c("--lowExp_col"), type="character", default = 'lowExp',
              help="Field in meta-data for low expression indicator (Default is lowExp)", metavar = "character"),
  make_option(c("--highExp_col"), type="character", default = 'highExp',
              help="Field in meta-data for high expression indicator (Default is highExp)", metavar = "character"),
  make_option(c("--CPUs"), type="numeric", default = 1,
              help="Number of threads to use (Default is 1)", metavar = "character"),
  make_option(c("--Memory"), type="numeric", default = NA,
              help="Allocate memory for Seurat parallelization in GB (Default is NA)", metavar = "character"),
  make_option(c("--SCTransform"), action="store_true", default = FALSE,
              help="Will use SCTransform in Seurat instead of standard log-normalization", metavar = "character"),
  make_option(c("--UseMagic"), action="store_true", default = FALSE,
              help="Will use Magic to Impute the data after normalization", metavar = "character"),
  make_option(c("--MagicCondaEnv"), type="character", default = NA,
              help="If --UseMagic is set it will this Conda env to fined the package [Must Be a Full Path!]", metavar = "character"),
  make_option(c("--SCT_ncells"), type="numeric", default = 5000,
              help="Number of subsampling cells used to build NB regression (Default is 5000)", metavar = "character"),
  make_option(c("--SCT_variable_features"), type="numeric", default = 3000,
              help="Use this many features as variable features after ranking by residual variance (Default is 3000)", metavar = "character"),
  make_option(c("--SCT_regress_vars"), type="character", default = NA,
              help="Variables to regress out in a second non-regularized linear regression (comma-separated)", metavar = "character"),
  make_option(c("--CellCycleScoring"), action="store_true", default = FALSE,
              help=" Perform CellCycleScoring to identify the cell-cycle stat of the cells (Default is False)", metavar = "character"),
  make_option(c("--CellCycle2regressOut"), action="store_true", default = FALSE,
              help=" Regress Out CellCycle stat of the cells (Default is False)", metavar = "character"),
  make_option(c("--CellCycle2regressOut_CC_Difference"), action="store_true", default = FALSE,
              help="Using the CC.Difference regression method (Default is FALSE)", metavar = "character"),
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
  plan("multicore", workers = opt$CPUs)
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
if (!is.na(opt$Remove_Clusters)) {
  writeLog(logfile, paste("Removing unwanted clusters..."))
  clusters2remove <- as.numeric(unlist(stringi::stri_split(str = opt$Remove_Clusters,fixed = ',')))
  if (isFALSE(all(clusters2remove %in% obj_seurat$seurat_clusters))) {
    writeLog(logfile, paste("ERROR: Cannot find the following identities in the object: ",opt$Remove_Clusters,sep=''))
    stop()
  }
  obj_seurat = subset(obj_seurat, idents = clusters2remove, invert = TRUE)
  writeLog(logfile, paste("Successfully removed the following identities from the object: ",opt$Remove_Clusters,sep=''))
}
if (opt$Exclude_Outliers) {
  writeLog(logfile, paste("Removing outlier clusters..."))
  if (opt$lowExp_col %in% colnames(obj_seurat@meta.data)) {
    rmLowClusters <- unique(obj_seurat$seurat_clusters[obj_seurat@meta.data[,opt$lowExp_col]==TRUE]) %>% 
      as.vector() %>% as.numeric()
    if (length(rmLowClusters)>0) {
      writeLog(logfile, paste("Removed low expression clusters: ",paste(rmLowClusters, collapse = ', ')))
      obj_seurat = subset(obj_seurat, idents = rmLowClusters, invert = TRUE)
    }
  } else {
    writeLog(logfile, paste("ERROR: Cannot find '",opt$lowExp_col,"' field in meta-data, low-expression clusters will not be removed.",sep=''))
  }
  if (opt$highExp_col %in% colnames(obj_seurat@meta.data)) {
    rmHighClusters <- unique(obj_seurat$seurat_clusters[obj_seurat@meta.data[,opt$highExp_col]==TRUE]) %>% 
      as.vector() %>% as.numeric()
    if (length(rmHighClusters)>0) {
      writeLog(logfile, paste("Removed high expression clusters: ",paste(rmHighClusters, collapse = ', ')))
      obj_seurat = subset(obj_seurat, idents = rmHighClusters, invert = TRUE)
    }
  } else {
    writeLog(logfile, paste("Cannot find '",opt$highExp_col,"' field in meta-data, high-expression clusters will not be removed.",sep=''))
  }
}

if (!is.na(opt$SCT_regress_vars)) {
   vars_to_regress = unlist(stringi::stri_split(str = opt$SCT_regress_vars,fixed = ','))
}else{
   vars_to_regress = c()
}

if (opt$CellCycleScoring){
    s.genes    = Seurat::cc.genes$s.genes
    g2m.genes  = Seurat::cc.genes$g2m.genes
    if (opt$CellCycle2regressOut){
        obj_seurat = CellCycleScoring(obj_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
        if (!opt$CellCycle2regressOut_CC_Difference){
            vars_to_regress = c(vars_to_regress,c("S.Score", "G2M.Score"))
            # obj_seurat = ScaleData(obj_seurat, vars.to.regress = vars_to_regress, features = rownames(obj_seurat))
        }else{
            obj_seurat$CC.Difference = obj_seurat$S.Score - obj_seurat$G2M.Score
            vars_to_regress = c(vars_to_regress,"CC.Difference")
            # obj_seurat               = ScaleData(obj_seurat, vars.to.regress = vars_to_regress, features = rownames(obj_seurat))
        }
    }else{
        obj_seurat = CellCycleScoring(obj_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    }
}



if (opt$SCTransform) {
  writeLog(logfile, paste("Using SCTransform based normalization instead of standard log-normalization",sep=''))
  writeLog(logfile, paste("Number of subsampling cells used to build NB regression: ",opt$SCT_ncells,sep=''))
  writeLog(logfile, paste("Number of variable features used after ranking by residual variance: ",opt$SCT_variable_features,sep=''))
  if (length(vars_to_regress)>0) {
    writeLog(logfile, paste("Variables to regress out in a second non-regularized linear regression: ",opt$SCT_regress_vars,sep=''))
    # vars_to_regress = unlist(stringi::stri_split(str = opt$SCT_regress_vars,fixed = ','))
    
    obj_seurat = Seurat::SCTransform(object = obj_seurat, ncells = opt$SCT_ncells, 
                                     variable.features.n = opt$SCT_variable_features, 
                                     vars.to.regress = vars_to_regress)
  } else {
    obj_seurat = Seurat::SCTransform(object = obj_seurat, ncells = opt$SCT_ncells, 
                                     variable.features.n = opt$SCT_variable_features)
  }
  
  # if (opt$UseMagic){
      # if (!is.na(opt$MagicCondaEnv)){
         # .libPaths(c(.libPaths(),paste(opt$MagicCondaEnv,"/lib/R/library",sep='')))
         # reticulate::use_condaenv(condaenv=opt$MagicCondaEnv, required =T)
      # }
      # library(reticulate)
      # library(Rmagic)
      # writeLog(logfile, paste("Imputing Data using Magic:"))
      # obj_seurat <- magic(obj_seurat,assay = "SCT")
      # obj_seurat@active.assay = "MAGIC_SCT"
      # # Identification of highly variable features (feature selection)
      # writeLog(logfile, paste("Identifying highly variable features..."))
      # obj_seurat <- FindVariableFeatures(obj_seurat, nfeatures = opt$VariableFeatures)
      # obj_seurat <- ScaleData(obj_seurat, features = rownames(obj_seurat))
  # }
  # Linear dimensional reduction
  obj_seurat <- RunPCA(object = obj_seurat, npcs = opt$JackStraw_dims, verbose = FALSE)
  dims_i = opt$JackStraw_dims
  
} else {
  # Normalizing the data
  writeLog(logfile, paste("Normalizing data..."))
  obj_seurat <- NormalizeData(obj_seurat)
  
  # Identification of highly variable features (feature selection)
  writeLog(logfile, paste("Identifying highly variable features..."))
  obj_seurat <- FindVariableFeatures(obj_seurat, nfeatures = opt$VariableFeatures)
  
  # Scaling the data
  writeLog(logfile, paste("Scaling data..."))
  # vars_to_regress = unlist(stringi::stri_split(str = opt$SCT_regress_vars,fixed = ','))
  obj_seurat <- ScaleData(obj_seurat, 
                          vars.to.regress = vars_to_regress,
                          features = rownames(obj_seurat))
  
  
  
  # if (opt$UseMagic){
      # if (!is.na(opt$MagicCondaEnv)){
         # .libPaths(c(.libPaths(),paste(opt$MagicCondaEnv,"/lib/R/library",sep='')))
         # reticulate::use_condaenv(condaenv=opt$MagicCondaEnv, required =T)
      # }
      # library(reticulate)
      # library(Rmagic)
      # writeLog(logfile, paste("Imputing Data using Magic:"))
      # obj_seurat <- magic(obj_seurat,assay = "RNA")
      # obj_seurat@active.assay = "MAGIC_RNA"
      # # Identification of highly variable features (feature selection)
      # writeLog(logfile, paste("Identifying highly variable features..."))
      # obj_seurat <- FindVariableFeatures(obj_seurat, nfeatures = opt$VariableFeatures)
      # obj_seurat <- ScaleData(obj_seurat, features = rownames(obj_seurat))
  # }
  # Determine the ‘dimensionality’ of the dataset
  writeLog(logfile, paste("Calculating dataset dimensionality..."))
  dims_seq <- seq.int(opt$JackStraw_dims, opt$JackStraw_dims+50, 10)
  i=1
  dimsFlag = FALSE
  
  while (isFALSE(dimsFlag)) {
    dims_i = dims_seq[i]
    writeLog(logfile, paste("Testing for ",dims_i," PCs to compute significance and score",sep=''))
    obj_i = obj_seurat
    
    # Linear dimensional reduction
    obj_i <- RunPCA(object = obj_i, npcs = dims_i, verbose = FALSE)
    
    # JackStraw
    obj_i <- JackStraw(obj_i, dims = dims_i)
    obj_i <- ScoreJackStraw(obj_i, dims = 1:dims_i, reduction = "pca")
    
    # Test for dimensions significance
    pAll <- JS(object = obj_i[["pca"]])
    nonSigDims = as.numeric(which(pAll$overall.p.values[,'Score'] > 0.05))
    if (length(nonSigDims) != 0) {
      SigDims = min(nonSigDims)-1
      writeLog(logfile, paste("Found ",SigDims," significant dimensions with p-value < 0.05",sep = ''))
      dimsFlag = TRUE
    } else {
      if (dims_i == dims_seq[length(dims_seq)]) {
        SigDims = dims_i
        writeLog(logfile, paste("All ",SigDims," dimensions with p-value < 0.05, continues with this amount of dimensions.",sep=''))
        dimsFlag = TRUE
      } else {
        writeLog(logfile, paste("All ",dims_i," dimensions with p-value < 0.05, Testing for higher dimensionality...",sep=''))
        i = i+1
        rm(obj_i)
      }
    }
  }
  obj_seurat = obj_i
}

# Plots
writeLog(logfile, paste("Plotting dimensionality figures..."))
if (!opt$SCTransform) {
  plt1 <- JackStrawPlot(obj_seurat, dims = 1:dims_i)
  jpeg(paste(opt$outDir,opt$Sample,"_jackstrawPlot.jpeg", sep = ""), 
       width = 1300, height = 750)
  print(plt1)
  dev.off()
}

plt2 <- ElbowPlot(obj_seurat, ndims = dims_i)
jpeg(paste(opt$outDir,opt$Sample,"_elbowPlot.jpeg", sep = ""), 
     width = 1000, height = 750)
print(plt2)
dev.off()


if (opt$UseMagic){
      if (!is.na(opt$MagicCondaEnv)){
         .libPaths(c(.libPaths(),paste(opt$MagicCondaEnv,"/lib/R/library",sep='')))
         reticulate::use_condaenv(condaenv=opt$MagicCondaEnv, required =T)
      }
      library(reticulate)
      library(Rmagic)
      writeLog(logfile, paste("Imputing Data using Magic:"))
      
      obj_seurat@active.assay = "RNA"
      
      writeLog(logfile, paste("Normalizing data For Magic ..."))
      obj_seurat <- NormalizeData(obj_seurat)
      
      # Identification of highly variable features (feature selection)
      # writeLog(logfile, paste("Identifying highly variable features..."))
      # obj_seurat <- FindVariableFeatures(obj_seurat, nfeatures = opt$VariableFeatures)
      
      # Scaling the data
      writeLog(logfile, paste("Scaling data..."))
      obj_seurat <- ScaleData(obj_seurat, 
                          vars.to.regress = vars_to_regress,
                          features = rownames(obj_seurat))
      obj_seurat <- magic(obj_seurat,assay = "RNA")
      # obj_seurat@active.assay = "MAGIC_RNA"
  }


# Save results
saveRDS(obj_seurat, file = paste(opt$outDir, opt$Sample, '.rds', sep=''))
writeLog(logfile, paste("Finished"))