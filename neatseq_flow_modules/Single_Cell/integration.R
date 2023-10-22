.libPaths()
library(optparse)
library(Seurat)
library(dplyr)
library(harmony)
library(future)
#library(cowplot)
library(ggpubr)
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--inputRDS"), type="character", default = NA,
              help="Vector of paths to Seurat object RDS files (comma-separated)", metavar = "character"),
  make_option(c("--ReIntegration"), type="character", default = NA,
              help="Path to Seurat integrated object for re-integration", metavar = "character"),
  make_option(c("--Samples"), type="character", default = NA,
              help="Vector of Samples names (comma-separated)", metavar = "character"),
  make_option(c("--ID"), type="character", default = NA,
              help="Project's ID to save integration results", metavar = "character"),
  make_option(c("--Subset_Samples"), type="character", default = NA,
              help="Vector of Samples to subset from all input data (comma-separated)", metavar = "character"),
  make_option(c("--VariableFeatures"), type="numeric", default = 3000,
              help="Number of features to select as top variable features (Default is 3000)", metavar = "character"),
  make_option(c("--JackStraw_dims"), type="numeric", default = 50,
              help="Number of PCs to compute significance for in JackStraw (Default is 50)", metavar = "character"),
  make_option(c("--CPUs"), type="numeric", default = 1,
              help="Number of threads to use (Default is 1)", metavar = "character"),
  make_option(c("--Memory"), type="numeric", default = NA,
              help="Allocate memory for Seurat parallelization in GB (Default is NA)", metavar = "character"),
  make_option(c("--Resolution"), type="numeric", default = 0.8,
              help="Value of the resolution parameter (Default is 0.8), use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities", metavar = "character"),
  make_option(c("--dims"), type="numeric", default = 50,
              help="Dimensions of reduction to use as input for clustering (Default is 50). Must include --overwrite_dims", metavar = "character"),
  make_option(c("--overwrite_dims"), action="store_true", default = FALSE,
              help="Cluster with manually number of dimensions (Set dims with --dims)", metavar = "character"),
  make_option(c("--SCTransform"), action="store_true", default = FALSE,
              help="Integrate datasets using Seurat with SCTransform", metavar = "character"),
  make_option(c("--UseMagic"), action="store_true", default = FALSE,
              help="Will use Magic to Impute the data after normalization", metavar = "character"),
  make_option(c("--MagicCondaEnv"), type="character", default = NA,
              help="If --UseMagic is set it will this Conda env to fined the package [Must Be a Full Path!]", metavar = "character"),
  make_option(c("-o", "--outDir"), type="character", default = NA,
              help="Path to the output directory", metavar = "character")
);

#Get user information
opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", 
                                    option_list=option_list,
                                    epilogue="\n\nAuthor:Gil Sorek");
opt = optparse::parse_args(opt_parser);

# Set log framework
logfile = paste(opt$outDir,opt$ID,'_log.txt', sep = "")
cat(paste('[',Sys.time(),']: Integration: ',opt$ID,sep=''), file=logfile, sep='\n')
writeLog <- function(logfile, msg) { cat(paste('(',strftime(Sys.time(), format = "%H:%M:%S"),'): ',msg,sep=''),file=logfile,append=TRUE,sep='\n') }
saveRDS(opt, paste(opt$outDir,'opt.rds',sep=''))

print("Input var:",quote = F)
print.data.frame(as.data.frame(x = unlist(opt),row.names = names(opt)),right = F,quote = F)

# Import data
if (!is.na(opt$inputRDS)) {
  writeLog(logfile, paste("Importing Seurat RDS objects...",sep=''))
  files = unlist(stringi::stri_split(str = opt$inputRDS,fixed = ','))
  if (is.na(opt$Samples)) {
    writeLog(logfile, paste("ERROR: Samples name must be specified [--Samples]"))
    stop()
  } else {
    writeLog(logfile, paste("Importing Sample names...",sep=''))
    samples = unlist(stringi::stri_split(str = opt$Samples,fixed = ','))
    if (length(files)!=length(samples)){
      writelog(logfile, paste("ERROR: The number of sample's names is not equal to the number of sample's files"))
      stop()
    }
  }
  if (!is.na(opt$Subset_Samples)) {
    writeLog(logfile, paste("Subsetting Samples...",sep=''))
    subset_samples <- unlist(stringi::stri_split(str = opt$Subset_Samples,fixed = ','))
    if (!all(subset_samples %in% samples)) {
      writeLog(logfile, paste("Cannot find all subset samples: ",paste(subset_samples,collapse = ","),sep=''))
      stop()
    }
    writeLog(logfile, paste("Using only the following sample(s) for integration: ",paste(subset_samples,collapse = ","),sep=''))
    files = files[which(samples %in% subset_samples)]
    samples = samples[which(samples %in% subset_samples)]
  }
} else {
  if (is.na(opt$ReIntegration)) {
    writeLog(logfile, paste("ERROR: Seurat object RDS files must be specified [--inputRDS]"))
    stop()
  }
}
if (!is.na(opt$ReIntegration)) {
  writeLog(logfile, paste("Importing Seurat Object for re-Integration...",sep=''))
  reinteg_obj = readRDS(opt$ReIntegration)
}
if (is.na(opt$outDir)) {
  writeLog(logfile, paste("ERROR: Output directory must be specified [--outDir]"))
  stop()
}
if (!is.na(opt$ID)) {
  dir_path = paste(opt$outDir,opt$ID,sep='')
  if (!dir.exists(dir_path))
    dir.create(dir_path)
} else {
  writeLog(logfile, paste("Project ID must be specified [--ID]"))
}
if (opt$overwrite_dims) {
  SigDims <- opt$dims
  writeLog(logfile, paste("Skipping dimensionality testing, will use ",SigDims,' dimensions',sep=''))
}

# Parallelization
if (opt$CPUs>1) {
  writeLog(logfile, paste("Using multicore parallelization with ",opt$CPUs,' threads',sep=''))
  plan("multicore", workers = opt$CPUs)
}
  
if (!is.na(opt$Memory)) {
    writeLog(logfile, paste("Allocating ",as.numeric(opt$Memory),'GB of memory',sep=''))
    options(future.globals.maxSize = as.numeric(opt$Memory) * 1024 ^ 3)
}


seurat_list = list()
# Load RDS files
if (!is.na(opt$inputRDS)) {
  writeLog(logfile, paste("Loading RDS files..."))
  for (sample in samples) {
    file = files[which(sample==samples)]
    if (file.exists(file)) {
      seurat_obj = readRDS(file)
      seurat_obj@meta.data[,stringr::str_detect(colnames(seurat_obj@meta.data),c('snn_res','seurat_clusters'))] <- NULL
      seurat_list <- c(seurat_list, seurat_obj)
      rm(seurat_obj)
    } else {
      writeLog(logfile, paste("ERROR: Couldn't find ",file,sep=''))
      stop()
    }
  }
}
if (!is.na(opt$ReIntegration)) {
  idents = unique(reinteg_obj$orig.ident)
  for (ident in idents) {
    seurat_obj = subset(reinteg_obj, subset = orig.ident == ident)
    seurat_list <- c(seurat_list, seurat_obj)
  }
}


# Integration type
if (opt$SCTransform) {
  writeLog(logfile, paste("Integration method: Seurat (using SCTransform)"))
  writeLog(logfile, paste("Selecting Integration features..."))
  integ_features = SelectIntegrationFeatures(object.list = seurat_list, nfeatures = opt$VariableFeatures)
  
  writeLog(logfile, paste("Preparing for SCT Integration..."))
  seurat_list = PrepSCTIntegration(object.list = seurat_list, anchor.features = integ_features, verbose = FALSE)
  
  writeLog(logfile, paste("Finding Integration anchors..."))
  integ_anchors = FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT",
                                         anchor.features = integ_features, dims = 1:SigDims, verbose = FALSE)
  
  writeLog(logfile, paste("Integrating data..."))
  integ_obj = IntegrateData(anchorset = integ_anchors, normalization.method = "SCT",
                            dims = 1:SigDims, verbose = FALSE)
  
  if (opt$UseMagic){
     if (!is.na(opt$MagicCondaEnv)){
         .libPaths(c(.libPaths(),paste(opt$MagicCondaEnv,"/lib/R/library",sep='')))
         reticulate::use_condaenv(condaenv=opt$MagicCondaEnv, required =T)
     }
     library(reticulate)
     library(Rmagic)
     writeLog(logfile, paste("Imputing Data using Magic:"))
     integ_obj <- magic(integ_obj)
     integ_obj@active.assay = "MAGIC_integrated"
  }
  
  # Identification of highly variable features (feature selection)
  writeLog(logfile, paste("Identifying highly variable features..."))
  integ_obj <- FindVariableFeatures(integ_obj, nfeatures = opt$VariableFeatures)
  
  writeLog(logfile, paste("Scaling data..."))
  integ_obj <- ScaleData(integ_obj, features = rownames(integ_obj))
  
  writeLog(logfile, paste("Running PCA..."))
  integ_obj <- RunPCA(object = integ_obj, npcs = SigDims, verbose = FALSE)
  
  # Clustering
  writeLog(logfile, paste("Clustering object..."))
  integ_obj <- integ_obj %>%
    RunUMAP(dims = 1:SigDims) %>% 
    FindNeighbors(dims = 1:SigDims) %>% 
    FindClusters(resolution = opt$Resolution)
  
  integ_obj <- RUNTSNE(integ_obj,dims = 1:SigDims)
  
} else {
  writeLog(logfile, paste("Integration method: Harmony"))
  
  # Create merged Seurat object
  writeLog(logfile, paste("Creating merged Seurat object..."))
  merge_obj <- merge(seurat_list[[1]], 
                     y = seurat_list[2:length(seurat_list)], 
                     #add.cell.ids = names(samples), 
                     project = opt$ID)
  
  # Normalizing the data
  writeLog(logfile, paste("Normalizing data..."))
  merge_obj <- NormalizeData(merge_obj)
  

  if (opt$UseMagic){
     if (!is.na(opt$MagicCondaEnv)){
         .libPaths(c(.libPaths(),paste(opt$MagicCondaEnv,"/lib/R/library",sep='')))
         reticulate::use_condaenv(condaenv=opt$MagicCondaEnv, required =T)
     }
     library(reticulate)
     library(Rmagic)
     writeLog(logfile, paste("Imputing Data using Magic:"))
     merge_obj <- magic(merge_obj)
     merge_obj@active.assay = "MAGIC_integrated"
  }
  
    # Identification of highly variable features (feature selection)
  writeLog(logfile, paste("Identifying highly variable features..."))
  merge_obj <- FindVariableFeatures(merge_obj, nfeatures = opt$VariableFeatures)
  
  # Scaling the data
  writeLog(logfile, paste("Scaling data..."))
  merge_obj <- ScaleData(merge_obj, features = rownames(merge_obj))
  
  
  
  if (isFALSE(opt$overwrite_dims)) {
    # Determine the ‘dimensionality’ of the dataset
    writeLog(logfile, paste("Calculating dataset dimensionality..."))
    dims_seq <- seq.int(opt$JackStraw_dims, opt$JackStraw_dims+50, 10)
    i=1
    dimsFlag = FALSE
    
    while (isFALSE(dimsFlag)) {
      dims_i = dims_seq[i]
      writeLog(logfile, paste("Testing for ",dims_i," PCs to compute significance and score",sep=''))
      obj_i = merge_obj
      
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
    
    # JackStraw Plots
    plt1 <- JackStrawPlot(obj_i, dims = 1:SigDims)
    pdf(paste(opt$outDir,opt$ID,"/",opt$ID,"_jackstrawPlot.pdf", sep = ""), 
         width = 20, height = 15)
    print(plt1)
    dev.off()
    plt2 <- ElbowPlot(obj_i, ndims = SigDims)
    pdf(paste(opt$outDir,opt$ID,"/",opt$ID,"_elbowPlot.pdf", sep = ""), 
         width = 20, height = 15)
    print(plt2)
    dev.off()
  } else {
    obj_i <- RunPCA(object = merge_obj, npcs = SigDims, verbose = FALSE)
  }
  
  # RunHaromny
  writeLog(logfile, paste("Running Harmony..."))
  integ_obj <- RunHarmony(obj_i, "orig.ident", plot_convergence = FALSE)
  
  # Embeddings
  writeLog(logfile, paste("Saving Embeddings..."))
  harmony_embeddings <- Embeddings(integ_obj, 'harmony')
  write.csv(harmony_embeddings, file = paste(opt$outDir,opt$ID,"/harmony_embeddings.csv", sep = ""))
  
  
  # Plots
  writeLog(logfile, paste("Plotting Harmony figures..."))
  p1 <- DimPlot(object = integ_obj, reduction = "harmony", pt.size = 0.7, group.by = "orig.ident")
  p2 <- VlnPlot(object = integ_obj, features = "harmony_1", group.by = "orig.ident", pt.size = 0)
  pdf(paste(opt$outDir,opt$ID,"/harmony_PCs.pdf", sep = ""), width = 20, height = 15)
  cowplot::plot_grid(p1,p2)
  dev.off()
  
  # Clustering
  writeLog(logfile, paste("Clustering object..."))
  integ_obj <- integ_obj %>%
    RunUMAP(reduction = "harmony", dims = 1:SigDims) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:SigDims) %>% 
    FindClusters(resolution = opt$Resolution)
    
    integ_obj <- RUNTSNE(integ_obj,reduction = "harmony",dims = 1:SigDims)
}


legend_pos = "right"
writeLog(logfile, paste("Plotting Integrated object figures..."))
p1 <- DimPlot(integ_obj, reduction = "umap", group.by = "orig.ident", pt.size = 1.5)+
              theme(legend.position = legend_pos)
p2 <- DimPlot(integ_obj, reduction = "umap", label = T, pt.size = 1.5, label.size = 8)+
  ggplot2::theme(legend.position = "none")
pdf(paste(opt$outDir,opt$ID,'/',opt$ID,'_umapBySample.pdf', sep = ""), width = 20, height = 15)
print(p1)
dev.off()
pdf(paste(opt$outDir,opt$ID,'/',opt$ID,'_umapByCluster.pdf', sep = ""), width = 20, height = 15)
print(p2)
dev.off()

p1 <- DimPlot(integ_obj, reduction = "tsne", group.by = "orig.ident", pt.size = 1.5)+
              theme(legend.position = legend_pos)
p2 <- DimPlot(integ_obj, reduction = "tsne", label = T, pt.size = 1.5, label.size = 8)+
  ggplot2::theme(legend.position = "none")
pdf(paste(opt$outDir,opt$ID,'/',opt$ID,'_tsneBySample.pdf', sep = ""), width = 20, height = 15)
print(p1)
dev.off()
pdf(paste(opt$outDir,opt$ID,'/',opt$ID,'_tsneByCluster.pdf', sep = ""), width = 20, height = 15)
print(p2)
dev.off()


# Clusters distribution per Sample
samples = unique(integ_obj$orig.ident)
samples = samples[order(samples)]
df = data.frame()
for (sample in samples) {
  sample_df = data.frame(table(integ_obj$seurat_clusters[integ_obj$orig.ident==sample]))
  colnames(sample_df) = c('Cluster','Cells')
  sample_df$Ratio = sample_df$Cells / length(integ_obj$seurat_clusters[integ_obj$orig.ident==sample])
  sample_df$Sample = sample
  df = rbind(df, sample_df)
}
samples_plot = ggplot(df, aes(x=Cluster,y=Ratio,fill=Sample)) +
  geom_col(position = position_dodge()) +
  ggtitle("Clusters distribution per Sample")
ggplot2::ggsave(filename = paste(opt$outDir,opt$ID,'/',opt$ID,'_ClustersDistPerSample.pdf', sep = ""),
                samples_plot, dpi=600, width=12, height=6,device='pdf')

samples_plot = ggplot(df)+ aes(x=Sample,y=Ratio,fill=Cluster) + 
               geom_bar(stat="identity", color = 'white', size = 0.5)+
               ggtitle("Clusters distribution per Sample")

if (length(unique( df$Cluster))<12){
    samples_plot = samples_plot + scale_fill_brewer(palette="Paired")
}
ggplot2::ggsave(filename = paste(opt$outDir,opt$ID,'/',opt$ID,'_PerSample_ClustersDist.pdf', sep = ""),
                samples_plot, dpi=600, width=12, height=6,device='pdf')
write.csv(x = df,file = paste(opt$outDir,opt$ID,'/',opt$ID,'ClustersDist.csv', sep = ""))

# Samples distribution per Cluster
clusters=as.numeric(as.vector(unique(integ_obj$seurat_clusters)))
clusters = clusters[order(clusters)]
df = data.frame()
for (cluster in clusters) {
  cluster_df = data.frame(table(integ_obj$orig.ident[integ_obj$seurat_clusters==cluster]))
  colnames(cluster_df) = c('Sample','Cells')
  cluster_df$Ratio = cluster_df$Cells / length(integ_obj$orig.ident[integ_obj$seurat_clusters==cluster])
  cluster_df$Cluster = cluster
  df = rbind(df, cluster_df)
}
clusters_plot = ggplot(df, aes(x=Cluster,y=Ratio,fill=Sample)) +
  geom_col(position = position_stack()) +
  geom_hline(yintercept = 0.5, lty = "dashed") +
  scale_x_continuous("Cluster", labels = as.character(clusters), breaks = clusters) +
  ggtitle("Samples distribution per Cluster")
ggplot2::ggsave(filename = paste(opt$outDir,opt$ID,'/',opt$ID,'_SamplesDistPerCluster.pdf', sep = ""),
                clusters_plot, dpi=600, width=12, height=6,device='pdf')

# Save results
writeLog(logfile, paste("Saving results..."))
saveRDS(integ_obj, file = paste(opt$outDir,opt$ID,'/integrated_obj.rds',sep=''))
writeLog(logfile, paste("Finished"))