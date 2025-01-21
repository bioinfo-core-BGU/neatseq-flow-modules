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
  make_option(c("--REFFRDS"), type="character", default = NA,
              help="Path to Seurat object RDS file That will be used as Reference for Cell-Type Annotation", metavar = "character"),
  make_option(c("--REFFCellTypeCol"), type="character", default = NA,
              help="The Reference MetaData Column name to use for Cell-Type Annotation", metavar = "character"),
  make_option(c("--dims"), type="numeric", default = 50,
              help="Dimensions of reduction to use as input for clustering (Default is 50). Must include --overwrite_dims", metavar = "character"),
  make_option(c("--Resolution"), type="numeric", default = 0.8,
              help="Value of the resolution parameter (Default is 0.8), use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities", metavar = "character"),
  make_option(c("--Sample"), type="character", default = NA,
              help="Sample's name", metavar = "character"),
  make_option(c("--reduction"), type="character", default = NA,
              help="Which dimensional reduction to use in FindNeighbors and RunUMAP (Default is PCA)", metavar = "character"),
  make_option(c("--Alpha"), type="numeric", default = 0.05,
              help="Significance level cutoff for JackStraw (Default is 0.05)", metavar = "character"),
  make_option(c("--human_adipose_markers"), action="store_true", default = FALSE,
              help="Plot human adipocytes markers FeaturePlot (ADIPOQ+AQP7)", metavar = "character"),
  make_option(c("--mouse_adipose_markers"), action="store_true", default = FALSE,
              help="Plot mouse adipocytes markers FeaturePlot (Adipoq+Aqp7)", metavar = "character"),
  make_option(c("--Show_Legend"), action="store_true", default = FALSE,
              help="Show legend in figures (Default is False)", metavar = "character"),
  make_option(c("--tSNE"), action="store_true", default = FALSE,
              help="Run 'RunTSNE' NOT 'RunUMAP' (Default is False)", metavar = "character"),
  make_option(c("--slingshot"), action="store_true", default = FALSE,
              help="Perform a slingshot Trajectory analysis (Default is False)", metavar = "character"),
  make_option(c("--outDir"), type="character", default = NA,
              help="Path to the output directory", metavar = "character"),
  make_option(c("--overwrite_dims"), action="store_true", default = FALSE,
              help="Cluster with manually number of dimensions (Set dims with --dims)", metavar = "character")
);

#Get user information
opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", 
                                    option_list=option_list,
                                    epilogue="\n\nAuthor:Gil Sorek");
opt = optparse::parse_args(opt_parser);

# Set log framework
logfile = paste(opt$outDir,opt$Sample,'_log.txt',sep='')
cat(paste('[',Sys.time(),']: clustering: ',opt$Sample,sep=''), file=logfile, sep='\n')
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

if (opt$tSNE){
    reduction = "tsne"
}else{
    reduction = "umap"
}


# if (opt$CellCycleScoring){
    # s.genes    = Seurat::cc.genes$s.genes
    # g2m.genes  = Seurat::cc.genes$g2m.genes
    # if (opt$CellCycle2regressOut){
        # obj_seurat = CellCycleScoring(obj_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
        # if (!opt$CellCycle2regressOut_CC_Difference){
            # obj_seurat = ScaleData(obj_seurat, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(obj_seurat))
        # }else{
            # obj_seurat$CC.Difference = obj_seurat$S.Score - obj_seurat$G2M.Score
            # obj_seurat               = ScaleData(obj_seurat, vars.to.regress = "CC.Difference", features = rownames(obj_seurat))
        # }
    # }else{
        # obj_seurat = CellCycleScoring(obj_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    # }
# }

# Set number of dimensions for clustering
if (obj_seurat@active.assay == "RNA") {
  writeLog(logfile, paste("Active assay detected: RNA"))
  if (opt$overwrite_dims) {
    SigDims = opt$dims
  } else {
    writeLog(logfile, paste("Significance level cutoff: ",opt$Alpha))
    pAll <- JS(object = obj_seurat[["pca"]])
    nonSigDims = as.numeric(which(pAll$overall.p.values[,'Score'] > opt$Alpha))
    if (length(nonSigDims) != 0) {
      SigDims = min(nonSigDims)-1
    } else {
      SigDims = dim(pAll$overall.p.values)[1]
    }
  }
} else {
  writeLog(logfile, paste("Active assay detected: ",obj_seurat@active.assay,sep=''))
  SigDims = opt$dims
}

writeLog(logfile, paste("Using ",SigDims," dimensions of reduction as input for clustering",sep=''))

# Reduction usage
if (is.na(opt$reduction)) {
  writeLog(logfile, paste("Using 'pca' dimensional reduction as input for FindNeighbors and ",reduction,sep=''))
  use_reducs = 'pca'
} else {
  if (opt$reduction %in% names(obj_seurat@reductions)) {
	writeLog(logfile, paste("Using '",opt$reduction,"' dimensional reduction as input for FindNeighbors and ",reduction,sep=''))
	use_reducs = opt$reduction
  } else {
	writeLog(logfile, paste("ERROR: Unable to find reduction '",opt$reduction,"'. Available reductions: ",paste(names(obj_seurat@reductions),collapse=','),sep=''))
	stop()
  }
}



# Cluster the Cells
writeLog(logfile, paste("Clustering cells..."))
obj_seurat <- FindNeighbors(obj_seurat, reduction = use_reducs, dims = 1:SigDims)
obj_seurat <- FindClusters(obj_seurat, resolution = opt$Resolution)
if (opt$tSNE){
    writeLog(logfile, paste("Running tSNE..."))
    obj_seurat <- RunTSNE(obj_seurat, reduction = use_reducs, dims = 1:SigDims)
}else{
    writeLog(logfile, paste("Running UMAP..."))
    obj_seurat <- RunUMAP(obj_seurat, reduction = use_reducs, dims = 1:SigDims)
}



if (!is.na(opt$REFFRDS)){
    Reff_seurat <- readRDS(opt$REFFRDS)
    
    if (is.na(opt$REFFCellTypeCol)){
        REFFCellTypeCol = Reff_seurat@active.ident
    }else{
        REFFCellTypeCol = opt$REFFCellTypeCol
    }


    anchors    <- FindTransferAnchors(reference = Reff_seurat,
                                  query   = obj_seurat,
                                  dims    = 1:SigDims,
                                  verbose = T)
    obj_seurat <- TransferData(anchorset = anchors,
                               query     = obj_seurat,
                               reference = Reff_seurat,
                               refdata   = REFFCellTypeCol,
                               verbose   = T)
    obj_seurat <- SetIdent(obj_seurat, value = "predicted.id")
    
    samples = unique(obj_seurat$orig.ident)
    samples = samples[order(samples)]
    df = data.frame()
    for (sample in samples) {
      sample_df = data.frame(table(obj_seurat@active.ident[obj_seurat$orig.ident==sample]))
      colnames(sample_df) = c('Cluster','Cells')
      sample_df$Ratio = sample_df$Cells / length(obj_seurat@active.ident[obj_seurat$orig.ident==sample])
      sample_df$Sample = sample
      df = rbind(df, sample_df)
    }
    write.csv(df, file=paste(opt$outDir,"Cell_Prop.csv", sep = ""))
}




# Plots
legend_pos = "none"
if (opt$Show_Legend)
  legend_pos = "right"
writeLog(logfile, paste("Plotting reduction figures..."))
plt1 <- DimPlot(obj_seurat, reduction = reduction, label = T, pt.size = 1.5, label.size = 8)+
  theme(legend.position = legend_pos)
jpeg(paste(opt$outDir,opt$Sample,"_",reduction,"ByCluster.jpeg", sep = ""), 
     width = 1500, height = 1250)
print(plt1)
dev.off()
if (length(unique(obj_seurat$orig.ident))>1) {
  plt2 <- DimPlot(obj_seurat, reduction = reduction, group.by = "orig.ident", pt.size = 1.5)
  jpeg(paste(opt$outDir,opt$Sample,"_",reduction,"BySample.jpeg", sep = ""), width = 1500, height = 1500)
  print(plt2)
  dev.off()
}

if ( "Phase" %in% colnames(obj_seurat@meta.data)) {
    writeLog(logfile, paste("Plotting Cell Cycle Scoring reduction figures..."))
    plt3 <- DimPlot(obj_seurat,
                    reduction = reduction,
                    group.by = "Phase",
                    label = T, pt.size = 1.5,
                    label.size = 8)+
      theme(legend.position = legend_pos)
    jpeg(paste(opt$outDir,opt$Sample,"_",reduction,"ByCellCycle.jpeg", sep = ""), 
         width = 1500, height = 1250)
    print(plt3)
    dev.off()
}

plt3 <- VlnPlot(obj_seurat, features = "nCount_RNA", pt.size = 0)+ 
  geom_boxplot(width=0.3)+
  theme(legend.position = legend_pos) +
  scale_y_continuous(trans = "log10")
plt4 <- VlnPlot(obj_seurat, features = "nFeature_RNA", pt.size = 0)+ 
  geom_boxplot(width=0.3)+
  theme(legend.position = legend_pos) +
  scale_y_continuous(trans = "log10")
jpeg(paste(opt$outDir,opt$Sample,"_violinPlots.jpeg", sep = ""), 
     width = 1800, height = 1400)
print(plt3+plt4)
dev.off()
if (opt$human_adipose_markers) {
  plt5 <- FeaturePlot(obj_seurat, features = c("AQP7", "ADIPOQ"))
  jpeg(paste(opt$outDir,opt$Sample,"_featurePlot_adipocytes.jpeg", sep = ""), 
       width = 2000, height = 1200)
  print(plt5)
  dev.off()
} else if (opt$mouse_adipose_markers) {
  plt5 <- FeaturePlot(obj_seurat, features = c("Aqp7", "Adipoq"))
  jpeg(paste(opt$outDir,opt$Sample,"_featurePlot_adipocytes.jpeg", sep = ""), 
       width = 2000, height = 1200)
  print(plt5)
  dev.off()
}

if (length(unique(obj_seurat$orig.ident))>1) {
# Clusters distribution per Sample
samples = unique(obj_seurat$orig.ident)
samples = samples[order(samples)]
df = data.frame()
for (sample in samples) {
  sample_df = data.frame(table(obj_seurat$seurat_clusters[obj_seurat$orig.ident==sample]))
  colnames(sample_df) = c('Cluster','Cells')
  sample_df$Ratio = sample_df$Cells / length(obj_seurat$seurat_clusters[obj_seurat$orig.ident==sample])
  sample_df$Sample = sample
  df = rbind(df, sample_df)
}
samples_plot = ggplot(df, aes(x=Cluster,y=Ratio,fill=Sample)) +
  geom_col(position = position_dodge()) +
  ggtitle("Clusters distribution per Sample")
ggplot2::ggsave(filename = paste(opt$outDir,opt$Sample,'_ClustersDistPerSample.jpeg', sep = ""),
                samples_plot, dpi=600, width=12, height=6)

# Samples distribution per Cluster
clusters=as.numeric(as.vector(unique(obj_seurat$seurat_clusters)))
clusters = clusters[order(clusters)]
df = data.frame()
for (cluster in clusters) {
  cluster_df = data.frame(table(obj_seurat$orig.ident[obj_seurat$seurat_clusters==cluster]))
  colnames(cluster_df) = c('Sample','Cells')
  cluster_df$Ratio = cluster_df$Cells / length(obj_seurat$orig.ident[obj_seurat$seurat_clusters==cluster])
  cluster_df$Cluster = cluster
  df = rbind(df, cluster_df)
}
clusters_plot = ggplot(df, aes(x=Cluster,y=Ratio,fill=Sample)) +
  geom_col(position = position_stack()) +
  geom_hline(yintercept = 0.5, lty = "dashed") +
  scale_x_continuous("Cluster", labels = as.character(clusters), breaks = clusters) +
  ggtitle("Samples distribution per Cluster")
ggplot2::ggsave(filename = paste(opt$outDir,opt$Sample,'_SamplesDistPerCluster.jpeg', sep = ""),
                clusters_plot, dpi=600, width=12, height=6)
}

if (opt$slingshot){
    library(slingshot)
    if (reduction == "tsne"){
        sce <- slingshot(as.matrix(obj_seurat@reductions$tsne@cell.embeddings),
                                omega = T, 
                                omega_scale = 1.5,
                                stretch=1.5,
                                clusterLabels = obj_seurat@active.ident)
    }else{
        sce <- slingshot(as.matrix(obj_seurat@reductions$umap@cell.embeddings),
                                omega = T, 
                                omega_scale = 1.5,
                                stretch=1.5,
                                clusterLabels = obj_seurat@active.ident)
    
    }
    Pseudotime = as.data.frame(slingPseudotime(sce))
    data       = SlingshotDataSet(sce)
    curves     = slingCurves(sce, as.df = TRUE)
    mst        = slingMST(sce, as.df = F)
    
    curves$x = curves[,1]
    curves$y = curves[,2]
    Pseudotime = Pseudotime[rownames(obj_seurat@meta.data),]
    obj_seurat = Seurat::AddMetaData(object =obj_seurat, Pseudotime)
  
    plt_Slingshot = Seurat::DimPlot(object =obj_seurat,reduction = reduction,label = T)+
                                    geom_path(data = curves %>% arrange(Order),
                                              arrow = arrow(length = unit(0.1, "inches")),
                                              aes(x=x,y=y,group = Lineage ),size=1, colour = "black")+
                                    facet_wrap(vars(Lineage))
    jpeg(paste(opt$outDir,opt$Sample,"_",reduction,"Slingshot.jpeg", sep = ""), 
    width = 1500, height = 1250)
    print(plt_Slingshot)
    dev.off()
    for (Lineage in colnames(Pseudotime)){
        Lineage_id = as.integer(strsplit(x = Lineage ,split  = "Lineage" )[[1]][2])
        plt_Lineage =Seurat::FeaturePlot(object =obj_seurat,reduction = reduction,features = Lineage)+
                            geom_path(data = curves[curves$Lineage==Lineage_id,] %>% arrange(Order),
                            arrow = arrow(length = unit(0.1, "inches")),
                            aes(x=x,y=y,group = Lineage ), size = 1)
        jpeg(paste(opt$outDir,opt$Sample,"_",reduction,"Slingshot_",Lineage,".jpeg", sep = ""), 
        width = 1500, height = 1250)
        print(plt_Lineage)
        dev.off()
    }
    
    

}

# Save results
saveRDS(obj_seurat, file = paste(opt$outDir, opt$Sample, '.rds', sep=''))
writeLog(logfile, paste("Finished"))