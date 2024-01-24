.libPaths()
library(optparse)
library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(ggplot2)
library(grid)
library(RColorBrewer)
library(patchwork)
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--inputRDS"), type="character", default = NA,
              help="Path to Seurat object RDS file", metavar = "character"),
  make_option(c("-s", "--Sample"), type="character", default = NA,
              help="Sample's name", metavar = "character"),
  make_option(c("-o", "--outDir"), type="character", default = NA,
              help="Path to the output directory", metavar = "character"),
  make_option(c("--Features_file"), type="character", default = NA,
              help="Path to file with location of feature files to plot", metavar = "character"),
  make_option(c("--UMAP"), action="store_true", default = FALSE,
              help="Plot UMAP", metavar = "character"),
  make_option(c("--TSNE"), action="store_true", default = FALSE,
              help="Plot TSNE", metavar = "character"),
  make_option(c("--RUNUMAP"), action="store_true", default = FALSE,
              help="Calculate  a New UMAP", metavar = "character"),
  make_option(c("--RUNTSNE"), action="store_true", default = FALSE,
              help="Calculate  a New tSNE", metavar = "character"),
  make_option(c("--reduction"), type="character", default = NA,
              help="Which dimensional reduction to use in FindNeighbors and RunUMAP (Default is PCA)", metavar = "character"),
  make_option(c("--Assay"), type="character", default = "RNA",
              help="Which Assay to use, the default is RNA. if set to 'active' will use the current active assay)", metavar = "character"),
  make_option(c("--UseMagic"), action="store_true", default = FALSE,
              help="Will use Magic to Impute the data after normalization, Will set --Assay to RNA", metavar = "character"),
  make_option(c("--MagicCondaEnv"), type="character", default = NA,
              help="If --UseMagic is set it will this Conda env to fined the package [Must Be a Full Path!]", metavar = "character"),
  make_option(c("--dims"), type="numeric", default = 50,
              help="Dimensions of reduction to use as input for clustering (Default is 50). Must include --overwrite_dims", metavar = "character"),
  make_option(c("--Heatmap"), action="store_true", default = FALSE,
              help="Plot input features heatmap", metavar = "character"),
  make_option(c("--Show_Heatmap_Legend"), action="store_true", default = FALSE,
              help="Show heatmap legend", metavar = "character"),
  make_option(c("--ViolinPlot"), action="store_true", default = FALSE,
              help="Plot input features violin-plot", metavar = "character"),
  make_option(c("--DotPlot"), action="store_true", default = FALSE,
              help="Plot input features dot-plot", metavar = "character"),
  make_option(c("--StackedViolinPlot"), action="store_true", default = FALSE,
              help="Plot input features stacked violin-plots", metavar = "character"),
  make_option(c("--FeaturePlot"), action="store_true", default = FALSE,
              help="Plot input features feature-plot", metavar = "character")
);

#Get user information
opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", 
                                    option_list=option_list,
                                    epilogue="\n\nAuthor:Gil Sorek");
opt = optparse::parse_args(opt_parser);

# Set log framework
logfile = paste(opt$outDir,opt$Sample,'_log.txt',sep='')
cat(paste('[',Sys.time(),']: Plotting: ',opt$Sample,sep=''), file=logfile, sep='\n')
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
} else {
  if (!dir.exists(paste(opt$outDir,opt$Sample,sep='')))
    dir.create(paste(opt$outDir,opt$Sample,sep=''))
}
if (is.na(opt$outDir)) {
  writeLog(logfile, paste("ERROR: Output directory must be specified [--outDir]"))
  stop()
}
if (!is.na(opt$Features_file)) {
  features_list = list()
  featureLocations <- read.table(opt$Features_file, sep = '\t', header = TRUE)
  if (isFALSE(all(colnames(featureLocations) == c('Source','Path')))) {
    writeLog(logfile, paste("ERROR: Invalid column names. Please use Source and Path"))
	stop()
  }
  for (src in unique(featureLocations$Source)) {
    file_path = featureLocations$Path[featureLocations$Source==src]
	writeLog(logfile, paste("Reading '",src,"' features list...",sep=''))
	if (file.exists(file_path)) {
	  features_list[[src]] = read.table(file_path, sep = '\t', header = TRUE)
	  if (isFALSE(all(colnames(features_list[[src]]) == c('Feature','CellType')))) {
        writeLog(logfile, paste("ERROR: Invalid column names. Please use Feature and CellType"))
	    stop()
      }
	} else {
	  writeLog(logfile, paste("ERROR: ",file_path," does not exists.",sep=''))
	  stop()
	}
  }
}



# Set active assay for visualization
if (opt$Assay!="active"){
    writeLog(logfile, paste(c("Setting active assay as ",opt$Assay,"for data visualization..."),sep=""))
    obj_seurat@active.assay = opt$Assay
}else{
    writeLog(logfile, paste("Using active assay for data visualization..."))
}


if (opt$UseMagic){

     if (!is.na(opt$MagicCondaEnv)){
         .libPaths(c(.libPaths(),paste(opt$MagicCondaEnv,"/lib/R/library",sep='')))
         reticulate::use_condaenv(condaenv=opt$MagicCondaEnv, required =T)
     }
     library(reticulate)
     library(Rmagic)
     writeLog(logfile, "Setting active assay as RNA for data visualization...")
     obj_seurat@active.assay = "RNA"
     writeLog(logfile, paste("Imputing Data using Magic:"))
     obj_seurat <- NormalizeData(object = obj_seurat)
     obj_seurat <- ScaleData(object = obj_seurat)
     obj_seurat <- magic(obj_seurat,assay = "RNA")
     obj_seurat@active.assay = "MAGIC_RNA"
}







# Plot Heatmap
if (opt$Heatmap) {
  writeLog(logfile, paste("Plotting Heatmap...",sep=''))
  if (is.na(opt$Features_file)) {
    writeLog(logfile, paste("ERROR: Please provide cell-type marker genes list [--Features_file]"))
	stop()
  }
  for (src in names(features_list)) {
  featureData = features_list[[src]]
  if (!dir.exists(paste(opt$outDir,opt$Sample,"/",src,sep='')))
    dir.create(paste(opt$outDir,opt$Sample,"/",src,sep=''))

  # Filter features present in expression matrix
  feat2use = featureData$Feature[featureData$Feature %in% rownames(obj_seurat)]
  feat2use = featureData[featureData$Feature %in% feat2use,]
  feat2use = feat2use[order(feat2use$CellType),]
  writeLog(logfile, paste("ATTENTION: The following markers were not found: ",
						paste(featureData$Feature[which(!(featureData$Feature %in% feat2use$Feature))],collapse=','),sep=''))
  
  # Determine row clusters
  row_clusters = feat2use$CellType %>% as.vector %>% as.character()
  names(row_clusters) = feat2use$Feature
  
  # Determine column clusters
  col_clusters = obj_seurat$seurat_clusters %>% as.vector() %>% as.numeric()
  names(col_clusters) = colnames(obj_seurat)
  col_clusters = col_clusters[order(col_clusters)]
  
  # Extract scaled and normalized counts
  counts = GetAssayData(object = obj_seurat, slot = "data")
  if (length(feat2use$Feature)>1) {
	  counts = counts[feat2use$Feature,]
	  counts = as.matrix(counts)
	  no_var=names(which(apply(counts,MARGIN = 2,FUN = sd)==0))
	  counts = counts[,!(colnames(counts) %in% no_var)]
	  seurat_clusters = obj_seurat$seurat_clusters[!(names(obj_seurat$seurat_clusters) %in% no_var)]
	  col_clusters = col_clusters[!(names(col_clusters) %in% no_var)]
	  
	  # Re-order clusters by cell expression
	  counts = t(counts)
	  New_clusters = c()
	  for (num in unique(col_clusters)) {
	    if (sum(seurat_clusters==num)>1) {
			mat = counts[seurat_clusters==num,]
			mat_novar = apply(mat,MARGIN = 2,sd) != 0
			if (sum(mat_novar)>1) {
			  mat = mat[,which(mat_novar)]
			  heat_map = pheatmap::pheatmap(mat = mat, cluster_rows = TRUE, show_colnames = FALSE, scale = "column",
								clustering_distance_rows = 'correlation', clustering_method = 'ward.D2',
								cluster_cols = FALSE, show_rownames = FALSE, silent = TRUE)
			  New_clusters = c(New_clusters, col_clusters[heat_map$tree_row$labels[heat_map$tree_row$order]])
			} else {
			  vec = rep(num, nrow(mat))
			  names(vec) = rownames(mat)
			  New_clusters = c(New_clusters, vec)
			}
		} else {
			vec = rep(num, sum(seurat_clusters==num))
			names(vec) = names(which(col_clusters==num))
			New_clusters = c(New_clusters, vec)
		}
	  }
	  col_clusters = New_clusters
	  rm(New_clusters,num,mat,heat_map)
	  counts = t(counts[names(col_clusters),])
	  counts_scaled = scale(counts)
  } else {
	  counts = counts[feat2use$Feature,]
	  counts_scaled = scale(counts)
	  counts_scaled = as.matrix(data.frame(t(counts_scaled), row.names=feat2use$Feature))
	  colnames(counts_scaled) <- names(counts)
  }
  
  
  # ComplexHeatmap
  col_fun = circlize::colorRamp2(c(-5,0,5), c('yellow','black','purple'))
  heat_map = ComplexHeatmap::Heatmap(matrix = counts_scaled, cluster_rows = FALSE, cluster_columns = FALSE,
                                   col = col_fun, name = "Heatmap", show_column_names = FALSE, row_title_rot = 0,
                                   row_names_side = "right", row_names_gp = gpar(fontsize = 12), border = TRUE, show_heatmap_legend = opt$Show_Heatmap_Legend,
                                   column_split = col_clusters, row_split = row_clusters, row_gap = grid::unit(3, "mm"), column_gap = grid::unit(4, "mm"),
                                   top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 1:length(unique(col_clusters))))),
                                   left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 1:length(unique(row_clusters))))))
  pdf(file = paste(opt$outDir,opt$Sample,"/",src,"/",opt$Sample,"_features_heatmap.pdf", sep = ""), width = 25, height = 15)
  draw(heat_map)
  dev.off()
  
  # Small Clusters only
  small_clusters = names(which(table(col_clusters)<300))
  col_smallclusters = col_clusters[which(col_clusters %in% small_clusters)]
  counts_scaled_smallclusters = counts_scaled[,names(col_smallclusters)]
  if (length(feat2use$Feature)==1)
      counts_scaled_smallclusters = as.matrix(data.frame(t(counts_scaled_smallclusters), row.names=feat2use$Feature))
  heat_map_small = ComplexHeatmap::Heatmap(matrix = counts_scaled_smallclusters, cluster_rows = FALSE, cluster_columns = FALSE,
                                   col = col_fun, name = "Heatmap", show_column_names = FALSE, row_title_rot = 0,
                                   row_names_side = "right", row_names_gp = gpar(fontsize = 12), border = TRUE, show_heatmap_legend = opt$Show_Heatmap_Legend,
                                   column_split = col_smallclusters, row_split = row_clusters, row_gap = grid::unit(3, "mm"), column_gap = grid::unit(4, "mm"),
                                   top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 1:length(unique(col_smallclusters))))),
                                   left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 1:length(unique(row_clusters))))))
  pdf(file = paste(opt$outDir,opt$Sample,"/",src,"/",opt$Sample,"_features_heatmap_small_clusters.pdf", sep = ""), width = 25, height = 15)
  draw(heat_map_small)
  dev.off()
  rm('heat_map')
  }
}

# Plot FeaturePlot
if (opt$FeaturePlot) {
  writeLog(logfile, paste("Plotting FeaturePlots...",sep=''))
  if (is.na(opt$Features_file)) {
    writeLog(logfile, paste("ERROR: Please provide cell-type marker genes list [--Features_file]"))
	stop()
  }
  for (src in names(features_list)) {
    featureData = features_list[[src]]
    if (!dir.exists(paste(opt$outDir,opt$Sample,"/",src,sep='')))
      dir.create(paste(opt$outDir,opt$Sample,"/",src,sep=''))
	if (!dir.exists(paste(opt$outDir,opt$Sample,"/",src,"/FeaturePlots/",sep='')))
      dir.create(paste(opt$outDir,opt$Sample,"/",src,"/FeaturePlots/",sep=''))
    celltypes <- unique(featureData$CellType)
    for (ct in celltypes) {
      feat2use = featureData[featureData$CellType==ct,'Feature']
      feat2use = feat2use[feat2use %in% rownames(obj_seurat)]
      if (length(feat2use)>0) {
        ct_plt <- FeaturePlot(obj_seurat, features = feat2use)
        pdf(file = paste(opt$outDir,opt$Sample,"/",src,"/FeaturePlots/",opt$Sample,'_',stringr::str_replace(ct,' ','_'),'_FeaturePlot.pdf',sep=''), 
             width = 25, height = 15)
        print(ct_plt)
        dev.off()
      } else {
        writeLog(logfile, paste("ATTENTION: None of the requested features were found for ",ct,sep=''))
      }
      rm('ct','ct_plt','feat2use')
    }
  }
}

# Plot ViolinPlot
if (opt$ViolinPlot) {
  writeLog(logfile, paste("Plotting ViolinPlots...",sep=''))
  if (is.na(opt$Features_file)) {
    writeLog(logfile, paste("ERROR: Please provide cell-type marker genes list [--Features_file]"))
	stop()
  }
  for (src in names(features_list)) {
    featureData = features_list[[src]]
    if (!dir.exists(paste(opt$outDir,opt$Sample,"/",src,sep='')))
      dir.create(paste(opt$outDir,opt$Sample,"/",src,sep=''))
	if (!dir.exists(paste(opt$outDir,opt$Sample,"/",src,"/ViolinPlots/",sep='')))
      dir.create(paste(opt$outDir,opt$Sample,"/",src,"/ViolinPlots/",sep=''))
    celltypes <- unique(featureData$CellType)
    for (ct in celltypes) {
      feat2use = featureData[featureData$CellType==ct,'Feature']
      feat2use = feat2use[feat2use %in% rownames(obj_seurat)]
      if (length(feat2use)>0) {
        ct_plt <- VlnPlot(obj_seurat, features = feat2use)
        pdf(paste(opt$outDir,opt$Sample,"/",src,"/ViolinPlots/",opt$Sample,'_',stringr::str_replace(ct,' ','_'),'_ViolinPlot.pdf',sep=''), 
             width = 25, height = 15)
        print(ct_plt)
        dev.off()
      } else {
        writeLog(logfile, paste("ATTENTION: None of the requested features were found for ",ct,sep=''))
      }
      rm('ct','ct_plt','feat2use')
    }
    All_VLN_plt <- VlnPlot(obj_seurat, features = featureData$Feature,sort=T,stack=T,same.y.lims=F,flip=T)
    pdf(paste(opt$outDir,opt$Sample,"/",src,"/ViolinPlots/",opt$Sample,'_Stacked_ViolinPlot.pdf',sep=''), 
              width = 25, height = 15)
    print(All_VLN_plt)
    dev.off()
  }
}

# Plot DotPlot
if (opt$DotPlot) {
  writeLog(logfile, paste("Plotting DotPlots...",sep=''))
  if (is.na(opt$Features_file)) {
    writeLog(logfile, paste("ERROR: Please provide cell-type marker genes list [--Features_file]"))
	stop()
  }
  for (src in names(features_list)) {
    featureData = features_list[[src]]
    if (!dir.exists(paste(opt$outDir,opt$Sample,"/",src,sep='')))
      dir.create(paste(opt$outDir,opt$Sample,"/",src,sep=''))
	if (!dir.exists(paste(opt$outDir,opt$Sample,"/",src,"/DotPlots/",sep='')))
      dir.create(paste(opt$outDir,opt$Sample,"/",src,"/DotPlots/",sep=''))
    celltypes <- unique(featureData$CellType)
    for (ct in celltypes) {
      feat2use = featureData[featureData$CellType==ct,'Feature']
      feat2use = feat2use[feat2use %in% rownames(obj_seurat)]
      if (length(feat2use)>0) {
        ct_plt <- DotPlot(obj_seurat, features = unique(feat2use))
        pdf(paste(opt$outDir,opt$Sample,"/",src,"/DotPlots/",opt$Sample,'_',stringr::str_replace(ct,' ','_'),'_DotPlot.pdf',sep=''), 
             width = 25, height = 15)
        print(ct_plt)
        dev.off()
      } else {
        writeLog(logfile, paste("ATTENTION: None of the requested features were found for ",ct,sep=''))
      }
      rm('ct','ct_plt','feat2use')
    }
  }
}

# StackedViolinPlot functions
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}
# Plot StackedViolinPlot
if (opt$StackedViolinPlot) {
  writeLog(logfile, paste("Plotting Stacked Violin-Plots...",sep=''))
  if (is.na(opt$Features_file)) {
    writeLog(logfile, paste("ERROR: Please provide cell-type marker genes list [--Features_file]"))
	stop()
  }
  for (src in names(features_list)) {
    featureData = features_list[[src]]
    if (!dir.exists(paste(opt$outDir,opt$Sample,"/",src,sep='')))
      dir.create(paste(opt$outDir,opt$Sample,"/",src,sep=''))
	feat2use = featureData$Feature[featureData$Feature %in% rownames(obj_seurat)]
	vlnplt = StackedVlnPlot(obj = obj_seurat, features = feat2use)
	pdf(file = paste(opt$outDir,opt$Sample,"/",src,"/",opt$Sample,"_StackedViolinPlot.pdf", sep = ""), width = 25, height = 15)
	print(vlnplt)
	dev.off()
  }
}

if (opt$RUNUMAP){
    writeLog(logfile, paste("Running UMAP..."))
	# Reduction usage
	if (is.na(opt$reduction)) {
	  writeLog(logfile, paste("Using 'pca' dimensional reduction as input for FindNeighbors and RunUMAP",sep=''))
	  use_reducs = 'pca'
	} else {
	  if (opt$reduction %in% names(obj_seurat@reductions)) {
		writeLog(logfile, paste("Using '",opt$reduction,"' dimensional reduction as input for FindNeighbors and RunUMAP",sep=''))
		use_reducs = opt$reduction
	  } else {
		writeLog(logfile, paste("ERROR: Unable to find reduction '",opt$reduction,"'. Available reductions: ",paste(names(obj_seurat@reductions),collapse=','),sep=''))
		stop()
	  }
	}
	obj_seurat <- RunUMAP(obj_seurat, reduction = use_reducs, dims = 1:opt$dims)
	opt$UMAP = TRUE
	
}


if (opt$UMAP) {
	  writeLog(logfile, paste("Plotting UMAP...",sep=''))
		# Reduction usage
		if ("umap" %in% names(obj_seurat@reductions)) {
			writeLog(logfile, paste("Using '",opt$reduction,"' dimensional reduction as input for FindNeighbors and RUNTSNE",sep=''))
			use_reducs = "umap"
		  } else {
			writeLog(logfile, paste("ERROR: Unable to find reduction '",opt$reduction,"'. Available reductions: ",paste(names(obj_seurat@reductions),collapse=','),sep=''))
			stop()
		  }
	# Plots
	# legend_pos = "none"
	# if (opt$Show_Legend)
	legend_pos = "right"
	writeLog(logfile, paste("Plotting reduction figures..."))
	plt1 <- DimPlot(obj_seurat, reduction = use_reducs, label = T, pt.size = 1.5, label.size = 8)+
	  theme(legend.position = legend_pos)
	pdf(paste(opt$outDir,opt$Sample,"_",use_reducs,"ByCluster.pdf", sep = ""),width = 20, height = 15)
	print(plt1)
	dev.off()
}




if (opt$RUNTSNE){
    writeLog(logfile, paste("Running tSNE..."))
	# Reduction usage
	if (is.na(opt$reduction)) {
	  writeLog(logfile, paste("Using 'pca' dimensional reduction as input for FindNeighbors and RUNTSNE",sep=''))
	  use_reducs = 'pca'
	} else {
	  if (opt$reduction %in% names(obj_seurat@reductions)) {
		writeLog(logfile, paste("Using '",opt$reduction,"' dimensional reduction as input for FindNeighbors and RUNTSNE",sep=''))
		use_reducs = opt$reduction
	  } else {
		writeLog(logfile, paste("ERROR: Unable to find reduction '",opt$reduction,"'. Available reductions: ",paste(names(obj_seurat@reductions),collapse=','),sep=''))
		stop()
	  }
	}
	obj_seurat <- RunTSNE(obj_seurat, reduction = use_reducs, dims = 1:opt$dims)
	opt$TSNE = TRUE
}


if (opt$TSNE) {
	  writeLog(logfile, paste("Plotting TSNE...",sep=''))
		# Reduction usage
		if ("tsne" %in% names(obj_seurat@reductions)) {
			writeLog(logfile, paste("Using '",opt$reduction,"' dimensional reduction as input for FindNeighbors and RUNTSNE",sep=''))
			use_reducs = "tsne"
		  } else {
			writeLog(logfile, paste("ERROR: Unable to find reduction '",opt$reduction,"'. Available reductions: ",paste(names(obj_seurat@reductions),collapse=','),sep=''))
			stop()
		  }
	# Plots
	# legend_pos = "none"
	# if (opt$Show_Legend)
	legend_pos = "right"
	writeLog(logfile, paste("Plotting reduction figures..."))
	plt1 <- DimPlot(obj_seurat, reduction = use_reducs, label = T, pt.size = 1.5, label.size = 8)+
	  theme(legend.position = legend_pos)
	pdf(paste(opt$outDir,opt$Sample,"_",use_reducs,"ByCluster.pdf", sep = ""),width = 20, height = 15)
	print(plt1)
	dev.off()
}




opt$ID = "Test"

# Clusters distribution per Sample
samples = unique(obj_seurat$orig.ident)
samples = samples[order(samples)]
df = data.frame()
for (sample in samples) {
  sample_df = data.frame(table(Idents(obj_seurat)[obj_seurat$orig.ident==sample]))
  colnames(sample_df) = c('Cluster','Cells')
  sample_df$Ratio = sample_df$Cells / length(Idents(obj_seurat)[obj_seurat$orig.ident==sample])
  sample_df$Sample = sample
  df = rbind(df, sample_df)
}
samples_plot = ggplot(df, aes(x=Cluster,y=Ratio,fill=Sample)) +
                      geom_col(position = position_dodge()) +
                      ggtitle("Clusters distribution per Sample")+
                      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
ggplot2::ggsave(filename = paste(opt$outDir,opt$ID,'_ClustersDistPerSample.pdf', sep = ""),
                samples_plot, dpi=600, width=12, height=6,device='pdf')

samples_plot = ggplot(df)+ aes(x=Sample,y=Ratio,fill=Cluster) + 
               geom_bar(stat="identity", color = 'white', size = 0.5)+
               ggtitle("Clusters distribution per Sample")
               

# if (length(unique( df$Cluster))<12){
    # samples_plot = samples_plot + scale_fill_brewer(palette="Paired")
# }

ggplot2::ggsave(filename = paste(opt$outDir,opt$ID,'_PerSample_ClustersDist.pdf', sep = ""),
                samples_plot, dpi=600, width=12, height=6,device='pdf')
write.csv(x = df,file = paste(opt$outDir,opt$ID,'ClustersDist.csv', sep = ""))

# Samples distribution per Cluster
# clusters=as.numeric(as.vector(unique(Idents(obj_seurat))))
clusters=unique(Idents(obj_seurat))
# clusters = clusters[order(clusters)]
df = data.frame()
for (cluster in clusters) {
  # cluster_df = data.frame(table(obj_seurat$orig.ident[obj_seurat$seurat_clusters==cluster]))
  cluster_df = data.frame(table(obj_seurat$orig.ident[Idents(obj_seurat)==cluster]))
  colnames(cluster_df) = c('Sample','Cells')
  cluster_df$Ratio = cluster_df$Cells / length(obj_seurat$orig.ident[Idents(obj_seurat)==cluster])
  cluster_df$Cluster = cluster
  df = rbind(df, cluster_df)
}

clusters_plot = ggplot(df, aes(x=Cluster,y=Ratio,fill=Sample)) +
                       geom_col(position = position_stack()) +
                       geom_hline(yintercept = 0.5, lty = "dashed") +
                       # scale_x_continuous("Cluster", labels = as.character(clusters), breaks = clusters) +
                       ggtitle("Samples distribution per Cluster")+
                       theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
ggplot2::ggsave(filename = paste(opt$outDir,opt$ID,'_SamplesDistPerCluster.pdf', sep = ""),
                clusters_plot, dpi=600, width=12, height=6,device='pdf')








writeLog(logfile, paste("Finished",sep=''))