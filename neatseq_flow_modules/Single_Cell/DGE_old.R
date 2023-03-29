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
  make_option(c("-s", "--Sample"), type="character", default = NA,
              help="Sample's name", metavar = "character"),
  make_option(c("--test_use"), type="character", default = 'wilcox',
              help="Denotes which test to use (Default is wilcox)", metavar = "character"),
  make_option(c("--logfc_threshold"), type="numeric", default = 0.25,
              help="Limit testing to genes which show, on average, at least X-fold difference (log-scale). Default is 0.25", metavar = "character"),
  make_option(c("--min_pct1"), type="numeric", default = 0.9,
              help="only test genes that are detected in a minimum fraction of min_pct1 cells in population 1 (Default is 0.9)", metavar = "character"),
  make_option(c("--min_pct2"), type="numeric", default = 0.1,
              help="only test genes that are detected in a minimum fraction of min_pct2 cells in population 2 (Default is 0.1)", metavar = "character"),
  make_option(c("--only_pos"), action="store_true", default = FALSE,
              help="Only return positive markers (Default is False)", metavar = "character"),
  make_option(c("--Human"), action="store_true", default = FALSE,
              help="Use Human Genome wide annotation", metavar = "character"),
  make_option(c("--Mouse"), action="store_true", default = FALSE,
              help="Use Mouse Genome wide annotation", metavar = "character"),
  make_option(c("--No_DGE_vs_all"), action="store_true", default = FALSE,
              help="Skip differential gene expression analysis for clusters vs all other cells", metavar = "character"),
  make_option(c("--pAdjustMethod"), type="character", default = 'fdr',
              help="one of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none' (Default is fdr)", metavar = "character"),
  make_option(c("--ont"), type="character", default = 'BP',
              help="One of 'BP', 'MF', and 'CC' subontologies, or 'ALL' for all three. (Default is BP)", metavar = "character"),
  make_option(c("--simplify_cutoff"), type="numeric", default = 0.6,
              help="Simplify similarity cutoff (Default is 0.6)", metavar = "character"),
  make_option(c("--ALPHA"), type="numeric", default = 0.05,
              help="Significant Level Cutoff (Default is 0.05)", metavar = "character"),
  make_option(c("--CPUs"), type="numeric", default = 1,
              help="Number of threads to use (Default is 1)", metavar = "character"),
  make_option(c("--Memory"), type="numeric", default = NA,
              help="Allocate memory for Seurat parallelization in GB (Default is NA)", metavar = "character"),
  make_option(c("--Subsample"), type="character", default = NA,
              help="Subsample number of cells per cluster in DoHeatmap (Default is all cells)", metavar = "character"),
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
cat(paste('[',Sys.time(),']: Differential Gene Expression: ',opt$Sample,sep=''), file=logfile, sep='\n')
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
if (!is.na(opt$Subsample)) {
  subsample = as.numeric(opt$Subsample)
  writeLog(logfile, paste("Using ",subsample," cells per cluster in DoHeatmap",sep=''))
}
if (opt$Human) {
  writeLog(logfile, paste("Loading Genome wide annotation for Human"))
  library(org.Hs.eg.db)
  OrgDb = org.Hs.eg.db
  opt$organism = 'hsa'
} else if (opt$Mouse) {
  writeLog(logfile, paste("Loading Genome wide annotation for Mouse"))
  library(org.Mm.eg.db)
  OrgDb = org.Mm.eg.db
  opt$organism = 'mmu'
} else {
  writeLog(logfile, paste("Genome wide annotation OrgDb object must be specified [--Human or --Mouse]"))
  stop
}

# Parallelization
if (opt$CPUs>1) {
  writeLog(logfile, paste("Using multicore parallelization with ",opt$CPUs,' threads',sep=''))
  plan("multicore", workers = opt$CPUs)
  if (!is.na(opt$Memory)) {
    writeLog(logfile, paste("Allocating ",as.numeric(opt$Memory),'GB of memory',sep=''))
    options(future.globals.maxSize = as.numeric(opt$Memory) * 1024 ^ 3)
  }
}

# Prepare for DGE analysis
clusters = unique(obj_seurat$seurat_clusters) %>% as.vector()
counts = GetAssayData(object = obj_seurat, slot = "counts")
background.genes = rownames(counts)[apply(counts>0,MARGIN = 1, FUN=any)] %>% as.character()
background_conv <- clusterProfiler::bitr(geneID=background.genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
background.genes_entrezid  <- background_conv$ENTREZID
if (!dir.exists(paste(opt$outDir,opt$Sample,sep='')))
  dir.create(paste(opt$outDir,opt$Sample,sep=''))
if (!opt$No_DGE_vs_all)
  if (!dir.exists(paste(opt$outDir,opt$Sample,'/FindMarkers_vs_all',sep='')))
    dir.create(paste(opt$outDir,opt$Sample,'/FindMarkers_vs_all',sep=''))
if (!dir.exists(paste(opt$outDir,opt$Sample,'/FindMarkers_vs_each_other',sep='')))
  dir.create(paste(opt$outDir,opt$Sample,'/FindMarkers_vs_each_other',sep=''))

# Prepare variables for heatmaps
# norm_counts = as.matrix(GetAssayData(object = obj_seurat, slot = "data"))
# #norm_counts = scale(norm_counts)
# no_var=names(which(apply(norm_counts,MARGIN = 2,FUN = sd)==0))
# norm_counts = norm_counts[,!(colnames(norm_counts) %in% no_var)]
# seurat_clusters = obj_seurat$seurat_clusters[!(names(obj_seurat$seurat_clusters) %in% no_var)]
# # Determine column clusters
# col_clusters = obj_seurat$seurat_clusters %>% as.vector() %>% as.numeric()
# names(col_clusters) = colnames(obj_seurat)
# col_clusters = col_clusters[order(col_clusters)]
# col_clusters = col_clusters[!(names(col_clusters) %in% no_var)]

### DGE analysis for clusters vs all other cells ###
if (!opt$No_DGE_vs_all) {
  writeLog(logfile, paste("Running DGE Analysis: Cluster vs. all"))
  DGE_clusters_vs_all <- list()
  sig_genes = list()
  for (cluster in clusters) {
    writeLog(logfile, paste("Started Cluster ",cluster,sep=''))
    file = paste(opt$outDir,opt$Sample,'/FindMarkers_vs_all/DGE_Cluster',cluster,'.tab',sep='')
    if (file.exists(file)) {
      writeLog(logfile, paste("Loaded FindMarkers results for cluster ",cluster,sep=''))
      DGE_clusters_vs_all[[cluster]] <- read.table(file, header = TRUE, row.names = 1)
    } else {
      DGE_clusters_vs_all[[cluster]] <- Seurat::FindMarkers(obj_seurat, ident.1 = cluster, test.use = opt$test_use,
                                                            logfc.threshold = opt$logfc_threshold, min.pct = 0,
                                                            only.pos = opt$only_pos, verbose = FALSE)
      writeLog(logfile, paste("Finished FindMarkers for cluster ",cluster,sep=''))
      write.table(DGE_clusters_vs_all[[cluster]], file = file, sep = '\t')
    }
    df <- DGE_clusters_vs_all[[cluster]]
    df = df[df$pct.1>=opt$min_pct1 & df$pct.2>=opt$min_pct2 & df$p_val_adj<opt$ALPHA,]
    writeLog(logfile, paste("Found ",dim(df)[1],' differentially expressed features with min_pct1 >= ',opt$min_pct1,' and min_pct2 >= ',opt$min_pct2,sep=''))
    if (dim(df)[1]>0) {
      # heatmap_counts = norm_counts[rownames(df),]
      # heatmap_counts = t(heatmap_counts)
      # New_clusters = c()
      # for (num in unique(col_clusters)) {
      #   mat = heatmap_counts[seurat_clusters==num,]
      #   mat = mat[,which(apply(mat,MARGIN = 2,sd)!=0)]
      #   heat_map = pheatmap::pheatmap(mat = mat, cluster_rows = TRUE, show_colnames = FALSE, scale = "column",
      #                                 clustering_distance_rows = 'correlation', clustering_method = 'ward.D2',
      #                                 cluster_cols = FALSE, show_rownames = FALSE, silent = TRUE)
      #   New_clusters = c(New_clusters, col_clusters[heat_map$tree_row$labels[heat_map$tree_row$order]])
      # }
      # col_clusters = New_clusters
      # rm(New_clusters,num,mat,heat_map)
      # heatmap_counts = t(heatmap_counts[names(col_clusters),])
      # no_var=names(which(apply(heatmap_counts,MARGIN = 2,FUN = sd)==0))
      # heatmap_counts = heatmap_counts[,!(colnames(heatmap_counts) %in% no_var)]
      # heatmap_counts = scale(heatmap_counts)
      # cols_cluster = col_clusters[!(names(col_clusters) %in% no_var)]
      # 
      # 
      # # Draw Heatmap
      # col_fun = circlize::colorRamp2(c(-5,0,5), c('blue','white','red'))
      # heat_map = ComplexHeatmap::Heatmap(matrix = heatmap_counts, cluster_rows = FALSE, cluster_columns = FALSE,
      #                                    col = col_fun, name = "Heatmap", show_column_names = FALSE, row_title_rot = 0,
      #                                    row_names_side = "right", row_names_gp = grid::gpar(fontsize = 12), border = TRUE, show_heatmap_legend = FALSE,
      #                                    column_split = cols_cluster, row_gap = grid::unit(3, "mm"), column_gap = grid::unit(4, "mm"),
      #                                    top_annotation = ComplexHeatmap::HeatmapAnnotation(foo = ComplexHeatmap::anno_block(gp = grid::gpar(fill = 1:length(unique(cols_cluster))))))
      # png(paste(opt$outDir,'FindMarkers_vs_all/Cluster',cluster,'_heatmap_vs_all2.jpeg',sep=''), width = 2000, height = 1000)
      # ComplexHeatmap::draw(heat_map)
      # dev.off()
      
      sig_genes[[cluster]] <- rownames(df)
	  if (!is.na(opt$Subsample)) {
	    heat_map = DoHeatmap(subset(obj_seurat, downsample=subsample), features = rownames(df)) + NoLegend()
	  } else {
	    heat_map = DoHeatmap(obj_seurat, features = rownames(df)) + NoLegend()
	  }
      jpeg(paste(opt$outDir,opt$Sample,'/FindMarkers_vs_all/Cluster',cluster,'_heatmap_vs_all.jpeg',sep=''),width = 3000, height = 1500)
      print(heat_map)
      dev.off()
    }
  }
  openxlsx::write.xlsx(DGE_clusters_vs_all, 
                       file=paste(opt$outDir,opt$Sample,'/FindMarkers_vs_all/DGE_FindMarkers_vs_all.xlsx',sep=''), 
                       row.names = TRUE)
  
  # GO Enrichment
  enrichGO_compClusters <- try(clusterProfiler::compareCluster(geneClusters = sig_genes, fun = "enrichGO", pAdjustMethod = opt$pAdjustMethod,
                                                               OrgDb = OrgDb, ont = opt$ont, universe = background.genes, keyType = "SYMBOL"),
                               silent = TRUE)
  if (!inherits(enrichGO_compClusters,"try-error")) {
    enrichGO_compClusters@compareClusterResult$GeneRatio=paste("'",enrichGO_compClusters@compareClusterResult$GeneRatio,sep='')
    write.csv(enrichGO_compClusters@compareClusterResult,
              file=paste(opt$outDir,opt$Sample,'/',opt$Sample,"_GO_compareClusters_vs_all_others.csv",sep=''),quote = FALSE,row.names = TRUE)
    enrichGO_simplify <- try(clusterProfiler::simplify(enrichGO_compClusters, cutoff = opt$simplify_cutoff),
                             silent = TRUE)
    if (!inherits(enrichGO_simplify,"try-error")) {
      enrichGO_simplify@compareClusterResult$GeneRatio=paste("'",enrichGO_simplify@compareClusterResult$GeneRatio,sep='')
      write.csv(enrichGO_simplify@compareClusterResult,
                file=paste(opt$outDir,opt$Sample,'/',opt$Sample,"_GO_simplify_vs_all_others.csv",sep=''),quote = FALSE,row.names = TRUE)
      
      if (dim(enrichGO_simplify@compareClusterResult)[1]<250) {
        fontsize=9
      } else {
        fontsize=4
      }
      enrichGO_simplify@compareClusterResult$GeneRatio=stringr::str_remove_all(enrichGO_simplify@compareClusterResult$GeneRatio,"'")
      enrichGO_dotplot <- enrichplot::dotplot(enrichGO_simplify, showCategory=10000, font.size=fontsize, label_format=100)
      ggplot2::ggsave(filename = paste(opt$outDir,opt$Sample,'/',opt$Sample,"_GO_simplify_vs_all_others_dotplot.pdf",sep=''),
                      enrichGO_dotplot, dpi = 600,device = "pdf",width = 30,height = 50, limitsize = FALSE)
    }
  }
  
  # Convert sig.genes to entrez gene ids
  sig_genes_entrez = list()
  for (cluster in names(sig_genes)) {
    conv = clusterProfiler::bitr(geneID=sig_genes[[cluster]], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
    sig_genes_entrez[[cluster]] = conv$ENTREZID
    rm(conv)
  }
  
  # KEGG Enrichment
  enrichKEGG_compClusters <- try(clusterProfiler::compareCluster(geneClusters = sig_genes_entrez, fun = "enrichKEGG", pAdjustMethod = opt$pAdjustMethod,
                                                                 organism = opt$organism, universe = background.genes_entrezid),
                                 silent = TRUE)
  if (!inherits(enrichKEGG_compClusters,"try-error")) {
    enrichKEGG_compClusters@compareClusterResult$GeneRatio=paste("'",enrichKEGG_compClusters@compareClusterResult$GeneRatio,sep='')
    write.csv(enrichKEGG_compClusters@compareClusterResult,
              file=paste(opt$outDir,opt$Sample,'/',opt$Sample,"_KEGG_compareClusters_vs_all_others.csv",sep=''),quote = FALSE,row.names = TRUE)
    if (dim(enrichKEGG_compClusters@compareClusterResult)[1]<250) {
      fontsize=9
    } else {
      fontsize=4
    }
    enrichKEGG_compClusters@compareClusterResult$GeneRatio=stringr::str_remove_all(enrichKEGG_compClusters@compareClusterResult$GeneRatio,"'")
    enrichKEGG_dotplot <- enrichplot::dotplot(enrichKEGG_compClusters, showCategory=10000, font.size=fontsize, label_format=100)
    ggplot2::ggsave(filename = paste(opt$outDir,opt$Sample,'/',opt$Sample,"_KEGG_compareClusters_vs_all_others.pdf",sep=''),
                    enrichKEGG_dotplot, dpi = 600,device = "pdf",width = 30,height = 50, limitsize = FALSE)
    rm(sig_genes,cluster,file,df,sig_genes_entrez,enrichGO_compClusters,enrichGO_simplify,enrichGO_dotplot,enrichKEGG_compClusters,
       enrichKEGG_dotplot)
  }
}



### DGE analysis for clusters vs each other ###
for (cluster1 in clusters) {
  print(paste("Running contrasts for cluster ",cluster1,sep=''))
  DGE_clusters_vs_clusters <- list()
  sig_genes = list()
  if (!dir.exists(paste(opt$outDir,opt$Sample,'/Cluster_',cluster1,sep='')))
    dir.create(paste(opt$outDir,opt$Sample,'/Cluster_',cluster1,sep=''))
  clusters2 = clusters[clusters != cluster1]
  for (cluster2 in clusters2) {
    comb <- paste(cluster1,cluster2,sep='vs')
    file = paste(opt$outDir,opt$Sample,'/FindMarkers_vs_each_other/DGE_',comb,'.tab',sep='')
    if (file.exists(file)) {
      print(paste("Loaded FindMarkers results for contrast ",comb,sep=''))
      DGE_clusters_vs_clusters[[comb]] <- read.table(paste(opt$outDir,opt$Sample,'/FindMarkers_vs_each_other/DGE_',comb,'.tab',sep=''), 
                                                     header = TRUE, row.names = 1)
    } else {
      DGE_clusters_vs_clusters[[comb]] <- Seurat::FindMarkers(obj_seurat, ident.1 = cluster1, ident.2 = cluster2,
                                                              test.use = opt$test_use, min.pct = 0,
                                                              logfc.threshold = opt$logfc_threshold,
                                                              only.pos = opt$only_pos, verbose = FALSE)
      print(paste("Finished FindMarkers for contrast ",comb,sep=''))
      write.table(DGE_clusters_vs_clusters[[comb]], file = file, sep = '\t')
    }
    df <- DGE_clusters_vs_clusters[[comb]]
    df = df[df$pct.1>=opt$min_pct1 & df$pct.2>=opt$min_pct2 & df$p_val_adj<opt$ALPHA,]
    print(paste("Found ",dim(df)[1],' differentially expressed features with min_pct1 >= ',opt$min_pct1,' and min_pct2 >= ',opt$min_pct2,sep=''))
    if (dim(df)[1]>0) {
      sig_genes[[comb]] <- rownames(df)
    }
  }
  if (!is.na(opt$Subsample)) {
    heat_map = DoHeatmap(subset(obj_seurat, downsample=subsample), features = unique(unlist(sig_genes))) + NoLegend()
  } else {
    heat_map = DoHeatmap(obj_seurat, features = unique(unlist(sig_genes))) + NoLegend()
  }
  jpeg(paste(opt$outDir,opt$Sample,'/Cluster_',cluster1,'/Cluster',cluster1,'_heatmap_vs_each_other.jpeg',sep=''),width = 3000, height = 1500)
  print(heat_map)
  dev.off()
  
  openxlsx::write.xlsx(DGE_clusters_vs_clusters,
                       file=paste(opt$outDir,opt$Sample,'/Cluster_',cluster1,'/DGE_FindMarkers_Cluster',cluster1,'_vs_each_other.xlsx',sep=''),
                       row.names = TRUE)
  
  # GO Enrichment
  enrichGO_compClusters <- try(clusterProfiler::compareCluster(geneClusters = sig_genes, fun = "enrichGO", pAdjustMethod = opt$pAdjustMethod,
                                                               OrgDb = OrgDb, ont = opt$ont, universe = background.genes, keyType = "SYMBOL"),
                               silent = TRUE)
  if (!inherits(enrichGO_compClusters,"try-error")) {
    enrichGO_compClusters@compareClusterResult$GeneRatio=paste("'",enrichGO_compClusters@compareClusterResult$GeneRatio,sep='')
    write.csv(enrichGO_compClusters@compareClusterResult,
              file=paste(opt$outDir,opt$Sample,'/Cluster_',cluster1,'/Cluster',cluster1,'_GO_compareClusters_vs_each_other.csv',sep=''),quote = FALSE,row.names = TRUE)
    enrichGO_simplify <- try(clusterProfiler::simplify(enrichGO_compClusters, cutoff = opt$simplify_cutoff),
                             silent = TRUE)
    if (!inherits(enrichGO_simplify,"try-error")) {
      enrichGO_simplify@compareClusterResult$GeneRatio=paste("'",enrichGO_simplify@compareClusterResult$GeneRatio,sep='')
      write.csv(enrichGO_simplify@compareClusterResult,
                file=paste(opt$outDir,opt$Sample,'/Cluster_',cluster1,'/Cluster',cluster1,'_GO_simplify_vs_each_other.csv',sep=''),quote = FALSE,row.names = TRUE)
      
      if (dim(enrichGO_simplify@compareClusterResult)[1]<250) {
        fontsize=9
      } else {
        fontsize=4
      }
      enrichGO_simplify@compareClusterResult$GeneRatio=stringr::str_remove_all(enrichGO_simplify@compareClusterResult$GeneRatio,"'")
      enrichGO_dotplot <- enrichplot::dotplot(enrichGO_simplify, showCategory=10000, font.size=fontsize, label_format=100)
      ggplot2::ggsave(filename = paste(opt$outDir,opt$Sample,'/Cluster_',cluster1,'/Cluster',cluster1,'_GO_simplify_vs_each_other_dotplot.pdf',sep=''),
                      enrichGO_dotplot, dpi = 600,device = "pdf",width = 30,height = 50, limitsize = FALSE)
    }
  }
  
  # Convert sig.genes to entrez gene ids
  sig_genes_entrez = list()
  for (contrast in names(sig_genes)) {
    conv = try(clusterProfiler::bitr(geneID=sig_genes[[contrast]], fromType = "SYMBOL", 
                                     toType = "ENTREZID", OrgDb = OrgDb),
               silent = TRUE)
    if (!inherits(conv,"try-error")) {
      sig_genes_entrez[[contrast]] = conv$ENTREZID
      rm(conv)
    }
  }
  
  # KEGG Enrichment
  enrichKEGG_compClusters <- try(clusterProfiler::compareCluster(geneClusters = sig_genes_entrez, fun = "enrichKEGG", pAdjustMethod = opt$pAdjustMethod,
                                                                 organism = opt$organism, universe = background.genes_entrezid),
                                 silent = TRUE)
  if (!inherits(enrichKEGG_compClusters,"try-error")) {
    enrichKEGG_compClusters@compareClusterResult$GeneRatio=paste("'",enrichKEGG_compClusters@compareClusterResult$GeneRatio,sep='')
    write.csv(enrichKEGG_compClusters@compareClusterResult,
              file=paste(opt$outDir,opt$Sample,'/Cluster_',cluster1,'/Cluster',cluster1,'_KEGG_compareClusters_vs_each_other.csv',sep=''),quote = FALSE,row.names = TRUE)
    if (dim(enrichKEGG_compClusters@compareClusterResult)[1]<250) {
      fontsize=9
    } else {
      fontsize=4
    }
    enrichKEGG_compClusters@compareClusterResult$GeneRatio=stringr::str_remove_all(enrichKEGG_compClusters@compareClusterResult$GeneRatio,"'")
    enrichKEGG_dotplot <- enrichplot::dotplot(enrichKEGG_compClusters, showCategory=10000, font.size=fontsize, label_format=100)
    ggplot2::ggsave(filename = paste(opt$outDir,opt$Sample,'/Cluster_',cluster1,'/Cluster',cluster1,'_KEGG_compareClusters_vs_each_other_dotplot.pdf',sep=''),
                    enrichKEGG_dotplot, dpi = 600,device = "pdf",width = 30,height = 50, limitsize = FALSE)
  }
}

