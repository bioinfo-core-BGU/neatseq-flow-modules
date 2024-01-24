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
  make_option(c("--Treatment_file"), type="character", default = NA,
              help="Path to file with location of the table showing treatment per sample", metavar = "character"),
  make_option(c("-A", "--ident1"), type="character", default = NA,
              help="identity 1 to use in a Pairwise differential gene expression analysis", metavar = "character"),
  make_option(c("-B", "--ident2"), type="character", default = NA,
              help="identity 2 to use in a Pairwise differential gene expression analysis", metavar = "character"),
  make_option(c("--test_use"), type="character", default = 'wilcox',
              help="Denotes which test to use (Default is wilcox)", metavar = "character"),
  make_option(c("--logfc_threshold"), type="numeric", default = 0.25,
              help="Limit testing to genes which show, on average, at least X-fold difference (log-scale). Default is 0.25", metavar = "character"),
  make_option(c("--min_pct"), type="numeric", default = 0.1,
              help="only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations  (Default is 0.1)", metavar = "character"),
  make_option(c("--min_pct1"), type="numeric", default = 0.9,
              help="only test genes that are detected in a minimum fraction of min_pct1 cells in population 1 (Default is 0.9)", metavar = "character"),
  make_option(c("--max_pct2"), type="numeric", default = 0.5,
              help="only test genes that are detected in a maximum fraction of max_pct2 cells in population 2 (Default is 0.5)", metavar = "character"),
  make_option(c("--Not_only_pos"), action="store_true", default = FALSE,
              help="Don't only return positive markers (Default is to return only positive markers)", metavar = "character"),
  make_option(c("--Human"), action="store_true", default = FALSE,
              help="Use Human Genome wide annotation", metavar = "character"),
  make_option(c("--Mouse"), action="store_true", default = FALSE,
              help="Use Mouse Genome wide annotation", metavar = "character"),
  make_option(c("--Force"), action="store_true", default = TRUE,
              help="Force running DGE analysis while over-writing previous results (Default is False)", metavar = "character"),
  make_option(c("--DGE_Pairwise"), action="store_true", default = FALSE,
              help="Run pair-wise differential gene expression analysis of clusters vs each other (Default is False)", metavar = "character"),
  make_option(c("--DGE_Clusters_vs_all"), action="store_true", default = TRUE,
              help="Run Clusters vs all other differential gene expression analysis of clusters vs each other (Default is False)", metavar = "character"),
  make_option(c("--DGE_Within_Clusters"), action="store_true", default = FALSE,
              help="Run Within Clusters differential gene expression analysis [use --ident1 to set the contrast and --DGE_GroupBy to set the source ] (Default is False)", metavar = "character"),
  make_option(c("--DGE_GroupBy"), type="character", default = 'orig.ident',
              help="Group by MetaData col", metavar = "character"),
  make_option(c("--pAdjustMethod"), type="character", default = 'fdr',
              help="one of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none' (Default is fdr)", metavar = "character"),
  make_option(c("--ont"), type="character", default = 'BP',
              help="One of 'BP', 'MF', and 'CC' subontologies, or 'ALL' for all three. (Default is BP)", metavar = "character"),
  make_option(c("--simplify_cutoff"), type="numeric", default = 0.7,
              help="Simplify similarity cutoff (Default is 0.7)", metavar = "character"),
  make_option(c("--Top_n"), type="numeric", default = 10,
              help="How many genes to plot as top marker genes in each cluster (Default is 10)", metavar = "character"),
  make_option(c("--ALPHA"), type="numeric", default = 0.05,
              help="Significant Level Cutoff (Default is 0.05)", metavar = "character"),
  make_option(c("--AllClustersCountCutoff"), type="numeric", default = 0,
              help="All Clusters Count Cutoff, will select genes expressed above the cutoff (Default is counts>0) in percent of cells [--AllClustersCellCutoff] ", metavar = "character"),
  make_option(c("--AllClustersCellCutoff"), type="numeric", default = 0.8,
              help="Percent of Cells Cutoff, will select genes expressed in percent of cells (Default is percent of cells>0.8). Expression cutoff set by --AllClustersCountCutoff ", metavar = "character"),
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


if (!is.na(opt$Treatment_file)) {
 # treatments_list = list()
  treatmentTable <- read.table(opt$Treatment_file, sep = '\t', header = TRUE)
  if (isFALSE(all(colnames(treatmentTable) == c(colnames(treatmentTable)[1],colnames(treatmentTable)[2])))) {
    writeLog(logfile, paste("ERROR: Invalid column names. Please use SampleID and Treatment"))
    stop()
  }
  }

ind= match(obj_seurat@meta.data$orig.ident, treatmentTable$SampleID)
obj_seurat@meta.data$treatment= treatmentTable$Treatment[ind]
#Subset only dent1 and ident2 samples:
ident1 <- opt$ident1
ident2 <- opt$ident2
obj_seurat <- obj_seurat[ ,grepl(ident1, obj_seurat$treatment) | grepl(ident2, obj_seurat$treatment)]


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
clusters = unique(obj_seurat@active.ident) %>% as.vector()
writeLog(logfile, paste("Found ",length(clusters)," Clusters: ",paste(clusters,collapse=','),sep=''))
counts = GetAssayData(object = obj_seurat, slot = "counts", assay = "RNA")
background.genes = rownames(counts)[apply(counts>0,MARGIN = 1, FUN=any)] %>% as.character()
background_conv <- clusterProfiler::bitr(geneID=background.genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
background.genes_entrezid  <- background_conv$ENTREZID
if (!dir.exists(paste(opt$outDir,opt$Sample,sep='')))
  dir.create(paste(opt$outDir,opt$Sample,sep=''))
if (!dir.exists(paste(opt$outDir,opt$Sample,'/FindMarkers_vs_all',sep='')))
    dir.create(paste(opt$outDir,opt$Sample,'/FindMarkers_vs_all',sep=''))
if (opt$DGE_Pairwise)
    if (!dir.exists(paste(opt$outDir,opt$Sample,'/FindMarkers_vs_each_other',sep='')))
        dir.create(paste(opt$outDir,opt$Sample,'/FindMarkers_vs_each_other',sep=''))

if (opt$DGE_Within_Clusters)
    if (!dir.exists(paste(opt$outDir,opt$Sample,'/FindMarkers_Within_Clusters',sep='')))
        dir.create(paste(opt$outDir,opt$Sample,'/FindMarkers_Within_Clusters',sep=''))




### DGE analysis for clusters vs all other cells ###
if (opt$DGE_Clusters_vs_all) {
    writeLog(logfile, paste("Running DGE Analysis: Cluster vs. all"))
    DGE_clusters_vs_all <- list()
    sig_genes = list()

    # FindAllMarkers
    file = paste(opt$outDir,opt$Sample,'/FindMarkers_vs_all/DGE_Markers.csv',sep='')
    if (any(opt$Force,!file.exists(file))) {
      writeLog(logfile, paste("Running FindAllMarkers..."))
      if (obj_seurat@active.assay == "SCT"){
          obj_seurat = PrepSCTFindMarkers(object = obj_seurat)
      }
      DGE_Markers <- FindAllMarkers(obj_seurat, only.pos = !opt$Not_only_pos, min.pct = opt$min_pct , logfc.threshold = opt$logfc_threshold,
                                    test.use = opt$test_use, verbose = FALSE)
      writeLog(logfile, paste("Finished FindAllMarkers, Saving results..."))
      writeLog(logfile, paste("Filtering Genes After DGE Analysis: min.pct1 >= ",opt$min_pct1,"& max.pct.2 <",opt$max_pct2," & ALPHA < ",opt$ALPHA,sep=''))
      if (opt$Not_only_pos) {
        DGE_Markers$Group = 'Up'
        DGE_Markers[DGE_Markers$avg_log2FC<0,'Group'] = 'Down'
        DGE_Markers_UP = DGE_Markers[DGE_Markers$Group == 'Up',]
        DGE_Markers_UP = DGE_Markers_UP[DGE_Markers_UP$pct.1>=opt$min_pct1 & DGE_Markers_UP$pct.2<opt$max_pct2 & DGE_Markers_UP$p_val_adj<opt$ALPHA,]
        
        DGE_Markers_Down = DGE_Markers[DGE_Markers$Group == 'Down',]
        DGE_Markers_Down = DGE_Markers_Down[DGE_Markers_Down$pct.1<opt$max_pct2 & DGE_Markers_Down$pct.2>=opt$min_pct1 & DGE_Markers_Down$p_val_adj<opt$ALPHA,]
        
        DGE_Markers = rbind(DGE_Markers_UP,DGE_Markers_Down)
      }else{
        DGE_Markers = DGE_Markers[DGE_Markers$pct.1>=opt$min_pct1 & DGE_Markers$pct.2<opt$max_pct2 & DGE_Markers$p_val_adj<opt$ALPHA,]
      }
      write.csv(DGE_Markers, file)
    } else {
      writeLog(logfile, paste("Loaded FindAllMarkers results from previous run"))
      DGE_Markers <- read.csv(file, header=T)
    }

    # Top Marker Genes
    writeLog(logfile, paste("Calculating the top ",opt$Top_n," marker genes in each cluster",sep=''))
    if ('Group' %in% colnames(DGE_Markers)) {
      Markers_Up = DGE_Markers[DGE_Markers$Group=='Up',]
      Top_Up = Markers_Up %>% group_by(cluster) %>% slice_head(n=opt$Top_n)
      Markers_Down = DGE_Markers[DGE_Markers$Group=='Down',]
      Top_Down = Markers_Down %>% group_by(cluster) %>% slice_head(n=opt$Top_n)
      DGE_Top_Markers = rbind(Top_Up, Top_Down)
    } else {
      DGE_Top_Markers <- DGE_Markers %>% group_by(cluster) %>% slice_head(n=opt$Top_n)
    }
    write.csv(DGE_Top_Markers, paste(opt$outDir,opt$Sample,'/FindMarkers_vs_all/DGE_Markers_top',opt$Top_n,'.csv',sep=''))
    writeLog(logfile, paste("Plotting top ",opt$Top_n," marker genes",sep=''))
    if (!is.na(opt$Subsample)) {
      heat_map = DoHeatmap(subset(obj_seurat, downsample=subsample), features = DGE_Top_Markers$gene) + NoLegend()
    } else {
      heat_map = DoHeatmap(obj_seurat, features = DGE_Top_Markers$gene) + NoLegend()
    }
    pdf(paste(opt$outDir,opt$Sample,'/FindMarkers_vs_all/DGE_Markers_top',opt$Top_n,'_heatmap.pdf',sep=''),width = 40, height = 15)
    print(heat_map)
    dev.off()

    # Enrichment Analysis
    writeLog(logfile, paste("Setting Genes for Enrichment Analysis: min.pct1 >= ",opt$min_pct1,"& max.pct.2 <",opt$max_pct2," & ALPHA < ",opt$ALPHA,sep=''))
    for (cluster in clusters) {
      df <- DGE_Markers[DGE_Markers$cluster==cluster,]
      if ('Group' %in% colnames(df))
        df = df[df$Group=='Up',]
      df = df[df$pct.1>=opt$min_pct1 & df$pct.2<opt$max_pct2 & df$p_val_adj<opt$ALPHA,]
      writeLog(logfile, paste("Cluster ",cluster,': ',dim(df)[1]," Differentially Expressed Genes",sep=''))
      write.csv(df, paste(opt$outDir,opt$Sample,'/FindMarkers_vs_all/DGE_Cluster',cluster,'.csv',sep=''), row.names=FALSE)
      if (dim(df)[1]>0) {
        sig_genes[[cluster]] <- df$gene
        if (!is.na(opt$Subsample)) {
          heat_map = DoHeatmap(subset(obj_seurat, downsample=subsample), features = df$gene) + NoLegend()
        } else {
          heat_map = DoHeatmap(obj_seurat, features = df$gene) + NoLegend()
        }
        pdf(paste(opt$outDir,opt$Sample,'/FindMarkers_vs_all/Cluster',cluster,'_heatmap_vs_all.pdf',sep=''),width = 40, height = 15)
        print(heat_map)
        dev.off()
      }
    }
    
    counts = GetAssayData(object = obj_seurat, slot = "count", assay = "RNA")
    flag =apply(X = counts,MARGIN=c(1),FUN=function(x) mean(x>opt$AllClustersCountCutoff)>opt$AllClustersCellCutoff)
    AllClustes_sig = row.names(counts[flag,])
    if (length(AllClustes_sig)>0){
        sig_genes[["All_Clustes"]] <- AllClustes_sig
        write.csv(as.data.frame(AllClustes_sig),file = paste(opt$outDir,opt$Sample,'/FindMarkers_vs_all/Genes_expressed_in_more_then',opt$AllClustersCellCutoff,'.csv',sep=''))
    }
    
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


if (opt$DGE_Within_Clusters) {
    writeLog(logfile, paste("Running Within Clusters-DGE Analysis"))
    opt$Not_only_pos = TRUE
    sig_genes = list()
    if (!is.na(opt$ident1)){
        DGE_GroupBy = opt$DGE_GroupBy
        if (DGE_GroupBy %in% colnames(obj_seurat@meta.data)){
            ident1 = opt$ident1
            ident2 = opt$ident2
            idents = c(ident1,ident2)
            if (!all(idents %in% unique(obj_seurat@meta.data[,DGE_GroupBy]))) {
                writeLog(logfile, paste("ERROR: Cannot find the following identities in the object: ",
                                paste(idents[which(!(idents %in% unique(obj_seurat@meta.data[,DGE_GroupBy])))],collapse=','),sep=''))
                stop()
            }
            
            for (cluster in clusters){
                dir_path = paste(opt$outDir,opt$Sample,'/FindMarkers_Within_Clusters/',ident1,'_Within_Cluster_',cluster,'/',sep='')
                if (!dir.exists(dir_path))
                    dir.create(dir_path)
                
                # FindMarkers
                file = paste(dir_path,'DGE_',ident1,"_VS_",ident2,'_Within_Cluster_',cluster,'_Markers.csv',sep='')
                if (any(opt$Force,!file.exists(file))) {
                    writeLog(logfile, paste("Running FindMarkers: ",ident1,'_Within_Cluster_',cluster,sep=''))
                    
                    # DGE_Markers <- FindMarkers(obj_seurat, ident.1 = ident1, subset.ident = cluster ,group.by = DGE_GroupBy, only.pos = !opt$Not_only_pos, min.pct = opt$min_pct, 
                    #                               logfc.threshold = opt$logfc_threshold, test.use = opt$test_use, verbose = FALSE)
                    cluster_obj         = obj_seurat[,obj_seurat@active.ident==cluster]
                    Idents(cluster_obj) = cluster_obj@meta.data[DGE_GroupBy]
                    
                    if (cluster_obj@active.assay == "SCT"){
                        cluster_obj = PrepSCTFindMarkers(object = cluster_obj)
                    }
                    cluster_counts = GetAssayData(object = cluster_obj, slot = "scale.data", assay = "RNA")
                    
                    DGE_Markers <- FindMarkers(cluster_obj, ident.1 = ident1,ident.2 = ident2, only.pos = !opt$Not_only_pos, min.pct = opt$min_pct, 
                                                  logfc.threshold = opt$logfc_threshold, test.use = opt$test_use, verbose = FALSE)
                    writeLog(logfile, paste("Finished FindMarkers, Saving results..."))
                    writeLog(logfile, paste("FindMarkers: Found ", dim(DGE_Markers)[1] ,"Significant Genes Before Filtration" ,ident1,' VS ',ident2," in Cluster ",cluster,sep=''))
                    if (dim(DGE_Markers)[1]>0){
                        DGE_Markers$gene = rownames(DGE_Markers)
                        
                        writeLog(logfile, paste("Filtering Genes After DGE Analysis: min.pct1 >= ",opt$min_pct1,"& max.pct.2 <",opt$max_pct2," & ALPHA < ",opt$ALPHA,sep=''))
                        
                        if (opt$Not_only_pos) {
                            DGE_Markers$Group = 'Up'
                            DGE_Markers[DGE_Markers$avg_log2FC<0,'Group'] = 'Down'
                            
                            DGE_Markers_UP = DGE_Markers[DGE_Markers$Group == 'Up',]
                            DGE_Markers_UP = DGE_Markers_UP[DGE_Markers_UP$pct.1>=opt$min_pct1 & DGE_Markers_UP$pct.2<opt$max_pct2 & DGE_Markers_UP$p_val_adj<opt$ALPHA,]
                            
                            DGE_Markers_Down = DGE_Markers[DGE_Markers$Group == 'Down',]
                            DGE_Markers_Down = DGE_Markers_Down[DGE_Markers_Down$pct.1<opt$max_pct2 & DGE_Markers_Down$pct.2>=opt$min_pct1 & DGE_Markers_Down$p_val_adj<opt$ALPHA,]
                            
                            DGE_Markers = rbind(DGE_Markers_UP,DGE_Markers_Down)
                            
                        }else{
                            DGE_Markers = DGE_Markers[DGE_Markers$pct.1>=opt$min_pct1 & DGE_Markers$pct.2<opt$max_pct2 & DGE_Markers$p_val_adj<opt$ALPHA,]
                        }
                        # saveRDS(cluster_obj, paste(opt$outDir,cluster,'_cluster_obj.rds',sep=''))
                        writeLog(logfile, paste("FindMarkers: Found ", dim(DGE_Markers)[1] ,"Significant Genes After Filtration" ,ident1,' VS ',ident2," in Cluster ",cluster,sep=''))
                        if (dim(DGE_Markers)[1]>0){
                            DGE_Markers = cbind(DGE_Markers, cluster_counts[ DGE_Markers$gene, rownames(cluster_obj@meta.data[cluster_obj@meta.data[DGE_GroupBy]==ident1,])])
                            DGE_Markers['space'] = ''
                            DGE_Markers = cbind(DGE_Markers, cluster_counts[ DGE_Markers$gene, rownames(cluster_obj@meta.data[cluster_obj@meta.data[DGE_GroupBy]==ident2,])])
                            write.csv(DGE_Markers, file)
                        }else{
                            writeLog(logfile, paste("No Significant Genes were Found After Filtration "))
                        }
                    }else{
                        writeLog(logfile, paste("FindMarkers Found No Significant Genes"))
                    }
                } else {
                  writeLog(logfile, paste("Loaded FindMarkers results from previous run"))
                  DGE_Markers <- read.csv(file, header=T)
                }
                if (dim(DGE_Markers)[1]>0){
                    # Top Marker Genes
                    writeLog(logfile, paste("Calculating the top ",opt$Top_n," marker genes",sep=''))
                    if ('Group' %in% colnames(DGE_Markers)) {
                      Markers_Up = DGE_Markers[DGE_Markers$Group=='Up',]
                      Top_Up = Markers_Up %>% slice_head(n=opt$Top_n)
                      Markers_Down = DGE_Markers[DGE_Markers$Group=='Down',]
                      Top_Down = Markers_Down %>% slice_head(n=opt$Top_n)
                      DGE_Top_Markers = rbind(Top_Up, Top_Down)
                    } else {
                      DGE_Top_Markers <- DGE_Markers %>% slice_head(n=opt$Top_n)
                    }
                    write.csv(DGE_Top_Markers, paste(dir_path,'DGE_',ident1,"_VS_",ident2,'_Within_Cluster_',cluster,'_Markers_top',opt$Top_n,'.csv',sep=''))
                    writeLog(logfile, paste("Plotting top ",opt$Top_n," marker genes",sep=''))
                    
                    # cluster_obj         = obj_seurat[,obj_seurat@active.ident==cluster]
                    # Idents(cluster_obj) = cluster_obj@meta.data[DGE_GroupBy]
                    
                    
                    
                    if (!is.na(opt$Subsample)) {
                      heat_map = DoHeatmap(subset(cluster_obj, downsample=subsample), features = DGE_Top_Markers$gene) + ggplot2::scale_fill_gradientn(colours=c("blue","white","red")) # + NoLegend()
                    } else {
                      heat_map = DoHeatmap(cluster_obj, features = DGE_Top_Markers$gene) + ggplot2::scale_fill_gradientn(colours=c("blue","white","red")) # + NoLegend()
                    }
                    pdf(paste(dir_path,'DGE_',ident1,"_VS_",ident2,'_Within_Cluster_',cluster,'_top',opt$Top_n,'_heatmap.pdf',sep=''),width = 40, height = 15)
                    print(heat_map)
                    dev.off()
                    
                    # Enrichment Analysis
                    writeLog(logfile, paste("Setting Genes for Enrichment Analysis: min.pct1 >= ",opt$min_pct1," & ALPHA < ",opt$ALPHA,sep=''))
                    df <- DGE_Markers
                    df = df[df$pct.1>=opt$min_pct1 & df$pct.2<opt$max_pct2 & df$p_val_adj<opt$ALPHA,]
                    writeLog(logfile, paste("Contrast ",ident1,' Within Cluster ',cluster,': ',dim(df)[1],"  Differentially Expressed Genes",sep=''))
                    if (dim(df)[1]>0) {
                        # Significant Genes Heatmap
                        if ('Group' %in% colnames(df)){
                            sig_genes[[paste(cluster,"_UP",sep='')]]   <- df[df$Group=='Up',]$gene
                            sig_genes[[paste(cluster,"_Down",sep='')]] <- df[df$Group=='Down',]$gene
                        } else {
                            sig_genes[[cluster]] <- df$gene
                        }
                        if (!is.na(opt$Subsample)) {
                          heat_map = DoHeatmap(subset(cluster_obj, downsample=subsample), features = df$gene) + ggplot2::scale_fill_gradientn(colours=c("blue","white","red"))# + NoLegend()
                        } else {
                          heat_map = DoHeatmap(cluster_obj, features = df$gene) + ggplot2::scale_fill_gradientn(colours=c("blue","white","red"))# + NoLegend()
                        }
                        pdf(paste(dir_path,'DGE_',ident1,"_VS_",ident2,'_Within_Cluster_',cluster,'_heatmap.pdf',sep=''),width = 40, height = 15)
                        print(heat_map)
                        dev.off()
                    }
                }else{
                    writeLog(logfile, paste("FindMarkers Found No Significant Genes"))
                }
            }
            # GO Enrichment
            enrichGO_compClusters <- try(clusterProfiler::compareCluster(geneClusters = sig_genes, fun = "enrichGO", pAdjustMethod = opt$pAdjustMethod,
                                                                       OrgDb = OrgDb, ont = opt$ont, universe = background.genes, keyType = "SYMBOL"),
                                       silent = TRUE)
            if (!inherits(enrichGO_compClusters,"try-error")) {
              enrichGO_compClusters@compareClusterResult$GeneRatio=paste("'",enrichGO_compClusters@compareClusterResult$GeneRatio,sep='')
              write.csv(enrichGO_compClusters@compareClusterResult,
                      file=paste(opt$outDir,opt$Sample,'/',opt$Sample,"_GO_compare_Within_Clusters.csv",sep=''),quote = FALSE,row.names = TRUE)
              enrichGO_simplify <- try(clusterProfiler::simplify(enrichGO_compClusters, cutoff = opt$simplify_cutoff),
                                     silent = TRUE)
              if (!inherits(enrichGO_simplify,"try-error")) {
                enrichGO_simplify@compareClusterResult$GeneRatio=paste("'",enrichGO_simplify@compareClusterResult$GeneRatio,sep='')
                write.csv(enrichGO_simplify@compareClusterResult,
                        file=paste(opt$outDir,opt$Sample,'/',opt$Sample,"_GO_simplify_Within_Clusters.csv",sep=''),quote = FALSE,row.names = TRUE)
                if (dim(enrichGO_simplify@compareClusterResult)[1]<250) {
                  fontsize=9
                } else {
                  fontsize=4
                }
                enrichGO_simplify@compareClusterResult$GeneRatio=stringr::str_remove_all(enrichGO_simplify@compareClusterResult$GeneRatio,"'")
                enrichGO_dotplot <- enrichplot::dotplot(enrichGO_simplify, showCategory=10000, font.size=fontsize, label_format=100)
                ggplot2::ggsave(filename = paste(opt$outDir,opt$Sample,'/',opt$Sample,"_GO_simplify_Within_Clusters_dotplot.pdf",sep=''),
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
                      file=paste(opt$outDir,opt$Sample,'/',opt$Sample,"_KEGG_compareClusters_Within_Clusters.csv",sep=''),quote = FALSE,row.names = TRUE)
              if (dim(enrichKEGG_compClusters@compareClusterResult)[1]<250) {
                fontsize=9
              } else {
                fontsize=4
              }
              enrichKEGG_compClusters@compareClusterResult$GeneRatio=stringr::str_remove_all(enrichKEGG_compClusters@compareClusterResult$GeneRatio,"'")
              enrichKEGG_dotplot <- enrichplot::dotplot(enrichKEGG_compClusters, showCategory=10000, font.size=fontsize, label_format=100)
              ggplot2::ggsave(filename = paste(opt$outDir,opt$Sample,'/',opt$Sample,"_KEGG_compare_Within_Clusters.pdf",sep=''),
                            enrichKEGG_dotplot, dpi = 600,device = "pdf",width = 30,height = 50, limitsize = FALSE)
              rm(sig_genes,cluster,file,df,sig_genes_entrez,enrichGO_compClusters,enrichGO_simplify,enrichGO_dotplot,enrichKEGG_compClusters,
                 enrichKEGG_dotplot)
            }
            
        }else{
            writeLog(logfile, paste("ERROR: Cannot find the following Group in the object: ",
                        paste(idents[which(!(DGE_GroupBy %in% colnames(obj_seurat@meta.data)))],collapse=','),sep=''))
            stop()
            
        }
    }
}



### DGE analysis for clusters vs each other ###
if (opt$DGE_Pairwise) {
    writeLog(logfile, paste("Running Pairwise-DGE Analysis: Cluster vs. other"))
    ident1 = unlist(stringi::stri_split(str = opt$ident1,fixed = ',')) 
    ident2 = unlist(stringi::stri_split(str = opt$ident2,fixed = ',')) 
    idents = c(ident1, ident2)
    sig_genes = list()
    
    if (isFALSE(all(idents %in% Idents(obj_seurat)))) {
        writeLog(logfile, paste("ERROR: Cannot find the following identities in the object: ",
                        paste(idents[which(!(idents %in% Idents(obj_seurat)))],collapse=','),sep=''))
        stop()
    }
    dir_path = paste(opt$outDir,opt$Sample,'/FindMarkers_vs_each_other/',paste(ident1,collapse = '.'),'_vs_',paste(ident2,collapse = '.'),'/',sep='')
    if (!dir.exists(dir_path))
        dir.create(dir_path)
    
    # FindMarkers
    file = paste(dir_path,'DGE_',paste(ident1,collapse = '.'),'_vs_',paste(ident2,collapse = '.'),'_Markers.csv',sep='')
    if (any(opt$Force,!file.exists(file))) {
      writeLog(logfile, paste("Running FindMarkers: ",paste(ident1,collapse = '.'),' vs. ',paste(ident2,collapse = '.'),sep=''))
      obj_seurat = PrepSCTFindMarkers(object = obj_seurat )
      DGE_Markers <- FindMarkers(obj_seurat, ident.1 = ident1, ident.2 = ident2, only.pos = !opt$Not_only_pos, min.pct = opt$min_pct, 
                                    logfc.threshold = opt$logfc_threshold, test.use = opt$test_use, verbose = FALSE)
      DGE_Markers$gene = rownames(DGE_Markers)
      writeLog(logfile, paste("Finished FindMarkers, Saving results..."))
      if (opt$Not_only_pos) {
        DGE_Markers$Group = 'Up'
        DGE_Markers[DGE_Markers$avg_log2FC<0,'Group'] = 'Down'
        
        DGE_Markers_UP = DGE_Markers[DGE_Markers$Group == 'Up',]
        DGE_Markers_UP = DGE_Markers_UP[DGE_Markers_UP$pct.1>=opt$min_pct1 & DGE_Markers_UP$pct.2<opt$max_pct2 & DGE_Markers_UP$p_val_adj<opt$ALPHA,]
        
        DGE_Markers_Down = DGE_Markers[DGE_Markers$Group == 'Down',]
        DGE_Markers_Down = DGE_Markers_Down[DGE_Markers_Down$pct.1<opt$max_pct2 & DGE_Markers_Down$pct.2>=opt$min_pct1 & DGE_Markers_Down$p_val_adj<opt$ALPHA,]
        
        DGE_Markers = rbind(DGE_Markers_UP,DGE_Markers_Down)
        
      }else{
        DGE_Markers = DGE_Markers[DGE_Markers$pct.1>=opt$min_pct1 & DGE_Markers$pct.2<opt$max_pct2 & DGE_Markers$p_val_adj<opt$ALPHA,]
      }
      write.csv(DGE_Markers, file)
    } else {
      writeLog(logfile, paste("Loaded FindMarkers results from previous run"))
      DGE_Markers <- read.csv(file, header=T)
    }
    
    # Top Marker Genes
    writeLog(logfile, paste("Calculating the top ",opt$Top_n," marker genes",sep=''))
    if ('Group' %in% colnames(DGE_Markers)) {
      Markers_Up = DGE_Markers[DGE_Markers$Group=='Up',]
      Top_Up = Markers_Up %>% slice_head(n=opt$Top_n)
      Markers_Down = DGE_Markers[DGE_Markers$Group=='Down',]
      Top_Down = Markers_Down %>% slice_head(n=opt$Top_n)
      DGE_Top_Markers = rbind(Top_Up, Top_Down)
    } else {
      DGE_Top_Markers <- DGE_Markers %>% slice_head(n=opt$Top_n)
    }
    write.csv(DGE_Top_Markers, paste(dir_path,'DGE_',paste(ident1,collapse = '.'),'_vs_',paste(ident2,collapse = '.'),'_Markers_top',opt$Top_n,'.csv',sep=''))
    writeLog(logfile, paste("Plotting top ",opt$Top_n," marker genes",sep=''))
    if (!is.na(opt$Subsample)) {
      heat_map = DoHeatmap(subset(obj_seurat, downsample=subsample), features = DGE_Top_Markers$gene) + NoLegend()
    } else {
      heat_map = DoHeatmap(obj_seurat, features = DGE_Top_Markers$gene) + NoLegend()
    }
    pdf(paste(dir_path,'DGE_',paste(ident1,collapse = '.'),'_vs_',paste(ident2,collapse = '.'),'_top',opt$Top_n,'_heatmap.pdf',sep=''),width = 40, height = 15)
    print(heat_map)
    dev.off()
    
    # Enrichment Analysis
    writeLog(logfile, paste("Setting Genes for Enrichment Analysis: ALPHA < ",opt$ALPHA,sep=''))
    df <- DGE_Markers
    if ('Group' %in% colnames(df))
        df = df[df$Group=='Up',]
    df = df[df$p_val_adj<opt$ALPHA,]
    writeLog(logfile, paste("Contrast ",paste(ident1,collapse = '.'),' vs. ',paste(ident2,collapse = '.'),': ',dim(df)[1]," Up-regulated Differentially Expressed Genes",sep=''))
    if (dim(df)[1]>0) {
        # Significant Genes Heatmap
        sig_genes <- df$gene
        if (!is.na(opt$Subsample)) {
          heat_map = DoHeatmap(subset(obj_seurat, downsample=subsample), features = sig_genes) + NoLegend()
        } else {
          heat_map = DoHeatmap(obj_seurat, features = sig_genes) + NoLegend()
        }
        pdf(paste(dir_path,'DGE_',paste(ident1,collapse = '.'),'_vs_',paste(ident2,collapse = '.'),'_Up_Regulated_heatmap.pdf',sep=''),width = 40, height = 15)
        print(heat_map)
        dev.off()
        
        # GO Enrichment
        enrichGO_compClusters <- try(clusterProfiler::compareCluster(geneClusters = sig_genes, fun = "enrichGO", pAdjustMethod = opt$pAdjustMethod,
                                                                   OrgDb = OrgDb, ont = opt$ont, universe = background.genes, keyType = "SYMBOL"),
                                   silent = TRUE)
        if (!inherits(enrichGO_compClusters,"try-error")) {
            enrichGO_compClusters@compareClusterResult$GeneRatio=paste("'",enrichGO_compClusters@compareClusterResult$GeneRatio,sep='')
            write.csv(enrichGO_compClusters@compareClusterResult,
                      file=paste(dir_path,'DGE_',paste(ident1,collapse = '.'),'_vs_',paste(ident2,collapse = '.'),'_GO_compareClusters.csv',sep=''),quote = FALSE,row.names = TRUE)
            enrichGO_simplify <- try(clusterProfiler::simplify(enrichGO_compClusters, cutoff = opt$simplify_cutoff),
                                     silent = TRUE)
            if (!inherits(enrichGO_simplify,"try-error")) {
              enrichGO_simplify@compareClusterResult$GeneRatio=paste("'",enrichGO_simplify@compareClusterResult$GeneRatio,sep='')
              write.csv(enrichGO_simplify@compareClusterResult,
                        file=paste(dir_path,'DGE_',paste(ident1,collapse = '.'),'_vs_',paste(ident2,collapse = '.'),'_GO_simplify.csv',sep=''),quote = FALSE,row.names = TRUE)
              
              if (dim(enrichGO_simplify@compareClusterResult)[1]<250) {
                fontsize=9
              } else {
                fontsize=4
              }
              enrichGO_simplify@compareClusterResult$GeneRatio=stringr::str_remove_all(enrichGO_simplify@compareClusterResult$GeneRatio,"'")
              enrichGO_dotplot <- enrichplot::dotplot(enrichGO_simplify, showCategory=10000, font.size=fontsize, label_format=100)
              ggplot2::ggsave(filename = paste(dir_path,'DGE_',paste(ident1,collapse = '.'),'_vs_',paste(ident2,collapse = '.'),'_GO_simplify_dotplot.pdf',sep=''),
                              enrichGO_dotplot, dpi = 600,device = "pdf",width = 30,height = 50, limitsize = FALSE)
            }
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

