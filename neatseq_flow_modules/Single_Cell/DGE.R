.libPaths()
library(optparse)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
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
  make_option(c("--DGE_Clusters_vs_all"), action="store_true", default = FALSE,
              help="Run Clusters vs all other differential gene expression analysis of clusters vs each other (Default is False)", metavar = "character"),
  make_option(c("--DGE_Within_Clusters"), action="store_true", default = FALSE,
              help="Run Within Clusters differential gene expression analysis [use --ident1 and --ident2 to set the contrast and --DGE_GroupBy to set the source ] (Default is False)", metavar = "character"),
  make_option(c("--CoExpTest"), action="store_true", default = FALSE,
              help="Run Within Clusters Co-Expression interaction Analysis [ONLY works when --DGE_Within_Clusters is set. Use --DGE_GroupBy to set the source ] (Default is False)", metavar = "character"),
  make_option(c("--CoExpMagic"), action="store_true", default = FALSE,
              help="Use Magic assay for the CoExpTest (Default is False: use RNA assay)", metavar = "character"),
  make_option(c("--CoExpRcutoff"), type="numeric", default = 0.8,
              help="Only consider correlations between genes that are higher than this value  (Default is 0.8)", metavar = "character"),
  make_option(c("--CoExpVARcutoff"), type="numeric", default = 0.4,
              help="Filter genes for correlations with Variance smaller than this value  (Default is 0.4)", metavar = "character"),
  make_option(c("--CoExpPlotPcutoff"), type="numeric", default = 1e-8,
              help="Only Plot genes with interaction P-value smaller than this value  (Default is 1e-8)", metavar = "character"),
  make_option(c("--CoExpTestUniform"), action="store_true", default = FALSE,
              help="Use only genes with uniformly distributed expression values", metavar = "character"),
  make_option(c("--Chisq"), action="store_true", default = FALSE,
              help="Run Within Clusters differential Cell Proportion Chisq TEST [ONLY works when --DGE_Within_Clusters is set. Use --ident1 and --ident2 to set the contrast and --DGE_GroupBy to set the source ] (Default is False)", metavar = "character"),
  make_option(c("--Chisq_only_top_genes"), action="store_true", default = FALSE,
              help="Use only --Top_n genes for the Enrichment Analysis  (Default is False)", metavar = "character"),
  make_option(c("--Chisq_minprop_cutoff"), type="numeric", default = 0,
              help="Only consider genes with proportion of Positive cells greater then the cutoff in ident1 OR ident2 (Default is 0)", metavar = "character"),
  make_option(c("--Chisq_diffprop_cutoff"), type="numeric", default = 0,
              help="Only consider genes with a difference in Positive cells proportion between ident1 and ident2 greater then the cutoff (Default is 0)", metavar = "character"),
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
  make_option(c("--StackedViolinPlot"), action="store_true", default = FALSE,
              help="Plot input features stacked violin-plots", metavar = "character"),
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
# obj_seurat <- SetIdent(obj_seurat, value= obj_seurat@meta.data$seurat_clusters)
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



if (!is.na(opt$Treatment_file)) {
    treatmentTable <- read.table(opt$Treatment_file, sep = '\t', header = TRUE)
    if (dim(treatmentTable)[1]<2) {
      writeLog(logfile, paste("ERROR: Invalid number of columns."))
      stop()
    }
    if (!(opt$DGE_GroupBy %in% colnames(treatmentTable))) {
      writeLog(logfile, paste("ERROR: Your DGE_GroupBy parameter:",opt$DGE_GroupBy," is not in the Treatment_file columns names"))
      stop()
    }
    ind= match(obj_seurat@meta.data$orig.ident, treatmentTable[[1]])
    obj_seurat@meta.data[opt$DGE_GroupBy]= treatmentTable[ind,opt$DGE_GroupBy]
    
}


Prefix  = ""
PosCell = "_PosCell"
NegCell = "_NegCell"
Proportion = "_Proportion"
Expression_Cutoff = 0
Ncell_Cutoff      = 50

CellProp = c()
DGE_GroupBy = opt$DGE_GroupBy
if (DGE_GroupBy %in% colnames(obj_seurat@meta.data)){
    for (cluster in clusters) {
        ClusterCellProp = data.frame()
        # cluster_obj     = obj_seurat[,obj_seurat@meta.data$seurat_clusters==cluster]
        # cluster_obj     = obj_seurat[,rownames(obj_seurat@meta.data[obj_seurat@meta.data$seurat_clusters==cluster,])]
        cluster_obj     = obj_seurat[, names(Idents(obj_seurat)[Idents(obj_seurat)==cluster])]
        cluster_counts  = as.matrix(GetAssayData(object = cluster_obj, slot = "count", assay = "RNA"))
        cells           = apply(X= cluster_counts>Expression_Cutoff,MARGIN=1,FUN=sum)
        Genes           = names(cells[cells>Ncell_Cutoff])
        # fetch_data      = FetchData(cluster_obj, vars = c(DGE_GroupBy,Genes))
        fetch_data      = FetchData(cluster_obj, vars = c(DGE_GroupBy))
        print(cluster)
        if (length(Genes)>0){
            for (ident in unique(obj_seurat@meta.data[,DGE_GroupBy])){
                Ident_Cells      = rownames(fetch_data)[fetch_data[DGE_GroupBy]==ident]
                Toatal_Cells     = length(Ident_Cells)
                if (Toatal_Cells>0){
                    ExpressedInCells = apply(X=as.data.frame(cluster_counts[Genes,Ident_Cells])>Expression_Cutoff,MARGIN=1,FUN=sum)
                    ExpressedInCells = as.data.frame(ExpressedInCells)
                    ClusterCellProp[Genes,paste(ident,PosCell,sep='')]    = ExpressedInCells
                    ClusterCellProp[Genes,paste(ident,NegCell,sep='')]    = Toatal_Cells-ExpressedInCells
                    ClusterCellProp[Genes,paste(ident,Proportion,sep='')] = ExpressedInCells/Toatal_Cells
                }else{
                    ClusterCellProp[Genes,paste(ident,PosCell,sep='')]    = Toatal_Cells
                    ClusterCellProp[Genes,paste(ident,NegCell,sep='')]    = Toatal_Cells
                    ClusterCellProp[Genes,paste(ident,Proportion,sep='')] = 0
                }
            }
            file = paste(opt$outDir,opt$Sample,'_Cluster',cluster,"_GeneProp.csv",sep='')
            write.csv(ClusterCellProp, file)
            CellProp[cluster] = list(ClusterCellProp)
        }
    }
}

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
          fontsize=12
        } else {
          fontsize=12
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
                                 silent = FALSE)
    if (!inherits(enrichKEGG_compClusters,"try-error")) {
        if (!is.null(enrichKEGG_compClusters)){
            enrichKEGG_compClusters@compareClusterResult$GeneRatio=paste("'",enrichKEGG_compClusters@compareClusterResult$GeneRatio,sep='')
            write.csv(enrichKEGG_compClusters@compareClusterResult,
                    file=paste(opt$outDir,opt$Sample,'/',opt$Sample,"_KEGG_compareClusters_vs_all_others.csv",sep=''),quote = FALSE,row.names = TRUE)
            if (dim(enrichKEGG_compClusters@compareClusterResult)[1]<250) {
              fontsize=12
            } else {
              fontsize=10
            }
            enrichKEGG_compClusters@compareClusterResult$GeneRatio=stringr::str_remove_all(enrichKEGG_compClusters@compareClusterResult$GeneRatio,"'")
            enrichKEGG_dotplot <- enrichplot::dotplot(enrichKEGG_compClusters, showCategory=10000, font.size=fontsize, label_format=100)
            ggplot2::ggsave(filename = paste(opt$outDir,opt$Sample,'/',opt$Sample,"_KEGG_compareClusters_vs_all_others.pdf",sep=''),
                          enrichKEGG_dotplot, dpi = 600,device = "pdf",width = 30,height = 50, limitsize = FALSE)
            rm(sig_genes,cluster,file,df,sig_genes_entrez,enrichGO_compClusters,enrichGO_simplify,enrichGO_dotplot,enrichKEGG_compClusters,
               enrichKEGG_dotplot)
        }
    }
}


if (opt$DGE_Within_Clusters) {
    writeLog(logfile, paste("Running Within Clusters-DGE Analysis"))
    opt$Not_only_pos = TRUE
    sig_genes = list()
    # if (!is.na(opt$Treatment_file)) {
        # treatmentTable <- read.table(opt$Treatment_file, sep = '\t', header = TRUE)
        # if (isFALSE(all(colnames(treatmentTable) == c(colnames(treatmentTable)[1],colnames(treatmentTable)[2])))) {
          # writeLog(logfile, paste("ERROR: Invalid column names. Please use SampleID and Treatment"))
          # stop()
        # }
        # ind= match(obj_seurat@meta.data$orig.ident, treatmentTable$SampleID)
        # obj_seurat@meta.data$treatment= treatmentTable$Treatment[ind]
        # #Subset only dent1 and ident2 samples:
        # ident1 <- opt$ident1
        # ident2 <- opt$ident2
        # obj_seurat <- obj_seurat[ ,grepl(ident1, obj_seurat$treatment) | grepl(ident2, obj_seurat$treatment)]
    # }
    if (!is.na(opt$ident1)){
        DGE_GroupBy = opt$DGE_GroupBy
        df = data.frame()
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
                cluster_obj = obj_seurat[, names(Idents(obj_seurat)[Idents(obj_seurat)==cluster])]
                Idents(cluster_obj) = cluster_obj@meta.data[,DGE_GroupBy]
                
                
                if (opt$CoExpTest){
                    dir_path = paste(opt$outDir,opt$Sample,'/FindMarkers_Within_Clusters/','Within_Cluster_',cluster,'/',sep='')
                    if (!dir.exists(dir_path))
                        dir.create(dir_path)
                    
                }else{
                    dir_path = paste(opt$outDir,opt$Sample,'/FindMarkers_Within_Clusters/',ident1,'_Within_Cluster_',cluster,'/',sep='')
                    if (!dir.exists(dir_path))
                        dir.create(dir_path)

                    # FindMarkers
                    Prefix = paste(ident1,"_VS_",ident2,sep='')
                    Prefix = paste('DGE_Within_Clusters_',Prefix,sep='')
                    file   = paste(dir_path,Prefix,cluster,'_Markers.csv',sep='')
                    writeLog(logfile, paste("Running FindMarkers: ",ident1,'_Within_Cluster_',cluster,sep=''))
                    DGE_Markers = data.frame()
                }
                if ("SCT" %in% names(cluster_obj@assays)){
                    assay2use   = "SCT"
                }else{
                    assay2use   = cluster_obj@active.assay
                }
                slot2use    = "scale.data"
                
                colours=c("blue","white","red")
                if (opt$Chisq){
                    if (cluster %in% names(CellProp)){
                        Prefix = paste(ident2,"_VS_",ident1,sep='')
                        Prefix = paste('Chisq_CellProp_Within_Cluster_',Prefix,sep='')
                        file = paste(dir_path,Prefix,cluster,'_Markers.csv',sep='')
                        ClusterCellProp = CellProp[[cluster]]
                        if (dim(ClusterCellProp)[1]>0){
                            ident1_PosCell          = which(colnames(ClusterCellProp) %in% paste(ident1,PosCell,sep=''))
                            ident2_PosCell          = which(colnames(ClusterCellProp) %in% paste(ident2,PosCell,sep=''))
                            ident1_NegCell          = which(colnames(ClusterCellProp) %in% paste(ident1,NegCell,sep=''))
                            ident2_NegCell          = which(colnames(ClusterCellProp) %in% paste(ident2,NegCell,sep=''))
                            
                            ident1_Proportion       = which(colnames(ClusterCellProp) %in% paste(ident1,Proportion,sep=''))
                            ident2_Proportion       = which(colnames(ClusterCellProp) %in% paste(ident2,Proportion,sep=''))
                            
                            gene_list = apply(X= ClusterCellProp,MARGIN=1,FUN= function(x) as.numeric( x[ident2_Proportion]) - as.numeric( x[ident1_Proportion]) )
                            # names(gene_list)=rownames(ClusterCellProp)
                            gene_list = gene_list[sort.list(gene_list,decreasing = T)]
                            saveRDS(gene_list,file=paste(dir_path,Prefix,cluster,'_GeneList',sep=''))
                            filter_rows             = apply(X= ClusterCellProp,MARGIN=1,FUN= function(x) sum(c(as.numeric(x[ident1_Proportion]),as.numeric(x[ident2_Proportion])) > opt$Chisq_minprop_cutoff)>0)
                            ClusterCellProp         = ClusterCellProp[filter_rows,]
                            
                            filter_rows             = apply(X= ClusterCellProp,MARGIN=1,FUN= function(x) abs(as.numeric(x[ident1_Proportion]) - as.numeric(x[ident2_Proportion]) )  > opt$Chisq_diffprop_cutoff )
                            ClusterCellProp         = ClusterCellProp[filter_rows,]
                            if (dim(ClusterCellProp)[1]>0) {
                                ClusterCellProp$p.value = apply(X= ClusterCellProp,MARGIN=1,FUN= function(x) chisq.test(x = c(x[ident2_PosCell]+1,x[ident2_NegCell]+1),
                                                                                              p = c(x[ident1_PosCell]+1,x[ident1_NegCell]+1),
                                                                                              rescale.p        = TRUE,
                                                                                              simulate.p.value = FALSE)$p.value)
                                
                                ClusterCellProp["Group"]           = 'Up'
                                ClusterCellProp["Proportion_Diff"] = apply(X= ClusterCellProp,MARGIN=1,FUN= function(x) as.numeric(x[ident2_Proportion]) - as.numeric(x[ident1_Proportion]))
                                # ClusterCellProp[apply(X= ClusterCellProp,MARGIN=1,FUN= function(x) as.numeric(x[ident2_Proportion]) < as.numeric(x[ident1_Proportion])),"Group"]  = 'Down'
                                ClusterCellProp[ClusterCellProp["Proportion_Diff"]<0,"Group"] = 'Down'
                                ClusterCellProp$p.adjust   = p.adjust(p =ClusterCellProp$p.value ,method = opt$pAdjustMethod)
                                ClusterCellProp            = ClusterCellProp[order(ClusterCellProp$p.adjust, decreasing = TRUE),]
                                ClusterCellProp            = ClusterCellProp[sort.list(ClusterCellProp$p.adjust),]
                                
                                ClusterCellProp[ClusterCellProp$p.adjust>opt$ALPHA,"Group"] = "Not_Significant"
                                write.csv(ClusterCellProp, file)
                                
                                volcanoPlot = ggplot2::ggplot(data = ClusterCellProp,mapping = aes(x=Proportion_Diff,y=-log10(p.adjust) ))+
                                                       geom_point(aes(colour = Group ),size = 1)+
                                                       scale_color_manual(values=c("black","blue","red"))+
                                                       xlim(-3, 3) +
                                                       theme_classic()
                                pdf(paste(dir_path,Prefix,cluster,'_volcanoPlot_Chisq.pdf',sep=''),width = 40, height = 15)
                                print(volcanoPlot)
                                dev.off()
                                
                                df        = ClusterCellProp[ClusterCellProp$p.adjust< opt$ALPHA,]
                                df$gene   = rownames(df)
                                df = df[order(abs(df$Proportion_Diff),decreasing = T), ]
                                UP_genes = df[df$Group=='Up',]
                                if (length(UP_genes$gene)>0){
                                    if (length(UP_genes$gene)>opt$Top_n){
                                        pheatmap::pheatmap(mat = ClusterCellProp[UP_genes$gene[1:opt$Top_n],stringi::stri_endswith(str = colnames(ClusterCellProp),fixed = Proportion)],
                                                       cluster_rows = F,cluster_cols = T,
                                                       filename = paste(dir_path,Prefix,cluster,"_Top_",opt$Top_n,'_UP_ChisqHeatMap.pdf',sep=''))
                                        if (opt$Chisq_only_top_genes){
                                            df = df[1:opt$Top_n,]
                                        }
                                    }else{
                                        pheatmap::pheatmap(mat = ClusterCellProp[UP_genes$gene,stringi::stri_endswith(str = colnames(ClusterCellProp),fixed = Proportion)],
                                                       cluster_rows = F,cluster_cols = T,
                                                       filename = paste(dir_path,Prefix,cluster,'_UP_ChisqHeatMap.pdf',sep=''))
                                    }
                                }
                                
                                Down_genes = df[df$Group=='Down',]
                                if (length(Down_genes$gene)>0){
                                    if (length(Down_genes$gene)>opt$Top_n){
                                        pheatmap::pheatmap(mat = ClusterCellProp[Down_genes$gene[1:opt$Top_n],stringi::stri_endswith(str = colnames(ClusterCellProp),fixed = Proportion)],
                                                       cluster_rows = F,cluster_cols = T,
                                                       filename = paste(dir_path,Prefix,cluster,"_Top_",opt$Top_n,'_Down_ChisqHeatMap.pdf',sep=''))
                                        if (opt$Chisq_only_top_genes){
                                            df = df[1:opt$Top_n,]
                                        }
                                    }else{
                                        pheatmap::pheatmap(mat = ClusterCellProp[Down_genes$gene,stringi::stri_endswith(str = colnames(ClusterCellProp),fixed = Proportion)],
                                                       cluster_rows = F,cluster_cols = T,
                                                       filename = paste(dir_path,Prefix,cluster,'_Down_ChisqHeatMap.pdf',sep=''))
                                    }
                                }
                                # df = df[order(df$Proportion_Diff,decreasing = T), ]
                                GOgse = try(clusterProfiler::gseGO(geneList=gene_list,
                                                                   pAdjustMethod = opt$pAdjustMethod,
                                                                   OrgDb = OrgDb,
                                                                   ont = opt$ont,
                                                                   keyType = "SYMBOL"),
                                          silent = F)
                                if (!inherits(GOgse,"try-error")) {
                                    if (dim(GOgse)[1]>1){
                                        write.csv(GOgse,
                                                  file=paste(dir_path,Prefix,cluster,'_GOgse.csv',sep=''))
                                        GOgse_plt = clusterProfiler::ridgeplot(GOgse,core_enrichment = TRUE)
                                        pdf(paste(dir_path,Prefix,cluster,'_GOgse.pdf',sep=''),width = 40, height = 15)
                                        print(GOgse_plt)
                                        dev.off()
                                    }
                                }
                                
                                conv = clusterProfiler::bitr(geneID=names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
                                names(gene_list) = conv$ENTREZID
                                saveRDS(gene_list,file=paste(dir_path,Prefix,cluster,'_ENTREZID_GeneList',sep=''))
                                
                                KEGGgse = try(clusterProfiler::gseKEGG(geneList = gene_list,
                                                                       pAdjustMethod = opt$pAdjustMethod,
                                                                       organism = opt$organism),
                                          silent = F)
                                if (!inherits(KEGGgse,"try-error")) {
                                    if (dim(KEGGgse)[1]>1){
                                        write.csv(KEGGgse,
                                                  file=paste(dir_path,Prefix,cluster,'_KEGGgse.csv',sep=''))
                                        KEGGgse_plt = clusterProfiler::ridgeplot(KEGGgse,core_enrichment = TRUE)
                                        pdf(paste(dir_path,Prefix,cluster,'_KEGGgse.pdf',sep=''),width = 40, height = 15)
                                        print(KEGGgse_plt)
                                        dev.off()
                                    }
                                }
                                
                            }else{
                                df       = data.frame()
                                writeLog(logfile, paste("No Significant Genes were Found After Filtration "))
                            }
                            slot2use  = "counts"
                            assay2use = "RNA"
                            colours=c("white","black")
                        }
                    }
                }else if (opt$CoExpTest){
                    CoExpTest     = data.frame()
                    if (opt$CoExpMagic & ("MAGIC_RNA" %in% names(cluster_obj@assays))){
                        slot2use = "data"
                        cluster_exp   = as.matrix(GetAssayData(object = cluster_obj, slot = slot2use, assay = "MAGIC_RNA"))   
                    }else{
                        slot2use    = "data"
                        # slot2use    = "scale.data"
                        cluster_exp   = as.matrix(GetAssayData(object = cluster_obj, slot = slot2use, assay = assay2use))
                    }
                    
                    fetch_data    = FetchData(cluster_obj, vars = c(DGE_GroupBy))
                    T_cluster_exp =  as.data.frame( t(cluster_exp) )
                    T_cluster_exp[rownames(fetch_data),DGE_GroupBy] = fetch_data[DGE_GroupBy]
                    Genes2Test    = c()
                    # save.image(paste(dir_path,"Cluster_",cluster,".rds",sep=''))
                    for (ident in unique(obj_seurat@meta.data[,DGE_GroupBy])){
                        Ident_Cells      = rownames(fetch_data)[fetch_data[DGE_GroupBy]==ident]
                        Toatal_Cells     = length(Ident_Cells)
                        if (Toatal_Cells>10){
                            Ident_cluster_exp = cluster_exp[ ,Ident_Cells]
                            genes2use_var  = apply(X= Ident_cluster_exp,MARGIN=1,FUN=var)>opt$CoExpVARcutoff
                            genes2use_var  = names(genes2use_var)[genes2use_var]
                            if (opt$CoExpTestUniform){
                                genes2use_unif = apply(X= Ident_cluster_exp,MARGIN=1,FUN=function(Values) ks.test(Values, "punif", min(Values), max(Values))$p.value)>0.05
                                genes2use_unif = names(genes2use_unif)[genes2use_unif]
                                genes2use      = intersect(genes2use_unif,genes2use_var)
                            }else{
                                genes2use = genes2use_var
                            }
                            
                            if (length(genes2use)>1){
                                Ident_cor  = cor(t(Ident_cluster_exp[genes2use,]))
                                Ident_unif = apply(X= Ident_cluster_exp[genes2use,],MARGIN=1,FUN=function(Values) ks.test(Values, "punif", min(Values), max(Values))$statistic[["D"]])
                                for (geneA in colnames(Ident_cor)){
                                    genes2test = row.names(Ident_cor)[abs(Ident_cor[,geneA])>opt$CoExpRcutoff]
                                    genes2test = genes2test[!(genes2test %in% geneA)]
                                    
                                    for (geneB in genes2test){
                                        if (cor.test(Ident_cluster_exp[geneA,],Ident_cluster_exp[geneB,])$p.value<0.05){
                                            ForLm = data.frame(geneA =T_cluster_exp[,geneA],geneB =T_cluster_exp[,geneB],Treatment = T_cluster_exp[,DGE_GroupBy] )
                                            ForLm = ForLm[apply(ForLm[,c('geneA','geneB')],MARGIN = 1,FUN = function(x) x[1]*x[2])>0,]
                                            if (length(unique(ForLm$Treatment))>1){
                                                Lm_Results = anova(lm(formula = geneA ~  geneB*Treatment,data = ForLm))
                                                pval = Lm_Results["geneB:Treatment","Pr(>F)"]
                                                CoExpTest_temp = data.frame(geneA = geneA ,geneB = geneB,geneA_D = Ident_unif[[geneA]],geneB_D = Ident_unif[[geneB]],pval = pval,ident = ident,R=Ident_cor[geneA,geneB])
                                                CoExpTest = rbind(CoExpTest,CoExpTest_temp)
                                                if (!is.na(pval)){
                                                    if (pval < opt$CoExpPlotPcutoff){
                                                        gg = ggplot(data = ForLm,mapping =  aes(x = geneA, y =geneB , color = Treatment)) +
                                                                    geom_point() +
                                                                    xlab(geneA) +
                                                                    ylab(geneB) +
                                                                    geom_smooth(method = "lm", se = FALSE)+ 
                                                                    ggplot2::ggtitle(label = paste("Cluster:",cluster,"pval:",pval,sep=""))
                                                        ggplot2::ggsave(filename = paste(dir_path,"Cluster_",cluster,"_",ident,"_",geneA,"_",geneB,".pdf",sep=''),
                                                                        gg, dpi = 600,device = "pdf",width = 10,height = 10, limitsize = FALSE)
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    
                                    Genes2Test = c(Genes2Test,genes2test)
                                }
                            }
                        }
                    }
                    CoExpTest$p.adjust   = p.adjust(p =CoExpTest$pval ,method = opt$pAdjustMethod)
                    write.csv(CoExpTest, paste(dir_path,"CoExpTest_Cluster_",cluster,".csv",sep=''))
                    

                }else{
                    if (any(opt$Force,!file.exists(file))) {
                        writeLog(logfile, paste("Running FindMarkers: ",ident1,'_Within_Cluster_',cluster,sep=''))
                        #cluster_obj = obj_seurat[,obj_seurat@active.ident==cluster]
                        
                        
                        if (cluster_obj@active.assay == "SCT"){
                            cluster_obj = PrepSCTFindMarkers(object = cluster_obj)
                        }
                        
                        if ("SCT" %in% names(cluster_obj@assays)){
                            cluster_counts = GetAssayData(object = cluster_obj, slot = slot2use, assay = "SCT")
                        }else{
                            cluster_counts = GetAssayData(object = cluster_obj, slot = slot2use, assay = "RNA")
                        }
                        DGE_Markers <- FindMarkers(cluster_obj, ident.1 = ident1,ident.2 = ident2, only.pos = !opt$Not_only_pos, min.pct = opt$min_pct, 
                                                      logfc.threshold = opt$logfc_threshold, test.use = opt$test_use, verbose = FALSE)
                        writeLog(logfile, paste("Finished FindMarkers, Saving results..."))
                        writeLog(logfile, paste("FindMarkers: Found ", dim(DGE_Markers)[1] ,"Significant Genes Before Filtration" ,ident1,' VS ',ident2," in Cluster ",cluster,sep=''))
                        
                        if (dim(DGE_Markers)[1]>0){
                            DGE_Markers$gene = rownames(DGE_Markers)
                            
                            writeLog(logfile, paste("Filtering Genes After DGE Analysis: min.pct1 >= ",opt$min_pct1,"& max.pct.2 <",opt$max_pct2," & ALPHA < ",opt$ALPHA,sep=''))
                            write.csv(DGE_Markers, file)
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
                            # save.image(file =paste(opt$outDir,"WorkSpace.RData",sep='') )
                            writeLog(logfile, paste("FindMarkers: Found ", dim(DGE_Markers)[1] ,"Significant Genes After Filtration" ,ident1,' VS ',ident2," in Cluster ",cluster,sep=''))
                            # if (dim(DGE_Markers)[1]>0){
                                # DGE_Markers = cbind(DGE_Markers, cluster_counts[ DGE_Markers$gene, rownames(cluster_obj@meta.data[cluster_obj@meta.data[DGE_GroupBy]==ident1,])])
                                # DGE_Markers['space'] = ''
                                # DGE_Markers = cbind(DGE_Markers, cluster_counts[ DGE_Markers$gene, rownames(cluster_obj@meta.data[cluster_obj@meta.data[DGE_GroupBy]==ident2,])])
                                # write.csv(DGE_Markers, file)
                            # }else{
                                # writeLog(logfile, paste("No Significant Genes were Found After Filtration "))
                            # }
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
                        
                                            
        #DGE_all <- FindMarkers(cluster_obj, ident.1 = ident1,ident.2 = ident2, only.pos = !opt$Not_only_pos, min.pct = 0, 
    #                    logfc.threshold = -Inf, test.use = opt$test_use, verbose = FALSE)
 #Volcano plot:
   #     volcano.plot <- ggplot(aes(DGE_all, aes(x = avg_log2FC, y = -log10(p_val_adj)))) +
 #geom_point(aes(color = ifelse(p_val_adj < 0.05, "Significant", "Not Significant")), alpha = 0.6) +
  #scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  #labs(x = "Log2(Fold Change)", y = "-log10(P-value)") +  geom_text_repel() + theme_classic()
  #
             #      pdf(paste(dir_path,'Volcano_',ident1,"_VS_",ident2,'_Within_Cluster_',cluster,'.pdf',sep=''),width = 40, height = 15)
        #            print(volcano.plot)
        #            dev.off() 
                    
                        if (!is.na(opt$Subsample)) {
                          violin_plot <- VlnPlot(subset(cluster_obj, downsample=subsample), features= DGE_Top_Markers$gene, log = TRUE) + ggplot2::theme(axis.text.y = element_text(size = 12))
                        } else {
                          violin_plot <- VlnPlot(cluster_obj, features= DGE_Top_Markers$gene, log = TRUE) + ggplot2::theme(axis.text.y = element_text(size = 12))
                        }
                        
                        pdf(paste(dir_path,'DGE_',ident1,"_VS_",ident2,'_Within_Cluster_',cluster,'_top',opt$Top_n,'_violinPlot.pdf',sep=''),width = 40, height = 15)
                        print(violin_plot)
                        dev.off()
                         if (!is.na(opt$Subsample)) {
                          heat_map = DoHeatmap(subset(cluster_obj, downsample=subsample), features = DGE_Top_Markers$gene, group.by= DGE_GroupBy,assay=assay2use,slot = slot2use)+ theme(axis.text.y = element_text(size = 12)) + ggplot2::scale_fill_gradientn(colours=colours) 
                         } else {
                          heat_map = DoHeatmap(cluster_obj, features = DGE_Top_Markers$gene, group.by= DGE_GroupBy,assay=assay2use,slot = slot2use) + theme(axis.text.y = element_text(size = 12)) + ggplot2::scale_fill_gradientn(colours=colours) # + NoLegend()
                         } 
                        pdf(paste(dir_path,'DGE_',ident1,"_VS_",ident2,'_Within_Cluster_',cluster,'_top',opt$Top_n,'_heatmap.pdf',sep=''),width = 40, height = 15)
                        print(heat_map)
                        dev.off()
                        
                        # Enrichment Analysis
                        writeLog(logfile, paste("Setting Genes for Enrichment Analysis: min.pct1 >= ",opt$min_pct1," & ALPHA < ",opt$ALPHA,sep=''))
                        df <- DGE_Markers
                        df = df[df$pct.1>=opt$min_pct1 & df$pct.2<opt$max_pct2 & df$p_val_adj<opt$ALPHA,]
                        writeLog(logfile, paste("Contrast ",ident1,' Within Cluster ',cluster,': ',dim(df)[1],"  Differentially Expressed Genes",sep=''))
                    }else{
                        df       = data.frame()
                        writeLog(logfile, paste("FindMarkers Found No Significant Genes"))
                    }
                }
                if (dim(df)[1]>0) {
                    # Significant Genes Heatmap
                    if ('Group' %in% colnames(df)){
                        sig_genes[[paste(cluster,"_UP",sep='')]]   <- df[df$Group=='Up',]$gene
                        sig_genes[[paste(cluster,"_Down",sep='')]] <- df[df$Group=='Down',]$gene
                    } else {
                        sig_genes[[cluster]] <- df$gene
                    }
                    if (!is.na(opt$Subsample)) {
                      heat_map     = DoHeatmap(subset(cluster_obj, downsample=subsample), features = df$gene, group.by= DGE_GroupBy,assay=assay2use,slot = slot2use)+ theme(axis.text.y = element_text(size = 12)) + ggplot2::scale_fill_gradientn(colours=colours) 
                      heat_map_top = DoHeatmap(subset(cluster_obj, downsample=subsample), features = df$gene[1:opt$Top_n], group.by= DGE_GroupBy,assay=assay2use,slot = slot2use)+ theme(axis.text.y = element_text(size = 12)) + ggplot2::scale_fill_gradientn(colours=colours) 
                    } else {
                      heat_map     = DoHeatmap(cluster_obj, features = df$gene, group.by= DGE_GroupBy,assay=assay2use,slot = slot2use) + theme(axis.text.y = element_text(size = 12))+ ggplot2::scale_fill_gradientn(colours=colours)  #+ NoLegend()
                      heat_map_top = DoHeatmap(cluster_obj, features = df$gene[1:opt$Top_n], group.by= DGE_GroupBy,assay=assay2use,slot = slot2use) + theme(axis.text.y = element_text(size = 12))+ ggplot2::scale_fill_gradientn(colours=colours)  #+ NoLegend()
                    }
                    pdf(paste(dir_path,Prefix,cluster,'_heatmap.pdf',sep=''),width = 40, height = 15)
                    print(heat_map)
                    dev.off()
                    
                    pdf(paste(dir_path,Prefix,cluster,"_Top",opt$Top_n,'_heatmap.pdf',sep=''),width = 40, height = 15)
                    print(heat_map_top)
                    dev.off()
                    
                }
            }
            
            if (length(sig_genes)>0){
                # GO Enrichment
                enrichGO_compClusters <- try(clusterProfiler::compareCluster(geneClusters = sig_genes, fun = "enrichGO", pAdjustMethod = opt$pAdjustMethod,
                                                                           OrgDb = OrgDb, ont = opt$ont, universe = background.genes, keyType = "SYMBOL"),
                                           silent = TRUE)
                if (!inherits(enrichGO_compClusters,"try-error")) {
                  enrichGO_compClusters@compareClusterResult$GeneRatio=paste("'",enrichGO_compClusters@compareClusterResult$GeneRatio,sep='')
                  write.csv(enrichGO_compClusters@compareClusterResult,
                          file=paste(opt$outDir,opt$Sample,'/',opt$Sample,"_",Prefix,"_GO_compare.csv",sep=''),quote = FALSE,row.names = TRUE)
                  enrichGO_simplify <- try(clusterProfiler::simplify(enrichGO_compClusters, cutoff = opt$simplify_cutoff),
                                         silent = TRUE)
                  if (!inherits(enrichGO_simplify,"try-error")) {
                    enrichGO_simplify@compareClusterResult$GeneRatio=paste("'",enrichGO_simplify@compareClusterResult$GeneRatio,sep='')
                    write.csv(enrichGO_simplify@compareClusterResult,
                            file=paste(opt$outDir,opt$Sample,'/',opt$Sample,"_",Prefix,"_GO_simplify.csv",sep=''),quote = FALSE,row.names = TRUE)
                    if (dim(enrichGO_simplify@compareClusterResult)[1]<250) {
                      fontsize=14
                    } else {
                      fontsize=14
                    }
                    enrichGO_simplify@compareClusterResult$GeneRatio=stringr::str_remove_all(enrichGO_simplify@compareClusterResult$GeneRatio,"'")
                    enrichGO_dotplot <- enrichplot::dotplot(enrichGO_simplify, showCategory=10000, font.size=fontsize, label_format=100)
                    ggplot2::ggsave(filename = paste(opt$outDir,opt$Sample,'/',opt$Sample,"_",Prefix,"_GO_simplify_dotplot.pdf",sep=''),
                                  enrichGO_dotplot, dpi = 600,device = "pdf",width = 30,height = 50, limitsize = FALSE)
                  }
                }
                
                # save.image(file =paste(opt$outDir,"WorkSpace.RData",sep='') )
                
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
                                             silent = FALSE)
                if (!inherits(enrichKEGG_compClusters,"try-error")) {
                    if (!is.null(enrichKEGG_compClusters)){
                        enrichKEGG_compClusters@compareClusterResult$GeneRatio=paste("'",enrichKEGG_compClusters@compareClusterResult$GeneRatio,sep='')
                        write.csv(enrichKEGG_compClusters@compareClusterResult,
                                file=paste(opt$outDir,opt$Sample,'/',opt$Sample,"_",Prefix,"_KEGG_compareClusters.csv",sep=''),quote = FALSE,row.names = TRUE)
                        if (dim(enrichKEGG_compClusters@compareClusterResult)[1]<250) {
                          fontsize=14
                        } else {
                          fontsize=14
                        }
                        enrichKEGG_compClusters@compareClusterResult$GeneRatio=stringr::str_remove_all(enrichKEGG_compClusters@compareClusterResult$GeneRatio,"'")
                        enrichKEGG_dotplot <- enrichplot::dotplot(enrichKEGG_compClusters, showCategory=10000, font.size=fontsize, label_format=100)
                        ggplot2::ggsave(filename = paste(opt$outDir,opt$Sample,'/',opt$Sample,Prefix,"_KEGG_compareClusters.pdf",sep=''),
                                      enrichKEGG_dotplot, dpi = 600,device = "pdf",width = 30,height = 50, limitsize = FALSE)
                        rm(sig_genes,cluster,file,df,sig_genes_entrez,enrichGO_compClusters,enrichGO_simplify,enrichGO_dotplot,enrichKEGG_compClusters,
                           enrichKEGG_dotplot)
                    }
                }
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
      heat_map = DoHeatmap(subset(obj_seurat, downsample=subsample, group.by= "orig.ident")+ theme(axis.text.y = element_text(size = 12)), features = DGE_Top_Markers$gene) + NoLegend()
    } else {
      heat_map = DoHeatmap(obj_seurat, features = DGE_Top_Markers$gene, group.by= "orig.ident")+ theme(axis.text.y = element_text(size = 12)) + NoLegend()
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
          heat_map = DoHeatmap(subset(obj_seurat, downsample=subsample), features = sig_genes, group.by= "orig.ident") + theme(axis.text.y = element_text(size = 12)) + NoLegend()
        } else {
          heat_map = DoHeatmap(obj_seurat, features = sig_genes, group.by= "orig.ident")+ theme(axis.text.y = element_text(size = 12)) + NoLegend()
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
                fontsize=14
              } else {
                fontsize=14
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
        if (!is.null(enrichKEGG_compClusters)){
            enrichKEGG_compClusters@compareClusterResult$GeneRatio=paste("'",enrichKEGG_compClusters@compareClusterResult$GeneRatio,sep='')
            write.csv(enrichKEGG_compClusters@compareClusterResult,
                      file=paste(opt$outDir,opt$Sample,'/Cluster_',cluster1,'/Cluster',cluster1,'_KEGG_compareClusters_vs_each_other.csv',sep=''),quote = FALSE,row.names = TRUE)
            if (dim(enrichKEGG_compClusters@compareClusterResult)[1]<250) {
              fontsize=12
            } else {
              fontsize=10
            }
            enrichKEGG_compClusters@compareClusterResult$GeneRatio=stringr::str_remove_all(enrichKEGG_compClusters@compareClusterResult$GeneRatio,"'")
            enrichKEGG_dotplot <- enrichplot::dotplot(enrichKEGG_compClusters, showCategory=10000, font.size=fontsize, label_format=100)
            ggplot2::ggsave(filename = paste(opt$outDir,opt$Sample,'/Cluster_',cluster1,'/Cluster',cluster1,'_KEGG_compareClusters_vs_each_other_dotplot.pdf',sep=''),
                            enrichKEGG_dotplot, dpi = 600,device = "pdf",width = 30,height = 50, limitsize = FALSE)
        }
      }
    }

