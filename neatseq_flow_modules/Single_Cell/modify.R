.libPaths()
library(optparse)
library(Seurat)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--inputRDS"), type="character", default = NA,
              help="Path to Seurat object RDS file", metavar = "character"),
  make_option(c("-s", "--Sample"), type="character", default = NA,
              help="Sample's name", metavar = "character"),
  make_option(c("-o", "--outDir"), type="character", default = NA,
              help="Path to the output directory", metavar = "character"),
  make_option(c("--DietSeurat"), action="store_true", default = FALSE,
              help="Run DietSeurat (keep counts and data)", metavar = "character"),
  make_option(c("--Keep_reducs"), type="character", default = NA,
              help="Which reductions to keep if running DietSeurat (Default is None)", metavar = "character"),
  make_option(c("--Keep_assays"), type="character", default = NA,
              help="Which assays to keep if running DietSeurat (Default is active assay)", metavar = "character"),
  make_option(c("--Keep_graphs"), type="character", default = NA,
              help="Which graphs to keep if running DietSeurat (Default is None)", metavar = "character"),
  make_option(c("--Subset_Clusters"), type="character", default = NA,
              help="Vector of clusters to keep (comma-separated)", metavar = "character"),
  make_option(c("--Remove_Clusters"), type="character", default = NA,
              help="Vector of clusters to remove (comma-separated)", metavar = "character"),
  make_option(c("--Set_Idents"), type="character", default = NA,
              help="Path to file with old and new identities (two-column matrix)", metavar = "character"),
  make_option(c("--Set_Ident_by_Marker"), type="character", default = NA,
              help="Set identity of cells based on marker expression using template: <MARKER>|<CUTOFF>|<Idents>|<New_Ident>", metavar = "character"),
  make_option(c("--Active_Assay"), type="character", default = "RNA",
              help="Set the Active_Assay for this analysis, default = RNA" , metavar = "character"),
  make_option(c("--CleanClusters"), action="store_true", default = FALSE,
              help="Remove clustering information from meta-data (Default is False)", metavar = "character"),
  make_option(c("--Pseudobulk"), action="store_true", default = FALSE,
              help="Create pseudobulk object from each sample (Default is False)", metavar = "character"),
  make_option(c("--group.by1"), type="character", default = "orig.ident",
              help="Create pseudobulk object from each sample and a second variable, i.e. cluster", metavar = "character"),
  make_option(c("--group.by2"), type="character", default = "seurat_clusters",
              help="Create pseudobulk object from each sample and a second variable, i.e. cluster", metavar = "character"),
  make_option(c("--Cell_cycle"), , action="store_true", default = FALSE,
              help="examine cell cycle variation in our data", metavar = "character"),
  make_option(c("--Cell_cycle_genes"), type="character", default = NA,
              help="Path to csv file with cel cycle associated genes", metavar = "character"),			  
  #mouse genes from https://github.com/hbc/tinyatlas/tree/master/cell_cycle
  make_option(c("--Remove_Doublets"), type="character", default = NA,
              help="Path to CSV file with classifications for each cell (DF.classification_homotypic column)", metavar = "character")
);

#Get user information
opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", 
                                    option_list=option_list,
                                    epilogue="\n\nAuthor:Gil Sorek");
opt = optparse::parse_args(opt_parser);


Doublets_colname = 'DF.classification_homotypic'

# Set log framework
logfile = paste(opt$outDir,opt$Sample,'_log.txt',sep='')
cat(paste('[',Sys.time(),']: Modify: ',opt$Sample,sep=''), file=logfile, sep='\n')
writeLog <- function(logfile, msg) { cat(paste('(',strftime(Sys.time(), format = "%H:%M:%S"),'): ',msg,sep=''),file=logfile,append=TRUE,sep='\n') }
saveRDS(opt, paste(opt$outDir,'opt.rds',sep=''))

print("Input var:",quote = F)
print.data.frame(as.data.frame(x = unlist(opt),row.names = names(opt)),right = F,quote = F)

# Import data
if (!is.na(opt$inputRDS)) {
  writeLog(logfile, paste("Importing Seurat RDS object...",sep=''))
  obj_seurat <- readRDS(opt$inputRDS)
  if ( opt$Active_Assay  %in% Seurat::Assays(obj_seurat)){
      Seurat::DefaultAssay(obj_seurat) <- opt$Active_Assay
    }else{
      writeLog(logfile, paste("ERROR: Unable to find Assay '",opt$Active_Assay,"'. Available Assays: ",paste(names(Seurat::Assays(obj_seurat)),collapse=','),sep=''))
      stop()
    }
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
if (!is.na(opt$Subset_Clusters)) {
  writeLog(logfile, paste("Subsetting clusters..."))
  clusters2keep <- unlist(stringi::stri_split(str = opt$Subset_Clusters,fixed = ','))
  if (isFALSE(all(clusters2keep %in% Idents(obj_seurat)))) {
    writeLog(logfile, paste("ERROR: Cannot find the following identities in the object: ",
								paste(clusters2keep[which(!(clusters2keep %in% Idents(obj_seurat)))],collapse=','),sep=''))
	writeLog(logfile, paste("The available idents are: ",paste(unique(Idents(obj_seurat)),collapse=','),sep=''))
    stop()
  }
  obj_seurat = subset(obj_seurat, idents = clusters2keep, invert = FALSE)
  writeLog(logfile, paste("Subsetting the following identities from the object: ",opt$Subset_Clusters,sep=''))
}
if (!is.na(opt$Remove_Clusters)) {
  writeLog(logfile, paste("Removing clusters..."))
  clusters2remove <- unlist(stringi::stri_split(str = opt$Remove_Clusters,fixed = ','))
  if (isFALSE(all(clusters2remove %in% Idents(obj_seurat)))) {
    writeLog(logfile, paste("ERROR: Cannot find the following identities in the object: ",
								paste(clusters2remove[which(!(clusters2remove %in% Idents(obj_seurat)))],collapse=','),sep=''))
	writeLog(logfile, paste("The available idents are: ",paste(unique(Idents(obj_seurat)),collapse=','),sep=''))
    stop()
  }
  obj_seurat = subset(obj_seurat, idents = clusters2remove, invert = TRUE)
  writeLog(logfile, paste("Removing the following identities from the object: ",opt$Remove_Clusters,sep=''))
}



# Remove Doublets
if (!is.na(opt$Remove_Doublets)) {
  writeLog(logfile, paste("Removing Doublets..."))
  files = unlist(stringi::stri_split(str = opt$Remove_Doublets,fixed = ','))
  Allcells2remove = c()
  for (file in files){
    if (file.exists(file)) {
        DF_classifications <- read.csv(file, header = TRUE, row.names = 1)
        if ( Doublets_colname %in% colnames(DF_classifications)) {
            Doublets_val_colname = colnames(DF_classifications)[which( Doublets_colname == colnames(DF_classifications))-1]
            temp_obj_seurat <- SeuratObject::AddMetaData(obj_seurat,DF_classifications)
            temp_obj_seurat@meta.data[which(is.na(temp_obj_seurat@meta.data[Doublets_colname])),Doublets_colname] = "Singlet"
            temp_obj_seurat@meta.data[temp_obj_seurat@meta.data[Doublets_colname]!="Doublet",Doublets_colname]    = "Singlet"
            ct_plt <- DimPlot(temp_obj_seurat, group.by = Doublets_colname)
            pdf(file = paste(opt$outDir,opt$Sample,tail(stringi::stri_split_fixed(str=file,pattern='/')[[1]],2)[1],'_DoubletsPlot.pdf',sep=''), 
                width = 25, height = 15)
            print(ct_plt)
            dev.off()
            
            temp_obj_seurat@meta.data[temp_obj_seurat@meta.data[Doublets_colname]!="Doublet",Doublets_val_colname] = 0
            ct_plt <- FeaturePlot(temp_obj_seurat, features = Doublets_val_colname)
            pdf(file = paste(opt$outDir,opt$Sample,tail(stringi::stri_split_fixed(str=file,pattern='/')[[1]],2)[1],'_Doublets_Value_Plot.pdf',sep=''), 
                width = 25, height = 15)
            print(ct_plt)
            dev.off()
            cells2remove <- row.names(DF_classifications[which(DF_classifications[Doublets_colname]=="Doublet"),])
            writeLog(logfile, paste("Removed ",length(cells2remove)," Doublets",sep=''))
            Allcells2remove = c(Allcells2remove,cells2remove)
        } else {
            writeLog(logfile, paste("ERROR: Missing column",Doublets_colname))
            stop()
        }
    } else {
        writeLog(logfile, paste("ERROR: ",file," does not exists.",sep=''))
        stop()
    }
 }
 if (length(Allcells2remove)>0){
    obj_seurat@meta.data[Doublets_colname] = "Singlet"
    obj_seurat@meta.data[rownames(obj_seurat@meta.data) %in% Allcells2remove,Doublets_colname] = "Doublet"
    ct_plt <- DimPlot(obj_seurat, group.by = Doublets_colname)
    pdf(file = paste(opt$outDir,opt$Sample,'_ALL_DoubletsPlot.pdf',sep=''), 
       width = 25, height = 15)
    print(ct_plt)
    dev.off()
    obj_seurat <- subset(obj_seurat, cells = unique(Allcells2remove), invert = TRUE)
    saveRDS(Allcells2remove, file = paste(opt$outDir, opt$Sample, '_DoubletsCells2Remove.rds', sep=''))
 }
}
# Rename Identities
if (sum(!is.na(opt$Set_Idents),!is.na(opt$Set_Ident_by_Marker))>1) {
  writeLog(logfile, paste("ERROR: Cannot rename identities using two methods, Please use either --Rename_Idents / --Set_Ident_by_Marker",sep=''))
  stop()
}
# Rename Idents based on meta-data matrix
if (!is.na(opt$Set_Idents)) {
  writeLog(logfile, paste("Rename Idents using meta-data matrix"))
  if (!file.exists(opt$Set_Idents)) {
    writeLog(logfile, paste("ERROR: Matrix file does not exists"))
	stop()
  }
  data <- read.table(opt$Set_Idents, header=TRUE)
  if (isFALSE(all(colnames(data) == c('ident.old','ident.new')))) {
    writeLog(logfile, paste("ERROR: Invalid column names. Please use ident.old and ident.new"))
	stop()
  }
  writeLog(logfile, paste("Loaded idents matrix, will rename ",length(data$ident.old)," identities",sep=''))
  idents2rename = data$ident.old
  if (isFALSE(all(idents2rename %in% Idents(obj_seurat)))) {
    writeLog(logfile, paste("ERROR: Cannot find the following identities in the object: ",
								paste(idents2rename[which(!(idents2rename %in% Idents(obj_seurat)))],collapse=','),sep=''))
	writeLog(logfile, paste("The available idents are: ",paste(unique(Idents(obj_seurat)),collapse=','),sep=''))
    stop()
  }
  fetch_idents = Idents(obj_seurat)
  for (i in 1:nrow(data)) {
	ident.old = paste(data$ident.old[i])
	ident.new = paste(data$ident.new[i])
	cells2rename = names(fetch_idents[fetch_idents==ident.old])
    Idents(obj_seurat, cells = cells2rename) = ident.new
	obj_seurat@meta.data$seurat_clusters <- Idents(obj_seurat)
  }
}
# Rename Idents based on Marker Expression (--Set_Ident_by_Marker)
if (!is.na(opt$Set_Ident_by_Marker)) {
  writeLog(logfile, paste("Rename Cells based on Marker Expression: ",opt$Set_Ident_by_Marker,sep=''))
  pattern_vector = unlist(stringi::stri_split(str = opt$Set_Ident_by_Marker,fixed = '|'))
  if (length(pattern_vector)!=4) {
    writeLog(logfile, paste("ERROR: --Set_Ident_by_Marker Invalid argument. Please use template: <MARKER>|<CUTOFF>|<Idents>|<New_Ident>",sep=''))
	stop()
  }
  marker = pattern_vector[1]
  if (marker %in% rownames(obj_seurat)) {
    writeLog(logfile, paste("Detected marker: ",marker,sep=''))
  } else {
    writeLog(logfile, paste("ERROR: Unable to find marker '",marker,"'",sep=''))
	writeLog(logfile, rownames(obj_seurat))
	stop()
  }
  threshold = pattern_vector[2]
  writeLog(logfile, paste("Rename cells with '",marker,"' expression above ",threshold,sep=''))
  old_idents = pattern_vector[3]
  new_ident = pattern_vector[4]
  idents2rename = unlist(stringi::stri_split(str = old_idents,fixed = ','))
  if (isFALSE(all(idents2rename %in% Idents(obj_seurat)))) {
    writeLog(logfile, paste("ERROR: Cannot find the following identities in the object: ",
								paste(idents2rename[which(!(idents2rename %in% Idents(obj_seurat)))],collapse=','),sep=''))
	writeLog(logfile, paste("The available idents are: ",paste(unique(Idents(obj_seurat)),collapse=','),sep=''))
    stop()
  }
  writeLog(logfile, paste("Rename cells to '",new_ident,"' within idents: ",old_idents,sep=''))
  fetch_data = FetchData(obj_seurat, vars = c('ident',marker))
  cells2rename = rownames(fetch_data[fetch_data$ident %in% idents2rename & fetch_data[[marker]]>threshold,])
  if (length(cells2rename)==0) {
    writeLog(logfile, paste("ERROR: No cells within idents '",old_idents,
	                        "' express '",marker,"' above ",threshold," Max:",max(fetch_data[fetch_data$ident %in% idents2rename,marker]),
	                        " Min:",min(fetch_data[fetch_data$ident %in% idents2rename,marker]),sep=''))
	stop()
  }
  writeLog(logfile, paste("Found ",length(cells2rename)," cells, renaming to '",new_ident,"'",sep=''))
  Idents(obj_seurat, cells=cells2rename) <- new_ident
  obj_seurat@meta.data$seurat_clusters <- Idents(obj_seurat)
}

# DietSeurat
if (opt$DietSeurat) {
  writeLog(logfile, paste("Running DietSeurat..."))
  if (is.na(opt$Keep_assays)) {
    writeLog(logfile, paste("DietSeurat: Active assay (",obj_seurat@active.assay,") will be kept",sep=''))
	keep_assay = NULL
  } else {
    if (opt$Keep_assays %in% names(obj_seurat@assays)) {
	  writeLog(logfile, paste("DietSeurat: '",opt$Keep_assays,"' assay will be kept",sep=''))
	  keep_assay = opt$Keep_assays
	} else {
	  writeLog(logfile, paste("ERROR: Unable to find assay '",opt$Keep_assays,"'. Available assays: ",paste(names(obj_seurat@assays),collapse=','),sep=''))
	  stop()
	}
  }
  if (is.na(opt$Keep_reducs)) {
    writeLog(logfile, paste("DietSeurat: No reductions will be kept",sep=''))
	keep_reducs = NULL
  } else {
    if (opt$Keep_reducs %in% names(obj_seurat@reductions)) {
	  writeLog(logfile, paste("DietSeurat: '",opt$Keep_reducs,"' reduction will be kept",sep=''))
	  keep_reducs = opt$Keep_reducs
	} else {
	  writeLog(logfile, paste("ERROR: Unable to find reduction '",opt$Keep_reducs,"'. Available reductions: ",paste(names(obj_seurat@reductions),collapse=','),sep=''))
	  stop()
	}
  }
  if (is.na(opt$Keep_graphs)) {
    writeLog(logfile, paste("DietSeurat: No graphs will be kept",sep=''))
	keep_graphs = NULL
  } else {
    if (opt$Keep_graphs %in% names(obj_seurat@graphs)) {
	  writeLog(logfile, paste("DietSeurat: '",opt$Keep_graphs,"' graph will be kept",sep=''))
	  keep_graphs = opt$Keep_graphs
	} else {
	  writeLog(logfile, paste("ERROR: Unable to find graph '",opt$Keep_graphs,"'. Available graphs: ",paste(names(obj_seurat@graphs),collapse=','),sep=''))
	  stop()
	}
  }
  obj_seurat <- DietSeurat(obj_seurat, counts = TRUE, data = TRUE, assays = keep_assay, dimreducs = keep_reducs, graphs = keep_graphs)
}

# Clean Clusters
if (opt$CleanClusters) {
  writeLog(logfile, paste("Cleaning previous clusters information..."))
  obj_seurat@meta.data[,stringr::str_detect(colnames(obj_seurat@meta.data),c('snn_res','seurat_clusters'))] <- NULL
}

# Save results
saveRDS(obj_seurat, file = paste(opt$outDir, opt$Sample, '.rds', sep=''))
writeLog(logfile, paste("Finished"))

#Create pseudobulk from Seurat object:
if (opt$Pseudobulk) {
  writeLog(logfile, paste("Creating pseudobulk object per group..."))
  
  library(tidyr)
  library(reshape2)
 #Use active.ident as cellType in metadata:
  df1 <- as.data.frame(Idents(obj_seurat))
  obj_seurat@meta.data$cellType <- df1[,1]
  counts <- GetAssayData(object = obj_seurat, assay = "RNA", layer = "counts")
  
  count <- as.data.frame(t(as.matrix(counts)))
		if(all(rownames(count)== rownames(obj_seurat@meta.data))){
  count$group1 <-  obj_seurat@meta.data[ ,names(obj_seurat@meta.data) %in% opt$group.by1]
  count$group2 <- obj_seurat@meta.data[ ,names(obj_seurat@meta.data) %in% opt$group.by2]
  count$group1_2 <- paste(count$group1, count$group2, sep ="_")  
  ps <- gather(count, "Gene","gene_sum", 1:(ncol(count)-3))
  pseudo <- dcast(ps, Gene ~ group1_2, fun.aggregate = sum, value.var = 'gene_sum')
  rownames(pseudo) <- pseudo$Gene
#Save results
write.table(pseudo[,-1],file = paste(opt$outDir, opt$Sample, '.txt', sep=''), sep = "\t", quote= FALSE)
}
}

if (opt$Cell_cycle) {
  writeLog(logfile, paste("Storing S and G2/M scores in object meta data..."))
if (opt$  Cell_cycle_genes){
	writeLog(logfile, paste("Using S and G2/M gene list from file"))
	}
#Assign Cell-Cycle Scores:
#This stores S and G2/M scores in object meta data, along with the predicted classification of each cell
#in either G2M, S or G1 phase
Cell_cycle_genes <- read.csv(opt$Cell_cycle_genes)
s.genes <-  Cell_cycle_genes[Cell_cycle_genes$phase == 'S',]
g2m.genes <- Cell_cycle_genes[Cell_cycle_genes$phase == 'G2/M',]
obj_seurat <- CellCycleScoring(obj_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
# Save results
saveRDS(obj_seurat, file = paste(opt$outDir, opt$Sample, '.rds', sep=''))
writeLog(logfile, paste("Finished"))
}