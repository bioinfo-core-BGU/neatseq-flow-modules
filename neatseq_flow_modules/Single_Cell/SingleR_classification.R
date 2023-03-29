.libPaths()
library(optparse)
#library(Seurat)
library(dplyr)
library(SingleR)
library(celldex)
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--inputRDS"), type="character", default = NA,
              help="Path to Seurat object RDS file", metavar = "character"),
  make_option(c("-s", "--Sample"), type="character", default = NA,
              help="Sample's name", metavar = "character"),
  make_option(c("-o", "--outDir"), type="character", default = NA,
              help="Path to the output directory", metavar = "character")
);

#Get user information
opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", 
                                    option_list=option_list,
                                    epilogue="\n\nAuthor:Gil Sorek");
opt = optparse::parse_args(opt_parser);

print("Input var:",quote = F)
print.data.frame(as.data.frame(x = unlist(opt),row.names = names(opt)),right = F,quote = F)

# Import data
if (!is.na(opt$inputRDS)) {
  obj_seurat <- readRDS(opt$inputRDS)
} else {
  stop("Seurat object RDS file must be specified [--inputRDS]")
}
if (is.na(opt$Sample)) {
  stop("Samples name must be specified [--Sample]")
}
if (is.na(opt$outDir)) {
  stop("Output directory must be specified [--outDir]")
}

# Prepare counts data
counts <- as.matrix(Seurat::GetAssayData(obj_seurat, slot = "counts"))

# Human Primary Cell Atlas classification
hpca.se <- HumanPrimaryCellAtlasData()
pred.hesc <- SingleR(test = counts, ref = hpca.se, assay.type.test=1, labels = hpca.se$label.main)
hpca.preds <- data.frame(row.names =row.names(pred.hesc),cell_Type.hesc= pred.hesc$pruned.labels)

# BluePrint classification
bp.se <- BlueprintEncodeData()
pred.Blueprint <- SingleR(test = counts, ref = bp.se, assay.type.test=1, labels = bp.se$label.main)
bp.preds <- data.frame(row.names =row.names(pred.Blueprint),cell_Type.Blueprint= pred.Blueprint$pruned.labels)
