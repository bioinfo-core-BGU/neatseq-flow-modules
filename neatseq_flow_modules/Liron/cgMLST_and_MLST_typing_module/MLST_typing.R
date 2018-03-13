# AUTHOR: Menachem Sklarz & Michal Gordon



# Check if required packages are installed:
if(!(all(c("magrittr") %in% installed.packages()))) {
    if (Sys.getenv("CONDA_PREFIX")!=""){
        install.packages("magrittr", dependencies=TRUE,repos = "http://cran.us.r-project.org")
    } else {
        cat("The 'magrittr' package is not installed. You must install it for this script to work!")
        }
}
library("magrittr")
if(!(all(c("plyr") %in% installed.packages()))) {
    if (Sys.getenv("CONDA_PREFIX")!=""){
        install.packages("plyr", dependencies=TRUE,repos = "http://cran.us.r-project.org")
    } else {
        cat("The 'plyr' package is not installed. You must install it for this script to work!")
        }
    
    
}
library("plyr")
if(!(all(c("optparse") %in% installed.packages()))) {
    if (Sys.getenv("CONDA_PREFIX")!=""){
        install.packages("optparse", dependencies=TRUE,repos = "http://cran.us.r-project.org")
    } else {
        cat("The 'optparse' package is not installed. You must install it for this script to work!")
        }
}
library("optparse")
library(tools)


args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-b", "--blast"), type="character", default=NULL, 
              help="Path to PARSED blast results", metavar="character"),
  make_option(c("-s", "--scheme"), type="character", default=NULL, 
              help="Path to tab-delimited MLST scheme WITH TITLE LINE.", metavar="character"),
  make_option(c("-n", "--num_of_genes"), type="numeric", default=7, 
              help="Number of genes  (default: 7)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Path to output file (default=<blast_input>.MLST)", metavar="character"),
  make_option(c("-c", "--Type_col_name"), type="character", default = NULL, 
              help="Columns in the scheme that are not locus names", metavar="character"),
  make_option(c("-F", "--Find_close_match"), type="character", default = "N", 
              help="whether to show the closest allele combination match? (Y/N)", metavar="character"),
  make_option(c("-I", "--ignore_missing_genes"), type="character", default = "N", 
              help="If allele combination was not found: Try ignoring missing genes when searching allele combination match? (Y/N)", metavar="character")
); 



opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", 
                                    option_list=option_list,
                                    epilogue="\n\nAuthor: Menachem Sklarz");
opt = optparse::parse_args(opt_parser);

opt$Type_col_name <- strsplit(x = opt$Type_col_name,
                               split = ",") %>% unlist

# REad mlst scheme
mlst <- read.table(opt$scheme,
                   sep="\t",
                   header  = T,
                   stringsAsFactors = F)
                   
ignore_missing_genes=FALSE
gene_names = names(mlst)

for (i in opt$Type_col_name){
    gene_names=gene_names[gene_names != i]
    } 
opt$num_of_genes=length(gene_names)

#seting the default in status new column
mlst[dim(mlst)[1]+1,"Status"]=''
mlst["Status"]="OK"
mlst["Percentage_of_missing_genes"]=0
#adding sample name

sample<- strsplit(x = opt$blast,
              split = "/")%>% unlist %>% tail(.,n=1) %>% strsplit(.,split=".blast.parsed") %>% unlist
mlst["Samples"]=sample


# Read BLAST data:
blast_raw <- read.delim(opt$blast,
                        he = T,
                        stringsAsFactors = F)
# Set Gene as rowname
rownames(blast_raw) <- blast_raw$Gene
blast_raw<-blast_raw[intersect(blast_raw$Gene,gene_names),]

if (length(blast_raw$Gene) < opt$num_of_genes){
 
  #sending erorr
  missing_genes <- setdiff(gene_names,blast_raw$Gene)
  write(sprintf("Genes %s not found",
                paste(missing_genes,
                      collapse = " ")), 
        file = stderr())
  
  mlst["Status"]=sprintf("Genes %s not found",
                         paste(missing_genes,
                               collapse = " "))
  mlst["Percentage_of_missing_genes"]=length(missing_genes)/opt$num_of_genes
  
  blast_raw[missing_genes,"Gene"]=missing_genes
  blast_raw[missing_genes,"pident"]=100
  blast_raw[missing_genes,"coverage"]=100
  blast_raw[missing_genes,"Number"]='N'
  
  # # Removing non-perfectly matching alleles
  # blast <- blast_raw[(blast_raw$pident==100 & blast_raw$coverage==100),]
  # # Finding genes in scheme not found perfectly in sample
  # non_existing_allel <- setdiff(blast_raw$Gene,blast$Gene)
  # if (length(non_existing_allel) > 0) {
    
    # write(sprintf("Genes %s not found perfectly",
                  # paste(non_existing_allel,
                        # collapse = " ")), 
          # file = stderr())
    # #informing the not perfect match in the status column 
    # mlst["Status"]=paste( mlst[dim(mlst)[1],"Status"] , sprintf("Genes %s not found perfectly",
                           # paste(non_existing_allel,
                                 # collapse = " ")),sep=" / ")
  # }
  # if (opt$Find_close_match=='Y'){
      # #adding new line in the scheme
      # mlst[dim(mlst)[1],blast_raw$Gene]=blast_raw[blast_raw$Gene,"Number"]
  # }else{
    # #adding new line in the scheme
    # mlst[dim(mlst)[1],blast$Gene]=blast[blast$Gene,"Number"]
    # #adding new alleles to the scheme
    # mlst[dim(mlst)[1],non_existing_allel]=sapply(X = blast_raw[non_existing_allel,"sseq"]  , FUN = function(x) paste("New_Allele=",stringi::stri_replace_all(str = x,replacement ='',regex = '-') ,sep = "")) 
    # }
  # write.table(x         = mlst[dim(mlst)[1],,drop=F] , 
              # file      = opt$output,
              # sep       = "\t",
              # row.names = F,
              # quote     = F,
              # col.names = T)
  
  } 
if (TRUE) {
    # Removing non-perfectly matching alleles
    blast <- blast_raw[(blast_raw$pident==100 & blast_raw$coverage==100),]
    
    
    # Finding genes in scheme not found perfectly in sample
    non_existing_allel <- setdiff(gene_names,blast$Gene)
    if (length(non_existing_allel) > 0) {
 
          write(sprintf("Genes %s not found perfectly",
                        paste(non_existing_allel,
                              collapse = " ")), 
                file = stderr())
        #informing the not perfect match in the status column 
        mlst["Status"]=sprintf("Genes %s not found perfectly",
                               paste(non_existing_allel,
                                     collapse = " "))
        
        if (opt$Find_close_match=='Y'){
            #going back before the identity cutoff 
            blast <- blast_raw
        }else{
          #going back before the identity cutoff 
          blast <- blast_raw
          #adding new alleles to the scheme
          blast[non_existing_allel,"Number"]=sapply(X = blast_raw[non_existing_allel,"sseq"]  , FUN = function(x) paste("New_Allele=",stringi::stri_replace_all(str = x,replacement ='',regex = '-') ,sep = "")) 
        }
    
    }


    
    if ((dim(mlst)[1]-1)>0){
        MLST_num <- sapply(X      = 1:(dim(mlst)[1]-1), 
                          FUN    = function(x) all(mlst[x,gene_names]  ==blast[gene_names,"Number"])  ) 
        
        if (is.na(sum(MLST_num))){
            MLST_num=0
        }
        if (sum(MLST_num) == 0) {
            temp=blast[gene_names,"Number"]
            temp[temp=='N']='0'
            MLST_num <- sapply(X      = 1:(dim(mlst)[1]-1), 
                          FUN    = function(x) all(mlst[x,gene_names]  ==temp)  ) 
            if (is.na(sum(MLST_num))){
                MLST_num=0
            }
        }
        if (sum(MLST_num) == 0) {
            temp=mlst
            temp[is.na(temp)]='N'
            MLST_num <- sapply(X      = 1:(dim(temp)[1]-1), 
                          FUN    = function(x) all(temp[x,gene_names]  == blast[gene_names,"Number"])  ) 
            if (is.na(sum(MLST_num))){
                MLST_num=0
            }
        }
        
        if (sum(MLST_num) == 0) {
            if (opt$ignore_missing_genes=='Y'){
                temp=blast[gene_names,"Number"]
                MLST_num <- sapply(X      = 1:(dim(mlst)[1]-1), 
                              FUN    = function(x) all(mlst[x,gene_names[temp!='N']]  ==blast[gene_names[temp!='N'],"Number"])  ) 
                if (is.na(sum(MLST_num))){
                    MLST_num=0
                }else{
                    ignore_missing_genes=TRUE
                }
            }
        }
    }else{
        MLST_num=0
    }
    
    

    
    if(sum(MLST_num) == 0) {
      
      write("Allele combination not found in scheme", file = stderr())
      mlst[dim(mlst)[1],blast$Gene]=blast[blast$Gene,"Number"]
      #adding new alleles to the scheme
      mlst[dim(mlst)[1],non_existing_allel]=sapply(X = blast_raw[non_existing_allel,"sseq"]  , FUN = function(x) paste("New_Allele=",stringi::stri_replace_all(str = x,replacement ='',regex = '-') ,sep = "")) 
      if (mlst[dim(mlst)[1],"Status"]=="OK") {
          mlst["Status"]="Allele combination not found in scheme"
      } else {
          mlst["Status"]=paste( mlst[dim(mlst)[1],"Status"] , sprintf("Allele combination not found in scheme"),
                                sep=" / ")                                                                      
      }
      
      
      mlst["Samples"]=sample
      write.table(x         = mlst[dim(mlst)[1],,drop=F] , 
                  file      = opt$output,
                  sep       = "\t",
                  row.names = F,
                  quote     = F,
                  col.names = T)
      
    } else if(sum(MLST_num) > 1) {
      write("Allele combination appears more than once in scheme", file = stderr())
      mlst[dim(mlst)[1],blast$Gene]=blast[blast$Gene,"Number"]
      if (mlst[dim(mlst)[1],"Status"]=="OK") {
        mlst["Status"]="Allele combination appears more than once in scheme"
      }else {
        mlst["Status"]=paste( mlst[dim(mlst)[1],"Status"] , sprintf("Allele combination appears more than once in scheme"),
                              sep=" / ")                                                                      
      }
      
      mlst["Samples"]=sample
      write.table(x         = mlst[dim(mlst)[1],,drop=F] , 
                  file      = opt$output,
                  sep       = "\t",
                  row.names = F,
                  quote     = F,
                  col.names = T)
    } else {
      if (mlst[dim(mlst)[1],"Status"]!="OK") {
        if (ignore_missing_genes){
            write("Found the closest allele combination ignoring missing genes", file = stderr())
            mlst["Status"]=paste( mlst[dim(mlst)[1],"Status"] , sprintf("Found the closest allele combination ignoring missing genes"),
                                  sep=" / ")  
        }
      }
      write.table(x         = mlst[which(MLST_num),,drop=F] , 
                  file      = opt$output,
                  sep       = "\t",
                  row.names = F,
                  quote     = F,
                  col.names = T)
    }
}