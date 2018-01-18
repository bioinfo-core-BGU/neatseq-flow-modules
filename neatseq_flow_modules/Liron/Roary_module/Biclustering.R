# Check if required packages are installed:
if(!(all(c("optparse") %in% installed.packages()))) {
    if (Sys.getenv("CONDA_PREFIX")!=""){
        install.packages("optparse", dependencies=TRUE,repos = "http://cran.us.r-project.org")
    } else {
        cat("The 'optparse' package is not installed. You must install it for this script to work!")
        }
}
library("optparse")



args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-a", "--Annotation"), type="character", default=NULL, 
              help="Path to VFDB Annotation file", metavar="character"), 
  make_option(c("-r", "--Roary_Results"), type="character", default=NULL, 
              help="Path to the Roary presence-absence csv Results file", metavar="character"), 
  make_option(c("-d", "--metadata"), type="character", default=NULL, 
              help="Path to a tabular MetaData file", metavar="character"), 
  make_option(c("-f", "--ID_field"), type="character", default="Samples", 
              help="Column name in the MetaData file of the Samples IDs corresponding to the Results file Samples IDs [default='Samples']", metavar="character"), 
  make_option(c("--Minimal_samples"), type="numeric", default = 2, 
              help="Minimal number of samples in a bicluster [default = 2]", metavar="character"), 
  make_option(c( "--Minimal_genes"), type="numeric", default = 2, 
              help="Minimal number of genes in a bicluster [default = 2]", metavar="character"), 
  make_option(c("-m", "--filter_min"), type="numeric", default = 0.05, 
              help="Filter genes with frequency lower than [default = 0.05]", metavar="character"), 
  make_option(c("-u", "--filter_max"), type="numeric", default = 0.95, 
              help="Filter genes with frequency higher than [default = 0.95]", metavar="character"), 
  make_option(c("-P", "--pvalueCutoff"), type="numeric", default = 0.05, 
              help="Enrichment pvalue Cutoff [default = 0.05]", metavar="character"), 
  make_option(c("-p", "--pAdjustMethod"), type="character", default='fdr', 
              help="Enrichment Multiple correction method [holm, hochberg, hommel, bonferroni, BH, BY, fdr, none] [default='fdr']", metavar="character"), 
  make_option(c("-g", "--thrgene"), type="numeric", default=3 , 
              help="The threshold parameters for the ISA, for features [default=3]", metavar="character"), 
  make_option(c("-c", "--thrcond"), type="numeric", default=2.75 , 
              help="The threshold parameters for the ISA, for samples [default=2.75]", metavar="character"), 
  make_option(c("-s", "--noseeds"), type="numeric", default=100 , 
              help="Number of seeds to run ISA from [default=100]", metavar="character"), 
  make_option(c("-o", "--output"), type="character", default=getwd(), 
              help="Path to output directory", metavar="character"), 
  make_option(c("-t", "--cols_to_use"), type="character", default = NULL, 
              help="Columns in the MetaData file to use", metavar="character")
); 



opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", 
                                    option_list=option_list,
                                    epilogue="\n\nAuthor:Liron Levin based on Eliad Levi script");
opt = optparse::parse_args(opt_parser);


print("Input var:",,quote = F)
print.data.frame(as.data.frame(x = unlist(opt),row.names = names(opt)),right = F,quote = F)
Annotation_file=opt$Annotation
Roary_Results_file=opt$Roary_Results
metadata_file=opt$metadata
pAdjustMethod=opt$pAdjustMethod
pvalueCutoff=opt$pvalueCutoff
metadata_sample_field=opt$ID_field
min.bicluster.samples =opt$Minimal_samples
min.bicluster.genes =opt$Minimal_genes
Roary_filter_min=opt$filter_min
Roary_filter_max=opt$filter_max
metadata_select_fileds=opt$cols_to_use
thr.gene =opt$thrgene
thr.cond =opt$thrcond
no.seeds =opt$noseeds
output=opt$output

if (is.null(Roary_Results_file)){
    print("No Roary file!!!!")
}else if (is.null(metadata_file)){
    print("No MetaData file!!!!")
} else if (is.null(Annotation_file)) {
    print("No Annotation file!!!!")
}else{



    # Check if required packages are installed:
    if(!(all(c("hgu95av2.db") %in% installed.packages()))) {
        if (Sys.getenv("CONDA_PREFIX")!=""){
            source("https://bioconductor.org/biocLite.R")
            biocLite("hgu95av2.db")
        } else {
            cat("The 'hgu95av2.db' package is not installed. You must install it for this script to work!")
        }
    }
    library("hgu95av2.db")
    if(!(all(c("clusterProfiler") %in% installed.packages()))) {
        if (Sys.getenv("CONDA_PREFIX")!=""){
            source("https://bioconductor.org/biocLite.R")
            biocLite("clusterProfiler")
        } else {
            cat("The 'clusterProfiler' package is not installed. You must install it for this script to work!")
        }
    }
    library("clusterProfiler")

    
    
    if(!(all(c("ExpressionView") %in% installed.packages()))) {
        if (Sys.getenv("CONDA_PREFIX")!=""){
            source("https://bioconductor.org/biocLite.R")
            biocLite("ExpressionView")
        } else {
            cat("The 'ExpressionView' package is not installed. You must install it for this script to work!")
        }
    }
    library("ExpressionView")
    
    if(!(all(c("eisa") %in% installed.packages()))) {
        if (Sys.getenv("CONDA_PREFIX")!=""){
            source("https://bioconductor.org/biocLite.R")
            biocLite("eisa",suppressUpdates=TRUE,ask=FALSE,siteRepos=c("http://cran.us.r-project.org"))
        } else {
            cat("The 'eisa' package is not installed. You must install it for this script to work!")
        }
        
    }
    library("eisa")
    




    if(!(all(c("openxlsx") %in% installed.packages()))) {
        if (Sys.getenv("CONDA_PREFIX")!=""){
            install.packages("openxlsx", dependencies=TRUE,repos = "http://cran.us.r-project.org")
        } else {
            cat("The 'openxlsx' package is not installed. You must install it for this script to work!")
        }
        
    }
    library("openxlsx")


    convert_agregate<-function(df,index,subject,sep){
      l1=apply(X = df,MARGIN = 1,FUN = function(x) {
        m=as.data.frame( x = stringi::stri_split(str = x[subject],regex = sep),col.names = c("v1"))
        m["index"]<-x[index]
        return(m[c("index","v1")])
      })
      return(do.call(what = "rbind",args = l1) )
    }

    cluster_enricher<-function(clusters,num,TERM2NAME,TERM2GENE,pAdjustMethod,pvalueCutoff){
      res_red=unique(TERM2GENE[,2])
      Genes=res_red[res_red  %in% unlist(eisa::getFeatureNames(clusters[[num]])) ]
      
      res=enricher(Genes, TERM2GENE=TERM2GENE, 
                   TERM2NAME=TERM2NAME,minGSSize = 0,maxGSSize = length(res_red),pAdjustMethod = pAdjustMethod,pvalueCutoff = pvalueCutoff,qvalueCutoff = pvalueCutoff) 
      return(res)
    }

    generat_urls<-function(allRes,gene2ko){
      temp_table=allRes@compareClusterResult
      temp_table$URL=apply(X = temp_table,MARGIN = 1,FUN = function(x) paste(c("http://www.kegg.jp/kegg-bin/show_pathway?",stringi::stri_replace_all(str = x["ID"],replacement = "",regex = "path:",collapse = ""),"/", paste(sapply( unlist(stringi::stri_split(str = x["geneID"],regex = "/")),FUN = function(x) gene2ko[gene2ko$V1==x,"V2"]),collapse = "+") ),collapse = ""))
      allRes@compareClusterResult=temp_table
      return(allRes)
    }

    clusters_enricher=function(clusters,TERM2NAME,TERM2GENE,file_name,Type,pAdjustMethod,pvalueCutoff,gene2ko=FALSE){
      allRes=list()
      cluster_names=list()
      count=1
      for (i in c(1:length(clusters))){
        temp<-cluster_enricher(clusters,i,TERM2NAME,TERM2GENE,pAdjustMethod,pvalueCutoff)
        if (length(temp)>0){
          allRes[count]<-temp
          cluster_names[count]=i
          count=count+1
          cat("In cluster ",i," : ",length(temp@result$ID)," enriched terms were found\n")
        }
      } 
      names(allRes)=cluster_names
      allRes=clusterProfiler::merge_result(enrichResultList =allRes)
      allRes@fun<-"enrichGO"
      if (Type=="KO"){
        allRes<-generat_urls(allRes,gene2ko)
      }
      write.table(x = allRes@compareClusterResult,file =file_name ,quote = FALSE,row.names = TRUE,sep = "\t")
      if (dim(allRes@compareClusterResult)[1]>0){
        if (Type=="GO"){
          allRes@fun<-"enrichGO"
          allRes=simplify(allRes, cutoff=0.7, by="p.adjust", select_fun=min)
          write.csv(x = allRes@compareClusterResult,file =paste("simplify",file_name,collapse = "_") ,quote = FALSE,row.names = TRUE)
        }
        
        DOSE::dotplot(allRes,x=allRes@compareClusterResult$Cluster,showCategory=1000)+
          ggplot2::ggsave(filename = paste(file_name,".pdf",collapse = ""),dpi = 600,device = "pdf",width = 20,height = 20)
      }
      return(allRes)
    }


    # Frames the bi-cluster values, according to rows and cols
    frameGrid = function(wb, sheet, rows, cols) {
      
      frameSegment = function(wb, sheet, rows, cols, side) {
        addStyle(wb=wb, sheet=sheet, gridExpand=FALSE, stack=TRUE,
                 style=createStyle(border = side, borderStyle = "thick"),
                 rows = rows, cols = cols)
      }
      
      frameSegment(wb=wb, sheet=sheet, rows=min(rows), cols=cols , "Top")
      frameSegment(wb=wb, sheet=sheet, rows=max(rows), cols=cols , "Bottom")
      frameSegment(wb=wb, sheet=sheet, rows=rows , cols=min(cols), "Left")
      frameSegment(wb=wb, sheet=sheet, rows=rows , cols=max(cols), "Right")
    }



    # Writes a sheet per bi-cluster with samples of subset, features of bi-cluster, ordered by bi-cluster.
    writeSheetPerBicluster = function(my.modules, my.eset, outputFile) {
      
      # Create workbook
      wb <- createWorkbook()
      numberStyle = createStyle(numFmt = "0")
      i=0
      new.modules=list()
      # For each bi-cluster
      for (count in 1:length(my.modules)) {
        module = my.modules[[count]]
        
        if ((getNoSamples(module) < min.bicluster.samples)|(getNoFeatures(module) < min.bicluster.genes)) next
        i=i+1
        print(sprintf("Processing bi-cluster %d / %d", i, length(my.modules)))
        
        new.modules[[i]]=module
        # Get module's features and samples
        module.features = unlist(getFeatureNames(module))
        module.samples  = unlist(getSampleNames(module))
        eset.samples    = colnames(my.eset)
        
        # Convert eSet to a matrix 
        my.eset.mat = exprs(my.eset)
        
        # Calculate Order of features by their mean in the bi-cluster (or if 0, in entire subset)
        module.feature.means     = apply(my.eset.mat[module.features ,                       module.samples  , drop=FALSE], 1, mean)
        non.module.feature.means = apply(my.eset.mat[module.features , setdiff(eset.samples, module.samples) , drop=FALSE], 1, mean)
        #means.ratio = pmax(module.feature.means, non.module.feature.means) / (pmin(module.feature.means, non.module.feature.means) + 1)
        #features.order = order(means.ratio, decreasing = TRUE)
        features.order = order(module.feature.means, non.module.feature.means, decreasing = TRUE)
        
        # Order features by their mean in the bi-cluster (or if 0, by the mean of the rest of the subset)
        module.features = module.features[features.order]
        
        # Set and\or sort the relevant features and samples for the biclusters
        features.to.use = c(varLabels(my.eset),
                            module.features)
        samples.to.use  = c(fvarLabels(my.eset),
                            module.samples,
                            setdiff(eset.samples, module.samples))
        
        # Merge matrix with feature metadata
        my.eset.df = merge(my.eset.mat, as.matrix(fData(my.eset)),
                           by="row.names", all.x=TRUE, all.y=FALSE, sort=FALSE)
        # Fix row names
        rownames(my.eset.df) = my.eset.df$Row.names
        my.eset.df$Row.names = NULL
        
        # Merge with samples metadata
        my.eset.df = merge(t(my.eset.df), as.matrix(pData(my.eset)),
                           by="row.names", all.x=TRUE, all.y=FALSE, sort=FALSE)
        # Fix row (transposed columns) names
        rownames(my.eset.df) = my.eset.df$Row.names
        my.eset.df$Row.names = NULL
        # Tranpose back
        my.eset.df = t(my.eset.df)
        
        # Get the subset of the data
        my.eset.df = my.eset.df[features.to.use, samples.to.use]
        
        # Prepare ordered subset to be written
        bicluster.rows       = 1 + length( varLabels(my.eset)) + (1 : length(module.features))
        bicluster.cols       = 1 + length(fvarLabels(my.eset)) + (1 : length(module.samples))
        numerical.data.rows  = bicluster.rows
        numerical.data.cols  = 1 + length(fvarLabels(my.eset)) + (1 : length(eset.samples))
        sample.metadata.rows = function(field) 1 + which(rownames(my.eset.df) == field)
        sample.metadata.cols = numerical.data.cols
        
        # Set sheet's name (e.g. '3_18s' means Bi-cluster #3, 18 samples)
        sheetName = sprintf("%d_%ds", i, getNoSamples(module))
        
        # Add sheet
        addWorksheet(wb, sheetName)
        # Write data + samples metadata + features metadata
        writeData(wb=wb, sheet=sheetName, x = my.eset.df,
                  colNames = TRUE, rowNames = TRUE,  withFilter = FALSE)
        # Overwrite data with numeric values (this is a workaround, disregard the overwrite warning)
        suppressWarnings(
          writeData(wb=wb, sheet=sheetName,
                    x = my.eset.mat[rownames(my.eset.df)[numerical.data.rows - 1],
                                    colnames(my.eset.df)[numerical.data.cols - 1]],
                    startRow = min(numerical.data.rows),
                    startCol = min(numerical.data.cols),
                    colNames = FALSE, rowNames = FALSE,  withFilter = FALSE)
        )
        suppressWarnings(
          writeData(wb=wb, sheet=sheetName,
                    x=paste(featureThreshold(module),
                            sampleThreshold(module) , sep="_"))
        )
        # Style the data
        addStyle(wb=wb, sheet=sheetName,  gridExpand=TRUE, stack=TRUE,
                 style=numberStyle, rows = numerical.data.rows, cols = numerical.data.cols)
        # Frame the bi-cluster
        frameGrid(wb=wb, sheet=sheetName, rows = bicluster.rows, cols = bicluster.cols)
        # Add coloring to sample metadata
        # for (coloring in sample.metadata.coloring)
        #   conditionalFormatting(wb=wb, sheet=sheetName,
        #                         rows = sample.metadata.rows(coloring$field),
        #                         cols = sample.metadata.cols,
        #                         rule = coloring$value, style = createStyle(bgFill = coloring$color),
        #                         type = "contains")
        #       conditionalFormatting(wb=wb, sheet=sheetName,
        #                             rows = sample.metadata.rows(coloring$field),
        #                             cols = sample.metadata.cols,
        #                             rule = paste("='", coloring$value, "'", sep=""), style = createStyle(bgFill = coloring$color),
        #                             type = "expression")
        # Add conditional formatting to data (red colorscale)
        conditionalFormatting(wb=wb, sheet=sheetName,
                              rows = numerical.data.rows, cols = numerical.data.cols,
                              rule = c(0, 1), style = c("white", "red"), type = "colourscale")
      }
      
      # Save workbook
      print("Saving...")
      saveWorkbook(wb=wb, file=outputFile, overwrite = TRUE)
      print("Done.")
      return(new.modules)
    }



    # Read files

    Roary_data =read.csv(Roary_Results_file,
                        header=TRUE , row.names=1, as.is=TRUE, check.names = FALSE)
    metadata =read.delim(metadata_file,
                         header=TRUE , as.is=TRUE, check.names = FALSE)
    
    row.names(metadata)=metadata[,metadata_sample_field]

    if (is.na(metadata_select_fileds)){
      metadata_select_fileds=colnames(metadata)
    }else{
      metadata_select_fileds=unlist(lapply(X =stringi::stri_split(metadata_select_fileds,regex ="," ),FUN = function(x)  stringr::str_trim(gsub('"','',gsub("'",'',x))))  )
    }
    
    metadata_select_fileds=metadata_select_fileds[metadata_select_fileds %in% colnames(metadata)]
    if (length(metadata_select_fileds)>0){
      metadata=subset.data.frame(x = metadata,select=metadata_select_fileds )
    }
    
    Sys.setenv(R_ZIPCMD="zip")

    Annotation_data = as.data.frame(read.delim(Annotation_file,
                                          header=TRUE , row.names=2, as.is=TRUE))

    Roary_matrix=Roary_data[,14:(dim(Roary_data)[2]-1)]
    Roary_matrix=as.matrix(Roary_matrix!="")
    Roary_matrix[Roary_matrix[,]]=1
    Roary_matrix.mean=apply(X =Roary_matrix,MARGIN = c(1),FUN =mean)
    Roary_matrix_filter=Roary_matrix[(Roary_matrix.mean>Roary_filter_min)&(Roary_matrix.mean<Roary_filter_max),]
    
    
    Roary_data$VFG=lapply(X =Roary_data$Annotation,FUN = function(x) unlist(stringr::str_extract_all(string = x,pattern = "VFG[0-9]+")))
    sub_Roary_data=subset.data.frame(x = Roary_data,subset = !S4Vectors::isEmpty(Roary_data$VFG),select = c('VFG','Annotation','No. isolates'))
    sub_Roary_data$'No. isolates'=sub_Roary_data$'No. isolates'/max(sub_Roary_data$'No. isolates')
    sub_Roary_data$vfclass=lapply(X =sub_Roary_data$VFG,FUN = function(x) unique(as.vector(na.exclude( Annotation_data[unlist(x),"vfclass"]))))

    Roary_matrix_filter=Roary_matrix_filter[row.names(Roary_matrix_filter) %in% row.names(sub_Roary_data),] 
    sub_Roary_data_filter=sub_Roary_data[row.names(Roary_matrix_filter),]

    # Slice meta data according to results table
    names=intersect(colnames(Roary_matrix_filter) , rownames(metadata))

    Roary_matrix_filter=Roary_matrix_filter[,names]
    metadata=subset.data.frame(x = metadata,rownames(metadata) %in%  names)
    metadata$LIRON.temp=""
    metadata=metadata[names,]
    metadata$LIRON.temp=NULL

    # Assert that the results table match the meta data in number and order of features and samples
    stopifnot(all(rownames(sub_Roary_data_filter) == rownames(Roary_matrix_filter)))
    stopifnot(all(rownames(metadata)  == colnames(Roary_matrix_filter)))

    # Create the proper annotation structures

    MetaFeature = new("AnnotatedDataFrame", data=sub_Roary_data_filter)
    MetaSample  = new("AnnotatedDataFrame", data=metadata)

    # Create the proper structure for data and meta-data all together
    dbESet = ExpressionSet(assayData   = Roary_matrix_filter,
                           featureData = MetaFeature,
                           phenoData   = MetaSample,
                           annotation  = "hgu95av2")    ## Without this "chip" argument, it doesn't work.


    # Run ISA with random seeds
    set.seed(1)
    dbModules.random = ISA(data = dbESet,
                           flist = NA, uniqueEntrez = FALSE,
                           thr.gene = thr.gene, thr.cond = thr.cond, no.seeds = no.seeds)
    # The next few lines help with adjusting the ISA thresholds
    getNoFeatures(dbModules.random)
    getNoSamples(dbModules.random)

    if (length(dbModules.random)>0){
        # Print randomly-seeded clusters to Excel, one sheet per bi-cluster
        dbModules.random=writeSheetPerBicluster(dbModules.random, dbESet, file.path(output,"Bicluster.xlsx"))


        sub_Roary_data_filter$Gene=row.names(sub_Roary_data_filter)
        sub_Roary_data_filter$VFG2=sapply(X = sub_Roary_data_filter$VFG,FUN = function(x) stringr::str_c(x,collapse = ","))# stringr::str_replace_all(string = as.character(x),pattern = '[c( \" )]',replacement = "") )
        sub_Roary_data_filter$vfclass2=sapply(X = sub_Roary_data_filter$vfclass,FUN = function(x) stringr::str_c(x,collapse = ","))# stringr::str_replace_all(string = as.character(x),pattern = '[c( \" )]',replacement = "") )
          
        mm=convert_agregate(sub_Roary_data_filter,"Gene","vfclass2",",")
        mm=mm[mm$v1!='',]
        Pathway2gene=mm[,c("v1","index")]
        Pathway2name=mm[,c("v1","v1")]
        
        lines=list()
        for (cluster in c(1:length(dbModules.random))){
           lines[cluster]=paste( c( as.character(cluster) ,unlist(eisa::getFeatureNames(dbModules.random[[cluster]])) ) ,collapse = "\t")
        }
        print(unlist(lines))
        writeLines(unlist(lines),file.path(output,"Bicluster_clusters"))

        print("Enrichment analysis... ")
        suppressWarnings(
            invisible(
                      clusters_enricher(dbModules.random,Pathway2name,Pathway2gene,file.path(output,"Bicluster_Enrichment"),"",pAdjustMethod,pvalueCutoff)
            )
        )
        print("Done")
    } else {
    write.table(x =data.frame() ,file =file.path(output,"No_clusters_were_found") ,quote = FALSE,row.names = TRUE,sep = "\t")
    print("No clusters were found!!!")
    
    }
}