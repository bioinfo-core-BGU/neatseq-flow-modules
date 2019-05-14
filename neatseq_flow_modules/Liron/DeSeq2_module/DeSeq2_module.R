library("optparse")
library("ggplot2")
library("pheatmap")
library("mclust")
library("factoextra")
library("cowplot")
library("gridExtra")

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-c", "--COUNT_DATA_FILE"), type="character", default=NA, 
              help="Path to Gene Level Count Data Files [comma sep]", metavar="character"),
  make_option(c("--COUNT_SOURCE"), type="character", default='Matrix', 
              help="The Source of the count Data Files ['Matrix'.'RSEM','HTSEQ']", metavar="character"),
  make_option(c("--SAMPLES"), type="character", default=NA, 
              help="If the COUNT_SOURCE is RSEM/HTSEQ it is possible to link the count files to samples by listing [comma sep] samples names in the same order as the files [sample1,sample2...]", metavar="character"),
  make_option(c("-s", "--SAMPLE_DATA_FILE"), type="character", default=NA, 
              help="Path to Samples Information File", metavar="character"),    
  make_option(c("-o","--outDir"), type="character", default=NA, 
              help="path to the output directory", metavar="character"),
  make_option(c("-t", "--GENE_ID_TYPE"), type="character", default=NA, 
              help="The Gene ID Type i.e 'ENSEMBL'[for Bioconductor] OR 'ensembl_gene_id'/'ensembl_transcript_id' [for ENSEMBL] DON'T USE IF IT IS DE-NOVO ASSEMBLY ", metavar="character"),
  make_option(c("-d", "--Annotation_db"), type="character", default=NA, 
              help="Bioconductor Annotation Data Base Name from https://bioconductor.org/packages/release/BiocViews.html#___OrgDb  ", metavar="character"),
  make_option(c("--Species"), type="character", default=NA, 
              help="Species Name to Retrieve Annotation Data from ENSEMBL", metavar="character"),
  make_option(c("--KEGG_Species"), type="character", default=NA, 
              help="Species Name to Retrieve Annotation Data from KEGG", metavar="character"),
  make_option(c("--KEGG_KAAS"), type="character", default=NA, 
              help="Gene to KO file from KEGG KAAS [first column gene id, second column KO number]", metavar="character"),
  make_option(c("--FILTER_SAMPLES"), action="store_true",default=FALSE, 
              help="Filter Samples with Low Number of expressed genes OR with Small Library size using 'scater' package ", metavar="character"),
  make_option(c("--FILTER_GENES"), action="store_true",default=FALSE, 
              help="Filter Low-Abundance Genes using 'scater' package ", metavar="character"),
  make_option(c( "--Trinotate"), type="character", default=NA, 
              help="Path to a Trinotate annotation file in which the first column is the genes names", metavar="character"),
  make_option(c( "--Rmarkdown"), type="character", default=NA, 
              help="Path to a Rmarkdown File [ usually 'DeSeq2_module.Rmd' ]", metavar="character"),
             
  make_option(c("-n", "--NORMALIZATION_TYPE"), type="character", default='VSD', 
              help="The Normalization Type To Use [VSD , RLOG] The Default is VSD", metavar="character"),
  make_option(c("--BLIND_NORM"), action="store_true",default=FALSE, 
              help="Perform Blind Normalization", metavar="character"),            
  make_option(c("--DESIGN"), type="character", default=NA, 
              help="The Main DeSeq Design", metavar="character"),                      
  make_option(c("--LRT"), type="character", default=NA, 
              help="The LRT DeSeq Design", metavar="character"),
  make_option(c("--ALPHA"), type="numeric", default = 0.05, 
              help="Significant Level Cutoff, The Default is 0.05", metavar="character"),
  make_option(c("--Post_statistical_ALPHA"), type="numeric", default = NA, 
              help="Post Statistical FDR P-value Filtering ", metavar="character"),
  make_option(c("--FoldChange"), type="numeric", default = 1, 
              help="Fold change Cutoff [testing for fold changes greater in absolute value], The Default is 1", metavar="character"),
  make_option(c("--Post_statistical_FoldChange"), type="numeric", default = NA, 
              help="Post Statistical Fold change Filtering ", metavar="character"),
  make_option(c("--CONTRAST"), type="character", default=NA,
              help="The DeSeq Contrast Design [Not For LTR]", metavar="character"),
  make_option(c("--modelMatrixType"), type="character", default="standard",
              help="How the DeSeq model matrix of the GLM formula is formed [standard or expanded] ,The Default is standard", metavar="character"),
  make_option(c("--removeBatchEffect"), type="character", default=NA,
              help="Will Remove Batch Effect from the Normalized counts data up to 2 [using the limma package and only one using the sva package] Batch Effect fields [from the Sample Data ] separated by , ", metavar="character"),
  make_option(c("--removeBatchEffect_method"), type="character", default="sva",
              help="The method to Remove Remove Batch Effect from the Normalized counts data using the limma or sva packages [sva is the default] ", metavar="character"),
  make_option(c("--collapseReplicates"), type="character", default=NA,
              help="Will collapse technical replicates using a Sample Data field indicating which samples are technical replicates", metavar="character"),
  
  
  make_option(c("--GENES_PLOT"), type="character", default=NA,
              help="Genes Id To Plot count Data [separated by ,] ", metavar="character"),
  make_option(c("--X_AXIS"), type="character", default=NA,
              help="The Filed In the Sample Data To Use as X Axis ", metavar="character"),
  make_option(c("--GROUP"), type="character", default=NA,
              help="The Filed In the Sample Data To Group By [can be two fields separated by ,]", metavar="character"),
  make_option(c("--SPLIT_BY"), type="character", default=NA,
              help="The Filed In the Sample Data To Split the Analysis By ", metavar="character"),
  make_option(c("--SPLIT_BY_CONTRAST"), action="store_true",default=FALSE, 
              help="Only use Samples found in the relevant contrast for Clustering and Enrichment Analysis", metavar="character"),
  
  make_option(c("--FUNcluster"), type="character", default='hclust',
              help='A clustering function including [kmeans,pam,clara,fanny,hclust,agnes,diana,click]. The default is hclust', metavar="character"),
  make_option(c("--hc_metric"), type="character", default='euclidean',
              help="Hierarchical clustering metric to be used for calculating dissimilarities between observations. The default is pearson", metavar="character"),
  make_option(c("--hc_method"), type="character", default='ward.D2',
              help="Hierarchical clustering agglomeration method to be used. The default is ward.D2 ", metavar="character"),
  make_option(c("--CLICK_PATH"), type="character", default=NA,
              help="The path to CLICK program (Shamir et al. 2000). If your using click cite: CLICK and Expander ( Shamir et al. 2000 and Ulitsky et al. 2010) ", metavar="character"),
  make_option(c("--CLICK_HOMOGENEITY"), type="numeric", default=0.5,
              help="The HOMOGENEITY [0-1] of clusters using CLICK program (Shamir et al. 2000). The default is 0.5 ", metavar="character"),
  
  
  make_option(c("--k.max"), type="numeric", default=20,
              help="The maximum number of clusters to consider, must be at least two. The default is 20", metavar="character"),
  make_option(c("--nboot"), type="numeric", default=10,
              help="Number of Monte Carlo (bootstrap) samples for determining the number of clusters [Not For Mclust]. The default is 10 ", metavar="character"),
  make_option(c("--stand"), action="store_true",default=FALSE, 
              help="The Data will be Standardized Before Clustering", metavar="character"), 
  make_option(c("--Mclust"), action="store_true",default=FALSE, 
              help="Use Mclust for determining the number of clusters", metavar="character"), 
  
  make_option(c("--Enriched_terms_overlap"), action="store_true",default=FALSE, 
              help="Test for genes overlap in enriched terms", metavar="character"),
  make_option(c("--USE_INPUT_GENES_AS_BACKGROUND"), action="store_true",default=FALSE, 
              help="Use The input Genes as the Background for Enrichment Analysis", metavar="character"),           
  make_option(c("--PCA_COLOR"), type="character", default=NA,
              help="The Filed In the Sample Data To Determine Color In The PCA Plot", metavar="character"),
  make_option(c("--PCA_SHAPE"), type="character", default=NA,
              help="The Filed In the Sample Data To Determine Shape In The PCA Plot", metavar="character"),
  make_option(c("--PCA_SIZE"), type="character", default='Library_sizes',
              help="The Filed In the Sample Data To Determine Size In The PCA Plot", metavar="character"),
              
  make_option(c("--only_clustering"), action="store_true",default=FALSE, 
              help="Don't Perform Differential Analysis!!!", metavar="character"), 
  make_option(c("--significant_genes"), type="character", default=NA,
              help="Use these genes as the set of significant genes [a comma separated list]", metavar="character")
); 


#Get user information
opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", 
                                    option_list=option_list,
                                    epilogue="\n\nAuthor:Liron Levin");
opt = optparse::parse_args(opt_parser);



print("Input var:",,quote = F)
print.data.frame(as.data.frame(x = unlist(opt),row.names = names(opt)),right = F,quote = F)

# Check if required packages are installed:
if(!(all(c('DESeq2') %in% installed.packages()))) {
  if (Sys.getenv("CONDA_PREFIX")!=""){
    source("http://bioconductor.org/biocLite.R")
    biocLite('DESeq2')
    library('DESeq2',character.only =T)
  }else{
    cat("The Bioconductor DESeq2 package is not installed. You must install it for this script to work!")
  }
}else{
  library('DESeq2',character.only =T)  
} 


#####################Functions################################

GetcountDataFromHTSeqCount <- function(sampleTable) # This function was modified from the original DESeq2 DESeqDataSetFromHTSeqCount function 
{
  l <- lapply( as.character( sampleTable[,2] ), function(fn) read.table( file.path( fn ), fill=TRUE ) )
  if( ! all( sapply( l, function(a) all( a$V1 == l[[1]]$V1 ) ) ) )
    stop( "Gene IDs (first column) differ between files." )
  # select last column of 'a', works even if htseq was run with '--additional-attr'
  tbl <- sapply( l, function(a) a[,ncol(a)] )
  colnames(tbl) <- sampleTable[,1]
  rownames(tbl) <- l[[1]]$V1
  rownames(sampleTable) <- sampleTable[,1]
  oldSpecialNames <- c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique")
  # either starts with two underscores
  # or is one of the old special names (htseq-count backward compatability)
  specialRows <- (substr(rownames(tbl),1,1) == "_") | rownames(tbl) %in% oldSpecialNames
  tbl <- tbl[ !specialRows, , drop=FALSE ]
 
  return(tbl)
}   




keggLink_retry <- function (Type,Data,RETRY_NUM=10){
      error_flag=TRUE
      res=c()
      try_count=0
      while ((error_flag)&&(try_count<RETRY_NUM) ){
          try_count = try_count+1
          res<-try(keggLink(Type,Data),silent = T)
          if (inherits(res,"try-error")){
            error_flag=FALSE
          }
      } 
      return(res)
}

plotPCA<-function (object, intgroup = "condition", ntop = 500, returnData = FALSE) { #This function was modified from the original DESeq2 plotPCA function
  if (ntop!='all'){
    rv <- genefilter::rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                       length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
  }else{
    pca <- prcomp(t(assay(object)))
  }
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }else{
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],PC3 = pca$x[, 3], group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:3]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
    geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 
                                                        100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 
                                                                                                            100), "% variance")) + coord_fixed()
}


cluster_enricher<-function(clusters,num,Type,organism,universe,keyType,Ontology="BP",pvalueCutoff=0.05,pAdjustMethod='fdr',qvalueCutoff=0.05){
  res_red=unique(names(clusters))
  Genes=res_red[res_red  %in% names(clusters[clusters %in% num])]
  if (Type=="GO"){
    res=enrichGO(gene      = Genes,
                 OrgDb         = organism,
                 universe      = universe,
                 keytype       = keyType,
                 ont           = Ontology,
                 pAdjustMethod = pAdjustMethod,
                 pvalueCutoff  = pvalueCutoff,
                 qvalueCutoff  = qvalueCutoff)
  }else{
    res=enrichKEGG(gene          = Genes, 
                   organism      = organism,
                   pAdjustMethod = pAdjustMethod,
                   pvalueCutoff  = pvalueCutoff,
                   qvalueCutoff  = qvalueCutoff,
                   universe      = universe,
                   keyType       = keyType,  
                   minGSSize     = 0,
                   maxGSSize     = length(res_red))
    
    
    
  }
  return(res)
}


clusters_enricher=function(outDir, clusters,Type,organism,universe,keyType,file_name,Ontology="BP",pvalueCutoff=0.05,pAdjustMethod='fdr',qvalueCutoff=0.05){
  allRes=list()
  cluster_names=list()
  count=1
  for (i in sort(unique(clusters))){
    temp=cluster_enricher(clusters,c(i),Type,organism,universe,keyType,Ontology,pvalueCutoff,pAdjustMethod,qvalueCutoff)
    if (length(temp)>0){
      allRes[count]<-temp
      cluster_names[count]=i
      count=count+1
    }
  } 
  names(allRes)=cluster_names
  allRes=clusterProfiler::merge_result(enrichResultList =allRes)
  allRes@fun<-"enrichGO"
  write.csv(x = allRes@compareClusterResult,
            file =file.path(outDir,file_name) ,
            quote = TRUE,
            row.names = TRUE)
  if (dim(allRes@compareClusterResult)[1]>0){
    if (Type=="GO"){
        allRes@fun<-"enrichGO"
        temp_allRes<-try(clusterProfiler::simplify(allRes, cutoff=0.7, by="p.adjust", select_fun=min),silent = T)
        if (!inherits(temp_allRes,"try-error")){
            allRes=temp_allRes
            file_name=paste("simplify",file_name,collapse = "_")
            write.csv(x = allRes@compareClusterResult,file =file.path(outDir,file_name) ,quote = FALSE,row.names = TRUE)
        }
    }else{
      allRes@fun<-"enrichKEGG"
    }
    if (dim(allRes@compareClusterResult)[1]>50){
        font.size=4
    }else{
        font.size=9
    }
    
    if (dim(allRes@compareClusterResult)[1]>0){
        DOSE::dotplot(allRes,x=allRes@compareClusterResult$Cluster,font.size =font.size,showCategory=1000)+
        ggplot2::ggsave(filename = file.path(outDir,paste(file_name,".pdf",collapse = "")),dpi = 600,device = "pdf",width = 20,height = 20)
    }
  }
  return(allRes)
}


General_Enrichment_Test<-function(clusters,num,TERM2NAME,TERM2GENE,pAdjustMethod='fdr',pvalueCutoff=0.05){
  res_red=unique(TERM2GENE[,2])
  Genes=res_red[res_red  %in% names(clusters[clusters %in% num])]
  
  res=enricher(Genes, TERM2GENE=TERM2GENE, 
               TERM2NAME=TERM2NAME,
               minGSSize     = 0,
               maxGSSize     = length(res_red),
               pAdjustMethod = pAdjustMethod,
               pvalueCutoff  = pvalueCutoff
  ) 
  
  return(res)
}


generat_urls<-function(allRes,gene2ko){
  temp_table=allRes@compareClusterResult
  temp_table$URL=apply(X = temp_table,MARGIN = 1,FUN = function(x) paste(c("http://www.kegg.jp/kegg-bin/show_pathway?",stringi::stri_replace_all(str = x["ID"],replacement = "",regex = "path:",collapse = ""),"/", paste(sapply( unlist(stringi::stri_split(str = x["geneID"],regex = "/")),FUN = function(x) gene2ko[gene2ko$V1==x,"V2"]),collapse = "+") ),collapse = ""))
  allRes@compareClusterResult=temp_table
  return(allRes)
}


Clusters_Enrichment_Test=function(outDir,clusters,TERM2NAME,TERM2GENE,file_name,Type,pAdjustMethod='fdr',pvalueCutoff=0.05,gene2ko=FALSE){
  allRes=list()
  cluster_names=list()
  count=1
  for (i in sort(unique(clusters))){
    temp=General_Enrichment_Test(clusters,c(i),TERM2NAME,TERM2GENE,pAdjustMethod,pvalueCutoff )
    if (length(temp)>0){
      allRes[count]<-temp
      cluster_names[count]=i
      count=count+1
    }
  } 
  names(allRes)=cluster_names
  allRes=clusterProfiler::merge_result(enrichResultList =allRes)
  if (Type=="KO"){
    allRes<-generat_urls(allRes,gene2ko)
  }
  write.csv(x = allRes@compareClusterResult,
            file =file.path(outDir,file_name) ,
            quote = TRUE,
            row.names = TRUE)
  if (Type=="GO"){
    allRes@fun<-"enrichGO"
    temp_allRes<-try(clusterProfiler::simplify(allRes, cutoff=0.7, by="p.adjust", select_fun=min),silent = T)
    if (!inherits(temp_allRes,"try-error")){
        allRes=temp_allRes
        file_name=paste("simplify",file_name,collapse = "_")
        write.csv(x = allRes@compareClusterResult,file =file.path(outDir,file_name) ,quote = FALSE,row.names = TRUE)
    }
  }else{
    allRes@fun<-"enrichKEGG"
  }
    if (dim(allRes@compareClusterResult)[1]>50){
        font.size=4
    }else{
        font.size=9
    }
  
  if (dim(allRes@compareClusterResult)[1]>0){
    DOSE::dotplot(allRes,showCategory=1000,font.size=font.size)+ 
    ggplot2::ggsave(filename = file.path(outDir,paste(file_name,".pdf",collapse = "")),dpi = 600,device = "pdf",width = 20,height = 20)
  }
  return(allRes)
}


plot_clusters<-function(clusters,vs_ano,color_group=c('Type','Time_int'),titles=c("Type","Time","EXP"),split_by=c(),smooth=T,X_AXIS_ORDER=NA){
    plot_list=list()
    count=1
    for (i in sort(unique(clusters))){
        Df=data.frame()
        genes=names(clusters[clusters==i])
        for (j in genes){
          if (length(split_by)>0){
            temp=subset.data.frame(x =vs_ano,,c(split_by,color_group,j))
            colnames(temp)=c("split_by","Type","Time","EXP")
            Df=rbind(Df,temp)
          }else{
            temp=subset.data.frame(x =vs_ano,,c(color_group,j))
            colnames(temp)=c("Type","Time","EXP")
            Df=rbind(Df,temp)
          }
        }
        if (color_group[1]==color_group[2]){
          Df$Type='Trend'
        }
        if (length(X_AXIS_ORDER)>1){
            Df$Time <-factor(Df$Time,levels=intersect(X_AXIS_ORDER,unique(Df$Time)))
        }
        if (length(split_by)>0){
            plot_list[count]=list(ggplot(data=Df, aes(x=Time, y=EXP, group=Type)) + 
                                    facet_wrap(~split_by, ncol=1,scales = "free_y",strip.position = c("right"))+
                                    theme(strip.text.y = element_text(size = 7, colour = "black"))+#,face="bold" ,angle = 90))+
                                    #theme(strip.background = element_blank())+
                                    ggtitle( paste(paste("Cluster ", i  ,sep="") ,length(genes) ,sep="\n")) +
                                    #theme(plot.title = element_text(colour = "black", size = 12)) + 
                                    theme(legend.position  ="none")+
                                    xlab(titles[2]) +
                                    ylab(titles[3]) +
                                    theme(axis.title =element_text(colour = "black", size = 10) )+
                                    theme(legend.title =element_text(colour = "black", size = 10) )+
                                    theme(legend.text =element_text(colour = "black", size = 8) )+
                                    theme(axis.text.y =   element_text(colour = "black", size = 6))+
                                    theme(axis.text.x =   element_text(colour = "black", size = 6))+
                                    theme(plot.margin = unit(c(0.5,0.2,0.5,0.2),"cm")  )+
                                    #scale_color_discrete(name=titles[1])+
                                    #scale_linetype_manual(name=titles[1],values=c(1,5))+
                                    scale_x_discrete(limits=unique(Df$Time))+
                                    theme(legend.key.width =unit(3,"line")) 
                                    #scale_y_continuous(breaks=c(seq(0,10,by=2)) )
                                )
        }else{
            plot_list[count]=list(ggplot(data=Df, aes(x=Time, y=EXP, group=Type)) + 
                                    ggtitle( paste(paste("Cluster ", i  ,sep="") ,length(genes) ,sep="\n")) +
                                    #theme(plot.title = element_text(colour = "black", size = 12)) + 
                                    theme(legend.position  ="none")+
                                    xlab(titles[2]) +
                                    ylab(titles[3]) +
                                    theme(axis.title =element_text(colour = "black", size = 10) )+
                                    theme(legend.title =element_text(colour = "black", size = 10) )+
                                    theme(legend.text =element_text(colour = "black", size = 8) )+
                                    theme(axis.text.y =   element_text(colour = "black", size = 6))+
                                    theme(axis.text.x =   element_text(colour = "black", size = 6))+
                                    theme(plot.margin = unit(c(0.5,0.2,0.5,0.2),"cm")  )+
                                    #scale_color_discrete(name=titles[1])+
                                    #scale_linetype_manual(name=titles[1],values=c(1,5))+
                                    scale_x_discrete(limits=unique(Df$Time))+
                                    theme(legend.key.width = unit(3,"line")) 
                                    #scale_y_continuous(breaks=c(seq(0,10,by=2)) )
                                )
        }

       
        if (smooth){
            plot_list[count]=list(plot_list[[count]]+
                                    stat_smooth(method="loess", fullrange=F, size=0.5 ,aes(linetype=Type,color=Type)) 
                                )
            gp2=plot_list[count]
            res<-try(get_legend(gp2[[1]] + theme(legend.position  ="right")+guides(color=guide_legend(title=titles[1]),linetype=guide_legend(title=titles[1]))) ,silent = T)
            if (inherits(res,"try-error")){
                plot_list[count]=list(plot_list[[count]]+
                                        stat_summary(fun.y=mean, geom="point", size = 3,aes(color=Type))+
                                        stat_summary(fun.data = "mean_se", geom = "errorbar", width = .3,size = 1,aes(color=Type),show.legend = F)
                                     )
            }
        }else{
            plot_list[count]=list(plot_list[[count]]+
                                    stat_summary(fun.y=mean, geom="line", size = 1.3,aes(linetype=Type,color=Type))+
                                    stat_summary(fun.data = "mean_se", geom = "errorbar", width = .3,size = 1,aes(color=Type),show.legend = F)
                                )
        }

        count=count+1  
    }
    
    gp2=plot_list[1]
    if (color_group[1]==color_group[2]){
        ml<-marrangeGrob(plot_list, nrow=2, ncol=3, top=NULL)
    }else{
        legend <- get_legend(gp2[[1]] + theme(legend.position  ="right")+guides(color=guide_legend(title=titles[1]),linetype=guide_legend(title=titles[1])))
        ml<-marrangeGrob(plot_list, nrow=2, ncol=3,right=legend, top=NULL)
        
    }
    return(ml)
}


convert_agregate<-function(df,index,subject,sep){
  l1=apply(X = df,MARGIN = 1,FUN = function(x) {
    m=as.data.frame( x = stringi::stri_split(str = x[subject],regex = sep),col.names = c("v1"))
    m["index"]<-x[index]
    return(m[c("index","v1")])
  })
  return(do.call(what = "rbind",args = l1) )
}


plot_sheard_genes<-function(allRes,outDir,file_name){
    if (length(allRes)>1){
        for (cluster in unique(allRes$Cluster)){
            TEMP = na.omit(allRes[allRes$Cluster==cluster,])
            genes = unique(unlist(stringi::stri_split(paste(TEMP$geneID,sep = ,collapse = '/'),fixed = '/')))
            
            if ((length(TEMP$ID)>1) &&((length(genes)>1))){
                mat  = matrix(0, nrow = length(TEMP$ID), ncol = length(genes))
                rownames(mat) = TEMP$ID
                colnames(mat) = genes


                for (x in TEMP$ID){
                  for (y in genes){
                    mat[x,y] = length( unlist(intersect( unlist(stringi::stri_split(TEMP[TEMP$ID==x,'geneID'],fixed = '/')),y )))
                    
                  }
                }

                rownames(mat) = TEMP$Description
                
                if ((length(TEMP$ID)>20) || ((length(genes)>200))){
                    heat_map = pheatmap(mat = mat,cluster_rows = T,
                             cluster_cols = T,
                             silent = T,
                             legend = F)
                
                }else{
                    if (length(unique(as.vector(mat)))==2 ){
                        color=colorRampPalette(c('white','pink'))(2)
                        mybreaks=NA
                    }else{
                        color=colorRampPalette(c('pink'))(1)
                        mybreaks=c(0,1)
                    }
                    heat_map = pheatmap(mat = mat,cluster_rows = T,
                            treeheight_col = 0,
                            treeheight_row = 0,
                            width = 5,
                            height =5,
                            cellheight = 10,
                            fontsize_row =5,
                            fontsize_col = 1,
                            main = paste('Cluster ',cluster),
                            border_color = 'white',
                            cluster_cols = T,
                            breaks = mybreaks,
                            filename = file.path(outDir,paste('Cluster',cluster,'_genes2term_',file_name,".pdf",collapse = "")),
                            color = color,
                            legend = F)
                         
                         }
                write.csv(x = mat[heat_map$tree_row$order,heat_map$tree_col$order],
                      quote = T,
                      row.names = T,
                      file = file.path(outDir,paste('Cluster',cluster,'_genes2term_',file_name,".csv",collapse = "")))
        
            
                mat  = matrix(0, nrow = length(TEMP$ID), ncol = length(TEMP$ID))
                rownames(mat) = TEMP$ID
                colnames(mat) = TEMP$ID

                for (x in TEMP$ID){
                  for (y in TEMP$ID){
                    mat[x,y] = (2*length( unlist(intersect( unlist(stringi::stri_split(TEMP[TEMP$ID==x,'geneID'],fixed = '/')),unlist(stringi::stri_split(TEMP[TEMP$ID==y,'geneID'],fixed = '/'))))))/(length(unlist(stringi::stri_split(TEMP[TEMP$ID==x,'geneID'],fixed = '/'))) +  length(unlist(stringi::stri_split(TEMP[TEMP$ID==y,'geneID'],fixed = '/'))))
                    
                  }
                }

                rownames(mat) = TEMP$Description
                colnames(mat) = TEMP$Description
                
                if (length(TEMP$ID)>20){
                    heat_map = pheatmap(mat = mat,cluster_rows = T,
                                 cluster_cols = T,
                                 silent = T,
                                 legend = T)
                
                }else{
                    if (length(unique(as.vector(mat)))>1 ){
                        color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu") ))( 100)
                        mybreaks=seq(0,1,0.01)
                    }else{
                        color=colorRampPalette(c('red'))(1)
                        mybreaks=c(0,1)
                    }
                
                    heat_map = pheatmap(mat = mat,cluster_rows = T,
                                width = 5,
                                height =5,
                                treeheight_col = 0,
                                treeheight_row = 0,
                                breaks = mybreaks,
                                fontsize_row =5,
                                fontsize_col = 5,
                                main = paste('Cluster ',cluster),
                                border_color = 'white',
                                cluster_cols = T,
                                color = color,
                                filename = file.path(outDir,paste('Cluster',cluster,'_term2term_',file_name,".pdf",collapse = "")),
                                legend = T)
                }
                write.csv(x = mat[heat_map$tree_row$order,heat_map$tree_col$order],
                          quote = T,
                          row.names = T,
                          file = file.path(outDir,paste('Cluster',cluster,'_term2term_',file_name,".csv",collapse = "")))
            
                
            }
        }
        
    }
}


Run_click<-function(mat,click_path,outDir,HOMOGENEITY = 0.5){
      writeLines(text =  paste(dim(mat)[1],dim(mat)[2],sep = ' ' ),con = file.path(outDir,'clickInput.orig'))
      write.table(x = mat,append = T,file = file.path(outDir,'clickInput.orig'),sep = '\t',col.names = F)
      lines=c()
      lines = c(lines,'DATA_TYPE')
      lines = c(lines,'FP ')
      lines = c(lines,'INPUT_FILES_PREFIX')
      lines = c(lines,file.path(outDir,'clickInput '))
      lines = c(lines,'OUTPUT_FILE_PREFIX')
      lines = c(lines,file.path(outDir,'clickOutput '))
      lines = c(lines,'SIMILARITY_TYPE')
      lines = c(lines,'CORRELATION ')
      lines = c(lines,'HOMOGENEITY')
      lines = c(lines,as.character(HOMOGENEITY))
      writeLines(text =  lines,con = file.path(outDir,'params.txt'))
      system(command = paste(click_path,file.path(outDir,'params.txt'),sep = ' ') )
      clusters = read.table(file = file.path(outDir,'clickOutput.res.sol'),sep = '\t',header = F,row.names = 1)
      clusters = setNames(unlist(as.list(clusters)), row.names(clusters))
      return(clusters)
}


#####################Functions-End################################

#####################Read count data and pars it##################

if (opt$COUNT_SOURCE=='RSEM'){
  print('Reading count data..')
  # Check if required packages are installed:
  if(!(all(c('tximport') %in% installed.packages()))) {
    if (Sys.getenv("CONDA_PREFIX")!=""){
      source("https://bioconductor.org/biocLite.R")
      biocLite('tximport')
      library('tximport',character.only =T)
    }else{
      cat("The Bioconductor tximport package is not installed. You must install it for this script to work!")
    }
  }else{
    library('tximport',character.only =T)  
  } 
  
  
  
  files   = unlist(stringi::stri_split(str = opt$COUNT_DATA_FILE,fixed = ','))
  if (!is.na(opt$SAMPLES)){
    samples = unlist(stringi::stri_split(str = opt$SAMPLES,fixed = ','))
    if (length(files)==length(samples)){
      names(files) = samples
    }else{
        stop( "Number of Sample's Names is not Equal the Number of Sample's Files" )
    }
  }else{
    stop( "Sample's names mast be indicated! [--SAMPLES] and in the same order as in the sample's files location [--COUNT_DATA_FILE]" )
  }
    
  txi.rsem <- tximport(files, type = "rsem")
  if (!is.na(opt$GENE_ID_TYPE)) {
      rownames(txi.rsem$abundance)=lapply(X =rownames(txi.rsem$abundance),FUN = function(x) unlist(stringi::stri_split(str = x,regex = "_"))[1] )
      rownames(txi.rsem$abundance)=lapply(X =rownames(txi.rsem$abundance),FUN = function(x) rev(unlist(stringi::stri_split(str = x,regex = ":")))[1] )
      
      rownames(txi.rsem$counts)=lapply(X =rownames(txi.rsem$counts),FUN = function(x) unlist(stringi::stri_split(str = x,regex = "_"))[1] )
      rownames(txi.rsem$counts)=lapply(X =rownames(txi.rsem$counts),FUN = function(x) rev(unlist(stringi::stri_split(str = x,regex = ":")))[1] )
      
      rownames(txi.rsem$length)=lapply(X =rownames(txi.rsem$length),FUN = function(x) unlist(stringi::stri_split(str = x,regex = "_"))[1] )
      rownames(txi.rsem$length)=lapply(X =rownames(txi.rsem$length),FUN = function(x) rev(unlist(stringi::stri_split(str = x,regex = ":")))[1] )
  }
  
  colnames(txi.rsem$abundance) = make.names(colnames(txi.rsem$abundance))
  colnames(txi.rsem$counts) = make.names(colnames(txi.rsem$counts))
  colnames(txi.rsem$length) = make.names(colnames(txi.rsem$length))
  
  countData = txi.rsem$counts
  
  print('Done reading count data')
}else if (opt$COUNT_SOURCE=='HTSEQ'){
    print('Reading count data..')
    files   = unlist(stringi::stri_split(str = opt$COUNT_DATA_FILE,fixed = ','))

    if (!is.na(opt$SAMPLES)){
        samples = unlist(stringi::stri_split(str = opt$SAMPLES,fixed = ','))
        if (length(files)!=length(samples)){
            stop( "Number of Sample's Names is not Equal the Number of Sample's Files" )
        }
    }else{
        stop( "Sample's names mast be indicated! [--SAMPLES] and in the same order as in the sample's files location [--COUNT_DATA_FILE]" )
    }
    sampleTable <- data.frame(sampleName = samples,
                              fileName   = files )
    countData = GetcountDataFromHTSeqCount(sampleTable)
    
    colnames(countData) = make.names(colnames(countData))
    
    print('Done reading count data')
    
    
}else{
  print('Reading count data..')
  countData <- as.matrix(read.csv(opt$COUNT_DATA_FILE,sep="\t",row.names=1))
  if (opt$only_clustering==FALSE){
    countData = apply(X =countData,MARGIN = c(1,2),FUN =round )
  }
  colnames(countData)=lapply(X =colnames(countData),FUN = function(x) unlist(stringi::stri_split(str = x,regex = ".genes.results"))[1] )
  if (!is.na(opt$GENE_ID_TYPE)) {
    rownames(countData)=lapply(X =rownames(countData),FUN = function(x) unlist(stringi::stri_split(str = x,regex = "_"))[1] )
    rownames(countData)=lapply(X =rownames(countData),FUN = function(x) rev(unlist(stringi::stri_split(str = x,regex = ":")))[1] )
  }
  colnames(countData) = make.names(colnames(countData))
  print('Done reading count data')
  
}


#Get Annotation
Annotation = NA
Annotation_flag = ''
if (!is.na(opt$Annotation_db) ){
  
  # Check if required packages are installed:
  if(!(all(c(opt$Annotation_db) %in% installed.packages()))) {
    if (Sys.getenv("CONDA_PREFIX")!=""){
      source("https://bioconductor.org/biocLite.R")
      biocLite(opt$Annotation_db)
      library(opt$Annotation_db,character.only =T)
    }else{
      cat("The Bioconductor Annotation Data Base package is not installed. You must install it for this script to work!")
    }
  }else{
    library(opt$Annotation_db,character.only =T)  
  } 
  
  if (paste("package:",opt$Annotation_db,sep='') %in% search()){
    Annotation=bitr(rownames(countData), fromType=opt$GENE_ID_TYPE, toType=c("ENTREZID",'UNIPROT',"GO","ONTOLOGY","GENENAME","PATH","SYMBOL","PFAM"), OrgDb=opt$Annotation_db)
    Annotation$kegg= sapply(X = Annotation$PATH,FUN = function(X) paste('path:map',x,sep='' ) )
    Annotation_flag = 'OrgDb'
  }
}

if (is.na(Annotation)) {
  if (!is.na(opt$Species)){
    # Check if required packages are installed:
    if(!(all(c('biomaRt') %in% installed.packages()))) {
      if (Sys.getenv("CONDA_PREFIX")!=""){
        source("https://bioconductor.org/biocLite.R")
        biocLite('biomaRt')
        library('biomaRt',character.only =T)
      }else{
        cat("The Bioconductor biomaRt Annotation package is not installed. You must install it for this script to work!")
      }
    }else{
      library('biomaRt',character.only =T)  
    } 
    
    if ("package:biomaRt" %in% search()){
      dataset=c()
      try_count=0
      while ((length(dataset) == 0)&&(try_count<10) ){
          try_count = try_count+1
          host = "www.ensembl.org"
          Test_host<-try(biomaRt::listMarts(host=host),silent = T)
          if (inherits(Test_host,"try-error")){
            host = "ensembl.org"
          }
          biomart = "ENSEMBL_MART_ENSEMBL"
          M=biomaRt::useMart(biomart = biomart, host = host)
          ensembl=biomaRt::listDatasets(M)
          dataset= ensembl$dataset[(sapply(X = ensembl$description,FUN = function(X) stringi::stri_startswith(fixed =  stringi::stri_trans_tolower(opt$Species),str =  stringi::stri_trans_tolower(X)    )))]
          if (length(dataset) == 0){
            list_dataset=unlist(stringi::stri_split(str = unlist(stringi::stri_split(str =  unlist( t(ensembl) ),regex = '\n \nTableSet\t')),regex = '\t'))
            dataset<-try( list_dataset[which(stringi::stri_startswith(str = stringi::stri_trans_tolower(list_dataset),fixed = stringi::stri_trans_tolower(opt$Species)))-1],silent = T) 
            if (inherits(Test_host,"try-error")){
                dataset=''
            }
          }
          if (length(dataset) == 0){
            host = "www.plants.ensembl.org"
            Test_host<-try(biomaRt::listMarts(host=host),silent = T)
            if (inherits(Test_host,"try-error")){
              host = "plants.ensembl.org"
            }
            biomart = 'plants_mart'
            M=biomaRt::useMart(biomart = biomart,host = host)
            plant=biomaRt::listDatasets(M)
            dataset= plant$dataset[(sapply(X = plant$description,FUN = function(X) stringi::stri_startswith(fixed =  stringi::stri_trans_tolower(opt$Species),str =  stringi::stri_trans_tolower(X)    )))]
          }
          if (length(dataset) == 0){
            list_dataset=unlist(stringi::stri_split(str = unlist(stringi::stri_split(str =  unlist( t(plant) ),regex = '\n \nTableSet\t')),regex = '\t'))
            dataset<-try( list_dataset[which(stringi::stri_startswith(str = stringi::stri_trans_tolower(list_dataset),fixed = stringi::stri_trans_tolower(opt$Species)))-1],silent = T) 
            if (inherits(Test_host,"try-error")){
                dataset=c()
            }
          }
          try_count2 =0
          ensembl <- try(biomaRt::useEnsembl(biomart=biomart, dataset=as.character(dataset),host = host),silent = T)
          while ((inherits(ensembl,"try-error")) &&( try_count2<3) ){
            try_count2 = try_count2+1
            ensembl <- try(biomaRt::useEnsembl(biomart=biomart, dataset=as.character(dataset),host = host),silent = T)
          }
          if (inherits(ensembl,"try-error")){
            dataset=c()
          }
        }
        
      if (length(dataset)==1){
        print('Found Species:')
        print(dataset)
        if ( opt$GENE_ID_TYPE %in% listAttributes(ensembl)$name){
          Total_Annotation = c()
          if (opt$USE_INPUT_GENES_AS_BACKGROUND){
            genes <- rownames(countData)
          } else {
            genes <- getBM(attributes = opt$GENE_ID_TYPE, mart = ensembl)
          }
          attributes=c( 'description',
                        'external_gene_name',
                        'kegg_enzyme',
                        'go_id',
                        'namespace_1003')
          if (!opt$GENE_ID_TYPE %in% attributes){
            attributes = c(opt$GENE_ID_TYPE,attributes)
          }
          for (subset_genes in split(genes[,opt$GENE_ID_TYPE], ceiling(seq_along(genes[,opt$GENE_ID_TYPE] )/100))){
              Annotation=getBM(values = subset_genes,
                               filters = opt$GENE_ID_TYPE,
                               attributes = attributes,
                               uniqueRows = TRUE,
                               mart = ensembl)
            Total_Annotation = rbind(Total_Annotation,Annotation)
            }
          Annotation = Total_Annotation
          Annotation$kegg = sapply(X = Annotation$kegg_enzyme,FUN = function(X) unlist(stringi::stri_split(str = X,fixed  = '+'))[1] )
          Annotation$kegg = sapply(X = Annotation$kegg,FUN = function(X) if (!is.na(X)) if (X!='')  paste('path:map',X,collapse = "",sep='' ) else NA else NA  )
          Annotation$ontology = Annotation$namespace_1003
          Annotation$namespace_1003 = NULL
          Kegg_Annotation = getBM(attributes=c(opt$GENE_ID_TYPE,
                                               'uniprotsptrembl'),
                                  mart = ensembl)
                                  
          colnames(Kegg_Annotation) = c(opt$GENE_ID_TYPE,"UNIPROT")
          Annotation_flag = 'ensembl'
        }else{
          
          cat("Your GENE ID TYPE is not recognized ")
        }
        
        
      }else{
        cat("Unable to find the Species Annotation in Ensembl, choose from:")
        print(ensembl$description)
        print(plant$description)
      }
      
    }
  }
  
}


# Get KEGG information:
KEGG_flag          = FALSE
ORGANISM_KEGG_flag = FALSE
KEGG_KAAS_flag     = FALSE
GO_flag            = FALSE
if(!(all(c('clusterProfiler') %in% installed.packages()))) {
  if (Sys.getenv("CONDA_PREFIX")!=""){
    source("https://bioconductor.org/biocLite.R")
    biocLite('clusterProfiler')
  }else{
    cat("The Bioconductor clusterProfiler package is not installed. You must install it for this script to work!")
  }
} 
library("clusterProfiler")

if ("package:clusterProfiler" %in% search()){
  
  if(!(all(c('KEGGREST') %in% installed.packages()))) {
    if (Sys.getenv("CONDA_PREFIX")!=""){
      source("https://bioconductor.org/biocLite.R")
      biocLite('KEGGREST')
    }else{
      cat("The Bioconductor KEGGREST package is not installed. You must install it for this script to work!")
    }
  } 
  library("KEGGREST") 
  if ("package:KEGGREST" %in% search()){
    if (!is.na(Annotation)) {      
      if (Annotation_flag == 'OrgDb'){
        x=as.data.frame(eval(parse(text=sub(".db$", "_dbInfo()", opt$Annotation_db))))
        ORGANISM=x[x["name"]=="ORGANISM","value"]
        ORGANISM=clusterProfiler::search_kegg_organism(ORGANISM)$kegg_code
        if (length(ORGANISM)>0){
          Kegg_Annotation=bitr(rownames(countData), fromType=opt$GENE_ID_TYPE, toType=c('UNIPROT'), OrgDb=opt$Annotation_db)
        }
        
      }else{
        ORGANISMs = clusterProfiler::search_kegg_organism('')
        if (!is.na(opt$KEGG_Species)){
            opt$Species = opt$KEGG_Species
        }
        ORGANISM  = ORGANISMs$kegg_code[(sapply(X = ORGANISMs$scientific_name,FUN = function(X) stringi::stri_startswith(fixed =  stringi::stri_replace_last(str = stringi::stri_trans_tolower(opt$Species), fixed = ' genes',replacement = ''),str =  stringi::stri_trans_tolower(X)    )))]
        if (length(ORGANISM)==0){
          ORGANISMs = na.omit(ORGANISMs)
          ORGANISM  = ORGANISMs$kegg_code[(sapply(X = ORGANISMs$common_name,FUN = function(X) stringi::stri_cmp_eq( stringi::stri_replace_last(str = stringi::stri_trans_tolower(opt$Species), fixed = ' genes',replacement = '') , stringi::stri_trans_tolower(X)    )))]
          if (length(ORGANISM)==0){
            ORGANISM  = ORGANISMs$kegg_code[(sapply(X = ORGANISMs$common_name,FUN = function(X) stringi::stri_startswith(fixed = stringi::stri_replace_last(str = stringi::stri_trans_tolower(opt$Species), fixed = ' genes',replacement = '') , str = stringi::stri_trans_tolower(X)    )))]
          }
        }
      }
      
      
      
      if (length(ORGANISM)==1){       
        print(ORGANISM)
        kegg=keggConv(ORGANISM,"uniprot") #'ncbi-geneid')
        uniprot =  unlist(lapply(X =names(kegg) ,FUN = function(x) rev(unlist(stringi::stri_split(str = x,regex = ":")))[1]))
        kegg=as.data.frame(kegg)
        kegg["uniprot"]=uniprot
        
        KO=keggLink(ORGANISM,"ko")
        ko=names(KO)
        KO=as.data.frame(KO )
        colnames(KO)="kegg"
        KO["ko"]=ko
        
        temp=merge.data.frame(kegg,KO,by = "kegg",all = T)
        
        Pathway=keggLink(ORGANISM,"pathway")
        pathway=names(Pathway)
        Pathway=as.data.frame(Pathway)
        colnames(Pathway)="kegg"
        Pathway["pathway"]=pathway
        
        temp=merge.data.frame(temp,Pathway,by = "kegg",all = T)
        Pathway_info=keggList("pathway",organism =ORGANISM )
        info=names(Pathway_info)
        Pathway_info=as.data.frame(Pathway_info)
        colnames(Pathway_info)="info"
        Pathway_info["pathway"]=info
        temp=merge.data.frame(temp,Pathway_info,by = "pathway",all.x = T)
        Kegg_Annotation=merge.data.frame(Kegg_Annotation,temp,by.x ="UNIPROT" ,by.y = "uniprot",all.x = T)
        
        Annotation$kegg = NULL 
        Annotation = merge.data.frame(Annotation,Kegg_Annotation,by.x =opt$GENE_ID_TYPE ,by.y = opt$GENE_ID_TYPE,all.x = T)
        Annotation$pathway_info = Annotation$info
        Annotation$info = NULL
        
        
        
        aggregate_by=Kegg_Annotation[opt$GENE_ID_TYPE]
        Kegg_Annotation[opt$GENE_ID_TYPE]=NULL
        Kegg_Annotation=aggregate.data.frame(x =Kegg_Annotation ,by = aggregate_by,FUN =function(x) paste(unique(na.omit(x)),sep = "|",collapse = "|"))
        
        ORGANISM_Pathway2gene=convert_agregate(Kegg_Annotation,opt$GENE_ID_TYPE,"pathway","[\\|]")
        ORGANISM_Pathway2gene=ORGANISM_Pathway2gene[ORGANISM_Pathway2gene$v1!='',]
        ORGANISM_Pathway2gene=ORGANISM_Pathway2gene[rev(colnames(ORGANISM_Pathway2gene))]
        ORGANISM_Pathway2name=Pathway_info[rev(colnames(Pathway_info))]
        ORGANISM_KEGG_flag  = TRUE
        
        write.table(ORGANISM_Pathway2name,
                file = file.path(opt$outDir,paste('ORGANISM_Pathway2name','tab', sep = '.')) ,
                quote = F,
                row.names = F,
                sep = "\t")
                
        write.table(ORGANISM_Pathway2gene,
                file = file.path(opt$outDir,paste('ORGANISM_Pathway2gene','tab', sep = '.')) ,
                quote = F,
                row.names = F,
                sep = "\t")
        
      }else{
        if ("kegg" %in% colnames(Annotation)){
          
          Pathway_info=keggList("pathway")
          info=names(Pathway_info)
          Pathway_info=as.data.frame(Pathway_info)
          colnames(Pathway_info)="info"
          Pathway_info["pathway"]=info
          
          
          Annotation = merge.data.frame(Annotation,Pathway_info,by.x ="kegg" ,by.y = "pathway",all.x = T)
          Annotation$pathway = Annotation$kegg
          Annotation$kegg = NULL
          Annotation$pathway_info = Annotation$info
          Annotation$info = NULL
          General_Kegg_Annotation = Annotation[c(opt$GENE_ID_TYPE,'pathway')]
          
          
          
          aggregate_by=General_Kegg_Annotation[opt$GENE_ID_TYPE]
          General_Kegg_Annotation[opt$GENE_ID_TYPE]=NULL
          General_Kegg_Annotation=aggregate.data.frame(x =General_Kegg_Annotation ,by = aggregate_by,FUN =function(x) paste(unique(na.omit(x)),sep = "|",collapse = "|"))
          
          Pathway2gene=convert_agregate(General_Kegg_Annotation,opt$GENE_ID_TYPE,"pathway","[\\|]")
          Pathway2gene=Pathway2gene[Pathway2gene$v1!='',]
          Pathway2gene=Pathway2gene[rev(colnames(Pathway2gene))]
          Pathway2name=Pathway_info[rev(colnames(Pathway_info))]
          KEGG_flag = TRUE
          
          write.table(Pathway2name,
                file = file.path(opt$outDir,paste('Pathway2name','tab', sep = '.')) ,
                quote = F,
                row.names = F,
                sep = "\t")
                
          write.table(Pathway2gene,
                file = file.path(opt$outDir,paste('Pathway2gene','tab', sep = '.')) ,
                quote = F,
                row.names = F,
                sep = "\t")
        }      
      }
    } else if (!is.na(opt$Trinotate)) {
        Annotation = read.csv(file = opt$Trinotate,sep = '\t')
        opt$GENE_ID_TYPE = colnames(Annotation)[1]
        
        
        Annotation$GO = apply(X = Annotation[,c('gene_ontology_blast','gene_ontology_pfam')],
                              FUN = function(x) unlist(paste(unique(x),collapse = ' ')),
                              MARGIN = 1)


        Annotation$GO_CC = sapply(X = Annotation$GO,
                                  FUN = function(X) stringi::stri_replace_all(str =  paste(unlist(unique(unlist(stringi::stri_extract(regex = 'GO:[0-9]+\\^cellular_component' ,
                                                                                                                                      str   = X ,
                                                                                                                                      mode = 'all')))),
                                                                                            collapse = '/'),
                                                                              fixed = '^cellular_component',
                                                                              replacement = '') )

        Annotation$GO_BP = sapply(X = Annotation$GO,
                                  FUN = function(X) stringi::stri_replace_all(str =  paste(unlist(unique(unlist(stringi::stri_extract(regex = 'GO:[0-9]+\\^biological_process' ,
                                                                                                                                      str   = X ,
                                                                                                                                      mode = 'all')))),
                                                                                           collapse = '/'),
                                                                              fixed = '^biological_process',
                                                                              replacement = '') )

        Annotation$GO_MF = sapply(X = Annotation$GO,
                                  FUN = function(X) stringi::stri_replace_all(str =  paste(unlist(unique(unlist(stringi::stri_extract(regex = 'GO:[0-9]+\\^molecular_function' ,
                                                                                                                                      str   = X ,
                                                                                                                                      mode = 'all')))),
                                                                                           collapse = '/'),
                                                                              fixed = '^molecular_function',
                                                                              replacement = '') )



        Annotation$go_id = sapply(X = Annotation$GO,
                                  FUN = function(X) stringi::stri_replace_all(str =  paste(unlist(unique(unlist(stringi::stri_extract(regex = 'GO:[0-9]+' ,
                                                                                                                                      str   = X ,
                                                                                                                                      mode = 'all')))),
                                                                                          collapse = '/'),
                                                                              fixed = '^',
                                                                              replacement = '') )

        Annotation$KO = sapply(X = Annotation$Kegg,
                                   FUN = function(X) stringi::stri_replace_all(str =  paste(unlist(stringi::stri_extract(regex = 'KO:K[0-9]+' ,
                                                                                                                         str   = X ,
                                                                                                                         mode = 'all')),
                                                                                            collapse = '/'),
                                                                               fixed = 'KO:',
                                                                               replacement = 'ko:') )
    }
    
    
    
    KEGG_KAAS_Data=NA
    if (!is.na(opt$KEGG_KAAS)) { 
      KEGG_KAAS_Data <- read.csv(opt$KEGG_KAAS,sep="\t",header = FALSE)
      colnames(KEGG_KAAS_Data) = c("Genes","KO")
      KEGG_KAAS_Data$KO = sapply(X =KEGG_KAAS_Data$KO, FUN = function(x) if (!is.na(x)) if (x!='')   paste( "ko:", x,collapse = "",sep = "") else NA else NA)
      if (!is.na(opt$GENE_ID_TYPE)) {
        KEGG_KAAS_Data$Genes = lapply(X =KEGG_KAAS_Data$Genes,FUN = function(x) unlist(stringi::stri_split(str = x,regex = "_"))[1] )
        KEGG_KAAS_Data$Genes = lapply(X =KEGG_KAAS_Data$Genes,FUN = function(x) rev(unlist(stringi::stri_split(str = x,regex = ":")))[1] )
      }
      
    }else if (!is.na(opt$Trinotate)) {
      KEGG_KAAS_Data <- Annotation[,c(opt$GENE_ID_TYPE,"KO")] 
      KEGG_KAAS_Data = convert_agregate(Annotation, opt$GENE_ID_TYPE ,"KO","/")
      colnames(KEGG_KAAS_Data) = c("Genes","KO")
      KEGG_KAAS_Data$KO = sapply(X =KEGG_KAAS_Data$KO, FUN = function(x) if (x=='')  NA else x)
      
    }
    
    if (!is.na(KEGG_KAAS_Data)) { 
      
      
      temp=unique(as.vector(na.omit(KEGG_KAAS_Data$KO)))
      temp=split(temp, ceiling(seq_along(temp)/50))
      
      pathway=unlist(sapply(X =temp,FUN = function(x) keggLink_retry("pathway",x)))
      
      pathway= data.frame(KO = names(pathway) , pathway = pathway)
      pathway$KO = sapply(X =pathway$KO, FUN = function(x)  paste( "ko:",  rev(unlist(stringi::stri_split(str = x,regex = ":")))[1] ,collapse = "",sep = "") )
      pathway = pathway[!stringi::stri_startswith(str = pathway$pathway ,fixed = 'path:ko'),]
      
      KEGG_KAAS_Data = merge.data.frame(KEGG_KAAS_Data,pathway,by.x ="KO" ,by.y = "KO",all.x = T)
      
      Pathway_info=keggList("pathway")
      info=names(Pathway_info)
      Pathway_info=as.data.frame(Pathway_info)
      colnames(Pathway_info)="pathway_info"
      Pathway_info["pathway"]=info
      KEGG_KAAS_Data = merge.data.frame(KEGG_KAAS_Data,Pathway_info,by = "pathway",all.x = T)
      
      KO_info=keggList("ko")
      info=names(KO_info)
      KO_info=as.data.frame(KO_info)
      colnames(KO_info)="ko_info"
      KO_info["KO"]=info
      
      KEGG_KAAS_Data = merge.data.frame(KEGG_KAAS_Data,KO_info,by = "KO",all.x = T)
      colnames(KEGG_KAAS_Data) = sapply(X =colnames(KEGG_KAAS_Data), FUN = function(x) if (x!='Genes')  paste( "KASS_", x,collapse = "",sep = "") else x)
      if (!is.na(Annotation)) {  
        Annotation = merge.data.frame(Annotation,KEGG_KAAS_Data,by.x =opt$GENE_ID_TYPE ,by.y = "Genes",all.x = T)
      }else{
        Annotation = data.frame(Genes = rownames(countData))
        Annotation = merge.data.frame(Annotation,KEGG_KAAS_Data,by = "Genes",all.x = T)
        opt$GENE_ID_TYPE = 'Genes'
      }
      
      KASS_Pathway2gene = na.omit(Annotation[c('KASS_pathway',opt$GENE_ID_TYPE)])
      KASS_Pathway2gene = KASS_Pathway2gene[!duplicated.data.frame(KASS_Pathway2gene),]
      
      KASS_Pathway2name = na.omit(Annotation[c('KASS_pathway','KASS_pathway_info')])
      KASS_Pathway2name = KASS_Pathway2name[!duplicated.data.frame(KASS_Pathway2name),]
      
      KEGG_KAAS_flag = TRUE
      
      write.table(KASS_Pathway2name,
                file = file.path(opt$outDir,paste('KASS_Pathway2name','tab', sep = '.')) ,
                quote = F,
                row.names = F,
                sep = "\t")
                
      write.table(KASS_Pathway2gene,
                file = file.path(opt$outDir,paste('KASS_Pathway2gene','tab', sep = '.')) ,
                quote = F,
                row.names = F,
                sep = "\t")
    }
    
  }else{
    cat("The Bioconductor KEGGREST package is not installed. You must install it for this script to work!")
  }
}

if (!is.na(Annotation)) { 
  aggregate_by=Annotation[opt$GENE_ID_TYPE]
  Annotation[opt$GENE_ID_TYPE]=NULL
  Annotation=aggregate.data.frame(x =Annotation ,
                                  by = aggregate_by,
                                  FUN =function(x) paste(unique(na.omit(x)),sep = "/",collapse = "/"))
  rownames(Annotation)=unlist(Annotation[opt$GENE_ID_TYPE])
  Annotation<-apply(X =Annotation,
                    MARGIN = c(1,2),
                    FUN =function(x) stringi::stri_replace_all(str = x,replacement = " ",regex = ",") )
  Annotation=as.data.frame(Annotation)
  
  
  write.table(Annotation,
                file = file.path(opt$outDir,paste('Annotation','tab', sep = '.')) ,
                quote = F,
                row.names = F,
                sep = "\t")
  
    # Get GO information:
    if ("package:clusterProfiler" %in% search()){      
        if ((Annotation_flag != 'OrgDb') & ("go_id" %in% colnames(Annotation))) {
          
          GO2gene           = convert_agregate(Annotation,opt$GENE_ID_TYPE,"go_id","/")
          GO2gene           = GO2gene[GO2gene$v1!='',]
          GO2gene           = sapply(GO2gene,FUN = function(x) stringi::stri_replace_all(str = x,replacement = "",regex = " "))
          GO2gene           = GO2gene[,c("v1","index")]
          #convert GO indirect to direct GOs
          GO2gene           = clusterProfiler::buildGOmap(GO2gene)
          # annotate the GO to TERMs
          Goterm            = merge.data.frame(x= clusterProfiler::go2term(GO2gene$GO), y=clusterProfiler::go2ont(GO2gene$GO) ,by='go_id')
          Goterm            = Goterm[!duplicated.data.frame(Goterm),]
          GO2name           = merge.data.frame(x = GO2gene,y = Goterm,by.y = "go_id",by.x ="GO" ,sort = FALSE)
          
          GO2name_MF        = GO2name[GO2name$Ontology=="MF",c('GO','Term')]
          GO2name_BP        = GO2name[GO2name$Ontology=="BP",c('GO','Term')]
          GO2name_MF        = GO2name_MF[!duplicated.data.frame(GO2name_MF),]
          GO2name_BP        = GO2name_BP[!duplicated.data.frame(GO2name_BP),]
          
          GO2gene_MF        = GO2name[GO2name$Ontology=="MF",c('GO','Gene')]
          GO2gene_BP        = GO2name[GO2name$Ontology=="BP",c('GO','Gene')]
          GO2gene_MF        = GO2gene_MF[!duplicated.data.frame(GO2gene_MF),]
          GO2gene_BP        = GO2gene_BP[!duplicated.data.frame(GO2gene_BP),]
          
          GO_flag           = TRUE
          write.table(GO2gene_MF,
                file = file.path(opt$outDir,paste('GO2gene_MF','tab', sep = '.')) ,
                quote = F,
                row.names = F,
                sep = "\t")
                
          write.table(GO2name_MF,
                file = file.path(opt$outDir,paste('GO2name_MF','tab', sep = '.')) ,
                quote = F,
                row.names = F,
                sep = "\t")
                
                
          write.table(GO2gene_BP,
                file = file.path(opt$outDir,paste('GO2gene_BP','tab', sep = '.')) ,
                quote = F,
                row.names = F,
                sep = "\t")
                
          write.table(GO2name_BP,
                file = file.path(opt$outDir,paste('GO2name_BP','tab', sep = '.')) ,
                quote = F,
                row.names = F,
                sep = "\t")
        }
    }
  
}
print('Reading samples data...')    
#Read Sample Data
colData <- read.csv(opt$SAMPLE_DATA_FILE, sep="\t", row.names=1)
rownames(colData) = make.names(rownames(colData))

used_samples = intersect(rownames(colData),colnames(countData))
colData_col <- colnames(colData)
countData <- countData[, used_samples]
colData   <- as.data.frame(colData[used_samples ,])
rownames(colData) <- used_samples
colnames(colData) <- colData_col
all(rownames(colData) == colnames(countData))
Library_sizes=apply(X = countData,MARGIN = 2,FUN = sum)
colData=merge(colData,as.data.frame(Library_sizes),by="row.names",sort=F)
rownames(colData)=colData$Row.names
colData$Row.names=NULL
countData <- countData[, rownames(colData)]
print('Done reading samples data')  


##############Filttering the data##############################

if ((opt$FILTER_SAMPLES) | (opt$FILTER_GENES)){
  
  # Check if required packages are installed:
  if(!(all(c('scater') %in% installed.packages()))) {
    if (Sys.getenv("CONDA_PREFIX")!=""){
      source("https://bioconductor.org/biocLite.R")
      biocLite('scater')
    }else{
      cat("The Bioconductor scater package is not installed. You must install it for this script to work!")
    }
  } 
  library(scater)
  
  
  sce <- newSCESet(countData=countData)
  dim(sce)
  sce <- calculateQCMetrics(sce)
  pdf(file = file.path(opt$outDir,"Pre_QA_hist.pdf"))
  par(mfrow=c(1,2))
  hist(sce$total_counts/1e6, xlab="Library sizes (millions)", main="", 
       breaks=20, col="grey80", ylab="Number of Samples")
  hist(sce$total_features, xlab="Number of expressed genes", main="", 
       breaks=20, col="grey80", ylab="Number of Samples")
  dev.off() 
  
  if (opt$FILTER_SAMPLES){
    print('Filtering Samples')
    libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
    feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE)
    sce <- sce[,!(libsize.drop | feature.drop )]
    print.data.frame(data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), Remaining=ncol(sce)))
    use_names=colnames(sce)
    countData=countData[,use_names]
    colData=colData[use_names,]
  }
  
  if (opt$FILTER_GENES){
    print('Filtering low-abundance genes')
    print('Before Filtering:')
    print(dim(sce))
    ave.counts <- rowMeans(counts(sce))
    keep <- ave.counts >= 1
    #print(sum(keep))
    sce <- sce[names(keep[keep==T]),]
    print('After Filtering:')
    print(dim(sce))
    countData=countData[names(keep[keep==T]),]
  }
  
  pdf(file = file.path(opt$outDir,"Post_QA_hist.pdf"))
  par(mfrow=c(1,2))
  hist(sce$total_counts/1e6, xlab="Library sizes (millions)", main="", 
       breaks=20, col="grey80", ylab="Number of Samples")
  hist(sce$total_features, xlab="Number of expressed genes", main="", 
       breaks=20, col="grey80", ylab="Number of Samples")
  
  dev.off()
}        
##################Filttering the data END######################

##################DeSeq- design and normalization ########################################

if (opt$only_clustering==FALSE){
        if (is.na(opt$DESIGN)){
            opt$DESIGN   = '~ 1'
            opt$BLIND    = T
            opt$CONTRAST = 'No_DESIGN'
            print('No DESIGN found, will only normalize and cluster all genes')
        }
        if (opt$COUNT_SOURCE=='RSEM'){
          zero_len           = unique( rownames(which(txi.rsem$length == 0,arr.ind = TRUE)))
          countData          = countData[!( rownames(countData) %in% zero_len ),]
          txi.rsem$counts    = countData
          txi.rsem$abundance = txi.rsem$abundance[row.names(countData),colnames(countData)]
          txi.rsem$length    = txi.rsem$length[row.names(countData),colnames(countData)]
          
          DESeqDataSet_base <- DESeqDataSetFromTximport(txi = txi.rsem,
                                                   colData = colData,
                                                   design  = as.formula(opt$DESIGN))
        }else{ 
          DESeqDataSet_base <- DESeqDataSetFromMatrix(countData = countData,
                                                 colData = colData,
                                                 design = as.formula(opt$DESIGN))
        }
        if (opt$modelMatrixType=='standard'){ 
          DESeqDataSet_base <- DESeq(DESeqDataSet_base,parallel=F)
        }else{
          DESeqDataSet_base <- DESeq(DESeqDataSet_base,
                                parallel=F,
                                betaPrior = TRUE,
                                modelMatrixType=opt$modelMatrixType)
        }

        test_count=list()

        if (!is.na(opt$CONTRAST)){
            test_count = unlist(stringi::stri_split(str = opt$CONTRAST ,fixed = "|"))
        }

        if (!is.na(opt$LRT)){
            test_count=c('LTR',test_count)
        }
        
        if (!is.na(opt$collapseReplicates)){
            if (opt$collapseReplicates %in% colnames(colData)){
                collapsed_DESeqDataSet <- try( collapseReplicates(DESeqDataSet_base,DESeqDataSet_base[[opt$collapseReplicates]] ) ,silent = T)
                if (inherits(combat_edata,"try-error")){
                    print('Failed to collapse Replicates !!!')
                }else{
                    collapsed_DESeqDataSet$Library_sizes=apply(X = as.data.frame(counts(collapsed_DESeqDataSet)),MARGIN = 2,FUN = sum)
                    DESeqDataSet_base = collapsed_DESeqDataSet
                    colData = as.data.frame(colData(collapsed_DESeqDataSet))
                    countData = as.data.frame(counts(collapsed_DESeqDataSet))
                }
            }else{
                 print('The field to be used for Replicates collapse is not found in the Sample Data !!!')
            }
        }
        
}else{
    DESeqDataSet_base=NA
    test_count=c('only_clustering')
}





base_out_dir=opt$outDir

Rmarkdown_flag=FALSE
if (!is.na(opt$Rmarkdown)){
    if((all(c('plotly','DT','clusterProfiler','rmarkdown') %in% installed.packages()))) {
      Rmarkdown_flag=TRUE
    } else {
      cat("In order to use the Rmarkdown option all of the following R packages need to be installed: 'plotly','DT','clusterProfiler','rmarkdown'")  
    }
}

save.image(file.path(base_out_dir,'Session.RSession'))

Original_colData   = colData
Original_countData = countData

for (test2do in test_count){
    print(test2do)
    if (Rmarkdown_flag){
        if ((test2do=='LTR')||(test2do=='No_DESIGN')||(test2do=='only_clustering') ){
            test_name=test2do
        }else{
            contrast_list = unlist(stringi::stri_split(str = test2do,regex = ','))
            test_name = paste(contrast_list[2],contrast_list[3],sep = '_vs_')
        }
        rmarkdown::render(input = opt$Rmarkdown,
                          output_file = file.path(base_out_dir,paste('Final_Report',test_name,'html', sep = '.')),
                          intermediates_dir = base_out_dir,
                          clean = TRUE,
                          params = list(
                                        opt = opt,
                                        test2do = test2do,
                                        DESeqDataSet_base = DESeqDataSet_base,
                                        base_out_dir = base_out_dir,
                                        colData = colData,
                                        countData = countData,
                                        Annotation = Annotation,
                                        Annotation_flag = Annotation_flag,
                                        KEGG_flag = KEGG_flag,
                                        ORGANISM_KEGG_flag = ORGANISM_KEGG_flag,
                                        KEGG_KAAS_flag = KEGG_KAAS_flag,
                                        GO_flag = GO_flag
                          ))
        if (test2do=='LTR'){ 
           opt$LRT=NA 
        }    
    } else {
        
        
        if (!is.na(opt$Post_statistical_ALPHA)){
            FDR_PVAL_CUTOFF = opt$Post_statistical_ALPHA
        }else{
            FDR_PVAL_CUTOFF = opt$ALPHA
        }
        
        if (!is.na(opt$Post_statistical_FoldChange)){
            FC_CUTOFF = opt$Post_statistical_FoldChange #in linear scale
            FC_CUTOFF_log2 = log2(opt$Post_statistical_FoldChange) 
        }else{
            FC_CUTOFF = opt$FoldChange #in linear scale
            FC_CUTOFF_log2 = log2(opt$FoldChange) 
        }
        
        
        DESeqDataSet_Results = NA
        Normalized_counts = NA
        Normalized_counts_list=c()
        DESeqDataSet=DESeqDataSet_base
        if (test2do=='LTR'){
          opt$outDir=file.path(base_out_dir,'LTR')
          dir.create(opt$outDir, showWarnings = FALSE)
          
          DESeqDataSet <- DESeq(DESeqDataSet, test="LRT",
                                reduced = as.formula(opt$LRT),
                                parallel=F)
          
          DESeqDataSet_Results <- results(DESeqDataSet,alpha = opt$ALPHA)
          
        }else if (test2do=='No_DESIGN') {
          opt$outDir=file.path(base_out_dir,'No_DESIGN')
          dir.create(opt$outDir, showWarnings = FALSE)
        }else if (test2do=='only_clustering'){
          opt$outDir=file.path(base_out_dir,'Only_Clustering')
          dir.create(opt$outDir, showWarnings = FALSE)
          Normalized_counts_assay = countData
        }else{
            opt$CONTRAST=test2do
            contrast_list = unlist(stringi::stri_split(str =  opt$CONTRAST,regex = ','))
            test_name = paste(contrast_list[2],contrast_list[3],sep = '_vs_')
            opt$outDir=file.path(base_out_dir,test_name)
            dir.create(opt$outDir, showWarnings = FALSE)
            
            DESeqDataSet_Results <- results(DESeqDataSet,
                                            contrast     = unlist(stringi::stri_split(str = opt$CONTRAST ,regex = ",")) ,
                                            alpha        = opt$ALPHA,
                                            lfcThreshold = log2(opt$FoldChange)   )
          
        }

        if (!is.na(DESeqDataSet_Results)){
            DESeqDataSet_Results_redOrdered <- as.matrix(DESeqDataSet_Results[order(DESeqDataSet_Results$padj),])
            summary(DESeqDataSet_Results)
        }

        # Normalization
        if (!is.na(DESeqDataSet)){
            if (opt$NORMALIZATION_TYPE=='RLOG'){
                Normalized_counts <- rlog(DESeqDataSet, blind=opt$BLIND)
            }else{
                Normalized_counts <- varianceStabilizingTransformation(DESeqDataSet, blind=opt$BLIND)
            }
            Normalized_counts_list=list(Normalized_counts)
            if (!is.na(opt$removeBatchEffect)){
                if (opt$removeBatchEffect_method=="sva"){
                    if(!(all(c('sva') %in% installed.packages()))) {
                        if (Sys.getenv("CONDA_PREFIX")!=""){
                            source("http://bioconductor.org/biocLite.R")
                            try(biocLite("sva", version = "3.8"),silent = F)
                        }else{
                          cat("The Bioconductor sva package is not installed. You must install it for this script to work!")
                        }
                    } 
                    try(library("sva"),silent = F)
                    if ("package:sva" %in% search()==FALSE){
                        cat(sprintf("%s","The Bioconductor sva package is not installed!! will use limma instead"))
                        opt$removeBatchEffect_method="limma"
                    }
                }
                BatchEffects = unlist(stringi::stri_split(str = opt$removeBatchEffect ,regex = ','))
                if (opt$removeBatchEffect_method=="sva"){
                        if (length(BatchEffects)>1){
                            cat(sprintf("%s","It is Not possible to use more then one batch effect using the sva package, will use only the first one!!"))
                        }
                        DESIGN_list = unlist(stringi::stri_split( stringi::stri_replace(str = opt$DESIGN,
                                                                                   regex =' |~',
                                                                                   mode = 'all',
                                                                                   replacement = ''),
                                                            fixed = '+'))
                        New_DESIGN =paste0('~',paste(DESIGN_list[DESIGN_list!=BatchEffects[1]],collapse ='+'))
                        
                        modcombat = model.matrix( as.formula(New_DESIGN) , data=Original_colData)
                        combat_edata<-try(ComBat(dat = assay(Normalized_counts),
                                              batch = Original_colData[,BatchEffects[1]],
                                              mod = modcombat,
                                              par.prior = T,
                                              mean.only=F),
                                          silent = T)
                        if (inherits(combat_edata,"try-error")){
                            cat(sprintf("%s","Failed to correct batch effect using parametric adjustment will try using limma instead!!\n"))
                           
                            assay(Normalized_counts) = limma::removeBatchEffect(assay(Normalized_counts),
                                                                                batch = Normalized_counts@colData[,BatchEffects[1]])
                            
                        }else{
                            assay(Normalized_counts) = combat_edata
                        }   
                        
                }else{
                
                    if (length(BatchEffects)>1){
                            assay(Normalized_counts) = limma::removeBatchEffect(assay(Normalized_counts),
                                                                                batch = Normalized_counts@colData[,BatchEffects[1]],
                                                                                batch2= Normalized_counts@colData[,BatchEffects[2]])
                    }else{
                            assay(Normalized_counts) = limma::removeBatchEffect(assay(Normalized_counts),
                                                                                batch = Normalized_counts@colData[,BatchEffects[1]])
                    }
                }
                Normalized_counts_list = c(Normalized_counts_list,Normalized_counts)
            }
            
            Normalized_counts_assay = assay(Normalized_counts)
        }



        ##################DeSeq-END####################################



        ##################DeSeq-Print Results As Is (Vered 14.6.2018)##
        if (!is.na(DESeqDataSet_Results)){
            Results_for_print = as.data.frame(DESeqDataSet_Results)

            if (!is.na(opt$LRT)){
                Results_for_print = Results_for_print[,c('pvalue','padj')]
                test_name = 'LRT'

            }else{
                if (!is.na(opt$CONTRAST)){
                    Results_for_print = Results_for_print[,c('log2FoldChange','pvalue','padj')]
                    contrast_list = unlist(stringi::stri_split(str =  opt$CONTRAST,regex = ','))
                    test_name = paste(contrast_list[2],contrast_list[3],sep = '.vs.')
                }
            }

            colnames(Results_for_print) = unlist(lapply(colnames(Results_for_print), function(x)  paste(x, test_name,sep = '.')) ) 
            write.table(Results_for_print,
                        file = file.path(opt$outDir,paste(test_name,'txt', sep = '.')) ,
                        quote = F,
                        sep = "\t",
                        col.names=NA)
        }
        ##################DeSeq-Print Results END######################

        ##################Print Normalized counts######################
        if (!is.na(Normalized_counts)){
            write.table(Normalized_counts_assay,
                        file = file.path(opt$outDir,paste(test_name,opt$NORMALIZATION_TYPE,'Normalized_counts.txt', sep = '_')) ,
                        quote = F,
                        sep = "\t")
        }
        ##################Print Normalized counts END#######################

        ##################Volcano plot and MA plot Vered##################

        if (!is.na(DESeqDataSet_Results)){
           
            if (!is.na(opt$LRT)){
              print('No Volcano plot for LTR')
            }else{
                if (!is.na(opt$CONTRAST)){
                    contrast_list = unlist(stringi::stri_split(str =  opt$CONTRAST,regex = ','))
                    test_name = paste(contrast_list[2],contrast_list[3],sep = ' vs ')
                    plotTitle = paste('Volcano Plot (', test_name, ')', sep="")
                    pdf(file.path(opt$outDir,paste(c(plotTitle,".pdf"),collapse ="")))
                    par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)

                    topT <- as.data.frame(DESeqDataSet_Results)

                    #Adjusted P values (FDR Q values)
                    with(topT, plot(log2FoldChange, -log10(padj), pch=20, main=plotTitle, cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~'FDR p-'~value)))

                    with(subset(topT, padj<FDR_PVAL_CUTOFF & abs(log2FoldChange)>FC_CUTOFF_log2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))

                    #Add lines for absolute adjusted p-value and FC cutoffs 
                    #abline(v=0, col="black", lty=3, lwd=1.0)  #adds a vertical line at 0
                    abline(v=-FC_CUTOFF_log2, col="black", lty=4, lwd=2.0)
                    abline(v=FC_CUTOFF_log2, col="black", lty=4, lwd=2.0)
                    if (length(topT$padj[topT$padj<FDR_PVAL_CUTOFF])>0){
                        abline(h=-log10(max(topT$padj[topT$padj<FDR_PVAL_CUTOFF], na.rm=TRUE)), col="black", lty=4, lwd=FC_CUTOFF_log2)
                    }
                    dev.off() 
                }
            }

            #plot MA

            plotTitle = paste('MA Plot (', test_name, ')', sep="")
            pdf(file.path(opt$outDir,paste(c(plotTitle,".pdf"),collapse ="")))
            plotMA(DESeqDataSet_Results, alpha=FDR_PVAL_CUTOFF, colNonSig = 'blue', ylim=c(-10,10), main=plotTitle)
            abline(h=c(-FC_CUTOFF_log2,FC_CUTOFF_log2), col="red")
            dev.off() 
        }
        
        ############################Original-PCA-START##############################################
        
        if (length(Normalized_counts_list)==2){
            PCA_Normalized_counts = Normalized_counts_list[[1]]
        }else if (length(Normalized_counts_list)==1){
            PCA_Normalized_counts = Normalized_counts_list[[1]]
        }else if (length(Normalized_counts_list)==0){
            PCA_Normalized_counts = NA
        }
        ##################PCA-TOP500####################################
        if (!is.na(PCA_Normalized_counts)){
            intgroup=c(opt$PCA_COLOR,opt$PCA_SHAPE,opt$PCA_SIZE)
            intgroup=intgroup[!is.na(intgroup)]
            PCA_data <- plotPCA(PCA_Normalized_counts, intgroup=make.names(intgroup), returnData=TRUE,ntop=500)
            percentVar <- round(100 * attr(PCA_data, "percentVar"))

            if (is.na(opt$PCA_COLOR)){
              opt$PCA_COLOR='name'
            }

            if (is.na(opt$PCA_SHAPE)){
                pca_res=ggplot(PCA_data, aes_string('PC1', 'PC2', color=make.names(opt$PCA_COLOR),size=make.names(opt$PCA_SIZE)) )+
                    geom_count()+
                    scale_shape_discrete(solid=F)+
                    guides(color=guide_legend(order = 1 ,title=opt$PCA_COLOR))+
                    guides(size=guide_legend(order  = 2,title=opt$PCA_SIZE))+
                    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
                    coord_fixed()
                ggsave(file.path(opt$outDir,paste(c('TOP500_PCA_pc1_pc2',".pdf"),collapse ="")),pca_res,dpi = 600, width = 6.99, height =6.99, units = "in")

                pca_res=ggplot(PCA_data, aes_string('PC1', 'PC3', color=make.names(opt$PCA_COLOR),size=make.names(opt$PCA_SIZE)) )+
                    geom_count()+
                    guides(color=guide_legend(order = 1 ,title=opt$PCA_COLOR))+
                    guides(size=guide_legend(order  = 2,title=opt$PCA_SIZE))+
                    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                    ylab(paste0("PC3: ",percentVar[3],"% variance")) +
                    coord_fixed()
                ggsave(file.path(opt$outDir,paste(c('TOP500_PCA_pc1_pc3',".pdf"),collapse ="")),pca_res,dpi = 600, width = 6.99, height =6.99, units = "in")
              
            }else{
                pca_res=ggplot(PCA_data, aes_string('PC1', 'PC2', color=make.names(opt$PCA_COLOR),shape=make.names(opt$PCA_SHAPE),size=make.names(opt$PCA_SIZE)) )+
                    geom_count()+
                    scale_shape_discrete(solid=F)+
                    guides(color=guide_legend(order = 1 ,title=opt$PCA_COLOR))+
                    guides(shape=guide_legend(order  = 2,title=opt$PCA_SHAPE))+
                    guides(size=guide_legend(order  = 2,title=opt$PCA_SIZE))+
                    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
                    coord_fixed()
                ggsave(file.path(opt$outDir,paste(c('TOP500_PCA_pc1_pc2',".pdf"),collapse ="")),pca_res,dpi = 600, width = 6.99, height =6.99, units = "in")

                pca_res=ggplot(PCA_data, aes_string('PC1', 'PC3', color=make.names(opt$PCA_COLOR),shape=make.names(opt$PCA_SHAPE),size=make.names(opt$PCA_SIZE)) )+
                    geom_count()+
                    guides(color=guide_legend(order = 1 ,title=opt$PCA_COLOR))+
                    guides(shape=guide_legend(order  = 2,title=opt$PCA_SHAPE))+
                    guides(size=guide_legend(order  = 2,title=opt$PCA_SIZE))+
                    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                    ylab(paste0("PC3: ",percentVar[3],"% variance")) +
                    coord_fixed()
                ggsave(file.path(opt$outDir,paste(c('TOP500_PCA_pc1_pc3',".pdf"),collapse ="")),pca_res,dpi = 600, width = 6.99, height =6.99, units = "in")
            }
        }
        ##################PCA-TOP500-END####################################

        ##################PCA-ALL####################################
        if (!is.na(PCA_Normalized_counts)){
            PCA_data <- plotPCA(PCA_Normalized_counts, intgroup=make.names(intgroup), returnData=TRUE,ntop='all')
            percentVar <- round(100 * attr(PCA_data, "percentVar"))

            if (is.na(opt$PCA_COLOR)){
              opt$PCA_COLOR='name'
            }

            if (is.na(opt$PCA_SHAPE)){
                pca_res=ggplot(PCA_data, aes_string('PC1', 'PC2', color=make.names(opt$PCA_COLOR),size=make.names(opt$PCA_SIZE)) )+
                    geom_count()+
                    scale_shape_discrete(solid=F)+
                    guides(color=guide_legend(order = 1 ,title=opt$PCA_COLOR))+
                    guides(size=guide_legend(order  = 2,title=opt$PCA_SIZE))+
                    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
                    coord_fixed()
                ggsave(file.path(opt$outDir,paste(c('ALL_PCA_pc1_pc2',".pdf"),collapse ="")),pca_res,dpi = 600, width = 6.99, height =6.99, units = "in")

                pca_res=ggplot(PCA_data, aes_string('PC1', 'PC3', color=make.names(opt$PCA_COLOR),size=make.names(opt$PCA_SIZE)) )+
                    geom_count()+
                    guides(color=guide_legend(order = 1 ,title=opt$PCA_COLOR))+
                    guides(size=guide_legend(order  = 2,title=opt$PCA_SIZE))+
                    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                    ylab(paste0("PC3: ",percentVar[3],"% variance")) +
                    coord_fixed()
                ggsave(file.path(opt$outDir,paste(c('ALL_PCA_pc1_pc3',".pdf"),collapse ="")),pca_res,dpi = 600, width = 6.99, height =6.99, units = "in")
              
            }else{
                pca_res=ggplot(PCA_data, aes_string('PC1', 'PC2', color=make.names(opt$PCA_COLOR),shape=make.names(opt$PCA_SHAPE),size=make.names(opt$PCA_SIZE)) )+
                    geom_count()+
                    scale_shape_discrete(solid=F)+
                    guides(color=guide_legend(order = 1 ,title=opt$PCA_COLOR))+
                    guides(shape=guide_legend(order  = 2,title=opt$PCA_SHAPE))+
                    guides(size=guide_legend(order  = 2,title=opt$PCA_SIZE))+
                    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
                    coord_fixed()
                ggsave(file.path(opt$outDir,paste(c('ALL_PCA_pc1_pc2',".pdf"),collapse ="")),pca_res,dpi = 600, width = 6.99, height =6.99, units = "in")
              
                pca_res=ggplot(PCA_data, aes_string('PC1', 'PC3', color=make.names(opt$PCA_COLOR),shape=make.names(opt$PCA_SHAPE),size=make.names(opt$PCA_SIZE)) )+
                    geom_count()+
                    guides(color=guide_legend(order = 1 ,title=opt$PCA_COLOR))+
                    guides(shape=guide_legend(order  = 2,title=opt$PCA_SHAPE))+
                    guides(size=guide_legend(order  = 2,title=opt$PCA_SIZE))+
                    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                    ylab(paste0("PC3: ",percentVar[3],"% variance")) +
                    coord_fixed()
                ggsave(file.path(opt$outDir,paste(c('ALL_PCA_pc1_pc3',".pdf"),collapse ="")),pca_res,dpi = 600, width = 6.99, height =6.99, units = "in")
            }
        }
        ##################PCA-END####################################
        
    ############################Original-PCA-ENDS##############################################


   ############################Batch-Effect-Corrected-PCA-START##############################################
    
    if (length(Normalized_counts_list)==2){
        PCA_Normalized_counts = Normalized_counts_list[[2]]
        extra_string="Batch_Effect_Corrected_"
    }        
        ##################PCA-TOP500####################################
        if (length(Normalized_counts_list)==2){
            intgroup=c(opt$PCA_COLOR,opt$PCA_SHAPE,opt$PCA_SIZE)
            intgroup=intgroup[!is.na(intgroup)]
            PCA_data <- plotPCA(PCA_Normalized_counts, intgroup=make.names(intgroup), returnData=TRUE,ntop=500)
            percentVar <- round(100 * attr(PCA_data, "percentVar"))

            if (is.na(opt$PCA_COLOR)){
              opt$PCA_COLOR='name'
            }

            if (is.na(opt$PCA_SHAPE)){
                pca_res=ggplot(PCA_data, aes_string('PC1', 'PC2', color=make.names(opt$PCA_COLOR),size=make.names(opt$PCA_SIZE)) )+
                    geom_count()+
                    scale_shape_discrete(solid=F)+
                    guides(color=guide_legend(order = 1 ,title=opt$PCA_COLOR))+
                    guides(size=guide_legend(order  = 2,title=opt$PCA_SIZE))+
                    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
                    coord_fixed()
                ggsave(file.path(opt$outDir,paste(c(extra_string,'TOP500_PCA_pc1_pc2',".pdf"),collapse ="")),pca_res,dpi = 600, width = 6.99, height =6.99, units = "in")

                pca_res=ggplot(PCA_data, aes_string('PC1', 'PC3', color=make.names(opt$PCA_COLOR),size=make.names(opt$PCA_SIZE)) )+
                    geom_count()+
                    guides(color=guide_legend(order = 1 ,title=opt$PCA_COLOR))+
                    guides(size=guide_legend(order  = 2,title=opt$PCA_SIZE))+
                    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                    ylab(paste0("PC3: ",percentVar[3],"% variance")) +
                    coord_fixed()
                ggsave(file.path(opt$outDir,paste(c(extra_string,'TOP500_PCA_pc1_pc3',".pdf"),collapse ="")),pca_res,dpi = 600, width = 6.99, height =6.99, units = "in")
              
            }else{
                pca_res=ggplot(PCA_data, aes_string('PC1', 'PC2', color=make.names(opt$PCA_COLOR),shape=make.names(opt$PCA_SHAPE),size=make.names(opt$PCA_SIZE)) )+
                    geom_count()+
                    scale_shape_discrete(solid=F)+
                    guides(color=guide_legend(order = 1 ,title=opt$PCA_COLOR))+
                    guides(shape=guide_legend(order  = 2,title=opt$PCA_SHAPE))+
                    guides(size=guide_legend(order  = 2,title=opt$PCA_SIZE))+
                    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
                    coord_fixed()
                ggsave(file.path(opt$outDir,paste(c(extra_string,'TOP500_PCA_pc1_pc2',".pdf"),collapse ="")),pca_res,dpi = 600, width = 6.99, height =6.99, units = "in")

                pca_res=ggplot(PCA_data, aes_string('PC1', 'PC3', color=make.names(opt$PCA_COLOR),shape=make.names(opt$PCA_SHAPE),size=make.names(opt$PCA_SIZE)) )+
                    geom_count()+
                    guides(color=guide_legend(order = 1 ,title=opt$PCA_COLOR))+
                    guides(shape=guide_legend(order  = 2,title=opt$PCA_SHAPE))+
                    guides(size=guide_legend(order  = 2,title=opt$PCA_SIZE))+
                    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                    ylab(paste0("PC3: ",percentVar[3],"% variance")) +
                    coord_fixed()
                ggsave(file.path(opt$outDir,paste(c(extra_string,'TOP500_PCA_pc1_pc3',".pdf"),collapse ="")),pca_res,dpi = 600, width = 6.99, height =6.99, units = "in")
            }
        }
        ##################PCA-TOP500-END####################################

        ##################PCA-ALL####################################
         if (length(Normalized_counts_list)==2){
            PCA_data <- plotPCA(PCA_Normalized_counts, intgroup=make.names(intgroup), returnData=TRUE,ntop='all')
            percentVar <- round(100 * attr(PCA_data, "percentVar"))

            if (is.na(opt$PCA_COLOR)){
              opt$PCA_COLOR='name'
            }

            if (is.na(opt$PCA_SHAPE)){
                pca_res=ggplot(PCA_data, aes_string('PC1', 'PC2', color=make.names(opt$PCA_COLOR),size=make.names(opt$PCA_SIZE)) )+
                    geom_count()+
                    scale_shape_discrete(solid=F)+
                    guides(color=guide_legend(order = 1 ,title=opt$PCA_COLOR))+
                    guides(size=guide_legend(order  = 2,title=opt$PCA_SIZE))+
                    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
                    coord_fixed()
                ggsave(file.path(opt$outDir,paste(c(extra_string,'ALL_PCA_pc1_pc2',".pdf"),collapse ="")),pca_res,dpi = 600, width = 6.99, height =6.99, units = "in")

                pca_res=ggplot(PCA_data, aes_string('PC1', 'PC3', color=make.names(opt$PCA_COLOR),size=make.names(opt$PCA_SIZE)) )+
                    geom_count()+
                    guides(color=guide_legend(order = 1 ,title=opt$PCA_COLOR))+
                    guides(size=guide_legend(order  = 2,title=opt$PCA_SIZE))+
                    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                    ylab(paste0("PC3: ",percentVar[3],"% variance")) +
                    coord_fixed()
                ggsave(file.path(opt$outDir,paste(c(extra_string,'ALL_PCA_pc1_pc3',".pdf"),collapse ="")),pca_res,dpi = 600, width = 6.99, height =6.99, units = "in")
              
            }else{
                pca_res=ggplot(PCA_data, aes_string('PC1', 'PC2', color=make.names(opt$PCA_COLOR),shape=make.names(opt$PCA_SHAPE),size=make.names(opt$PCA_SIZE)) )+
                    geom_count()+
                    scale_shape_discrete(solid=F)+
                    guides(color=guide_legend(order = 1 ,title=opt$PCA_COLOR))+
                    guides(shape=guide_legend(order  = 2,title=opt$PCA_SHAPE))+
                    guides(size=guide_legend(order  = 2,title=opt$PCA_SIZE))+
                    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
                    coord_fixed()
                ggsave(file.path(opt$outDir,paste(c(extra_string,'ALL_PCA_pc1_pc2',".pdf"),collapse ="")),pca_res,dpi = 600, width = 6.99, height =6.99, units = "in")
              
                pca_res=ggplot(PCA_data, aes_string('PC1', 'PC3', color=make.names(opt$PCA_COLOR),shape=make.names(opt$PCA_SHAPE),size=make.names(opt$PCA_SIZE)) )+
                    geom_count()+
                    guides(color=guide_legend(order = 1 ,title=opt$PCA_COLOR))+
                    guides(shape=guide_legend(order  = 2,title=opt$PCA_SHAPE))+
                    guides(size=guide_legend(order  = 2,title=opt$PCA_SIZE))+
                    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                    ylab(paste0("PC3: ",percentVar[3],"% variance")) +
                    coord_fixed()
                ggsave(file.path(opt$outDir,paste(c(extra_string,'ALL_PCA_pc1_pc3',".pdf"),collapse ="")),pca_res,dpi = 600, width = 6.99, height =6.99, units = "in")
            }
        }
        ##################PCA-END####################################
        
    ############################Batch-Effect-Corrected-PCA-ENDS##############################################

        

        ##################Plot Counts################################    
        if (!is.na(Normalized_counts)){
            if (!is.na(opt$GENES_PLOT)&(!is.na(opt$X_AXIS)   )){
                intgroup=c(opt$X_AXIS,opt$GROUP)
                intgroup=intgroup[!is.na(intgroup)]
                X_AXIS=opt$X_AXIS
                if (is.na(opt$GROUP)){
                    GROUP=opt$X_AXIS
                }else{
                    GROUP=opt$GROUP
                }
                if (is.na(opt$SPLIT_BY)){
                    for (gene in unlist(stringi::stri_split(str =  opt$GENES_PLOT,regex = ","))){
                        data=plotCounts(DESeqDataSet, gene=gene,transform=F, intgroup=make.names(intgroup),returnData=TRUE)
                        p=ggplot(data, aes_string(make.names(X_AXIS), 'count',color=make.names(GROUP), group=make.names(GROUP))) +
                        geom_point(size=4)+
                        scale_shape_discrete(solid=F)+
                        geom_smooth(method="gam", fullrange=F, size=0.5)
                        ggsave(file.path(opt$outDir,paste(c('Count_plot_',gene,".pdf"),collapse ="")),p,dpi = 600, width = 6.99, height =6.99, units = "in")
                    }
                }else{
                    for (gene in unlist(stringi::stri_split(str =  opt$GENES_PLOT,regex = ","))){
                        data=plotCounts(DESeqDataSet, gene=gene,transform=F, intgroup=make.names(intgroup),returnData=TRUE)
                        p=ggplot(data, aes_string(make.names(X_AXIS), 'count',color=make.names(GROUP), group=make.names(GROUP))) +
                        facet_wrap(~opt$SPLIT_BY, ncol=1,scales = "free_y",strip.position = c("right"))+
                        geom_point(size=4)+
                        scale_shape_discrete(solid=F)+
                        stat_smooth(method="auto", fullrange=F, size=0.5)
                        ggsave(file.path(opt$outDir,paste(c('Count_plot_',gene,".pdf"),collapse ="")),p,dpi = 600, width = 6.99, height =6.99, units = "in")
                    }
                }
            }
        }
        ##################Plot Counts-END#############################    

        ##################Clustering##################################    
        if (!is.na(opt$significant_genes)){ 
            sig.genes = unlist(stringi::stri_split(str = opt$significant_genes,regex = ','))
            Normalized_counts_assay = Normalized_counts_assay[sig.genes,]
        }else if (!is.na(DESeqDataSet_Results)){
            #sig.genes=rownames(DESeqDataSet_Results_redOrdered[(DESeqDataSet_Results_redOrdered[,"padj"]<opt$ALPHA)&(!is.na(DESeqDataSet_Results_redOrdered[,"padj"])),])
            up_reg    = rownames(DESeqDataSet_Results_redOrdered[(DESeqDataSet_Results_redOrdered[,"log2FoldChange"]>FC_CUTOFF_log2)&(DESeqDataSet_Results_redOrdered[,"padj"]<FDR_PVAL_CUTOFF)&(!is.na(DESeqDataSet_Results_redOrdered[,"padj"])),])
            down_reg  = rownames(DESeqDataSet_Results_redOrdered[(DESeqDataSet_Results_redOrdered[,"log2FoldChange"]<FC_CUTOFF_log2)&(DESeqDataSet_Results_redOrdered[,"padj"]<FDR_PVAL_CUTOFF)&(!is.na(DESeqDataSet_Results_redOrdered[,"padj"])),])
            sig.genes = c(up_reg,down_reg)
        }else{
            sig.genes=rownames(countData)
        }

        if (!is.na(opt$X_AXIS)){

            if (length(sig.genes)>0){

                X_AXIS=opt$X_AXIS
                if (is.na(opt$GROUP)){
                    GROUP=X_AXIS
                }else{
                    GROUP=opt$GROUP
                }
                
                if (is.na(opt$SPLIT_BY)){
                    if (!is.na(opt$CONTRAST)){
                        if (opt$SPLIT_BY_CONTRAST){
                            contrast_list_split = unlist(stringi::stri_split(str =  opt$CONTRAST,regex = ','))
                            if (length(contrast_list_split)==3){
                                colData = Original_colData[(Original_colData[contrast_list_split[1]]==contrast_list_split[2]) | (Original_colData[contrast_list_split[1]]==contrast_list_split[3]),]
                                
                                if (length(unique(colData[,X_AXIS]))==1){
                                    if (length(unique(colData[,GROUP]))>1){
                                        temp_X_AXIS = X_AXIS
                                        X_AXIS      = GROUP
                                        GROUP       = temp_X_AXIS
                                        
                                    }
                                }
                            }
                        }
                    }
                }
                
                
                if (X_AXIS!=GROUP){
                    if (is.na(opt$SPLIT_BY)){
                        Normalized_counts_Annotated=merge(colData[c(X_AXIS,GROUP)],t(Normalized_counts_assay),all.x = T  ,by="row.names",sort=F)
                        SPLIT_BY=opt$SPLIT_BY
                        Normalized_counts_Annotated_mean=aggregate.data.frame(x = Normalized_counts_Annotated,by = c(Normalized_counts_Annotated[GROUP],Normalized_counts_Annotated[X_AXIS]),FUN = mean)
                        rownames(Normalized_counts_Annotated_mean)=paste(unlist(Normalized_counts_Annotated_mean[GROUP]),unlist(Normalized_counts_Annotated_mean[X_AXIS]))
                        Normalized_counts_Annotated_mean=Normalized_counts_Annotated_mean[order(rownames(Normalized_counts_Annotated_mean)),]
                        Normalized_counts_Annotated_mean$Row.names=NULL
                        Normalized_counts_Annotated_mean[unique(c(X_AXIS,GROUP))]=NULL
                        Normalized_counts_Annotated_mean[unique(c(X_AXIS,GROUP))]=NULL
                        
                        Heatmap_Normalized_counts_Annotated=Normalized_counts_Annotated
                        rownames(Heatmap_Normalized_counts_Annotated)=make.unique(paste(unlist(Heatmap_Normalized_counts_Annotated[GROUP]),
                                                                                        unlist(Heatmap_Normalized_counts_Annotated[X_AXIS]),
                                                                                        unlist(Heatmap_Normalized_counts_Annotated['Row.names'])),
                                                                                 sep=' ')
                        Heatmap_Normalized_counts_Annotated=Heatmap_Normalized_counts_Annotated[order(rownames(Heatmap_Normalized_counts_Annotated)),]
                        Heatmap_Normalized_counts_Annotated[unique(c(X_AXIS,GROUP))]=NULL
                        Heatmap_Normalized_counts_Annotated[unique(c(X_AXIS,GROUP))]=NULL

                    }else{
                        Normalized_counts_Annotated=merge(colData[unique(c(X_AXIS,GROUP,opt$SPLIT_BY))],t(Normalized_counts_assay),all.x = T  ,by="row.names",sort=F)
                        SPLIT_BY=unlist(as.list(unique(Normalized_counts_Annotated[opt$SPLIT_BY])))
                        Normalized_counts_Annotated_mean=aggregate.data.frame(x = Normalized_counts_Annotated,by = c(Normalized_counts_Annotated[GROUP],Normalized_counts_Annotated[X_AXIS]),FUN = mean)
                        rownames(Normalized_counts_Annotated_mean)=paste(unlist(Normalized_counts_Annotated_mean[GROUP]),unlist(Normalized_counts_Annotated_mean[X_AXIS]))
                        Normalized_counts_Annotated_mean=Normalized_counts_Annotated_mean[order(rownames(Normalized_counts_Annotated_mean)),]
                        SPLIT_BY_col=Normalized_counts_Annotated_mean[opt$SPLIT_BY]
                        Normalized_counts_Annotated_mean$Row.names=NULL
                        Normalized_counts_Annotated_mean[unique(c(X_AXIS,GROUP,opt$SPLIT_BY))]=NULL
                        Normalized_counts_Annotated_mean[unique(c(X_AXIS,GROUP,opt$SPLIT_BY))]=NULL
                        
                        Heatmap_Normalized_counts_Annotated=Normalized_counts_Annotated
                        rownames(Heatmap_Normalized_counts_Annotated)=make.unique(paste(unlist(Heatmap_Normalized_counts_Annotated[GROUP]),
                                                                                        unlist(Heatmap_Normalized_counts_Annotated[X_AXIS]),
                                                                                        unlist(Heatmap_Normalized_counts_Annotated['Row.names'])),
                                                                                 sep=' ')
                        Heatmap_Normalized_counts_Annotated=Heatmap_Normalized_counts_Annotated[order(rownames(Heatmap_Normalized_counts_Annotated)),]
                        Heatmap_Normalized_counts_Annotated[unique(c(X_AXIS,GROUP,opt$SPLIT_BY))]=NULL
                        Heatmap_Normalized_counts_Annotated[unique(c(X_AXIS,GROUP,opt$SPLIT_BY))]=NULL
                    }

                }else{
                    if (is.na(opt$SPLIT_BY)){
                        Normalized_counts_Annotated=merge(colData[c(X_AXIS)],t(Normalized_counts_assay),all.x = T ,by="row.names",sort=F)
                        SPLIT_BY=opt$SPLIT_BY
                        Normalized_counts_Annotated_mean=aggregate.data.frame(x = Normalized_counts_Annotated,by = Normalized_counts_Annotated[X_AXIS],FUN = mean)
                        rownames(Normalized_counts_Annotated_mean)=unlist(Normalized_counts_Annotated_mean[X_AXIS])
                        Normalized_counts_Annotated_mean=Normalized_counts_Annotated_mean[order(rownames(Normalized_counts_Annotated_mean)),]
                        Normalized_counts_Annotated_mean[X_AXIS]=NULL
                        Normalized_counts_Annotated_mean[X_AXIS]=NULL
                        Normalized_counts_Annotated_mean$Row.names=NULL
                        
                        Heatmap_Normalized_counts_Annotated=Normalized_counts_Annotated
                        rownames(Heatmap_Normalized_counts_Annotated)=make.unique(as.character(paste(unlist(Heatmap_Normalized_counts_Annotated[X_AXIS]),
                                                                                                     unlist(Heatmap_Normalized_counts_Annotated['Row.names']))) ,
                                                                                 sep=' ')
                        Heatmap_Normalized_counts_Annotated=Heatmap_Normalized_counts_Annotated[order(rownames(Heatmap_Normalized_counts_Annotated)),]
                        Heatmap_Normalized_counts_Annotated[X_AXIS]=NULL
                        Heatmap_Normalized_counts_Annotated[X_AXIS]=NULL

                      
                    }else{
                        Normalized_counts_Annotated=merge(colData[unique(c(X_AXIS,opt$SPLIT_BY))],t(Normalized_counts_assay),all.x = T ,by="row.names",sort=F)
                        SPLIT_BY=unlist(as.list(unique(Normalized_counts_Annotated[opt$SPLIT_BY])))
                        Normalized_counts_Annotated_mean=aggregate.data.frame(x = Normalized_counts_Annotated,by = Normalized_counts_Annotated[X_AXIS],FUN = mean)
                        rownames(Normalized_counts_Annotated_mean)=unlist(Normalized_counts_Annotated_mean[X_AXIS])
                        Normalized_counts_Annotated_mean=Normalized_counts_Annotated_mean[order(rownames(Normalized_counts_Annotated_mean)),]
                        SPLIT_BY_col=Normalized_counts_Annotated_mean[opt$SPLIT_BY]
                        Normalized_counts_Annotated_mean[unique(c(X_AXIS,opt$SPLIT_BY))]=NULL
                        Normalized_counts_Annotated_mean[unique(c(X_AXIS,opt$SPLIT_BY))]=NULL
                        Normalized_counts_Annotated_mean$Row.names=NULL
                        
                        Heatmap_Normalized_counts_Annotated=Normalized_counts_Annotated
                        rownames(Heatmap_Normalized_counts_Annotated)=make.unique(as.character( paste( unlist(Heatmap_Normalized_counts_Annotated[X_AXIS]),
                                                                                                       unlist(Heatmap_Normalized_counts_Annotated['Row.names']))),
                                                                                 sep=' ')
                        Heatmap_Normalized_counts_Annotated=Heatmap_Normalized_counts_Annotated[order(rownames(Heatmap_Normalized_counts_Annotated)),]
                        Heatmap_Normalized_counts_Annotated[unique(c(X_AXIS,opt$SPLIT_BY))]=NULL
                        Heatmap_Normalized_counts_Annotated[unique(c(X_AXIS,opt$SPLIT_BY))]=NULL
                    }
                    
                    
                }

                Heatmap_Normalized_counts_Annotated=t(Heatmap_Normalized_counts_Annotated)
                
                for (SPLIT in SPLIT_BY){ 
                    if (is.na(SPLIT)){
                        SPLIT_Normalized_counts_Annotated_mean=t(Normalized_counts_Annotated_mean)
                        Heatmap_Normalized_counts_Annotated_mean=SPLIT_Normalized_counts_Annotated_mean
                        SPLIT=''
                        SPLIT_Normalized_counts_Annotated=Normalized_counts_Annotated
                    }else{
                        SPLIT_Normalized_counts_Annotated_mean=Normalized_counts_Annotated_mean
                        SPLIT_Normalized_counts_Annotated_mean=SPLIT_Normalized_counts_Annotated_mean[SPLIT_BY_col==SPLIT, ]
                        SPLIT_Normalized_counts_Annotated_mean=t(SPLIT_Normalized_counts_Annotated_mean)
                        Heatmap_Normalized_counts_Annotated_mean=t(Normalized_counts_Annotated_mean)
                        SPLIT_Normalized_counts_Annotated=Normalized_counts_Annotated[Normalized_counts_Annotated[opt$SPLIT_BY]==SPLIT, ]
                    }   


                    SPLIT_Normalized_counts_Annotated_mean=SPLIT_Normalized_counts_Annotated_mean[sig.genes,]

                    no_var=names(which( apply(SPLIT_Normalized_counts_Annotated_mean,MARGIN = 1,FUN = sd)==0))
                    if (length(no_var)>0){
                      SPLIT_Normalized_counts_Annotated_mean=SPLIT_Normalized_counts_Annotated_mean[!(row.names(SPLIT_Normalized_counts_Annotated_mean)  %in% no_var),]
                    }
                    
                    if (opt$stand){
                      SPLIT_Normalized_counts_Annotated_mean = t(scale(t(SPLIT_Normalized_counts_Annotated_mean)))
                      SPLIT_Normalized_counts_Annotated_mean = apply(SPLIT_Normalized_counts_Annotated_mean,
                                                                     MARGIN = c(1,2),
                                                                     FUN = function(x) round(x,digits = 5))
                    }
                    
                    if (sum(duplicated(SPLIT_Normalized_counts_Annotated_mean)==F)>1){
                        if (nrow(SPLIT_Normalized_counts_Annotated_mean)>2) {
                           
                            
                            if (nrow(SPLIT_Normalized_counts_Annotated_mean)<=opt$k.max) {
                                k.max=nrow(SPLIT_Normalized_counts_Annotated_mean)-1
                            }else{
                                k.max=opt$k.max
                            }
                            if ((opt$FUNcluster=='click')&(!is.na(opt$CLICK_PATH))){
                                clusters = Run_click(mat =as.matrix(SPLIT_Normalized_counts_Annotated_mean),click_path = opt$CLICK_PATH,outDir = opt$outDir,HOMOGENEITY = opt$CLICK_HOMOGENEITY)
                            }else{
                                if (opt$Mclust){
                                    Fit_Mclust <- Mclust( as.matrix(SPLIT_Normalized_counts_Annotated_mean) ,G=1:k.max,modelNames = mclust.options("emModelNames"))
                                    Number_Of_Clusters=Fit_Mclust$G
                                    Clustering <- eclust(as.matrix(SPLIT_Normalized_counts_Annotated_mean),
                                                       stand =F, 
                                                       FUNcluster= opt$FUNcluster,
                                                       hc_metric =opt$hc_metric, 
                                                       graph = FALSE,
                                                       hc_method = opt$hc_method, 
                                                       k=Number_Of_Clusters ) 
                                    clusters=Clustering$cluster
                                }else{

                                    Clustering <- eclust(as.matrix(SPLIT_Normalized_counts_Annotated_mean),stand =F,
                                                       FUNcluster= opt$FUNcluster,
                                                       hc_metric =opt$hc_metric,
                                                       graph = FALSE,
                                                       hc_method = opt$hc_method,
                                                       k.max =  k.max,
                                                       nboot = opt$nboot,
                                                       gap_maxSE = list(method= "Tibs2001SEmax", SE.factor = 1) )
                                    clusters=Clustering$cluster
                                }
                            }
                        }else if (nrow(SPLIT_Normalized_counts_Annotated_mean)==1){
                            clusters=c(1)
                            names(clusters)=rownames((SPLIT_Normalized_counts_Annotated_mean))
                        }else if (nrow(SPLIT_Normalized_counts_Annotated_mean)==2){
                            Clustering <- eclust(as.matrix(SPLIT_Normalized_counts_Annotated_mean),stand = F,
                                                       FUNcluster= 'kmeans',
                                                       graph = FALSE,
                                                       k.max = 2,
                                                       nboot = opt$nboot,
                                                       gap_maxSE = list(method= "Tibs2001SEmax", SE.factor = 1) )
                            clusters=Clustering$cluster
                            sprintf("%s",'Since only 2 significant genes were detected, the clustering was performed using the kmeans method')
                            
                        }else{
                            clusters=c()
                        }
                    
                    }else{
                        clusters=rep(1,nrow(SPLIT_Normalized_counts_Annotated_mean) )
                        names(clusters)=rownames((SPLIT_Normalized_counts_Annotated_mean))
                    }
                    if (length(no_var)>0){
                      if (length(clusters)>0){
                        temp=rep(max(clusters)+1,length(no_var) )
                      }else{
                        temp=rep(1,length(no_var) )
                      }
                      names(temp)=no_var
                      clusters=c(clusters,temp)
                    }

                    if ((opt$hc_metric=='spearman') | ((opt$hc_metric=='pearson') )){
                      hc_metric='correlation'
                    }else{
                      hc_metric=opt$hc_metric
                    }
                    
                    New_clusters=c()
                    for (num in sort(unique(clusters))){
                        res_red=unique(names(clusters))
                        Genes=res_red[res_red  %in% names(clusters[clusters %in% num])]
                        if (length(Genes)>1){
                            if (opt$stand){
                                heat_map=pheatmap(mat = as.matrix(Heatmap_Normalized_counts_Annotated_mean[Genes,]),
                                                cluster_rows = T,show_rownames=F,
                                                clustering_distance_rows = hc_metric ,
                                                clustering_method=opt$hc_method,
                                                cluster_cols = F,silent = T,scale = "row")
                            }else{
                                heat_map=pheatmap(mat = as.matrix(Heatmap_Normalized_counts_Annotated_mean[Genes,]),
                                                cluster_rows = T,show_rownames=F,
                                                clustering_distance_rows = hc_metric ,
                                                clustering_method=opt$hc_method,
                                                cluster_cols = F,silent = T)
                            }
                            New_clusters=c(New_clusters,clusters[heat_map$tree_row$labels[heat_map$tree_row$order]])
                        }else{
                            New_clusters=c(New_clusters,clusters[Genes])
                        }
                    }

                    clusters=New_clusters

                    intgroup=c(opt$GROUP,opt$X_AXIS)
                    intgroup=intgroup[!is.na(intgroup)]
                    Clustering_Plot=plot_clusters(clusters,SPLIT_Normalized_counts_Annotated,c(GROUP,X_AXIS),c(GROUP,X_AXIS,"Normalized counts"),X_AXIS_ORDER=unique(colData[X_AXIS]))
                    ggsave(file.path(opt$outDir,paste('Clusters_',SPLIT,".pdf",sep="")),Clustering_Plot,dpi = 600,  width = 6.99, height =4.99, units = "in")




                    genes=sig.genes
                    for (Heatmap_Type in c('Mean','Raw')){
                    
                        if (Heatmap_Type=='Mean'){
                            mat = as.matrix(Heatmap_Normalized_counts_Annotated_mean[genes,])
                        }else{
                            mat = apply(X =Heatmap_Normalized_counts_Annotated[genes,] ,MARGIN = c(1,2),FUN = as.numeric)
                            
                            New_clusters=c()
                            for (num in sort(unique(clusters))){
                                res_red=unique(names(clusters))
                                Genes=res_red[res_red  %in% names(clusters[clusters %in% num])]
                                if (length(Genes)>1){
                                    if (opt$stand){
                                        heat_map=pheatmap(mat = mat[Genes,],
                                                        cluster_rows = T,show_rownames=F,
                                                        clustering_distance_rows = hc_metric ,
                                                        clustering_method=opt$hc_method,
                                                        cluster_cols = F,silent = T,scale = "row")
                                    }else{
                                        heat_map=pheatmap(mat = mat[Genes,],
                                                        cluster_rows = T,show_rownames=F,
                                                        clustering_distance_rows = hc_metric ,
                                                        clustering_method=opt$hc_method,
                                                        cluster_cols = F,silent = T)
                                    }
                                    New_clusters=c(New_clusters,clusters[heat_map$tree_row$labels[heat_map$tree_row$order]])
                                }
                            }

                            clusters=New_clusters

                            
                        }
                        mat=mat[names(clusters),]
                        
                        write.table(mat,
                                    file = file.path(opt$outDir,paste('Significant_Genes_',Heatmap_Type,opt$NORMALIZATION_TYPE,'Normalized_counts.txt', sep = '_')) ,
                                    quote = F,
                                    sep = "\t")
                        
                        
                        
                        if (X_AXIS!=GROUP){
                            split=sapply(X =colnames(mat),FUN = function(x) unlist(strsplit(x =x,split = " "))[1])
                            #colnames(split)='split'

                            annotation_col=data.frame(row.names =colnames(mat) ,
                                                    X_AXIS=sapply(X =colnames(mat),FUN = function(x) unlist(strsplit(x =x,split = " "))[2] ),
                                                    GROUP=sapply(X =colnames(mat),FUN = function(x) unlist(strsplit(x =x,split = " "))[1]) )
                            colnames(annotation_col)=c(X_AXIS,GROUP)
                            annotation_col[X_AXIS]  <- factor(annotation_col[,X_AXIS] ,levels=intersect(unique(colData[,X_AXIS]),unique(annotation_col[,X_AXIS]) ))
                            annotation_col[GROUP]   <- factor(annotation_col[,GROUP] ,levels=intersect(unique(colData[,GROUP]),unique(annotation_col[,GROUP]) ))
                            new_order = order(annotation_col[,GROUP],annotation_col[,X_AXIS])
                            annotation_col = annotation_col[new_order,,drop=FALSE]
                            mat = mat[,new_order,drop=FALSE]

                            heat_map=pheatmap(mat =mat ,
                                            cluster_rows= F,show_rownames=F,
                                            cluster_cols=F,scale = "row",
                                            annotation_row=data.frame("Clusters"=factor(clusters) ),
                                            #annotation_names_row=F,
                                            border_color=NA,
                                            silent = T,
                                            annotation_col=annotation_col ,
                                            show_colnames = T,
                                            gaps_row = cumsum(aggregate.data.frame(x = as.data.frame(clusters),by =list(clusters),FUN = length )[,"clusters"] ),
                                            gaps_col = c(cumsum(rev(aggregate.data.frame(x = annotation_col[X_AXIS],by =annotation_col[c(X_AXIS,GROUP)],FUN = length))[,1]),
                                                         rep(cumsum(rev(aggregate.data.frame(x = annotation_col[GROUP],by =annotation_col[c(GROUP)],FUN = length))[,1]),2)
                                                         )
                                            
                            )
                        
                        
                        }else{
                            annotation_col=data.frame(row.names =colnames(mat) ,
                                                    X_AXIS=sapply(X =colnames(mat),FUN = function(x) unlist(strsplit(x =x,split = " "))[1] ))
                            colnames(annotation_col)=c(X_AXIS)
                            annotation_col[X_AXIS]  <- factor(annotation_col[,X_AXIS] ,levels=intersect(unique(colData[,X_AXIS]),unique(annotation_col[,X_AXIS]) ))
                            new_order = order(annotation_col[,X_AXIS])
                            annotation_col = annotation_col[new_order,,drop=FALSE]
                            mat = mat[,new_order,drop=FALSE]
                            if (opt$stand){
                                heat_map=pheatmap(mat = mat ,
                                                  cluster_rows= F,show_rownames=F,
                                                  cluster_cols=F,scale = "row",
                                                  silent = T,
                                                  annotation_row=data.frame("Clusters"=factor(clusters) ),
                                                  #annotation_names_row=T,
                                                  annotation_col=annotation_col,
                                                  show_colnames = T,
                                                  gaps_row = cumsum(aggregate.data.frame(x = as.data.frame(clusters),by =list(clusters),FUN = length )[,"clusters"] ),
                                                  gaps_col = cumsum(aggregate.data.frame(x = annotation_col[X_AXIS],by =annotation_col[X_AXIS],FUN = length)[,-1]),
                                                  border_color=NA
                                )
                            }else{
                                heat_map=pheatmap(mat =mat ,
                                                  cluster_rows= F,show_rownames=F,
                                                  cluster_cols=F,
                                                  silent = T,
                                                  annotation_row=data.frame("Clusters"=factor(clusters) ),
                                                  #annotation_names_row=T,
                                                  annotation_col=annotation_col,
                                                  show_colnames = T,
                                                  gaps_row = cumsum(aggregate.data.frame(x = as.data.frame(clusters),by =list(clusters),FUN = length )[,"clusters"] ),
                                                  gaps_col = cumsum(aggregate.data.frame(x = annotation_col[X_AXIS],by =annotation_col[X_AXIS],FUN = length)[,-1]),
                                                  border_color=NA
                                )
                            }
                        }
                        ggsave(file.path(opt$outDir,paste('Clustering_heatmap_',SPLIT,'_',Heatmap_Type,".pdf",sep="")),heat_map$gtable,dpi = 600, width = 15, height = 15, units = "cm",scale = 1.5)

                        mat2=t(scale(t(mat)))
                        mat3=cbind(clusters[rownames(mat2)],mat2)
                        if (!is.na(DESeqDataSet_Results)){
                            mat3=cbind(mat3,DESeqDataSet_Results_redOrdered[rownames(mat3),c('log2FoldChange',"padj")])
                        }
                        if (!is.na(Annotation)) {
                          mat3=merge.data.frame(mat3,Annotation,by="row.names",all.x = T,sort = F)
                          rownames(mat3)=mat3$Row.names
                          mat3$Row.names=NULL
                        }
                        colnames(mat3)[colnames(mat3)=="V1"]="Clusters"
                        colnames(mat3)[colnames(mat3)==""]="Clusters"
                        mat3=mat3[rownames(mat2),]
                        write.csv(x = mat3,
                                  file = file.path(opt$outDir, paste('Clustering_heatmap_',SPLIT,'_',Heatmap_Type,".csv",sep="") ),
                                  quote = TRUE,
                                  row.names = TRUE)
                    }



                    ##################Clustering-End##################################    

                    ##################GO/KEGG STAR####################################
                    # Go Enrichment:
                    if (Annotation_flag == 'OrgDb'){
                        gene2entrez = bitr(names(clusters), fromType=opt$GENE_ID_TYPE, toType=c("ENTREZID"), OrgDb=opt$Annotation_db)
                        gene2entrez = gene2entrez[!duplicated(gene2entrez[opt$GENE_ID_TYPE]),]
                        gene2entrez$CLUSTERS=sapply(X =gene2entrez[opt$GENE_ID_TYPE],FUN = function(x)  unlist(clusters[x]))

                        clusters_for_GO= gene2entrez$CLUSTERS
                        names(clusters_for_GO)=gene2entrez$ENTREZID

                        allRes <- clusters_enricher(opt$outDir,clusters_for_GO,"GO",opt$Annotation_db,Annotation$ENTREZID,"ENTREZID",paste('Clusters_GO_Enrichment_',SPLIT,".csv",sep=""))
                        if (opt$Enriched_terms_overlap){
                          plot_sheard_genes(allRes@compareClusterResult,opt$outDir,'Clusters_GO_Enrichment')
                        }
                        One_cluster=clusters_for_GO
                        One_cluster[]=1
                        allRes <- clusters_enricher(opt$outDir,One_cluster,"GO",opt$Annotation_db,Annotation$ENTREZID,"ENTREZID",paste('One_Cluster_GO_Enrichment_',SPLIT,".csv",sep=""))
                        if (opt$Enriched_terms_overlap){
                          plot_sheard_genes(allRes@compareClusterResult,opt$outDir,'One_Cluster_GO_Enrichment')
                        }
                    }else{
                        if (GO_flag){
                            allRes <- Clusters_Enrichment_Test(opt$outDir,clusters,GO2name_BP,GO2gene_BP,paste('Clusters_GO_Enrichment_',SPLIT,".csv",sep=""),"GO",pAdjustMethod='fdr',pvalueCutoff=0.05,gene2ko=FALSE)
                            if (opt$Enriched_terms_overlap){
                                plot_sheard_genes(allRes@compareClusterResult,opt$outDir,'Clusters_GO_Enrichment')
                            }
                            One_cluster=clusters
                            One_cluster[]=1
                            allRes <- Clusters_Enrichment_Test(opt$outDir,One_cluster,GO2name_BP,GO2gene_BP,paste('One_Cluster_GO_Enrichment_',SPLIT,".csv",sep=""),"GO",pAdjustMethod='fdr',pvalueCutoff=0.05,gene2ko=FALSE)
                            if (opt$Enriched_terms_overlap){
                                plot_sheard_genes(allRes@compareClusterResult,opt$outDir,'One_Cluster_GO_Enrichment')
                            }
                        }
                    }
                    if (KEGG_flag){
                        allRes <- Clusters_Enrichment_Test(opt$outDir,clusters,Pathway2name,Pathway2gene,paste('KEGG_Clusters_Enrichment_Analysis_',SPLIT,".csv",sep=""),"KEGG",pAdjustMethod='fdr',pvalueCutoff=0.05,gene2ko=FALSE)
                        if (opt$Enriched_terms_overlap){
                            plot_sheard_genes(allRes@compareClusterResult,opt$outDir,'KEGG_Clusters_Enrichment_Analysis_')
                        }
                        One_cluster=clusters
                        One_cluster[]=1
                        allRes <- Clusters_Enrichment_Test(opt$outDir,One_cluster,Pathway2name,Pathway2gene,paste('KEGG_One_Clusters_Enrichment_Analysis_',SPLIT,".csv",sep=""),"KEGG",pAdjustMethod='fdr',pvalueCutoff=0.05,gene2ko=FALSE)
                        if (opt$Enriched_terms_overlap){
                            plot_sheard_genes(allRes@compareClusterResult,opt$outDir,'KEGG_One_Clusters_Enrichment_Analysis_')
                        }
                    }
                    if (ORGANISM_KEGG_flag){
                        allRes  <- Clusters_Enrichment_Test(opt$outDir,clusters,ORGANISM_Pathway2name,ORGANISM_Pathway2gene,paste(ORGANISM,'_KEGG_Clusters_Enrichment_Analysis_',SPLIT,".csv",sep=""),"KEGG",pAdjustMethod='fdr',pvalueCutoff=0.05,gene2ko=FALSE)
                        if (opt$Enriched_terms_overlap){
                            plot_sheard_genes(allRes@compareClusterResult,opt$outDir,'KEGG_Clusters_Enrichment_Analysis_')
                        }
                        One_cluster=clusters
                        One_cluster[]=1
                        allRes <- Clusters_Enrichment_Test(opt$outDir,One_cluster,ORGANISM_Pathway2name,ORGANISM_Pathway2gene,paste(ORGANISM,'_KEGG_One_Clusters_Enrichment_Analysis_',SPLIT,".csv",sep=""),"KEGG",pAdjustMethod='fdr',pvalueCutoff=0.05,gene2ko=FALSE)
                        if (opt$Enriched_terms_overlap){
                            plot_sheard_genes(allRes@compareClusterResult,opt$outDir,'KEGG_One_Clusters_Enrichment_Analysis_')
                        }
                    }
                    if (KEGG_KAAS_flag){
                        allRes <- Clusters_Enrichment_Test(opt$outDir,clusters,KASS_Pathway2name,KASS_Pathway2gene,paste('KASS_KEGG_Clusters_Enrichment_Analysis_',SPLIT,".csv",sep=""),"KEGG",pAdjustMethod='fdr',pvalueCutoff=0.05,gene2ko=FALSE)
                        if (opt$Enriched_terms_overlap){
                            plot_sheard_genes(allRes@compareClusterResult,opt$outDir,'KASS_KEGG_Clusters_Enrichment_Analysis_')
                        }
                        One_cluster=clusters
                        One_cluster[]=1
                        allRes <- Clusters_Enrichment_Test(opt$outDir,One_cluster,KASS_Pathway2name,KASS_Pathway2gene,paste('KAAS_KEGG_One_Clusters_Enrichment_Analysis_',SPLIT,".csv",sep=""),"KEGG",pAdjustMethod='fdr',pvalueCutoff=0.05,gene2ko=FALSE)
                        if (opt$Enriched_terms_overlap){
                            plot_sheard_genes(allRes@compareClusterResult,opt$outDir,'KAAS_KEGG_One_Clusters_Enrichment_Analysis_')
                        }
                    }        
                    ##################GO/KEGG-End####################################    
                }
            }else{
                print('No significant genes were detected')
            }
        }else{
            print('The "X_AXIS" parameter must be given for clustering analysis')
        }
        colData   = Original_colData   
        countData = Original_countData 
    }
    if (test2do=='LTR'){ 
       opt$LRT=NA 
    }
    
}




