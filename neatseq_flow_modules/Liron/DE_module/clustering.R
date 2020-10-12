

library("optparse")
library("ggplot2")
library("pheatmap")
library("mclust")
library("factoextra")
library("cowplot")
library("gridExtra")
library("dplyr")
library("stringr")
library("openxlsx")
library("RColorBrewer")
library("colorspace")

theme_set(theme_half_open())

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-c", "--COUNT_DATA_FILE"), type="character", default=NA, 
              help="Path to Gene Level Count Data Files [comma sep]", metavar="character"),
  make_option(c("-s", "--SAMPLE_DATA_FILE"), type="character", default=NA, 
              help="Path to Samples Information File", metavar="character"),    
  make_option(c("-o","--outDir"), type="character", default=NA, 
              help="path to the output directory", metavar="character"),
  make_option(c("--X_AXIS"), type="character", default=NA,
              help="The Filed In the Sample Data To Use as X Axis ", metavar="character"),
  make_option(c("--Y_AXIS_TITLE"), type="character", default='',
              help="The Title to be used for the Y Axis ", metavar="character"),
  make_option(c("--GROUP"), type="character", default=NA,
              help="The Filed In the Sample Data To Group By [can be two fields separated by ,]", metavar="character"),
  make_option(c("--FUNcluster"), type="character", default='hclust',
              help='A clustering function including [kmeans,pam,clara,fanny,hclust,agnes,diana,click,network]. The default is hclust', metavar="character"),
  make_option(c("--Network_corr_method"), type="character", default='pearson',
              help='The correlatin method to use in network clustering. The default is pearson', metavar="character"),
  make_option(c("--Network_clustring_method"), type="character", default='louvain',
              help='The clustering method to use in network clustering. options are [louvain, greedy, infomap, label_prop, spinglass, walktrap] The default is louvain', metavar="character"),
  make_option(c("--Network_r_cutoff"), type="numeric", default=0.85,
              help="The Correlation R cutoff to use in Network clustering. The default is 0.85", metavar="character"),
  make_option(c("--hc_metric"), type="character", default='pearson',
              help="Hierarchical clustering metric to be used for calculating dissimilarities between observations. The default is pearson", metavar="character"),
  make_option(c("--hc_method"), type="character", default='ward.D2',
              help="Hierarchical clustering agglomeration method to be used. The default is ward.D2 ", metavar="character"),
  make_option(c("--CLICK_PATH"), type="character", default=NA,
              help="The path to CLICK program (Shamir et al. 2000). If your using click cite: CLICK and Expander ( Shamir et al. 2000 and Ulitsky et al. 2010) ", metavar="character"),
  make_option(c("--CLICK_HOMOGENEITY"), type="numeric", default=0.5,
              help="The HOMOGENEITY [0-1] of clusters using CLICK program (Shamir et al. 2000). The default is 0.5 ", metavar="character"),
  make_option(c("--smooth"), type="numeric", default=0.5,
              help="The smoothness level of clusters plots. If set to 0 no smoothness will apply, the default is 0.5 ", metavar="character"),
  make_option(c("--k.max"), type="numeric", default=20,
              help="The maximum number of clusters to consider, must be at least two. The default is 20", metavar="character"),
  make_option(c("--k.fixed"), type="numeric", default=NULL,
              help="The maximum number of clusters to consider, must be at least two. The default is 20", metavar="character"),
  make_option(c("--nboot"), type="numeric", default=10,
              help="Number of Monte Carlo (bootstrap) samples for determining the number of clusters [Not For Mclust]. The default is 10 ", metavar="character"),
  make_option(c("--gap_maxSE_method"), type="character", default='firstSEmax',
              help="eclust parameter for determining the number of clusters [Not For Mclust]. The default is firstSEmax ", metavar="character"),
  make_option(c("--gap_maxSE_SE.factor"), type="numeric", default=1,
              help="eclust parameter for determining the number of clusters [Not For Mclust] higher number is less sensitive. The default is 1 ", metavar="character"),
  make_option(c("--stand"), action="store_true",default=FALSE, 
              help="The Data will be Standardized Before Clustering", metavar="character"), 
  make_option(c("--Mclust"), action="store_true",default=FALSE, 
              help="Use Mclust for determining the number of clusters", metavar="character")
); 


#Get user information
opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", 
                                    option_list=option_list,
                                    epilogue="\n\nAuthor:Liron Levin");
opt = optparse::parse_args(opt_parser);


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

Prepair_Normalized_counts_for_clustering <- function(Normalized_counts_assay,colData,opt){
  if (!is.na(opt$X_AXIS)){    
    X_AXIS=opt$X_AXIS
    if (is.na(opt$GROUP)){
      GROUP=X_AXIS
    }else{
      GROUP=opt$GROUP
    }
    if (X_AXIS!=GROUP){
      Normalized_counts_Annotated                               = merge(colData[c(X_AXIS,GROUP)],t(Normalized_counts_assay),all.x = T ,by="row.names",sort=F)
      Normalized_counts_Annotated_mean                          = aggregate.data.frame(x = Normalized_counts_Annotated,by = c(Normalized_counts_Annotated[GROUP],Normalized_counts_Annotated[X_AXIS]),FUN = mean)
      rownames(Normalized_counts_Annotated_mean)                = paste(unlist(Normalized_counts_Annotated_mean[GROUP]),unlist(Normalized_counts_Annotated_mean[X_AXIS]))
      Normalized_counts_Annotated_mean                          = Normalized_counts_Annotated_mean[order(rownames(Normalized_counts_Annotated_mean)),]
      Normalized_counts_Annotated_mean$Row.names                = NULL
      Normalized_counts_Annotated_mean[unique(c(X_AXIS,GROUP))] = NULL
      Normalized_counts_Annotated_mean[unique(c(X_AXIS,GROUP))] = NULL
      
      Heatmap_Normalized_counts_Annotated                       = Normalized_counts_Annotated
      rownames(Heatmap_Normalized_counts_Annotated)             = make.unique(paste(unlist(Heatmap_Normalized_counts_Annotated[GROUP]),
                                                                      unlist(Heatmap_Normalized_counts_Annotated[X_AXIS]),
                                                                      unlist(Heatmap_Normalized_counts_Annotated['Row.names'])),
                                                                sep=' ')
      
      Heatmap_Normalized_counts_Annotated                          = Heatmap_Normalized_counts_Annotated[order(rownames(Heatmap_Normalized_counts_Annotated)),]
      Heatmap_Normalized_counts_Annotated[unique(c(X_AXIS,GROUP))] = NULL
      Heatmap_Normalized_counts_Annotated[unique(c(X_AXIS,GROUP))] = NULL
      #row_names=Heatmap_Normalized_counts_Annotated$Row.names
      Heatmap_Normalized_counts_Annotated$Row.names=NULL
      #row.names(Heatmap_Normalized_counts_Annotated)=row_names
      
    }else{
      Normalized_counts_Annotated                   = merge(colData[c(X_AXIS)],t(Normalized_counts_assay),all.x = T,by="row.names",sort=F)
      Normalized_counts_Annotated_mean              = aggregate.data.frame(x = Normalized_counts_Annotated,by = Normalized_counts_Annotated[X_AXIS],FUN = mean)
      rownames(Normalized_counts_Annotated_mean)    = unlist(Normalized_counts_Annotated_mean[X_AXIS])
      Normalized_counts_Annotated_mean              = Normalized_counts_Annotated_mean[order(rownames(Normalized_counts_Annotated_mean)),]
      Normalized_counts_Annotated_mean[X_AXIS]      = NULL
      Normalized_counts_Annotated_mean[X_AXIS]      = NULL
      Normalized_counts_Annotated_mean$Row.names    = NULL
      
      Heatmap_Normalized_counts_Annotated           = Normalized_counts_Annotated
      rownames(Heatmap_Normalized_counts_Annotated) = make.unique(as.character(paste(unlist(Heatmap_Normalized_counts_Annotated[X_AXIS]),
                                                                                   unlist(Heatmap_Normalized_counts_Annotated['Row.names']))) ,
                                                                sep=' ')
      Heatmap_Normalized_counts_Annotated           = Heatmap_Normalized_counts_Annotated[order(rownames(Heatmap_Normalized_counts_Annotated)),]
      Heatmap_Normalized_counts_Annotated[X_AXIS]   = NULL
      Heatmap_Normalized_counts_Annotated[X_AXIS]   = NULL
      #row_names=Heatmap_Normalized_counts_Annotated$Row.names
      Heatmap_Normalized_counts_Annotated$Row.names = NULL
      #row.names(Heatmap_Normalized_counts_Annotated)=row_names
    }
    
    if (opt$stand){
      Normalized_counts_Annotated         = merge(colData[c(X_AXIS,GROUP)],scale(t(Normalized_counts_assay)),all.x = T ,by="row.names",sort=F)
      Heatmap_Normalized_counts_Annotated = scale(Heatmap_Normalized_counts_Annotated) 
      Normalized_counts_Annotated_mean    = scale(Normalized_counts_Annotated_mean)
    }
    
    row.names(Normalized_counts_Annotated)          = Normalized_counts_Annotated$Row.names
    Normalized_counts_Annotated$Row.names           = NULL
    Heatmap_Normalized_counts_Annotated             = t(Heatmap_Normalized_counts_Annotated)
    Normalized_counts_Annotated_mean                = t(Normalized_counts_Annotated_mean)
    return(list(Normalized_counts_Annotated_mean,Normalized_counts_Annotated,Heatmap_Normalized_counts_Annotated))
  }
}

sub_cluster  <- function(Normalized_counts_assay,opt){
  
  no_var=names(which( apply(Normalized_counts_assay,MARGIN = 1,FUN = sd)==0))
  if (length(no_var)>0){
    Normalized_counts_assay=Normalized_counts_assay[!(row.names(Normalized_counts_assay)  %in% no_var),]
  }
  if (dim(Normalized_counts_assay)[1]>0){
    New_order= c()
    if ((opt$hc_metric=='spearman') | ((opt$hc_metric=='pearson') )){
      hc_metric='correlation'
    }else{
      hc_metric=opt$hc_metric
    }
    if (opt$stand){
      heat_map=pheatmap(mat = as.matrix(Normalized_counts_assay),
                        cluster_rows = T,show_rownames=F,
                        clustering_distance_rows = hc_metric ,
                        clustering_method=opt$hc_method,
                        cluster_cols = F,silent = T,scale = "row")
    }else{
      heat_map=pheatmap(mat = as.matrix(Normalized_counts_assay),
                        cluster_rows = T,show_rownames=F,
                        clustering_distance_rows = hc_metric ,
                        clustering_method=opt$hc_method,
                        cluster_cols = F,silent = T)
    }
    New_order =heat_map$tree_row$labels[heat_map$tree_row$order]
  }
  New_order = c(New_order,no_var)
  return(New_order)
}

plot_clusters<-function(clusters,vs_ano,color_group=c('Type','Time_int'),titles=c("Type","Time","Normalized counts"),split_by=c(),smooth=0.5,X_AXIS_ORDER=NA,show_lines=F){
  plot_list=list()
  count=1
  for (i in sort(unique(clusters))){
    Df=data.frame()
    genes=names(clusters[clusters==i])
    for (j in genes){
      if (length(split_by)>0){
        temp=subset.data.frame(x =vs_ano,,c(split_by,color_group,j))
        colnames(temp)=c("split_by","Type","Time","EXP")
        temp['Gene'] = j
        Df=rbind(Df,temp)
      }else{
        temp=subset.data.frame(x =vs_ano,,c(color_group,j))
        colnames(temp)=c("Type","Time","EXP")
        temp['Gene'] = j
        Df=rbind(Df,temp)
      }
    }
    if (color_group[1]==color_group[2]){
      Df$Type='Trend'
    }else{
      Df$Type = as.character(Df$Type)
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
                              theme(plot.title = element_text(hjust = 0.5)) + 
                              theme(legend.position  ="none")+
                              xlab(titles[2]) +
                              ylab(titles[3]) +
                              theme(axis.title   = element_text(colour = "black", size = 10) )+
                              theme(legend.title = element_text(colour = "black", size = 10) )+
                              theme(legend.text  = element_text(colour = "black", size = 8) )+
                              theme(axis.text.y  = element_text(colour = "black", size = 6))+
                              theme(axis.text.x  = element_text(colour = "black", size = 6))+
                              theme(plot.margin  = unit(c(0.5,0.2,0.5,0.2),"cm")  )+
                              #scale_color_discrete(name=titles[1])+
                              #scale_linetype_manual(name=titles[1],values=c(1,5))+
                              scale_x_discrete(limits=unique(Df$Time))+
                              theme(legend.key.width =unit(3,"line")) 
                            #scale_y_continuous(breaks=c(seq(0,10,by=2)) )
      )
    }else{
      plot_list[count]=list(ggplot(data=Df, aes(x=Time, y=EXP, group=Type)) + 
                              ggtitle( paste(paste("Cluster ", i  ,sep="") ,length(genes) ,sep="\n")) +
                              theme(plot.title = element_text(hjust = 0.5)) + 
                              theme(legend.position  ="none")+
                              xlab(titles[2]) +
                              ylab(titles[3]) +
                              theme(axis.title   = element_text(colour = "black", size = 10) )+
                              theme(legend.title = element_text(colour = "black", size = 10) )+
                              theme(legend.text  = element_text(colour = "black", size = 8) )+
                              theme(axis.text.y  = element_text(colour = "black", size = 6))+
                              theme(axis.text.x  = element_text(colour = "black", size = 6))+
                              theme(plot.margin  = unit(c(0.5,0.2,0.5,0.2),"cm")  )+
                              #scale_color_discrete(name=titles[1])+
                              #scale_linetype_manual(name=titles[1],values=c(1,5))+
                              scale_x_discrete(limits=unique(Df$Time))+
                              theme(legend.key.width = unit(3,"line")) 
                            #scale_y_continuous(breaks=c(seq(0,10,by=2)) )
      )
    }
    
    if (show_lines){
      plot_list[count]=list(plot_list[[count]]+
                              ggplot2::geom_line(aes(colour = Type,group=Gene),size=0.1,alpha=0.1)+
                              scale_color_hue(l=40, c=35)
      )
    }
    
    if (smooth>0){
      if (dim(Df)[1]>20000){
        plot_list[count]=list(plot_list[[count]]+
                                stat_summary(fun.y=mean,geom="line",aes(linetype=Type,color=Type))+
                                stat_summary(fun.data = "mean_se", geom = "errorbar", width = .3,size = 1,aes(color=Type),show.legend = F)
                                
        )
      }else{
        plot_list[count]=list(plot_list[[count]]+
                                stat_smooth(method="loess", span = smooth, fullrange=F, size=0.5,alpha=1.0 ,aes(linetype=Type,color=Type))
                               
        )
      }
      
      gp2=plot_list[count]
      res<-try(get_legend(gp2[[1]] + theme(legend.position  ="right")+guides(color=guide_legend(title=titles[1]),linetype=guide_legend(title=titles[1]))) ,silent = T)
      if (inherits(res,"try-error")){
        plot_list[count]=list(plot_list[[count]]+
                                stat_summary(fun.y=mean, geom="point", size = .3,aes(color=Type))+
                                stat_summary(fun.data = "mean_se", geom = "errorbar", width = .3,size = 1,aes(color=Type),show.legend = F)
        )
      }
    }else{
      plot_list[count]=list(plot_list[[count]]+
                              stat_summary(fun.y=mean, geom="line", size = 0.3,aes(linetype=Type,color=Type))+
                              stat_summary(fun.data = "mean_se", geom = "errorbar", width = .3,size = 1,aes(color=Type),show.legend = F)
      )
    }
    
    count=count+1  
  }
  
  gp2=plot_list[1]
  if (color_group[1]==color_group[2]){
    ml<-marrangeGrob(plot_list, nrow=2, ncol=3, top=NULL,as.table = F,layout_matrix = matrix(1:6, 2, 3, TRUE))
  }else{
    legend <- get_legend(gp2[[1]] + theme(legend.position  ="right")+guides(color=guide_legend(title=titles[1]),linetype=guide_legend(title=titles[1])))
    ml<-marrangeGrob(plot_list, nrow=2, ncol=3,right=legend, top=NULL,as.table = F,layout_matrix = matrix(1:6, 2, 3, TRUE))
    
  }
  return(ml)
}

cluster_Normalized_counts <- function(SPLIT_Normalized_counts_Annotated_mean,SPLIT_Normalized_counts_Annotated,colData,opt) {
  X_AXIS=opt$X_AXIS
  if (is.na(opt$GROUP)){
    GROUP=X_AXIS
  }else{
    GROUP=opt$GROUP
  }
  if (length(unique(colData[,X_AXIS]))==1){
    if (length(unique(colData[,GROUP]))>1){
      temp_X_AXIS = X_AXIS
      X_AXIS      = GROUP
      GROUP       = temp_X_AXIS
      
    }
  }
  
  
  
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
      
      if (opt$FUNcluster=='network'){
        
        clusters = Network_clustring(SPLIT_Normalized_counts_Annotated_mean,opt)
      
      }else if ((opt$FUNcluster=='click')&(!is.na(opt$CLICK_PATH))){
        clusters = Run_click(mat =as.matrix(SPLIT_Normalized_counts_Annotated_mean),click_path = opt$CLICK_PATH,outDir = opt$outDir,HOMOGENEITY = opt$CLICK_HOMOGENEITY)
      }else{
        if (opt$Mclust){
          Fit_Mclust <- Mclust( as.matrix(SPLIT_Normalized_counts_Annotated_mean) ,G=1:k.max,modelNames = mclust.options("emModelNames"))
          Number_Of_Clusters=Fit_Mclust$G
          Clustering <- eclust(as.matrix(SPLIT_Normalized_counts_Annotated_mean),
                               stand      = F, 
                               FUNcluster = opt$FUNcluster,
                               hc_metric  = opt$hc_metric, 
                               graph      = FALSE,
                               hc_method  = opt$hc_method, 
                               k          = Number_Of_Clusters ) 
          clusters=Clustering$cluster
        }else{
          
          Clustering <- eclust(as.matrix(SPLIT_Normalized_counts_Annotated_mean),stand =F,
                               FUNcluster = opt$FUNcluster,
                               hc_metric  = opt$hc_metric,
                               graph      = FALSE,
                               hc_method  = opt$hc_method,
                               k.max      = k.max,
                               k          = opt$k.fixed, 
                               nboot      = opt$nboot,
                               gap_maxSE  = list(method= opt$gap_maxSE_method, SE.factor = opt$gap_maxSE_SE.factor) )
          
          clusters=Clustering$cluster
        }
      }
    }else if (nrow(SPLIT_Normalized_counts_Annotated_mean)<3){
      clusters=c(1)
      names(clusters)=rownames((SPLIT_Normalized_counts_Annotated_mean))
      sprintf("%s",'Since only 2 significant genes were detected, No clustering was performed')
    }else if (nrow(SPLIT_Normalized_counts_Annotated_mean)==2){
      Clustering <- eclust(as.matrix(SPLIT_Normalized_counts_Annotated_mean),stand = F,
                           FUNcluster= 'kmeans',
                           graph = FALSE,
                           k.max = 2,
                           nboot = opt$nboot,
                           gap_maxSE = list(method= opt$gap_maxSE_method, SE.factor = opt$gap_maxSE_SE.factor) )
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
        heat_map=pheatmap(mat = as.matrix(SPLIT_Normalized_counts_Annotated_mean[Genes,]),
                          cluster_rows = T,show_rownames=F,
                          clustering_distance_rows = hc_metric ,
                          clustering_method=opt$hc_method,
                          cluster_cols = F,silent = T,scale = "row")
      }else{
        heat_map=pheatmap(mat = as.matrix(SPLIT_Normalized_counts_Annotated_mean[Genes,]),
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
  Clustering_Plot=plot_clusters(clusters,SPLIT_Normalized_counts_Annotated,c(GROUP,X_AXIS),c(GROUP,X_AXIS,opt$Y_AXIS_TITLE),X_AXIS_ORDER=unique(colData[X_AXIS]),smooth = opt$smooth)
  ggsave(file.path(opt$outDir,'Clusters.pdf'),Clustering_Plot,dpi = 600,  width = 6.99, height =4.99, units = "in")
  return(clusters)
}

Heatmap_cluster <- function(mat,clusters,Heatmap_Type,opt,genes=c()){
  X_AXIS=opt$X_AXIS
  if (is.na(opt$GROUP)){
    GROUP=X_AXIS
  }else{
    GROUP=opt$GROUP
  }
  if (length(unique(colData[,X_AXIS]))==1){
    if (length(unique(colData[,GROUP]))>1){
      temp_X_AXIS = X_AXIS
      X_AXIS      = GROUP
      GROUP       = temp_X_AXIS
      
    }
   } 
  
  if ((opt$hc_metric=='spearman') | ((opt$hc_metric=='pearson') )){
    hc_metric='correlation'
  }else{
    hc_metric=opt$hc_metric
  }
  
   if (length(genes)==0){
      genes=row.names(mat)
    }
      
    if (Heatmap_Type=='Mean'){
      mat = as.matrix(mat[genes,])
    }else{
      mat = apply(X =mat[genes,] ,MARGIN = c(1,2),FUN = as.numeric)
      
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
      ann_colors = list()
      num_of_X_AXIS = length(unique(colData[,X_AXIS]))
      ann_colors[X_AXIS]=list(colorRampPalette(rev(brewer.pal(n = 11, name ="Set3")))(num_of_X_AXIS))
      names(ann_colors[[X_AXIS]]) <- c(as.character(unique(colData[,X_AXIS])))
      
      num_of_GROUP = length(unique(colData[,GROUP]))
      ann_colors[GROUP]=list(colorRampPalette(rev(brewer.pal(n = 11, name ="Set1")))(num_of_GROUP))
      names(ann_colors[[GROUP]]) <- c(as.character(unique(colData[,GROUP])))
      
      num_of_CLUSTERS = length(unique(clusters))
      ann_colors['Clusters']=list(colorRampPalette(rev(brewer.pal(n = 11, name ="Spectral")))(num_of_CLUSTERS))
      names(ann_colors[['Clusters']]) <- c(as.character(unique(clusters)))
      
      
      heat_map=pheatmap(mat =mat ,
                        annotation_colors = ann_colors,
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
      ann_colors = list()
      num_of_X_AXIS = length(unique(colData[,X_AXIS]))
      ann_colors[X_AXIS]=list(colorRampPalette(rev(brewer.pal(n = 11, name ="Set3")))(num_of_X_AXIS))
      names(ann_colors[[X_AXIS]]) <- c(as.character(unique(colData[,X_AXIS])))
      
      num_of_CLUSTERS = length(unique(clusters))
      ann_colors['Clusters']=list(colorRampPalette(rev(brewer.pal(n = 11, name ="Spectral")))(num_of_CLUSTERS))
      names(ann_colors[['Clusters']]) <- c(as.character(unique(clusters)))
      
      
      
      if (opt$stand){
        heat_map=pheatmap(mat = mat ,
                          annotation_colors = ann_colors,
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
                          annotation_colors = ann_colors,
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
    ggsave(file.path(opt$outDir,paste('Clustering_heatmap_',Heatmap_Type,".pdf",sep="")),heat_map$gtable,dpi = 600, width = 15, height = 15, units = "cm",scale = 1.5)
    
    mat2=t(scale(t(mat)))
    mat3=cbind(clusters[rownames(mat2)],mat2)

    colnames(mat3)[colnames(mat3)=="V1"]="Clusters"
    colnames(mat3)[colnames(mat3)==""]="Clusters"
    mat3=mat3[rownames(mat2),]
    write.csv(x = mat3,
              file = file.path(opt$outDir, paste('Clustering_heatmap_',Heatmap_Type,".csv",sep="") ),
              quote = TRUE,
              row.names = TRUE)
  }
  
Network_clustring <- function(mat,opt){
  corr_mat <- cor(t(mat),method = opt$Network_corr_method)
  corr_mat[corr_mat<opt$Network_r_cutoff]<-0
  network <- igraph::graph_from_adjacency_matrix( corr_mat, weighted=T, mode="undirected", diag=F)
  CL = switch(   
    opt$Network_clustring_method,   
    "louvain"    = igraph::cluster_louvain(network),   
    "greedy"     = igraph::cluster_fast_greedy(network),   
    "infomap"    = igraph::cluster_infomap(network),   
    "label_prop" = igraph::cluster_label_prop(network), 
    "spinglass"  = igraph::cluster_spinglass(network), 
    "walktrap"   = igraph::cluster_walktrap(network) 
  )   
  
  # CL   = igraph::cluster_louvain(network)
  # CL   = igraph::cluster_fast_greedy(network)
  ## CL   = igraph::cluster_leading_eigen(network)
  ## CL   = igraph::cluster_optimal(network)
  ## CL   = igraph::cluster_edge_betweenness(network)
  # CL   = igraph::cluster_infomap(network)
  # CL   = igraph::cluster_label_prop(network)
  ## CL   = igraph::cluster_spinglass(network)
  # CL   = igraph::cluster_walktrap(network)

  
  print(length(CL))
  #sizes_data= data.frame(row.names =names(sizes(CL)),sizes=as.vector(sizes(CL)) )
  #membership_data=data.frame(row.names =names(membership(CL)),membership=as.vector(membership(CL)))
  #Non_singular_clusters=as.integer(row.names(sizes_data)[sizes_data$sizes>1])
  
  Clusters         = as.vector(igraph::membership(CL))
  names(Clusters)  = names(igraph::membership(CL))   
  
  return(Clusters)
  
}
  
countData <- as.matrix(read.csv(opt$COUNT_DATA_FILE,sep="\t",row.names=1))
colnames(countData) = make.names(colnames(countData))
colData <- read.csv(opt$SAMPLE_DATA_FILE, sep="\t", row.names=1)
rownames(colData) = make.names(rownames(colData))

used_samples = intersect(rownames(colData),colnames(countData))
colData_col <- colnames(colData)
countData <- countData[, used_samples]
colData   <- as.data.frame(colData[used_samples ,])
rownames(colData) <- used_samples
colnames(colData) <- colData_col
all(rownames(colData) == colnames(countData))
countData <- countData[, rownames(colData)]

#colData = as.data.frame(apply(colData, MARGIN = c(1,2), FUN = function(x) as.character(str_trim(x))))


matrix_data = Prepair_Normalized_counts_for_clustering(countData,colData,opt)

matrix_data_mean     = matrix_data[[1]]
matrix_data_raw      = matrix_data[[2]]
matrix_data_heatmap  = matrix_data[[3]]


opt$stand = F

clusters = cluster_Normalized_counts(matrix_data_mean,matrix_data_raw,colData,opt )

opt$stand = T

Heatmap_cluster(matrix_data_mean,clusters,'Mean',opt)

Heatmap_cluster(matrix_data_heatmap,clusters,'Raw',opt)
