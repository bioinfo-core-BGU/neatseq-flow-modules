# Check if required packages are installed:
if(!(all(c("optparse") %in% installed.packages()))) {
    if (Sys.getenv("CONDA_PREFIX")!=""){
        install.packages("optparse", dependencies=TRUE,repos = "http://cran.us.r-project.org")
    } else {
        cat("The 'optparse' package is not installed. You must install it for this script to work!")
        }
}
library("optparse")
# Check if required packages are installed:
if(!(all(c("ape") %in% installed.packages()))) {
    if (Sys.getenv("CONDA_PREFIX")!=""){
        install.packages("ape", dependencies=TRUE,repos = "http://cran.us.r-project.org")
    } else {
        cat("The 'ape' package is not installed. You must install it for this script to work!")
        }
}
library("ape")
# Check if required packages are installed:
if(!(all(c("ggtree") %in% installed.packages()))) {
    if (Sys.getenv("CONDA_PREFIX")!=""){
        source("https://bioconductor.org/biocLite.R")
        biocLite("ggtree")
    } else {
        cat("The 'ggtree' package is not installed. You must install it for this script to work!")
    }
} 
library("ggtree")

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option("--tree", type="character", default=NA, 
              help="Path to Tree file", metavar="character"),
  make_option("--layout", type="character", default='rectangular', 
              help="Tree layout [fan or rectangular]", metavar="character"),
  make_option("--Meta_Data", type="character", default=NA, 
              help="Path to tab-delimited Meta Data file with header line.", metavar="character"),
  make_option("--M_Excel", type="character", default=NA, 
              help="If the Meta_Data input is an Excel file indicate the sheet name to use", metavar="character"),
  make_option("--ID_field", type="character", default="Samples",  
              help="Column name in the Meta Data file for IDs found in the tips of the tree", metavar="character"),
  make_option("--output", type="character", default=file.path(getwd(),"TREE.PDF"), 
              help="Path to output file in PDF format", metavar="character"),
  make_option("--cols_to_use", type="character", default = NA, 
              help="Columns in the Meta Data file to use and the order from the center up", metavar="character"),
  make_option("--open.angle", type="numeric", default = 20, 
              help="Tree open angle", metavar="character"),
  make_option("--branch.length", action="store_true",default=FALSE, 
              help="Don't use branch length [cladogram]", metavar="character"),
  make_option("--connect.tip", action="store_true",default=FALSE, 
              help="Connect the tip to its label?", metavar="character"),
  make_option("--pre_spacer", type="numeric", default=0.05, 
              help="Space before the label text [default=0.05] ", metavar="character"),
  make_option("--post_spacer", type="numeric", default=0.01, 
              help="Space after the label text [default=0.01]", metavar="character"),
  make_option("--OTU", type="character",default=NA,
              help="Column name in the Meta Data file to use as OTU annotation", metavar="character"),
  make_option("--labels",action="store_true",default=FALSE,
              help="Use branch length labels", metavar="character"),
  make_option( "--Tip_labels",action="store_true",default=FALSE,
              help="Show tip labels", metavar="character"),
  make_option("--heatmap", type="character", default=NA, 
              help="Path to Data file to generate a heatmap", metavar="character"),
  make_option("--H_Excel", type="character", default=NA, 
              help="If the heatmap input is an Excel file indicate the sheet name to use", metavar="character"),
  make_option( "--heatmap_count_by_sep",type="character", default=NA, 
               help="Count the sep in each cell to generate the heatmap", metavar="character"),
  make_option( "--heatmap_cell_border", type="character", default="white", 
              help="Color of heatmap cell border [default='white']", metavar="character"),
  make_option("--heatmap_lowest_value", type="character", default="white", 
              help="Color of heatmap lowest value [default='white'] ", metavar="character"),
  make_option("--heatmap_highest_value", type="character", default="red", 
              help="Color of heatmap highest value [default='red']", metavar="character"),
  make_option("--cols_to_use_heatmap", type="character", default = NA, 
              help="Columns in the heatmap Data file to use and the order from the center up", metavar="character"),
  make_option("--heatmap_variable", action="store_true",default=FALSE, 
              help="Use only variable columns in the heatmap", metavar="character"),
  make_option("--tree_by_heatmap", action="store_true",default=FALSE, 
              help="Generate tree using Hierarchical Clustering of the heatmap", metavar="character"),
  make_option("--heatmap_HC_dist", type="character", default = 'maximum', 
              help="The heatmap Hierarchical Clustering dist method", metavar="character"),
  make_option("--heatmap_HC_agg", type="character", default = 'average', 
              help="The heatmap Hierarchical Clustering agglomeration method", metavar="character"),
  make_option("--ID_heatmap_field", type="character", default=NA,  
              help="Column name for IDs found in the tips of the tree in the heatmap Data file", metavar="character")
  ); 




opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", 
                                    option_list=option_list,
                                    epilogue="\n\nAuthor:Liron Levin");
opt = optparse::parse_args(opt_parser);


if (is.na(opt$heatmap)==F){
    
    if (is.na(opt$H_Excel)){
      More_Data<- read.csv(opt$heatmap,check.names = F, sep="\t",row.names = opt$ID_heatmap_field)
    }else {
      if(!(all(c("openxlsx") %in% installed.packages()))) {
        if (Sys.getenv("CONDA_PREFIX")!=""){
          install.packages("openxlsx", dependencies=TRUE,repos = "http://cran.us.r-project.org")
        } else {
          cat("The 'openxlsx' package is not installed. You must install it for this script to work!")
        }
        
      }
      library("openxlsx")
      
      wb_heatmap=loadWorkbook(opt$heatmap)
      if (is.na(opt$ID_heatmap_field)){ 
        More_Data=read.xlsx(wb_heatmap, sheet =opt$H_Excel , startRow = 1, colNames = TRUE,
                            rowNames = TRUE, detectDates = FALSE, skipEmptyRows = TRUE,
                            skipEmptyCols = TRUE)
      }else{
        More_Data=read.xlsx(wb_heatmap, sheet =opt$H_Excel , startRow = 1, colNames = TRUE,
                            rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE,
                            skipEmptyCols = TRUE)
        rownames(More_Data)=More_Data[,opt$ID_heatmap_field]
        More_Data=subset.data.frame(More_Data,,opt$ID_heatmap_field !=colnames(More_Data))
        
      }
      
    }
    
    if (is.na(opt$cols_to_use_heatmap)){
      opt$cols_to_use_heatmap=colnames(More_Data)
    }
    
    colnames(More_Data)=unlist(lapply(X =colnames(More_Data),FUN = function(x)  stringr::str_trim(gsub('"','',gsub("'",'',x))))  )
    More_Data=More_Data[,opt$cols_to_use_heatmap]
    if (is.na(opt$heatmap_count_by_sep)==FALSE){
      More_Data=as.data.frame(row.names=rownames(More_Data),x=apply(X = More_Data,MARGIN = c(1,2),FUN =function(x) stringr::str_count(x,opt$heatmap_count_by_sep)+1) )
      More_Data=as.data.frame(row.names=rownames(More_Data),x=apply(X = More_Data,MARGIN = c(1,2),FUN =function(x)  if (is.na(x)) 0 else x) )
    } else if (all(apply(X = More_Data,MARGIN = c(1,2),FUN =function(x) is.numeric(x)))==FALSE){
      temp_More_Data=as.data.frame(row.names=rownames(More_Data),x=apply(X = More_Data,MARGIN = c(1,2),FUN =function(x) as.numeric(x) ))
      if (all(apply(X = temp_More_Data,MARGIN = c(1,2),FUN =function(x) is.na(x)==FALSE))){
        More_Data=temp_More_Data
      }else{
        More_Data=as.data.frame(row.names=rownames(More_Data),x=apply(X = More_Data,MARGIN = c(1,2),FUN =function(x) if (length(x)>0) 1 else 0 ) )
      }
    }
    
    h=hclust(dist(t.data.frame(More_Data),method=opt$heatmap_HC_dist),method=opt$heatmap_HC_agg )
    More_Data=More_Data[,h$order]
    
    rownames(More_Data)=unlist(lapply(X =rownames(More_Data),FUN = function(x)  stringr::str_trim(gsub(' ','',gsub('"','',gsub("'",'',x))))  ))
    
    if (opt$heatmap_variable){
      More_Data=subset.data.frame(x=More_Data,,apply(X = More_Data,MARGIN = 2,FUN =function(x) length(unique(x))>1 ))
    }
    
    tree=as.phylo.hclust( hclust(dist(More_Data, method=opt$heatmap_HC_dist), method=opt$heatmap_HC_agg ))
}



if (opt$tree_by_heatmap==FALSE){
  tree <- read.tree(opt$tree)
}
  



if (is.na(opt$M_Excel)){
    Meta_Data<- read.csv(opt$Meta_Data, sep="\t")
}else {
    if(!(all(c("openxlsx") %in% installed.packages()))) {
        if (Sys.getenv("CONDA_PREFIX")!=""){
            install.packages("openxlsx", dependencies=TRUE,repos = "http://cran.us.r-project.org")
        } else {
            cat("The 'openxlsx' package is not installed. You must install it for this script to work!")
        }
        
    }
    library("openxlsx")
    
    wb=loadWorkbook(opt$Meta_Data)
    Meta_Data=read.xlsx(wb, sheet =opt$M_Excel , startRow = 1, colNames = TRUE,
                        rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE,
                        skipEmptyCols = TRUE)
}
Meta_Data=as.data.frame(apply(X = Meta_Data,MARGIN = 2,FUN =function(x) as.character(gsub(" ", "-", stringr::str_trim(x)) )))
Meta_Data[Meta_Data==""]=NA



tree$tip.label=unlist(lapply(X =tree$tip.label,FUN = function(x)  stringr::str_trim(gsub('"','',gsub("'",'',x)))))
#tree=root(tree,outgroup ="",resolve.root=T,edgelabel = TRUE)
Meta_Data=Meta_Data[c(opt$ID_field, setdiff(names(Meta_Data), opt$ID_field))]
Meta_Data[,opt$ID_field]=sapply(X =Meta_Data[,opt$ID_field],FUN = function(x)   stringr::str_trim(gsub('"','',gsub("'",'',x))) ) 
Meta_Data=subset.data.frame(Meta_Data,Meta_Data[,opt$ID_field] %in%  tree$tip.label,)
colnames(Meta_Data)=make.names(unlist(lapply(X =colnames(Meta_Data),FUN = function(x)  stringr::str_trim(gsub('"','',gsub("'",'',x))))  ))

print('step1')


if (is.na(opt$cols_to_use)){
  opt$cols_to_use=colnames(Meta_Data)
  if (!is.na(opt$OTU)){
    if (opt$OTU %in% colnames(Meta_Data)){
      opt$cols_to_use=c(opt$OTU, setdiff(opt$cols_to_use, opt$OTU))
    }else{
      opt$OTU=NA
    }
  }
  
}else{
  opt$cols_to_use=opt$cols_to_use[opt$cols_to_use %in% colnames(Meta_Data)]
}

if (opt$ID_field  %in%  opt$cols_to_use){
  opt$cols_to_use=opt$cols_to_use[opt$ID_field !=opt$cols_to_use]
} 
names2write=make.names(unlist(lapply(X =stringi::stri_split(opt$cols_to_use,regex ="," ),FUN = function(x)  stringr::str_trim(gsub('"','',gsub("'",'',x))))  ))

print('step2')

OTU=names2write[1]
if (opt$branch.length){
  p=ggtree(tree,ladderize = T, layout=opt$layout,open.angle=opt$open.angle,branch.length="none")  
}else{
  p=ggtree(tree,ladderize = T, layout=opt$layout,open.angle=opt$open.angle)  
}


print('step3')


p <- p %<+% Meta_Data
p$data$j<-0.05

if (opt$layout!='fan'){
  p$data$angle=0
  angle=0
  jump=0.02
}else{
  jump=0.04
  angle=335
}

print('step4')
pre_spacer=opt$pre_spacer
post_spacer=opt$post_spacer 
size=3
p$data$xx[p$data$isTip]=max(p$data$x)
if(opt$connect.tip){
  p$data$x[p$data$isTip]=max(p$data$x)
}
data=p$data
if (opt$Tip_labels){
  p =p+geom_tiplab2(data=p$data,aes_string(x="xx+(j*xx)", label="label",angle="angle"),color="black",size=size,align=F) 
  p$data$j=p$data$j+(size*jump*max(nchar(as.character(p$data$label)),na.rm = T))+post_spacer
}
print('step5')

for (i in names2write){ 
  if (all(is.na(p$data[, i]))==FALSE){
      print(i)
      p$data$color=as.character(as.integer(p$data[,i])^3)
      if (!is.na(opt$OTU)){
        if (i==make.names(opt$OTU)){
          p=groupOTU(p,split(p$data$label,as.character(p$data[,i]))) 
          dic=p$data$color[p$data$isTip]
          names(dic)=p$data$group[p$data$isTip]
          p$data$group=sapply(X = p$data$group,FUN = function(x) ifelse(is.na(dic[as.character(x)]),x,dic[as.character(x)]))
          p=p+aes(color=group)
        }
      }
      
      p$data$str=gsub('\\.',' ',i)#stringi::stri_trans_totitle(i))
      p =p+geom_tippoint(data=p$data,aes_string(x="xx+(j*xx)",color="color",shape="str"),size=size,na.rm=T) 
      p$data$j=p$data$j+size*jump+pre_spacer
      if (opt$layout!='fan'){
        p=p+annotate('text',family = "Helvetica",x=max(p$data$x)+(max(p$data$x)*p$data$j), y =max(p$data$y)+1, size=size*1.2,angle=angle,label=stringi::stri_pad_left(str = p$data$str,use_length = T,width = 2.4*nchar(p$data$str)) ) 
      }else{
        p=p+annotate('text',family = "Helvetica",x=max(p$data$x)+(max(p$data$x)*p$data$j), y =min(p$data$y), size=size*1.2,angle=angle,label=stringi::stri_pad_left(str = p$data$str,use_length = T,width = 2.4*nchar(p$data$str)+15) ) 
      }
      
      p =p+geom_tiplab2(data=p$data,aes_string(x="xx+(j*xx)",color="color", label=i,angle="angle"),size=size ,align=F) 
      if (opt$layout!='fan'){
        p$data$j=p$data$j+(size*jump*max(nchar(c(as.character(p$data[,i]),p$data$str)),na.rm = T))+post_spacer
      }else{
        p$data$j=p$data$j+(size*jump*max(nchar(as.character(p$data[,i])),na.rm = T))+post_spacer
      }
  }
}
print('step7')
if (opt$labels){
  p=p+ geom_label(size=size*0.4,aes(label=as.integer(branch.length)))  
}

p =p+geom_tippoint(data=p$data,aes_string(x="xx+(j*xx)",shape="str"),color='white',size=0,na.rm=T) 
p=p+scale_shape_manual(values = c(11:24))+scale_color_discrete(direction=-1,na.value=NA)+guides(color="none",shape=guide_legend(title="Annotation",title.position = "top"  ))+ theme(legend.position = "bottom",legend.box.margin=margin(0, 30, 0, 0))


if (is.na(opt$heatmap)==F){

    if (opt$layout!='fan'){
      if (length(More_Data)>0){
        p2=gheatmap( p,More_Data,colnames_position='bottom',colnames_offset_y=0.4, width=size*0.5, hjust='left', colnames_angle=-15, font.size=size*0.3,offset = max(p$data$x)*max(p$data$j),color =opt$heatmap_cell_border,low = opt$heatmap_lowest_value,high = opt$heatmap_highest_value)
        ggsave(opt$output,p2,device='pdf',dpi = 600,width=6*length(names2write)+(length(More_Data)*size*0.05),height=0.5*length(tree$tip.label)+10,units='cm',limitsize = FALSE)
      }else{
        ggsave(opt$output,p,device='pdf',dpi = 600,width=6*length(names2write),height=0.5*length(tree$tip.label)+10,units='cm',limitsize = FALSE)
      }
     
    }else{
      if (length(More_Data)>0){
        p2=gheatmap( p,More_Data,colnames_offset_y=0.4, width=size*1, hjust='left', colnames_angle=-50, font.size=size*0.7,offset = max(p$data$x)*max(p$data$j),color =opt$heatmap_cell_border,low = opt$heatmap_lowest_value,high = opt$heatmap_highest_value)
        ggsave(opt$output,p2,device='pdf',dpi = 600,width=7*length(names2write)+(length(More_Data)*size*0.7),height=7*length(names2write)+(length(More_Data)*size*0.7),units='cm',limitsize = FALSE)
      }else{
        ggsave(opt$output,p,device='pdf',dpi = 600,width=7*length(names2write)+(length(More_Data)*size*0.7),height=8*length(names2write)+(length(More_Data)*size*0.7),units='cm',limitsize = FALSE)
      }
      
    }
    
}else{
    if (opt$layout!='fan'){
      ggsave(opt$output,p,device='pdf',dpi = 600,width=6*length(names2write),height=0.5*length(tree$tip.label)+10,units='cm',limitsize = FALSE)
    }else{
      ggsave(opt$output,p,device='pdf',dpi = 600,width=7*length(names2write)+(length(More_Data)*size*0.7),height=8*length(names2write)+(length(More_Data)*size*0.7),units='cm',limitsize = FALSE)
    }
  }



