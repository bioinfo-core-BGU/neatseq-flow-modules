library("optparse")
library("KEGGREST")
library("pheatmap")
library(KEGGgraph)
library(doParallel)
library('igraph')
library('visNetwork')
library('NetPathMiner')
library('pathview')

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("NetPathMiner")
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("pathview")

option_list = list(
  make_option(c("-D", "--GFF_DIR"), type="character", default=NA, 
              help="Path to a Directory of GFF Files ", metavar="character"),
  make_option(c("-P", "--PATTERN"), type="character", default='*.gff', 
              help="Pattern to search for GFF Files ", metavar="character"),
  make_option(c("-M", "--MATRIX"), type="character", default=NA, 
              help="Path to a TAB Delimited KO/Sample Matrix ", metavar="character"),
  make_option(c("-O", "--OUT_DIR"), type="character", default='*.gff', 
              help="Pattern to search for GFF Files ", metavar="character"),
  make_option(c("--GENES_CUTOFF"), type="numeric", default = 1000, 
              help="The Minimal Number of Genes in a GFF File, Files With Less Will be Discarded , The Default is 1000", metavar="character"),
  make_option(c("--Minimal_PATHWAY_SIZE"), type="numeric", default = 40, 
              help="The Minimal Number of Genes Observed in a Sample that Belong to Specific Pathway, Pathways With Less Will be Discarded , The Default is 40", metavar="character"),
  make_option(c("-T", "--TYPE"), type="character", default='metabolic', 
              help="The Type of Analysis:  metabolic or signaling", metavar="character"),
  make_option(c("--UNIPROT_ID_TYPE"), action="store_true",default=FALSE, 
              help="The Gene Homology IDs are of Uniprot", metavar="character"),
  make_option(c("-C","--CPUs"), type="numeric", default = 5, 
              help="The  Number of CPUs to Use , The Default is 5", metavar="character"),
  make_option(c("--PATHWAY_Reference"), type="character", default=NA, 
              help="A KEGG Pathay to use as a reference for sub-setting all other available pathways ", metavar="character"),
  make_option(c("--PATHWAYS_2_USE"), type="character", default=NA,
              help="List of KEGG Pathways to use [a comma separated list of IDs]", metavar="character")
); 


#Get user information
opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", 
                                    option_list=option_list,
                                    epilogue="\n\nAuthor:Liron Levin");
opt = optparse::parse_args(opt_parser);

Create_Graph_metabolic <- function(edges,pathway_data,original_path,overlap_edges){
        library('visNetwork')
        Nodes=unique(as.vector(as.matrix(edges[,c('from','to')])))
        nodes = apply(pathway_data['to'] , MARGIN = c(1),
                         FUN =  function(x)  (as.integer(x[1]) %in% Nodes ) )  
        nodes=pathway_data[nodes,c('to','to_name','to_reactions')]
        colnames(nodes)=c('id','label','shape')
        nodes[nodes$shape==FALSE,'shape'] = 'star'
        nodes[nodes$shape==TRUE,'shape']  = 'square'
        nodes   = nodes[duplicated(nodes$id)==FALSE,]
        nodes_from=edges[,c('from','from_name','from_reactions')]
        colnames(nodes_from)=c('id','label','shape')
        nodes_from[nodes_from$shape,'shape'] = 'square'
        nodes_from[nodes_from$shape==FALSE,'shape'] = 'star'
        nodes_to=edges[,c('to','to_name','to_reactions')]
        colnames(nodes_to)=c('id','label','shape')
        nodes_to[nodes_to$shape,'shape'] = 'square'
        nodes_to[nodes_to$shape==FALSE,'shape'] = 'star'
        nodes   = rbind(nodes_to,nodes_from)
        nodes   = nodes[duplicated(nodes$id)==FALSE,]
        row.names(nodes)= as.character(nodes$id)
        nodes[,'color'] = 'red'

        original = unique(as.vector(as.matrix(original_path[,c('from','to')])))
        nodes[as.character(original),'color'] = 'blue' 

        overlap = unique(as.vector(as.matrix(overlap_edges[,c('from','to')])))
        nodes[as.character(overlap),'color'] = 'purple'

        nodes = nodes[duplicated.data.frame(nodes)==F,]
        network <- visNetwork(nodes, edges, width = "100%")
        return(network)
}

covareg_pathway_metabolic <- function(pathway_data,ko_list1_original,node_ko_occupancy_cutoff=0.7){
  ko_list1 = ko_list1_original
  compound_pathway_data = pathway_data
  compound_pathway_data = compound_pathway_data[!duplicated(compound_pathway_data),]
  EDGES_of_ko_list1=NA
  if (nrow(compound_pathway_data)>0){
    
    FROM1 = sapply(compound_pathway_data$from_name, 
                   function(x)  mean(unlist(stringi::stri_split(x,fixed    = ' ')) %in%  ko_list1 )>node_ko_occupancy_cutoff )
    
    ko_list1_compound = compound_pathway_data[FROM1,'to_name']
    
    TO1 = apply(compound_pathway_data[c('from_name','to_name')] , MARGIN = c(1),
                FUN =  function(x)  (all(unlist(stringi::stri_split(x[1],fixed    = ' ')) %in%  ko_list1_compound )) & (mean(unlist(stringi::stri_split(x[2],fixed    = ' ')) %in%  ko_list1 )> node_ko_occupancy_cutoff )    )
    
    
    
    EDGES_of_ko_list1 = compound_pathway_data[(FROM1)|(TO1), ]
  }
    return(EDGES_of_ko_list1)
  }

covareg_pathway <- function(pathway_data,ko_list1_original){
  EDGES_of_ko_list1=NA
  ko_list1 = ko_list1_original
  compound_pathway_data = na.omit(pathway_data)
  compound_pathway_data = compound_pathway_data[!duplicated(compound_pathway_data),]

  if (nrow(compound_pathway_data)>0){
    
    
    FROM1 = sapply(compound_pathway_data$from_name, 
                   function(x)  any(unlist(stringi::stri_split(x,fixed    = ' ')) %in%  ko_list1 ))
    
    TO1   = sapply(compound_pathway_data$to_name, 
                   function(x)  any(unlist(stringi::stri_split(x,fixed    = ' ')) %in%  ko_list1 ))
    
    
    EDGES_of_ko_list1 = compound_pathway_data[(FROM1)&(TO1), ]
  }
  return(EDGES_of_ko_list1)
}

copaere_pathwats <- function(pathway_data,ko_list1_original,ko_list2_original){
  list1_path =  data.frame()
  ko_list1 = ko_list1_original
  ko_list2 = ko_list2_original
  original_path=c()
  merge_path=c()
  Cooperation_index=NA
  compound_pathway_data = na.omit(pathway_data)
  compound_pathway_data = compound_pathway_data[!duplicated(compound_pathway_data),]
  first_time=T
  increasing=T
  if (nrow(compound_pathway_data)>0){
    
    # Nodes   = unique( c(compound_pathway_data$from,compound_pathway_data$to))
    # overlap = sapply(Nodes, 
    #                function(x)  any(unlist(stringi::stri_split(x,fixed    = ' ')) %in%  ko_list1 ))
    # overlap = sapply(Nodes[overlap], 
    #                function(x)  any(unlist(stringi::stri_split(x,fixed    = ' ')) %in%  ko_list2 ))
    # overlap = sum(overlap)/length(Nodes)
    
    
    FROM1 = sapply(compound_pathway_data$from_name, 
                        function(x)  any(unlist(stringi::stri_split(x,fixed    = ' ')) %in%  ko_list1 ))
    
    TO1   = sapply(compound_pathway_data$to_name, 
                        function(x)  any(unlist(stringi::stri_split(x,fixed    = ' ')) %in%  ko_list1 ))
    
    
    EDGES_of_ko_list1 = compound_pathway_data[(FROM1)&(TO1), ]
      
    
    FROM2 = sapply(compound_pathway_data$from_name, 
                  function(x)  any(unlist(stringi::stri_split(x,fixed    = ' ')) %in%  ko_list2 ))
    
    TO2   = sapply(compound_pathway_data$to_name, 
                  function(x)  any(unlist(stringi::stri_split(x,fixed    = ' ')) %in%  ko_list2 ))
    
    
    EDGES_of_ko_list2 = compound_pathway_data[(FROM2)&(TO2), ]  
    overlap_edges     = rbind(EDGES_of_ko_list1 ,EDGES_of_ko_list2)
    overlap           = sum(duplicated(overlap_edges ))
    overlap_edges     = overlap_edges[duplicated(overlap_edges ),]
    
      while((length(ko_list2)&(increasing))>0){
          
          from_list1 = sapply(compound_pathway_data$from_name, 
                              function(x)  any(unlist(stringi::stri_split(x,fixed    = ' ')) %in%  ko_list1 ))
          
          to_list1   = sapply(compound_pathway_data$to_name, 
                              function(x)  any(unlist(stringi::stri_split(x,fixed    = ' ')) %in%  ko_list1 ))
          
          
          list1_path=compound_pathway_data[(from_list1)&(to_list1), ]
          
          if (first_time){
            original_path       = list1_path
            original_path_edges = nrow(list1_path)
            first_time          = F
          }

          
          
          from_list2 = sapply(compound_pathway_data$from_name, 
                              function(x)  any(unlist(stringi::stri_split(x,fixed    = ' ')) %in%  ko_list2 ))
          
          to_list2   = sapply(compound_pathway_data$to_name, 
                              function(x)  any(unlist(stringi::stri_split(x,fixed    = ' ')) %in%  ko_list2 ))
          
          list1_path_extra_from=compound_pathway_data[(from_list1)&(to_list2), ]
          list1_path_extra_to=compound_pathway_data[(from_list2)&(to_list1), ]
          
          new_list1_path   = rbind(list1_path,list1_path_extra_from,list1_path_extra_to)
          new_list1_path   = new_list1_path[!duplicated(new_list1_path),]
          if ( nrow(new_list1_path)  >nrow(list1_path) ){
             list1_path=new_list1_path
             new_ko_list1 = c(as.vector(unlist(stringi::stri_split(list1_path$from_name,fixed    = ' '))) ,
                              as.vector(unlist(stringi::stri_split(list1_path$to_name,fixed    = ' '))))
             new_ko_list1 = ko_list2[ko_list2 %in% new_ko_list1]
             ko_list2     = ko_list2[!(ko_list2 %in% new_ko_list1)]
             ko_list1     = unique(c(ko_list1 ,new_ko_list1))
          }else{
            increasing=F
          }
    
      }

    merge_path  = list1_path
    edges       = nrow(list1_path)
    new_edges   = edges-original_path_edges
    total_edges = nrow(compound_pathway_data)
    # Cooperation_index = (edges+1)/(original_path_edges+1)
    # Cooperation_index = (new_edges)/total_edges
    Cooperation_index = (new_edges+1)/(overlap+1)
   
    
  }
return(list(original_path,merge_path,Cooperation_index,overlap,new_edges,total_edges,overlap_edges))
}

old_copaere_pathwats_metabolic <- function(pathway_data,ko_list1_original,ko_list2_original,node_ko_occupancy_cutoff=0.7){
  list1_path =  data.frame()
  ko_list1 = ko_list1_original
  ko_list2 = ko_list2_original
  original_path=c()
  merge_path=c()
  Cooperation_index=NA
  compound_pathway_data = pathway_data
  compound_pathway_data = compound_pathway_data[!duplicated(compound_pathway_data),]
  first_time=T
  increasing=T
  if (nrow(compound_pathway_data)>0){
    
    FROM1 = sapply(compound_pathway_data$from_name, 
                   function(x)  mean(unlist(stringi::stri_split(x,fixed    = ' ')) %in%  ko_list1 )>node_ko_occupancy_cutoff )
    
    ko_list1_compound = compound_pathway_data[FROM1,'to_name']
    
    TO1 = apply(compound_pathway_data[c('from_name','to_name')] , MARGIN = c(1),
                FUN =  function(x)  (all(unlist(stringi::stri_split(x[1],fixed    = ' ')) %in%  ko_list1_compound )) & (mean(unlist(stringi::stri_split(x[2],fixed    = ' ')) %in%  ko_list1 )> node_ko_occupancy_cutoff )    )
    
    
    
    EDGES_of_ko_list1 = compound_pathway_data[(FROM1)|(TO1), ]
    
    FROM2 = sapply(compound_pathway_data$from_name, 
                   function(x)  mean(unlist(stringi::stri_split(x,fixed    = ' ')) %in%  ko_list2 )>node_ko_occupancy_cutoff)
    
    ko_list2_compound = compound_pathway_data[FROM2,'to_name']
    
    TO2 = apply(compound_pathway_data[c('from_name','to_name')] , MARGIN = c(1),
                FUN =  function(x)  (all(unlist(stringi::stri_split(x[1],fixed    = ' ')) %in%  ko_list2_compound )) & (mean(unlist(stringi::stri_split(x[2],fixed    = ' ')) %in%  ko_list2 )>node_ko_occupancy_cutoff)    )
    
    
    
    EDGES_of_ko_list2 = compound_pathway_data[(FROM2)|(TO2), ]
    
    overlap_edges     = rbind(EDGES_of_ko_list1 ,EDGES_of_ko_list2)
    overlap           = sum(duplicated(overlap_edges ))
    overlap_edges     = overlap_edges[duplicated(overlap_edges ),]
    
    while((length(ko_list2)&(increasing))>0){
      
      
      from_list1 = sapply(compound_pathway_data$from_name, 
                     function(x)  mean(unlist(stringi::stri_split(x,fixed    = ' ')) %in%  ko_list1 )>node_ko_occupancy_cutoff )
      
      ko_list1_compound = compound_pathway_data[from_list1,'to_name']
      
      to_list1 = apply(compound_pathway_data[c('from_name','to_name')] , MARGIN = c(1),
                  FUN =  function(x)  (all(unlist(stringi::stri_split(x[1],fixed    = ' ')) %in%  ko_list1_compound )) & (mean(unlist(stringi::stri_split(x[2],fixed    = ' ')) %in%  ko_list1 )> node_ko_occupancy_cutoff )    )
      
      list1_path = compound_pathway_data[(from_list1)|(to_list1), ]
      
      
      
      
      from_list2 = sapply(compound_pathway_data$from_name, 
                          function(x)  mean(unlist(stringi::stri_split(x,fixed    = ' ')) %in%  ko_list2 )>node_ko_occupancy_cutoff )
      
      ko_list2_compound = compound_pathway_data[from_list2,'to_name']
      
      
      list1_path_extra_from_list2_compound = apply(compound_pathway_data[c('from_name','to_name')] , MARGIN = c(1),
                                                    FUN =  function(x)  (all(unlist(stringi::stri_split(x[1],fixed    = ' ')) %in%  ko_list2_compound )) & (mean(unlist(stringi::stri_split(x[2],fixed    = ' ')) %in%  ko_list1 )> node_ko_occupancy_cutoff )    )
      
      subset_compound_pathway_data=compound_pathway_data[from_list2,]
      
      list1_path_extra_from_list2_genes   = sapply(subset_compound_pathway_data$to, function(x) x %in% compound_pathway_data[list1_path_extra_from_list2_compound,'from'])
      if (length(list1_path_extra_from_list2_genes)>0){
          list1_path_extra = subset_compound_pathway_data[list1_path_extra_from_list2_genes,]
      }else{
        list1_path_extra=c()
      }
      list1_path_extra_to_list2_genes     = apply(compound_pathway_data[c('from_name','to_name')] , MARGIN = c(1),
                                                   FUN =  function(x)  (all(unlist(stringi::stri_split(x[1],fixed    = ' ')) %in%  ko_list1_compound )) & (mean(unlist(stringi::stri_split(x[2],fixed    = ' ')) %in%  ko_list2 )> node_ko_occupancy_cutoff )    )
      
     
      list1_path_extra =  rbind(list1_path_extra , compound_pathway_data[(list1_path_extra_to_list2_genes)|(list1_path_extra_from_list2_compound),])
      
      
      if (first_time){
        original_path       = list1_path
        original_path_edges = nrow(list1_path)
        first_time          = F
      }
      
      
      
      new_list1_path   = rbind(list1_path,list1_path_extra)
      new_list1_path   = new_list1_path[!duplicated(new_list1_path),]
      if ( nrow(new_list1_path)  > nrow(list1_path) ){
        list1_path=new_list1_path
        new_ko_list1 = c(as.vector(unlist(stringi::stri_split(list1_path$from_name[list1_path$from_reactions],fixed    = ' '))) ,
                         as.vector(unlist(stringi::stri_split(list1_path$to_name[!list1_path$from_reactions],fixed    = ' '))))
        new_ko_list1 = ko_list2[ko_list2 %in% new_ko_list1]
        ko_list2     = ko_list2[!(ko_list2 %in% new_ko_list1)]
        ko_list1     = unique(c(ko_list1 ,new_ko_list1))
      }else{
        increasing=F
      }
      
    }
    merge_path  = list1_path
    edges       = nrow(list1_path)
    new_edges   = edges-original_path_edges
    total_edges = nrow(compound_pathway_data)
    # Cooperation_index = (edges+1)/(original_path_edges+1)
     Cooperation_index = (new_edges)/total_edges
    # Cooperation_index = (new_edges+1)/(overlap+1)
    
    
  }
  return(list(original_path,merge_path,Cooperation_index,overlap,new_edges,total_edges,overlap_edges))
}


copaere_pathwats_metabolic     <- function(pathway_data,ko_list1_original,ko_list2_original,node_ko_occupancy_cutoff=0.7){
  list1_path =  data.frame()
  ko_list1 = ko_list1_original
  ko_list2 = ko_list2_original
  original_path=c()
  merge_path=c()
  Cooperation_index=NA
  compound_pathway_data = pathway_data
  compound_pathway_data = compound_pathway_data[!duplicated(compound_pathway_data),]
  first_time=T
  increasing=T
  if (nrow(compound_pathway_data)>0){
    
    FROM1 = sapply(compound_pathway_data$from_name, 
                   function(x)  mean(unlist(stringi::stri_split(x,fixed    = ' ')) %in%  ko_list1 )>node_ko_occupancy_cutoff )
    
    ko_list1_compound = compound_pathway_data[FROM1,'to_name']
    
    TO1 = apply(compound_pathway_data[c('from_name','to_name')] , MARGIN = c(1),
                FUN =  function(x)  (all(unlist(stringi::stri_split(x[1],fixed    = ' ')) %in%  ko_list1_compound )) & (mean(unlist(stringi::stri_split(x[2],fixed    = ' ')) %in%  ko_list1 )> node_ko_occupancy_cutoff )    )
    
    
    
    EDGES_of_ko_list1 = compound_pathway_data[(FROM1)|(TO1), ]
    
    FROM2 = sapply(compound_pathway_data$from_name, 
                   function(x)  mean(unlist(stringi::stri_split(x,fixed    = ' ')) %in%  ko_list2 )>node_ko_occupancy_cutoff)
    
    ko_list2_compound = compound_pathway_data[FROM2,'to_name']
    
    TO2 = apply(compound_pathway_data[c('from_name','to_name')] , MARGIN = c(1),
                FUN =  function(x)  (all(unlist(stringi::stri_split(x[1],fixed    = ' ')) %in%  ko_list2_compound )) & (mean(unlist(stringi::stri_split(x[2],fixed    = ' ')) %in%  ko_list2 )>node_ko_occupancy_cutoff)    )
    
    
    
    EDGES_of_ko_list2 = compound_pathway_data[(FROM2)|(TO2), ]
    
    overlap_edges       = rbind(EDGES_of_ko_list1 ,EDGES_of_ko_list2)
    overlap             = sum(duplicated(overlap_edges ))
    overlap_edges       = overlap_edges[duplicated(overlap_edges ),]
    original_path       = EDGES_of_ko_list1
    original_path_edges = nrow(EDGES_of_ko_list1)
    
    while((length(ko_list2)>0)&(increasing)){
      
      Nodes = unique(as.vector(as.matrix(EDGES_of_ko_list1[,c('from','to')])))
      
      
      compound = apply(compound_pathway_data[c('from','from_reactions')] , MARGIN = c(1),
                       FUN =  function(x)  (as.integer(x[1]) %in% Nodes )  & (x[2]) )    
      
      
      list1_path = rbind(EDGES_of_ko_list1,compound_pathway_data[compound,] )
      
      ko_list1_compound = compound_pathway_data[compound,'to_name']   
      
      compounds = unique(c(ko_list2_compound,ko_list1_compound))
      
     
      
      to_list2 = apply(compound_pathway_data[c('from','from_name','to_name')] , MARGIN = c(1),
                       FUN =  function(x)  (as.integer(x[1]) %in% Nodes )  & (all(unlist(stringi::stri_split(x[2],fixed    = ' ')) %in%   compounds )) & (mean(unlist(stringi::stri_split(x[3],fixed    = ' ')) %in%  ko_list2 )> node_ko_occupancy_cutoff )    ) 
      list1_path = rbind(list1_path,compound_pathway_data[to_list2, ])
      
      Nodes = unique(as.vector(as.matrix(list1_path[,c('from','to')])))
      
      to_list2 = apply(compound_pathway_data[c('to','from_name','to_name')] , MARGIN = c(1),
                       FUN =  function(x)  (as.integer(x[1]) %in% Nodes )  & (all(unlist(stringi::stri_split(x[2],fixed    = ' ')) %in%  compounds )) & (mean(unlist(stringi::stri_split(x[3],fixed    = ' ')) %in%  ko_list2 )> node_ko_occupancy_cutoff )    ) 
      list1_path = rbind(list1_path,compound_pathway_data[to_list2, ])
      
      
      # 
      # to_list2 = apply(compound_pathway_data[c('from_name','to_name')] , MARGIN = c(1),
      #                  FUN =  function(x)  (all(unlist(stringi::stri_split(x[1],fixed    = ' ')) %in%  ko_list1_compound )) & (mean(unlist(stringi::stri_split(x[2],fixed    = ' ')) %in%  ko_list2 )> node_ko_occupancy_cutoff )    ) 
      # list1_path = rbind(EDGES_of_ko_list1,compound_pathway_data[to_list2, ])
      # 
      
      list1_path = list1_path[duplicated.data.frame(list1_path)==F,]
      
      ko_list1_compound = list1_path[list1_path$from_reactions,'to_name']

      to_list1 = apply(compound_pathway_data[c('from_name','to_name')] , MARGIN = c(1),
                       FUN =  function(x)  (all(unlist(stringi::stri_split(x[1],fixed    = ' ')) %in%  ko_list1_compound )) & (mean(unlist(stringi::stri_split(x[2],fixed    = ' ')) %in%  ko_list1 )> node_ko_occupancy_cutoff )    )
      list1_path = rbind(list1_path,compound_pathway_data[to_list1, ])
      list1_path = list1_path[duplicated.data.frame(list1_path)==F,]
      
      from_list2 = sapply(compound_pathway_data$from_name, 
                          function(x)  mean(unlist(stringi::stri_split(x,fixed    = ' ')) %in%  ko_list2 )>node_ko_occupancy_cutoff )
      
      ko_list2_compound = compound_pathway_data[from_list2,'to_name']
      

      all_nodes=unique(c(list1_path$from,list1_path$to))
      list1_path_extra_from_list2_compound = apply(compound_pathway_data[c('from_name','to')] , MARGIN = c(1),
                                                   FUN =  function(x)  (all(unlist(stringi::stri_split(x[1],fixed    = ' ')) %in%  ko_list2_compound )) & (x[2] %in% all_nodes  )    )
      list1_path_extra_from_list2_compound  = compound_pathway_data[list1_path_extra_from_list2_compound,]
      
      #list1_path_extra_from_list2_genes     = compound_pathway_data[(from_list2) & (compound_pathway_data$to %in% list1_path_extra$from),]
      
      
      new_list1_path   = rbind(list1_path,list1_path_extra_from_list2_compound)#,list1_path_extra_from_list2_genes)
      new_list1_path   = new_list1_path[!duplicated(new_list1_path),]
      if ( nrow(new_list1_path)  > nrow(EDGES_of_ko_list1) ){
        EDGES_of_ko_list1=new_list1_path
      }else{
        increasing=F
      }
      print( nrow(new_list1_path))
    }
    
    if (nrow(list1_path)==0){
        merge_path        = EDGES_of_ko_list1
        edges             = nrow(EDGES_of_ko_list1)
        new_edges         = 0
        total_edges       = nrow(compound_pathway_data)
        Cooperation_index = (new_edges)/total_edges
    }else{
        merge_path  = list1_path
        edges       = nrow(list1_path)
        new_edges   = edges-original_path_edges
        total_edges = nrow(compound_pathway_data)
        # Cooperation_index = (edges+1)/(original_path_edges+1)
        Cooperation_index = (new_edges)/total_edges
        # Cooperation_index = (new_edges+1)/(overlap+1)
    }
    
  }
  return(list(original_path,merge_path,Cooperation_index,overlap,new_edges,total_edges,overlap_edges))
}


read_GFF<- function(FILE,Uniprot_ID_type=F,Genes_cutoff=100,cpus=2){
  library(doParallel)
  run=TRUE
  while (run){
      cl2 <- try(makeCluster(cpus),silent = T)
      if (!inherits(cl2,"try-error")){
        run=FALSE
      }else{
        run=TRUE
      }
  }
  run=TRUE
  while (run){
      test <- try(registerDoParallel(cl2),silent = T)
      if (!inherits(test,"try-error")){
        run=FALSE
      }else{
        run=TRUE
      }
  }
  
  getDoParWorkers()
  print(FILE)
  GFF=read.delim(FILE,sep = '\t',header = F,comment.char = '#')
  if (length(colnames(GFF))>=9){
    GFF=GFF[GFF['V3']=="CDS",]
    if (dim(GFF)[1]>0){
      # GFF=GFF[sapply(GFF[,'V9'], function(x) stringi::stri_count(str = x,fixed = "similar to AA sequence")>0),]
      INDEX <- foreach(x=1:dim(GFF)[1] , .combine=c) %dopar% {
        result = stringi::stri_count(str = GFF[x,'V9'],fixed = "similar to AA sequence")>0
      }
      GFF=GFF[INDEX,]
      
      if (dim(GFF)[1]>0){
 
        # GFF["Kegg"] = sapply(GFF[,'V9'], function(x) unlist(stringi::stri_split(str =  unlist(stringi::stri_split(str = x,regex  = "similar to AA sequence:\\w+:"))[2],fixed = ";"))[1]) 
        Kegg <- foreach(x=1:dim(GFF)[1] , .combine=c) %dopar% {
          result = unlist(stringi::stri_split(str =  unlist(stringi::stri_split(str = GFF[x,'V9'],regex  = "similar to AA sequence:\\w+:"))[2],fixed = ";"))[1]
        }
        GFF["Kegg"] = Kegg
        
        
        # GFF["ID"]   = sapply(GFF[,'V9'], function(x) unlist(stringi::stri_split(str =  unlist(stringi::stri_split(str = x,regex  = "ID="))[2],fixed = ";"))[1]) 
        ID <- foreach(x=1:dim(GFF)[1] , .combine=c) %dopar% {
          result = unlist(stringi::stri_split(str =  unlist(stringi::stri_split(str = GFF[x,'V9'],regex  = "ID="))[2],fixed = ";"))[1]
        }
        GFF["ID"] = ID
        
        GFF=GFF[c('ID',"Kegg")]

        #colnames(GFF)=c('ContigID',"Kegg")
        if (Uniprot_ID_type){
          GFF[,"Kegg"]=paste('up:',GFF[,"Kegg"],sep = '')
          Steps=c('genes','ko')
        }else{
          Steps=c('ko')
        }
        for (Step in Steps){
          Keggs=GFF[,"Kegg"]
          Keggs= Keggs[c(duplicated(Keggs)==F)]
          Total_ko=data.frame()
          if (length(Keggs)>=Genes_cutoff){
            subset_Keggs=split(Keggs, ceiling(seq_along(Keggs )/100))
            REST<- foreach(x=1:length(subset_Keggs) , .combine=rbind) %dopar% {
              Run = TRUE
              while (Run){
                  if (Step=='genes'){
                    REST=try(KEGGREST::keggConv( Step,subset_Keggs[[x]]),silent = T)
                  }else{
                    REST=try(KEGGREST::keggLink( Step,subset_Keggs[[x]]),silent = T)
                  }
                 if (!inherits(REST,"try-error")){
                   Run = FALSE
                 }
              }
              REST=as.matrix(REST)
              
            }
            
            if (dim(REST)[1]>0){
              if (Step=='genes'){
                Total_ko=data.frame('Uniprot'=row.names(REST),'Kegg'=as.vector(REST))
                Total_ko=Total_ko[c(!duplicated(Total_ko$Uniprot)),]
              }else{
                Total_ko=data.frame('Kegg'=row.names(REST),'KO'=as.vector(REST))
              }
            }
            
            
            if (dim(Total_ko)[1]>0){
              if (Step=='genes'){
                GFF["Uniprot"]=GFF["Kegg"]
                GFF["Kegg"]=NULL
                GFF=merge.data.frame(GFF,Total_ko,by = 'Uniprot',sort = F,all.y = T)
                
              }else{
                GFF=merge.data.frame(GFF,Total_ko,by= 'Kegg',sort = F,all.y = T)
              }
            }
          }
        }
        stopCluster(cl2) 
        if (dim(GFF)[2]==3){
            GFF['Bin']=stringi::stri_replace_last(str = basename(FILE) ,replacement = '',fixed = '.gff')
            return(GFF)
        }else{
            return(matrix(nrow = 0,ncol = 3))
        }
        
      }
    }
  }
  
  stopCluster(cl2)  
  return(matrix(nrow = 0,ncol = 3))
}

path                 = opt$GFF_DIR
pattern              = opt$PATTERN
base_output_dir      = opt$OUT_DIR
Genes_cutoff         = opt$GENES_CUTOFF
Type                 = opt$TYPE
Uniprot_ID_type      = opt$UNIPROT_ID_TYPE
cpus                 = opt$CPUs
KOMatrix             = opt$MATRIX
PATHWAYS_2_USE       = opt$PATHWAYS_2_USE
Minimal_PATHWAY_SIZE = opt$Minimal_PATHWAY_SIZE
PATHWAY_Reference    = opt$PATHWAY_Reference

# path                 = NA
# pattern              = '*.gff'
# base_output_dir      = 'C://Users//levinl//Documents//PathAnalayser//New4' 
# Genes_cutoff         = 1000
# Type                 = 'metabolic'
# Uniprot_ID_type      = FALSE
# cpus                 = 5
# KOMatrix             = 'C://Users//levinl//Documents//PathAnalayser//New//Test.Asian_Elephant_2'
# PATHWAYS_2_USE       = 'path:map01200'#NA
# #PATHWAY_Reference    = 'path:map01100'
# PATHWAY_Reference    = NA
# Minimal_PATHWAY_SIZE = 50 

xml_dir = file.path(base_output_dir,'XML') 
dir.create(xml_dir, showWarnings = FALSE)

RUN=TRUE
while (RUN){
    ko2path=try(KEGGREST::keggLink('pathway','ko'),silent = T)
    if (!inherits(ko2path,"try-error")){
        RUN = FALSE
    }
}

ko2path_data_frame = data.frame(KO=names(ko2path),pathway=ko2path)

if (Type=='metabolic'){
    RUN=TRUE
    while (RUN){
        reaction2ko = try(KEGGREST::keggLink('ko','reaction'),silent = T)
        if (!inherits(reaction2ko,"try-error")){
            RUN = FALSE
        }
    }
  
  reaction2ko = data.frame(reaction=names(reaction2ko),ko=as.vector(reaction2ko))
}


    pivot_table=c()
   
    sample_name=''#basename(path)
    output_dir=base_output_dir # file.path(base_output_dir,sample_name) 
    if (file.exists(file.path(output_dir,paste(sample_name,'P_A_MATRIX.RData',sep = '',collapse = '')))){
        load( file.path(output_dir,paste(sample_name,'P_A_MATRIX.RData',sep = '',collapse = '')))
        #pathways=as.vector(na.omit(unique(as.vector(ko2path[row.names(pivot_table)]))))
    }else{
      if ( file.exists(file.path(KOMatrix))) {
        pivot_table = read.delim(file =KOMatrix,sep = '\t' ,row.names = 1)
        row.names(pivot_table) = paste('ko:',row.names(pivot_table),sep = '')
      }else{
          
          file.names <- dir(path, pattern = pattern)
          dir.create(output_dir, showWarnings = FALSE)
          run=TRUE
          while (run){
            cl <- try(makeCluster(round(cpus/2)),silent = T)
            if (!inherits(cl,"try-error")){
              run=FALSE
            }else{
              run=TRUE
            }
          }
          run=TRUE
          while (run){
            test <- try(registerDoParallel(cl),silent = T)
            if (!inherits(test,"try-error")){
              run=FALSE
            }else{
              run=TRUE
            }
          }
          # cl <- makeCluster(round(cpus/2))
          # registerDoParallel(cl)
          getDoParWorkers()
          print(file.names)
          ko_matrix <- foreach(x=1:length(file.names) , .combine=rbind) %dopar% {
                          GFF = read_GFF( FILE = file.path(path,file.names[x]),Uniprot_ID_type=Uniprot_ID_type,Genes_cutoff=Genes_cutoff)
                          if (dim(GFF)[1]>0){
                            result = GFF
                          } else {
                            result = matrix(nrow = 0,ncol = 3)
                          }
                       }
            print(2)
            stopCluster(cl)
          # ko_matrix=data.frame()
          # for (FILE in file.names){
            # GFF=read_GFF(file.path(path,FILE),Uniprot_ID_type=Uniprot_ID_type,Genes_cutoff=Genes_cutoff)
            # print(dim(GFF))
            # if (dim(GFF)[1]>0){
                # ko_matrix=rbind(ko_matrix,GFF)
            # }
          # }
          
          # cl <- makeCluster(round(cpus))
          # registerDoParallel(cl)
          # getDoParWorkers()
          
          write.table(x = ko_matrix,
                      file = file.path(output_dir,'GFF.tab'),
                      sep='\t')
          
          ko_matrix['Count']=1
          pivot_table=reshape2::dcast(ko_matrix,KO ~ Bin,value.var = 'Count' ,fun.aggregate = sum)
          rownames(pivot_table)=pivot_table$KO
          pivot_table$KO=NULL
        }
       if (length(pivot_table)>0){
          
          pivot_table_bin=pivot_table
          pivot_table_bin[pivot_table_bin>0]=1
          pivot_table_bin['CORE']=apply(pivot_table_bin, MARGIN = 1,FUN = mean)
          pivot_table_bin = pivot_table_bin[order(pivot_table_bin$CORE,decreasing = T),]
          pivot_table_bin = pivot_table_bin[pivot_table_bin['CORE']>0,]
          pivot_table_bin$CORE=NULL
          
          
          #pathways = as.vector(na.omit(unique(ko2path_data_frame[ko2path_data_frame$KO %in% row.names(pivot_table_bin),'pathway'])))
          
          #pathways=as.vector(na.omit(unique(as.vector(ko2path[row.names(pivot_table_bin)]))))
          
          pivot_table = pivot_table_bin
          
          mat=apply(pivot_table_bin, MARGIN = 1,FUN = function(x) if (mean(x)==1) {x+1} else {x} )
          pheatmap(mat =mat,cluster_rows = T,cluster_cols = F ,color = c('white','black','blue'),filename = file.path(output_dir,paste(sample_name,'P_A_MATRIX.pdf',sep = '',collapse = '')))
          save(pivot_table,file = file.path(output_dir,paste(sample_name,'P_A_MATRIX.RData',sep = '',collapse = '')))
          
      }
    }
    #pathways = c('path:map01100','path:map01120',pathways)
    if (!is.na(PATHWAYS_2_USE)){
        pathways = unlist(stringi::stri_split(str =  PATHWAYS_2_USE,regex = ','))
    }
    
    # pivot_table_pathways =  row.names(pivot_table) %in% ko2path_data_frame[ko2path_data_frame$pathway == pathway,'KO']
    
    original_pivot_table     = pivot_table
    original_base_output_dir = output_dir
    
    pathway.id_list = c()
    
    pathways = as.vector(na.omit(unique(ko2path_data_frame[ko2path_data_frame$KO %in% row.names(pivot_table),'pathway'])))
    for (pathway in pathways){
      print(pathway)
      pivot_table     = original_pivot_table
      output_dir      = original_base_output_dir
      pivot_table     = pivot_table[row.names(pivot_table) %in% ko2path_data_frame[ko2path_data_frame$pathway == pathway,'KO'],]
      pathway_data    = data.frame()
      if (is.na(PATHWAY_Reference)){
        pathway.id = stringi::stri_replace(str = pathway,regex = 'path:((map)|(ko))',replacement = '')
      }else{
        pathway.id = stringi::stri_replace(str = PATHWAY_Reference,regex = 'path:((map)|(ko))',replacement = '')
      }
      
      print(dim(pivot_table)[1])
      if ((dim(pivot_table)[1]>=Minimal_PATHWAY_SIZE)&((pathway.id %in% pathway.id_list)==F)){
        
            pathway.id_list = c(pathway.id_list,pathway.id)

            xml.file = paste('ko',pathway.id,'.xml',collapse = '', sep = '')
            
            met.file = paste(Type,'_',pathway.id,'.RData',sep = '',collapse = '')
            if (!(met.file %in% dir(xml_dir))){
              if (!(xml.file %in% dir(xml_dir))){
                  while (!(xml.file %in% dir(xml_dir))){
                      pathview_test = try(pathview::download.kegg(pathway.id = pathway.id,
                                                                  species = 'ko',
                                                                  kegg.dir = xml_dir),silent = T)
                  }
                  
              }  
              xml.file = file.path(xml_dir,xml.file)
            
              #pathway_data=KEGGgraph::parseKGML2DataFrame(xml.file,genesOnly=FALSE,reactions=T)
              pathway_igraph=try(NetPathMiner::KGML2igraph(xml.file,
                                               parse.as = Type,
                                               expand.complexes =F),silent = T)
              if (!inherits(pathway_igraph,"try-error")){
                pathway_data = igraph::as_long_data_frame(pathway_igraph)
                #pathway_data = igraph::as_data_frame(pathway_igraph)
                
                if (Type =='signaling'){
                  pathway_data = pathway_data[sapply(pathway_data$attr, function(x) names(x[1])=='miriam.kegg.compound' ),]
                  save(pathway_data,file = file.path(xml_dir,met.file))
                }else{
                  
                  pathway_data_from_reaction = pathway_data[pathway_data$from_reactions,]
                  reversible = sapply(c(1:dim(pathway_data_from_reaction)[1]),FUN = function(x)  pathway_data_from_reaction$from_attr[[x]]$reversible )
                  pathway_data_from_reaction_reversible = pathway_data_from_reaction[reversible,]
                  
                  
                  pathway_data_to_reactions  = pathway_data[pathway_data$to_reactions,]
                  reversible = sapply(c(1:dim(pathway_data_to_reactions)[1]),FUN = function(x)  pathway_data_to_reactions$to_attr[[x]]$reversible )
                  pathway_data_to_reactions_reversible = pathway_data_to_reactions[reversible,]
                  
                  pathway_data = pathway_data[,c('from','to','from_name','to_name','from_reactions','to_reactions')]
                  
                  pathway_data_to_reactions_reversible  = pathway_data_to_reactions_reversible[,c('to','from','to_name','from_name','to_reactions','from_reactions')]
                  pathway_data_from_reaction_reversible = pathway_data_from_reaction_reversible[,c('to','from','to_name','from_name','to_reactions','from_reactions')]
                  
                  pathway_data_reversible = rbind(pathway_data_to_reactions_reversible,pathway_data_from_reaction_reversible)
                  colnames(pathway_data_reversible) = c('from','to','from_name','to_name','from_reactions','to_reactions')
                  pathway_data = rbind(pathway_data,pathway_data_reversible)
                  pathway_data = pathway_data[!duplicated(pathway_data),]
                  
                  
                  pathway_data$from_name = sapply(pathway_data$from_name,
                                                  function(x)  {ifelse( stringi::stri_startswith(str = x,fixed = 'rn:'), 
                                                                        paste(unique(unlist(lapply(unlist(stringi::stri_split(x,fixed    = ' ')),
                                                                                                   function(y)  as.vector(reaction2ko[reaction2ko$reaction==y,'ko'])   ))),
                                                                              collapse = ' '),
                                                                        x)}
                                                  
                  )
                  
                  pathway_data$to_name = sapply(pathway_data$to_name,
                                                function(x)  {ifelse( stringi::stri_startswith(str = x,fixed = 'rn:'), 
                                                                      paste(unique(unlist(lapply(unlist(stringi::stri_split(x,fixed    = ' ')),
                                                                                                 function(y) as.vector(reaction2ko[reaction2ko$reaction==y,'ko']) ))),
                                                                            collapse = ' '),
                                                                      x)}
                                                
                  )
                  
                  
                  pathway_data = pathway_data[!apply(pathway_data[c('from_name','to_name')],MARGIN = c(1),FUN = function(x) (x[1]=="")|(x[2]=="") ),]
                  compounds_created = pathway_data[sapply(pathway_data$from_name, function(x) stringi::stri_startswith(str = x,fixed = 'ko:') ),'to_name']
                  pathway_data      = pathway_data[sapply(pathway_data$from_name, function(x) (stringi::stri_startswith(str = x,fixed = 'ko:'))|(x %in% compounds_created) ),]
                  
                  pathway_data_from_reaction = pathway_data[sapply(pathway_data$from_name, function(x) stringi::stri_startswith(str = x,fixed = 'ko:')),]
                  pathway_data_from_compound = pathway_data[!sapply(pathway_data$from_name, function(x) stringi::stri_startswith(str = x,fixed = 'ko:')),]
                  
                  
                  aggregated = aggregate.data.frame(pathway_data_from_compound[c('to','from_name')],by = list(pathway_data_from_compound$to) ,FUN = function(x) paste(unique(x),collapse = ' '))
                  pathway_data_from_compound$from_name = sapply(pathway_data_from_compound$to, function(x) aggregated[aggregated$to==x,'from_name'])
                  pathway_data = rbind(pathway_data_from_reaction,pathway_data_from_compound)
                  
                  save(pathway_data,file = file.path(xml_dir,met.file))
                  
                }
              }
            }else{
                load(file.path(xml_dir,met.file)  )
            }
            
            if (nrow(pathway_data)>0){ 
      
                # pathway_data = na.omit(pathway_data)
                Cooperation_index_matrix            = matrix(nrow = ncol(pivot_table),ncol = ncol(pivot_table))
                row.names(Cooperation_index_matrix) = colnames(pivot_table)
                colnames(Cooperation_index_matrix)  = colnames(pivot_table)
                
                Cooperation_index_matrix2            = matrix(nrow = ncol(pivot_table),ncol = ncol(pivot_table))
                row.names(Cooperation_index_matrix2) = colnames(pivot_table)
                colnames(Cooperation_index_matrix2)  = colnames(pivot_table)
                
                graph_data_frame          = data.frame()
                edge_data_frame           = data.frame()
                New_Edges_data_frame      = data.frame()
                Overlap_data_frame        = data.frame()
                original_data_frame       = data.frame()
                Merge_Edges_data_frame    = data.frame()
                Overlap_Edges_data_frame  = data.frame()
                Original_Edges_data_frame = data.frame()
                pathway_info =list()
                Bin_covarege=data.frame()
      
                  for (bin in colnames(pivot_table)){
                      ko_list=c()
                      ko_list=row.names(pivot_table[pivot_table[,bin]>=1,])
                      ko_list=unique(ko_list)
                      if (Type =='signaling'){
                        Bin_covarege[pathway,bin]=nrow(covareg_pathway(pathway_data,ko_list))
                      }else{
                        Bin_covarege[pathway,bin]=nrow(covareg_pathway_metabolic(pathway_data,ko_list))
                      }
                  }
                  
                
                
                    run=TRUE
                    while (run){
                        cl <- try(makeCluster(round(cpus)),silent = T)
                        if (!inherits(cl,"try-error")){
                          run=FALSE
                        }else{
                          run=TRUE
                        }
                    }
                    run=TRUE
                    while (run){
                        test <- try(registerDoParallel(cl),silent = T)
                        if (!inherits(test,"try-error")){
                          run=FALSE
                        }else{
                          run=TRUE
                        }
                    }
                  getDoParWorkers()
      
                if (nrow(pathway_data)>0){ 
                  output_dir = file.path(output_dir,stringi::stri_replace(str = pathway,fixed = 'path:',replacement = '')) 
                  dir.create(output_dir, showWarnings = FALSE)    
                  pivot_table_colnames = colnames(pivot_table)
                  if (!file.exists( file.path(output_dir,paste(pathway.id,"_Total_Data",'.RData',sep = '',collapse = '')) )){
                  #Cooperation_data_Total = foreach(i=1:length(pivot_table_colnames) ,.combine=rbind) %do% {
                                              #library(doParallel)
                                          Cooperation_data_Total=c()
                                          for(i in 1:length(pivot_table_colnames) ){
                                              print(pivot_table_colnames[i])
                                              file = file.path(output_dir,paste(pathway.id,pivot_table_colnames[i],'.RData',sep = '',collapse = ''))
                                              if (!file.exists(file)){
                                                  Cooperation_data_line = foreach(j=1:length(pivot_table_colnames) ,.combine=c) %dopar% {
                                                                              ko_list1=row.names(pivot_table[pivot_table[,pivot_table_colnames[i]]>=1,])
                                                                              ko_list2=row.names(pivot_table[pivot_table[,pivot_table_colnames[j]]>=1,])
                                                                              
                                                                              if (Type =='signaling'){
                                                                                Cooperation_data = copaere_pathwats(pathway_data,ko_list1,ko_list2)  
                                                                              }else{
                                                                                Cooperation_data = copaere_pathwats_metabolic(pathway_data,ko_list1,ko_list2,0)
                                                                              }
                                                                              Cooperation_data = list(Cooperation_data)
                                                                          }
                                                  names(Cooperation_data_line)=pivot_table_colnames
                                                  save(Cooperation_data_line,file = file)
                                              }
                                              Cooperation_data_Total[pivot_table_colnames[i]]=file
                                          }
                                          
                      save(Cooperation_data_Total,file = file.path(output_dir,paste(pathway.id,"_Total_Data",'.RData',sep = '',collapse = '')))
                  }else{
                      load(file.path(output_dir,paste(pathway.id,"_Total_Data",'.RData',sep = '',collapse = '')))
                  }
                  
                  
                  # for (i   in pivot_table_colnames){
                    # file_name           = Cooperation_data_Total[ which(pivot_table_colnames==i) ]
                    # print(file_name)
                    # load(file_name)
                    # for (j in pivot_table_colnames){ 
                      # Cooperation_data   = Cooperation_data_line[[which(pivot_table_colnames==j)]]
                      # #original_path      = Cooperation_data[[1]]
                      # #merge_path         = Cooperation_data[[2]]
                      # Cooperation_index  = Cooperation_data[[3]]
                      # overlap            = Cooperation_data[[4]]
                      # new_edges          = Cooperation_data[[5]]
                      # #total_edges        = Cooperation_data[[6]]
                      # #overlap_edges      = Cooperation_data[[7]]
                      
                      # if (i==j){
                        # Cooperation_index_matrix[i,j]      = NA
                        # New_Edges_data_frame[i,j]          = NA
                        # Overlap_data_frame[i,j]            = NA
                        # #original_data_frame[i,j]           = NA
                      # }else{
                        # Cooperation_index_matrix[i,j]      = Cooperation_index
                        # New_Edges_data_frame[i,j]          = new_edges
                        # Overlap_data_frame[i,j]            = overlap
                        # #original_data_frame[i,j]           = nrow(original_path)
                        
                        # if (!is.na(Cooperation_index_matrix[j,i])){
                          # Cooperation_index_matrix2[i,j] = Cooperation_index_matrix[i,j]*Cooperation_index_matrix[j,i]
                        # }
                      # }
                      # # if (dim(merge_path)[1]>0){
                          
                          
                          
                          
                        # # # original_list_A           = original_path[c('to','to_name')]
                        # # # colnames(original_list_A) = c('ID','Name')
                        # # # if (Type =='signaling'){
                          # # # original_list_B           = original_path[,c('from','from_name')]
                        # # # }else{
                          # # # original_list_B           = original_path[original_path$from_reactions,c('from','from_name')]
                        # # # }
                        # # # colnames(original_list_B) = c('ID','Name')
                        # # # original_list             = rbind(original_list_A,original_list_B)
                        # # # original_list             = original_list[!duplicated(original_list),]
                        # # # Name_list                 = original_list$ID
                        # # # original_list             = as.list(original_list$Name)
                        # # # names(original_list)      = Name_list
                        
                        
                        # # # extra_list_A  = merge_path[c('to','to_name')]
                        # # # colnames(extra_list_A) = c('ID','Name')
                        # # # if (Type =='signaling'){
                          # # # extra_list_B  = merge_path[,c('from','from_name')]
                        # # # }else{
                          # # # extra_list_B  = merge_path[merge_path$from_reactions,c('from','from_name')]
                        # # # }
                        # # # colnames(extra_list_B) = c('ID','Name')
                        # # # extra_list    = rbind(extra_list_A,extra_list_B)
                        # # # extra_list    = extra_list[!duplicated(extra_list),]
                        # # # Name_list     = extra_list$ID
                        # # # extra_list    = as.list(extra_list$Name)
                        # # # names(extra_list) = Name_list
                        
                        
                        # # # overlap_list_A  = overlap_edges[c('to','to_name')]
                        # # # colnames(overlap_list_A) = c('ID','Name')
                        # # # if (Type =='signaling'){
                          # # # overlap_list_B  = overlap_edges[,c('from','from_name')]
                        # # # }else{
                          # # # overlap_list_B  = overlap_edges[overlap_edges$from_reactions,c('from','from_name')]
                        # # # }
                        # # # colnames(overlap_list_B) = c('ID','Name')
                        # # # overlap_list    = rbind(overlap_list_A,overlap_list_B)
                        # # # overlap_list    = overlap_list[!duplicated(overlap_list),]
                        # # # Name_list       = overlap_list$ID
                        # # # overlap_list    = as.list(overlap_list$Name)
                        # # # names(overlap_list) = Name_list
                        
                        
                        
                        # # # g=igraph::graph_from_edgelist( as.matrix(apply(X = merge_path[,c('from','to')],
                                                                       # # # MARGIN = c(1,2),
                                                                       # # # FUN = as.character)),
                                                       # # # directed = T)
                        
                        
                        
                        # # # colors = extra_list[V(g)$name]
                        
                        # # # colors[names(extra_list)]    = 'red'
                        
                        # # # colors[names(original_list)] = 'blue'
                        
                        # # # colors[names(overlap_list)]  = 'purple'
                        
                        # # # V(g)$color = as.character(colors[V(g)$name])
                        
                       
                        
                        # # # size = extra_list[V(g)$name]
                        # # # size = lapply(unlist(size), function(x) if ( stringi::stri_startswith(str = x,fixed = 'ko:')) {15} else {4} )
                        
                        
                        # # # V(g)$size = unlist(size)
                        
                        # # # V(g)$name = as.character(extra_list[V(g)$name])
                        
                        # # # save(g,file = file.path(output_dir,paste(pathway.id,"_",i,"_vs_",j,'.RData',sep = '',collapse = ''),fsep = '\\'))
                        
                        # # #graph_data_frame[i,j]              = list(list(g))
                        
                        
                        # # # id<-tkplot(g,vertex.label.dist=1.5)
                        # # # canvas <- tk_canvas(id)
                        # # # tcltk::tkpostscript(canvas, file=".eps")
                        # # # tk_close(id)
                        # # # names(ko2)=stringi::stri_replace(str = extra_list,fixed = 'ko:',replacement = '')
                        # # # names(ko1)=stringi::stri_replace(str = original_list,fixed = 'ko:',replacement = '')
                        
                        # # # tt=pathview::pathview(gene.data = c(ko1,ko2),pathway.id = pathway.id ,species = 'ko',kegg.native = T,map.null = T,split.group = T,expand.node = F)
                        
                        
                        
                        # # Merge_Edges_data_frame[i,j]        = list(list(merge_path))
                      # # }else{
                        # # #graph_data_frame[i,j]              = list(list(NULL))
                        # # Merge_Edges_data_frame[i,j]        = list(list(NULL))
                      # # }
                      
                      # # if (dim(overlap_edges)[1]>0){
                          # # Overlap_Edges_data_frame[i,j]     = list(list(overlap_edges))
                      # # }else{
                          # # Overlap_Edges_data_frame[i,j]     = list(list(NULL))
                      # # }
                    
                      # # if (dim(original_path)[1]>0){
                          # # Original_Edges_data_frame[i,j]     = list(list(original_path))
                      # # }else{
                          # # Original_Edges_data_frame[i,j]     = list(list(NULL))
                      # # }
                    # }
                  # }
                  
                  
                  Graph_dir = file.path(output_dir,'Graphs') 
                  dir.create(Graph_dir, showWarnings = FALSE)
                  
                  Graph_data_frame      = foreach(i=1:length(Cooperation_data_Total) ,.combine=rbind) %dopar% {
                                            I = names(Cooperation_data_Total)[i]
                                            file_name           = Cooperation_data_Total[I]
                                            print(file_name)
                                            load(file_name)
                                            Graph_data = data.frame()
                                            for (j in names(Cooperation_data_line)){ 
                                                Cooperation_data   = Cooperation_data_line[[j]]
                                                original_path      = Cooperation_data[[1]]
                                                merge_path         = Cooperation_data[[2]]
                                                #Cooperation_index  = Cooperation_data[[3]]
                                                #overlap            = Cooperation_data[[4]]
                                                #new_edges          = Cooperation_data[[5]]
                                                #total_edges        = Cooperation_data[[6]]
                                                overlap_edges      = Cooperation_data[[7]]
                                                if (I==j){
                                                      #Cooperation_index_matrix[i,j]      = NA
                                                      Graph_data[I,j]          = NA
                                                      #Overlap_data_frame[i,j]            = NA
                                                      #original_data_frame[i,j]           = NA
                                                }else{
                                                      #Cooperation_index_matrix[i,j]      = Cooperation_index
                                                      graph <- try(list(Create_Graph_metabolic(merge_path,pathway_data,original_path,overlap_edges)),silent = T)
                                                      if (!inherits(graph,"try-error")){
                                                        Graph_data[I,j] =  graph
                                                        # visNetwork::visSave(graph = visNetwork(Graph_data[[I,j]]$nodes, Graph_data[[I,j]]$edges, width = 1800,height = 1800),
                                                                              # selfcontained = TRUE,
                                                                              # file = file.path(Graph_dir,paste(sample_name,pathway.id,'_',I,'__',j,'.html',sep = '',collapse = '')) )
                                                        # visNetwork::visExport(graph = visNetwork(Graph_data[[I,j]]$nodes, Graph_data[[I,j]]$edges, width = 1800,height = 1800),
                                                                              # type = 'pdf',
                                                                              # name = paste(sample_name,pathway.id,'_',I,'__',j,sep = '',collapse = '')  )
                                                        # save(graph,file = file.path(Graph_dir,paste(sample_name,pathway.id,'_',I,'__',j,'.Graphs.RData',sep = '',collapse = '')))
                                                      }else{
                                                        Graph_data[I,j] = NA
                                                      }
                                                      
                                                      #Overlap_data_frame[i,j]            = overlap
                                                      #original_data_frame[i,j]           = nrow(original_path)
                                                      
                                                    # if (!is.na(Cooperation_index_matrix[j,i])){
                                                        # Cooperation_index_matrix2[i,j] = Cooperation_index_matrix[i,j]*Cooperation_index_matrix[j,i]
                                                    # }
                                                }
                                            }
                                            Graph_data
                                        }
                  save(Graph_data_frame,file = file.path(Graph_dir,paste(sample_name,pathway.id,'Graphs.RData'                  ,sep = '',collapse = '')))
                  
                  New_Edges_data_frame = foreach(i=1:length(Cooperation_data_Total) ,.combine=rbind) %dopar% {
                                            I                     = names(Cooperation_data_Total)[i]
                                            file_name             = Cooperation_data_Total[I]
                                            print(file_name)
                                            load(file_name)
                                            New_Edges_data = data.frame()
                                            for (j in names(Cooperation_data_line)){ 
                                                Cooperation_data   = Cooperation_data_line[[j]]
                                                #original_path      = Cooperation_data[[1]]
                                                #merge_path         = Cooperation_data[[2]]
                                                #Cooperation_index  = Cooperation_data[[3]]
                                                #overlap            = Cooperation_data[[4]]
                                                new_edges          = Cooperation_data[[5]]
                                                #total_edges        = Cooperation_data[[6]]
                                                #overlap_edges      = Cooperation_data[[7]]
                                                if (I==j){
                                                      #Cooperation_index_matrix[i,j]      = NA
                                                      New_Edges_data[I,j]          = NA
                                                      #Overlap_data_frame[i,j]            = NA
                                                      #original_data_frame[i,j]           = NA
                                                }else{
                                                      #Cooperation_index_matrix[i,j]      = Cooperation_index
                                                      New_Edges_data[I,j]          = new_edges
                                                      #Overlap_data_frame[i,j]            = overlap
                                                      #original_data_frame[i,j]           = nrow(original_path)
                                                      
                                                    # if (!is.na(Cooperation_index_matrix[j,i])){
                                                        # Cooperation_index_matrix2[i,j] = Cooperation_index_matrix[i,j]*Cooperation_index_matrix[j,i]
                                                    # }
                                                }
                                            }
                                            New_Edges_data
                                        }
                  write.table(x = New_Edges_data_frame,
                            file = file.path(output_dir,paste(sample_name,pathway.id,'_New_Edges.tab',sep = '',collapse = '')),
                            sep = '\t')
                            
                  try(pheatmap(mat =New_Edges_data_frame ,# fontsize = 1,fontsize_number = 2,
                           cluster_rows = F,
                           cluster_cols = F,
                           display_numbers = F,
                           main = pathway,
                           #breaks = mybreaks,
                           #color = color,
                           silent = T,
                           filename = file.path(output_dir,paste(Type,'_',pathway.id,'_heatmap_New_Edges.pdf',sep = '',collapse = ''))
                  ),silent = T   )       
                            
                  Overlap_data_frame = foreach(i=1:length(Cooperation_data_Total) ,.combine=rbind) %dopar% {
                                          I                     = names(Cooperation_data_Total)[i]
                                          file_name             = Cooperation_data_Total[I]
                                          print(file_name)
                                          load(file_name)
                                          New_data = data.frame()
                                          for (j in names(Cooperation_data_line)){ 
                                              Cooperation_data   = Cooperation_data_line[[j]]
                                              #original_path      = Cooperation_data[[1]]
                                              #merge_path         = Cooperation_data[[2]]
                                              #Cooperation_index  = Cooperation_data[[3]]
                                              overlap            = Cooperation_data[[4]]
                                              #new_edges          = Cooperation_data[[5]]
                                              #total_edges        = Cooperation_data[[6]]
                                              #overlap_edges      = Cooperation_data[[7]]
                                              if (I==j){
                                                    #Cooperation_index_matrix[i,j]      = NA
                                                    New_data[I,j]          = NA
                                                    #Overlap_data_frame[i,j]            = NA
                                                    #original_data_frame[i,j]           = NA
                                              }else{
                                                    #Cooperation_index_matrix[i,j]      = Cooperation_index
                                                    New_data[I,j]          = overlap
                                                    #Overlap_data_frame[i,j]            = overlap
                                                    #original_data_frame[i,j]           = nrow(original_path)
                                                    
                                                  # if (!is.na(Cooperation_index_matrix[j,i])){
                                                      # Cooperation_index_matrix2[i,j] = Cooperation_index_matrix[i,j]*Cooperation_index_matrix[j,i]
                                                  # }
                                              }
                                          }
                                          New_data
                                      }
                  write.table(x = Overlap_data_frame,
                            file = file.path(output_dir,paste(sample_name,pathway.id,'_Overlap.tab',sep = '',collapse = '')),
                            sep = '\t')      
                  
                  try(pheatmap(mat =Overlap_data_frame ,# fontsize = 1,fontsize_number = 2,
                           cluster_rows = F,
                           cluster_cols = F,
                           display_numbers = F,
                           main = pathway,
                           #breaks = mybreaks,
                           #color = color,
                           silent = T,
                           filename = file.path(output_dir,paste(Type,'_',pathway.id,'_heatmap_Overlap.pdf',sep = '',collapse = ''))
                      ) ,silent = T)
                  
                  Original_data_frame = foreach(i=1:length(Cooperation_data_Total) ,.combine=rbind) %dopar% {
                                          I = names(Cooperation_data_Total)[i]
                                          file_name           = Cooperation_data_Total[I]
                                          print(file_name)
                                          load(file_name)
                                          New_data = data.frame()
                                          for (j in names(Cooperation_data_line)){ 
                                              Cooperation_data   = Cooperation_data_line[[j]]
                                              original_path      = Cooperation_data[[1]]
                                              #merge_path         = Cooperation_data[[2]]
                                              #Cooperation_index  = Cooperation_data[[3]]
                                              #overlap            = Cooperation_data[[4]]
                                              #new_edges          = Cooperation_data[[5]]
                                              #total_edges        = Cooperation_data[[6]]
                                              #overlap_edges      = Cooperation_data[[7]]
                                              if (I==j){
                                                    #Cooperation_index_matrix[i,j]      = NA
                                                    New_data[I,j]          = NA
                                                    #Overlap_data_frame[i,j]            = NA
                                                    #original_data_frame[i,j]           = NA
                                              }else{
                                                    #Cooperation_index_matrix[i,j]      = Cooperation_index
                                                    New_data[I,j]          = dim(original_path)[1]
                                                    #Overlap_data_frame[i,j]            = overlap
                                                    #original_data_frame[i,j]           = nrow(original_path)
                                                    
                                                  # if (!is.na(Cooperation_index_matrix[j,i])){
                                                      # Cooperation_index_matrix2[i,j] = Cooperation_index_matrix[i,j]*Cooperation_index_matrix[j,i]
                                                  # }
                                              }
                                          }
                                          New_data
                                      }
                  write.table(x = t(Original_data_frame),
                            file = file.path(output_dir,paste(sample_name,pathway.id,'_Original.tab',sep = '',collapse = '')),
                            sep = '\t')          
      
                  try(pheatmap(mat =t(Original_data_frame) ,# fontsize = 1,fontsize_number = 2,
                           cluster_rows = F,
                           cluster_cols = F,
                           display_numbers = F,
                           main = pathway,
                           #breaks = mybreaks,
                           #color = color,
                           silent = T,
                           filename = file.path(output_dir,paste(Type,'_',pathway.id,'_heatmap_Original.pdf',sep = '',collapse = ''))
                     ),silent = T)          
                  # write.table(x = Overlap_data_frame,
                            # file = file.path(output_dir,paste(sample_name,pathway.id,'_Overlap_Edges.tab',sep = '',collapse = '')),
                            # sep = '\t')
                  # write.table(x = Bin_covarege,
                            # file = file.path(output_dir,paste(sample_name,pathway.id,'Bin_covarege.tab',sep = '',collapse = '')),
                            # sep = '\t')
                  # pathway_info['Graphs']                 = list(graph_data_frame)
                  # pathway_info['Original_Edges_data']    = list(Original_Edges_data_frame)
                  # pathway_info['Merge_Edges_Data']       = list(Merge_Edges_data_frame)
                  # pathway_info['Overlap_Edges_Data']     = list(Overlap_Edges_data_frame)
                  # pathway_info['New_Edges']              = list(New_Edges_data_frame)
                  # pathway_info['Overlap']                = list(Overlap_data_frame)
                  # pathway_info['Cooperation_index_Full'] = list(Cooperation_index_matrix)
                  # pathway_info['Cooperation_index']      = list(Cooperation_index_matrix2)
                  # pathway_info['original']               = list(original_data_frame)
                  # pathway_info['Bin_covarege']           = list(Bin_covarege)
                  # pathway_info['Total_Edges']            = total_edges
                  # save(pathway_info,file = file.path(output_dir,paste(sample_name,pathway.id,'.RData',sep = '',collapse = '')))
                  #total_pathway_info[pathway] = list(pathway_info)
                  
                  # save(graph_data_frame            ,file = file.path(output_dir,paste(sample_name,pathway.id,'Graphs.RData'                  ,sep = '',collapse = '')))
                  
                  # save(Original_Edges_data_frame   ,file = file.path(output_dir,paste(sample_name,pathway.id,'Original_Edges_data.RData'     ,sep = '',collapse = '')))
                  
                  # save(Merge_Edges_data_frame      ,file = file.path(output_dir,paste(sample_name,pathway.id,'Merge_Edges_Data.RData'        ,sep = '',collapse = '')))
                  
                  # save(Overlap_Edges_data_frame    ,file = file.path(output_dir,paste(sample_name,pathway.id,'Overlap_Edges_Data.RData'      ,sep = '',collapse = '')))
                  
                  # save(New_Edges_data_frame        ,file = file.path(output_dir,paste(sample_name,pathway.id,'New_Edges.RData'               ,sep = '',collapse = '')))
                  
                  # save(Overlap_data_frame          ,file = file.path(output_dir,paste(sample_name,pathway.id,'Overlap.RData'                 ,sep = '',collapse = '')))
                  
                  # save(Cooperation_index_matrix    ,file = file.path(output_dir,paste(sample_name,pathway.id,'Cooperation_index_Full.RData'  ,sep = '',collapse = '')))
                  
                  # save(Cooperation_index_matrix2   ,file = file.path(output_dir,paste(sample_name,pathway.id,'Cooperation_index.RData'       ,sep = '',collapse = '')))
                  
                  # save(total_edges                 ,file = file.path(output_dir,paste(sample_name,pathway.id,'Total_Edges.RData'             ,sep = '',collapse = '')))
                  
                  # save(original_data_frame         ,file = file.path(output_dir,paste(sample_name,pathway.id,'Original.RData'                ,sep = '',collapse = '')))
                  
                  # save(Bin_covarege                ,file = file.path(output_dir,paste(sample_name,pathway.id,'Bin_covarege.RData'            ,sep = '',collapse = '')))
                  
                  
                  
                  if (length(na.omit(unique(as.vector(Cooperation_index_matrix))))>1){
                    #color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu") ))( 100)
                    #mybreaks=seq(0,1,0.01)
                    pheatmap(mat =Cooperation_index_matrix ,fontsize = 1,fontsize_number = 2,
                             cluster_rows = F,
                             cluster_cols = F,
                             display_numbers = F,
                             main = pathway,
                             #breaks = mybreaks,
                             #color = color,
                             silent = T,
                             filename = file.path(output_dir,paste(Type,'_',pathway.id,'_heatmap_full.pdf',sep = '',collapse = ''))
                              )
                  }
                  if (length(na.omit(unique(as.vector(Cooperation_index_matrix2))))>1){
                    pheatmap(mat =Cooperation_index_matrix2 ,fontsize = 1 ,fontsize_number = 2,
                             cluster_rows = F,
                             cluster_cols = F,
                             display_numbers = F,
                             main = pathway,
                             
                             #breaks = mybreaks,
                             #color = color,
                             silent = T,
                             filename = file.path(output_dir,paste(Type,'_',pathway.id,'_heatmap.pdf',sep = '',collapse = ''))
                    )
                  }
                }  
             }
     }
    }
    
    stopCluster(cl)
    #save(total_pathway_info,file = file.path(output_dir,paste(Type,'_',sample_name,'.RData',sep = '',collapse = ''),fsep = '\\'))
    