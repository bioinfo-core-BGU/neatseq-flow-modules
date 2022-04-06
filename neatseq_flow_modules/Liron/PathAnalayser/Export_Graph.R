library('visNetwork')
library("optparse")


option_list = list(
  make_option(c("-G", "--Graph_Object"), type="character", default=NA, 
              help="Path to Graph Object File [PathwayIDGraphs.RData]", metavar="character"),
  make_option(c("-P","--Print_Genomes"), action="store_true",default=FALSE, 
              help="Print All Genomes Names in the Graph Object", metavar="character"),
  make_option(c("-I", "--Genome_I"), type="character", default=NA, 
              help="Name of Genome 1 ", metavar="character"),
  make_option(c("-J", "--Genome_J"), type="character", default=NA, 
              help="Name of Genome 2", metavar="character"),
  make_option(c("-O", "--OUT_DIR"), type="character", default=NA, 
              help="Path to Output Directory ", metavar="character"),
  make_option(c("--width"), type="numeric", default = 1800, 
              help="Width of the ouput Figure, The Default is 1800", metavar="character"),
  make_option(c("--height"), type="numeric", default = 1800, 
              help="height of the ouput Figure, The Default is 1800", metavar="character")
); 


#Get user information
opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", 
                                    option_list=option_list,
                                    epilogue="\n\nAuthor:Liron Levin");
opt = optparse::parse_args(opt_parser);
load(opt$Graph_Object)
if (opt$Print_Genomes){
    print(colnames(Graph_data_frame))
}else{
    if (opt$Genome_I != opt$Genome_J ){
        if (opt$Genome_I %in% row.names(Graph_data_frame)){
            if (opt$Genome_J %in% colnames(Graph_data_frame)){
                visNetwork::visSave(graph = visNetwork::visNetwork(Graph_data_frame[[opt$Genome_I,opt$Genome_J]]$nodes, Graph_data_frame[[opt$Genome_I,opt$Genome_J]]$edges, width = opt$width, height = opt$height),
                                    selfcontained = TRUE,
                                    file = file.path(opt$OUT_DIR,paste(opt$Genome_I,'__',opt$Genome_J,'.html',sep = '',collapse = '')) )
            }
        }
    }
}