

library(rvest)
library(stringr)
library(tidyverse)
library(yaml)



# Get Plugin-Method list:

address <- str_interp("https://docs.qiime2.org/2019.1/plugins/available/")
feattable <- read_html(address)
# html_structure(feattable)

max_meths = feattable %>%
    html_nodes("a.reference.internal") %>%
    html_text() %>% length
method_index <- data.frame(plugin=character(length = max_meths),
                           method=character(length = max_meths),
                           type=character(length = max_meths),stringsAsFactors = F)
i<-0
for(methodl in feattable %>%
    html_nodes("a.reference.internal") %>%
    html_text()) {
    # If the text contains "Getting started", reached the end of the list. Break from loop
    if(str_detect(string = methodl,pattern = "Getting started")) {
        break
    }
    # If the text contains "Plugin", extract the plugin name:
    if(str_detect(string = methodl,pattern = "Plugin")) {
        # print(methodl)
        plugin = str_extract(string = methodl,pattern = "^([^\\:]+)")
        # print(plugin)
        next
        # print("in 1")
    }
    # If the text is one of Visualizers, Methods, Pipelines, save as type
    if(str_detect(string = methodl,pattern = "Visualizers")|
       str_detect(string = methodl,pattern = "Methods")|
       str_detect(string = methodl,pattern = "Pipelines")) {
        typel=methodl
        next
        # print("in 2")
    }
    # Else, extract method name and store
    i <- i+1
    methodname = str_extract(string = methodl,pattern = "^([^\\:]+)")
    cat(sprintf("%s\t%s\n",plugin,methodname))
    method_index[i,] <- c(plugin,methodname,typel)

} 
method_index <- method_index[1:i,]

method_index %>% View

plugin_index <- method_index

# Extracting argument info from plugin method page

# plugin <- "feature-classifier"
# p_method <- "fit-classifier-naive-bayes"
# plugin <- "fragment-insertion"
# p_method <- "sepp"
# plugin <- "diversity"
# p_method <- "core-metrics-phylogenetic"
# 
# plugin_index <- method_index[method_index$plugin %in% c("vsearch","demux", "dada2","longitudinal","emperor","diversity") & 
#                                  method_index$method %in% c("dereplicate-sequences", "emp-single", "denoise-single","first-differences","biplot","beta-rarefaction"),]
# 
# plugin <- plugin_index$plugin[1]; print(plugin)
# p_method <- plugin_index$method[1]; print(p_method)

list4yaml <-
    lapply(X = unique(plugin_index$plugin),
       FUN = function(plugin) {
           t1 <- lapply(X = unique(plugin_index[plugin_index$plugin==plugin,"method"]),
                        FUN = function(p_method) {

                            address <- str_interp("https://docs.qiime2.org/2019.1/plugins/available/${plugin}/${p_method}/")
                            print(address)     
                            feattable <- read_html(address)
                            
                            helpcode <- feattable %>%
                                html_node("div pre") %>%
                                html_text() 
                            helpcode %>% str_split("\\n")
                            plugin_method <- 
                                helpcode %>% 
                                str_match_all(pattern = "qiime ([\\-\\w]+) ([\\-\\w]+)")  %>% 
                                unlist %>% 
                                '['(-1)
                            if (plugin_method[1]!=plugin | plugin_method[2]!=p_method) {
                                stop("Wrong page!!")
                            }
                            helplines <-
                                data.frame(raw_lines=(helpcode %>% 
                                                          str_split(pattern = "\n") %>% unlist)) %>% 
                                as_tibble %>% 
                                mutate(newarg = str_detect(string = raw_lines,
                                                           pattern = "^\\s+\\-\\-") %>% cumsum) %>% 
                                group_by(newarg) %>% 
                                summarise(arg_txt=paste0(raw_lines,collapse = " NL ") %>% 
                                              str_replace(pattern = "\\%",replacement = " NL ")) %>%   # Remove sections beginning with '%', see emperor biplot
                                ungroup %>% 
                                mutate(arg_type=unlist(str_match(string = arg_txt,
                                                                 pattern = "\\-\\-(\\w)\\-")[,2])) %>% 
                                mutate(required = str_detect(string = arg_txt,
                                                             pattern = "required")) %>% 
                                mutate(optional = str_detect(string = arg_txt,
                                                             pattern = "optional")) %>% 
                                # mutate(art_type = str_match(string = arg_txt,pattern = "ARTIFACT PATH (\\w+(\\[.*?\\])?)")[,2]) %>%
                                mutate(art_type = str_match(string = arg_txt,pattern = "ARTIFACT PATH (.*?) NL")[,2]) %>%
                                mutate(arg_name = str_match(string = arg_txt,pattern = "^\\s+([\\-\\w]+)\\s")[,2])
                            
                            method_list <- list()
                            
                            # Required inputs:
                            seclines <- helplines %>% filter(arg_type %in% c("i") & required)
                            # method_list$inputs <- list()
                            # secline <- seclines[1,,drop=F]
                            method_list$inputs <- apply(seclines,1,function(secline) {
                                # Special treatments: 
                                # Split into artifact types. Extreme cases:
                                # 1. "SampleData[PairedEndSequencesWithQuality | SequencesWithQuality]"
                                # 2. "EMPPairedEndSequences | EMPSingleEndSequences | RawSequences"
                                # 3. "SampleData[JoinedSequencesWithQuality] | SampleData[SequencesWithQuality] | SampleData[Sequences]"
                                t6 <- str_match_all(string = secline["art_type"],pattern = "(\\w+)(\\[.*?\\])?\\%?")[[1]][,1]
                                
                                # For each type in list t6, expand if list is inside [] (example 1 above)
                                art_type <- lapply(t6,function(x) {
                                    x[!is.na(str_match(string = x,pattern = "(\\w+)(\\[.*?\\|.*?\\])?")[,3])] <-
                                        x[!is.na(str_match(string = x,pattern = "(\\w+)(\\[.*?\\|.*?\\])?")[,3])] %>% 
                                        lapply(FUN = function(x) {
                                            # cat(sprintf("%s %s\n",plugin,p_method))
                                            t2 <- str_match_all(string = x,pattern = "\\w+") %>% unlist
                                            # t2 <- str_match(string = x,pattern = "(\\w+)\\[(.*?)\\s*\\|\\s*(.*?)\\]")
                                            lapply(t2[-1], 
                                                   FUN = function(x) {paste(t2[1],"[",x,"]",sep="")}) %>% unlist %>% return    
                                        })
                                    return(x)
                                }) %>% unlist#(recursive = F)

                                return(list(art_type))
                            }) %>% unlist(recursive = F)
                            if(!is.null(method_list$inputs)) names(method_list$inputs)<-seclines$arg_name
                            rm(seclines)
                            # optional inputs:
                            seclines <- helplines %>% filter(arg_type %in% c("i") & !required)
                            # secline <- seclines[1,,drop=F]
                            method_list$optional_inputs <- apply(seclines,1,function(secline) {
                                t6 <- str_match_all(string = secline["art_type"],pattern = "(\\w+)(\\[.*?\\])?")[[1]][,1]
                                # For each type in list t6, expand if list is inside [] (example 1 above)
                                art_type <- lapply(t6,function(x) {
                                    x[!is.na(str_match(string = x,pattern = "(\\w+)(\\[.*?\\|.*?\\])?")[,3])] <-
                                        x[!is.na(str_match(string = x,pattern = "(\\w+)(\\[.*?\\|.*?\\])?")[,3])] %>% 
                                        lapply(FUN = function(x) {
                                            # cat(sprintf("%s %s\n",plugin,p_method))
                                            t2 <- str_match_all(string = x,pattern = "\\w+") %>% unlist
                                            # t2 <- str_match(string = x,pattern = "(\\w+)\\[(.*?)\\s*\\|\\s*(.*?)\\]")
                                            lapply(t2[-1], 
                                                   FUN = function(x) {paste(t2[1],"[",x,"]",sep="")}) %>% unlist %>% return    
                                        })
                                    return(x)
                                }) %>% unlist #%>% list
                                
                                return(list(art_type))
                            }) %>% unlist(recursive = F)
                            if(!is.null(method_list$optional_inputs)) names(method_list$optional_inputs)<-seclines$arg_name
                            # seclines <- helplines %>% filter(arg_type=="o")
                            rm(seclines)
                            seclines <- helplines %>% filter(arg_type %in% c("o"))
                            seclines$art_type[!is.na(str_match(string = seclines$art_type,pattern = "(\\w+)(\\[.*?\\|.*?\\])?")[,3])] <-
                                seclines$art_type[!is.na(str_match(string = seclines$art_type,pattern = "(\\w+)(\\[.*?\\|.*?\\])?")[,3])] %>% 
                                lapply(FUN = function(x) {
                                    t2 <- str_match_all(string = x,pattern = "\\w+") %>% unlist
                                    # t2 <- str_match(string = x,pattern = "(\\w+)\\[(.*?)\\s*\\|\\s*(.*?)\\]")
                                    lapply(t2[-1], 
                                           FUN = function(x) {paste(t2[1],"[",x,"]",sep="")}) %>% unlist %>% return    
                                })
                            seclines$art_type[is.na(seclines$art_type)] <- "Visualization"
                            inputs <- c(seclines$art_type) %>% trimws
                            names(inputs) <- seclines$arg_name
                            inputs <- as.list(inputs)
                            method_list$outputs <- inputs
                            
                            rm(seclines)
                            seclines <- helplines %>% filter(arg_type %in% c("p", "m") & required)
                            method_list$required <- c(seclines$arg_name)
                            
                            rm(seclines)
                            seclines <- helplines %>% filter(arg_type %in% c("p", "m") & optional)
                            method_list$optional <- c(seclines$arg_name)
                            
                            # Remove empty elements:
                            method_list <- method_list[lapply(names(method_list),
                                                              FUN = function(x)
                                                                  length(method_list[[x]])!=0) %>%
                                                           unlist]

                            return(method_list)
                        })
           names(t1) <- unique(plugin_index[plugin_index$plugin==plugin,"method"])
           return (t1)
       })
names(list4yaml) <- unique(plugin_index$plugin)

sink("qiime2_arguments_index.yaml")
list4yaml %>% as.yaml %>% cat
sink()



