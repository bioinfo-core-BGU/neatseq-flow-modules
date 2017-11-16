
library(magrittr)
library(hash)

# The purpose of this script is two fold:
# A. Merge KRAKEN output files into a biom table that can be converted to biom format and used as input for phyloseq or QIIME
# B. Convert the ncbi taxonomy identifiers into a taxonomy table that can be read by phyloseq | QIIME

# Getting a list of arguments as follows:
# First argument: The name of the output file prefix. Will create <prefix>.count.tsv and <prefix>.taxonomy.tsv
# Second argument: The suffix of the report files to remove from the filenames to obtain the sample name
# The rest of the argument list: files with the above mentioned suffix to use in the analysis.

library("optparse")
library(tools)

args = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
    make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="An existing dir in which to create the output files. These will be of the form kraken_biom.* [index_file location]", metavar="character"),
    make_option(c("-d", "--db"), type="character", default=NULL, 
              help="path to database used for kraken", metavar="character"),
  make_option(c("--redo_hashes"), default=FALSE, action = "store_true", 
              help="recalculate the taxonomy lookup hashes", metavar="character")

              
); 

opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", option_list=option_list,epilogue="\n\nAuthor: Menachem Sklarz");
opt = optparse::parse_args(opt_parser);

errors_exist = FALSE
if (is.null(opt$file)){
    
    errors_message = "Please supply an index file of kraken reports (-f)\n"
    errors_exist = TRUE
}
if (is.null(opt$db)){
    
    errors_message = paste0(errors_message,"Please supply the database used by kraken (--db)\n")
    errors_exist = TRUE
}

if(errors_exist) {
    print_help(opt_parser)
    cat(errors_message)
    stop()
}


if (is.null(opt$outdir)){
    cat("Using ",dirname(normalizePath(opt$file))," as output directory!\n")
    opt$outdir <- dirname(normalizePath(opt$file))
}



out.dir <- opt$outdir
sample_file <- opt$file


##########
# # For debuggiing:
if(FALSE) {
    sample_file = "/gpfs0/bioinfo/projects/SPICE/20.Mekorot/01.Metagenomics/03.May2017_workflow/data/kraken/kraken_k21/kraken_files_index.txt"
    out.dir = "/gpfs0/bioinfo/projects/SPICE/20.Mekorot/01.Metagenomics/03.May2017_workflow/r_temp"
    opt$db = "/gpfs0/bioinfo/databases/kraken_databases/stdDB_k21"
    opt$redo_hashes = TRUE
}
#########

out.dir <- sub("/$","",out.dir, perl=T)
if (!file.exists(out.dir)) {
    print_help(opt_parser)
    cat("Please supply a legitimate output directory in --outdir parameter\n")
    stop()
}
out.prefix <- sprintf("%s/%s",out.dir,"kraken_biom")

if (!file.exists(sample_file)) {
    cat ("Please supply a legitimate file with an index of kraken output files in --file parameter")
    # print_help(opt_parser)
    stop()
}
sample_file <- read.delim(sample_file,he=T)
filelist <- sample_file$kraken_report
sample_names <- sample_file[,1]

print(sample_file)

##########
# # For debugging:
# filelist <- h(filelist)
# sample_names <- h(sample_names)
##########


merge_kraken_files <- function(filelist,sample_names) {
    # Read all files (set the for-loop accordingly) and merge them by taxid
    count=0
    # for (rep_file in dir(pattern="*report")) {
    print("Reading and merging report files")

    for (i in 1:length(filelist)) {
        rep_file <- as.character(filelist[i])
        sample_name <- sample_names[i]
        sprintf("Reading file %i of %i: %s\n",i,length(filelist),rep_file) %>% cat
        # print(paste("\tReading file:",rep_file,sep=" "))
        if(count==0) {
            krak_mer = read.delim(rep_file,he=F)
            names(krak_mer) <- paste(c("Perc_reads","Num_reads_clade","Num_reads","rank_code","taxid","sci_name"),sample_name,sep="_")
            names(krak_mer)[5] <- "taxid"
            count = count+1;
        } else {
            kraken_rep = read.delim(rep_file,he=F)
            names(kraken_rep) <- paste(c("Perc_reads","Num_reads_clade","Num_reads","rank_code","taxid","sci_name"),sample_name,sep="_")
            names(kraken_rep)[5] <- "taxid"
            krak_mer <- merge(krak_mer,kraken_rep,by="taxid",all=T)
        }

    }
    return(krak_mer)
}
print("Preparing count table")

krak_mer <- merge_kraken_files(filelist,sample_names);

# Editting count table
# (make sure this is general enough)
# count_table <- krak_mer[,c("taxid",paste("Num_reads_clade",sample_names,sep="_"))]
# Using Num_reads rather than Num_reads_clade because they are later summed up for each rank (in phyloseq and the like)
count_table <- krak_mer[,c("taxid",paste("Num_reads",sample_names,sep="_"))]
count_table[is.na(count_table)] <- 0
# Converting integers to floats. This is possibly a requirement by "biom convert"
count_table[,-1] <- as.data.frame(apply(count_table[,-1],2,function(x) formatC(x,digits=1,format="f")))
    taxonomy <- krak_mer[,c("taxid",paste("sci_name",sample_names,sep="_"))]
    taxonomy <- data.frame(taxid=taxonomy[,1],sci_name=apply(taxonomy[,-1],1,function(x) {x[which(!is.na(x))[1]]}))
count_table$taxonomy <- taxonomy$sci_name

taxrank <- krak_mer[,c("taxid",paste("rank_code",sample_names,sep="_"))]
taxrank <- data.frame(taxid=taxrank[,1],sci_name=apply(taxrank[,-1],1,function(x) {x[which(!is.na(x))[1]]}))

# Write count table
writeLines(text=paste(c("#OTU ID",sample_names,"taxonomy"),collapse="\t"),con=paste(out.prefix,".count.tsv",sep=""))
write.table(count_table,file=paste(out.prefix,".count.tsv",sep=""),sep="\t",col.names=F,row.names=F,append=T,quote=F)


# The command for converting the file into a json based biom file (the files shouldn't be so big as to require hdf5)
sprintf("# Use the following commad to convert the tsv file into biom format:\n") %>% cat
sprintf("qiimeVE biom convert -i %s -o %s  --table-type='OTU table' --to-json --process-obs-metadata taxonomy\n\n",
        paste(out.prefix,".count.tsv",sep=""),
        paste(out.prefix,".count.biom",sep="")) %>% cat


############
# Function for creating or loading taxonomy hashes
load_hashes <- function(path) {
    path = paste0(path,"/taxonomy")
    # Check if hash file exists in "path". If it does, load the file. If not, create it
    if(file.exists(paste(path,"hashes.R",sep="/")) & !opt$redo_hashes) {
        print("Loading hashes:  Takes a few minutes. Please be patient...")
        system.time(load(paste(path,"hashes.R",sep="/"))) %>% print
    } else {
        print("hashes file does not exist. Creating the hashes from original files. Might take time...")
        print("Getting taxonomy rank info from NCBI. This might take a while")
        print("Note that taxonomy files from kraken db are being used. This might produce an error if the files are moved")
        # Reading NCBI taxonomy files from kraken db:
        # Nodes file:
        print("Reading NCBI nodes.dmp file")
        tax_node <- read.table(paste(path,"nodes.dmp",sep="/"),he=F,sep="|")
        tax_node <- tax_node[,1:3]
        names(tax_node) <- c("taxid","parent","level")
        tax_node$level <- gsub("\t","",tax_node$level)
        tax_node <- rbind(c(0,0,"no rank"),tax_node)

        # Names file
        print("Reading NCBI names.dmp file")
        tax_name <- read.table(paste(path,"names.dmp",sep="/"),he=F,sep="|",quote = "",comment.char = "")
        tax_name <- tax_name[,c(1,2,4)]
        names(tax_name) <- c("taxid","name","type")
        tax_name$name <- gsub("\t","",tax_name$name)
        tax_name$type <- gsub("\t","",tax_name$type)
        tax_name <- tax_name[(tax_name$type=="scientific name"),-3]
        tax_name <- rbind(c(0,"unclassified"),tax_name)
        rownames(tax_name) <- tax_name$taxid

        # Using hash for fast searches:
        print("preparing hash tables for NCBI taxonomy data:")
        print("1 of 3")
        tax_node_hash_parent <- hash(tax_node$taxid,tax_node$parent)
        print("2 of 3")
        tax_node_hash_level <- hash(tax_node$taxid,tax_node$level)
        print("3 of 3")
        tax_name_hash <- hash(tax_name$taxid,tax_name$name)
        tax_name_hash[["-"]] <- "unclassified"

        rm(tax_node)
        rm(tax_name)
        print("Saving hashes: Takes a few minutes. Please be patient...")
        system.time(save(list=c("tax_node_hash_parent","tax_node_hash_level","tax_name_hash"),file=paste(path,"hashes.R",sep="/"))) %>% print
    }
    return(list(hash_parent=tax_node_hash_parent,
                hash_level=tax_node_hash_level,
                hash_name=tax_name_hash))
}

# Get the hashes from load_hashes function:
hash_list <- load_hashes(path = opt$db);
# Retrieve the hashes from the returned list:
tax_node_hash_parent <- hash_list$hash_parent
tax_node_hash_level <- hash_list$hash_level
tax_name_hash <- hash_list$hash_name
rm(hash_list)


###############################
# New method for creating taxonomy table:

# Function for finding relevant parent of taxid (relevant = belongs to one of: "kingdom","phylum","class","order","family","genus","species")
# Returns a list of structure: list(parent=parent,parent_lev=level) containing the taxid of the parent and it's level
# Example: find_parent(1706371)
find_parent <- function(taxid) {
    if(taxid==0) return (list(0,"unclassified"))
    if(taxid==1) return (list(1,"root"))
    taxid <- as.character(taxid)
    level <- NA;
    while(is.na(level)) {
        parent <- tax_node_hash_parent[[taxid]]
        level <- tax_node_hash_level[[as.character(parent)]]
        # if(level=="superkingdom") level <- "kingdom";
        if(parent==1) return(list(parent=1,parent_lev="no rank"))
        if (!(level %in% c("superkingdom","kingdom","phylum","class","order","family","genus","species"))) {
#             print(level)
#             print(taxid)
#             print(parent)
            taxid <- as.character(parent);
            level <- NA     # Continue loop until relevant taxonomy found
        }
        
    }
    return(list(parent=parent,parent_lev=level))
}


# Level name index: (fast conversion level name to one-letter code)
level_ind <- c("u","r","d","d","p","c","o","f","g","s")
names(level_ind) <- c("unclassified","root","superkingdom","kingdom","phylum","class","order","family","genus","species")
# Note kingdom and superkingdom are both coded by "d"

existing_taxids <- krak_mer[,"taxid"] %in% keys(tax_node_hash_level)
if(any(existing_taxids==FALSE)) stop("Not all taxids in the report could be identified! Make sure you are using the correct database. Pass --db with the same db used for kraken.")


# krak_mer <- krak_mer[existing_taxids,]
# Creating taxonomy structure:
taxonomy <- cbind(taxid=krak_mer[,"taxid"],
                  level=values(tax_node_hash_level,
                               keys=krak_mer[,"taxid"])
                  ) %>% as.data.frame
taxonomy$level <- gsub("superkingdom","kingdom",taxonomy$level) # Convert all "superkingdom" to "kingdom"
taxonomy <- taxonomy[,c("taxid","level")]
# Loop over all lines and replace all non-standard levels with the closest taxid above the level with a standard level
# e.g taxid 323
# Didn't manage with apply for some reason. Looping. Not too bad...
t1 <- data.frame(taxid = numeric(dim(taxonomy)[1]),
                 level = character(dim(taxonomy)[1]),
                 name = character(dim(taxonomy)[1]))
for(i in 1:dim(taxonomy)[1]) {
    # if((i %% 200) == 0) cat(i)

    x <- taxonomy[i,]
    if(x[1] == 0) {
        t1[i,1:2] <- c(0,"unclassified")
    } else if(x[2] %in% c("superkingdom","kingdom","phylum","class","order","family","genus","species")) {
        t1[i,1:2] <- x[1:2]
    } else {
        t1[i,1:2] <- find_parent(x[1])%>%unlist
    }
    if((i %% 200) == 0) sprintf("%i: %s\n",i, tax_name_hash[[t1$taxid[i]]]) %>% cat
    t1$name[i] <- tax_name_hash[[t1$taxid[i]]]
}   
# See what this does with following command: 
# cbind(t1,taxonomy)%>%head(100)
taxonomy <- t1
taxonomy$level <- gsub("superkingdom","kingdom",taxonomy$level)
rm(t1)

taxonomy$rank_code <- level_ind[taxonomy$level]
table(taxonomy$level,taxonomy$rank_code)

##### --------------------

# Creating columns in taxonomy to recieve taxonomy info:
taxonomy[,c("u","r","d","p","c","o","f","g","s")] <- NA
# Loop over lines and store taxid in correct column in taxonomy based on the rank_code
print("Creating taxonomy file:")
for (i in 1:length(taxonomy$rank_code)) {
    if((i %% 200) == 0) sprintf("Line %i of %i",i,length(taxonomy$rank_code))
    # Ignore this kind of level:
    if(is.na(taxonomy$rank_code[i]) | taxonomy$rank_code[i] == "-") {
        next;
    }
    taxonomy[i,taxonomy$rank_code[i]] <- taxonomy$taxid[i]
}
# Loop over lines again and for each level fill in all the levels above it in one of two ways:
# 1. The slow track: Use NCBI nodes.dmp file to find the taxids parents. 
# 2. The fast track: If the taxonomy has been identified, use it rather than the slow track
taxa_hash <- hash()
for (i in 1:length(taxonomy$rank_code)) {
    if((i %% 200) == 0) sprintf("# %i Level: %s\n",i,taxonomy[i,"level"]) %>% cat
    # Ignore this kind of level:
    if(taxonomy$rank_code[i] %in% c(NA,"-","u","r")) {
        # print("Skipping -ud & NA")
        next;
    }
    # j is the column number with known taxid. The mission is to find taxids for all columns to the left of j.
    j <- which(c("d","p","c","o","f","g","s") == taxonomy$rank_code[i])
    # Looping backwards from column j to column 1 (=right to left)
    for (k in j:1) {
        if(is.na(taxonomy[i,c("d","p","c","o","f","g","s")[k]])) next;
        parent <- find_parent(taxonomy[i,c("d","p","c","o","f","g","s")[k]])
        if (parent$parent_lev == "no rank") parent$parent_lev <- "root"
        taxonomy[i,level_ind[(parent$parent_lev)]] <- parent$parent

    }
    # Building taxa_hash to use in fast track, which isn't implemented at the moment...
#     if(is.null(taxa_hash[[as.character(taxonomy[i,c("d","p","c","o","f","g","s")][j])]])) {
#         taxa_hash[[as.character(taxonomy$taxid[i])]] <- taxonomy[i,c("d","p","c","o","f","g","s")[1:j]]
#     }
    
}

print("Done...")
print("Finding taxid names:")

# # Helper function. Given a vector of taxids returns the NCBI names of those taxids
# get_tax_name <- function(x) {
  # x[is.na(x)] <- 0;
  # return(tax_name[unlist(x) %>% as.character,"name"]);
# }

# # Apply get_tax_name to each line in taxonomy (Note the t() and conversion to dataframe at end of command.
# final_tax <- apply(taxonomy[,c("d","p","c","o","f","g","s")],1,get_tax_name) %>% t %>% data.frame

# Converting NA's into "-" for all taxonomy level columns:
taxonomy[,c("r","d","p","c","o","f","g","s")][is.na(taxonomy[,c("r","d","p","c","o","f","g","s")])] <- "-"

############# 
# Critical. Converting all level numbers into level names using hash::values
# Also transposing and converting to df
final_tax <- apply(taxonomy[,c("r","d","p","c","o","f","g","s")],
                    MARGIN = 1,
                    FUN = function(x) {
                        # print(x)
                        hash::values(tax_name_hash,keys=x)
                    }) %>% t %>% data.frame

# Add an "unclassified" column
# final_tax <- cbind(unclassified="-",final_tax)
names(final_tax) <- names(level_ind)[c(-1,-3)]
# Add taxid as first column in final_tax:
final_tax <- cbind(taxid = taxonomy$taxid,
                    taxonomy = taxonomy$name,
                    final_tax)

# Sum up duplicate taxa names ()
cols_to_ignore <- which(count_table %>%names %in% c("taxid","taxonomy"))

t1 <- merge(count_table,final_tax,all.x=F,all.y=T)
t1 <- t1[,names(count_table)]
count_table_sum <- apply(count_table[,-(cols_to_ignore)],
                         MARGIN = 2,
                         FUN = function(x) tapply(as.numeric(x),
                                                  INDEX=final_tax$taxid,
                                                  FUN = sum)
                        ) %>% as.data.frame

final_tax_sum <- apply(final_tax[,-1],
                       MARGIN = 2,
                       FUN = function(x) tapply(x,
                                                INDEX=final_tax$taxid,
                                                FUN = max)
                        ) %>% as.data.frame

## This is a fix!!
# final_tax_sum$root[final_tax_sum$taxonomy=="root"] <- "root"

# Apply after discarding taxonomy and species level:
final_tax_sum$taxonomyQIIME <- apply(final_tax_sum[,setdiff(names(final_tax_sum), 
                                                            c("taxonomy","root", "species"))],  
                                     1,
                                     function(x) {
                                             x[x=="unclassified"] = ""
                                             paste(paste(substring(names(x),
                                                                   first = 1,
                                                                   last = 1),
                                                         x,
                                                         sep="__"),
                                                   collapse=";") %>% return 
                                     })

# Add taxid and taxonomy to count_table_sum 
count_table_sum$taxonomy <- final_tax_sum$taxonomyQIIME
count_table_sum$taxid <- row.names(count_table_sum)
# Change order so that "taxid" is first and "taxonomy" is last
taxid_col <- which(names(count_table_sum) == "taxid")
taxonomy_col <- which(names(count_table_sum) == "taxonomy")

count_table_sum <- count_table_sum[,c(taxid_col,
                                      setdiff(1:(dim(count_table_sum)[2]),
                                              c(taxonomy_col,
                                                taxid_col)),  # All other columns. This is a good place to put them in a user defined order!
                                      taxonomy_col)] 

      
# Write taxonomy table and updated count tables
print("Writing taxonomy file:")
writeLines(text = paste(c(names(final_tax_sum)),
                        collapse="\t"),
           con = paste(out.prefix,
                       ".taxonomy.sum.tsv",
                       sep=""))
write.table(final_tax_sum,
            file = paste(out.prefix,".taxonomy.sum.tsv",sep=""),
            sep="\t",
            col.names=F,
            row.names=F,
            append=T,
            quote=F)

# Write count table
print("Writing updated count file:")
writeLines(text=paste(c("#OTU ID",
                        sample_names,
                        "taxonomy"),
                      collapse="\t"),
           con = paste(out.prefix,".count.sum.tsv",sep=""))
write.table(format(count_table_sum,
                   digits = 10,
                   nsmall = 1),
            file = paste(out.prefix,".count.sum.tsv",sep=""),
            sep="\t",
            col.names=F,
            row.names=F,
            append=T,
            quote=F)


# The command for converting the file into a json based biom file (the files shouldn't be so big as to require hdf5)
sprintf("# Use the following commad to convert the tsv file into biom format:");
sprintf("# Note - you might have to set up QIIME environmental variables. Don't ask...");
conv_CMD <- sprintf("biom convert -i %s -o %s  --table-type='OTU table' --to-json --process-obs-metadata taxonomy\n\n",
                    paste(out.prefix,".count.sum.tsv",sep=""),
                    paste(out.prefix,".count.sum.biom",sep=""))
cat(conv_CMD)
system(conv_CMD)


conv_CMD <- sprintf("biom summarize-table -i %s -o %s \n\n",
                    paste(out.prefix,".count.sum.biom",sep=""),
                    paste(out.prefix,".count.sum.biom.summary",sep=""))
cat(conv_CMD)
system(conv_CMD)






# Uploading to phyloseq
print("Printing phyloseq commands to....")

# Creates code to load the data files into phyloseq
sprintf("
library(phyloseq)
mybiom <- import_biom(\"%s\")
final_tax_sum <- read.table(file=\"%s\",sep=\"\\t\",stringsAsFactors=F, quote = \"\",he=T)
taxid_ind <- read.table(file=\"%s\",
                        sep=\"\\t\",quote = \"\",
                        stringsAsFactors=F,
                        colClasses=c(NA,rep(\"NULL\",times=dim(mybiom@otu_table)[2]),NA))


final_tax_mat <- as.matrix(final_tax_sum)
rownames(final_tax_mat) <- taxid_ind[,1]
final_tax_mat <- final_tax_mat[,-1]
tax_table(mybiom) <- tax_table(final_tax_mat)
\n",
paste(out.prefix,".count.sum.biom",sep=""),
paste(out.prefix,".taxonomy.sum.tsv",sep=""),
paste(out.prefix,".count.sum.tsv",sep="")
) -> phyloseq_cmd

writeLines(phyloseq_cmd,sprintf("%s_phyloseq_script.R",out.prefix))



# t1 <- mybiom
# t1 = filter_taxa(t1, function(x) mean(x) > 100, TRUE)
# # Remove very large = unclassified
# # t1 = filter_taxa(t1, function(x) mean(x) < 1000, TRUE)
# # t1 = transform_sample_counts(t1, log)

# toplot <- tax_glom(t1, taxrank="phylum", NArm=TRUE)

# jpeg("phyla.jpg",width=15,height=7,units="in",res=200)
# plot_bar(toplot,fill="phylum")
# dev.off()



