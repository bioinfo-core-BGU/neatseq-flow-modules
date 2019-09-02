if (any(lapply(X = c("Biostrings","stringr", "tibble", "tidyverse", "optparse", "ggplot2"), 
               FUN = function(x) !(x %in% installed.packages()))))     {
    cat("Make sure the following packages are installed:\n
        install.packages('stringr')\n
        install.packages('tibble')\n
        install.packages('tidyverse')\n
        install.packages('Biostrings')\n
        install.packages('optparse')\n
        install.packages('ggplot2')\n\n")
    stop()
}
    

library(Biostrings)
    
library(stringr)
library(tibble)
library(tidyverse)
library(ggplot2)


library(optparse)
args = commandArgs(trailingOnly=TRUE)


option_list = list(
    make_option(c("-c", "--counts"), 
                type = "character", 
                default = NULL, 
                help = "Path to counts table, one of the three tables created by abundance_estimates_to_matrix.pl", 
                metavar = "character"),
    make_option(c("-g", "--grouping"), 
                type = "character", 
                default = NULL, 
                help = "Path to table of sample groupings. The First column is the samples names (identical to those in --counts header) and a Grouping column ( indicate the name of this column in --grouping_Field  ).", 
                metavar = "character"),
    make_option(c("-d", "--grouping_Field"), 
                type = "character", 
                default = NULL, 
                help = "Name of the Field (Column Name) in the Grouping File (--grouping) to Group the Samples", 
                metavar = "character"),
    make_option(c("-f", "--FASTA"), 
                type = "character", 
                default = NULL, 
                help = "Path to trinity fasta file.", 
                metavar = "character"),
    make_option(c("-m", "--min_count"), 
                type = "numeric", 
                default = NULL, 
                help = "Minimum value of counts table to consider the gene existant in sample.", 
                metavar = "numeric"),
    make_option(c("-s", "--min_sample"), 
                type = "integer", 
                default = 2, 
                help = "Minimum samples per group to consider existant in group", 
                metavar = "integer"),
    make_option(c("-r", "--min_groups"), 
                type = "integer", 
                default = 1, 
                help = "Minimum number of groups to consider the gene existant.", 
                metavar = "integer"),
    make_option(c("-o", "--output"), 
                type = "character", 
                default = NULL, 
                help = "Name of output fasta file", 
                metavar = "character"),
    make_option(c("--gene_map"), 
                type = "character", 
                default = NULL, 
                help = "Path to gene_map file. If passed, will return representative transcript per gene.", 
                metavar = "character"),
    make_option(c("-p", "--plot"), 
                default = FALSE, 
                action = "store_true",
                help = "Plot an histogram of the counts data over all samples and genes"),
    make_option(c("-t", "--transrate"), 
                type = "character", 
                default = NULL, 
                help = "Path to TransRate good contigs fasta file", 
                metavar = "character")
); 



opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", 
                                    description = "",
                                    option_list=option_list,
                                    epilogue="\n\nAuthor: Menachem Sklarz");
opt = optparse::parse_args(opt_parser);




if (any(lapply(X = c("counts", "grouping", "FASTA", "min_count"), 
               FUN = function(x) !(x %in% names(opt)))))     {
    opt_parser %>% print_help
    stop(cat("Please pass values for parameters --counts, --grouping, --FASTA and --min_count\n\n"))
}

# Check fasta file exists:
if(!file.exists(opt$FASTA)) {
    stop(sprintf("FASTA file %s does not exist!\n",opt$FASTA))
}

if (!("output" %in% names(opt))) {
    opt$output = sprintf("%s.filtered.fasta",opt$FASTA)
    cat("No --output passed. Using %s\n\n", opt$output)
}


cat(sprintf("\n[%s] \n\tFinding genes with: \n\t\tmore than %f (counts or normalized counts) \n\t\tin at least %s samples \n\t\tin at least %s treatments. \n\tUsing files: \n\t\tcounts: %s\n\t\tgrouping: %s\n\t\tfasta: %s\n",
            date(),
            opt$min_count, 
            opt$min_sample,opt$min_groups,opt$counts,opt$grouping,opt$FASTA))

if ("gene_map" %in% names(opt)) {
    cat(sprintf("\tAlso, using gene to trans map: %s\n",opt$gene_map))
} 
cat(sprintf("\tOutputting to %s\n",opt$output))

cat(sprintf("\n[%s] Reading table %s\n",date(),opt$counts))
FPKM_tab <- read.delim(opt$counts, he=T)
names(FPKM_tab)[1] <- "Transcript"

if ("gene_map" %in% names(opt)) {
    cat(sprintf("\n[%s] Reading gene map %s\n",date(),opt$gene_map))
    gene_map <- read.delim(opt$gene_map, he=F)
    names(gene_map) <- c("Gene","Transcript")
} 


if (opt$plot) {
    cat(sprintf("Writing plot file to %s\n",sprintf("%s.hist.png",opt$output)))
    png(filename = sprintf("%s.hist.png",opt$output))
    mydf <- FPKM_tab[,-1] %>%
        as.matrix %>%
        as.vector %>%
        data.frame(counts=.)
    
    myplot <- ggplot(mydf,aes(counts)) +
        geom_histogram(bins=150) +
        scale_x_log10( breaks = c(1,10,100,1000)) 
    print(myplot)
        
    
    # hist(log10(FPKM_tab_long$FPKM), 
         # breaks = 100, 
         # freq = FALSE, 
         # main = "Histogram of counts, pooled over samples and genes", 
         # xlab = "log10 of counts", 
         # ylab = "density")
    dev.off()
    cat(sprintf("\n[%s] Created histogram of counts. Stopping\n",date()))
    # stop()

}


cat(sprintf("\n[%s] Reading table %s\n",date(),opt$grouping))
sample_division <- read.delim(opt$grouping,
                              header = T,
                              stringsAsFactors = F,
                              comment.char = "#")
if (opt$grouping_Field %in% colnames(sample_division)){                     
    sample_division = sample_division[,c(colnames(sample_division)[1],opt$grouping_Field ) ]
} else{
    stop(cat("Your grouping_Field could not be found in the grouping file"))
}

colnames(sample_division) = c('SampleID','Treatment')

sample_division$SampleID <- str_replace(string = sample_division$SampleID,
                                        pattern = "^(\\d)",
                                        replacement = "X\\1")
## Testing partial coverage of samples in sample_division:
# sample_division[1:29,] -> sample_division

# Convert to tibble
FPKM_tab <- as.tibble(FPKM_tab)
cat(sprintf("\n[%s] Step 1: Calculate counts >= %f\n",date(),opt$min_count))

FPKM_tab_long <- 
    FPKM_tab                                                 %>% 
    # Convert to long format
    gather(sample, FPKM, -Transcript)                        %>% 
    # Add treatemnt field from sample_division
    inner_join(sample_division,by = c("sample" = "SampleID")) %>% 
    # Add field indicating whther higher than cutoff
    mutate(passCO = FPKM>=opt$min_count)   
    

       
cat(sprintf("\n[%s] Step 2: Calculate samples >= %d\n",date(),opt$min_sample))
FPKM_passCO <- FPKM_tab_long                                 %>% 
    group_by(Transcript,Treatment)                           %>% 
    summarise(treatPass = sum(passCO) >= opt$min_sample)     %>% 
    group_by(Transcript)                                     %>% 
    summarise(genePassed = sum(treatPass))


# Create fasta file index
# Doing this here so that the seqlength can be used for transcript selection when gene_map is passed
cat(sprintf("\n[%s] Reading original FASTA into index.\n", date()))
fastaindex <- fasta.index(filepath =  opt$FASTA)
# Extract gene name from fasta header (desc)
fastaindex$gene_name <-
    fastaindex$desc %>% 
    str_extract(string = .,
                pattern = ".*?\\s") %>% 
    str_trim(string = .,side = "right")


if ("gene_map" %in% names(opt)) {
    
    # Getting sequence lengths from fastaindex
    cat(sprintf("\n[%s] Extracting sequence lengths from FASTA index.\n", date()))
    seq_len_index <-
        fastaindex %>% 
        as_tibble %>% 
        select(gene_name,seqlength) %>% 
        rename(Transcript=gene_name)
    
    cat(sprintf("\n[%s] Finding one transcript per gene passing filters\n", date()))
    passed_genes <- 
        FPKM_passCO %>% 
        # Keep only transcripts which pass min_group filter
        filter(genePassed>=opt$min_groups) %>% 
        # Add reads per sample
        left_join(FPKM_tab_long) %>%
        # Calcultate sum of reads per Transcript (over samples)
        group_by(Transcript) %>%
        summarise(sum_reads=sum(FPKM)) %>% 
        # Add gene names from gene_map
        left_join(gene_map)  %>% 
        group_by(Gene) %>% 
        # Keep only highest expressed transcript per gene 
        filter(sum_reads==max(sum_reads)) %>% 
        # If tied, keep longest transcript
        left_join(seq_len_index) %>% 
        filter(seqlength==max(seqlength)) %>% 
        # If still tied, keep transcript with lowest name
        summarise(Transcript=min(as.character(Transcript))) %>% 
        ungroup() %>% 
        # Convert to character vector
        .$Transcript
    cat(sprintf("\n[%s] Step 3: Calculate treatments >= %d\n",date(),opt$min_groups))
    cat(sprintf("\n[%s] Number of genes with transcripts that passed filter: %d\n",date(),length(passed_genes)))
    
} else {
    passed_genes <- as.character(FPKM_passCO$Transcript[FPKM_passCO$genePassed>=opt$min_groups])
    cat(sprintf("\n[%s] Step 3: Calculate treatments >= %d\n",date(),opt$min_groups))
    cat(sprintf("\n[%s] Number of transcripts that passed filter: %d\n",date(),length(passed_genes)))
}

# Get intersect of selected genes and TransRate good genes:
if ("transrate" %in% names(opt)) {
    cat(sprintf("\n[%s] Reading transrate 'good' FASTA into index.\n", date()))
    fastaindex <- fasta.index(filepath =  opt$transrate)
    cat(sprintf("\n[%s] Intersecting transrate 'good' contigs with passed contigs.\n", date()))
    passed_genes <- intersect(passed_genes,fastaindex$desc)
}



# For each blast result:

t1 <- fastaindex[fastaindex$gene_name %in% passed_genes,]
newseq <- readBStringSet(t1)
cat(sprintf("\n[%s] Writing new fasta to %s\nNumber of sequences: %s\n\n",
            date(),
            opt$output,
            length(newseq)))
writeXStringSet(newseq, 
                filepath = opt$output, 
                compress=FALSE, 
                format="fasta")

cat(sprintf("\n[%s] Writing new counts table to %s\n",
            date(),
            paste0(opt$output,".tab")))
FPKM_tab %>% 
    filter(Transcript %in% passed_genes) %>% 
    write.table(file = paste0(opt$output,".tab"),quote = F,sep = "\t",row.names = F)


len_dist <- rbind(data.frame(type="original",length=fastaindex$seqlength),
                  data.frame(type="filtered",length=t1$seqlength))

cat(sprintf("\n[%s] Creating length distribution plot and writing to %s\n",
            date(),
            paste0(opt$output,".png")))

ggplot(len_dist, aes(x=length, fill=type)) +
    scale_x_log10() +
    ggtitle("Distribution of lengths before and after filtering") +
    xlab("log10 of length") +
    ylab("Density") +
    geom_histogram(aes(y=0.001*..density..),binwidth=0.01)+
#    geom_vline(xintercept=c(290,695))
ggsave(paste0(opt$output,".png"), 
       plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 15, height = 10, units = "cm",
       dpi = 300, limitsize = TRUE)

cat("Medians of lengths before and after filtering:")
print(tapply(len_dist$length,
            INDEX = len_dist$type,
            FUN = median))
