export DBDIR="/gpfs0/bioinfo/users/sklarz/Metagenomics_4_readthedocs/DBdir"

export PERL5LIB=""

# Metaphlan
metaphlan2.py \
 --input_type fastq \
 --bowtie2_exe bowtie2 \
 --bowtie2db $DBDIR/MetaPhlAn

# Kraken
mkdir -p $DBDIR/kraken2
kraken2-build \
    --standard \
    --threads 10 \
    --db $DBDIR/kraken2

# Krona
mkdir -p $DBDIR/krona
ktUpdateTaxonomy.sh $DBDIR/krona/taxonomy

# Kaiju
mkdir -p $DBDIR/kaiju
cd $DBDIR/kaiju
kaiju-makedb -s progenomes -t 10
kaiju-makedb -s nr_euk -t 10
cd -

# HUMAnN2
mkdir -p databases/HUMAnN2
humann2_databases --download chocophlan full  $DBDIR/HUMAnN2
humann2_databases --download uniref uniref90_diamond  $DBDIR/HUMAnN2/uniref90
humann2_databases --download uniref uniref50_diamond  $DBDIR/HUMAnN2/uniref50

humann2_config --update database_folders nucleotide $DBDIR/HUMAnN2/chocophlan
humann2_config --update database_folders protein $DBDIR/HUMAnN2/uniref90
