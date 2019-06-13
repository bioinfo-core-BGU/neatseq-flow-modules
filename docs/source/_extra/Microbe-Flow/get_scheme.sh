#!/bin/bash

echo 'AUTHOR: Liron Levin'
echo 'Usege: sh get_scheme.sh [ORGANISM_NAME] [SCHEME_TYPE]'

ORGANISM=$1
BASE_URL=http://rest.pubmlst.org/
SCHEM_TYPE=$2

DB_URL=$(curl -s   $BASE_URL | awk 'BEGIN { RS="[}]"; } {print ;}'| grep -Ei "$ORGANISM" |grep -Ei "sequence/profile" | grep -Ei -o "http[ / \. _ : A-z 0-9 & =]+")

if [ -n $DB_URL ]
then

    SCHEM_URL=$DB_URL"/schemes"

    SCHEM_COUNT=$(curl -s $SCHEM_URL | awk 'BEGIN { RS="[}]"; } {print ;}'| grep -Ei '"description".+'$SCHEM_TYPE |grep -Ei "scheme" | grep -Eic -o "http[ / \. _ : A-z 0-9 & =]+")
    if [ $SCHEM_COUNT -gt 1 ]
    then
        SCHEM_COUNT=$(curl -s  $SCHEM_URL | awk 'BEGIN { RS="[}]"; } {print ;}'| grep -Ei '"description".+"'$SCHEM_TYPE'"' |grep -Ei "scheme" | grep -Eic -o "http[ / \. _ : A-z 0-9 & =]+")
        SCHEM_URL=$(curl -s  $SCHEM_URL | awk 'BEGIN { RS="[}]"; } {print ;}'| grep -Ei '"description".+"'$SCHEM_TYPE'"' |grep -Ei "scheme" | grep -Ei -o "http[ / \. _ : A-z 0-9 & =]+")        
    else
        SCHEM_COUNT=$(curl -s  $SCHEM_URL | awk 'BEGIN { RS="[}]"; } {print ;}'| grep -Ei '"description".+'$SCHEM_TYPE |grep -Ei "scheme" | grep -Eic -o "http[ / \. _ : A-z 0-9 & =]+")
        SCHEM_URL=$(curl -s  $SCHEM_URL | awk 'BEGIN { RS="[}]"; } {print ;}'| grep -Ei '"description".+'$SCHEM_TYPE |grep -Ei "scheme" | grep -Ei -o "http[ / \. _ : A-z 0-9 & =]+")
        
    fi

    if [ $SCHEM_COUNT -eq 1 ]
    then
        echo "Found " $(curl -s -S  $SCHEM_URL | grep -Eio 'description":"[ a-Z \. : , 0-9]+' | sed 's/^description":"//') "scheme"
        SCHEM_FILE=$(echo $ORGANISM"_"$SCHEM_TYPE"_scheme.tsv" | sed 's/ //')
        if [ $(curl -s -S  $SCHEM_URL | grep -Eioc "/profiles_csv") -gt 0 ]
        then 
            echo "Downloading profiles..."
            curl  -s -S  $SCHEM_URL"/profiles_csv" >  $SCHEM_FILE 
        else
            echo 'No profiles found for this scheme!!!'
            echo 'Generates profiles less scheme file'
            echo -e $SCHEM_TYPE"\t\c" > $SCHEM_FILE
            curl -s -S  $SCHEM_URL | awk 'BEGIN { RS="[}]"; } {print ;}' | grep -Ei -o "loci/\w+" | sed 's/^loci\///' | tr '\n' '\t' >>  $SCHEM_FILE 
        fi
        if [ -s $SCHEM_FILE ]
        then 
            LOCI_DATA=$(curl -s -S  $SCHEM_URL)
            LOCI_NUM=$(echo  $LOCI_DATA | awk 'BEGIN { RS="[}]"; } {print ;}'  | grep -Ei -o "http[ / \. _ : A-z 0-9 & =]+" | grep -Eic "loci")
            loci_URL=$(echo  $LOCI_DATA | awk 'BEGIN { RS="[}]"; } {print ;}'  | grep -Ei -o "http[ / \. _ : A-z 0-9 & =]+" | grep -Ei "loci"|  awk '{print $0"/alleles_fasta";}' )
            echo "Found "$LOCI_NUM " Locus "
            if [  $LOCI_NUM  -gt 0 ]
            then 
                echo "Downloading Alleles..."
                ALLELE_FILE=$(echo $ORGANISM"_"$SCHEM_TYPE"_FASTA.fasta" | sed 's/ //')
                curl $loci_URL > $ALLELE_FILE
                if [ -s $ALLELE_FILE ]
                then 
                    ALLELE_MAP_FILE=$(echo $ORGANISM"_"$SCHEM_TYPE"_Alleles.tab" | sed 's/ //')
                    #grep '>' $ALLELE_FILE | sed 's/^>//' | awk 'BEGIN { FS="_"; print "Allele\tGene\tNumber"; } {print $1"_"$2"\t"$1"\t"$2;}' > $ALLELE_MAP_FILE
                    grep '>' $ALLELE_FILE | sed 's/^>//' | awk 'BEGIN { FS="_" ; print "Allele\tGene\tNumber"; } { temp=sprintf("%s",$0); temp_NF=sprintf("%s",$NF); sub("_"$NF,"",$0);} { print temp"\t"$0"\t"temp_NF;}' > $ALLELE_MAP_FILE                     
                else
                    echo "Could not download the alleles sequences of scheme " $SCHEM_TYPE " for " $ORGANISM
                fi
                
            else
                echo "No locus data for " $SCHEM_TYPE " scheme found for " $ORGANISM
            fi    
        else
            echo "Could not download the " $SCHEM_TYPE " scheme profiles for " $ORGANISM
        fi
    else
        echo "No scheme " $SCHEM_TYPE " found for " $ORGANISM
    fi
else
    echo ' Could not fined the organism database '
fi
# echo $DB_URL
# echo $SCHEM_URL
