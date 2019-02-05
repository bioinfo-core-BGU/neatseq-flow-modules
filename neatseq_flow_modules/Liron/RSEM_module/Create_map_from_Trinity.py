import os, re 
import argparse
from Bio import SeqIO
import sys
import pandas as pd

# parser = argparse.ArgumentParser(description='Create map file from Trinity')

# parser.add_argument('-i', type=str,dest='input',
                    # help='Input FASTA File')
# parser.add_argument('-o', dest='Output', type=str,default=os.path.join(os.getcwd()),
                    # help='Output file')
# args = parser.parse_args()

args=pd
args.input=sys.argv[1]
args.Output=sys.argv[2]
record_dict = SeqIO.index(args.input, "fasta")
TEMP_list=list(record_dict)
lines=[re.split("_i.+$",x)[0]+"\t"+x+"\n" for x in TEMP_list]
h_file=open(args.Output,"w")
h_file.writelines(lines)
h_file.close