from __future__ import division
import os, re
import pandas as pd
import argparse
from functools import partial
from multiprocessing import Pool
import time
import sys

parser = argparse.ArgumentParser(description='Cog file maker for Gecko3')

parser.add_argument('-D', type=str,dest='DIR',default=os.getcwd(),
                    help='Directory of GFF files')
parser.add_argument('-C', dest='clusters',default="clustered_proteins",
                    help='A Cluster file from Roary program [clustered_proteins]')
parser.add_argument('-B', dest='Bicluster',
                    help='A Bicluster clusters file from the Biclustering analysis')
parser.add_argument('-P', dest='processes',default=1,type=int,
                    help='The number of CPUs to use')
parser.add_argument('-o', type=str,dest='out',default="Gecko_cog.cog",
                    help='output cog file')

args = parser.parse_args()

def run_on_GFF(gff_file,args,clusters):
    cog_file=[]
    Data=pd.read_table(os.path.join(args.DIR,gff_file),comment='#',header=None,engine="python")
    Data=Data.loc[~Data[2].isnull(),]
    for contig in set(Data[0]):
        cog_file.append(gff_file.strip(".gff")+ " Contig_" + contig)        
        contig_Data=Data.loc[Data[0]==contig,].copy()
        contig_Data=contig_Data.loc[contig_Data[2]=="CDS",].copy()
        contig_Data["locus_tag"]=map(lambda x: re.findall("ID=[A-z 0-9 _ . ,]+",x)[0].strip('ID='),contig_Data[8])
        contig_Data["product"]=map(lambda x: re.findall("product=\S+",x)[0].strip('product='),contig_Data[8]) 
        contig_Data["NAME"]=map(lambda x: re.findall("Name=[A-z 0-9 _ . ,]+",x)[0].strip('Name=') if "Name=" in x else re.findall("ID=[A-z 0-9 _ . ,]+",x)[0].strip('ID='),contig_Data[8])
        
        cog_file.append(str(len(contig_Data))+ " proteins" )
        cog_file=cog_file+map(lambda x: clusters.setdefault(contig_Data.loc[x,"locus_tag"],'0').replace("group_","")+"\t"+contig_Data.loc[x,6]+"\t?\t"+contig_Data.loc[x,"NAME"]+"\t" +
                              contig_Data.loc[x,"product"] +"\t" +contig_Data.loc[x,"locus_tag"]+"\t" +contig_Data.loc[x,"product"]  ,contig_Data.index)
        cog_file.append("")
    return cog_file

print "Reading Clusters file..."

with open(args.clusters,"rb") as h_file:
    clustered_proteins=h_file.readlines()
    h_file.close()
clusters={}
map(lambda x: map(lambda y: clusters.setdefault(y,clusters[y]+","+x.split(": ")[0]) if clusters.has_key(y) else  clusters.setdefault(y,x.split(": ")[0])   ,x.split(": ")[1].split("\t"))       ,clustered_proteins)
print "Done Reading Clusters file!"

cog_file=[]
if args.Bicluster!=None:
    Bicluster_h=open(args.Bicluster,"rb")
    Bicluster=Bicluster_h.readlines()
    Bicluster_h.close()
    Bicluster=map(lambda x: x.strip().split("\t"),Bicluster)
    for cluster in Bicluster:
        cog_file.append("Reference_clusters"+ " Contig_" + str(cluster[0] ))
        Data=cluster[1:]
        cog_file.append(str(len(Data))+ " proteins" )
        cog_file=cog_file+map(lambda x: str(x).replace("group_","")+"\t+\t?\t"+"Reference_clusters_"+str(x)+"\t?\tReference_clusters_"+str(x)+"\t?" ,Data)
        cog_file.append("")


files=filter(lambda x:len(re.findall(".gff$",x)) ,os.listdir(args.DIR))
Bar_wide=20
pool = Pool(processes=args.processes)
num_tasks=len(files)
if args.processes>num_tasks:
    args.processes=num_tasks
g=pool.imap_unordered(partial(run_on_GFF,args=args,clusters=clusters ), files,int(num_tasks/args.processes))
print "START:"
sys.stdout.write("\r[{}] {:.0f}%".format("#" * 0 + "-" * (Bar_wide - 0),0))
sys.stdout.flush() 
for i, m in enumerate(g, 1):
    progress=i/num_tasks
    Done=int(Bar_wide*progress)
    cog_file.extend(m)
    sys.stdout.write("\r[{}] {:.0f}%".format("#" * Done + "-" * (Bar_wide - Done),round(progress*100,0))) 
    sys.stdout.flush()
pool.close()
print ""
cog_file=map(lambda x: x+"\n",cog_file )
with open(args.out,"wb") as h_file:
    h_file.writelines(cog_file)
    h_file.close()