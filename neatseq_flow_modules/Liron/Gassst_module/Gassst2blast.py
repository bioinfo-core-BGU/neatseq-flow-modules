import os, re 
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Pars GASSST to blast output')

parser.add_argument('-i', type=str,dest='input',
                    help='Input GASSST output File')
parser.add_argument('-o', dest='Output', type=str,default=os.path.join(os.getcwd(),"out.txt"),
                    help='Output file')
args = parser.parse_args()
f=open(args.input,"r")
lines=f.readlines()
f.close()
data=pd.DataFrame()
BANKs=[x for x in lines if x.startswith("BANK")]
QUERYs=[x for x in lines if x.startswith("QUERY")]
BANKs_fields=[[x for x in re.split("\s",y) if x!=''] for y in BANKs]
QUERYs_fields=[[x for x in re.split("\s",y) if x!=''] for y in QUERYs]
gaps=[float(z[0].replace('gap(s):','')) for z in [y for y in [re.findall('gap\(s\):[0-9 ]+',x) for x in lines] if len(y)>0]]
mismatches=[float(z[0].replace('# mismatche(s):','')) for z in [y for y in [re.findall('# mismatche\(s\):[0-9 ]+',x) for x in lines] if len(y)>0]]
data["qseqid"]=[x[-1] for x in QUERYs_fields]
data["sseqid"]=[x[-1] for x in BANKs_fields]
data["qlen"]=[abs(int(x[3]) -int(x[1]) )+1 for x in QUERYs_fields]
data["slen"]=[abs(int(x[3]) -int(x[1]) )+1 for x in BANKs_fields]
data["qstart"]=[int(x[1]) for x in QUERYs_fields]
data["qend"]=[int(x[3]) for x in QUERYs_fields]
data["sstart"]=[int(x[1]) for x in BANKs_fields]
data["send"]=[int(x[3]) for x in BANKs_fields]
data["length"]=[len(x[2]) for x in BANKs_fields]                  
data["evalue"]=[float(z[0].replace('e-value:','')) for z in [y for y in [re.findall('e-value:[0-9 e \- \+ \.]+',x) for x in lines] if len(y)>0]]
data["pident"]=list(map(lambda gap,mis,alen: 100*(1-((gap+mis)/alen)),gaps,mismatches,data["length"]))
data["sseq"]=[x[2].upper() for x in BANKs_fields]   
data.to_csv(args.Output, sep='\t',header=False,index=False,float_format="%g")