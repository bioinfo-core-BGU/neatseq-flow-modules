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
BANKs=filter(lambda x:x.startswith("BANK"),lines)
QUERYs=filter(lambda x:x.startswith("QUERY"),lines)
BANKs_fields=map(lambda y: filter(lambda x:x!='' ,re.split("\s",y) ) ,BANKs)
QUERYs_fields=map(lambda y:filter(lambda x:x!='' ,re.split("\s",y) ),QUERYs)
gaps=map(lambda z: float(z[0].replace('gap(s):','')),filter(lambda y: len(y)>0 ,map(lambda x: re.findall('gap\(s\):[0-9 ]+',x),lines )))
mismatches=map(lambda z: float(z[0].replace('# mismatche(s):','')),filter(lambda y: len(y)>0 ,map(lambda x: re.findall('# mismatche\(s\):[0-9 ]+',x),lines )))
data["qseqid"]=map(lambda x: x[-1], QUERYs_fields)
data["sseqid"]=map(lambda x: x[-1], BANKs_fields)
data["qlen"]=map(lambda x: abs(int(x[3]) -int(x[1]) )+1, QUERYs_fields)
data["slen"]=map(lambda x: abs(int(x[3]) -int(x[1]) )+1, BANKs_fields)
data["qstart"]=map(lambda x: int(x[1]), QUERYs_fields)
data["qend"]=map(lambda x: int(x[3]), QUERYs_fields)
data["sstart"]=map(lambda x: int(x[1]), BANKs_fields)
data["send"]=map(lambda x: int(x[3]), BANKs_fields)
data["length"]=map(lambda x: len(x[2]), BANKs_fields)                  
data["evalue"]=map(lambda z: float(z[0].replace('e-value:','')),filter(lambda y: len(y)>0 ,map(lambda x: re.findall('e-value:[0-9 e \- \+ \.]+',x),lines )))
data["pident"]=map(lambda gap,mis,alen: 100*(1-((gap+mis)/alen)),gaps,mismatches,data["length"])
data["sseq"]=map(lambda x: x[2].upper(), BANKs_fields)   
data.to_csv(args.Output, sep='\t',header=False,index=False,float_format="%g")