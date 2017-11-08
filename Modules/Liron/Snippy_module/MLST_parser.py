import os, re 
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Pars MLST')
parser.add_argument('-M', type=str,
                    help='MetaData file')
parser.add_argument('-F',  type=str,
                    help='Merged MLST typing file')
parser.add_argument('-O' , type=str, default=os.getcwd(),
                    help='Output file directory')
parser.add_argument('-C', type=float, default=0.95,
                    help='Percentage of identified allele cutoff to consider sample [0.0 - 1.0]')
parser.add_argument('--S_MetaData', type=str, default="Sample",
                    help='samples ID field in the metadata file')
parser.add_argument('--S_Merged', type=str, default="Sample",
                    help='samples ID field in the Merged file')
parser.add_argument('--Non_allelic', nargs='+', type=str, default=["Sample",'Status','Percentage_of_missing_genes'],
                    help='Non allelic fields in the Merged file')
parser.add_argument('--Fields', nargs='+', type=str, default=['Status','Percentage_of_missing_genes'],
                    help='Fields to move to the metadata file')
parser.add_argument('--Cut', action='store_true', default=False,
                    help='Use only samples with metadata information')
parser.add_argument('--FASTA', action='store_true', default=False,
                    help='The input is a FASTA file')
args = parser.parse_args()
Fields=[]
if args.Fields != None:
    for field in args.Fields:
        if len(field.split(","))>1:
            for field_t in field.split(","):
                Fields=Fields+[field_t]
        else:
            Fields=Fields+[field]
    args.Fields=Fields

Fields=[]
if args.Non_allelic != None:
    for field in args.Non_allelic:
        if len(field.split(","))>1:
            for field_t in field.split(","):
                Fields=Fields+[field_t]
        else:
            Fields=Fields+[field]
    args.Non_allelic=Fields

flag=0
if args.FASTA:
    from Bio import AlignIO
    msa=AlignIO.read(args.F,"fasta")
    data=pd.DataFrame.from_records(msa)
    data.index=map(lambda x:msa[int(x)].id,list(data.index))
    data=data.drop(data.columns[data.apply(lambda x: len(re.findall("[- N]",x.sum().upper()))>0 ,axis=0)],axis=1)
    temp_data=data
else:
    temp_data = pd.read_csv(args.F, sep='\t',index_col=False)
    temp_data=temp_data.set_index(args.S_Merged,drop=False).copy()

for j in temp_data.index:
    for i in temp_data.columns:
        if type(temp_data.loc[j,i])==str:
            if temp_data.loc[j,i].startswith("New_Allele="):
                temp_data.loc[temp_data.loc[:,i]==temp_data.loc[j,i],i]=i+"_"+str(j)



temp_data.index=map(lambda x:str(x),temp_data.index)
if (args.M != None)&(args.Cut):
    MetaData = pd.read_csv(args.M , sep='\t',index_col=False)
    MetaData=MetaData.set_index(args.S_MetaData,drop=False).copy()
    MetaData.index=map(lambda x:str(x),MetaData.index)
    flag=1
    temp_data=temp_data.loc[map(lambda x:x in MetaData.index,temp_data.index),].copy()
args.Non_allelic.extend([args.S_Merged])
args.Non_allelic.extend(args.Fields)
if None in args.Non_allelic:
    args.Non_allelic.remove(None)
args.Non_allelic=set(args.Non_allelic)

def cut_rows(temp_data,cutoff,Non_allelic_rows):
    drop=[x for x in temp_data.columns if x not in Non_allelic_rows]
    temp_data=temp_data[drop]
    stay=list()
    for row in temp_data.index:
        if (float(temp_data.ix[row].count())/ float(temp_data.shape[1]))>=cutoff:
            stay.append(row)
        else:
            print "The Sample %s has lower percentage of identified allele (%%s) than the cutoff" % row % (float(temp_data.ix[row].count())/ float(temp_data.shape[1]))
    return  temp_data.ix[stay].copy()

def cut_col(temp_data,Non_allelic):
    stay=list()
    for col in temp_data.columns:
        if col not in Non_allelic:
            if temp_data[col].count()==temp_data.shape[0]:
                stay.append(col)
        else:
            stay.append(col)
    return temp_data[stay].copy()

def drop(data,fields,op=1):
    if op==1:
        for i in fields:
            if i in data.columns:
                data=data.drop(i,axis=1).copy()
    else:
        for i in fields:
            if i in data:
                data=data.drop(i).copy()
    return data

new_temp_data=cut_col(cut_rows(temp_data, args.C ,args.Non_allelic),args.Non_allelic)
m=new_temp_data.columns
m=drop(m,args.Non_allelic,0).copy()
g=new_temp_data.groupby(list(m[:]))
new_temp_data["Index"]=''
num=1
for i in g:
    new_temp_data.loc[i[1].index,"Index"]=str(num)
    num=num+1

if args.M != None:
    if flag!=1:
        MetaData = pd.read_csv(args.M , sep='\t',index_col=False)
        MetaData=MetaData.set_index(args.S_MetaData,drop=False).copy()
        MetaData.index=map(lambda x:str(x),MetaData.index)
    MetaData=MetaData.join(new_temp_data["Index"])
    MetaData=MetaData.loc[~MetaData["Index"].isnull(),:]
    if args.Fields != None:
        for field in args.Fields:
            if field in temp_data.columns:
                MetaData=MetaData.join(temp_data[field],lsuffix='_Old')
    MetaData=MetaData.set_index("Index").copy()
else:
    MetaData=new_temp_data[["Index"]].copy()
    if args.Fields != None:
        for field in args.Fields:
            if field in temp_data.columns:
                MetaData=MetaData.join(temp_data[field],lsuffix='_Old')

MetaData.to_csv(os.path.join(args.O,'New_MetaData.tab'), sep='\t',float_format="%g",index=True) 
new_temp_data=new_temp_data.set_index("Index").copy()

#new_temp_data=drop(new_temp_data,args.Non_allelic).copy()
new_temp_data.to_csv(os.path.join(args.O, 'New_Merged_cut.tab'), sep='\t',index=True,float_format="%g")

