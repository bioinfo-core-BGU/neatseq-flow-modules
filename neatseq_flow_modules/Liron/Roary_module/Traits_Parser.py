import os, re
import argparse
import pandas as pd


parser = argparse.ArgumentParser(description='Traits Parser')
parser.add_argument('-M', type=str,
                    help='MetaData file')
parser.add_argument('-P', type=str,
                    help='Optional gene presence/absence file to use only shared samples')
parser.add_argument('-O' , type=str, default=os.getcwd(),
                    help='Output file directory')
parser.add_argument('--S_MetaData', type=str, default="Sample",
                    help='samples ID field in the metadata file')
parser.add_argument('--Fields_val', nargs='+', type=str, default=[],
                    help=" Pairs of field and operator + value to convert to boolean traits: field_name1/op_value1 .. field_nameN/op_valueN \n ,\
                    Example: ""field_1/>=val_1<val_2""    ""feild_2/=='str_val'"" ")
args = parser.parse_args()

op={}
op['==']="eq"
op['<']="lt"
op['<=']="le"
op['=<']="le"
op['!=']="ne"
op['=!']="ne"
op['>']="gt"
op['>=']="ge"
op['=>']="ge"
op['&']="and"
op['|']="or"

MetaData = pd.read_csv(args.M , sep='\t',index_col=False)
MetaData.loc[:,args.S_MetaData]=map(lambda x:str(x),MetaData.loc[:,args.S_MetaData])
temp_MetaData=MetaData[[args.S_MetaData]].copy()
temp_MetaData=temp_MetaData.set_index(args.S_MetaData).copy()
print args.Fields_val
MetaData=MetaData.set_index(args.S_MetaData).copy()
Final_Fields=[]
for Field_val in args.Fields_val:
    if len(Field_val.split("&"))==2:
        full_Field_val=Field_val
        Filter_val,Field_val=Field_val.split("&")
        if len(Field_val.split("/"))==2:
            field,val=Field_val.split("/")
            if field in MetaData.columns:
                temp_MetaData["Temp"]=MetaData[field]
                Field_val=reduce(lambda x, y: x.replace(y, op[y]), op, full_Field_val.replace("/",'_'))
                MetaData[Field_val]=0
                MetaData.loc[temp_MetaData.eval("Temp"+val.replace('/','')),Field_val]=1
                if len(Filter_val.split("/"))==2:
                    Filter,val_Filter=Filter_val.split("/")
                    if Filter in MetaData.columns:
                        temp_MetaData["Temp"]=MetaData[Filter]
                        MetaData.loc[temp_MetaData.eval("Temp"+val_Filter.replace('/','')),Field_val]="NA"
                Final_Fields.append(Field_val)
    else:
        if len(Field_val.split("/"))==2:
            field,val=Field_val.split("/")
            if field in MetaData.columns:
                temp_MetaData["Temp"]=MetaData[field]
                Field_val=reduce(lambda x, y: x.replace(y, op[y]), op, Field_val.replace("/",'_'))
                MetaData[Field_val]=0
                MetaData.loc[temp_MetaData.eval("Temp"+val.replace('/','')),Field_val]=1
                Final_Fields.append(Field_val)
new_MetaData=MetaData[Final_Fields].copy() 
if args.P is not None:
    genes_file= pd.read_csv(args.P, sep=',',index_col=False, low_memory=False)
    if "Inference" in genes_file.columns:
        genes_file.drop(["Inference"],axis=1,inplace=True)
    #genes_file.rename(columns=lambda x: filter(lambda y: y in x ,new_MetaData.index)[0] if len(filter(lambda y: y in x ,new_MetaData.index))>0 else x, inplace=True)    
    #col=map(lambda x: filter(lambda y: x.startswith(y+'_') ,new_MetaData.index) if  len(filter(lambda y: x.startswith(y+'_') ,new_MetaData.index))==1 else x  ,genes_file.columns)
    genes_file.rename(columns=lambda x: filter(lambda y: x.startswith(str(y)+'_') ,new_MetaData.index)[0] if  len(filter(lambda y: x.startswith(str(y)+'_') ,new_MetaData.index))==1 else x  , inplace=True)
    samples_names=genes_file.columns[14:]
    shared_samples= set(samples_names) & set(new_MetaData.index)
    new_MetaData=new_MetaData.ix[shared_samples].copy()
    temp_col=list(genes_file.columns[:14])
    temp_col.extend(list(shared_samples))
    new_genes_file=genes_file[temp_col].copy()
    new_genes_file.to_csv(os.path.join(args.O, args.P.split(os.sep)[-1]), sep=',',float_format="%g",index=False)
       
new_MetaData=new_MetaData.rename_axis('').copy()
new_MetaData.to_csv(os.path.join(args.O,'Traits_file.csv'), sep=',',index=True,float_format="%g")