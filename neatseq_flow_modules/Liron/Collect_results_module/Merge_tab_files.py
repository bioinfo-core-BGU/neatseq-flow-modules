import os, re 
import argparse
import pandas as pd
STRING_TYPES = (basestring, unicode, bytes)
parser = argparse.ArgumentParser(description='Merge tabular files')
parser.add_argument('-D', type=str,dest='directory', nargs='+',
                    help='Location to search')
parser.add_argument('-R', dest='Regular', type=str,
                    help='Regular Expression')
parser.add_argument('-O', dest='Output', type=str,default=os.path.join(os.getcwd(), 'Merged_files.tab'),
                    help='Output file directory')
parser.add_argument('--add_samples_names', dest='samples_names', action='store_true',default=False,
                    help='infer and add samples names from file directory to "Samples" column')
parser.add_argument('--pivot', dest='pivot', nargs='+', type=str, default=[],
                    help='Convert to pivot table by [index columns values]')
parser.add_argument('--Merge_by', type=str,dest='merge_by', nargs='+', default=None,
                    help='Merge files by common column (append is the default)')
parser.add_argument('--MetaData', type=str,dest='MetaData',
                    help='Use external MetaData file as the base for merging')
parser.add_argument('--header', dest='header', action='store_true',default=False,
                    help='dont use a header row, use integers')
parser.add_argument('--split_by', type=str,dest='split_by',
                    help='Split columns before pivot')
parser.add_argument('--Excel', dest='Excel', action='store_true',default=False,
                    help='Collect all results to one Excel file')
parser.add_argument('--sep', dest='sep', type=str,default="\t",
                    help='Columns separator for input file')
parser.add_argument('-T', dest='Trans', action='store_true',default=False,
                    help='write Transpose output')
parser.add_argument('--ignore_shared_col', dest='ignore', action='store_true',default=False,
                    help='Ignore shared columns when merging using the --Merge_by option')
args = parser.parse_args()

if args.header:
    args.pivot=map(lambda x:int(x)  if str.isdigit(x) else x ,args.pivot)
    header=None
else:
    header='infer'

directory=args.directory
Query=args.Regular
Output=args.Output

def find_files(directory, Query):
    file_name=list()
    #print("Searching..." + Query)
    for root, dirs, files in os.walk(directory):
        for name in files:
            if len(re.findall(Query,name))>0:
                file_name.append(os.path.join(root, name))
                #print name
                
    return file_name
index=False
files=[]
for Dir in directory:
    files += find_files(Dir, Query)

if args.Excel:
    if len(files)>0:
        print str(len(files))+" files were found"
        if os.path.exists(Output):
            from openpyxl import load_workbook
            book = load_workbook(Output)
            writer = pd.ExcelWriter(Output, engine='openpyxl')
            writer.book = book
            writer.sheets = dict((ws.title, ws) for ws in book.worksheets)
        else:
            writer = pd.ExcelWriter(Output)
        flag=0
        for file_name in files:
            print file_name
            if os.stat(file_name).st_size > 0:
                temp_data = pd.read_table(file_name, sep=args.sep,header=header,low_memory=False)
                if temp_data.shape[0]==0:
                    print file_name +" is empty!!!!"
                else:
                    temp_data=temp_data.applymap(lambda x:unicode(re.sub('[\000-\010]|[\013-\014]|[\016-\037]', "",x), errors='ignore') if isinstance(x, STRING_TYPES) else x)
                    temp_data.to_excel(writer, sheet_name=re.sub(Query,"",os.path.split(file_name)[-1]), engine='openpyxl')
            else:
                print file_name +" is empty!!!!"
        writer.save()

else:
    if args.MetaData!=None:
        files=files+[args.MetaData]
    if len(files)>0:
        print str(len(files))+" files were found"
        flag=0
        for file_name in files:
            if os.stat(file_name).st_size > 0:
                print file_name
                temp_data = pd.read_table(file_name, sep=args.sep,header=header)
                if temp_data.shape[0]==0:
                    print file_name +" is empty!!!!"
                if args.samples_names:
                    temp_data["Samples"]=os.path.basename(os.path.dirname(file_name))
                if args.merge_by!=None:
                    for merge_by in args.merge_by:
                        if merge_by in temp_data.columns:
                            temp_data.rename(columns=lambda x: args.merge_by[0] if  x==merge_by else x, inplace=True)
                
                if flag==0:
                    Data=temp_data.copy()
                    flag=1
                else:
                    try:
                        if args.merge_by!=None:
                            if args.MetaData==file_name:
                                temp_temp_data=Data.copy()
                                Data=temp_data.copy()
                                temp_data=temp_temp_data.copy()
                            if args.ignore:
                                Data=Data.merge(temp_data, on=args.merge_by[0] ,how='outer',suffixes=('_REMOVE_X', '_REMOVE_Y'))
                                
                                Data=Data[filter(lambda x: ((str(x).endswith('_REMOVE_X'))|(str(x).endswith('_REMOVE_Y')))==False, Data.columns )].copy()
                            else:
                                Data=Data.merge(temp_data, on=args.merge_by[0] ,how='outer',suffixes=('', "["+re.sub(Query,"",os.path.split(file_name)[-1]) +"]"))
                        else:
                            Data=Data.merge(temp_data, on=None ,how='outer',suffixes=('', "["+re.sub(Query,"",os.path.split(file_name)[-1]) +"]"))
                    except:
                         print "File %s could not be merged" %file_name
                    # if args.merge_by!=None:
                        # index=True
            else:
                print file_name +" is empty!!!!"
        if len(args.pivot)==3:
            if args.split_by!=None:
                Data=(Data.drop(args.pivot[1], axis=1)
                .join(
                Data[args.pivot[1]]
                .str
                .split(args.split_by,expand=True)
                .stack()
                .str
                .strip()
                .reset_index(drop=True, level=1)
                .rename(args.pivot[1])
                 ))#.reset_index(drop=True, level=0)
            index=True
            if len(Data)>0:
                Data=Data.groupby([args.pivot[0],args.pivot[1]])[args.pivot[2]].apply(list).reset_index()
                Data[args.pivot[2]]=map(lambda x: str(x).replace("'",'').replace('"','').replace("[ ","").replace("]","").replace("[","").replace(" "," "),Data[args.pivot[2]])
                Data=Data.pivot(index=args.pivot[0], columns=args.pivot[1], values=args.pivot[2]).copy()
                Data.columns=map(lambda x:x.strip(" "),Data.columns)
                if args.Trans:
                    Data.T.to_csv(Output ,sep='\t',index=index,float_format="%g")
                else:
                    Data.to_csv(Output ,sep='\t',index=index,float_format="%g")
        else:
            Data.columns=map(lambda x:x.strip(" "),Data.columns)
            if args.Trans:
                Data.T.to_csv(Output ,sep='\t',index=index,float_format="%g")
            else:
                Data.to_csv(Output ,sep='\t',index=index,float_format="%g")