import os, re 
import sys
import pandas as pd

files= sys.argv[1:]
flag=0
for file_name in files:
	temp_data = pd.read_csv(file_name, sep='\t',index_col=0)
	sample=file_name.split(os.sep)[-1]
	path=file_name.rstrip(sample).rstrip(os.sep)
	if sample.endswith('.genes.results'):
		sample=sample.rstrip('.genes.results')
		prefix='GeneMat_'
	else:
		sample=sample.rstrip('.isoforms.results')
		prefix='IsoMat_'
	if flag==0:
		Data_counts=temp_data["expected_count"].copy()
		Data_TPM=temp_data["TPM"].copy()
		Data_FPKM=temp_data["FPKM"].copy()
		Data_counts.name=sample
		Data_TPM.name=sample
		Data_FPKM.name=sample
		flag=1
	else:
		temp_counts=temp_data["expected_count"].copy()
		temp_TPM=temp_data["TPM"].copy()
		temp_FPKM=temp_data["FPKM"].copy()
		temp_counts.name=sample
		temp_TPM.name=sample
		temp_FPKM.name=sample
		Data_counts=pd.concat([Data_counts,temp_counts],axis=1)
		Data_TPM=pd.concat([Data_TPM,temp_TPM],axis=1)
		Data_FPKM=pd.concat([Data_FPKM,temp_FPKM],axis=1)
if flag!=0:
	with pd.option_context('display.max_rows', None, 'display.max_columns', None):
		print(Data_counts.to_csv(sep="\t",index=True))
	if len(path)>0:
		Data_TPM.to_csv(os.sep.join([path,prefix+"TPM"]) ,sep='\t',index=True)
		Data_FPKM.to_csv(os.sep.join([path,prefix+"FPKM"]) ,sep='\t',index=True)
	else:
		Data_TPM.to_csv(prefix+"TPM" ,sep='\t',index=True)
		Data_FPKM.to_csv(prefix+"FPKM" ,sep='\t',index=True)