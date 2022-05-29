from datetime import date,datetime,timedelta
from hashlib import blake2b
from pathlib import Path

#import pdb; pdb.set_trace()
import numpy as np
import pandas as pd
import sys

#Ask for input file, read
inputfile=str(sys.argv[1])
df=pd.read_csv(inputfile,sep='\t')

#Print column names
for col in df.columns:
	print(col)
	
#Define MRN column
mrnCol=input("Enter name of MRN Column: ")
	
#Define Report ID column
reportCol=input("Enter name of Report ID Column: ")

#Define Columns to be hashed
columns_to_hash=[]
columns_to_hash=[item for item in input("Enter names of columns to be hashed: ").split()]

#Drop blank MRNs
if mrnCol in df.columns:
	df[mrnCol]=df[mrnCol].apply('{:0>8}'.format)
	df=df[df[mrnCol].apply(len)==8]
	df.dropna(subset=[mrnCol],inplace=True)

#Drop blank Report IDs
if reportCol in df.columns:
	df.dropna(subset=[str(reportCol)],inplace=True)
	
#Define hashing function
def hash_data(data):
  hasher=blake2b(digest_size=8)
  hasher.update(data.encode())
  return(hasher.hexdigest())

#Define row hashing with specific identifiers
def hash_column(row):
  original=str(row[column])
  ######Set new seed#####
  hashed_original=hash_data(original)
  site_code='UCSD'
  if column==mrnCol:
  	identifier='P'
  elif column==reportCol:
  	identifier='S'
  else:
  	identifier='X'
  row[column+'_Hashed']=site_code+identifier+hashed_original
  return(row)
	
func1=lambda row: hash_column(row)
func1.__name__='HashColumn'

#Apply hashing function
for column in columns_to_hash:
	if column in df.columns:
		df=df.apply(func1,axis=1)
	else:
		print(column+" not in data,skipping")
	df.fillna('',inplace=True)

#Create file save name
saveName=input("Enter save file name: ")
saveDest=input("Enter save destination: ")
filename=saveDest+saveName+".txt"

#Write file
df.to_csv(filename,sep='\t',index=False)
print('Hashing complete')


