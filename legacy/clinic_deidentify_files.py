from datetime import datetime,timedelta
from hashlib import blake2b
from pathlib import Path

import numpy as np
import pandas as pd

#Upper limit-200000 lines table

#Read in clinic table, pad MRN to 8 digits and remove MRNs that have more/less than 8
df=pd.read_csv("/home/scp010/Clinic-Data/clinicFiles/20210829/xax.txt",sep='\t') #correct version
df['MRN']=df['MRN'].apply('{:0>8}'.format)
df=df[df['MRN'].apply(len)==8]

#Define Hash function
def hash_data(data):
  hasher=blake2b(digest_size=8)
  hasher.update(data.encode())
  return(hasher.hexdigest())
  
#Hash MRN
def generate_synthetic_mrn(row):
  mrn=str(row['MRN'])
  ######Set new seed#####
  mrn_hash=hash_data(mrn)
  site_code='UCSD'
  row['Patient_ID']=site_code+'P'+mrn_hash
  return(row)
  
func1=lambda row: generate_synthetic_mrn(row)
func1.__name__='HashMRN'

#apply hash function, re order columns
hashdf=df.apply(func1,axis=1)
hashdf=hashdf.drop(columns=['MRN'])
cols=list(hashdf)
cols=[cols[-1]]+cols[:-1]
hashdf=hashdf[cols]

#write hashed file
hashdf.to_csv('~/Clinic-Data/clinicFiles/20210829_hashed/xax.txt',sep='\t',index=False) #correct version
