from datetime import datetime,timedelta
from hashlib import blake2b
from pathlib import Path

import numpy as np
import pandas as pd

#mrn_new column is incorrect, mrn column is correct

#Read in clinic table, pad MRN to 8 digits and remove MRNs that have more/less than 8
df=pd.read_csv("~/Clinic-Data/clinicFiles/MRN_new.txt",sep='\t') #correct version
df['mrn']=df['mrn'].apply('{:0>8}'.format)
df=df[df['mrn'].apply(len)==8]
df['mrn_new']=df['mrn_new'].apply('{:0>8}'.format)
df=df[df['mrn_new'].apply(len)==8]

#Define Hash function
def hash_data(data):
  hasher=blake2b(digest_size=8)
  hasher.update(data.encode())
  return(hasher.hexdigest())
  
#Hash MRN
def generate_synthetic_mrn(row):
  mrn=str(row['mrn'])
  ######Set new seed#####
  mrn_hash=hash_data(mrn)
  site_code='UCSD'
  row['Patient_ID_wrong']=site_code+'P'+mrn_hash
  return(row)
  
#Hash MRN_new
def generate_synthetic_mrn_new(row):
  mrn=str(row['mrn_new'])
  ######Set new seed#####
  mrn_hash=hash_data(mrn)
  site_code='UCSD'
  row['Patient_ID_real']=site_code+'P'+mrn_hash
  return(row)
  
func1=lambda row: generate_synthetic_mrn(row)
func1.__name__='Hashmrn'

func2=lambda row: generate_synthetic_mrn_new(row)
func2.__name__='Hashmrn_new'

#apply hash function, re order columns
hashdf=df.apply(func1,axis=1).apply(func2,axis=1)
#hashdf=hashdf.drop(columns=['MRN'])
#cols=list(hashdf)
#cols=[cols[-1]]+cols[:-1]
#hashdf=hashdf[cols]

#write hashed file
hashdf.to_csv('~/Clinic-Data/clinicFiles/20210829_hashed/hashed_mrn_new.txt',sep='\t',index=False) #correct version
