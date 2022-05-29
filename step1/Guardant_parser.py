import pandas as pd
import json
import os
from os import path

inputDir='/home/scp010/Guardant/retrospective/'

def readFile(file):
	with open(file) as f:
		data=json.load(f)
	
	mrn=data['patient']['mrn']
	accessionId=data['accessionId']
	alteration=data['alterationDetected']
	
	df=pd.DataFrame(columns=('MRN','Accession_ID','Alteration'))

	df=df.append({'MRN':mrn,'Accession_ID':accessionId,'Alteration':alteration},ignore_index=True)
	
	if path.exists('Patients.txt'):
		df.to_csv('/home/scp010/Guardant/Patients.txt',header=False,mode='a+',sep='\t',index=False)
	else:
		df.to_csv('/home/scp010/Guardant/Patients.txt',header=True,sep='\t',index=False)
		
for file in os.listdir(inputDir):
  if file.endswith(".json") and "CANCELLED" not in file:
  	line=str(os.path.join(inputDir,file))
  	try:
  		readFile(line)
  		
  	except:
  		print(str(file))

