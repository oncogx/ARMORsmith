import pandas as pd
import xml.etree.ElementTree as ET
import os
from os import path
#import logging
#import pdb; pdb.set_trace()

#logging.basicConfig(filename='error.log',level=logging.DEBUG)

inputdir='/home/scp010/Caris/xml/'

def readMRN(file):
  tree=ET.parse(file)
  root=tree.getroot()
  df_pat=pd.DataFrame(columns=('MRN','Caris_Report_ID','Gender','Test_Type','Order_Date','Received_Date','DOB','ICD','Diagnosis','Pathological_Diagnosis','Site','Lineage','Sublineage','Organization','Sample','Specimen_Type','Specimen_Accession','Specimen_Site','Specimen_Date'))
  df_var=pd.DataFrame()
  
  columns=['biomarkerName','gene','hgvsCodingChange','hgvsProteinChange','chromosome','genomeBuild','exon','genomicSource','alterationDetails','molecularConsequence','alleleFrequencyInformation','mutationBurdenCall','mutationBurdenScore','mutationBurdenUnit','msiCall','genomicCoordinates','copyNumberType','copyNumber','fusionISOForm','gene1','exon1','transcriptID1','gene2','exon2','transcriptID2','genomicBreakpoint','LOHpercentage','readInformation','platformTechnology']
  genomicAlteration=['referenceNucleotide','observedNucleotide','transcriptStartPosition','transcriptStopPosition','transcriptID','transcriptIDsource']
  labLabel=['analysisConfigurationName','analysisConfigurationVersion','analysisPipelineName','analysisPipelineVersion','NGSPanelName','NGSPanelVersion']
  
  mrn=root.findtext('.//mrn')
  if mrn!="":
    reportID=root.findtext('.//labReportID')
    orderDate=root.findtext('.//orderedDate')
    receiveDate=root.findtext('.//receivedDate')
    dob=root.findtext('.//dob')
    gender=root.findtext('.//gender')
    icd=root.findtext('.//icd_code')
    diagnosis=root.findtext('.//diagnosis')
    pathologic=root.findtext('.//pathologicDiagnosis')
    site=root.findtext('.//primarySite')
    lineage=root.findtext('.//lineage')
    sublineage=root.findtext('.//subLineage')
    organization=root.findtext('.//healthcareOrganization/name')
    specimenID=root.findtext('.//specimenInformation/tumorSpecimenInformation/specimenID')
    specimenType=root.findtext('.//specimenInformation/tumorSpecimenInformation/specimenType')
    specimenAccession=root.findtext('.//specimenInformation/tumorSpecimenInformation/specimenAccessionID')
    specimenSite=root.findtext('.//specimenInformation/tumorSpecimenInformation/specimenSite')
    specimenDate=root.findtext('.//specimenInformation/tumorSpecimenInformation/specimenCollectionDate')
    for child in root.findall('.//tests'):
      testType=child.findtext('testMethodology')
      df_pat=df_pat.append({'MRN':mrn,'Caris_Report_ID':reportID,'Gender':gender,'Test_Type':testType,'Order_Date':orderDate,'Received_Date':receiveDate,'DOB':dob,'ICD':icd,'Diagnosis':diagnosis,'Pathological_Diagnosis':pathologic,'Site':site,'Lineage':lineage,'Sublineage':sublineage,'Organization':organization,'Sample':specimenID,'Specimen_Type':specimenType,'Specimen_Accession':specimenAccession,'Specimen_Site':specimenSite,'Specimen_Date':specimenDate},ignore_index=True)
      if testType=="NGS":
        for ngs_category in child.findall('.//testResults/'):
          for ngs_mutations in child.findall('.//testResults/'+ngs_category.tag):
            ngs_result=ngs_mutations.findtext('result_group')
            if ngs_result!="Normal" and ngs_result!="No Result":
              ngs_dictionary={}
              for column in columns:
                ngs_dictionary[column]=ngs_mutations.findtext(column)
                if ngs_mutations.findtext('alterationDetails') is not None:
                  for ngs_alterations in genomicAlteration:
                    ngs_dictionary[ngs_alterations]=ngs_mutations.findtext('./alterationDetails/transcriptAlterationDetails/'+ngs_alterations)
                if ngs_mutations.findtext('alleleFrequencyInformation') is not None:
                  ngs_dictionary['alleleFrequency']=ngs_mutations.findtext('./alleleFrequencyInformation/alleleFrequency')
                if ngs_mutations.findtext('readInformation') is not None:
                  ngs_dictionary['readDepth']=ngs_mutations.findtext('./readInformation/readDepth')
              for ngs_label in labLabel:
                ngs_dictionary[ngs_label]=ngs_mutations.findtext('./labSpecific/'+ngs_label)
              ngs_dictionary.update(MRN=root.findtext('.//mrn'),Caris_Report_ID=reportID,Test_Name=child.findtext('testName'),Test_Code=child.findtext('testCode'),Platform_Technology=child.findtext('platformTechnology'),Test_Type=testType,Category=ngs_category.tag,Result=ngs_result)
              for ngs_key in ngs_dictionary:
                if ngs_dictionary[ngs_key] is None:
                  ngs_dictionary[ngs_key]=''
              df_var=df_var.append(ngs_dictionary,ignore_index=True)
              df_var=df_var.fillna('')
              df_var=df_var.drop_duplicates()
      if testType=="CNA-NGS" or testType=="CNA-Seq" or testType=="RNA-Seq" or testType=="Seq":
        for cna_category in child.findall('.//testResults/'):
          for cna_mutations in child.findall('.//testResults/'+cna_category.tag):
            cna_result=cna_mutations.findtext('result_group')
            if cna_result!="Normal" and cna_result!="No Result":
              cna_dictionary={}
              for column in columns:
                cna_dictionary[column]=cna_mutations.findtext(column)
                if cna_mutations.findtext('alterationDetails') is not None:
                  for cna_alterations in genomicAlteration:
                    cna_dictionary[cna_alterations]=cna_mutations.findtext('./alterationDetails/transcriptAlterationDetails/'+cna_alterations)
                if cna_mutations.findtext('alleleFrequencyInformation') is not None:
                  cna_dictionary['alleleFrequency']=cna_mutations.findtext('./alleleFrequencyInformation/alleleFrequency')
                if cna_mutations.findtext('readInformation') is not None:
                  cna_dictionary['readDepth']=cna_mutations.findtext('./readInformation/readDepth')
                for cna_label in labLabel:
                  cna_dictionary[cna_label]=cna_mutations.findtext('./labSpecific/'+cna_label)
              cna_dictionary.update(MRN=root.findtext('.//mrn'),Caris_Report_ID=reportID,Test_Name=child.findtext('testName'),Test_Type=testType,Category=cna_category.tag,Result=cna_result)
              for cna_key in cna_dictionary:
                if cna_dictionary[cna_key] is None:
                  cna_dictionary[cna_key]=''
              df_var=df_var.append(cna_dictionary,ignore_index=True)
              df_var=df_var.fillna('')
              df_var=df_var.drop_duplicates()
  if path.exists('Patients.txt'):
    df_pat.to_csv('/home/scp010/Caris/Patients.txt',header=False,mode='a+',sep='\t',index=False)
  else:
    df_pat.to_csv('/home/scp010/Caris/Patients.txt',header=True,sep='\t',index=False)
  
  if not df_var.empty:
    if path.exists('Mutations.txt'):
      df_var_orig=pd.read_csv('/home/scp010/Caris/Mutations.txt',sep="\t")
      df_var_orig.append(df_var).to_csv('/home/scp010/Caris/Mutations.txt',sep='\t',header=True,index=False)
    else:
      df_var.to_csv('/home/scp010/Caris/Mutations.txt',header=True,sep='\t',index=False)
      
  #if path.exists('Mutations'):
    #df_var.to_csv('/home/scp010/Caris/Mutations',header=False,mode='a+',sep='\t',index=False)
    #df_var_orig=pd.read_csv('/home/scp010/Caris/Mutations',sep="\t")
    #df_var_orig.append(df_var).to_csv('/home/scp010/Caris/Mutations',sep='\t',header=True,index=False)
  #else:
    #df_var.to_csv('/home/scp010/Caris/Mutations',header=True,sep='\t',index=False)
    
    
for file in os.listdir(inputdir):
  if file.endswith(".xml"):
    line = str(os.path.join(inputdir, file))
    try:
      readMRN(line)
      
    except:
      print(str(file))

