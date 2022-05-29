library(plyr)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(splitstackshape)

mutations=read.table("~/JSON/mutations.tsv",sep="\t",colClasses="character")
patients=read.table("~/JSON/patients.tsv",sep="\t",colClasses="character")

mutations[mutations==""]=NA
patients[patients==""]=NA

substrRight=function(x,n){
  substr(x,nchar(x)-n+1,nchar(x))
}

filterDash=function(x){
  if(grepl('-',x)&nchar(x)>1){
    x=x
  } else {
    x=NA
  }
}

mutations=mutations%>%
  unite(DoB,c(DoB,dateOfBirth),na.rm=TRUE,remove=TRUE)%>%
  unite(emr_id,c(emr_id,emrId),na.rm=TRUE,remove=TRUE)%>%
  unite(tempusOrder_id,c(tempusOrder_id_name,tempusOrderId_name),na.rm=TRUE,remove=TRUE)%>%
  #unite(accessionId,c(accessionId,accessionId_name),na.rm=TRUE,remove=TRUE)%>%
  rename(accessionId=accessionId_name)%>%
  select(-c(firstName,lastName,physician_name,tempusOrder_id,tempusId,name_name))%>%
  drop_na(emr_id)%>%
  filter(emr_id!="")%>%
  #rowwise()%>%
  mutate(DoB=if_else(DoB=="",as.character(NA),DoB))%>%
  mutate(DoB=as.Date(DoB))%>%  
  mutate(diagnosisDate=case_when(
    grepl("/",diagnosisDate)~as.Date(diagnosisDate,"%m/%d/%Y"),
    !is.na(diagnosisDate)~lubridate::ymd(diagnosisDate,truncated=2L),
    TRUE~as.Date(NA)
  ))%>%
  group_by(accessionId)%>%
  fill(referenceGenome,.direction="downup")%>%
  unique()

patients=patients%>%
  unite(DoB,c(DoB,dateOfBirth),na.rm=TRUE,remove=TRUE)%>%
  unite(emr_id,c(emr_id,emrId),na.rm=TRUE,remove=TRUE)%>%
  unite(tempusOrder_id,c(tempusOrder_id,tempusOrder_id_name,tempusOrderId_name),na.rm=TRUE,remove=TRUE)%>%
  unite(accessionId,c(accessionId,accessionId_name),na.rm=TRUE,remove=TRUE)%>%
  unite(bioInfPipeline,c(bioInfPipeline,bioInfPipeline_name),na.rm=TRUE,remove=TRUE)%>%
  unite(code,c(code,code_name),na.rm=TRUE,remove=TRUE)%>%
  unite(description,c(description,description_name),na.rm=TRUE,remove=TRUE)%>%
  unite(details,c(details,details_name),na.rm=TRUE,remove=TRUE)%>%
  unite(institution,c(institution,institution_name),na.rm=TRUE,remove=TRUE)%>%
  unite(msiStatus,c(msiStatus,microsatelliteInstability.status),na.rm=TRUE,remove=TRUE)%>%
  unite(name,c(name,name_name),na.rm=TRUE,remove=TRUE)%>%
  unite(notes,c(notes,notes_name),na.rm=TRUE,remove=TRUE)%>%
  unite(physician,c(physician,physician_name),na.rm=TRUE,remove=TRUE)%>%
  unite(reportId,c(reportId,reportId_name),na.rm=TRUE,remove=TRUE)%>%
  unite(reportStatus,c(reportStatus,reportStatus_name),na.rm=TRUE,remove=TRUE)%>%
  unite(reportType,c(reportType,reportType_name),na.rm=TRUE,remove=TRUE)%>%
  unite(signing_pathologist,c(signing_pathologist_name,signingPathologist,signingPathologist_name,signing_pathologist),na.rm=TRUE,remove=TRUE)%>%
  unite(signout_date,c(signoutDate,signout_date,signout_date_name,signoutDate_name),na.rm=TRUE,remove=TRUE)%>%
  #Get BlockID
  unite(blockID,c(institutionData.caseId,institutionData.blockId),sep="-",na.rm=TRUE,remove=TRUE)%>%
  cSplit(.,'blockID',' ',type.convert ="as.character")%>%
  rowwise()%>%
  mutate_at(vars(contains('blockID')),filterDash)%>%
  unite(blockID,contains('blockID'),na.rm=TRUE)%>%
  select(-c(firstName,lastName,tempusId,physician,tempusOrder_id,name,signing_pathologist,bioInfPipeline,notes,reportStatus,reportId,
            details,reportType,microsatelliteInstability.therapies.association,
            microsatelliteInstability.therapies.drugClass,microsatelliteInstability.therapies.evidenceType,
            microsatelliteInstability.therapies.fdaApproved,microsatelliteInstability.therapies.isOnLabel,microsatelliteInstability.therapies.pubMedId,
            microsatelliteInstability.therapies.status,microsatelliteInstability.therapies.therapy,microsatelliteInstability.therapies.tissue,
            microsatelliteInstability.therapies.url))%>%
  drop_na(emr_id)%>%
  filter(emr_id!="")
  #rowwise()%>%

patients[patients==""]=NA

patients=patients%>%
  mutate(DoB=as.Date(DoB))%>%
  mutate(signout_date=as.Date(signout_date))%>%
  mutate(collectionDate=as.Date(collectionDate))%>%
  mutate(receiptDate=as.Date(receiptDate))%>%
  mutate(diagnosisDate=case_when(
    grepl("/",diagnosisDate)~as.Date(diagnosisDate,"%m/%d/%Y"),
    !is.na(diagnosisDate)~lubridate::ymd(diagnosisDate,truncated=2L),
    TRUE~as.Date(NA)
  ))%>%
  filter(sampleCategory!="normal")%>%
  filter(code!="MMR"&code!="PD-L1")%>%
  unique()

expMut=mutations%>%
  mutate(Age_at_Diagnosis=ifelse(!is.na(diagnosisDate),floor(time_length(difftime(diagnosisDate,DoB),"years")),NA))%>%
  dplyr::rename(MRN=emr_id)%>%
  select(-diagnosisDate,-DoB)%>%
  #filter(Age_at_Diagnosis>=18|is.na(Age_at_Diagnosis))%>%
  unique()

expMut$MRN=substrRight(expMut$MRN,8)

expPat=patients%>%
  mutate(Age_at_Report=ifelse(!is.na(signout_date),floor(time_length(difftime(signout_date,DoB),"years")),NA))%>%
  mutate(Days_since_Diagnosis=ifelse(!is.na(diagnosisDate),abs(floor(time_length(difftime(collectionDate,diagnosisDate),"days"))),NA))%>%
  dplyr::rename(MRN=emr_id,Sample=blockID)%>%
  select(-c(diagnosisDate,DoB,collectionDate,receiptDate))%>%
  #filter(Age_at_Report>=18|is.na(Age_at_Report))%>%
  unique()

expPat$MRN=substrRight(expPat$MRN,8)

write.table(expMut,"~/JSON/TempusVariants.tsv",sep='\t')
write.table(expPat,"~/JSON/TempusPatients.tsv",sep='\t')
