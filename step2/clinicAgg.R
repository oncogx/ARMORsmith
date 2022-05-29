library(plyr)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)

#Set file path, extract csv files into list of dataframes, select relevant columns
setwd("~/Clinic-Data/clinicFiles/20210829/") #correct location
filenames=list.files(path="~/Clinic-Data/clinicFiles/20210829/")
filenames=paste("~/Clinic-Data/clinicFiles/20210829/",filenames,sep="")
data=lapply(filenames,read.csv,sep="\t",colClasses="character")

#edit columnsKeep and columnsDrop for relevant columns, check spelling of columns
columnsKeep=c("MRN","person_id","birth_datetime","gender","race","ethnicity","visit_start_datetime",
              "visit_end_datetime","EHR_death_datetime","EDRS_death_datetime","vocabulary_id",
              "condition_source_value","condition_start_datetime")

#use person_id if any files don't have MRN

columnsDrop=c("gender_concept_id","race_concept_id","ethnicity_concept_id","atc_ancestor_vocabulary_id")

data=lapply(data,select,matches(columnsKeep))
data=lapply(data,select,-c(matches(columnsDrop)))
data=data[lapply(data,length)>0]
data=lapply(data,unique)

#loop through list, merge listed dataframes into large dataframe, remove unnecessary rows
for (i in 1:length(data)){
  if (exists("clinic")){
    if (exists("MRN",data[[i]])&exists("MRN",clinic)&exists("person_id",data[[i]])&exists("person_id",clinic)){
      #addition=unique(data[[i]])
      clinic=clinic%>%
        merge(unique(data[[i]]),all=TRUE)%>%
        unique()
    }else if (exists("MRN",data[[i]])&exists("MRN",clinic)){
      #addition=unique(data[[i]])
      clinic=clinic%>%
        merge(unique(data[[i]]),all=TRUE)%>%
        unique()
    }else if (exists("person_id",data[[i]])&exists("person_id",clinic)){
      clinic=clinic%>%
        merge(unique(data[[i]]),all=TRUE)%>%
        unique()
    }else{
      #addition=unique(data[[i]])
      clinic=clinic%>%
        merge(unique(data[[i]]),all=TRUE)%>%
        unique()
    }
  }else{
    clinic=unique(data[[i]])
  }
}

clinic[clinic=="NULL"]=NA
#clinic=clinic[!grepl("affected",clinic$MRN),]
#clinic=clinic[!grepl("affected",clinic$person_id),]
clinic=clinic[!grepl("RxNorm|None",clinic$vocabulary_id),]

#Data cleanup/prep for deidentify, may need to drop some columns
clinic$MRN=str_pad(clinic$MRN,8,pad="0")

#write identifiable version to file, clear environment--file too large to go directly
write.table(clinic,"~/Clinic-Data/clinicFiles/20210829_mid.csv",sep=",")

#Prepare output table, manipulate dates/times,add metastatic
clinic=read.csv("~/Clinic-Data/clinicFiles/20210829_mid.csv",sep=",",colClasses = "character")
oncotree=read.csv("~/ICD-Oncotree",sep=",",colClasses="character")
icd=read.csv("~/ICD10_codes",sep="\t",colClasses="character")

#Trim clinic ICD codes
clinic=clinic%>%
  #Remove spaces from ICD when 2 or more exist
  mutate(condition_source_value=gsub(" ","",condition_source_value,fixed=TRUE))%>%
  #split multiple ICD codes into multiple rows per patient
  mutate(CODE=strsplit(condition_source_value,","))%>%
  unnest(CODE)%>%
  #trim codes for Oncotree match
  mutate(ICD10=substr(.$CODE,1,3))%>%
  #Drop original ICD code columns
  select(-vocabulary_id,-condition_source_value)%>%
  unique()

#Read in oncotree organ level
oncotree=oncotree%>%
  select(ICD10,Oncotree.level.1)%>%
  unique()

#Create clinic file with relevant dates
outClinic=clinic%>%
  group_by(MRN)%>%
  #fill observations
  fill(names(.),.direction="downup")%>%
  ungroup()%>%
  #convert column type to dates
  {if(exists("birth_datetime",.)) mutate(.,birth_datetime=as.Date(.$birth_datetime)) else .}%>%
  {if(exists("condition_start_datetime",.)) mutate(.,condition_start_datetime=as.Date(.$condition_start_datetime)) else .}%>%
  {if(exists("EHR_death_datetime",.)) mutate(.,EHR_death_datetime=as.Date(.$EHR_death_datetime)) else .}%>%
  {if(exists("EDRS_death_datetime",.)) mutate(.,EDRS_death_datetime=as.Date(.$EDRS_death_datetime)) else .}%>%
  {if(exists("visit_start_datetime",.)) mutate(.,visit_start_datetime=as.Date(.$visit_start_datetime)) else .}%>%
  {if(exists("visit_end_datetime",.)) mutate(.,visit_end_datetime=as.Date(.$visit_end_datetime)) else .}%>%
  #select which death date is earliest when having EHR and EDRS columns, with consideration to NA
  {if(exists("EHR_death_datetime",.)&exists("EDRS_death_datetime",.)) 
    mutate(.,death_date=case_when(
      .$EHR_death_datetime<.$EDRS_death_datetime~.$EHR_death_datetime,
      .$EDRS_death_datetime<.$EHR_death_datetime~.$EDRS_death_datetime,
      .$EHR_death_datetime==.$EDRS_death_datetime~.$EDRS_death_datetime,
      is.na(.$EHR_death_datetime)~.$EDRS_death_datetime,
      is.na(.$EDRS_death_datetime)~.$EHR_death_datetime
    )) else .}%>%
  #select death date if only one of EHR and EDRS
  {if(exists("EHR_death_datetime",.)&!exists("EDRS_death_datetime",.)) rename(.,death_date=EHR_death_datetime) else .}%>%
  {if(exists("EDRS_death_datetime",.)&!exists("EHR_death_datetime",.)) rename(.,death_date=EDRS_death_datetime) else .}%>%
  select(-EDRS_death_datetime,-EHR_death_datetime)%>%
  #Merge oncotree
  left_join(oncotree,by="ICD10")%>%
  #Merge metastatic
  left_join(icd,by="CODE")%>%
  #Drop code and ICD10
  select(-CODE,-ICD10)%>%
  #Find oncotree date
  mutate(Oncotree_Date=if_else(!is.na(Oncotree.level.1),condition_start_datetime,as.Date(NA)))%>%
  #Find metastatic date
  mutate(Metastatic_Date=if_else(!is.na(Metastatic),condition_start_datetime,as.Date(NA)))%>%
  select(-condition_start_datetime)%>%
  #Start grouped transform
  group_by(MRN,Oncotree.level.1)%>%
  #Find earliest diagnosis date for each oncotree code
  mutate(Oncotree_Date=if_else(Oncotree_Date==min(Oncotree_Date),Oncotree_Date,min(Oncotree_Date)))%>%
  ungroup()%>%
  #Find earliest metastatic date for each type of metastasis
  group_by(MRN,Metastatic)%>%
  mutate(Metastatic_Date=if_else(Metastatic_Date==min(Metastatic_Date),Metastatic_Date,min(Metastatic_Date)))%>%
  #Reset grouping
  ungroup()%>%
  group_by(MRN)%>%
  #Fill columns by MRN
  fill(names(.),.direction="updown")%>%
  #Find earliest death date
  mutate(death_date=min(death_date))%>%
  #Find last visit date
  mutate(last_visit=if_else(visit_start_datetime>visit_end_datetime,visit_start_datetime,visit_end_datetime))%>%
  select(-visit_start_datetime,-visit_end_datetime)%>%
  unique()%>%
  #Determine survival status, need to consider NAs
  mutate(status=case_when(
    #check conditions
    !is.na(death_date)&!is.na(last_visit)&death_date<last_visit~"Deceased",
    !is.na(death_date)~"Deceased",
    is.na(death_date)&!is.na(last_visit)~"Alive",
    TRUE~"Unknown"
  ))%>%
  unique()%>%
  #Check slice
  #slice(which.min(Earliest_Oncotree))%>%
  ungroup()%>%
  unique()

write.table(outClinic,"~/Clinic-Data/20210829_metastatic.csv",sep=",",quote=FALSE,row.names=FALSE) #correct version
