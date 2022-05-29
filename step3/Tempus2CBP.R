library(tidyverse)
library(lubridate)
library(reshape2)
library(data.table)

#Import Tempus data, minimal QC to fill empty gene categories
mutationsTempus=read.csv("~/JSON/HashedTempus-Variants_08-29-2021.txt",sep="\t",colClasses="character")
patTempus=read.csv("~/JSON/HashedTempus-Patients_08-29-2021.txt",sep="\t",colClasses="character")
master=read.csv("~/Metadata.20210829.txt",sep="\t",colClasses="character")

mutationsTempus[mutationsTempus==""]=NA
patTempus[patTempus==""]=NA

mutationsTempus=mutationsTempus%>%
  drop_na(MRN_Hashed)%>%
  mutate(gene5=case_when(
    is.na(gene5)&is.na(gene5display)~gene,
    is.na(gene5)~gene5display,
    TRUE~gene5
  ))%>%
  mutate_all(~replace(.,str_detect(.,"EGFR"),"EGFR"))%>%
  mutate(hgncId=gsub("\\.0","",.$hgncId))%>%
  mutate(entrezId=gsub("\\.0","",.$entrezId))%>%
  unique()

patTempus=patTempus%>%
  drop_na(MRN_Hashed)%>%
  unique()

#####Write patient samples table
clinic=master%>%
  mutate(birth_datetime=as.Date(birth_datetime))%>%
  mutate(Oncotree_Date=as.Date(Oncotree_Date))%>%
  mutate(death_date=as.Date(death_date))%>%
  #mutate(Age=floor(time_length(difftime(Earliest_Diagnosis,birth_datetime),"years")))%>%
  unique()

merged=patTempus%>%
  select(MRN_Hashed,sex,Age_at_Report,tumorMutationalBurden,msiStatus)%>%
  left_join(clinic,by="MRN_Hashed")%>%
  #SEX for CBP
  mutate(sex=toupper(sex))%>%
  mutate(SEX=case_when(
    sex==gender~sex,
    is.na(sex)~gender,
    is.na(gender)~sex,
    sex!=gender~"UNKNOWN"
    ))%>%
  mutate(SEX=if_else(SEX=="MALE"|SEX=="FEMALE",SEX,"UNKNOWN"))%>%
  #AGE for CBP
  mutate(age=as.integer(floor(time_length(difftime(Anchor_Date,birth_datetime),"years"))))%>%
  mutate(Age_at_Report=as.integer(Age_at_Report))%>%
  mutate(AGE=case_when(
    age==Age_at_Report~age,
    is.na(age)~Age_at_Report,
    is.na(Age_at_Report)~age,
    age!=Age_at_Report~age
  ))%>%
  mutate(AGE=if_else(as.numeric(AGE)>90,90,as.numeric(AGE)))%>%
  #SURVIVAL for CBP
  mutate(Anchor_Date=as.Date(Anchor_Date),last_visit=as.Date(last_visit))%>%
  mutate(OS_MONTHS=case_when(
    status=="Alive"~floor(time_length(difftime(last_visit,Anchor_Date),"months")),
    status=="Deceased"~floor(time_length(difftime(death_date,Anchor_Date),"months")),
    TRUE~as.numeric(NA)
  ))%>%
  #mutate(OS_MONTHS=if_else(status=="Deceased",floor(time_length(difftime(death_date,signout_date),"months")),as.numeric(NA)))%>%  
  #STATUS for CBP
  mutate(OS_STATUS=case_when(
    status=="Alive"~"LIVING",
    status=="Deceased"~"DECEASED",
    is.na(status)~"LIVING"
  ))%>%
  #METASTATIC FOR CBP
  mutate(METASTATIC=case_when(
    Metastatic_Date<=Anchor_Date&!is.na(Metastatic)~"Met",
    Metastatic_Date>=Anchor_Date&!is.na(Metastatic)~"Non-Met",
    is.na(Metastatic)~"Non-Met",
    TRUE~"Unknown"
  ))%>%
  #ONCOTREE FOR CBP
  mutate(ONCOTREE=if_else(Oncotree_Date<=Anchor_Date&!is.na(Oncotree.level.1),Oncotree.level.1,NA_character_))%>%
  group_by(MRN_Hashed)%>%
  #select earliest Oncotree
  mutate(ONCOTREE=if_else(Oncotree_Date==min(Oncotree_Date),ONCOTREE,NA_character_))%>%
  rename(PATIENT_ID=MRN_Hashed,RACE=race,ETHNICITY=ethnicity)%>%
  unique()

cbioPat=merged%>%
  select(PATIENT_ID,OS_STATUS,OS_MONTHS,SEX,RACE,ETHNICITY,AGE,ONCOTREE,Oncotree_Date,METASTATIC,Metastatic_Date)%>%
  group_by(PATIENT_ID)%>%
  fill(names(.),.direction="downup")%>%
  #Convert metastatic, replacing unknowns with Met if one Met is before Anchor Date
  mutate(METASTATIC=case_when(
    any(METASTATIC=="Met")~"Met",
    TRUE~METASTATIC))%>%
  #Convert to multiple organ, replace Oncotree with Multiple_Organ if more than one organ before Anchor Date
  mutate(ONCOTREE=case_when(
    n_distinct(ONCOTREE)>1~"Multiple_Organ",
    TRUE~ONCOTREE
  ))%>%
  #Fix sex discrepancy for some patients
  mutate(SEX=case_when(
    any(SEX=="UNKNOWN")~"UNKNOWN",
    TRUE~SEX
  ))%>%
  filter(AGE>=18)%>%
  #Select youngest age
  replace_na(list(AGE=0))%>%
  slice(which.min(AGE))%>%
  mutate(AGE=na_if(AGE,0))%>%
  unique()

write.table(cbioPat,"~/JSON/cBioPortal/study/data_clinical_patient.txt",sep="\t",row.names = FALSE,quote=FALSE)

#####Write clinical samples table-Patient_ID and Accession_ID

tempDiagnosis=patTempus%>%
  select(MRN_Hashed,accessionId_Hashed,diagnosis,code,signout_date,sampleSite,tumorMutationalBurden,msiStatus,description)%>%
  unique()

anchor=master%>%
  select(MRN_Hashed,accessionId_Hashed,Anchor_Date)%>%
  unique()

cbioSample=mutationsTempus%>%
  select(MRN_Hashed,accessionId_Hashed)%>%
  drop_na(MRN_Hashed)%>%
  drop_na(accessionId_Hashed)%>%
  merge(tempDiagnosis,by=c("MRN_Hashed","accessionId_Hashed"),all=TRUE)%>%
  left_join(anchor,by=c("MRN_Hashed","accessionId_Hashed"))%>%
  #create report year
  mutate(REPORT_YEAR=year(as.Date(signout_date)))%>%
  #create days to report
  mutate(DAYS_TO_REPORT=as.numeric(as.Date(signout_date)-as.Date(Anchor_Date)))%>%
  rename(SAMPLE_ID=accessionId_Hashed,PATIENT_ID=MRN_Hashed,CANCER_TYPE=diagnosis,ASSAY=code,TMB_SCORE=tumorMutationalBurden,
         MICRO_INSTABILITY=msiStatus,SPECIMEN_SITE=sampleSite)%>%
  select(-signout_date,-Anchor_Date)%>%
  #Fix microsatellite instability
  mutate(MICRO_INSTABILITY=case_when(
    MICRO_INSTABILITY=="high"~"MSI-H",
    MICRO_INSTABILITY=="stable"~"MSS",
    MICRO_INSTABILITY=="undetermined"~"unknown",
    TRUE~"unknown"
  ))%>%
  #Add laboratory
  mutate(LABORATORY="Tempus")%>%
  #Add RNA
  mutate(RNA=if_else(grepl("RNA",description),"RNA",NA_character_))%>%
  select(-description)%>%
  unique()%>%
  group_by(SAMPLE_ID)%>%
  slice(which.min(DAYS_TO_REPORT),na.rm=TRUE)%>%
  unique()

write.table(cbioSample,"~/JSON/cBioPortal/study/data_clinical_samples.txt",sep="\t",row.names=FALSE,quote=FALSE)

#####Extended Mutations Table
mutationsExt=mutationsTempus%>%
  filter(is.na(mutation_type)|mutation_type!="Fusion Variants",is.na(variantType)|variantType!="fusion",is.na(variantDescription)|!grepl("copy number",variantDescription,ignore.case=TRUE),is.na(mutation_type)|!grepl("copy number",mutation_type,ignore.case=TRUE))%>%
  unique()

mutationsExt=mutationsExt%>%
  mutate(Variant_Classification=case_when(
    grepl("\\?",HGVS.p)~"Unknown",
    grepl("TERT",gene,ignore.case=TRUE)&grepl("variant",variantDescription,ignore.case=TRUE)&grepl("promoter",variantDescription,ignore.case=TRUE)~"5'Flank",
    #Del_ins first to avoid miscategorizing
    grepl("frameshift",variantDescription,ignore.case=TRUE)&grepl("delins",HGVS.c,ignore.case=TRUE)~"Del_Ins",#need to evaluate del_Ins cases
    grepl("frameshift",variantDescription,ignore.case=TRUE)&grepl("dup",HGVS.c,ignore.case=TRUE)~"Frame_Shift_Ins",
    grepl("frameshift",variantDescription,ignore.case=TRUE)&grepl("ins",HGVS.c,ignore.case=TRUE)~"Frame_Shift_Ins",
    grepl("frameshift",variantDescription,ignore.case=TRUE)&grepl("del",HGVS.c,ignore.case=TRUE)~"Frame_Shift_Del",
    grepl("inframe deletion",variantDescription,ignore.case=TRUE)~"In_Frame_Del",
    grepl("inframe insertion",variantDescription,ignore.case=TRUE)~"In_Frame_Ins",
    grepl("Missense variant",variantDescription,ignore.case=TRUE)~"Missense_Mutation",
    grepl("Splice region",variantDescription,ignore.case=TRUE)~"Splice_Region",
    grepl("Start",variantDescription,ignore.case=TRUE)~"Translation_Start_Site",
    grepl("Stop gain",variantDescription,ignore.case=TRUE)~"Nonsense_Mutation",
    grepl("Stop loss",variantDescription,ignore.case=TRUE)~"Nonstop_Mutation",
    grepl("Variant",variantDescription,ignore.case=TRUE)&grepl("\\-",HGVS.c)~"5'UTR",
    grepl("variant",variantDescription,ignore.case=TRUE)&grepl("\\>",HGVS.c)~"Silent",
    grepl("variant",variantDescription,ignore.case=TRUE)~"Variant",
    TRUE~"Unknown" #Need to check what it is
  ))

#Process Del_Ins
delIns=mutationsExt%>%
  filter(Variant_Classification=="Del_Ins")%>%
  unique()

# clean up HGVS string
delIns$HGVS.c.mod<-gsub("[()]","",delIns$HGVS.c)
delIns$HGVS.c.mod<-gsub("c\\.","",delIns$HGVS.c.mod)

# split HGVS string
delIns<-delIns%>%
  separate(HGVS.c.mod,into=c("loc","size"),sep="delins",remove = F)

#correct intron coordinate
delIns_nointron<-delIns%>%
  filter(!grepl("[+-]",loc))%>%
  separate(loc,into=c("beg","end"),sep="_")%>%
  mutate(len5=0)%>%
  mutate(len3=0)

delIns_intron5<-delIns%>%
  filter(grepl("-",loc))%>%
  separate(loc,into=c("beg1","end"),sep="_")%>%
  separate(beg1,into=c("beg","len5"),sep="-")%>%
  mutate(len3=0)

delIns_intron3<-delIns%>%
  filter(grepl("\\+",loc))%>%
  separate(loc,into=c("beg","end1"),sep="_")%>%
  separate(end1,into=c("end","len3"),sep="\\+")%>%
  mutate(len5=0)

delIns=rbind(delIns_nointron,delIns_intron5,delIns_intron3)

# fill in blank and convert to numeric
delIns<-delIns%>%
  mutate(end=ifelse(is.na(end),beg,end))
delIns$beg<-as.numeric(delIns$beg)
delIns$end<-as.numeric(delIns$end)
delIns$len3<-as.numeric(delIns$len3)
delIns$len5<-as.numeric(delIns$len5)

#calculate sizes
delIns<-delIns%>%
  mutate(size1=end-beg+1+len5+len3)
delIns<-delIns%>%
  mutate(size2=ifelse(grepl("[0-9]+",size),as.numeric(size),nchar(size)))

#compare sizes
delIns<-delIns%>%
  mutate(Variant_Classification=ifelse(size2>size1,"Frame_Shift_Ins","Frame_Shift_Del"))%>%
  select(-beg,-end,-HGVS.c.mod,-len3,-len5,-size,-size1,-size2)%>%
  unique()

#Merge tables together
mutationsExt=rbind(mutationsExt,delIns)

#Split & organize nucleotide alteration
mutationsExt=mutationsExt%>%
  rename(Hugo_Symbol=gene,HGVSp_Short=HGVS.p)%>%
  mutate(HGVSp_Short=ifelse(is.na(HGVSp_Short),"NA_splice",HGVSp_Short))%>%
  rename(Entrez_Gene_Id=entrezId,Tumor_Sample_Barcode=accessionId_Hashed,Variant_Type=variantType)%>%
  mutate(Variant_Type=if_else(is.na(Variant_Type),"UNK",Variant_Type))%>%
  mutate(Center="Tempus")%>%
  #Reference Genome: NCBI Build
  filter(!is.na(nucleotideAlteration))%>%
  separate(nucleotideAlteration,into=c("Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele1"),sep=":",remove=TRUE)%>%
  mutate(Start_Position=as.integer(Start_Position))%>%
  mutate(End_Position=Start_Position-nchar(Reference_Allele)+nchar(Tumor_Seq_Allele1))%>%
  mutate(Matched_Norm_Sample_Barcode=Tumor_Sample_Barcode)%>%
  mutate(allelicFraction=as.numeric(allelicFraction))%>%
  mutate(coverage=as.numeric(coverage))%>%
  mutate(coverage=replace_na(coverage,1000))%>%
  mutate(t_alt_count=round((coverage*allelicFraction)/100))%>%
  mutate(t_ref_count=round(coverage-t_alt_count))%>%
  mutate(Strand="",Tumor_Seq_Allele2="",dbSNP_RS="",dbSNP_Val_Status="",Match_Norm_Seq_Allele1="",Match_Norm_Seq_Allele2="",Tumor_Validation_Allele1="",Tumor_Validation_Allele2="",Match_Norm_Validation_Allele1="",
         Match_Norm_Validation_Allele2="",Verification_Status="Unknown",Validation_Status="Unknown",Mutation_Status="Somatic",Sequencing_Phase="",Sequence_Source="Exome",Validation_Method="",Score="",BAM_File="",n_alt_count="",n_ref_count="",
         Tumor_Sample_UUID="",Matched_Norm_Sample_UUID="",Comments="")%>%
  rename(cDNA_Change=HGVS.c)%>%
  rename(Transcript=transcript)%>%
  rename(Sequencer=code_name)%>%
  rename(NCBI_Build=referenceGenome)%>%
  unique()

mutationsExt$NCBI_Build[grepl("GRCh37",mutationsExt$NCBI_Build)]="GRCh37"

cBioMutExt=mutationsExt%>%
  select(Hugo_Symbol,Entrez_Gene_Id,Center,NCBI_Build,Chromosome,Start_Position,End_Position,Strand,Variant_Classification,Variant_Type,Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2,dbSNP_RS,dbSNP_Val_Status,Tumor_Sample_Barcode,Matched_Norm_Sample_Barcode,
         Match_Norm_Seq_Allele1,Match_Norm_Seq_Allele2,Tumor_Validation_Allele1,Tumor_Validation_Allele2,Match_Norm_Validation_Allele1,Match_Norm_Validation_Allele2,Verification_Status,Validation_Status,Mutation_Status,Sequencing_Phase,Sequence_Source,Validation_Method,Score,BAM_File,Sequencer,
         Tumor_Sample_UUID,Matched_Norm_Sample_UUID,HGVSp_Short,t_alt_count,t_ref_count,n_alt_count,n_ref_count,cDNA_Change,Comments,Transcript)%>%
  unique()

write.table(cBioMutExt,"~/JSON/cBioPortal/mutations/data_mutations_extended.txt",sep="\t",row.names=FALSE,quote=FALSE)
  
#####Copy Number Alterations Table

copyNum=mutationsTempus%>%
  #filter(grepl("copy number",variantDescription,ignore.case=TRUE)|grepl("copy number",mutation_type,ignore.case=TRUE))%>%
  select(accessionId_Hashed,variantDescription,mutation_type,gene,display,hgncId,entrezId,variantDescription,variantType,copyNumber)%>%
  mutate(status=case_when(
    (grepl("copy number",variantDescription,ignore.case=TRUE)|grepl("copy number",mutation_type,ignore.case=TRUE))&variantType=="deletion"~-2,
    (grepl("copy number",variantDescription,ignore.case=TRUE)|grepl("copy number",mutation_type,ignore.case=TRUE))&grepl("copy number loss",variantDescription,ignore.case=TRUE)&variantType=="CNV"~-1,
    (grepl("copy number",variantDescription,ignore.case=TRUE)|grepl("copy number",mutation_type,ignore.case=TRUE))&grepl("copy number gain",variantDescription,ignore.case=TRUE)&variantType=="CNV"~1,
    (grepl("copy number",variantDescription,ignore.case=TRUE)|grepl("copy number",mutation_type,ignore.case=TRUE))&grepl("copy number gain",variantDescription,ignore.case=TRUE)&variantType=="amplification"~2,
    ))%>%
  group_by(accessionId_Hashed,gene,entrezId)%>%
  fill(status,.direction="downup")%>%
  unique()

cBioCopyNum=copyNum%>%
  select(accessionId_Hashed,gene,entrezId,status)%>%
  unique()%>%
  filter(!is.na(status))%>%
  #add_count(gene,entrezId,accessionId_Hashed)
  reshape2::dcast(gene+entrezId~accessionId_Hashed,value.var="status")%>%
  mutate_all(~replace(., is.na(.), 0))%>%
  rename(Hugo_Symbol=gene,Entrez_Gene_Id=entrezId)

write.table(cBioCopyNum,"~/JSON/cBioPortal/mutations/data_CNA.txt",sep="\t",row.names=FALSE,quote=FALSE)

#####Fusions Table

#QC for fusions-some EGFR listed in gene/display columns and have fusion type, not listed in structuralVariant column
fusions=mutationsTempus%>%
  filter(mutation_type=="Fusion Variants"|variantType=="fusion")%>%
  select(accessionId_Hashed,description_name,code_name,gene5,gene5entrezId,gene3,gene3entrezId,structuralVariant)%>%
  unique()

fusions=fusions%>%
  #replace NA in gene5,gene3 with structural Variant
  mutate(gene5=ifelse(grepl("EGFR",structuralVariant),"EGFR",gene5),
         gene3=ifelse(grepl("EGFR",structuralVariant),"EGFR",gene3)
  )%>%
  #replace NA in gene3 with corresponding gene 5
  mutate(gene3=ifelse(is.na(gene3)&!is.na(gene5),gene5,gene3))%>%
  #assign entrezID for EGFR to gene3/gene5 entrezID
  mutate(gene5entrezId=ifelse(grepl("EGFR",gene5),"1956",gene5entrezId),
         gene3entrezId=ifelse(grepl("EGFR",gene3),"1956",gene3entrezId)
  )%>%
  mutate(DNA_support="yes")%>%
  mutate(RNA_support=ifelse(grepl("RNA",description_name),"yes","no"))%>%
  mutate(Fusion=case_when(
    gene5==gene3~"EGFR-intragenic",
    TRUE~paste(fusions$gene5,fusions$gene3,sep="-")
  ))%>%
  select(accessionId_Hashed,code_name,gene5,gene5entrezId,gene3,gene3entrezId,DNA_support,RNA_support,Fusion)%>%
  group_by(gene5,gene3)%>%
  fill(c(1:ncol(fusions)),.direction=c("downup"))%>%
  unique()

cbioFus=fusions%>%
  rename(gene5.hgnc=gene5,gene5.entrezId=gene5entrezId,gene3.hgnc=gene3,gene3.entrezId=gene3entrezId)%>%
  group_by(accessionId_Hashed)%>%
  reshape(direction='long',varying=list(c('gene5.hgnc','gene3.hgnc'),c('gene5.entrezId','gene3.entrezId')),timevar='var',times=c('gene5','gene3'),v.names=c('Hugo_Symbol','Entrez_Gene_Id'),idvar=c('accessionId_Hashed','code_name','DNA_support','RNA_support','Fusion'))%>%
  mutate(Frame="unknown")%>%
  mutate(Center="Tempus")%>%
  rename(Tumor_Sample_Barcode=accessionId_Hashed,Method=code_name)%>%
  select(Hugo_Symbol,Entrez_Gene_Id,Center,Tumor_Sample_Barcode,Fusion,DNA_support,RNA_support,Method,Frame)

write.table(cbioFus,"~/JSON/cBioPortal/mutations/data_fusions.txt",sep="\t",row.names=FALSE,quote=FALSE)








