library(tidyverse)
library(lubridate)
library(reshape2)
library(data.table)

mut=read.csv("~/ARMOR/HashedFM-Variants_08-29-2021.txt",sep="\t",colClasses="character")
pat=read.csv("~/ARMOR/HashedFM-Patients_08-29-2021.txt",sep="\t",colClasses="character")
master=read.csv("~/Metadata.20210829.txt",sep="\t",colClasses="character")

master=master%>%
	filter(!is.na(MRN))%>%
	unique()

mut[mut==""]=NA
pat[pat==""]=NA

mut=mut[,colSums(is.na(mut))<nrow(mut)]
pat=pat[,colSums(is.na(pat))<nrow(pat)]

mut=mut%>%
	filter(Category!="nonhuman")%>%
	unique()

pat=pat%>%
	filter(Quality!="Fail")%>%
	unique()

###Patients
masterPat=master%>%
	select(MRN_Hashed,birth_datetime,Anchor_Date,gender,race,ethnicity,death_date,last_visit,status,
				 Oncotree.level.1,Oncotree_Date,Metastatic,Metastatic_Date)%>%
	mutate(birth_datetime=as.Date(birth_datetime),Oncotree_Date=as.Date(Oncotree_Date),
				 death_date=as.Date(death_date),last_visit=as.Date(last_visit),Anchor_Date=as.Date(Anchor_Date),
				 Metastatic_Date=as.Date(Metastatic_Date))%>%
	unique()

cbpPat=masterPat%>%
	left_join(pat,by="MRN_Hashed")%>%
	filter(!is.na(TRF_Hashed)&TRF_Hashed%in%pat$TRF_Hashed)%>%
	select(-c("MRN","TRF","Sample","Short.Variants","Rearrangements","Copy.Number.Alterations",
						"Non.Human.Content","First.Name","Last.Name","Collection.Date","Disease","Disease.Ontology",
						"Purity.Assessment","Specimen.Site","Quality","Test.Type","Tumor.Mutation.Burden",
						"Tumor.Mutation.Burden.Score","Microsatellite.Instability","TRF_Hashed"))%>%
	#Standardize DOB
	mutate(DOB=as.Date(DOB))%>%
	mutate(DOB=if_else(DOB==birth_datetime,DOB,birth_datetime))%>%
	select(-birth_datetime)%>%
	#Standardize gender
	mutate(Gender=toupper(Gender),gender=toupper(gender))%>%
	mutate(Gender=if_else(gender==Gender,Gender,gender))%>%
	select(-gender)%>%
	#OS_STATUS
	mutate(OS_STATUS=case_when(
		status=="Alive"~"LIVING",
		status=="Deceased"~"DECEASED"
	))%>%
	select(-status)%>%
	#OS_MONTHS
	mutate(OS_MONTHS=case_when(
		OS_STATUS=="LIVING"~floor(time_length(difftime(last_visit,Anchor_Date),"months")),
		OS_STATUS=="DECEASED"~floor(time_length(difftime(death_date,Anchor_Date),"months")),
		TRUE~as.numeric(NA)
	))%>%
	select(-last_visit,-death_date)%>%
	#AGE
	mutate(AGE=floor(time_length(difftime(Anchor_Date,DOB),"years")))%>%
	mutate(AGE=if_else(AGE>=90,90,AGE))%>%
	select(-DOB,-Anchor_Date)%>%
	#METASTATIC,ONCOTREE
	group_by(MRN_Hashed)%>%
	fill(names(.),.direction="downup")%>%
	#Convert metastatic, replacing unknowns with Met if one Met is before Anchor Date
	mutate(METASTATIC=case_when(
		any(Metastatic=="Met")~"Met",
		TRUE~Metastatic))%>%
	#Convert to multiple organ, replace Oncotree with Multiple_Organ if more than one organ before Anchor Date
	mutate(ONCOTREE=case_when(
		n_distinct(Oncotree.level.1)>1~"Multiple_Organ",
		TRUE~Oncotree.level.1
	))%>%
	#RACE,ETHNICITY,PATIENT_ID,GENDER
	rename(RACE=race,ETHNICITY=ethnicity,PATIENT_ID=MRN_Hashed,SEX=Gender)%>%
	ungroup()%>%
	select(PATIENT_ID,OS_STATUS,OS_MONTHS,SEX,RACE,ETHNICITY,AGE,ONCOTREE,Oncotree_Date,
				 METASTATIC,Metastatic_Date)%>%
	unique()

write.table(cbpPat,"~/ARMOR/cBioPortal/study/data_clinical_patient.txt",sep="\t",quote=FALSE,row.names=FALSE)
	
###Samples
cbpSample=masterPat%>%
	left_join(pat,by="MRN_Hashed")%>%
	filter(!is.na(TRF_Hashed)&TRF_Hashed%in%pat$TRF_Hashed)%>%
	select(MRN_Hashed,Anchor_Date,Record.Date,Disease.Ontology,Specimen.Site,Test.Type,
				 Tumor.Mutation.Burden,Tumor.Mutation.Burden.Score,Microsatellite.Instability,TRF_Hashed)%>%
	#REPORT_YEAR, DAYS_TO_REPORT
	mutate(REPORT_YEAR=year(as.Date(Record.Date)))%>%
	mutate(DAYS_TO_REPORT=as.numeric(as.Date(Record.Date)-as.Date(Anchor_Date)))%>%
	select(-Anchor_Date,-Record.Date)%>%
	#SPECIMEN_TYPE
	mutate(SPECIMEN_TYPE=ifelse(grepl("Blood",Specimen.Site),"Peripheral Blood","Other"))%>%
	#LABORATORY
	mutate(LABORATORY="Foundation Medicine")%>%
	#TMB
	mutate(TMB=if_else(is.na(Tumor.Mutation.Burden),"Unknown",toupper(Tumor.Mutation.Burden)))%>%
	#MICRO_INSTABILITY
	mutate(MICRO_INSTABILITY=if_else(is.na(Microsatellite.Instability),"unknown",Microsatellite.Instability))%>%
	#TMB_SCORE,SPECIMEN_SITE,ASSAY,CANCER_TYPE,PATIENT_ID,SAMPLE_ID	
	rename(TMB_SCORE=Tumor.Mutation.Burden.Score,SPECIMEN_SITE=Specimen.Site,ASSAY=Test.Type,
				 CANCER_TYPE=Disease.Ontology,PATIENT_ID=MRN_Hashed,SAMPLE_ID=TRF_Hashed)%>%
	#Reorder
	select(PATIENT_ID,SAMPLE_ID,REPORT_YEAR,DAYS_TO_REPORT,CANCER_TYPE,SPECIMEN_SITE,SPECIMEN_TYPE,ASSAY,TMB_SCORE,
				 MICRO_INSTABILITY,LABORATORY)%>%
	unique()

write.table(cbpSample,"~/ARMOR/cBioPortal/study/data_clinical_sample.txt",sep="\t",row.names=FALSE,quote=FALSE)

###Mutations-Need to determine Variant Classification, Variant Type
entrezId=read.csv("~/Homo_sapiens.gene_info",sep="\t",colClasses="character")

entrezId=entrezId%>%
  select(GeneID,Symbol)%>%
  rename(Entrez_ID=GeneID,Hugo_Symbol=Symbol)%>%
  unique()

assays=pat%>%
  select(TRF_Hashed,Test.Type)%>%
  unique()

mutations=mut%>%
  left_join(entrezId,by=c("Gene"="Hugo_Symbol"))%>%
  left_join(assays,by=("TRF_Hashed"))%>%
  filter(TRF_Hashed%in%pat$TRF_Hashed)%>%
  filter(Category=="shortvariant"|Category=="rearrangement")%>%
  filter(Functional.Effect!="amplification",Functional.Effect!="fusion",Functional.Effect!="deletion")%>%
  select(-MRN,-TRF,-Sample,-Status,-Equivocal,-Subclonal,-MRN_Hashed)%>%
  unique()

mutations=unique(mutations[,colSums(is.na(mutations))<nrow(mutations)])

cbpMut=mutations%>%
  #Reference_Allele,Tumor_Seq_Allele1
  separate(CDS.Effect,into=c("Reference_Allele","Tumor_Seq_Allele1"),sep=">",remove = FALSE)%>%
  mutate(Tumor_Seq_Allele1=if_else(grepl(">",CDS.Effect),toupper(Tumor_Seq_Allele1),toupper(CDS.Effect)))%>%
  mutate(Reference_Allele=if_else(grepl(">",CDS.Effect),toupper(Reference_Allele),NA_character_))

cbpMut$Reference_Allele=gsub("[^A-Za-z]","",cbpMut$Reference_Allele)
cbpMut$Tumor_Seq_Allele1=gsub("[^A-Za-z]","",cbpMut$Tumor_Seq_Allele1)

#del end position=start position+length(deletion)-1
#ins end position=start position+1

cbpMut=cbpMut%>%
  #Chromosome,Start_Position,End_Position
  separate(Position,into=c("Chromosome","Start_Position"),sep=":",remove=TRUE)%>%
  mutate(End_Position=case_when(
    grepl("INS",Tumor_Seq_Allele1)~as.numeric(Start_Position)+1,
    grepl("DEL",Tumor_Seq_Allele1)~as.numeric(Start_Position)+nchar(Tumor_Seq_Allele1)-4,
    nchar(Reference_Allele)!=nchar(Tumor_Seq_Allele1)~as.numeric(Start_Position)+nchar(Reference_Allele)-1,
    nchar(Reference_Allele)==nchar(Tumor_Seq_Allele1)~as.numeric(Start_Position)+nchar(Reference_Allele)-1,
    TRUE~1000
  ))%>%
  #Variant_Classification
  mutate(Variant_Classification=case_when(
    Functional.Effect=="missense"~"Missense_Mutation",
    Functional.Effect=="nonsense"~"Nonsense_Mutation",
    Functional.Effect=="splice"~"Splice_Site",
    Functional.Effect=="promoter"~"5'Flank",
    Functional.Effect=="nonframeshift"&grepl("INS",Tumor_Seq_Allele1)~"In_Frame_Ins",
    Functional.Effect=="nonframeshift"&grepl("DEL",Tumor_Seq_Allele1)~"In_Frame_Del",
    Functional.Effect=="frameshift"&grepl("INS",Tumor_Seq_Allele1)~"Frame_Shift_Ins",
    Functional.Effect=="frameshift"&grepl("DEL",Tumor_Seq_Allele1)~"Frame_Shift_Del",
    Functional.Effect=="unknown"~"Unknown",
    TRUE~"Check"
    #truncation,duplication,rearrangement
  ))

cbpMut$Reference_Allele=gsub("INS|DEL","",cbpMut$Reference_Allele)
cbpMut$Tumor_Seq_Allele1=gsub("INS|DEL","",cbpMut$Tumor_Seq_Allele1)

cbpMut=cbpMut%>%
  #Variant_Classification
  mutate(Variant_Classification=case_when(
    Functional.Effect=="frameshift"&nchar(Reference_Allele)>nchar(Tumor_Seq_Allele1)~"Frame_Shift_Del",
    Functional.Effect=="frameshift"&nchar(Reference_Allele)<nchar(Tumor_Seq_Allele1)~"Frame_Shift_Ins",
    Functional.Effect=="nonframeshift"&nchar(Reference_Allele)>nchar(Tumor_Seq_Allele1)~"In_Frame_Del",
    Functional.Effect=="nonframeshift"&nchar(Reference_Allele)<nchar(Tumor_Seq_Allele1)~"In_Frame_Ins",
    Functional.Effect=="nonframeshift"&nchar(Reference_Allele)==nchar(Tumor_Seq_Allele1)&grepl("\\*",Protein.Effect)~"Nonsense_Mutation",
    Functional.Effect=="nonframeshift"&nchar(Reference_Allele)==nchar(Tumor_Seq_Allele1)&!grepl("\\*",Protein.Effect)~"Missense_Mutation",
    TRUE~Variant_Classification
  ))%>%
  #Variant_Type
  mutate(Variant_Type=case_when(
    Variant_Classification=="In_Frame_Ins"|Variant_Classification=="Frame_Shift_Ins"~"INS",
    Variant_Classification=="In_Frame_Del"|Variant_Classification=="Frame_Shift_Del"~"DEL",
    nchar(Reference_Allele)==nchar(Tumor_Seq_Allele1)&nchar(Reference_Allele)==1~"SNP",
    nchar(Reference_Allele)==nchar(Tumor_Seq_Allele1)&nchar(Reference_Allele)==2~"DNP",
    nchar(Reference_Allele)==nchar(Tumor_Seq_Allele1)&nchar(Reference_Allele)==3~"TNP",
    nchar(Reference_Allele)==nchar(Tumor_Seq_Allele1)&nchar(Reference_Allele)>3~"ONP",
    TRUE~"UNK"
  ))
  







###Fusions
entrezId=read.csv("~/Homo_sapiens.gene_info",sep="\t",colClasses="character")

entrezId=entrezId%>%
	select(GeneID,Symbol)%>%
	rename(Entrez_ID=GeneID,Hugo_Symbol=Symbol)%>%
	unique()

assays=pat%>%
	select(TRF_Hashed,Test.Type)%>%
	unique()

fusions=mut%>%
	left_join(entrezId,by=c("Gene"="Hugo_Symbol"))%>%
	left_join(assays,by=("TRF_Hashed"))%>%
	filter(TRF_Hashed%in%cbpSample$SAMPLE_ID)%>%
	filter(Functional.Effect=="fusion"|grepl("fusion",.$Description))%>%
	select(TRF_Hashed,Gene,Description,inFrame,Other.Gene,Entrez_ID,Test.Type)%>%
	unique()

fusions=unique(fusions[,colSums(is.na(fusions))<nrow(fusions)])

cbpFusions=fusions%>%
	mutate(Description=if_else(Description=="Rearrangement"|Description=="rearrangement"|
														 	Description=="Fusion"|Description=="fusion",
														 paste(Gene,Other.Gene,sep="-"),Description))%>%
	select(-Other.Gene)%>%
	#Frame
	mutate(Frame=case_when(
		inFrame=="No"~"frameshift",
		inFrame=="Yes"~"in-frame",
		TRUE~inFrame
	))%>%
	#DNA_support,RNA_support,Center
	mutate(DNA_support="yes",RNA_support="unknown",Center="Foundation Medicine")%>%
	#Method,Tumor_Sample_Barcode,Hugo_Symbol,Entrez_Gene_Id,Fusion
	rename(Method=Test.Type,Tumor_Sample_Barcode=TRF_Hashed,Hugo_Symbol=Gene,Entrez_Gene_Id=Entrez_ID,
				 Fusion=Description)%>%
	select(Hugo_Symbol,Entrez_Gene_Id,Center,Tumor_Sample_Barcode,Fusion,DNA_support,RNA_support,Method,Frame)%>%
	filter(Tumor_Sample_Barcode%in%pat$TRF_Hashed)%>%
  unique()

write.table(cbpFusions,"~/ARMOR/cBioPortal/mutations/data_fusions.txt",sep="\t",quote=FALSE,row.names=FALSE)

###CNA
entrezId=read.csv("~/Homo_sapiens.gene_info",sep="\t",colClasses="character")

entrezId=entrezId%>%
	select(GeneID,Symbol)%>%
	rename(Entrez_ID=GeneID,Hugo_Symbol=Symbol)%>%
	unique()

cna=mut%>%
	left_join(entrezId,by=c("Gene"="Hugo_Symbol"))%>%
	filter(Category=="copynumber")%>%
	unique()

cbpCNA=cna%>%
	select(TRF_Hashed,Gene,Entrez_ID,Functional.Effect,Copy.Number)%>%
	mutate(Copy.Number=as.numeric(Copy.Number))%>%
	mutate(Rank=case_when(
		Functional.Effect=="amplification"&Copy.Number<=170~1,
		Functional.Effect=="amplification"&Copy.Number>170~2,
		Functional.Effect=="loss"~-1,
		TRUE~1000
	))%>%
	select(-Copy.Number,-Functional.Effect)%>%
  filter(TRF_Hashed%in%pat$TRF_Hashed)%>%
	reshape2::dcast(Gene+Entrez_ID~TRF_Hashed,value.var = "Rank")%>%
	mutate_all(~replace(., is.na(.), 0))%>%
	rename(Hugo_Symbol=Gene,Entrez_Gene_Id=Entrez_ID)%>%
	unique()

write.table(cbpCNA,"~/ARMOR/cBioPortal/mutations/data_CNA.txt",sep="\t",quote=FALSE,row.names=FALSE)
	
	



