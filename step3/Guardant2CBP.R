library(tidyverse)
library(lubridate)
library(reshape2)
library(data.table)

mut=read.csv("~/Guardant/HashedGuardant-Mutations_08-29-2021.txt",sep="\t",colClasses="character")
pat=read.csv("~/Guardant/HashedGuardant-Patients_08-29-2021.txt",sep="\t",colClasses="character")
master=read.csv("~/Metadata.20210829.txt",sep="\t",colClasses="character")

master=master%>%
	filter(!is.na(MRN))

mut[mut==""]=NA
pat[pat==""]=NA

mut=mut[,colSums(is.na(mut))<nrow(mut)]
pat=pat[,colSums(is.na(pat))<nrow(pat)]

###Patients table
masterPat=master%>%
	select(MRN_Hashed,birth_datetime,Anchor_Date,gender,race,ethnicity,death_date,last_visit,status,
				 Oncotree.level.1,Oncotree_Date,Metastatic,Metastatic_Date)%>%
	mutate(birth_datetime=as.Date(birth_datetime),Oncotree_Date=as.Date(Oncotree_Date),
				 death_date=as.Date(death_date),last_visit=as.Date(last_visit),Anchor_Date=as.Date(Anchor_Date),
				 Metastatic_Date=as.Date(Metastatic_Date))%>%
	unique()

cbpPat=masterPat%>%
	left_join(pat,by="MRN_Hashed")%>%
	filter(!is.na(Guardant_Report_ID))%>%
	select(-c("MRN","Guardant_Report_ID","Guardant_Report_ID_Hashed","Test_Type","Spec_Type","Sample_ID",
						"Disease","Alteration"))%>%
	#Standardize DOB
	mutate(DOB=as.Date(DOB))%>%
	mutate(DOB=if_else(DOB==birth_datetime,DOB,birth_datetime))%>%
	select(-birth_datetime)%>%
	#Standardize gender
	mutate(Gender=toupper(Gender))%>%
	mutate(Gender=if_else(gender==Gender,Gender,gender))%>%
	select(-gender)%>%
	#OS_Status
	mutate(OS_STATUS=case_when(
		status=="Alive"~"LIVING",
		status=="Deceased"~"DECEASED"
	))%>%
	select(-status)%>%
	#OS_SURVIVAL
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

write.table(cbpPat,"~/Guardant/cBioPortal/study/data_clinical_patient.txt",sep="\t",quote=FALSE,row.names=FALSE)

###Samples
tmbMsi=mut%>%
	filter(type=="tmb"|type=="msi")%>%
	select(MRN_Hashed,Guardant_Report_ID_Hashed,name,type,value,description)%>%
	#Select detected variants
	filter(value!="NOT DETECTED")%>%
	filter(value!="Not Evaluable")%>%
	#Remove extras
	mutate(value=gsub(" mut/MB","",value))%>%
	mutate(name=gsub("MSI-","",name))%>%
	#Micro_instability
	mutate(MICRO_INSTABILITY=if_else(type=="msi",name,NA_character_))%>%
	#TMB_Score
	mutate(TMB_SCORE=if_else(type=="tmb",value,NA_character_))%>%
	select(-name,-type,-value,-description)%>%
	#Fill values
	group_by(MRN_Hashed,Guardant_Report_ID_Hashed)%>%
	fill(names(.),.direction="downup")%>%
	#Fix Micro_instability
	mutate(MICRO_INSTABILITY=case_when(
		MICRO_INSTABILITY=="High"~"MSI-H",
		TRUE~"unknown"
	))%>%
	unique()

assayName=pat%>%
	select(MRN_Hashed,Guardant_Report_ID_Hashed,Test_Type)%>%
	mutate(Test_Type=if_else(Test_Type=="Guardant360"|Test_Type=="Guardant 360","Guardant360",Test_Type))%>%
	rename(ASSAY=Test_Type)%>%
	unique()

masterSample=master%>%
	select(MRN_Hashed,Guardant_Report_ID_Hashed,Anchor_Date)%>%
	filter(!is.na(Guardant_Report_ID_Hashed))%>%
	unique()

cbpSample=masterSample%>%
	left_join(pat,by=c("MRN_Hashed","Guardant_Report_ID_Hashed"))%>%
	left_join(tmbMsi,by=c("MRN_Hashed","Guardant_Report_ID_Hashed"))%>%
	left_join(assayName,by=c("MRN_Hashed","Guardant_Report_ID_Hashed"))%>%
	select(-MRN,-Guardant_Report_ID,-Sample_ID,-Alteration,-DOB,-Gender,-Test_Type)%>%
	#REPORT_YEAR, DAYS_TO_REPORT
	mutate(REPORT_YEAR=year(as.Date(Report_Date)))%>%
	mutate(DAYS_TO_REPORT=as.numeric(as.Date(Report_Date)-as.Date(Anchor_Date)))%>%
	select(-Anchor_Date,-Report_Date)%>%
	#SPECIMEN_TYPE
	mutate(Spec_Type=ifelse(grepl("Tissue",Spec_Type),"Other","Peripheral Blood"))%>%
	#SPECIMEN_SITE
	mutate(SPECIMEN_SITE=ifelse(grepl("Blood",Spec_Type),"Blood","Unknown"))%>%
	#Laboratory
	mutate(LABORATORY="Guardant")%>%
	#PATIENT_ID,SAMPLE_ID,SPECIMEN_SITE,ASSAY,CANCER_TYPE
	rename(PATIENT_ID=MRN_Hashed,SAMPLE_ID=Guardant_Report_ID_Hashed,SPECIMEN_TYPE=Spec_Type,
				 CANCER_TYPE=Disease)%>%
	#Reorder
	select(PATIENT_ID,SAMPLE_ID,REPORT_YEAR,DAYS_TO_REPORT,CANCER_TYPE,SPECIMEN_SITE,SPECIMEN_TYPE,ASSAY,TMB_SCORE,
				 MICRO_INSTABILITY,LABORATORY)%>%
	unique()

write.table(cbpSample,"~/Guardant/cBioPortal/study/data_clinical_sample.txt",sep="\t",quote=FALSE,row.names=FALSE)

###Mutations
entrezId=read.csv("~/Homo_sapiens.gene_info",sep="\t",colClasses="character")

entrezId=entrezId%>%
	select(GeneID,Symbol)%>%
	rename(Entrez_ID=GeneID,Hugo_Symbol=Symbol)%>%
	unique()

assays=pat%>%
	select(Guardant_Report_ID_Hashed,Test_Type)%>%
	unique()

mutationsExt=mut%>%
	left_join(entrezId,by=c("gene"="Hugo_Symbol"))%>%
	left_join(assays,by=c("Guardant_Report_ID_Hashed"))%>%
	filter(type!="tmb",type!="fusion",type!="cnv",type!="msi",method!="amplification",method!="fusion")%>%
	select(-MRN,-Guardant_Report_ID,-cfdnaPercentage,-detected,-value,-MRN_Hashed,-exon,-description)%>%
	unique()
	
mutationsExt=mutationsExt[,colSums(is.na(mutationsExt))<nrow(mutationsExt)]

cbpMut=mutationsExt%>%
	separate(nucleotideMutation,into=c("Reference_Allele","Tumor_Seq_Allele1"),sep=">",remove = FALSE)%>%
	#Set positions-End position only depends on length of refAllele
	mutate(length=if_else(is.na(length),as.character(nchar(Reference_Allele)),length))%>%
	mutate(End_Position=as.numeric(position)+as.numeric(length))%>%
	mutate(Start_Position=as.numeric(position))%>%
	select(-length)%>%
	#Create Variant_Classification
	mutate(Variant_Classification=case_when(
		grepl("fs",aminoAcidMutation)&nchar(Reference_Allele)>nchar(Tumor_Seq_Allele1)~"Frame_Shift_Del",
		grepl("fs",aminoAcidMutation)&nchar(Reference_Allele)<nchar(Tumor_Seq_Allele1)~"Frame_Shift_Ins",
		grepl("splice",name,ignore.case=TRUE)|grepl("splice",reportingCategory,ignore.case=TRUE)|grepl("splice",spliceEffect,ignore.case=TRUE)~"Splice_Site",
		grepl("promotor",name,ignore.case=TRUE)|grepl("promoter",reportingCategory,ignore.case=TRUE)~"Translation_Start_Site",
		method=="substitution"&grepl("\\*",name)~"Nonsense_Mutation",
		method=="substitution"&nchar(Reference_Allele)==nchar(Tumor_Seq_Allele1)~"Missense_Mutation",
		method=="insertion"&grepl("fs",name)&!grepl("\\*",name)~"Nonsense_Mutation",
		method=="deletion"&grepl("fs",name)&!grepl("\\*",name)~"Nonsense_Mutation",
		#Nonsense takes priority over indels
		method=="insertion"&(nchar(End_Position)-1)%%3==0~"In_Frame_Ins",
		method=="deletion"&(nchar(End_Position)-1)%%3==0~"In_Frame_Del",
		method=="insertion"~"Frame_Shift_Ins",
		method=="deletion"~"Frame_Shift_Del",
		method=="indel"&nchar(Reference_Allele)==nchar(Tumor_Seq_Allele1)~"Missense_Mutation",
		TRUE~"Check"
	))%>%
	#Variant_Type
	mutate(Variant_Type=case_when(
		nchar(Reference_Allele)==nchar(Tumor_Seq_Allele1)&nchar(Tumor_Seq_Allele1)==1~"SNP",
		nchar(Reference_Allele)==nchar(Tumor_Seq_Allele1)&nchar(Tumor_Seq_Allele1)==2~"DNP",
		nchar(Reference_Allele)==nchar(Tumor_Seq_Allele1)&nchar(Tumor_Seq_Allele1)==3~"TNP",
		nchar(Reference_Allele)==nchar(Tumor_Seq_Allele1)&nchar(Tumor_Seq_Allele1)>3~"ONP",
		nchar(Reference_Allele)>nchar(Tumor_Seq_Allele1)~"DEL",
		nchar(Reference_Allele)<nchar(Tumor_Seq_Allele1)~"INS",
		TRUE~"Check"
	))%>%
	#cDNA_Change
	mutate(cDNA_Change=gsub("c.","",cdna))%>%
	#HGVSp_Short
	mutate(HGVSp_Short=case_when(
		Variant_Classification=="Splice_Site"&is.na(aminoAcidMutation)~"NA_splice",
		Variant_Classification=="Translation_Start_Site"&is.na(aminoAcidMutation)~"NA_promoter",
		name==aminoAcidMutation~name,
		name!=aminoAcidMutation~aminoAcidMutation,
		TRUE~"Check"
	))%>%
	#Center,NCBI_Build,Strand,Verification_Status,Validation_Status,Mutation_Status,Sequence_Source,t_ref_count,t_alt_count
	mutate(Center="Guardant",NCBI_Build="GRCh37",Strand="",Verification_Status="unknown",
				 Validation_Status="unknown",Mutation_Status="Somatic",Sequence_Source="Exome",t_ref_count="",
				 t_alt_count="")%>%
	#Hugo_Symbol,Entrez_Gene_Id,Chromosome,Tumor_Sample_Barcode,Transcript,Sequencer
	rename(Hugo_Symbol=gene,Entrez_Gene_Id=Entrez_ID,Chromosome=chromosome,
				 Tumor_Sample_Barcode=Guardant_Report_ID_Hashed,Transcript=transcriptId,Sequencer=Test_Type)%>%
	select(Hugo_Symbol,Entrez_Gene_Id,Center,NCBI_Build,Chromosome,Start_Position,End_Position,Strand,
				 Variant_Classification,Variant_Type,Reference_Allele,Tumor_Seq_Allele1,Tumor_Sample_Barcode,
				 Verification_Status,Validation_Status,Mutation_Status,Sequence_Source,Sequencer,t_ref_count,
				 t_alt_count,cDNA_Change,HGVSp_Short,Transcript)%>%
	unique()

cbpMut[is.na(cbpMut)]=""

write.table(cbpMut,"~/Guardant/cBioPortal/mutations/data_mutations_extended.txt",sep="\t",row.names=FALSE,quote=FALSE)

###Fusions
entrezId=read.csv("~/Homo_sapiens.gene_info",sep="\t",colClasses="character")

entrezId=entrezId%>%
	select(GeneID,Symbol)%>%
	rename(Entrez_ID=GeneID,Hugo_Symbol=Symbol)%>%
	unique()

assays=pat%>%
	select(Guardant_Report_ID_Hashed,Test_Type)%>%
	unique()

fusions=mut%>%
	left_join(entrezId,by=c("gene"="Hugo_Symbol"))%>%
	left_join(assays,by=c("Guardant_Report_ID_Hashed"))%>%
	filter(type=="fusion")%>%
	select(-MRN,-MRN_Hashed,-Guardant_Report_ID)%>%
	unique()

fusions=unique(fusions[,colSums(is.na(fusions))<nrow(fusions)])

cbpFusions=fusions%>%
	select(gene,chromosome,downstreamGene,downstreamChromosome,Entrez_ID,Guardant_Report_ID_Hashed,Test_Type)%>%
	#Create Fusion column
	mutate(chromosome=paste("chr",chromosome))%>%
	mutate(downstreamChromosome=gsub(".0","",downstreamChromosome))%>%
	mutate(downstreamChromosome=if_else(downstreamChromosome!="",paste("chr",downstreamChromosome),
																			downstreamChromosome))%>%
	unite(Gene,c(gene,downstreamGene),sep="-",na.rm=TRUE,remove=FALSE)%>%
	mutate(Chromosome=if_else(downstreamChromosome!="",paste(chromosome,downstreamChromosome,sep="-"),
														paste(chromosome,"chrN",sep="-")))%>%
	select(-chromosome,-downstreamChromosome)%>%
	unite(Fusion,c(Gene,Chromosome),sep=" ",na.rm=TRUE,remove=TRUE)%>%
	#add DNA/RNA support
	mutate(DNA_support="yes",RNA_support="no",Frame="unknown")%>%
	#Center
	mutate(Center="Guardant")%>%
	#Rename
	rename(Hugo_Symbol=gene,Entrez_Gene_Id=Entrez_ID,Method=Test_Type,Tumor_Sample_Barcode=Guardant_Report_ID_Hashed)%>%
	select(Hugo_Symbol,Entrez_Gene_Id,Center,Tumor_Sample_Barcode,Fusion,DNA_support,RNA_support,Method,Frame)%>%
	unique()

write.table(cbpFusions,"~/Guardant/cBioPortal/mutations/data_fusions.txt",sep="\t",quote=FALSE,row.names=FALSE)

###CNA
entrezId=read.csv("~/Homo_sapiens.gene_info",sep="\t",colClasses="character")

entrezId=entrezId%>%
	select(GeneID,Symbol)%>%
	rename(Entrez_ID=GeneID,Hugo_Symbol=Symbol)%>%
	unique()

cna=mut%>%
	left_join(entrezId,by=c("gene"="Hugo_Symbol"))%>%
	filter(type=="cnv")%>%
	unique()
	
cna=unique(cna[,colSums(is.na(cna))<nrow(cna)])

cbpCNA=cna%>%
	select(Guardant_Report_ID_Hashed,gene,name,value,amplification,copyNumber,Entrez_ID)%>%
	filter(!is.na(amplification))%>%
	mutate(Rank=case_when(
		name=="Amplification"&amplification=="high"~2,
		name=="Amplification"&amplification=="medium"~1,
		name=="Amplification"&amplification=="low"~1
	))%>%
	select(gene,Entrez_ID,Guardant_Report_ID_Hashed,Rank)%>%
	reshape2::dcast(gene+Entrez_ID~Guardant_Report_ID_Hashed,value.var = "Rank")%>%
	mutate_all(~replace(., is.na(.), 0))%>%
	rename(Hugo_Symbol=gene,Entrez_Gene_Id=Entrez_ID)%>%
	unique()

write.table(cbpCNA,"~/Guardant/cBioPortal/mutations/data_CNA.txt",sep="\t",quote=FALSE,row.names=FALSE)



