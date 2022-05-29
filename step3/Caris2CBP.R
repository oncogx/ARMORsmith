library(tidyverse)
library(lubridate)
library(reshape2)
library(data.table)

mut=read.csv("~/Caris/General/HashedCaris-Mutations_08-29-2021.txt",sep="\t",colClasses="character")
pat=read.csv("~/Caris/General/HashedCaris-Patients_08-29-2021.txt",sep="\t",colClasses="character")
master=read.csv("~/Metadata.20210829.txt",sep="\t",colClasses="character")

master=master%>%
	filter(!is.na(MRN))

mut[mut==""]=NA
pat[pat==""]=NA

mut=mut[,colSums(is.na(mut))<nrow(mut)]
pat=pat[,colSums(is.na(pat))<nrow(pat)]

mut=mut%>%
	select(-alleleFrequencyInformation,-alterationDetails,-readInformation)%>%
	unique()

#####Patients Table
masterPat=master%>%
	select(MRN_Hashed,birth_datetime,Anchor_Date,gender,race,ethnicity,death_date,last_visit,status,
				 Oncotree.level.1,Oncotree_Date,Metastatic,Metastatic_Date)%>%
	mutate(birth_datetime=as.Date(birth_datetime),Oncotree_Date=as.Date(Oncotree_Date),
				 death_date=as.Date(death_date),last_visit=as.Date(last_visit),Anchor_Date=as.Date(Anchor_Date),
				 Metastatic_Date=as.Date(Metastatic_Date))%>%
	unique()

cbpPat=masterPat%>%
	left_join(pat,by="MRN_Hashed")%>%
	filter(Test_Type %in% c("NGS","CNA-NGS","CNA-Seq","RNA-Seq","Seq"))%>%
	select(-c("MRN","Caris_Report_ID","Test_Type","Order_Date","Received_Date","ICD","Diagnosis",
	"Pathological_Diagnosis","Site","Lineage","Sublineage","Organization","Sample","Specimen_Type",
	"Specimen_Accession","Specimen_Site","Specimen_Date","Caris_Report_ID_Hashed"))%>%
	#standardize birth_date
	mutate(DOB=as.Date(DOB))%>%
	mutate(DOB=if_else(DOB==birth_datetime,DOB,birth_datetime))%>%
	select(-birth_datetime)%>%
	#standardize gender
	mutate(Gender=toupper(Gender))%>%
	mutate(Gender=if_else(gender==Gender,Gender,gender))%>%
	select(-gender)%>%
	#OS_STATUS
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

write.table(cbpPat,"~/Caris/cBioPortal/study/data_clinical_patient.txt",sep="\t",row.names=FALSE,quote=FALSE)

#####Samples Table
tmbMsiLOH=mut%>%
	filter(Category=="tumorMutationBurden"|Category=="microsatelliteInstability"|
				 	Category=="genomicLevelHeterozygosity")%>%
	select(MRN_Hashed,Caris_Report_ID_Hashed,Category,Result,msiCall,mutationBurdenScore,
				 mutationBurdenCall,LOHpercentage,Result)%>%
	#mutationBurdenCall-categorical
	mutate(mutationBurdenScore=gsub(" per Mb","",mutationBurdenScore))%>%
	rename(TMB_SCORE=mutationBurdenScore,MICRO_INSTABILITY=msiCall,TMB=mutationBurdenCall)%>%
	#TMB_VALUE=mutationBurdenCall-categorical ranking
	mutate(LOSS_OF_HETEROZYGOSITY=if_else(Category=="genomicLevelHeterozygosity",Result,NA_character_))%>%
	select(-Category,-Result)%>%
	group_by(MRN_Hashed,Caris_Report_ID_Hashed)%>%
	fill(names(.),.direction="downup")%>%
	mutate(MICRO_INSTABILITY=case_when(
		MICRO_INSTABILITY=="High"~"MSI-H",
		TRUE~"unknown"
	))%>%
	unique()

makeShortName<-mut%>% 
	select(Caris_Report_ID_Hashed,Test_Code,Test_Name,Platform_Technology,analysisConfigurationName,
				 analysisPipelineName)%>% 
	unique()%>% 
	group_by(Test_Code,Test_Name,Platform_Technology,analysisConfigurationName,analysisPipelineName)%>%
	mutate(shortName=ifelse(grepl("600",analysisConfigurationName),"DNA600","other"))%>% 
	mutate(shortName=ifelse(grepl("720",analysisConfigurationName),"DNA720",shortName))%>% 
	mutate(shortName=ifelse(grepl("Exome",analysisConfigurationName),"WXS",shortName))%>% 
	mutate(shortName=ifelse(grepl("WTS",analysisConfigurationName),"RNA-Seq",shortName))%>% 
	mutate(shortName=ifelse(grepl("Archer",analysisConfigurationName),"other_Fusion",shortName))%>% 
	mutate(shortName=ifelse(grepl("AMBRY",analysisPipelineName),"germline_panel",shortName))%>% 
	mutate(shortName=ifelse(grepl("TruSeq|BRCA",Test_Name),"DNA_other",shortName)) 

reportTypes<-mut%>% 
	select(Caris_Report_ID_Hashed,Test_Code,Test_Name,Platform_Technology,analysisConfigurationName,
				 analysisPipelineName)%>% 
	unique() 

Combo_agg<-left_join(reportTypes,makeShortName) %>% 
	select(Caris_Report_ID_Hashed,shortName) %>% 
	unique() %>% 
	group_by(Caris_Report_ID_Hashed) %>% 
	summarize(combo=paste(shortName,collapse=":"))

mut=mut%>%
	left_join(Combo_agg,by="Caris_Report_ID_Hashed")%>%
	unique()

masterSample=master%>%
	select(MRN_Hashed,Caris_Report_ID_Hashed,Anchor_Date)%>%
	filter(!is.na(Caris_Report_ID_Hashed))%>%
	unique()

cbpSample=masterSample%>%
	left_join(pat,by=c("MRN_Hashed","Caris_Report_ID_Hashed"))%>%
	left_join(Combo_agg,by="Caris_Report_ID_Hashed")%>%
	filter(Test_Type %in% c("NGS","CNA-NGS","CNA-Seq","RNA-Seq","Seq"))%>%
	select(-Test_Type)%>%
	left_join(tmbMsiLOH,by=c("MRN_Hashed","Caris_Report_ID_Hashed"))%>%
	select(-MRN,-Caris_Report_ID,-Received_Date,-DOB,-Gender,-ICD,-Organization,-Sample,
				 -Specimen_Accession,-Specimen_Date)%>%
	#REPORT_YEAR,DAYS_TO_REPORT
	mutate(REPORT_YEAR=year(as.Date(Order_Date)))%>%
	mutate(DAYS_TO_REPORT=as.numeric(as.Date(Order_Date)-as.Date(Anchor_Date)))%>%
	select(-Anchor_Date,-Order_Date)%>%
	#SPECIMEN_TYPE
	mutate(SPECIMEN_TYPE=ifelse(grepl("Tissue",Specimen_Type),"Other","Blood"))%>%
	select(-Specimen_Type)%>%
	#PATIENT_ID,SAMPLE_ID,SPECIMEN_SITE,ASSAY,CANCER_TYPE,LOH_PERCENTAGE
	rename(PATIENT_ID=MRN_Hashed,SAMPLE_ID=Caris_Report_ID_Hashed,SPECIMEN_SITE=Specimen_Site,
				 ASSAY=combo,CANCER_TYPE=Diagnosis,LOH_PERCENTAGE=LOHpercentage)%>%
	#LABORATORY
	mutate(LABORATORY="Caris")%>%
	#RNA
	mutate(RNA=if_else(grepl("RNA",ASSAY),"RNA",NA_character_))%>%
	select(PATIENT_ID,SAMPLE_ID,REPORT_YEAR,DAYS_TO_REPORT,CANCER_TYPE,SPECIMEN_SITE,SPECIMEN_TYPE,
				 ASSAY,TMB_SCORE,TMB,MICRO_INSTABILITY,LOSS_OF_HETEROZYGOSITY,LOH_PERCENTAGE,LABORATORY,RNA)%>%
	unique()

write.table(cbpSample,"~/Caris/cBioPortal/study/data_clinical_sample.txt",sep="\t",quote=FALSE,row.names=FALSE)

#####Mutations
entrezId=read.csv("~/Homo_sapiens.gene_info",sep="\t",colClasses="character")

entrezId=entrezId%>%
	select(GeneID,Symbol)%>%
	rename(Entrez_ID=GeneID,Hugo_Symbol=Symbol)%>%
	unique()

mutationsExt=mut%>%
	left_join(entrezId,by=c("gene"="Hugo_Symbol"))%>%
	filter(Category=="genomicAlteration")%>%
	unique()

mutationsExt=mutationsExt[,colSums(is.na(mutationsExt))<nrow(mutationsExt)]

mutationsExt=mutationsExt%>%
	select(-MRN,-Caris_Report_ID,-Category,-Result,-biomarkerName,-MRN_Hashed,-NGSPanelName,
				 -NGSPanelVersion,-Platform_Technology,-Test_Code,-Test_Name,-Test_Type,-analysisConfigurationName,
				 -analysisConfigurationVersion,-analysisPipelineName,-analysisPipelineVersion,)%>%
	#NCBI_Build
	mutate(genomeBuild=case_when(
		grepl("GRCh38",genomeBuild)~"GRCh38",
		grepl("GRCh37",genomeBuild)~"GRCh37",
		TRUE~"GRCh37"
	))%>%
	rename(NCBI_Build=genomeBuild)%>%
	filter(!is.na(NCBI_Build))%>%
	#Chromosome
	mutate(chromosome=gsub("chr","",chromosome))%>%
	rename(Chromosome=chromosome)%>%
	#cDNA_Change
	mutate(hgvsCodingChange=gsub("c.","",hgvsCodingChange))%>%
	rename(cDNA_Change=hgvsCodingChange)%>%
	#HGVSp_Short
	mutate(hgvsProteinChange=gsub("p.","",hgvsProteinChange))%>%
	rename(HGVSp_Short=hgvsProteinChange)%>%
	#Rename Hugo_Symbol, Mutation_Status, Entrez_Gene_Id, Tumor_Sample_Barcode, Reference_Allele, 
	#Transcript, Tumor_Seq_Allele1, combo
	rename(Hugo_Symbol=gene,Entrez_Gene_Id=Entrez_ID,Mutation_Status=genomicSource,
				 Tumor_Sample_Barcode=Caris_Report_ID_Hashed,Transcript=transcriptID,
				 Tumor_Seq_Allele1=observedNucleotide,Reference_Allele=referenceNucleotide,Sequencer=combo)%>%
	#Fill NAs in Tumor_Seq_Allele1,Reference_Allele with ""
	mutate(Tumor_Seq_Allele1=if_else(is.na(Tumor_Seq_Allele1),"",Tumor_Seq_Allele1))%>%
	mutate(Reference_Allele=if_else(is.na(Reference_Allele),"",Reference_Allele))%>%
	#Start_Position, End_Position
	mutate(transcriptStartPosition=as.integer(transcriptStartPosition),
				 transcriptStopPosition=as.integer(transcriptStopPosition))%>%
	rename(Start_Position=transcriptStartPosition,End_Position=transcriptStopPosition)%>%
	#Variant Classification
	mutate(Variant_Classification=case_when(
		molecularConsequence=="Nonsense"~"Nonsense_Mutation",
		molecularConsequence=="Missense"~"Missense_Mutation",
		molecularConsequence=="Silent"~"Silent",
		molecularConsequence=="Splicing"~"Splice_Site",
		molecularConsequence=="Noncoding"~"5'UTR",
		molecularConsequence=="Promoter"~"5'Flank",
		molecularConsequence=="CODON_DELETION"~"In_Frame_Del",
		molecularConsequence=="CODON_INSERTION"~"In_Frame_Ins",
		molecularConsequence=="CODON_CHANGE_PLUS_CODON_DELETION"~"In_Frame_Del",
		molecularConsequence=="Frameshift"&(nchar(Tumor_Seq_Allele1)>nchar(Reference_Allele))~"Frame_Shift_Ins",
		molecularConsequence=="Frameshift"&(nchar(Tumor_Seq_Allele1)<nchar(Reference_Allele))~"Frame_Shift_Del",
		molecularConsequence=="Frameshift"&is.na(Tumor_Seq_Allele1)&is.na(Reference_Allele)~"Unknown",
		TRUE~"Check"
	))%>%
	mutate(Variant_Classification=case_when(
		#For NAs in molecularConsequence and full stop in protein
		Variant_Classification=="Check"&grepl("fs|\\*",HGVSp_Short)~"Nonsense_Mutation",
		#For NAs in molecularConsequence and missense mutation in protein
			#Check if first and last character in HGVS are proteins and if they are different
		Variant_Classification=="Check"&str_detect(substring(HGVSp_Short,1,1),"^[:upper:]+$")&
			str_detect(str_sub(HGVSp_Short,-1,-1),"^[:upper:]+$")&
			substring(HGVSp_Short,1,1)!=str_sub(HGVSp_Short,-1,-1)~"Missense_Mutation",
		Variant_Classification=="Check"&grepl("del",HGVSp_Short)~"Splice_Site",
		TRUE~as.character(Variant_Classification)
	))%>%
	#Variant_Type
	mutate(Variant_Type=case_when(
		nchar(Tumor_Seq_Allele1)==1&nchar(Reference_Allele)==1~"SNP",
		nchar(Tumor_Seq_Allele1)==2&nchar(Reference_Allele)==2~"DNP",
		nchar(Tumor_Seq_Allele1)==3&nchar(Reference_Allele)==3~"TNP",
		Tumor_Seq_Allele1!=""&Reference_Allele!=""&nchar(Tumor_Seq_Allele1)==nchar(Reference_Allele)~"ONP",
		nchar(Tumor_Seq_Allele1)>nchar(Reference_Allele)~"INS",
		nchar(Tumor_Seq_Allele1)<nchar(Reference_Allele)~"DEL",
		TRUE~"UNK"
	))%>%
	#t_alt_count, t_ref_count
	mutate(t_alt_count=round(as.numeric(readDepth)*as.numeric(alleleFrequency)/100))%>%
	mutate(t_ref_count=round(as.numeric(readDepth)-t_alt_count))%>%
	#Center, Validation_Status, Verification_Status, Sequence_Source, Strand
	mutate(Center="Caris",Validation_Status="unknown",Verification_Status="unknown",
				 Sequence_Source="Exome",Strand="")%>%
	select(Hugo_Symbol,Entrez_Gene_Id,Center,NCBI_Build,Chromosome,Start_Position,End_Position,Strand,
				 Variant_Classification,Variant_Type,Reference_Allele,Tumor_Seq_Allele1,Tumor_Sample_Barcode,
				 Verification_Status,Validation_Status,Mutation_Status,Sequence_Source,Sequencer,t_ref_count,
				 t_alt_count,cDNA_Change,HGVSp_Short,Transcript)%>%
	unique()

write.table(mutationsExt,"~/Caris/cBioPortal/mutations/data_mutations_extended.txt",sep="\t",quote=FALSE,row.names=FALSE)

#####Fusions
fusions=mut%>%
	left_join(entrezId,by=c("gene"="Hugo_Symbol"))%>%
	filter(Category=="translocation")%>%
	unique()

cbpFusions=fusions%>%
	select(gene,Entrez_ID,Caris_Report_ID_Hashed,fusionISOForm,combo,Test_Type)%>%
	rename(Hugo_Symbol=gene,Entrez_Gene_Id=Entrez_ID,Tumor_Sample_Barcode=Caris_Report_ID_Hashed,
				 Fusion=fusionISOForm,Method=combo)%>%
	mutate(Center="Caris",Frame="unknown")%>%
	mutate(DNA_support=if_else(grepl("DNA",Method),"yes","no"))%>%
	mutate(RNA_support=if_else(grepl("RNA",Method)|grepl("RNA",Test_Type),"yes","no"))%>%
	select(Hugo_Symbol,Entrez_Gene_Id,Center,Tumor_Sample_Barcode,Fusion,DNA_support,RNA_support,
				 Method,Frame)%>%
	unique()

write.table(cbpFusions,"~/Caris/cBioPortal/mutations/data_fusions.txt",sep="\t",row.names=FALSE,quote=FALSE)

#####CNA
cna=mut%>%
	left_join(entrezId,by=c("gene"="Hugo_Symbol"))%>%
	filter(Category=="copyNumberAlteration")%>%
	unique()

cna=cna%>%
	select(Caris_Report_ID_Hashed,gene,Entrez_ID,copyNumberType)%>%
	mutate(copyNumberType=case_when(
		copyNumberType=="Amplified"~"2",
		copyNumberType=="Intermediate"~"1",
		TRUE~"100"
	))%>%
	reshape2::dcast(gene+Entrez_ID~Caris_Report_ID_Hashed,value.var = "copyNumberType")%>%
	mutate_all(~replace(., is.na(.), 0))%>%
	rename(Hugo_Symbol=gene,Entrez_Gene_Id=Entrez_ID)%>%
	unique()

write.table(cna,"~/Caris/cBioPortal/mutations/data_CNA.txt",sep="\t",row.names=FALSE,quote=FALSE)


