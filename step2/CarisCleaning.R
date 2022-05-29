library(plyr)
library(dplyr)
library(stringr)

#switch mrn_new to mrn--mrn_new column is incorrect
pat=read.csv("~/Caris/Patients.txt",sep="\t",colClasses="character")
mut=read.csv("~/Caris/Mutations.txt",sep="\t",colClasses="character")

pat[pat==""]=NA
mut[mut==""]=NA

pat=pat[,colSums(is.na(pat))<nrow(pat)]
mut=mut[,colSums(is.na(mut))<nrow(mut)]

pat$MRN=str_pad(pat$MRN,8,pad="0")
mut$MRN=str_pad(mut$MRN,8,pad="0")

pat=pat%>%
	filter(nchar(MRN)==8)%>%
	unique()

mut=mut%>%
	filter(nchar(MRN)==8)%>%
	unique()

#mrn_new is correct, switch other MRNs
mrn=read.csv("~/Clinic-Data/clinicFiles/MRN_new.txt",sep="\t",colClasses="character")
mrn=mrn%>%
	dplyr::rename(MRN=mrn)

pat=pat%>%
	merge(mrn,by="MRN",all.x=TRUE)%>%
	mutate(MRN=if_else(!is.na(mrn_new),mrn_new,MRN))%>%
	select(-mrn_new)%>%
	unique()

mut=mut%>%
	merge(mrn,by="MRN",all.x=TRUE)%>%
	mutate(MRN=if_else(!is.na(mrn_new),mrn_new,MRN))%>%
	select(-mrn_new)%>%
	unique()

write.table(pat,"~/Caris/Patients.txt",sep="\t",row.names = FALSE)
write.table(mut,"~/Caris/Mutations.txt",sep="\t",row.names=FALSE)



