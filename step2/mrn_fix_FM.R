library(plyr)
library(dplyr)
library(stringr)

#mrn_new column is correct
pat=read.csv("~/ARMOR/data/test_processed/Summary",sep="\t",colClasses="character")
variant=read.csv("~/ARMOR/data/test_processed/Variants",sep="\t",colClasses="character")

pat$MRN=str_pad(pat$MRN,8,pad="0")
variant$MRN=str_pad(variant$MRN,8,pad="0")

mrn=read.csv("~/Clinic-Data/clinicFiles/MRN_new.txt",sep="\t",colClasses="character")
mrn=mrn%>%
	rename(MRN=mrn)

pat=pat%>%
	merge(mrn,by="MRN",all.x=TRUE)%>%
	mutate(MRN=if_else(!is.na(mrn_new),mrn_new,MRN))%>%
	select(-mrn_new)%>%
	unique()

variant=variant%>%
	merge(mrn,by="MRN",all.x=TRUE)%>%
	mutate(MRN=if_else(!is.na(mrn_new),mrn_new,MRN))%>%
	select(-mrn_new)%>%
	unique()

write.table(pat,"~/ARMOR/data/test_processed/Summary",sep="\t",row.names=FALSE)
write.table(variant,"~/ARMOR/data/test_processed/Variants",sep="\t",row.names=FALSE)