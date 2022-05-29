library(plyr)
library(dplyr)
library(tidyr)
library(lubridate)

FMpat=read.csv("~/ARMOR/HashedFM-Patients_08-29-2021.txt",sep="\t",colClasses = "character")
FMvar=read.csv("~/ARMOR/HashedFM-Variants_08-29-2021.txt",sep="\t",colClasses="character")
Tempus=read.csv("~/JSON/HashedTempus-Patients_08-29-2021.txt",sep="\t",colClasses = "character")
clinic=read.csv("~/Clinic-Data/20210829_metastatic.csv",sep=",",colClasses="character")
mrn=read.csv("~/Clinic-Data/clinicFiles/20210829_hashed/hashed_mrn_new.txt",sep="\t",colClasses="character")
Caris=read.csv("~/Caris/General/HashedCaris_Patients_08-29-2021.txt",sep="\t",colClasses="character")
guardant=read.csv("~/Guardant/HashedGuardant-Patients_08-29-2021.txt",sep="\t",colClasses="character")

FMpat[FMpat==""]=NA
FMvar[FMvar==""]=NA
Tempus[Tempus==""]=NA
clinic[clinic==""]=NA
mrn[mrn==""]=NA
Caris[Caris==""]=NA
guardant[guardant==""]=NA

FMpat=FMpat%>%
	select(MRN,MRN_Hashed,TRF,TRF_Hashed,Record.Date)%>%
	mutate(Anchor_Date=Record.Date)%>%
	unique()

FMvar=FMvar%>%
	select(MRN,MRN_Hashed,TRF,TRF_Hashed,Sample)%>%
	unique()

FM=FMpat%>%
	merge(FMvar,by=c("MRN","MRN_Hashed","TRF","TRF_Hashed"),all=TRUE)%>%
	group_by(MRN,MRN_Hashed,TRF,TRF_Hashed)%>%
	fill(names(.),.direction="downup")%>%
	unique()

Tempus=Tempus%>%
	select(MRN,MRN_Hashed,accessionId,accessionId_Hashed,signout_date,Sample)%>%
	mutate(Anchor_Date=signout_date)%>%
	group_by(MRN,MRN_Hashed,accessionId,accessionId_Hashed)%>%
	fill(names(.),.direction="downup")%>%
	unique()

Caris=Caris%>%
	select(MRN,MRN_Hashed,Caris_Report_ID,Caris_Report_ID_Hashed,Order_Date,Sample)%>%
	mutate(Anchor_Date=Order_Date)%>%
	group_by(MRN,MRN_Hashed,Caris_Report_ID,Caris_Report_ID_Hashed)%>%
	fill(names(.),.direction="downup")%>%
	unique()

guardant=guardant%>%
	select(MRN,MRN_Hashed,Guardant_Report_ID,Guardant_Report_ID_Hashed,Report_Date,Sample_ID)%>%
	mutate(Anchor_Date=Report_Date)%>%
	rename(Sample=Sample_ID)%>%
	group_by(MRN,MRN_Hashed,Guardant_Report_ID,Guardant_Report_ID_Hashed)%>%
	fill(names(.),.direction="downup")%>%
	unique()

data=plyr::rbind.fill(FM,Tempus,Caris,guardant)

mrn=mrn%>%
	select(mrn_new,Patient_ID_wrong)%>%
	dplyr::rename(MRN=mrn_new)%>%
	unique()

merged=clinic%>%
	left_join(data,by="MRN")%>%
	left_join(mrn,by="MRN")%>%
	unite(Sample_ID,c(accessionId,Caris_Report_ID,Guardant_Report_ID,TRF),na.rm=TRUE,remove=FALSE)%>%
	group_by(MRN)%>%
	mutate(Anchor_Date=min(as.Date(Anchor_Date)))%>%
	unique()

write.table(merged,"~/Metadata.20210829.txt",sep="\t",quote=FALSE,row.names=FALSE)
write.table(merged,"~/Clinic-Data/Metadata.20210829.txt",sep="\t",quote=FALSE,row.names=FALSE)
	
