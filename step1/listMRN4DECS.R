library(plyr)
library(dplyr)
library(stringr)

#####Foundation
patients=read.csv("~/ARMOR/data/test_processed/Summary",sep="\t",colClasses = "character")
variants=read.csv("~/ARMOR/data/test_processed/Variants",sep="\t",colClasses="character")

mrnFMp=patients%>%
  select(MRN)%>%
  filter(MRN!="")%>%
  mutate(MRN=str_pad(.$MRN,8,pad="0"))%>%
  mutate(MRN=str_trunc(.$MRN,8,"left"))%>%
  unique()

mrnFMv=variants%>%
  select(MRN)%>%
  mutate(MRN=str_pad(.$MRN,8,pad="0"))%>%
  mutate(MRN=str_trunc(.$MRN,8,"left"))%>%
  unique()

mrnFM=rbind(mrnFMp,mrnFMv)

mrnFM=unique(mrnFM)

#####Caris
patients=read.csv("~/Caris/Patients.txt",sep="\t",colClasses = "character")
mutations=read.csv("~/Caris/Mutations.txt",sep="\t",colClasses="character")

mrnCp=patients%>%
  select(MRN)%>%
  filter(MRN!="")%>%
  mutate(MRN=str_pad(.$MRN,8,pad="0"))%>%
  mutate(MRN=str_trunc(.$MRN,8,"left"))%>%
  unique()

mrnCv=mutations%>%
  select(MRN)%>%
  filter(MRN!="")%>%
  mutate(MRN=str_pad(.$MRN,8,pad="0"))%>%
  mutate(MRN=str_trunc(.$MRN,8,"left"))%>%
  unique()

mrnC=rbind(mrnCp,mrnCv)

MRN=rbind(mrnFM,mrnC)
MRN=unique(MRN)

#####Tempus
patients=read.csv("~/JSON/TempusPatients.csv",sep=",",colClasses="character")
variants=read.csv("~/JSON/TempusVariants.csv",sep=",",colClasses = "character")

mrnTp=patients%>%
  select(MRN)%>%
  filter(MRN!="")%>%
  mutate(MRN=str_pad(.$MRN,8,pad="0"))%>%
  mutate(MRN=str_trunc(.$MRN,8,"left"))%>%
  unique()

mrnTv=variants%>%
  select(MRN)%>%
  filter(MRN!="")%>%
  mutate(MRN=str_pad(.$MRN,8,pad="0"))%>%
  mutate(MRN=str_trunc(.$MRN,8,"left"))%>%
  unique()

mrnT=rbind(mrnTp,mrnTv)
mrnT=unique(mrnT)

MRN=rbind(MRN,mrnT)
MRN=unique(MRN)

###Guardant
patients=read.csv("~/Guardant/Patients.txt",sep="\t",colClasses="character")

mrnG=patients%>%
  select(MRN)%>%
  filter(MRN!="")%>%
  mutate(MRN=str_pad(.$MRN,8,pad="0"))%>%
  mutate(MRN=str_trunc(.$MRN,8,"left"))%>%
  unique()

MRN=rbind(MRN,mrnG)
MRN=unique(MRN)


write.table(MRN,"~/ARMOR_MRN_20210829.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)


