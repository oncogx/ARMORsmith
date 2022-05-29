library(dplyr)

aggPatients=function(path){
  setwd(path)
  file.list=list.files(path)
  for (file in file.list){
    if(!exists("patients")){
      load(file)
      patients=outPatients
    }
    
    if (exists("patients")){
      load(file)
      patients=rbind.fill(patients,outPatients)
      rm(outPatients)
    }
  }
  #needs tab separator to be robust against commas in fields
  write.table(write.table(patients,"~/JSON/patients.tsv",sep="\t"))
}