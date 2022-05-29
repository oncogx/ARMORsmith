library(dplyr)

aggMutations=function(path){
  setwd(path)
  file.list=list.files(path)
  for (file in file.list){
    if(!exists("mutations")){
      load(file)
      mutations=outMutations
    }
    
    if (exists("mutations")){
      load(file)
      mutations=rbind.fill(mutations,outMutations)
      rm(outMutations)
    }
  }
  #needs tab separator to be robust against commas in fields
  write.table(write.table(mutations,"~/JSON/mutations.tsv",sep="\t"))
}