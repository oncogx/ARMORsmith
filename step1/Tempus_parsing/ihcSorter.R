library(jsonlite)
library(stringr)
library(dplyr)

ihcSort=function(filename){
  name=read_json(filename)
  title=str_remove(basename(filename),".json")
  options(max.print=5000)
  if (name$report$workflow$reportType!="DNA"){
    cat(sprintf("%s report is not DNA, moved\n",title))
    
    if (grepl("prospective",getwd())){
      from=paste("~/JSON/prospective/",title,".json",sep='')
    }else if (grepl("retrospective",getwd())){
      from=paste("~/JSON/retrospective/",title,".json",sep='')
    }
    #from=paste("~/JSON/retrospective/",title,".json",sep='')
    to=paste("~/JSON/nonDNA/",title,".json",sep='')
    file.rename(from=from,to=to)
    #print("Y")
  }
  #else{
    
  #}

}

