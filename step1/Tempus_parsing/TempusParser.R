library(jsonlite)
library(data.table)
library(plyr)
library(dplyr)
library(stringr)
library(tidyr)

tempusParse=function(filename){
  parsefile=read_json(filename)
  title=str_remove(basename(filename),".json")
  print(title)
  patCount=0
  mutCount=0
  schema=parsefile$metadata$schemaVersion
 
        #Extract report info
      if (exists("report",parsefile)&&length(parsefile$report)!=0){
        reportData=rbindlist(list(parsefile$report,parsefile$report$workflow),fill=TRUE)
        reportData$workflow=NULL
        reportData[is.na(reportData)]=""
        reportData=reportData%>%
          unique()%>%
          summarize_all(list(name=~(trimws(paste(.,collapse='')))))
        #assign(paste("reportData_",title,sep=''),reportData,envir=.GlobalEnv)
        if (exists("reportData")){
          patCount=patCount+1
        }
      }else{
        #print("No report info")
      }
      
      #Extract patient info  
      if(exists("patient",parsefile)&&length(parsefile$patient)!=0){
        patientData=as.data.frame(t(unlist(parsefile$patient)))
        #assign(paste("patientData_",title,sep=''),patientData,envir=.GlobalEnv)
        if (exists("patientData")){
          patCount=patCount+1
        }
      }else{
        #print("No patient info")
      }
      
      #Extract order info
      if (exists("report",parsefile)&&length(parsefile$report)!=0){
        orderData=rbindlist(list(parsefile$order,parsefile$order$test),fill=TRUE)
        orderData$test=NULL
        orderData[is.na(orderData)]=""
        orderData=orderData%>%
          unique()%>%
          #summarize_all(funs(trimws(paste(.,collapse=''))))
          summarize_all(list(name=~(trimws(paste(.,collapse='')))))
        if (exists("orderData")){
          patCount=patCount+1
        }
        #assign(paste("orderData_",title,sep=''),orderData,envir=.GlobalEnv)
      }else{
        #print("No order info")
      }
      
      #Extract specimens info
      if (exists("specimens",parsefile)&&length(parsefile$specimens)!=0){
        specimens=list()
        for(i in 1:length(parsefile$specimens)){
          specimens[[i]]=as.data.frame(t(unlist(parsefile$specimens[[i]])))
        }
        specimensData=do.call(rbind.fill,specimens)%>%
          unique()
        #assign(paste("specimensData_",title,sep=''),specimensData,envir=.GlobalEnv)
        if (exists("specimensData")){
          patCount=patCount+1
        }
      }else{
        #print("No specimens info")
      }
      
      #Extract   mutations info
      if (exists("results",parsefile)&&length(parsefile$results)!=0){
        #Tumor Mutational Burden
        if (length(parsefile$results$tumorMutationalBurden!=0)){
          tmbData=as.data.frame(t(unlist(c(parsefile$results[1],parsefile$results[2],parsefile$results[3]))))
          #assign(paste("tmbData_",title,sep=''),tmbdata,envir=.GlobalEnv)
          if (exists("tmbData")){
            patCount=patCount+1
          }
        }else{
          #print("No TMB")
        }
        #Somatic Potentially Actionable Mutations
        if (exists("somaticPotentiallyActionableMutations",parsefile$results)&&length(parsefile$results$somaticPotentiallyActionableMutations)!=0){
          #Remove Therapies
          for (i in 1:length(parsefile$results$somaticPotentiallyActionableMutations)){
            for (j in 1:length(parsefile$results$somaticPotentiallyActionableMutations[[i]]$variants)){
              parsefile$results$somaticPotentiallyActionableMutations[[i]]$variants[[j]]$therapies=NULL
            }
          }
          mutations=list()
          variants=list()
          #Extract Somatic Potentially Actionable Mutations data
          #Per Entry
          for(i in 1:length(parsefile$results$somaticPotentiallyActionableMutations)){
            #Per Variants in Entry
            for(j in 1:length(parsefile$results$somaticPotentiallyActionableMutations[[i]]$variants)){
              variants[[j]]=rbindlist(list(parsefile$results$somaticPotentiallyActionableMutations[[i]]$variants[j],parsefile$results$somaticPotentiallyActionableMutations[[i]]$variants[[j]]),fill=TRUE)
            }
            #Per Variants in Entry
            for (j in 1:length(variants)){
              variantsData=do.call(rbind,variants[j])
              mutationsList=parsefile$results$somaticPotentiallyActionableMutations[[i]]
              mutationsList$variants=NULL
              mutationTable=rbindlist(list(mutationsList,variantsData),fill=TRUE)
              mutations=append(mutations,list(rbindlist(list(mutationsList,variantsData),fill=TRUE)))
            }
          }
          #Build SPAMS table
          SPAMsData=do.call(rbind,mutations)
          SPAMsData$mutationEffect=NULL
          SPAMsData=SPAMsData%>%
            unique()%>%
            fill(c(1:4),.direction=c("down"))%>%
            fill(c(5:ncol(SPAMsData)),.direction=c("up"))%>%
            unique()%>%
            mutate(mutation_type="Somatic Potentially Actionable Mutation")
          if (exists("SPAMsData")){
            mutCount=mutCount+1
          }
        }else{
          #print("No SPAMS")
        }
        
        #Somatic Potentially Actionable Copy Number Variants
        if (exists("somaticPotentiallyActionableCopyNumberVariants",parsefile$results)&&length(parsefile$results$somaticPotentiallyActionableCopyNumberVariants)!=0){
          #Remove Therapies
          for (i in 1:length(parsefile$results$somaticPotentiallyActionableCopyNumberVariants)){
            parsefile$results$somaticPotentiallyActionableCopyNumberVariants[[i]]$therapies=NULL
          } 
          #Extract Somatic Potentially Actionable Copy Number Variants Data
          variants=list()
          for(i in 1:length(parsefile$results$somaticPotentiallyActionableCopyNumberVariants)){
            variants[[i]]=as.data.frame(t(unlist(parsefile$results$somaticPotentiallyActionableCopyNumberVariants[[i]])))
          }
          SPACNVsData=do.call(rbind,variants)%>%
            mutate(mutation_type="Somatic Potentially Actionable Copy Number Variants")
          rm(variants)
          if (exists("SPACNVsData")){
            mutCount=mutCount+1
          }
        }else{
          #print("No SPACNVs")
        }
        
        #Somatic Biologically Relevant Variants
        if (exists("somaticBiologicallyRelevantVariants",parsefile$results)&&length(parsefile$results$somaticBiologicallyRelevantVariants)!=0){
          #Extract Data
          variants=list()
          for(i in 1:length(parsefile$results$somaticBiologicallyRelevantVariants)){
            variants[[i]]=as.data.frame(t(unlist(parsefile$results$somaticBiologicallyRelevantVariants[[i]])))
          }
          SBRVsData=do.call(rbind,variants)%>%
            mutate(mutation_type="Somatic Biologically Relevant Variants")
          rm(variants)
          if (exists("SBRVsData")){
            mutCount=mutCount+1
          }
        }else{
          #print("No SBRVs")
        }
        
        #Somatic Variants of Unknown Significance
        if (exists("somaticVariantsOfUnknownSignificance",parsefile$results)&&length(parsefile$results$somaticVariantsOfUnknownSignificance)!=0){
          #Extract Data
          variants=list()
          for(i in 1:length(parsefile$results$somaticVariantsOfUnknownSignificance)){
            variants[[i]]=as.data.frame(t(unlist(parsefile$results$somaticVariantsOfUnknownSignificance[[i]])))
          }
          SVUSsData=do.call(rbind,variants)%>%
            mutate(mutation_type="Somatic Variants of Unknown Significance")
          rm(variants)
          if (exists("SVUSsData")){
            mutCount=mutCount+1
          }
        }else{
          #print("No SVUSs")
        }
        
        #Fusion Variants
        if (exists("fusionVariants",parsefile$results)&&length(parsefile$results$fusionVariants)!=0){
          #Remove Therapies
          for (i in 1:length(parsefile$results$fusionVariants)){
            parsefile$results$fusionVariants[[i]]$therapies=NULL
          }
          #Extract Data
          variants=list()
          for(i in 1:length(parsefile$results$fusionVariants)){
            variants[[i]]=as.data.frame(t(unlist(parsefile$results$fusionVariants[[i]])))
          }
          FVsData=do.call(rbind,variants)%>%
            mutate(mutation_type="Fusion Variants")%>%
            dplyr::rename(gene5display=gene5Display)%>%
            dplyr::rename(gene3display=gene3Display)
          rm(variants)
          if (exists("FVsData")){
            mutCount=mutCount+1
          }
        }else{
          #print("No FVs")
        }
        
        #Inherited Relevant Variants
        if (exists("inheritedRelevantVariants",parsefile$results)&&length(parsefile$results$inheritedRelevantVariants)!=0){
          #Extract Data
          if (length(parsefile$results$inheritiedIncidentalFindings)!=0){
            variants=list()
            for(i in 1:length(parsefile$results$inheritedRelevantVariants)){
              variants[[i]]=as.data.frame(t(unlist(parsefile$results$inheritedRelevantVariants[[i]])))
            }
            IRVsData=do.call(rbind,variants)%>%
              mutate(mutation_type="Inherited Relevant Variants")
            rm(variants)
            if (exists("IRVsData")){
              mutCount=mutCount+1
            }
          }else{
            #print("No IRVs")
          }
        }else{
          #print("No IRVs")
        }
        
        #Inherited Incidental Findings
        if (exists("inheritedIncidentalFindings",parsefile$results)&&length(parsefile$results$inheritedIncidentalFindings)!=0){
          #Extract Data
          if (length(parsefile$results$inheritiedIncidentalFindings)!=0){
            variants=list()
            for(i in 1:length(parsefile$results$inheritedIncidentalFindings)){
              variants[[i]]=as.data.frame(t(unlist(parsefile$results$inheritedIncidentalFindings[[i]])))
            }
            IIFsData=do.call(rbind,variants)%>%
              mutate(mutation_type="Inherited Incidental Findings")
            rm(variants)
            if (exists("IIFsData")){
              mutCount=mutCount+1
            }
          }else{
            #print("No IIFs")
          }
        }else{
          #print("No IIFs")
        }
        
        #Inherited Variants of Unknown Significance
        if (exists("inheritedVariantsOfUnknownSignificance",parsefile$results)&&length(parsefile$results$inheritedVariantsOfUnknownSignificance)!=0){
          #Extract Data
          if (length(parsefile$results$inheritiedIncidentalFindings)!=0){
            variants=list()
            for(i in 1:length(parsefile$results$inheritedVariantsOfUnknownSignificance)){
              variants[[i]]=as.data.frame(t(unlist(parsefile$results$inheritedVariantsOfUnknownSignificance[[i]])))
            }
            IVUSsData=do.call(rbind,variants)%>%
              mutate(mutation_type="Inherited Variants of Unknown Significance")
            rm(variants)
            if (exists("IVUSsData")){
              mutCount=mutCount+1
            }
          }else{
            #print("No IVUSs")
          }
        }else{
          #print("No IVUSs")
        }
        
        #Merge and Output data tables
        
        if (mutCount!=0){
          outMutations=rbind.fill(get0("SPAMsData"),get0("SPACNVsData"),get0("SBRVsData"),get0("SVUSsData"),get0("FVsData"),get0("IRVsData"),get0("IIFsData"),get0("IVUSsData"))%>%
            select(mutation_type,everything())
          outMutations=merge(orderData,outMutations)
          outMutations=merge(patientData,outMutations)
          
          i=sapply(outMutations,is.factor)
          outMutations[i]=lapply(outMutations[i],as.character)
          outMutations[is.na(outMutations)]=""
          outMutations=as.data.table(outMutations)
          
          if (grepl("prospect",getwd())){
            outMutations$TempusSource="Prospective"
          }else if(grepl("retro",getwd())){
            outMutations$TempusSource="Retrospective"
          }else{
            outMutations$TempusSource="Unknown"
          }
          
          outMutations=outMutations%>%
            mutate(Version=schema)
        
          mutationFilename=paste("~/JSON/parsedFiles/mutations/mut_",title,".rda",sep='')
          save(outMutations,file=mutationFilename,version=3)
          #write.table(outMutations,mutationFilename,sep=",",append=TRUE,row.names=FALSE)
          #assign(paste("mergedMutations_",title,sep=''),mergedMutations,envir=.GlobalEnv)
          #assign(paste("patientMutations_",title,sep=''),merge(patientData,mergedMutations),envir=.GlobalEnv)
        }else{
          print("FLAG")
        }
        
      }else{
        print("FLAG")
      } 
      
      if (patCount!=0){
        outPatients=rbind.fill(get0("patientData"),get0("orderData"),get0("reportData"),get0("tmbData"),get0("specimensData"))
        outPatients[outPatients==""]=NA
        outPatients=outPatients%>%
          fill(c(1:ncol(outPatients)),.direction=c("downup"))%>%
          unique()
        
        i=sapply(outPatients,is.factor)
        outPatients[i]=lapply(outPatients[i],as.character)
        outPatients[is.na(outPatients)]=""
        
        outPatients=as.data.table(outPatients)
        
        if (grepl("prospect",getwd())){
          outPatients$TempusSource="Prospective"
        }else if(grepl("retro",getwd())){
          outPatients$TempusSource="Retrospective"
        }else{
          outPatients$TempusSource="Unknown"
        }
        
        outPatients=outPatients%>%
          mutate(Version=schema)
        
        patientFilename=paste("~/JSON/parsedFiles/patients/pat_",title,".rda",sep='')
        save(outPatients,file=patientFilename,version=3)
        #write.table(outPatients,patientFilename,sep=",",append=TRUE,row.names=FALSE)
        
      }else{
        print("FLAG")
      }
      
}


