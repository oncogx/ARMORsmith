# ARMORsmith

ARMORsmith is a suite of scripts to assemble the data from the UC San Diego Aggregated Registry for Molecular Oncology Research (ARMOR). 

This work is still in progress. The workflow is not fully automatized and depends highly on local infrastructure and source data. The overall workflow is represented in Figure 1 and consist of the following major steps:

* Step 1: Parses the molecular reports and generate list of Medical Record Numbers for CTRI DECS request (clinical data)
* Step 2: Parses the metadata tables, aggregate them and de-identifies. Generate the master key table. Fixes the data MRN (typos and alternate MRNs) and de-identifies. 
* Step 3: Reformat the data table into cbioportal flat files, including harmonizing the mutation description and effects. 
* Step 4: Merge the cbioportal formatted tables into the flat files requeired for cbioportal ingestion
