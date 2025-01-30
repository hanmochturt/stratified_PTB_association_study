# this is the pipeline for the PTB phecodes project

library(DBI)
library(odbc)
library(glue)
library(lubridate)
library(dplyr)
library(tidyr)

# GLOBALS AND PRESETS ####
from_scratch <<- 0  # indicates whether to run the whole program from scratch
first_recorded <<- 0
primagravid <<- 0
all_pregnancies <<- 1
source_icd <<- 1
source_snomed <<- 0


file_select = function(myfile) {
  # this function picks the most recent file by datestring in the filename
  files = list.files(data_dir, pattern = paste0("[0-9]{4}-[0-9]{2}-[0-9]{2}_", myfile))
  maxdate = max(do.call(c, lapply(files, function(x) ymd(gsub("[a-zA-Z._]", "", x)))))
  return(paste0(maxdate, '_', myfile, '.csv'))
}


## directories ####
base_dir <<- 'Z:\\Birth_IRB1722929/projects/ptb_diagnoses/'
defs_dir <<- paste0(base_dir, 'defs/')
data_dir<<- paste0(base_dir, 'data/')
#intermediate_dir <<- paste0(base_dir, 'intermediate/')
image_dir <<- paste0(base_dir, 'images/')
output_dir <<- paste0(base_dir, 'output/')


csv_write = function(df, dfname) {
  write.csv(df, file = paste0(data_dir, Sys.Date(), '_', dfname, '.csv'), row.names = F)
}

# READ IN COHORT FILES ####
births = read.csv(paste0(data_dir, file_select('births')))
conditions_raw = read.csv(paste0(data_dir, file_select('conditions_raw')))
conditions_dates = read.csv(paste0(data_dir, file_select('conditions_dates')))

# read map
mymap = read.csv(paste0(defs_dir, 'pdb_cols_rename.csv'))

# RENAME PDB COLUMNS

rename_fn = function(df, mapping) {
  cnew = c()
  for (i in colnames(df)){
    if(i %in% mapping$pdb_name){
      tmp = mapping[mapping$pdb_name == i,  ]$new_name
    } else {
      tmp = i
    }
    cnew = c(cnew, tmp)
  }
  colnames(df) = cnew
  return(df)
}

births = rename_fn(births, mymap)
conditions_raw = rename_fn(conditions_raw, mymap)
conditions_dates = rename_fn(conditions_dates, mymap)


#################################
# WRITE FILES TO BE ANALYZED ####
#################################
csv_write(conditions_raw, 'conditions_raw')
csv_write(conds, 'conditions_dates')
csv_write(df, 'births')

