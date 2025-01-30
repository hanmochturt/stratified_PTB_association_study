#===============#
# SETTING UP ####
#===============#

set.seed(123)

## LIBRARIES ####
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(splines)
library(forcats)
library(lubridate)
library(table1)
library(htmltools)
library(ggplot2)

## SELECT WHAT TO PROCESS ####
calculate_crude <<- 0
calculate_adjusted <<- 0
calculate_crude_subtype <<- 0
calculate_adjusted_subtype <<- 0

## DIRECTORIES ####
base_dir <<- 'Z:\\Birth_IRB1722929/projects/ptb_diagnoses/'
defs_dir <<- paste0(base_dir, 'defs/')
data_dir<<- paste0(base_dir, 'data/')
intermediate_dir <<- paste0(base_dir, 'intermediate/')
image_dir <<- paste0(base_dir, 'images/')
output_dir <<- paste0(base_dir, 'output/')

## FUNCTIONS ####
file_select = function(myfile) {  # picks the most recent file by datestring in the filename
  files = list.files(data_dir, pattern = paste0("[0-9]{4}-[0-9]{2}-[0-9]{2}_", myfile))
  maxdate = max(do.call(c, lapply(files, function(x) ymd(gsub("[a-zA-Z._]", "", x)))))
  return(paste0(maxdate, '_', myfile))
}

image_name = function(img_name) {return(paste0(image_dir, Sys.Date(), '_', img_name, '.png'))}

csv_write = function(df, dfname) {
  write.csv(df, file = paste0(output_dir, Sys.Date(), '_', dfname, '.csv'), row.names = F)
}


phe_defs = read.csv(paste0(defs_dir, 'phecode_definitions1.2.csv'), colClasses = c('phecode'='character'))
conds = read.csv(paste0(data_dir, file_select('conditions_binary.csv')))
births = read.csv(paste0(data_dir, file_select('cohort_complete.csv')))

termbirths = births %>% filter(preterm == 0)
termids= termbirths %>% select(person_id_mom, person_id_inf)

conds = inner_join(conds, termids)


# how many people only have one condition?
cond_sums = rowSums(conds[, 3:ncol(conds)])
cond_sums = cbind(conds[,1:2], cond_sums)
one_cond = cond_sums[cond_sums$cond_sums == 1, ]
nrow(one_cond)  # 1730 people with only one condition

# what are the most common conditions?
phe_sums = colSums(conds[, 3:ncol(conds)])
phe_total = as.data.frame(phe_sums)
phe_total$phestring = row.names(phe_total)
phe_total$phecode = as.character(str_sub(phe_total$phestring, 2))
phe_total = merge(phe_defs, phe_total, by = 'phecode')
phe_total$phe_sums = as.numeric(phe_total$phe_sums)
phe_order <- phe_total[order(-phe_total$phe_sums),]
View(phe_order)

