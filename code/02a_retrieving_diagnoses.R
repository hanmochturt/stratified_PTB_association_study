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
all_pregnancies <<- 0
source_icd <<- 1
source_snomed <<- 0


file_select = function(myfile) {
  # this function picks the most recent file by datestring in the filename
  files = list.files(data_dir, pattern = paste0("[0-9]{4}-[0-9]{2}-[0-9]{2}_", myfile))
  maxdate = max(do.call(c, lapply(files, function(x) ymd(gsub("[a-zA-Z._]", "", x)))))
  return(paste0(maxdate, '_', myfile, '.csv'))
}

rerun = function(myfile) {
  if (from_scratch == 1 | !file.exists(myfile)) {
    return(1)
  } else {
    return(0)
  }

}

conditions_to_home_db <- function(df){

  con <- dbConnect(odbc::odbc(),
                   Driver = "ODBC Driver 17 for SQL Server",
                   #Database = "omop_deid",
                   Server = "QcDidDwDb001.ucsfmedicalcenter.org",
                   Trusted_Connection = "yes",
                   #Port = 1433
  )

  pids = unique(list(df$person_id_mom)[[1]])

  dbRemoveTable(con_home_dir, 'ptbd_condition_occurrence', fail_if_missing = F)
  sqlstring = glue_sql(
    "SELECT CO.person_id, CO.condition_occurrence_id, CO.condition_concept_id, CO.condition_start_date, CO.condition_end_date,
      CO.condition_type_concept_id, CO.condition_source_value
    INTO home_jcostello.dbo.ptbd_condition_occurrence
    FROM omop_deid.omop.condition_occurrence CO
    WHERE CO.person_id in ({pids*})", .con = con)
  dbExecute(con, sqlstring)
}

## directories ####
base_dir <<- 'Z:\\Birth_IRB1722929/projects/ptb_diagnoses/'
defs_dir <<- paste0(base_dir, 'defs/')
data_dir<<- paste0(base_dir, 'data/')
#intermediate_dir <<- paste0(base_dir, 'intermediate/')
image_dir <<- paste0(base_dir, 'images/')
output_dir <<- paste0(base_dir, 'output/')

image_name = function(img_name) {return(paste0(image_dir, Sys.Date(), '_', img_name, '.png'))}
csv_write = function(df, dfname) {
  write.csv(df, file = paste0(data_dir, Sys.Date(), '_', dfname, '.csv'), row.names = F)
}

# READ IN COHORT ####
cohort_full = read.csv(paste0(data_dir, file_select('pdk_from_pdb')))


df = cohort_full
df  = df %>% rename(DELDATE = DELDATE_ehr)
df <- df %>% dplyr::select(-c('PDK_inf', "X", 'LIVEBIRTH')) %>%
  filter(!is.na(PDK_mom) & !is.na(person_id_mom) & !is.na(person_id_inf) & !is.na(GAWEEKS) &
           SINGLETON == "1: Yes" & !is.na(DELDATE)) %>% distinct()

# filter pregnancies
if (first_recorded == 1){
  filter_df = df %>% dplyr::select('person_id_mom', 'DELDATE') %>%
    group_by(person_id_mom) %>%
    summarise(DELDATE = min(DELDATE))
  df = merge(filter_df, df, by =c('person_id_mom', 'DELDATE'))
} else if (primagravid == 1) {
  filter_df = df %>% filter(PARA == 0 & GRAVIDA == 1) %>% dplyr::select(person_id_mom, DELDATE)
  df = merge(filter_df, df, by =c('person_id_mom', 'DELDATE'))
}


#  remove records with multiple accounts for the same delivery
#   these are all singleton deliveries
non_duplicated = df %>% group_by(PDK_mom, DELDATE) %>%
  summarize(countid=n()) %>% filter(countid == 1) %>% dplyr::select(-countid)

df <- merge(df, non_duplicated, by = c('PDK_mom', 'DELDATE'))

# remove records where the same infant id is assigned to multiple maternal ids
non_duplicated = df %>% group_by(person_id_inf) %>% 
  summarise(countid = n()) %>% filter(countid == 1) %>% dplyr::select(-countid)

df = merge(df, non_duplicated, by = 'person_id_inf')

# only consider diagnoses before conception
df$preg_start = ymd(df$DELDATE) - weeks(df$GAWEEKS)

con_home_dir = con <- dbConnect(odbc::odbc(),
                                Driver = "ODBC Driver 18 for SQL Server",
                                Database = "home_jcostello",
                                Server = "QcDidDwDb001.ucsfmedicalcenter.org",
                                Trusted_Connection = "yes",
                                TrustServerCertificate="yes",
                                #Port = 1433
)


dbWriteTable(con_home_dir, "ptbd_pregnancy", df, overwrite= T)


#####################
# GET CONDITIONS ####
#####################


if (from_scratch) {  # this requires a large database query - run only if necessary
  conditions_to_home_db(df)
}

conditions = dbGetQuery(con_home_dir, 'SELECT o.*, c.* FROM home_jcostello.dbo.ptbd_condition_occurrence o
                        INNER JOIN OMOP_DEID.omop.concept c
                        ON o.condition_source_value = concept_code')

conditions = merge(conditions, df, by.x = 'person_id', by.y = 'person_id_mom')
conditions = conditions[conditions$condition_start_date < conditions$preg_start, ]


if (all_pregnancies == 1) { # date filtering on conditions
  # identify earliest recorded pregnancy per person
  prevstart = df %>% dplyr::select('person_id_mom', 'person_id_inf', 'DELDATE') %>%
    group_by(person_id_mom) %>%
    arrange(DELDATE) %>%
    mutate(row_rank = 1:n(), prev_row = row_rank - 1) %>%
    left_join(., ., by = c('prev_row' = 'row_rank', 'person_id_mom')) %>%
    mutate(prev_start_date = case_when(!is.na(DELDATE.y) ~ ymd(DELDATE.y), TRUE ~ ymd('1900-1-1')),
           past_cutoff = prev_start_date + months(6),
           person_id_inf = person_id_inf.x) %>%
    dplyr::select(-c('row_rank', 'prev_row', 'person_id_inf.y', 'DELDATE.y', 'prev_row.y', 'DELDATE.x',
                     'person_id_inf.x'))
  conditions = merge(conditions, prevstart, by.x = c('person_id', 'person_id_inf'),
                     by.y = c('person_id_mom', 'person_id_inf'))

  conditions = conditions[conditions$condition_start_date >= conditions$past_cutoff, ]

  cond_unique = conditions %>%  # consider conditions as binary, rather than count data
    dplyr::select(person_id, person_id_inf, condition_concept_id, condition_source_value) %>%
    distinct()

  cond_per_person = conditions %>% group_by(person_id, person_id_inf, condition_concept_id, condition_source_value) %>%
    summarize(earliest_start = min(condition_start_date),
              latest_start = max(condition_start_date))

} else {
  cond_unique = conditions %>%  # consider conditions as binary, rather than count data
    dplyr::select(person_id, condition_concept_id, condition_source_value) %>%
    distinct()

  cond_per_person = conditions %>% group_by(person_id, condition_concept_id, condition_source_value) %>%
    summarize(earliest_start = min(condition_start_date),
              latest_start = max(condition_start_date))
}



########################################
# CREATE CONDITIONS TO PHECODES MAP ####
########################################

conditions_to_map = conditions %>% dplyr::select(condition_concept_id, condition_source_value) %>%
  distinct()

sqlstring = "select * FROM OMOP_DEID.omop.concept
  WHERE vocabulary_id IN ('ICD10CM', 'ICD9CM', 'ICD9', 'ICD10')"
cond_concepts = dbGetQuery(con, sqlstring)


phe_defs = read.csv(paste0(defs_dir, 'phecode_definitions1.2.csv'), colClasses = c('phecode'='character'))
#phe_defs = phe_defs %>% filter(category_number != 12)
icd9_to_phecode = read.csv(paste0(defs_dir, 'phecode_icd9_rolled.csv'),  colClasses = c('PheCode'='character'))
icd10_to_phecode = read.csv(paste0(defs_dir, 'Phecode_map_v1_2_icd10cm_beta.csv'),  colClasses = c('phecode'='character'))

cond_concepts = cond_concepts %>% dplyr::select(concept_id, concept_code, vocabulary_id)

concepts9 = cond_concepts %>% filter(vocabulary_id %in% c('ICD9', 'ICD9CM'))
concepts10 = cond_concepts %>% filter(vocabulary_id %in% c('ICD10', 'ICD10CM'))

phecols = c('ICD', 'ICD_str', 'phecode', 'phenotype', 'exlude_codes', 'exclude_pheno', 'rollup',
            'leaf', 'concept_id', 'vocabulary_id')

icd9map = left_join(icd9_to_phecode, concepts9, by = c("ICD9" = "concept_code" )) %>%
  dplyr::select(-Ignore.Bool)
colnames(icd9map) = phecols

icd10map = left_join(icd10_to_phecode, concepts10, by = c('icd10cm' = 'concept_code'))
icd10map = icd10map[, c(1:6, 8, 7, 9, 10)]
colnames(icd10map) = phecols

phecode_map = rbind(icd10map, icd9map)

phecodes_map_simple = phecode_map %>% dplyr::select(ICD, phecode) %>%
  mutate(phecode = case_when(nchar(phecode) == 6 ~ substr(phecode, 1, 5),
                             TRUE ~ phecode))

phecodes_map = merge(phecode_map, phe_defs, by = 'phecode')
phecodes_map_simple = merge(phecodes_map_simple, phe_defs, by = 'phecode')

mapped = merge(conditions_to_map, phecodes_map_simple, by.x = 'condition_source_value', by.y = 'ICD')
mapped_full = merge(conditions_to_map, phecode_map, by.x = 'condition_source_value', by.y = 'ICD')


#################################
# MAP CONDITIONS TO PHECODES ####
#################################

conds = merge(cond_per_person, df, by.x = c('person_id', 'person_id_inf'), by.y = c('person_id_mom', 'person_id_inf'))
conds = merge(conds, mapped, by = 'condition_source_value')

conds$phestring = paste0('P', conds$phecode)

conds_binary = conds %>%
  dplyr::select(person_id, person_id_inf, phestring) %>%
  distinct() %>% mutate(presence = 1) %>%
  pivot_wider(names_from = phestring, values_from = presence, values_fill = 0)

conds_f = merge(cond_per_person, df, by.x = c('person_id', 'person_id_inf'), by.y = c('person_id_mom', 'person_id_inf'))
conds_f = merge(conds_f, mapped_full, by = 'condition_source_value')
conds_f$phestring = paste0('P', conds_f$phecode)

conds_binary_f = conds_f %>%
  dplyr::select(person_id, person_id_inf, phestring) %>%
  distinct() %>% mutate(presence = 1) %>%
  pivot_wider(names_from = phestring, values_from = presence, values_fill = 0)

conditions_raw = merge(conditions, mapped, by = 'condition_source_value')
conditions_raw = merge(conditions_raw, df[c('person_id_inf', 'person_id_mom')], 
                       by.x = c('person_id', 'person_id_inf'), by.y = c('person_id_mom', 'person_id_inf'))


#################################
# WRITE FILES TO BE ANALYZED ####
#################################
csv_write(conditions_raw, 'conditions_raw')
csv_write(conds_binary, 'conditions_binary')
csv_write(conds_binary_f, 'conditions_binary_fullphecode')
csv_write(conds, 'conditions_dates')
csv_write(df, 'births')

