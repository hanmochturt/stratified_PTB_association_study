library(dplyr)
library(DBI)
library(odbc)
library(glue)
library(lubridate)
library(tidyr)
library(stringr)
library(ggplot2)
library(forcats)
library(splines)

rawdf = read.csv("Z://Birth_IRB1722929/projects/ptb_diagnoses/data/2023-07-18_conditions_raw.csv")
cohort = read.csv('Z://Birth_IRB1722929/projects/ptb_diagnoses/data/2023-07-18_births.csv')
phe_defs= read.csv('Z://Birth_IRB1722929/projects/ptb_diagnoses/defs/phecode_definitions1.2.csv', colClasses = c('phecode'='character'))


cohort = cohort %>% mutate(ptb = case_when(GAWEEKS < 37 ~ 1, TRUE ~ 0),
                               edu = case_when(MATEDUC < 12 ~ 'lt12',
                                               MATEDUC > 12 ~ 'gt12',
                                               MATEDUC == 12 ~ 'eq12',
                                               MATEDUC == '17:more than college' ~ 'gt12',
                                               TRUE ~'unknown'),
                               insurance = case_when(INSPRIVATE == '1:Yes' ~ 'private',
                                                     INSPRIVATE == '0:No' ~ 'public',
                                                     TRUE ~ 'unknown'))

cohort$race9c = factor(cohort$race9c)
cohort$edu = factor(cohort$edu)
cohort$insurance = factor(cohort$insurance)

lbw = rawdf %>% filter(phecode == '288.1')


cond_occur = list(unique(lbw$condition_occurrence_id))[[1]]

pids = list(unique(rawdf$person_id))[[1]]


con_home_dir = con <- dbConnect(odbc::odbc(),
                                Driver = "ODBC Driver 17 for SQL Server",
                                Server = "QcDidDwDb001.ucsfmedicalcenter.org",
                                Trusted_Connection = "yes",
                                #Port = 1433
)

leukocyte_ct_list = c(3000905, 3030232, 3032393, 3003282, 3008511,  3010813)

sqlstring = glue::glue_sql(
  "SELECT  M.person_id, M.measurement_concept_id, M.value_as_number, M.unit_concept_id, M.measurement_date
  FROM OMOP_DEID.omop.measurement M
  WHERE M.measurement_concept_id IN ({leukocyte_ct_list*})
  AND M.person_id IN ({pids*})
  ", .con = con)
tmp = dbGetQuery(con, sqlstring)

length(unique(tmp$person_id))
labs = tmp %>% filter(!is.na(value_as_number) & 
                        !is.na(unit_concept_id) & 
                        unit_concept_id != '8792')  # 8792 units are Kelvin per microliter??
range(labs$value_as_number)
length(unique(labs$person_id))

# merge labs with pregnancy data to filter to labs prior to pregnancy
df = merge(cohort, labs, by.x = 'person_id_mom', by.y = 'person_id')

cutoff_dates = rawdf %>% dplyr::select(person_id, person_id_inf, past_cutoff) %>% distinct()
df = merge(df, cutoff_dates, by.x = c('person_id_mom', 'person_id_inf'), by. = c('person_id', 'person_id_inf'))
df = df %>% filter(measurement_date <= preg_start & measurement_date >= past_cutoff)
length(unique(df$person_id_mom))


df = df %>% mutate(P288.1 = case_when(person_id_inf %in% lbw$person_id_inf ~ 1,
                                      TRUE ~ 0))

df_birth = df %>% group_by(person_id_mom, person_id_inf) %>% summarize(meanval = mean(value_as_number),
                                                                       stdval = sd(value_as_number),
                                                                       minval = min(value_as_number),
                                                                       maxval = max(value_as_number),
                                                                       nmeas =  n())


df_preg = df %>% select(person_id_mom, person_id_inf, PRETERM, GAWEEKS, P288.1, ptb, edu, insurance, 
                        race9c, INSPRIVATE, MATEDUC, MATAGE) %>% distinct()

df_birth = merge(df_birth, df_preg, by = c('person_id_mom', 'person_id_inf'))

df_birth = df_birth %>% mutate( wbc_cat = case_when(meanval < 4.5 ~ 'low',
                                                    meanval > 11 ~ 'high',                                                    
                                                    TRUE ~ 'normal'),
                                wbc_extreme = case_when(minval < 4.5 & maxval > 11 ~ 'high_low',
                                                        minval < 4.5 ~ 'low',
                                                        maxval > 11 ~ 'high',
                                                        TRUE ~ 'normal'))
df_birth$wbc_cat = fct_relevel(df_birth$wbc_cat, c('normal', 'low', 'high'))
df_birth$wbc_extreme = fct_relevel(df_birth$wbc_extreme, c('normal', 'low', 'high', 'high_low'))

all_cont = glm(ptb ~ ns(meanval) + insurance + edu + race9c + MATAGE, data = df_birth, family = 'binomial')
summary(all_cont)

all_cont_crude = glm(ptb ~ ns(meanval), data = df_birth, family = 'binomial')
summary(all_cont_crude)


all_cat = glm(ptb ~ wbc_cat, data = df_birth, family = 'binomial')
summary(all_cat)

df_low = df_birth %>% filter(P288.1 == 1)
mean(df_low$minval)
sd(df_low$minval)

all_extreme = glm(ptb ~ wbc_extreme + insurance + edu + race9c + MATAGE, data = df_birth, family = 'binomial')
summary(all_extreme)
# alternative model examining extremes

