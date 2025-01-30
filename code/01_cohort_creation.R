pdb_write = pdb %>% dplyr::select( MATUNUM, NBUNUM, DELDATE, MATAGE, RACE, race9c, MATDOB, WOUTCOME,
                                  GAWEEKS, GRAVIDA, PARA, LIVEBIRTH, MATEDUC,
                                  DELDATE, INSURANCE, INSPRIVATE, INSMEDICAL, INSMEDICARE, INSKAISER,
                                  INSOTHER, INSURANCESPEC, PRETERM, R_HYPERTENSION, preeclampsia, DIABETES, 
                                  MATDEATH, VAGDEL, ASTHMA, ANXIETY, LUNG, BIPOLAR, HELLP, PREVPRETERM, NBDEL, SINGLETON)


mapdf = mapdf %>% filter(!is.na(PDK_mom) & !is.na(PDK_inf) & PDK_mom != "" & PDK_inf != "")

df_write = merge(mapdf, pdb_write, by = c('MATUNUM', "NBUNUM"))

df_write = df_write  %>% filter(!is.na(PDK_mom) & !is.na(PDK_inf) & PDK_mom != "" & PDK_inf != "")

# Recovering birthdates for mom and baby from Limited Data Set (LDS)
mdob = df_write[is.na(df_write$MATDOB),]

epic_mom = unique(list(mdob$EpicID_mom)[[1]])

# nicknames for columns are used instead of their real names because their real
# names may be protected by Epic EHR company patents, preventing us from publicly 
# publishing them. If you are a UCSF affiliate with CDW access and would like 
# the real string that we used, please email us to request access
sql_query = "SELECT PatientId, Key, DateOfBirth
    FROM CDW_NEW.LDS.Patient
    WHERE PatientId IN ({id_mom*})" 

sqlstring = glue::glue_sql(sql_query, .con = con)
mom_lds = dbGetQuery(con, sqlstring)
mom_lds = mom_lds %>% distinct()

df_write = left_join(df_write, mom_lds, by = c('id_mom' = 'PatientId'))
df_write = df_write %>% mutate(MATDOB  = case_when(is.na(MATDOB) ~ DateOfBirth,
                                                   TRUE ~ MATDOB)) %>%
  dplyr::select(-c( DateOfBirth, Key))

library(lubridate)
df_write$daydiff = ymd(df_write$MATDOB) - ymd(df_write$dob_mom)
df_write$DELDATE_ehr = ymd(df_write$DELDATE) - df_write$daydiff

df_write = df_write %>% dplyr::select(-c( daydiff, MATDOB, DELDATE, MATUNUM, NBUNUM, pdb_mrn, pdb_nbu, new_mrn, new_nbu, mrn_to_map, nbu_to_map))

write.csv(df_write, file = paste0(paste0(fp_base, 'projects/ptb_diagnoses/data/'),  Sys.Date(), '_pdk_from_pdb', '.csv'))
