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

## SELECT WHAT COHORT ####
diag_req <<- 1

## SELECT WHAT TO PROCESS ####
calculate_crude <<- 1  # overall odds
calculate_adjusted <<- 0
calculate_crude_subtype <<- 1  # spont and indicated odds
calculate_adjusted_subtype <<- 0

## DIRECTORIES ####
base_dir <<- 'Z:\\Birth_IRB1722929/projects/ptb_diagnoses/'
defs_dir <<- paste0(base_dir, 'defs/')
data_dir<<- paste0(base_dir, 'data/')
intermediate_dir <<- paste0(base_dir, 'intermediate/')

if (diag_req == 1){
  image_dir <<- paste0(base_dir, 'images/01_revisions/')
  output_dir <<- paste0(base_dir, 'output/01_revisions/')
} else {
  image_dir <<- paste0(base_dir, 'images/01_revisions/nodiagreq/')
  output_dir <<- paste0(base_dir, 'output/01_revisions/nodiagreq/')
}

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


crude_model = function(){
  cr_model = mapply(function(v1, v2) broom::tidy(glm(as.formula(paste0("preterm ~ ", v1)), 
                                                      na.action = na.omit, family = 'binomial', 
                                                      data = df[, colnames(df) %in% c(colnames(df)[1:5], v2)]))[1,5], 
                     colnames(df)[2:ncol(df)], colnames(df)[2:ncol(df)])
  return(data.frame(adj_model))
}

adjusted_model = function(df) {
  adj_model = mapply(function(v1, v2) broom::tidy(glm(as.formula(paste("preterm ~ MATEDUC + INSPRIVATE + ns(MATAGE, df=3)", v1, sep = " + ")), 
                                              na.action = na.omit, family = 'binomial', 
                                              data = df[, colnames(df) %in% c(colnames(df)[1:5], v2)]))[10,5], 
             colnames(df)[6:ncol(df)], colnames(df)[6:ncol(df)])
  return(data.frame(adj_model))
}


odds_ratios_calculate = function(df, remove_columns, cohort_name, adjustment){

  df <- df %>% dplyr::select(-all_of(remove_columns))

  if (adjustment == 'crude') {
    or_tmp = vector('list', length = ncol(df) - 1)
    for(i in 1:length(colnames(df[,-1]))) {  # this step is slow
      print(i)
      
      count_pt = nrow(df[df[,i+1] == 1 & df$preterm == 1,])
      count_t = nrow(df[df[,i+1] == 1 & df$preterm == 0,])
      
      mymodel = glm(preterm ~ df[[i + 1]], data = df, family = 'binomial')
      
      tmp = cbind(exp(coef(mymodel)), exp(confint.default(mymodel)), coef(summary(mymodel))[,4])[2,]
      names(tmp) = c('or', 'or_lower', 'or_upper', 'pval')
      
      tmpct = c('count_preterm' = round(count_pt,0), 'count_term' = round(count_t,0))
      tmp = c(tmp, tmpct)
      
      or_tmp[[i]] = tmp
    }
    ordf = do.call('rbind', or_tmp)
    phestring = colnames(df[-1])
    ordf = as.data.frame(cbind(ordf, phestring))
  } else {
    or_tmp = vector('list', length = ncol(df) - 5)
    upperlim = length(colnames(df)) - 5
    for(i in 1:upperlim) {
      print(i)
      count_pt = nrow(df[df[,i+5] %in% c(1) & df$preterm == 1,])
      count_t = nrow(df[df[,i+5] %in% c(1) & df$preterm == 0,])
      
      mymodel = glm(preterm ~ df[[i + 5]] + INSPRIVATE + MATEDUC + ns(MATAGE, df=3), data = df, family = 'binomial')
      tmp = cbind(exp(coef(mymodel)), exp(confint.default(mymodel)), coef(summary(mymodel))[,4])[2,]
      names(tmp) = c('or', 'or_lower', 'or_upper', 'pval')
      
      tmpct = c('count_preterm' = round(count_pt,0), 'count_term' = round(count_t,0))
      tmp = c(tmp, tmpct)
      
      or_tmp[[i]] = tmp
    }
    ordf = do.call('rbind', or_tmp)
    phestring = colnames(df[-seq(1:5)])
    ordf = as.data.frame(cbind(ordf, phestring))
  }

  ordf <- ordf %>% mutate(phecode = str_sub(phestring, 2))
  ordf = merge(ordf, phe_defs, by = 'phecode')
  ordf = ordf %>% mutate(pval = as.numeric(pval),
                         or = as.numeric(or),
                         or_lower = as.numeric(or_lower),
                         or_upper = as.numeric(or_upper),
                         pval_minus_log = -log10(pval))

  ordf = ordf[ordf$pval_minus_log != Inf,]

  ordf = ordf[order(ordf$pval_minus_log),]


  pval_bh = p.adjust(ordf$pval, method = 'BH', n = nrow(ordf))
  ordf = cbind(ordf, pval_bh)
  ordf$padj_minus_log = -log10(ordf$pval_bh)

  csv_write(ordf, paste0(cohort_name, '_', adjustment, '_', 'odds_ratios'))

  return(ordf)
}


permute_drop_test = function(df, remove_columns, cohort_name, adjustment){
  df <- df %>% dplyr::select(-all_of(remove_columns))
  threshold = 50
  tdf = df[, 6:ncol(df)]
  tdf = tdf[, colSums(tdf) > 1]
  
  df = cbind(df[,1:5], tdf)
  
  drop_list = vector('list', length = ncol(df) - 5)
  for (i in 6:ncol(df)) {
    positives = which(df[, i] == 1)
    if(sum(df[, i]) <= threshold) {
      tmp = rep(positives, length.out = threshold)
    } else{
      tmp = sample(positives, size = threshold, replace = FALSE)
    }
    drop_list[[i-5]] = tmp
  }
  
  dropdf = data.frame(drop_list)
  colnames(dropdf) = sapply(seq(6, ncol(df)), function(x, y) paste0(y, x), y = 'col')
  
  ptm = proc.time()
  pvals = vector('list', length = threshold)
  for (i in 1:threshold) {
    df_na = data.frame(df)
    for (j in 6:ncol(df_na)){
      df_na[dropdf[i, j-5], j] = NA
    }
    pvals[[i]] = adjusted_model(df_na)
  }
  b = proc.time() - ptm
  pvaldf = do.call(rbind, pvals)
  #return(pvaldf)
  colnames(pvaldf) = sapply(colnames(pvaldf), function(x) gsub('.p.value', '', x))
  #pvaldf = pvaldf %>% dplyr::select(-c('early', 'late'))
  csv_write(pvaldf, paste0(cohort_name, '_', adjustment, '_droptest_50'))
  
  pvdf = pvaldf
  cnames = colnames(pvdf)
  cnames = lapply(cnames, function(x) gsub('.p.value', '', x))
  colnames(pvdf) = cnames
  
  pertbh = vector('list', length = nrow(pvdf))  # BH adjusted pvalues
  for (i in 1:nrow(pvdf)){
    pval_bh = p.adjust(pvdf[i,], method = 'BH', n = ncol(pvdf))
    pertbh[[i]] = pval_bh
  }
  
  pertbh = data.frame(pertbh)
  colnames(pertbh) = letters[1:ncol(pertbh)]
  
  pertbh = data.frame(t(pertbh))
  
  pertrow = vector('list', length = ncol(pertbh))
  for (i in 1:ncol(pertbh)){
    pertrow[[i]] = sum(pertbh[,i] < .05)
  }
  
  pertbh = rbind(pertbh, pertrow)
  pertbh = data.frame(t(pertbh))
  pertbh = pertbh %>% rownames_to_column(var = 'phestring')
  pertbh$sig_pct = pertbh$X1/threshold * 100
  pertbh$phecode = gsub('P', '', pertbh$phestring)
  
  pertpltdf = merge(pertbh, phe_defs, by = 'phecode')  # add in phecode descriptions
  
  #pertpltdf = pertpltdf %>% filter(category_number != 12)
  
  pertpltdf = pertpltdf %>% filter(X1 > 1)
  
  # create factor to order plot
  pertpltdf = pertpltdf[order(pertpltdf$sig_pct),]
  rownames(pertpltdf) = NULL  # renumber rows
  myorder = pertpltdf$phenotype  
  pertpltdf$phenotype = as.factor(pertpltdf$phenotype)
  pertpltdf$phenotype = fct_relevel(pertpltdf$phenotype, myorder)  # relevel factor by row order
  csv_write(pertpltdf, paste0(cohort_name, '_', adjustment, '_droptest_50'))
  
  prtplt = ggplot(aes(x = sig_pct, y= phenotype), data = pertpltdf) + 
    geom_point() +
    #geom_vline(xintercept = 95) + 
    ggtitle('Indicated drop test') + 
    theme_bw() 
  
  prtplt
  ggsave(prtplt, width = 12, height = 6, 
         filename = image_name(paste0(cohort_name, '_', adjustment, '_droptest_50_summary')))
}


#====================#
# DATA PROCESSING ####
#====================#

## read in data ####

phe_defs = read.csv(paste0(defs_dir, 'phecode_definitions1.2.csv'), colClasses = c('phecode'='character'))
conds = read.csv(paste0(data_dir, file_select('conditions_binary.csv')))
births = read.csv(paste0(data_dir, file_select('births.csv')))
conds_dates = read.csv(paste0(data_dir, file_select('conditions_dates.csv')))

remove_cols = paste0('P', phe_defs[phe_defs$category_number == 12, ]$phecode)

conds_dates = conds_dates %>% dplyr::select(-condition_concept_id.y) %>% distinct()

# define preterm
births = births %>% dplyr::select(-RACE) %>% mutate(preterm = case_when(GAWEEKS < 37 ~ 1, TRUE ~ 0))
births_columns = colnames(births)[colnames(births) != 'person_id' & colnames(births) != 'preterm']

if (diag_req == 1){
# this restricts to a condition existing prior to pregnancy
  cohort = merge(births, conds, by.x = c('person_id_mom', 'person_id_inf'), by.y = c('person_id', 'person_id_inf'))
} else {
  #if no diagnosis is required prior to pregnancy
  cohort = left_join(births, conds, by = c('person_id_mom' = 'person_id', 'person_id_inf' = 'person_id_inf'))
}

cohort$MATEDUC = as.numeric(substr(cohort$MATEDUC, 1, 2))

cohort = cohort %>% filter(!is.na(MATAGE)) %>%
  mutate(MATEDUC = case_when(MATEDUC > 12 ~ 'gt12',
                                               MATEDUC == 12 ~ 'e12',
                                               MATEDUC < 12 ~ 'lt12',
                                               TRUE ~ 'unknown'))

cohort$INSPRIVATE = as.factor(cohort$INSPRIVATE)
cohort$MATEDUC = as.factor(cohort$MATEDUC)

cohort = cohort %>% mutate(race9c = case_when(race9c == '5:Single race-Native' ~ '8:Other race',  # no category with only one
                                              TRUE ~ race9c))
cohort$race9c = as.factor(cohort$race9c)

cohort = cohort %>% mutate(early = case_when(GAWEEKS < 32 ~ 1, TRUE ~ 0),
                           late = case_when(GAWEEKS < 37 & GAWEEKS >= 32 ~ 1, TRUE ~ 0))

cohort = cohort %>% mutate(INSPRIVATE = fct_na_value_to_level(INSPRIVATE, 'unknown'),
                           MATEDUC = fct_na_value_to_level(MATEDUC, 'unknown'))
cohort[is.na(cohort)] <- 0

## TABLE 1 ####

tabdf = cohort %>% mutate(Delivery = case_when(preterm == 0 ~ "Term", TRUE ~ "Preterm")) %>%
  mutate(Delivery = factor(Delivery),
         Race = factor(race9c),
         Maternal.Age = MATAGE,
         Private.Insurance = INSPRIVATE,
         Maternal.Education = case_when(MATEDUC == 'e12' ~ '12 years',
                                        MATEDUC == 'lt12' ~ '<12 years',
                                        TRUE ~ '>12 years'))



mytab = table1(~ Race + Maternal.Age + Private.Insurance + Maternal.Education | Delivery, data = tabdf)
save_html(mytab, file=paste0(output_dir,'table1.html'), libdir = 'lib', lang = 'en', background = 'white')

tabdf2 = cohort %>% mutate(Delivery = case_when(PRETERM %in% c('4:medically indicated', '5:Termination Iatrogenic') ~ "Indicated",
                                                PRETERM %in% c("1:spontaneous",  "2:PPROM", "6:PTL with TOCO and TERM ") ~ 'Spontaneous')) %>%
  filter(Delivery %in% c('Indicated', 'Spontaneous')) %>%
  mutate(Delivery = factor(Delivery),
         Race = factor(race9c),
         Maternal.Age = MATAGE,
         Private.Insurance = INSPRIVATE,
         Maternal.Education = case_when(MATEDUC == 'e12' ~ '12 years',
                                        MATEDUC == 'lt12' ~ '<12 years',
                                        TRUE ~ '>12 years'))

mytab = table1(~ Race + Maternal.Age + Private.Insurance + Maternal.Education | Delivery, data = tabdf2)
save_html(mytab, file=paste0(output_dir,'table1a.html'), libdir = 'lib', lang = 'en', background = 'white')
cohort = cohort %>% select(-any_of(remove_cols))

indic_df = cohort %>% filter(PRETERM %in% c('0:No','4:medically indicated', '5:Termination Iatrogenic'))
spont_df = cohort %>% filter(PRETERM %in% c('0:No', "1:spontaneous",  "2:PPROM", "6:PTL with TOCO and TERM "))

adj_col_remove = births_columns[!(births_columns %in% c('MATEDUC', 'INSPRIVATE', 'race9c', 'MATAGE'))]

# write cohort files
csv_write(cohort, 'cohort_complete')
csv_write(indic_df, 'cohort_indicated')
csv_write(spont_df, 'cohort_spontaneous')

#=================#
# MAIN ANALYIS ####
#=================#

## ODDS RATIOS ####

if(calculate_crude == 1){
  crude_odds = odds_ratios_calculate(cohort, births_columns, 'overall', 'crude')
}


if(calculate_adjusted == 1){
  adj_odds = odds_ratios_calculate(cohort, adj_col_remove, 'overall', 'adj')
}

if(calculate_crude_subtype == 1){
  indic_crude = odds_ratios_calculate(indic_df, births_columns, 'indicated', 'crude')
  spont_crude = odds_ratios_calculate(spont_df, births_columns, 'spontaneous', 'crude')
}

if(calculate_adjusted_subtype == 1){
  indic_adj = odds_ratios_calculate(indic_df, adj_col_remove, 'indicated', 'adj')
  spont_adj = odds_ratios_calculate(spont_df, adj_col_remove, 'spontaneous', 'adj')
}


#========================#
# PERMUTATION TESTING ####
#========================#

permute_drop_test(spont_df, adj_col_remove, 'spontaneous', 'adj')
permute_drop_test(indic_df, adj_col_remove, 'indicated', 'adj')
permute_drop_test(cohort, adj_col_remove, 'overall', 'adj')

permute_drop_test(spont_df, births_columns, 'spontaneous', 'crude')
permute_drop_test(indic_df, births_columns, 'indicated', 'crude')
permute_drop_test(cohort, births_columns, 'overall', 'crude')

## DROPPING INDIVIDUALS ####
cohort_name = 'indicated'
odf = indic_df %>% dplyr::select(-all_of(adj_col_remove))
df = data.frame(odf)
threshold = 50

tdf = df[, 6:ncol(df)]
tdf = tdf[, colSums(tdf) > 1]

df = cbind(df[,1:5], tdf)
cols = ncol(df)

drop_list = vector('list', length = cols - 5)
for (i in 6:cols) {
  positives = which(df[, i] == 1)
  if(sum(df[, i]) <= threshold) {
    tmp = rep(positives, length.out = threshold)
  } else{
    tmp = sample(positives, size = threshold, replace = FALSE)
  }
  drop_list[[i-5]] = tmp
}

dropdf = data.frame(drop_list)
colnames(dropdf) = sapply(seq(6, cols), function(x, y) paste0(y, x), y = 'col')

ptm = proc.time()
pvals = vector('list', length = threshold)
for (i in 1:threshold) {
  df_na = data.frame(df)
  for (j in 6:cols){
    df_na[dropdf[i, j-5], j] = NA
  }
  pvals[[i]] = adjusted_model(df_na)
}
b = proc.time() - ptm
b
pvaldf = do.call(rbind, pvals)
colnames(pvaldf) = sapply(colnames(pvaldf), function(x) gsub('.p.value', '', x))
#pvaldf = pvaldf %>% dplyr::select(-c('early', 'late'))
csv_write(pvaldf, paste0(cohort_name, '_droptest_50'))

pvdf = pvaldf
cnames = colnames(pvdf)
cnames = lapply(cnames, function(x) gsub('.p.value', '', x))
colnames(pvdf) = cnames

pertbh = vector('list', length = nrow(pvdf))  # BH adjusted pvalues
for (i in 1:nrow(pvdf)){
  pval_bh = p.adjust(pvdf[i,], method = 'BH', n = ncol(pvdf))
  pertbh[[i]] = pval_bh
}

pertbh = data.frame(pertbh)
colnames(pertbh) = letters[1:ncol(pertbh)]

pertbh = data.frame(t(pertbh))

pertrow = vector('list', length = ncol(pertbh))
for (i in 1:ncol(pertbh)){
  pertrow[[i]] = sum(pertbh[,i] < .05)
}

pertbh = rbind(pertbh, pertrow)
pertbh = data.frame(t(pertbh))
pertbh = pertbh %>% rownames_to_column(var = 'phestring')
pertbh$sig_pct = pertbh$X1/threshold * 100
pertbh$phecode = gsub('P', '', pertbh$phestring)

pertpltdf = merge(pertbh, phe_defs, by = 'phecode')  # add in phecode descriptions

#pertpltdf = pertpltdf %>% filter(category_number != 12)

pertpltdf = pertpltdf %>% filter(X1 > 1)

# create factor to order plot
pertpltdf = pertpltdf[order(pertpltdf$sig_pct),]
rownames(pertpltdf) = NULL  # renumber rows
myorder = pertpltdf$phenotype  
pertpltdf$phenotype = as.factor(pertpltdf$phenotype)
pertpltdf$phenotype = fct_relevel(pertpltdf$phenotype, myorder)  # relevel factor by row order
csv_write(pertpltdf, paste0(cohort_name, '_drop_50_summary'))

prtplt = ggplot(aes(x = sig_pct, y= phenotype), data = pertpltdf) + 
  geom_point() +
  #geom_vline(xintercept = 95) + 
  ggtitle('Indicated drop test') + 
  theme_bw() 

prtplt
ggsave(prtplt, width = 12, height = 6, 
       filename = image_name(paste0(cohort_name, '_drop_test')))


# ARCHIVE ####

# random seed testing for significance

myseeds = round(runif(2, 100, 2^16), 0)
sigtest_seed_df = vector('list', length = length(myseeds))
idf = indic_df %>% dplyr::select(-all_of(adj_col_remove))
for (i in 1:length(myseeds)){
  set.seed(myseeds[[i]])
  print(myseeds[[i]])
  sigtest_seed_df[[i]] = sapply(colnames(idf)[6:ncol(idf)],
                                function(x, d) 
                                  broom::tidy(glm(as.formula(paste("preterm ~ race9c + MATEDUC + INSPRIVATE + ns(MATAGE, df=3)", x, sep = " + ")), family = 'binomial', data = d))[18,4],
                                d = idf)
  set.seed(123)
}

dfs = rbindlist(sigtest_seed_df)
csv_write(dfs, 'sig_test_seed_change.csv')

ulst <- lapply(dfs, unique)
k <- lengths(ulst)
sum(k)

myvec = sapply(colnames(idf)[56:59],
       function(x, d) broom::tidy(glm(as.formula(paste("preterm ~ race9c + MATEDUC + INSPRIVATE + ns(MATAGE, df=3)", x, sep = " + ")), family = 'binomial', data = d))[18,4],
       d = idf)

test_remove_row = function(df, cohort_name, max_iter){
  for (tmpcol in 1:ncol(df)){
  
  }
}

idf = indic_df %>% dplyr::select(-all_of(adj_col_remove))
df = idf
max_iter = 10
myseeds = sample(2^16, ncol(df))
remove_vec = vector('list', length = ncol(df))  # empty vector to save output
for (i in 6:ncol(df)){
  print(i)
  set.seed(myseeds[[i - 5]])
  
  tmpname = colnames(df)[[i]]
  tmpsum = sum(df[, i])
  condition_rownums = which(df[, i] == 1)
  
  n_remove = min(tmpsum, max_iter)
  pvec = vector('list', length = max_iter)  # empty vector to save output
  if (tmpsum < max_iter & tmpsum > 0){
    for (j in 1:n_remove){
      newdf = df[-condition_rownums[[j]],]
      pvec[j] = as.numeric(broom::tidy(glm(as.formula(paste("preterm ~ race9c + MATEDUC + INSPRIVATE + ns(MATAGE, df=3)", tmpname, sep = " + ")), family = 'binomial', data = newdf))[18,5])
    }
  } else if (tmpsum >= max_iter) {
    remove_rownum = sample(condition_rownums, n_remove)
    for (j in 1:n_remove) {
      newdf = df[-condition_rownums[[j]],]
      pvec[j] = as.numeric(broom::tidy(glm(as.formula(paste("preterm ~ race9c + MATEDUC + INSPRIVATE + ns(MATAGE, df=3)", tmpname, sep = " + ")), family = 'binomial', data = newdf))[18,5])
    }
  }
  remove_vec[[i-5]] = pvec
}

# process output
rdf = rbindlist(remove_vec)  # to datatable
rdf = data.frame(rdf)  # to dataframe

write.csv(rdf,"Z://Birth_IRB1722929/projects/ptb_diagnoses/output/removal_overall_test.csv", row.names = F)
rdf = read.csv(paste0(data_dir, 'removal_indic_test.csv'))
rdf = rdf[,-1]

bh = vector('list', length = ncol(rdf))  # BH adjusted pvalues
for (i in 1:ncol(rdf)){
  pval_bh = p.adjust(rdf[,i], method = 'BH', n = nrow(rdf))
  bh[[i]] = pval_bh
}

pbh = data.frame(bh)
colnames(pbh) = letters[1:ncol(pbh)]

namesdf = df[6:ncol(df)]
namesdf = namesdf[,colSums(namesdf) > 0]

row.names(pbh) = colnames(namesdf) # phecodes as column names

total_sig = rowSums(pbh < .05, na.rm = T)
notna = rowSums(!is.na(pbh))



pbh = cbind(pbh, total_sig)
pbh = cbind(pbh, notna)
pbh = pbh[pbh$notna > 0, ]
pbh$percent_sig = pbh$total_sig/pbh$notna * 100


pbh$phestr = row.names(pbh)
pbh$phecode = str_sub(pbh$phestr, 2)
pbh = merge(pbh, phe_defs, by = 'phecode')

csv_write(pbh, 'bootstrap_indicated')

# plot results
#pbh = pbh[pbh$total_sig > 0 & pbh$category_number !=12,]
pbh = pbh[order(pbh$total_sig),]
rownames(pbh) = NULL  # renumber rows
myorder = as.vector(pbh$phenotype)
pbh$phenotype = as.factor(pbh$phenotype)
pbh$phenotype = fct_relevel(pbh$phenotype, myorder)  # relevel factor by row order


prtplt = ggplot(aes(x = total_sig, y= phenotype), data = pbh) + 
  geom_point() +
  ggtitle('Bootstrap, indicated, total') + 
  xlab('total # of significant iterations (max=10)') + 
  theme_bw() 

prtplt

ggsave(prtplt, width = 12, height = 6, 
       filename = image_name('bootstrap_indicated_total_sig'))


#pbh_pct = pbh[pbh$percent_sig > 0 & pbh$category_number !=12,]
pbh_pct = pbh_pct[order(pbh_pct$percent_sig),]
rownames(pbh_pct) = NULL  # renumber rows
myorder = as.vector(pbh_pct$phenotype)
pbh_pct$phenotype = as.factor(pbh_pct$phenotype)
pbh_pct$phenotype = fct_relevel(pbh_pct$phenotype, myorder)  # relevel factor by row order


prtplt = ggplot(aes(x = percent_sig, y= phenotype), data = pbh_pct) + 
  geom_point() +
  ggtitle('Bootstrap, indicated, percent') + 
  xlab('total % of significant iterations') + 
  theme_bw() 

prtplt

ggsave(prtplt, width = 12, height = 6, 
       filename = image_name('bootstrap_indicated_pct_sig'))


# random perturbation testing for significance
test_perturbations = function(df, remove_cols, cohort_name, iterations, pct, minrows){
  name_string = paste('preturbations', cohort_name, 
                      as.character(pct), as.character(minrows), 
                      as.character(iterations), sep = '_')
  
  df = df %>% dplyr::select(-all_of(remove_cols))
  myseeds = round(runif(iterations, 100, 2^16), 0)  # random seed for each iteration
  
  mysums = sapply(seq(6, ncol(df), 1),  # how many rows to change for each diagnosis
                  function(x) max(ceiling(sum(df[,x])*pct/100), minrows))

  
  pert = vector('list', length = length(myseeds))  # empty vector to save output
  
  for (k in 1:length(myseeds)){  # create a new perturbed dataframe
    set.seed(myseeds[[k]])

    for (i in seq(6, ncol(df), 1)){  # alter each column of the dataframe
      nchange = mysums[[i-5]]
      rchange = sample.int(nrow(ndf), nchange)
      
      for (j in rchange){
        ndf[j, i] = ndf[j, i] + 1 %% 2
      }
      
    }
    
    sigtest_perturb = sapply(colnames(idf)[6:ncol(df)],  # create model for each column of the dataframe
                             function(x, d) 
                               broom::tidy(glm(as.formula(paste("preterm ~ race9c + MATEDUC + INSPRIVATE + ns(MATAGE, df=3)", x, sep = " + ")), family = 'binomial', data = d))[18,5],
                             d = ndf)
  
    pert[[k]] = data.frame(sigtest_perturb)
  }
  set.seed(123)
  
  # process the results
  
  pertdf = rbindlist(pert)
  
  colnames(pertdf) = sapply(colnames(pertdf), function(x) gsub('.p.value', '', x))
  pertdf = pertdf %>% dplyr::select(-c(early, late))
  
  csv_write(pertdf, name_string)
  
  pertbh = vector('list', length = nrow(pertdf))  # BH adjusted pvalues
  for (i in 1:nrow(pertdf)){
    pval_bh = p.adjust(pertdf[i,], method = 'BH', n = nrow(pdf))
    pertbh[[i]] = pval_bh
  }
  
  pertbh = data.frame(pertbh)
  colnames(pertbh) = letters[1:ncol(pertbh)]
  
  pertbh = data.frame(t(pertbh))
  
  pertrow = vector('list', length = ncol(pertbh))
  for (i in 1:ncol(pertbh)){
    pertrow[[i]] = sum(pertbh[,i] < .05)
  }
  
  pertbh = rbind(pertbh, pertrow)
  pertbh = data.frame(t(pertbh))
  pertbh = pertbh %>% rownames_to_column(var = 'phestring')
  pertbh$sig_pct = pertbh$X1/25 * 100
  pertbh$phecode = gsub('P', '', pertbh$phestring)
  
  pertpltdf = merge(pertbh, phe_defs, by = 'phecode')  # add in phecode descriptions
  
  #pertbh = pertbh %>% filter(category_number != 12)
  
  pertpltdf = pertpltdf %>% filter(X1 > 1)
  
  # create factor to order plot
  pertpltdf = pertpltdf[order(pertpltdf$sig_pct),]
  rownames(pertpltdf) = NULL  # renumber rows
  myorder = pertpltdf$phenotype  
  pertpltdf$phenotype = as.factor(pertpltdf$phenotype)
  pertpltdf$phenotype = fct_relevel(pertpltdf$phenotype, myorder)  # relevel factor by row order
  
  
  prtplt = ggplot(aes(x = sig_pct, y= phenotype), data = pertpltdf) + 
    geom_point() +
    geom_vline(xintercept = 95) + 
    ggtitle('Overall perturbation test') + 
    theme_bw() 
  
  prtplt
  
  ggsave(prtplt, width = 12, height = 6, 
         filename = image_name(name_string))
  
}


# preparing data for drops
odf = cohort %>% dplyr::select(-all_of(adj_col_remove))
df = data.frame(odf)
threshold = 50
tdf = df[, 6:ncol(df)]
tdf = tdf[, colSums(tdf) > 1]

df = cbind(df[,1:5], tdf)

drop_list = vector('list', length = ncol(df) - 5)
for (i in 6:ncol(df)) {
  positives = which(df[, i] == 1)
  if(sum(df[, i]) <= threshold) {
    tmp = rep(positives, length.out = threshold)
  } else{
    tmp = sample(positives, size = threshold, replace = FALSE)
  }
  drop_list[[i-5]] = tmp
}

dropdf = data.frame(drop_list)
colnames(dropdf) = sapply(seq(6, ncol(df)), function(x, y) paste0(y, x), y = 'col')

ptm = proc.time()
pvals = vector('list', length = threshold)
for (i in 1:threshold) {
  df_na = data.frame(df)
  for (j in 6:ncol(df_na)){
    df_na[dropdf[i, j-5], j] = NA
  }
  pvals[[i]] = adjusted_model(df_na)
}
b = proc.time() - ptm
pvaldf = do.call(rbind, pvals)
colnames(pvaldf) = sapply(colnames(pvaldf), function(x) gsub('.p.value', '', x))
pvaldf = pvaldf %>% dplyr::select(-c('early', 'late'))
csv_write(pvaldf, 'overall_droptest_50')

pvaldf = read.csv('Z://Birth_IRB1722929/projects/ptb_diagnoses/output/2023-08-26_indicated_droptest_50.csv')

#pvdf = data.frame(t(pvaldf))
pvdf = pvaldf
cnames = colnames(pvdf)
cnames = lapply(cnames, function(x) gsub('.p.value', '', x))
colnames(pvdf) = cnames

pertbh = vector('list', length = nrow(pvdf))  # BH adjusted pvalues
for (i in 1:nrow(pvdf)){
  pval_bh = p.adjust(pvdf[i,], method = 'BH', n = ncol(pvdf))
  pertbh[[i]] = pval_bh
}

pertbh = data.frame(pertbh)
colnames(pertbh) = letters[1:ncol(pertbh)]

pertbh = data.frame(t(pertbh))

pertrow = vector('list', length = ncol(pertbh))
for (i in 1:ncol(pertbh)){
  pertrow[[i]] = sum(pertbh[,i] < .05)
}

pertbh = rbind(pertbh, pertrow)
pertbh = data.frame(t(pertbh))
pertbh = pertbh %>% rownames_to_column(var = 'phestring')
pertbh$sig_pct = pertbh$X1/threshold * 100
pertbh$phecode = gsub('P', '', pertbh$phestring)

pertpltdf = merge(pertbh, phe_defs, by = 'phecode')  # add in phecode descriptions

#pertpltdf = pertpltdf %>% filter(category_number != 12)

pertpltdf = pertpltdf %>% filter(X1 > 1)

# create factor to order plot
pertpltdf = pertpltdf[order(pertpltdf$sig_pct),]
rownames(pertpltdf) = NULL  # renumber rows
myorder = pertpltdf$phenotype  
pertpltdf$phenotype = as.factor(pertpltdf$phenotype)
pertpltdf$phenotype = fct_relevel(pertpltdf$phenotype, myorder)  # relevel factor by row order
csv_write(pertpltdf, 'indicated_drop_50_summary')

prtplt = ggplot(aes(x = sig_pct, y= phenotype), data = pertpltdf) + 
  geom_point() +
  #geom_vline(xintercept = 95) + 
  ggtitle('Indicated drop test') + 
  theme_bw() 

prtplt
ggsave(prtplt, width = 12, height = 6, 
                        filename = image_name('indicated_drop_test'))
