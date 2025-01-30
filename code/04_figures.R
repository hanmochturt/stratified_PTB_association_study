diag_req <<- 1  # flip this binary variable depending on whether you require a diagnosis prior to pregnancy


base_dir = '//ars-data-01.sde.net.ucsf.edu/MyResearchShared/sirotam1_shared/Birth_IRB1722929/projects/ptb_diagnoses/'
defs_dir <<- paste0(base_dir, 'defs/')

if (diag_req == 1) {
  data_dir <<- paste0(base_dir, 'output/01_revisions/')
  image_dir <<- paste0(base_dir, 'images/01_revisions/')
} else {
  data_dir <<- paste0(base_dir, 'output/01_revisions/nodiagreq/')
  image_dir <<- paste0(base_dir, 'images/01_revisions/nodiagreq/')
}


library(tidyverse)
library(forcats)
library(lubridate)
library(ggrepel)
library(ggplot2)
library(cowplot)

file_select = function(myfile) {
  # this function picks the most recent file by datestring in the filename
  files = list.files(data_dir, pattern = paste0("[0-9]{4}-[0-9]{2}-[0-9]{2}_", myfile))
  maxdate = max(do.call(c, lapply(files, function(x) ymd(str_sub(x, 1, 10)))))
  return(paste0(maxdate, '_', myfile, '.csv'))
}

image_name = function(img_name) {
  if(diag_req == 1){
    return(paste0(image_dir, Sys.Date(), '_', img_name, '.png'))
  } else {
   return(paste0(image_dir, Sys.Date(), '_nodiagreq_', img_name, '.png'))
  }
}

csv_write = function(df, dfname) {
  if (diag_req == 1) {
    write.csv(df, file = paste0(data_dir, Sys.Date(), '_', dfname, '.csv'), row.names = F)
  } else {
    write.csv(df, file = paste0(data_dir, Sys.Date(), '_nodiagreq_', dfname, '.csv'), row.names = F)
  }
}


category_rename = function(df) {
  # any renaming or reordering of categories should be done here
  df = df %>% mutate(category = case_when(category == "NULL" ~ "uncategorized",
                                          TRUE ~ category))
  return(df)
}

publish_file = function(bootdf, oddsdf, outname){
  bootdf = bootdf %>% dplyr::select(phecode, sig_pct)
  
  odds_print = left_join(bootdf, oddsdf, by = 'phecode')
  odds_print$n = odds_print$count_preterm + odds_print$count_term
  odds_print$robust_pct = odds_print$sig_pct
  
  odds_print = odds_print %>% filter(n >= 10)
  
  odds_publish = odds_print %>% 
    select(phecode, phenotype, n, or, or_lower, or_upper, pval_bh, robust_pct)
  odds_publish <- odds_publish[order(odds_publish$pval_bh),]

  odds_publish = odds_publish %>% mutate(pval_bh = signif(pval_bh, digits = 3),
                                         or = round(or, 2),
                                         or_lower = round(or_lower, 2),
                                         or_upper = round(or_upper, 2))
  
  csv_write(odds_publish, paste0('publish_', outname))
  
  return(odds_print)
  
}

spont_file = function(df, outname){
  df$n = df$count_preterm + df$count_term
  df = df %>% filter(n > 10)
  
  df_publish = df %>%
    select(phecode, phenotype, n, or, or_lower, or_upper, pval_bh)
  df_publish <- df_publish[order(df_publish$pval_bh),]
  df_publish = df_publish %>% filter(pval_bh <= .5)
  
  df_publish = df_publish %>% mutate(pval_bh = signif(pval_bh, digits = 3),
                                         or = round(or, 2),
                                         or_lower = round(or_lower, 2),
                                         or_upper = round(or_upper, 2))
  
  csv_write(df_publish, paste0('publish_', outname))
}


manhattan_plot = function(df, imgname, label_pval){
  pos <- position_jitter(width = 0.3, seed = 2)
  
  df = df %>% filter(category_number != 12)
  df = category_rename(df)
  
  catvars = df %>% dplyr::select(category_number, category) %>% distinct()
  catvars = catvars[order(catvars$category_number), ]
  catorder = catvars$category
  
  catcols = color_list[catorder]
  
  myplot = ggplot() +
    geom_jitter(data = filter(df, pval_bh > label_pval),
                mapping = aes(x = category, y = padj_minus_log, color = category)) +
    geom_jitter(data = filter(df, pval_bh <=label_pval),
                mapping = aes(x = category, y = padj_minus_log, color = category),
                position = pos) +
    geom_text_repel(data = filter(df, pval_bh <= label_pval),
                    mapping = aes(x=category, y=padj_minus_log, label = phenotype, color = category),
                    max.overlaps = 100,
                    box.padding = .35,
                    position = pos,
                    size = 4) +
    geom_hline(yintercept = -log10(.05), linetype = 2) +
    scale_color_manual(values = catcols) + 
    ylab('-log10(pvalue)') +
    xlab('Phecode category') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none",
          text = element_text(size = 19),
          plot.margin = margin(0,0,20,0))
  
  fname = image_name(paste0(imgname, '_manhattan'))
  
  ggsave(myplot, filename = fname, width = 12, height = 7)
  return(myplot)
  
}

manhattan_plot_fct = function(df, imgname, label_pval){
  pos <- position_jitter(width = 0.3, seed = 2)
  
  df = df %>% filter(category_number != 12)
  df = category_rename(df)
  
  catvars = df %>% dplyr::select(category_number, category) %>% distinct()
  catvars = catvars[order(catvars$category_number), ]
  catorder = catvars$category
  
  catcols = color_list[catorder]
  
  myplot = ggplot() +
    geom_jitter(data = filter(df, pval_bh > label_pval),
                mapping = aes(x = category, y = padj_minus_log, color = category)) +
    geom_jitter(data = filter(df, pval_bh <=label_pval),
                mapping = aes(x = category, y = padj_minus_log, color = category),
                position = pos) +
    geom_text_repel(data = filter(df, pval_bh <= label_pval),
                    mapping = aes(x=category, y=padj_minus_log, label = phenotype, color = category),
                    max.overlaps = 100,
                    box.padding = .35,
                    position = pos,
                    size = 4) +
    geom_hline(yintercept = -log10(.05), linetype = 2) +
    scale_color_manual(values = catcols) + 
    ylab('-log10(pvalue)') +
    xlab('Phecode category') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none",
          text = element_text(size = 19),
          plot.margin = margin(0,0,20,0)) + 
    facet_grid(subtype ~ .)
  
  fname = image_name(paste0(imgname, '_manhattan'))
  
  ggsave(myplot, filename = fname, width = 10, height = 12)
  
}

volcano_plot = function(df, imgname, label_pval){
  
  df = df %>% filter(category_number != 12)
  
  volc_plt = ggplot() +
    geom_point(aes(x = or, y = padj_minus_log, color = category), data = df[df$pval_bh < label_pval,]) +
    geom_text_repel(data = filter(df, pval_bh < label_pval),
                    mapping = aes(x=or, y=padj_minus_log, label = phenotype, color = category),
                    max.overlaps = 100,
                    box.padding = .35,
                    size = 5) +
    geom_point(aes(x = or, y = padj_minus_log), data = df[df$pval_bh >= label_pval, ]) +
    scale_x_log10() +
    geom_vline(xintercept = 1) +
    geom_hline(yintercept = -log10(.05), linetype = 2) +
    xlab('Odds Ratio') +
    ylab('-log(pvalue)') +
    theme_bw() +
    theme(text = element_text(size = 15), legend.position = "none")
  volc_plt
  
  fname = image_name(paste0(imgname, '_volcano'))
  
  ggsave(volc_plt, filename = fname)
  
}


forest_plot = function(df, imgname, ttl, legendx, legendy, ht){
  
  df = df %>% mutate(phenotype = case_when(phecode == '614.1' ~ 'Pelvic peritoneal adhesions',
                                           phecode == '580.3' ~ 'Nephritis and nephropathy',
                                           phecode == '793.2' ~ 'Intrathoracic organ abnormalities',
                                           phecode == '587' ~ 'Kidney replaced by transplant',
                                           phecode == '598.9' ~ 'Nonspecific urine findings',
                                           TRUE ~ phenotype))
  
  # group by category
  df = df[order(df$category_number), ]
  catvars = df %>% dplyr::select(category_number, category) %>% distinct()
  catvars = catvars[order(catvars$category_number), ]
  catorder = catvars$category
  
  df$category = as.factor(df$category)
  df$category = fct_relevel(df$category, catorder)
  
  ordervars = df$phenotype
  df$phenotype = as.factor(df$phenotype)
  df$phenotype = fct_relevel(df$phenotype, ordervars)
  
  df = df[df$pval_bh < .05, ]
  catcols = color_list[catorder]
  
  forestplt = ggplot(aes(x = or, xmin = or_lower, xmax = or_upper, y = phenotype, color = category), data = df) +
    geom_pointrange() +
    guides(colour = guide_legend(reverse=T)) + 
    geom_vline(xintercept = 1, lty = 2) +
    theme_bw() +
    scale_color_manual(values = catcols) + 
    scale_x_log10() +
    xlab('Odds Ratio') + 
    ggtitle(ttl) + 
    theme(legend.position = c(legendx, legendy),
      text = element_text(size = 18),
      axis.title.y = element_blank(),
      legend.text=element_text(size=8),
      legend.title = element_text(size=8),
      legend.margin = margin(0, 0, 0, 0), 
      plot.margin = margin(20, 0, 0,0))
  forestplt
  
  fname = image_name(paste0(imgname, '_forest'))
  
  ggsave(forestplt, filename = fname, width = 8, height = ht)
  return(forestplt)
}

#===============#
# READ FILES ####
#===============#
phe_defs = read.csv(paste0(defs_dir, 'phecode_definitions1.2.csv'), colClasses = c('phecode'='character'))
phe_defs = category_rename(phe_defs)
categories = lapply(list(unique(phe_defs$category)), sort)  # alphabetic order of category names

# read adjusted
adj_odds = read.csv(paste0(data_dir, file_select('overall_adj_odds_ratios')))
bootstrap_all = read.csv(paste0(data_dir, file_select('overall_drop_50_summary')))
adj_odds_ind = read.csv(paste0(data_dir, file_select('indicated_adj_odds_ratios')))
bootstrap_ind = read.csv(paste0(data_dir, file_select('indicated_drop_50_summary')))
adj_odds_spont = read.csv(paste0(data_dir, file_select('spontaneous_adj_odds_ratios')))


# read crude
crude_odds = read.csv(paste0(data_dir, file_select('overall_crude_odds_ratios')))
bootstrap_all_crude = read.csv(paste0(data_dir, file_select('overall_crude_drop_50_summary')))
crude_odds_ind = read.csv(paste0(data_dir, file_select('indicated_crude_odds_ratios')))
bootstrap_ind_crude = read.csv(paste0(data_dir, file_select('indicated_crude_drop_50_summary')))
crude_odds_spont = read.csv(paste0(data_dir, file_select('spontaneous_crude_odds_ratios')))

#=======================================#
# DEIDENTIFY OUTPUT FILES FOR UPLOAD ####
#=======================================#


# print adjusted odds files
adj_all = publish_file(bootstrap_all, adj_odds, 'overall_adjusted')
adj_ind = publish_file(bootstrap_ind, adj_odds_ind, 'indicated_adjusted')

# print crude odds files
crude_all = publish_file(bootstrap_all_crude, crude_odds, 'overall_crude')
crude_ind = publish_file(bootstrap_ind_crude, crude_odds_ind, 'indicated_crude')

# print spontaneous files (no significant results)
spont_file( adj_odds_spont, 'spontaneous_adjusted')
spont_file(crude_odds_spont, 'spontaneous_crude')

#==================#
# DEFINE COLORS ####
#==================#
# http://mkweb.bcgsc.ca/colorblind/palettes.mhtml#24-color-palette-for-colorbliness

color_list = list('circulatory system' = '#003D30',
              'congenital anomalies' = '#FF8735',
              'dermatologic' = '#810D49',
              'digestive' = '#DE0D2E',
              'endocrine/metabolic' = '#0079FA',
              'genitourinary' = '#6B069F',
              'hematopoietic' = '#FFB935',
              'infectious diseases' = '#009175',
              'injuries & poisonings' = '#D80D7B',
              'mental disorders' = '#00489E',
              'musculoskeletal' = '#B20725',
              'neoplasms' = '#00B408',
              'neurological' = '#FF4235',
              'respiratory' = '#00E5F8',
              'sense organs' = '#B40AFC',
              'symptoms' = '#005A01',
              'uncategorized' = '#808080')



#==================#
# ADJUSTED plots ####
#==================#
## overall plots ####
manhattan_plot(adj_odds, 'adj_overall', .002)

adj_all_sig = adj_all %>% filter(robust_pct == 100)


if (diag_req == 1) {
  forest_plot(adj_all_sig, 'adj_overall', 'Overall', .8, .22, 10)
} else {
  forest_plot(adj_all_sig, 'adj_overall', 'Overall', .8, .55, 8)
}

## indicated plots ####
adj_ind_sig = adj_ind %>% filter(robust_pct == 100)

manhattan_plot(adj_odds_ind, 'ind_adj', .0001)

if (diag_req == 1) {
  forest_plot(adj_ind_sig, 'ind_adj', 'Indicated subtype', .84, .89, 12) 
} else {
  forest_plot(adj_ind_sig, 'ind_adj', 'Indicated subtype', .84, .15, 10)
}




#====================#
# CRUDE plots ####
#====================#

bootstrap_crude = bootstrap_all_crude %>% filter(sig_pct == 100) %>%
  dplyr::select(phestring)
crude_odds_all_robust = merge(crude_odds, bootstrap_crude, by = 'phestring')

manhattan_plot(crude_odds_ind, 'crude_indicated', .000001)
manhattan_plot(crude_odds, 'crude_overall', .00005)

#======================#
# SPONTANEOUS plots ####
#======================#

manhattan_plot(adj_odds_spont, 'spont_adj', .1)


#================================#
# FACET PLOT: INDIC AND SPONT ####
#================================#
adj_odds_ind$subtype = 'Indicated'
adj_odds_spont$subtype = 'Spontaneous'

adj_odds_subtype = rbind(adj_odds_ind, adj_odds_spont)

manhattan_plot_fct(adj_odds_subtype, 'facet', .0001)

#==================#
# COMBINED PLOT ####
#==================#
panela = manhattan_plot(adj_odds, 'adj_overall', .002)
panelb =  forest_plot(adj_odds_robust, 'adj_overall')
pg = plot_grid(panela, panelb, nrow = 2, labels = "AUTO")
pg
save_plot(pg, filename = image_name('overall_panels'),  base_height = 16, base_width = 8 , dpi = 600)










