
library(DBI)
library(odbc)
library(glue)
library(lubridate)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(forcats)

fp_base <<- '//ars-data-01.sde.net.ucsf.edu/MyResearchShared/sirotam1_shared/Birth_IRB1722929/projects/ptb_diagnoses/'

fp_data <<- paste0(fp_base, 'data/')
fp_defs <<- paste0(fp_base, 'defs/')
fp_out <<- paste0(fp_base, 'output/')
fp_img <<- paste0(fp_base, 'images/')

mydata = read.csv(paste0(fp_data, '2023-07-18_conditions_binary.csv'))
phedefs = read.csv(paste0(fp_defs, 'phecode_definitions1.2.csv'), colClasses = c('phecode'='character'))

wbc = phedefs %>% filter(phenotype == 'Decreased white blood cell count', ) %>%
  select(phecode) %>% distinct()
wbc = wbc$phecode[1]
mycol = paste0('P', wbc)

cohort = mydata %>% filter(P288.1 == 1)  # all the people with a diagnosis of decreased white blood cell count

# look at what other diagnoses they have
cohort <- cohort %>% select(where(~ any(. != 0))) %>% select(-P288.1)
phesums = colSums(cohort[,-c(1,2)])
freqdf = data.frame(phecol = names(phesums), counts = phesums)

freqdf <- freqdf %>% mutate(phecode = str_sub(phecol, 2))
freqdf = merge(freqdf, phedefs, by= 'phecode')
freqdf = freqdf[order(freqdf$counts),]

pltdf = freqdf[freqdf$counts >= 10,]  # filter to 10 or more people
myorder = pltdf$phenotype
pltdf$phenotype = factor(pltdf$phenotype)
pltdf$phenotype = fct_relevel(pltdf$phenotype, myorder)

barplt = ggplot(aes(x = counts, y = phenotype), data = pltdf) + 
  geom_col() + 
  xlab('Number of people with phecode') + 
  ggtitle('Phecodes present in those with decreased white blood cell count') +
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        plot.title.position = "plot")
barplt
ggsave(barplt, file = paste0(fp_img, Sys.Date(), '_decreased_white_blood_cell_phecodes.png'))




