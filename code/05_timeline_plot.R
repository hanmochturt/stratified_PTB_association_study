library(vistime)
library(ggplot2)
#library(ggchicklet)

timeline_data <- data.frame(event = c('Conditions 1', "Pregnancy 1", 'Conditions 2', "Pregnancy 2"),
                            start = as.Date(c('2013-01-01', "2014-01-01", "2015-03-01", '2015-12-31')), 
                            end   = as.Date(c("2014-01-01", "2014-09-01", "2015-12-31", '2016-09-01')),
                            color = c('#a6cee3','#1f78b4','#b2df8a','#33a02c'),
                            fontcolor = c('black', 'white', 'black', 'white'),
                            group = "My Events")


mycolors = c('#a6cee3','#b2df8a', '#1f78b4','#33a02c')

myplot = ggplot(aes(y = group, xmin = start, xmax = end, color = event), data = timeline_data) + 
  geom_linerange(linewidth = 20) + 
  scale_color_manual(values = mycolors) + 
  geom_label(x = timeline_data$start + 90, label = timeline_data$event, size = 6) + 
  geom_linerange(aes(y = group, xmin = as.Date('2014-09-01'), xmax = as.Date('2015-03-01')), 
                 data = timeline_data,
                 linewidth = 5, 
                 color = '#fb9a99') + 
  geom_label(x = as.Date('2014-12-01'), label = 'Excluded\nconditions', color = '#fb9a99', size = 6) + 
  scale_x_date(date_breaks = "3 month", date_labels =  "%b %Y") + 
  theme_bw() + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 15),
        axis.title.y=element_blank(),legend.position="none",
        axis.ticks = element_blank(),
        panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())

myplot


ggsave(myplot, width = 15, height = 2, file = 'Z:/Birth_IRB1722929/projects/ptb_diagnoses/images/timeline.png')
