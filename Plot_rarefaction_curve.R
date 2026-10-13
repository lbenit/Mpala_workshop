########################################################################
# Make a rarefaction curve
# November 6, 2024
# Lorena Benitez

library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(tidyr)
library(iNEXT)
###############################
# Make species matrix
plots<-read.csv('plot_data_w_chave14.csv')
s<-read.csv('stem_data_w_model_height.csv')
species<-as.data.frame(table(s$plot_id,s$species_name_clean))


species<-species%>%
  rename(plot_id=Var1)%>%
  left_join(plots,by='plot_id')%>%#select(-contains('p_'))%>%
  relocate(veg_type, .after=plot_id)%>%
  mutate(Var2=gsub(' ','.',Var2))

species<-species[,1:4]

# Collect by veg type
species1<-species%>%
  group_by(Var2,veg_type)%>%
  summarise(n=sum(Freq))


species2<-species1%>%
  pivot_wider(names_from = veg_type,values_from = n)%>%
  rename(species=Var2)

input<-as.data.frame(species2)
row.names(input)<-input$species

#Remove species column
input<-input[,-1]

#####
D_abund <- iNEXT (input,q=0,datatype = 'abundance', nboot=100)
#, size=c(18,121,37,193,84,50)
plot(D_abund)
mtext('Species Richness', side=2, line=3)
D_abund

ggiNEXT(D_abund, type=2)

#####################################
# PLot 1
plot<-ggiNEXT(D_abund, type=1)

plot1<-plot+scale_fill_viridis_d(option='turbo', na.value = "transparent")+
  scale_color_viridis_d(option='turbo', na.value = "transparent")+
  theme_minimal()+xlim(0,6000)+
  annotate(geom="text", x=20, y=150, label="A",size=14,
           color="black")+
  ylab('Species Richness')+
  theme(
    legend.text = element_text(size=20),
    legend.title = element_text(size=30),
    legend.key.size = unit(1, 'cm'),
    axis.title.y=element_text(size=20),
    axis.title.x=element_text(size=20),
    axis.text = element_text(size=15)
  )
plot1



##############################################################
#Ignore
##############################################################################
# Generate based on incidence frequency
plot_n<-as.data.frame(table(plots$veg_type))


###########
# Species

inc_freq<-species%>%mutate(Freq=case_when(Freq>=1~1,
                                          TRUE~0))

sum_frq<-inc_freq%>%group_by(veg_type,Var2)%>%
  summarise(freq=sum(Freq))%>%
  pivot_wider(names_from = Var2,values_from = freq)


sum_frq<-sum_frq[,-1]

# COmbine
input3<-cbind(plot_n,sum_frq)
rownames(input3)<-input3[,1]
input3<-input3[,-1]

library(purrr)
row_list <- input3 |> split(seq_len(nrow(input3))) |> map(as.numeric)
names(row_list)<-row.names(input3)



#####
D_freq <- iNEXT (row_list,q=0,datatype = 'incidence_freq', nboot=100)
plot(D_freq)

D_freq


#####################################
# PLot 1
plot2<-ggiNEXT(D_freq, type=1)

plot2<-plot2+scale_fill_viridis_d(option='turbo', na.value = "transparent")+
  scale_color_viridis_d(option='turbo', na.value = "transparent")+
  annotate(geom="text", x=20, y=150, label="B",size=14,
           color="black")+
  theme_minimal()+
  ylab('Species Richness')+
  theme(
    legend.text = element_text(size=20),
    legend.title = element_text(size=30),
    legend.key.size = unit(1, 'cm'),
    axis.title.y=element_text(size=20),
    axis.title.x=element_text(size=20),
    axis.text = element_text(size=15)
  )

plot2



plots<-(plot1 / plot2)

final_plot <- plots +
  plot_layout(guides = "collect")

final_plot

ggsave('figure_3.png', final_plot, bg='white',
       width = 8,
       height = 11,
       units = 'in',
       dpi = 300)
