####
# Indicator species
library(indicspecies)

# indicator species analysis
ind <- multipatt(
  x = bio_hel,
  cluster = nms_sites$cluster4 ,
  func = "r.g"
)

summary(ind)
ind_table <- ind$sign

####################################
# Reorganize data for saving
plot_assignments <- nms_sites %>%
  mutate(
    `Floristic group` = recode(cluster4,
                               `1` = "Mopane",
                               `2` = "Baikiaea",
                               `3` = "Miombo",
                               `4` = "Burkea"))

plot_assignments$plot_id<-row.names(plot_assignments)
row.names(plot_assignments)<-NULL

############
# PLot with new colors
library(patchwork)
nmds12<-ggplot() +
  geom_point(data = plot_assignments,
             aes(NMDS1, NMDS2, color = `Floristic group`),
             size = 3) +
  theme_classic()+
  theme(#legend.position = 'none',
    #axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    legend.text = element_text(size=20),
    legend.title = element_text(size=30),
    legend.key.size = unit(1, 'cm'),
    axis.title=element_text(size=20),
    axis.text= element_text(size=18))
nmds23<-ggplot() +
  geom_point(data = plot_assignments,
             aes(NMDS1, NMDS3, color = `Floristic group`),
             size = 3) +
  theme_classic()+
  theme(#legend.position = 'none',
    #axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    legend.text = element_text(size=20),
    legend.title = element_text(size=30),
    legend.key.size = unit(1, 'cm'),
    axis.title=element_text(size=20),
    axis.text= element_text(size=18))

nmds_plot<-(nmds12/nmds23)+plot_layout(guides = 'collect')+ plot_annotation(tag_levels = 'A')
nmds_plot


###############################################################
# Save data

bi$plot_id<-rownames(bi)
bi_info<-bi%>%select('Brachystegia','Burkea','Baikiaea','Colophospermum','Guibourtia','Combretum','Erythrophleum','plot_id')
join<-left_join(plot_assignments,bi_info, by='plot_id')




write.csv(bi,file='Zim_2026/basal_area_genus_percentage.csv')
write.csv(plot_assignments,file='Zim_2026/NMDS_clustering.csv' )