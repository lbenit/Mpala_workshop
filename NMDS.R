#Load packages
library(dplyr)
library(ggplot2)
library(sf)
library(stringr)
library(vegan)
library(mobr)
library(tibble)
library(tidyr)

#############################
# Read in data
############################

plots = read.csv('C:/Users/s2437177/Documents/Github/ZFC_datset_cleaning/all_ZWE_plots_data.csv')
stems = read.csv('C:/Users/s2437177/Documents/Github/ZFC_datset_cleaning/all_ZWE_stem_data.csv')

plot_data<-plots%>%filter(census_date>=2012)%>%
  mutate(census_id=paste0(plot_id,'_',census_date))%>%
  filter(!census_id%in%c(
  'ZPF_10_2015','ZPF_9_2015','ZPF_8_2015','ZPF_7_2015','ZPF_6_2012',
  'ZPF_5_2012','ZPF_4_2012','ZPF_11_2015'))%>%
  filter(min_diam_thresh<=10)%>%
  filter(n_stems>0)%>%
  filter(agb_ge10>0)%>%
  distinct()

plots_basic<-plot_data%>%select(plot_id, longitude_of_centre,latitude_of_centre)


#Remove stems <10 dbh
stems<-stems%>%filter(diam>=10)%>%mutate(census_id=paste0(plot_id,'_',census_date))%>%
  filter(census_id%in%plot_data$census_id)

plots_basic<-plots_basic%>%filter(plot_id %in%stems$plot_id)

# Remove plots with non-native species
#remove<-c('Pinus', 'Eucalyptus', 'Elaeis', 'Lantana', 
#          'Jacaranda', 'Agave', 'Tectona', 'Cercidoideae')
#stems_filtered <- stems %>% group_by(plot_id) %>% filter(!any(genus_clean %in% remove)) %>% ungroup()


#Relative abundance
abundance<-stems%>% 
  group_by(plot_id,genus_clean)%>% 
  summarise(n=n())%>%pivot_wider(names_from = genus_clean,values_from = n)%>%
  mutate(
    across(everything(), ~replace_na(.x, 0))
  )%>%
  mutate(
    across(everything(), ~abs(.x))
  )%>%column_to_rownames("plot_id")

#Remove rare species
yn<-abundance
yn[yn > 0] <- 1 
keep<-as.data.frame(colSums(yn))%>%filter(colSums(yn)>1)# Change number to adjust threshold for inclusion
keep<-rownames(keep)
keep
keep<-keep[-8]#remove indet


#Prep for analysis
abund<-abundance[,keep]
#ab<-decostand(abund, method = "total")




######################################################

###########################
# Get basal area instead
# Get basal area instead
df_composition = stems %>%
  group_by(plot_id) %>%
  summarise(
    ba_total    = sum(ba, na.rm=T))

bio<-stems%>% 
  left_join(df_composition,by='plot_id')%>%
  group_by(plot_id,genus_clean)%>% 
  summarise(ba=sum(ba),
            ba_per=ba/ba_total)%>%
  select(-ba)%>%
  unique()


reframe<-bio%>%
  pivot_wider(names_from = genus_clean,values_from = ba_per)%>%
  filter_all(any_vars(!is.na(.)))%>% 
  mutate(
    across(everything(), ~replace_na(.x, 0))
  )%>%
  mutate(
    across(everything(), ~abs(.x))
  )%>%column_to_rownames("plot_id")%>%
  mutate(sum=rowSums(.))


summary(reframe$sum)
hist(reframe$sum)
#Filter species to remove singletons
bi<-reframe[,keep]

bi<-bi%>%
  select_if(function(x) !all(is.na(x)))%>%
  mutate(across(everything(), ~replace_na(., 0)))%>% 
  mutate(sum=rowSums(.))%>%
  filter(sum!=0)%>%
  select(-sum)#Very important because if included it will not converge

#Drop columns with no values
tail(colnames(bi))
tail(rownames(bi))


bio_hel <- decostand(bi, method = "hellinger")

##################################################################
library(vegan)

############################################

#Analysis
bio_distmat <- 
  vegdist(bio_hel, method = "bray")
bio_DistanceMatrix <- 
  as.matrix(bio_distmat, labels = T)

# Running NMDS in vegan (metaMDS)
bio_NMS <-
  metaMDS(bio_distmat,
          k = 3,
          maxit = 999, 
          trymax = 50)


sppscores(bio_NMS) <- bio_hel
plot(bio_NMS)
orditorp(bio_NMS, "species") 
orditorp(bio_NMS,'sites')

###############################################
# Remove outliers

scores_NMS <- scores(bio_NMS, display = "sites")

d <- dist(scores_NMS)
centroid <- colMeans(scores_NMS)
dist_to_centroid <- apply(scores_NMS, 1, function(x) sqrt(sum((x - centroid)^2)))

# Identify outliers (e.g., top 5% farthest)
cutoff <- quantile(dist_to_centroid, 0.95)
outliers <- names(dist_to_centroid[dist_to_centroid > cutoff])
outliers

bi_clean <- bi[!rownames(bi) %in% outliers, ]

bio_hel <- decostand(bi_clean, method = "hellinger")

#Analysis
bio_distmat <- 
  vegdist(bio_hel, method = "bray")
bio_DistanceMatrix <- 
  as.matrix(bio_distmat, labels = T)

# Running NMDS in vegan (metaMDS)
bio_NMS <-
  metaMDS(bio_distmat,
          k = 3,
          maxit = 999, 
          trymax = 200,
          try=20)


sppscores(bio_NMS) <- bio_hel
plot(bio_NMS)
orditorp(bio_NMS, "species") 
orditorp(bio_NMS,'sites')

####
# Stats
bio_NMS
bio_NMS$stress
stressplot(bio_NMS)
gof<-goodness(bio_NMS)
mean(gof)

# 4. Convergence info
bio_NMS$iters        # number of iterations
bio_NMS$converged    # TRUE/FALSE
bio_NMS$engine       # algorithm used


###############CLustering


#####################################
# K means
nms_scores <- as.data.frame(scores(bio_NMS, display = "sites"))

wss <- sapply(1:10, function(k)
  kmeans(nms_scores, centers = k, nstart = 20)$tot.withinss)

plot(1:10, wss, type = "b",
     xlab = "Number of clusters (k)",
     ylab = "Within-cluster sum of squares")


set.seed(123)   # for reproducibility
km <- kmeans(nms_scores, centers = 4)
nms_scores$cluster4 <- as.factor(km$cluster)



nms_species <- as.data.frame(scores(bio_NMS, display = "species"))

nms_species$species <- rownames(nms_species)

nms_sites <- as.data.frame(scores(bio_NMS, display = "sites"))
nms_sites$cluster <- as.factor(km$cluster)   # if using k-means
nms_sites$cluster4 <- as.factor(km$cluster)

ggplot(nms_scores, aes(NMDS1, NMDS2, color = cluster4)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "NMDS with k-means Clusters")


ggplot() +
  geom_point(data = nms_sites,
             aes(NMDS1, NMDS2, color = cluster4),
             size = 3) +
  # species arrows
  # geom_segment(data = nms_species,
  #              aes(x = 0, y = 0,
  #                  xend = NMDS1, yend = NMDS2),
  #              arrow = arrow(length = unit(0.25, "cm")),
  #              color = "black") +
  
  # species labels
  geom_text(data = nms_species,
            aes(NMDS1, NMDS2, label = species),
            color = "black", size = 3, vjust = -0.5) +
  
  theme_minimal() 

nms_scores$plot_id<-row.names(nms_scores)

write.csv(nms_scores,file='Zim_2026/nmds_scores.csv')
write.csv(nms_species,file='Zim_2026/nmds_species.csv')


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



