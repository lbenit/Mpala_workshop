# Zimbabwe Plot Data Figures
# Five panel figure with diversity metrics
#### Also creates summary tables
# By: Lorena Benitez
# April 2026
# Plot data from process_clusters

#Load packages
library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(seosawr)
library(tidyr)
library(tibble)
library(vegan)
library(viridis)

#############################
# Read in data
############################
plots<-read.csv('single_census_plot.csv') 

stems<-read.csv('single_census_stems.csv')


#############################
# Clean data
############################

p<-plots%>%
  filter(plot_radius==20)%>% # Only keep plots that are the same size
  select(plot_id,n_stems_ge10, plot_area, plot_cluster)

s<-stems%>%
  filter(plot_id%in% p$plot_id)%>% #Filter to remove stems that are in plots that were removed because of plot size differences
  filter(diam_adj>=10)%>% #Only keep stems greater than or equal to 10 as some plots include smaller stems
  select(plot_id,stem_id,species_name_clean,diam_adj)


######################################
#Create abundance matrix
#######################################
abundance<-s%>% 
  group_by(plot_id,species_name_clean)%>% 
  summarise(n=n())%>%pivot_wider(names_from = species_name_clean,values_from = n)%>%
  mutate(
    across(everything(), ~replace_na(.x, 0))
  )%>%
  mutate(
    across(everything(), ~abs(.x))
  )%>%column_to_rownames("plot_id")

comm_matrix<-as.matrix(abundance)


#######################################
# #Calculate diversity metrics from stems
#######################################


# Wrapper function for community diversity metrics
calc_diversity <- function(comm_matrix) {
  # Species richness (number of species per site)
  richness <- specnumber(comm_matrix)
  
  # Shannon diversity
  shannon <- diversity(comm_matrix, index = "shannon")
  
  # Simpson diversity (1 - D)
  simpson <- diversity(comm_matrix, index = "simpson")
  
  # Inverse Simpson diversity (1 / D)
  inv_simpson <- diversity(comm_matrix, index = "invsimpson")
  
  # Pielou's evenness (J' = H' / ln(Richness))
  evenness <- shannon / log(richness)
  
  
  # Combine into a clean data frame
  results <- data.frame(
    plot_id = rownames(comm_matrix),
    richness = richness,
    shannon = shannon,
    simpson = simpson,
    invsimpson = inv_simpson,
    evenness = evenness
  )
  
  return(results)
}

diversity_data<-calc_diversity(comm_matrix)


#################################################
# Calculate rarefied richness
###############################################

# Rarified species richness
filtered_matrix <- comm_matrix[rowSums(comm_matrix) >= 10, ]# Remove plots with fewer than 10 stems
rarefied<-rarefy(filtered_matrix,sample=10)

rare_data<-data.frame(
  plot_id = rownames(filtered_matrix),
  rarefied = rarefied
)


#Join diversity data with rarefied richness
biodiversity_data<-left_join(diversity_data,rare_data,by='plot_id')%>%
                  left_join(p,by='plot_id')#combine with plot metadata


##############################
#Summarise by site
#############################
diversity_summary<-biodiversity_data%>%
  group_by(plot_cluster)%>%
   summarise(plot_n=n(),
            mean_richness=mean(richness),
            sd_richness=sd(richness),
            se_richness=sd_richness/sqrt(plot_n),
            mean_shannon=mean(shannon,na.rm=TRUE),
            sd_shannon=sd(shannon,na.rm=TRUE),
            se_shannon=sd_shannon/sqrt(plot_n),
            mean_simpson=mean(simpson,na.rm=TRUE),
            sd_simpson=sd(simpson),
            se_simpson=sd_simpson/sqrt(plot_n),
            mean_invsimpson=mean(invsimpson,na.rm=TRUE),
            sd_invsimpson=sd(invsimpson),
            se_simpson=sd_simpson/sqrt(plot_n),
            mean_even=mean(evenness,na.rm=TRUE),
            sd_even=sd(evenness,na.rm=TRUE),
            se_even=sd_even/sqrt(plot_n)
  )#%>% mutate(across(2:13, round, 3))

# Summarise rarefied richness on a si te level
rarefied_summary<-biodiversity_data%>%
  filter(n_stems_ge10>=10)%>%
  group_by(plot_cluster)%>%
  summarise(plot_n=n(),
            mean_rarified=mean(rarefied,na.rm=TRUE),
            sd_rarified=sd(rarefied,na.rm=TRUE),
            se_rarified=sd_rarified/sqrt(plot_n),
            mean_even=mean(evenness,na.rm=TRUE),
            sd_even=sd(evenness,na.rm=TRUE),
            se_even=sd_even/sqrt(plot_n))

# Join together and save
final_diversity_metrics<-left_join(diversity_summary,rarefied_summary, by='plot_cluster')

write.csv(final_diversity_metrics,'diversity_metrics_by_site.csv')



############################################
# Make plots for each metric
###########################################

# Make boxplot plot for rarefied richness
plot_rare <- ggplot(data=biodiversity_data, aes(y=rarefied, x=plot_cluster,fill=plot_cluster)) + 
  geom_boxplot()+scale_fill_viridis_d(option='turbo', na.value = "transparent")+
  ylab("Rarified Richness")+ 
  xlab("Plot cluster")+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.x=element_text(size=16),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y=element_text(size=20),
        axis.text.y=element_text(size=16))

# View 
plot_rare


# Plot evenness
plot_even <- ggplot(data=biodiversity_data, aes(y=evenness, x=plot_cluster,fill=plot_cluster)) + 
  geom_boxplot()+scale_fill_viridis_d(option='turbo', na.value = "transparent")+
  ylab("Evenness")+ 
  xlab("Floristic group")+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.x=element_text(size=16),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y=element_text(size=20),
        axis.text.y=element_text(size=16))

plot_even

# Plot species richness
plot_richness <- ggplot(data=biodiversity_data, aes(y=richness, x=plot_cluster,fill=plot_cluster)) + 
  geom_boxplot()+scale_fill_viridis_d(option='turbo', na.value = "transparent")+
  ylab("Species Richness")+ 
  xlab("Floristic group")+
  theme_classic()+
  theme(
        axis.text.x=element_text(size=16),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y=element_text(size=20),
        axis.text.y=element_text(size=16))

plot_richness


# Use patchwork to make it into one figures
final_plot<-(plot_rare+plot_even)/(plot_richness+ guide_area()+plot_layout(guides = 'collect'))+plot_annotation(tag_levels = 'A')


final_plot 


# Save final plot
ggsave('diversity_metrics.png', final_plot, bg='transparent',
       width = 9,
       height = 11,
       units = 'in',
       dpi = 300)
