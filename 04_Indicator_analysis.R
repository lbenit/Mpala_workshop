# Identifying indicator genera
# Author: Lorena Benitez
# Date: October 2026

###########################################
# Description

# This script is used to identify the indicator species within each cluster.
# It allows for the clusters to be identified as ecological communities with
# indicator statistics for each species that describes how influential it is
# in each community. 


#############################
# Load packages
############################

# Install new package
install.packages('indicspecies')

library(indicspecies)
library(dplyr)


#############################
# Read in data
############################
# Read in plot metadata
plots<-read.csv('Data/single_census_plot.csv') 

# Load in community matrix
bio_hel<-read.csv('Data/hellinger_transformed_data.csv')%>%
  select(-X)


####################################
# Select which ordination and clustering method you want to use
# Important: Only select one!!!!!
#####################################

# Option 1: PCA
sites<-read.csv('Data/pca_site_scores.csv',row.names = NULL )%>%
  select(-X)
species<-read.csv('Data/pca_species_scores.csv',row.names = NULL )%>%
  select(-X)





# Option 2: NMDS
sites<-read.csv('Data/nmds_scores.csv',row.names = NULL )%>%
  select(-X)
species<-read.csv('Data/nmds_species.csv',row.names = NULL )%>%
  select(-X)




#Option 3: UMAP
sites<-read.csv('Data/umap_plot_scores.csv',row.names = NULL )%>%
  select(-X)
species<-read.csv('Data/umap_species_scores.csv',row.names = NULL )%>%
  select(-X)









#############################################
# Identify indicators
############################################

# Look at this R vignette if you want to learn more about indicspecies package
# https://cran.r-project.org/web/packages/indicspecies/vignettes/IndicatorSpeciesAnalysis.html


# Indicator species analysis
ind <- multipatt(
  x = bio_hel,
  cluster = sites$cluster,
  func = "r.g" # Change this to alter the calculation method
)

# Options for func are as follows:
# IndVal = Tests the association of species presence with groups
# Indval.g = Same as IndVal, but corrects for a different number of sites in each group
# r = Test the correlation between species abundance and groups
# r.g = Same as r, but corrects for a different number of sites in each group


# Look at summary of scores by cluster
summary(ind)




#*****************************
#* Q1: Look at the summary of indicator scores by cluster. 
#*     Do the top indicators make sense based on what you know 
#*     about southern African tree communities?














# Save indicator scores
ind_table <- ind$sign

write.csv(ind_table, 'indicator_statistics.csv') #change name based on method



####################################
# Reorganize data for saving
###################################

# Assign your clusters to an ecologically meaningful name
# I usually look at the the first few indicators
plot_assignments <- sites %>%
  mutate(
    `Floristic group` = recode(cluster, # Add numbers based on number of clusters you have
                               
                               `1` = "", # Fill in with your description based on cluster 1 indicators
                               `2` = "", # Fill in with your description based on cluster 2 indicators
                               `3` = "", # Fill in with your description based on cluster 3 indicators
                               `4` = ""))# Fill in with your description based on cluster 4 indicators

# Get plot_ids
plot_assignments$plot_id<-row.names(plot_assignments)
row.names(plot_assignments)<-NULL

#Save
write.csv(plot_assignments, 'plot_assignments.csv') #change file name based on selection


#################################
# Plot ordination space with assigned cluster
##################################

ggplot() +
  geom_point(data = plot_assignments,
             aes(x=UMAP1, #Change based on the technique you chose
                 y=UMAP2, #Change based on the technique you chose
                 color = `Floristic group`),
             size = 3) +
  theme_classic()+
  theme(
    legend.text = element_text(size=20),
    legend.title = element_text(size=30),
    legend.key.size = unit(1, 'cm'),
    axis.title=element_text(size=20),
    axis.text= element_text(size=18))


#*****************************
#* Q2: Do the communities make sense when plotted?







###################################################################
# Optional
##################################################################

# Identify the indicators using another method and compare them to the first.
# Do they result in similar communities? If not, think about why that might be. 


# Compare the clusters to the plot_cluster column in plots. 
# Are the clusters found at all sites or are they geographically separated?



