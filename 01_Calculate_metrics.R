# Introduction to calculating biodiversity metrics
# Author: Lorena Benitez
# Date: October 2026

###########################################
# Description

# This script starts with taking SEOSAW plot and stem data and cleaning it to
# prepare for calculating biodiversity metrics: species richness, rarified richness, 
# shannon diversity, simpson, inverse simpson, and evenness. 

# The scripts shows how to go from the SEOSAW format to a species abundance matrix 
# which can be used with the vegan package to analyse biodiversity. 

# Later on, it covers how to use ggplot and patchwork to make boxplots to summarise 
# the plot level metrics between sites. 



#############################
# Load packages
############################

# If you don't have any of these packages, use install.packages("_name_") 

library(dplyr)
library(ggplot2)
library(tibble)
library(viridis)
library(tidyr)
library(stringr)

# New package
install.packages('vegan')
library(vegan)


#############################
# Read in data
############################
plots<-read.csv('Data/single_census_plot.csv') 
stems<-read.csv('Data/single_census_stems.csv')


#############################
# Clean data
############################

# For this tutorial, we want to restrict the plots to ones that are the same size
# and make sure that we are only examining tress with dbh >10cm.
# We also only want to look at trees that are alive. 

# There are a lot of columns in both data frames and we won't need all of them. 
# The following code removes plots and columns that we will not be analyzing. 

# Clean plots
p<-plots%>%
  filter(plot_radius==20)%>% # Only keep plots that are the same size
  select(plot_id,n_stems_ge10, plot_area, plot_cluster)# Keep only a few useful columns

# Clean stems
s<-stems%>%
  filter(plot_id%in% p$plot_id)%>% #Filter to remove stems that are in plots that were removed because of plot size differences
  filter(diam_adj>=10)%>% #Only keep stems greater than or equal to 10 as some plots include smaller stems
 # We use diam_adj instead of diam to standarise between trees that have different points of measurement
  filter(stem_status=='a')%>% # Only keep trees that are alive
  select(plot_id,stem_id,species_name_clean,diam_adj)%>%# Keep only a few useful columns
  distinct() # This function is great for removing duplicates. It checks that there are no repeated rows.


#********
#*Q1a: How many plots are we using for this analysis? How many were removed? 
#*Q1b: How many stems are we using for this analysis? How many were removed? 





#############################
# Check taxonomy
############################

# It is good practice to chekc the species names in your data before proceeding
# Misspellings may mean there are duplicates of the same species
# We do not want to count one species as two, so we need to make sure that everything is spelled correctly
# Check species names with the following code
sort(unique(s$species_name_clean))



# The good news about this dataset is that it is pretty clean already, 
# there are no mispellings that lead to a species appearing twice in the dataset. 
# However, there are some individuals that are not determined to species. 
# e.g. "Vitex indet". (indet stands for indeterminate).
# We do not know if Vitex indet in one plot is the same as Vitex indet in another plot. 
# This is a common issue and common practice is to remove indets


# We can filter out rows that have indet in the species column with the following code
s_no_indets <- s %>% filter(!str_detect(species_name_clean, "indet"))

#Compare how many rows are removed
nrow(s)
nrow(s_no_indets)


# For this analysis, we are going to keep the stems not identified to species as we are calculating 
# diversity on a plot level and unidentified species are likely to be different from identified species 
# on the same plot







######################################
#Create abundance matrix
#######################################

# Reformat so that each row is a plot and each column is a species
abundance<-s%>% 
  group_by(plot_id,species_name_clean)%>% # Focus on plots and species data
  summarise(n=n())%>% # Count number of stems per species per plot
  pivot_wider(names_from = species_name_clean,values_from = n)%>% # Change format from long to wide
  mutate(across(everything(), ~replace_na(.x, 0)))%>% # Replace NA with zeros
  column_to_rownames("plot_id")# move plot_id to row names so that entire data frame is numbers

# Remove unidentified genus and species
abundance<-abundance%>%
  select(-`Indet indet`) # SEOSAW uses indet to mean that the species was not identified


# Convert to a matrix
# The vegan package needs the data in this format
comm_matrix<-as.matrix(abundance)


#********
#*Q2: How many species are found across all of the plots?





#######################################
# #Calculate diversity metrics from stems
#######################################
# The vegan package provides separate functions to calculate different 
# diversity metrics. I have written this function to calculate all of 
# them at the same time and put them into a data frame together. 

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

# Use function to calculate metrics
diversity_data<-calc_diversity(comm_matrix)


#################################################
# Calculate rarefied richness
###############################################
# This has to be calculated separately as not all plots have enough stems
# to have rarefied richness calculated for them
# We remove the plots with few 

# Remove plots with too few stems
filtered_matrix <- comm_matrix[rowSums(comm_matrix) >= 10, ]# Remove plots with fewer than 10 stems

# Rarefied species richness
rarefied<-rarefy(filtered_matrix,sample=10)

rare_data<-data.frame(
  plot_id = rownames(filtered_matrix),
  rarefied = rarefied
)


#Join diversity data with rarefied richness
biodiversity_data<-left_join(diversity_data,rare_data,by='plot_id')%>%
                  left_join(p,by='plot_id')#combine with plot metadata


##############################
# Summarise by site
#############################

# We often want to summarise the biodiversity across an area by taking the mean
# from many plots. In the following code, we use plot_cluster to get the means per
# site. You might use this technique to look at biodiversity between different hbabitat
# types or between areas with different human activites or management strategies.

diversity_summary<-biodiversity_data%>%
        group_by(plot_cluster)%>% #Group by site
        summarise(plot_n=n(),# Calculate number of plots being summarised per site
            mean_richness=mean(richness), # Calculate mean
            sd_richness=sd(richness), # Calculate standard deviation
            mean_shannon=mean(shannon,na.rm=TRUE),
            sd_shannon=sd(shannon,na.rm=TRUE),
            mean_simpson=mean(simpson,na.rm=TRUE),
            sd_simpson=sd(simpson),
            mean_invsimpson=mean(invsimpson,na.rm=TRUE),
            sd_invsimpson=sd(invsimpson),
            mean_even=mean(evenness,na.rm=TRUE),
            sd_even=sd(evenness,na.rm=TRUE),
            mean_rarified=mean(rarefied,na.rm=TRUE),
            sd_rarified=sd(rarefied,na.rm=TRUE))

# Save
write.csv(diversity_summary,'Data/diversity_metrics_by_site.csv')

#********
#*Q3: How would you calculate the mean across all sites?
#* 
#* Bonus question: How would you rewrite this using across()?






############################################
# Make plots for each metric
###########################################

# Make boxplot plot for rarefied richness
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

# View
plot_richness



# Make boxplot plot for evenness
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


# Plot rarefied richness


# Plot Simpson


# Plot Inverse Simpson



#*****************************
#* Q4: How would you use facet() to make a plot of all metrics?







###################################################################
# Optional: Use patchwork to make a complex plot
##################################################################

# Patchwork is a package that helps you make complex plots 
# It has good documentation here to help you make publication ready plots
# https://patchwork.data-imaginist.com/

# Install patchwork
install.packages('patchwork')
library(patchwork)

# Use patchwork to make different plots into one figure
# Play around with it as you see fit!
final_plot<-(plot_richness+plot_even)/(plot_rare)+
                                     plot_layout(guides = 'collect')+ # Combine keys into one
                                    plot_annotation(tag_levels = 'A') # Add leters for each panel


final_plot 


# Save final plot using ggplot
ggsave('diversity_metrics.png', final_plot, bg='transparent',
       width = 9,
       height = 11,
       units = 'in',
       dpi = 300)
