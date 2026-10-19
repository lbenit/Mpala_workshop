# Rarefaction
# Author: Lorena Benitez
# Date: October 2026

###########################################
# Description

# This script shows how to create rarefaction curves using two packages: vegan 
# and iNext. Rarefaction is useful for understanding how sampling affects
# observed species richness. 
# We rarefy to 10 stems in this script, removing plots with fewer stems. 
# It is important to remember how many stems you are rarefying to in this analysis
# for some of the questions.

#############################
# Load packages
############################

# If you don't have any of these packages, use install.packages("_name_") 

library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(vegan)

# New package
install.packages('iNext')
library(iNEXT)

#############################
# Read in data
############################
plots<-read.csv('Data/single_census_plot.csv') 
stems<-read.csv('Data/single_census_stems.csv')


#############################
# Clean data
############################

# This script starts the same as the previous one.


# Clean plots
p<-plots%>%
  filter(plot_radius==20)%>% # Only keep plots that are the same size
  filter(n_stems_ge10>=10)%>% # Only keep plots with more than 10 stems 
  select(plot_id,n_stems_ge10, plot_area, plot_cluster)# Keep only a few useful columns

# Clean stems
s<-stems%>%
  filter(plot_id%in% p$plot_id)%>% #Filter to remove stems that are in plots that were removed because of plot size differences
  filter(diam_adj>=10)%>% #Only keep stems greater than or equal to 10 as some plots include smaller stems
  filter(stem_status=='a')%>% # Only keep trees that are alive
  select(plot_id,stem_id,species_name_clean,diam_adj)%>%# Keep only a few useful columns
  distinct() # This function is great for removing duplicates. It checks that there are no repeated rows.



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

# Remove unidentified species
abundance<-abundance%>%
  select(-`Indet indet`) # SEOSAW uses indet to mean that the species was not identified


# Convert to a matrix
comm_matrix<-as.matrix(abundance)

# Remove plots with fewer than 10 individuals now that indet is removed
filtered_matrix <- comm_matrix[rowSums(comm_matrix) >= 10, ]# Remove plots with fewer than 10 stems



######################################
# Generate plot level rarefaction curves
#######################################
rarecurve(filtered_matrix, sample=10)

# The resulting plot has many horizontal lines



#*****************
#*Q1: The resulting plot has many horizontal lines.
#* Why is this happening and what does it mean?
#* How would you go about removing them?







#########################################################################
# Use iNext to look at site level rarefaction curves
#########################################################################

# Description of iNext

# iNext is an R package that is specifically designed for rarefaction analysis. 
# It has several advantages over vegan including:
# - It can extrapolate to predict the diversity with more plots/stems
# - It allows easy analysis of site level diversity and comparison
# - It allows calculation of rarefaction using species richness, shannon diveristy,
#     and Simpson diversity, while vegan only uses species richness
# - It automatically generates confidence intervals
# - It automatically uses ggplot, so making nice figures is a lot easier

# You can take a look at this tutorial for more information of iNext
# https://johnsonhsieh.github.io/iNEXT/inst/doc/Introduction.html



###############################
# For this tutorial, we are using iNext to compare the sampling and diversity
# between our three sites.
# iNext needs a different data fomrat than vegan and we need to retain the 
# information about which plots belong to which sites


# First we keep infomration about the sites
plot_cluster<-p%>%select(plot_id,plot_cluster)


# Then, we simply the stem data in a long format
# This gives us a data frame of plots by species like we had before
species<-as.data.frame(table(s$plot_id,s$species_name_clean))


# Modify species data and add plot cluster assignments
species<-species%>%
  rename(plot_id=Var1, species=Var2)%>%
  left_join(plot_cluster,by='plot_id')%>%
  relocate(plot_cluster, .after=plot_id)


# Group data by site
species_long<-species%>%
  group_by(species, plot_cluster)%>%
  summarise(n=sum(Freq))

# Change shape to wide format
species_wide<-species_long%>%
  pivot_wider(names_from = plot_cluster,values_from = n)%>%
  column_to_rownames('species')


##############################
# Caclulate rarefaction
##############################
D_abund <- iNEXT(species_wide,
                 q=1, # Change to alter diversity metric 
                 datatype = 'abundance', 
                 nboot=100)

plot(D_abund)

#*****************
#*Q2: iNext allows you to examine diversity using three metrics.
#*    Change q in the code above to compare different metrics. 
#*    1 is richness, 2 is shannon, and 3 is simpson
#*    How do they differ?












######################################
# Visualisation 
######################################

# iNext allows for the data to be visualised using three methods
# 1. Sample-size-based R/E curve (type=1). 
# This curve plots diversity estimates with confidence intervals (if se=TRUE) as a 
# function of sample size up to double the reference sample size, by default, 
# or a user‐specified endpoint.

# 2. Sample completeness curve (type=2) with confidence intervals. 
# This curve plots the sample coverage with respect to sample size for the same 
# range described in (1).

# 3. Coverage-based R/E curve (type=3). 
# This curve plots the diversity estimates with confidence intervals as a function 
# of sample coverage up to the maximum coverage obtained from the maximum size described in (1).



ggiNEXT(D_abund, type=1)
ggiNEXT(D_abund, type=2)
ggiNEXT(D_abund, type=3)


#*****************
#*Q3: What do these plots mean and how do they differ?







#####################################
# Example: Modify iNext plot using ggplot
#####################################

# Specify plot type
plot<-ggiNEXT(D_abund, type=1)

# Change colors and annotate
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


# View plot
plot1





