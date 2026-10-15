# Ordination and clustering techniques
# Author: Lorena Benitez
# Date: October 2026

###########################################
# Description

# This script 


#############################
# Load packages
############################

# If you don't have any of these packages, use install.packages("_name_") 

library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(vegan)
library(stringr)

# New package
install.packages('factoextra')
library(factoextra)

install.packages('umap')
library(umap)


#############################
# Read in data
############################
plots<-read.csv('single_census_plot.csv') 
stems<-read.csv('single_census_stems.csv')


#############################
# Clean data- same as last two scripts
############################

# Clean plots
p<-plots%>%
  filter(plot_radius==20)%>% # Only keep plots that are the same size
  select(plot_id,n_stems_ge10, plot_area, plot_cluster)# Keep only a few useful columns

# Clean stems
s<-stems%>%
  filter(plot_id%in% p$plot_id)%>% #Filter to remove stems that are in plots that were removed because of plot size differences
  filter(diam_adj>=10)%>% #Only keep stems greater than or equal to 10 as some plots include smaller stems
  filter(stem_status=='a')%>% # Only keep trees that are alive
  select(plot_id,stem_id,species_name_clean,diam_adj)%>%# Keep only a few useful columns
  distinct() # This function is great for removing duplicates. It checks that there are no repeated rows.


######################################
#Create abundance matrix- same as last two scripts
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



######################################
# Prepare data for ordination: Really important!
######################################

#Remove rare species
# This is a good idea because rare species make the ordination space unstable. 
# You don't have to do this, but 
yn<-abundance
yn[yn > 0] <- 1 
keep<-as.data.frame(colSums(yn))%>%filter(colSums(yn)>1)# Change number to adjust threshold for inclusion
keep<-rownames(keep)
keep

# Remove rare species from abundance dataframe
abund<-abundance[,keep]


# Convert to a matrix
# The vegan package needs the data in this format
comm_matrix<-as.matrix(abund)


# Filter out rows equal to zero
# If you don't do this, the analysis will fail
filtered_matrix <- comm_matrix[rowSums(comm_matrix) > 0, ]# Remove plots with fewer than 10 stems


# Convert values using Hellinger method
# Hellinger...
bio_hel <- decostand(filtered_matrix, method = "hellinger")







##################################################################
# Ordination technique 1: Principal Components Analysis (PCA)
##################################################################

# Run PCA
pca_result <- prcomp(bio_hel, center = TRUE, scale. = FALSE)

# Look at results
summary(pca_result)

# Assess how much variation is described by the first few components
fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50))


#*********************
#*Q1: How many components does it take to represent 50% of the variation in 
#*    the plot data?
#*    Hint: Look at the Cumulative proportion from the summary. 



#########################
# Clustering with k-means
#########################


# Plot first two components
fviz_pca_ind(pca_result,
             geom.ind = "point", # Show points only
             title = "PCA of Community Matrix (Hellinger Transformed)",
             axes=c(1,2))# Specify which components to plot


# Identify ideal number of clusters
fviz_nbclust(bio_hel, kmeans, method = "wss") + 
  labs(subtitle = "Elbow Method")
# Look at the plot and determine at which number of clusters the curve
# goes from being vertical to being horizontal
# This is called the elbow method


# Specify the number of clusters
kmeans_res <- kmeans(pca_result$x[, 1:2], 
                     centers = 10, # Change this to alter the number of clusters
                     nstart = 25)

# Plot first two components with 
fviz_pca_ind(pca_result,
             col.ind = as.factor(kmeans_res$cluster), # Color by cluster
             palette = "jco",
             addEllipses = TRUE,
             ellipse.type = "convex",
             legend.title = "Clusters",
             title = "PCA with K-Means Clusters",
             axes=c(1,2)# Specify which components to view
             )



#*Q2: Test out a few different K numbers. Which visually looks the best?
#*    There is not a correct answer with this one, just a thought exercise.










##################################################################
# Ordination technique 2:  Non-metric Mulitdimensional Scaling (NMDS)
##################################################################

# Calculate distances between plots 
# vegan needs the data in this format
bio_distmat <- 
  vegdist(bio_hel, method = "bray")


# Run NMDS in vegan (metaMDS)
bio_NMS <-
  metaMDS(bio_distmat,
          k = 3,
          maxit = 999, 
          trymax = 100)

#######
# Test fit of model
########
# Stats
bio_NMS
stressplot(bio_NMS)



###################
# Plot NMDS scores
##################
# Specify species
sppscores(bio_NMS) <- bio_hel

# View results
plot(bio_NMS)

# Add labels for species
orditorp(bio_NMS, "species") 

# Add labels for plots
#orditorp(bio_NMS,'sites')




#######################################
# Clustering with k-means using base r
########################################
# Create dataframe with site scores
nms_scores <- as.data.frame(scores(bio_NMS, display = "sites"))

# Test 
wss <- sapply(1:10, function(k)
  kmeans(nms_scores, centers = k, nstart = 20)$tot.withinss)

# Elbow method with following plot to determine the best k number
plot(1:10, wss, type = "b",
     xlab = "Number of clusters (k)",
     ylab = "Within-cluster sum of squares")

# Specify number of clsuters (k) and assign plots to clusters
set.seed(123)   # for reproducibility
km <- kmeans(nms_scores, centers = 6)
nms_scores$cluster <- as.factor(km$cluster)
nms_scores$plot_id<-row.names(nms_scores) #Make plot_id a column


# Create dataframe with scores for species
nms_species <- as.data.frame(scores(bio_NMS, display = "species"))
nms_species$species <- rownames(nms_species)


# Create dataframe with scores for sites
nms_sites <- as.data.frame(scores(bio_NMS, display = "sites"))
nms_sites$cluster <- as.factor(km$cluster)


# Plot first two dimensions using ggplot
ggplot(nms_scores, aes(NMDS1, NMDS2, color = cluster)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "NMDS with k-means Clusters")


# Add species names to plot
ggplot() +
  geom_point(data = nms_sites,
             aes(NMDS1, NMDS2, color = cluster),
             size = 2) +
  # add species labels
  geom_text(data = nms_species,
            aes(NMDS1, NMDS2, label = species),
            color = "black", size = 3, vjust = -0.5) +
  theme_minimal() 


# Save these for further analyses
write.csv(nms_scores,file='nmds_scores.csv')
write.csv(nms_species,file='nmds_species.csv')



#****************************
#*Q3: How does the NMDS result compare to the PCA? Which one works better for the data?









##################################################################
# Ordination technique 3:  Uniform Manifold Approximation and Projection (UMAP)
##################################################################

# UMAP is a relatively new ordination technique. Like NMDS<,it is a non-linear method,
# but it has the advantage of being much more computationally efficient.

# UMAP was originally designed for python. The python package has really good documentation,
# which I recommend looking at to understand how UMAP works and the different functions available.
# https://umap-learn.readthedocs.io/en/latest/


# The umap r package is a basic version of the python package. It doesn't have all of the functions,
# but the basics are good. 
# Tutorial using umap r package: https://cran.r-project.org/web/packages/umap/vignettes/umap.html


##################
# Run umap using distance matrix calculated earlier
##################
# Convert to matrix
dist_matrix<-as.matrix(bio_distmat)

# Run umap
umap_result <- umap(
  dist_matrix, 
  input = "dist", 
  n_neighbors = 15, 
  min_dist = 0.1, 
  n_components = 2 # Number of dimensions to reduce to: use either 2 or 3
)


# Extract results to dataframe
umap_data <- as.data.frame(umap_result$layout)
colnames(umap_data) <- c("UMAP1", "UMAP2")
umap_coords$plot_id <- rownames(umap_data)

# Plot umap
ggplot(umap_coords, aes(x = UMAP1, y = UMAP2)) +
  geom_point(color = "darkgreen", size = 2) +
  labs(title = "UMAP Projection (Bray-Curtis Distance)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


#################
# Cluster plots using hierarchical clustering
#################
# K-means clustering performs poorly umap analysis
# Instead, we use hierarhcical clustering to cluster the plots


# Select only scores
umap_coords<-umap_data[,1:2]


#Hierarchical clustering of umap scores
hc <- hclust(dist(umap_coords))

# PLot the clustering tree
plot(hc)
rect.hclust(hc, k = 5, border = "red")
rect.hclust(hc, k = 6, border = "blue")
rect.hclust(hc, k = 7, border = "cyan")


# Determine where to cut the trees
clusters <- cutree(hc, k = 4)
df$cluster <- as.factor(clusters)



#####
# Visualize in 3D
plot_ly(umap_data, x = ~x, y = ~y, 
        type = "scatter3d", mode = "markers",color= ~cluster,colors = c("red",'yellow','orange',
                                                                        'purple','pink', 'darkred','brown','deeppink',
                                                                        'darkmagenta','coral1'),
        marker = list(size = 4))

##############
# Get genus scores
species_scores <- wascores(umap_data[,1:2], bio_hel)
spscores <- data.frame(
  UMAP1 = species_scores[,1],
  UMAP2 = species_scores[,2],
  species=colnames(bio_hel)
)



# PLotting
library(ggrepel)
ggplot(df, aes(x = x, y = y, color = cluster)) +
  geom_point(size = 1.5, alpha = 0.5) +
  
  scale_color_viridis_d(option = "turbo") +
  guides(color = guide_legend(override.aes = list(size = 8))) +
  geom_text_repel(data = important3, aes(x = x, y = y, label=genus),
                  color = "black", size = 5)+
  theme_bw() +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Clusters"   # legend title
  ) +
  theme(
    legend.position = "right",   # move legend if needed
    axis.title = element_text(size = 50),                  # axis titles
    axis.text = element_text(size = 40),                   # axis tick labels
    legend.text = element_text(size = 30),                 # legend labels
    legend.title = element_text(size = 40)  ,
    legend.key.size = unit(3, "lines")
    
  )



