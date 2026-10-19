# Ordination and clustering techniques
# Author: Lorena Benitez
# Date: October 2026

###########################################
# Description

# Ordination is used to reduce the complexity of ecological 
# data down to fewer dimensions. In this case, we are
# reducing the species abundances to identify ecological communities. Once we 
# have reduced the data down to a few dimensions, we use clustering to identify 
# the community (aka clusters).

# This script compares three ordination techniques and two clustering methods.

# 1. Principal Components Analysis (PCA):
# A linear method of dimension reduction based on the distance of plots
# from one another.


# 2. Non-metric Multidimensional Scaling (NMDS):
# Is a ranked-based, non-linear method. This means that it doesn't rely on the 
# distance between plots as much as which plots are closer to one another. 
# You need a dissimilarity matrix to calculate NMDS.
# It determines is plot A is closer to plot B or C in the original high dimensional space.
# And if this can be preserved in fewer dimensions (usually 2 or 3). 



# 3. Uniform Manifold Approximation and Projection (UMAP):
# UMAP is a non-linear dimension reduction technique. It is similar to PCA, 
# but instead of fitting data to straight line relationships, it allows for 
# more complex representations of data patterns. 



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
plots<-read.csv('Data/single_census_plot.csv') 
stems<-read.csv('Data/single_census_stems.csv')


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
# Hellinger is particularly good for data with a lot of zeros
# It divides all values by the row sum and takes the square root of the values
bio_hel <- decostand(filtered_matrix, method = "hellinger")


# Save file as we will need it in the next tutorial
write.csv(bio_hel,'Data/hellinger_transformed_data.csv')




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


# Make dataframe with plot scores
pca_site<-as.data.frame(pca_result$x)




# Identify ideal number of clusters
fviz_nbclust(pca_site, kmeans, method = "wss", k.max=20) + 
  labs(subtitle = "Elbow Method")
# Look at the plot and determine at which number of clusters the curve
# goes from being vertical to being horizontal
# This is called the elbow method


# Specify the number of clusters
kmeans_res <- kmeans(pca_result$x, 
                     centers = 3, # Change this to alter the number of clusters
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






####
# Save files with plot and species scores
# Prep species data
pca_species<-as.data.frame(pca_result$rotation)
pca_species$species<-rownames(pca_species)
write.csv(pca_species,'Data/pca_species_scores.csv')


# Prep plot data

# Add cluster assignments with your chosen K
pca_site$cluster<-as.factor(kmeans_res$cluster)

pca_site$plot_id<-rownames(pca_site)
write.csv(pca_site,'Data/pca_site_scores.csv')












##################################################################
# Ordination technique 2:  Non-metric Mulitdimensional Scaling (NMDS)
##################################################################

# Calculate distances between plots 
# vegan needs the data in this format
bio_distmat <- 
  vegdist(bio_hel, method = "bray")


# Run NMDS in vegan (metaMDS)
# This may take a few minutes to run
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
  geom_point(data = nms_scores,
             aes(NMDS1, NMDS2, color = cluster),
             size = 2) +
  # add species labels
  geom_text(data = nms_species,
            aes(NMDS1, NMDS2, label = species),
            color = "black", size = 3, vjust = -0.5) +
  theme_minimal() 


# Save these for further analyses
write.csv(nms_scores,file='Data/nmds_scores.csv')
write.csv(nms_species,file='Data/nmds_species.csv')



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
  n_neighbors = 100, 
  min_dist = 0.5, 
  n_components = 2 # Number of dimensions to reduce to: use either 2 or 3
)


# Extract results to dataframe
umap_data <- as.data.frame(umap_result$layout)
colnames(umap_data) <- c("UMAP1", "UMAP2")
umap_coords$plot_id <- rownames(umap_data)

# Plot umap
ggplot(umap_data, aes(x = UMAP1, y = UMAP2)) +
  geom_point(color = "darkgreen", size = 2) +
  labs(title = "UMAP Projection (Bray-Curtis Distance)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


#################
# Cluster plots using hierarchical clustering
#################
# K-means clustering performs poorly umap analysis
# Instead, we use hierarchical clustering to cluster the plots


# Select umap scores scores
umap_coords<-umap_data[,1:2]


#Hierarchical clustering of umap scores
dist_umap<-dist(umap_coords)# first compute the distance between points
hc <- hclust(dist_umap, method='ward.D2')# Then cluster the distances

# Plot the clustering tree
plot(hc)

# add lines to show where breaks for four and six are
rect.hclust(hc, k = 4, border = "red")
rect.hclust(hc, k = 6, border = "blue")


##########
# Determine the optimal number of clusters using the average silhouette scores
############

# In this section, we are calculating the mean silhouette width of each cluster
# for clusters numbered 2 to 15. The higher mean silhouette width is, the better 
# defined and separated the clusters are. 

# Initialize vector to store scores
sil_width_hc <- c()

# Loop through k = 2 to 15 (adjust max if you expect more clusters)
max_k <- 15 #chnage number to adjust how many clusters to test

# Run loop to get scores for each number of clusters
for(i in 2:max_k){
  # Cut the tree to get cluster assignments
  clusters_i <- cutree(hc, k = i)
  
  # Calculate silhouette width using the UMAP distance matrix
  sil <- cluster::silhouette(clusters_i, dist_umap)
  
  # Store the mean silhouette width
  sil_width_hc[i] <- mean(sil[, 3])
}

# Remove the first element (NA from index 1 since we started loop at 2)
sil_width_hc <- sil_width_hc[!is.na(sil_width_hc)]


# Plot the silhouette score results and use the elbow method to determine the 
# best number of clusters
plot(2:max_k, sil_width_hc,
     type = "b", pch = 19,
     xlab = "Number of clusters (k)",
     ylab = "Average Silhouette Width",
     main = "Optimal Clusters- use elbow method")



#****************************
#*Q4: What is the best number of clusters according to the silhouette scores?
#* Hint: Which numbers of clusters has the highest average silhouette width?











# Determine where to cut the trees i.e. how many clusters do you want?
# Input the number of clusters you have selected
clusters <- cutree(hc, k = ) # Fill in k =  with your choice

# Make sure  that cluster is a factor
umap_data$cluster <- as.factor(clusters)
umap_data$plot_id<-rownames(umap_data)


# Plot UMAP space with clusters
ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 2) +
  guides(color = guide_legend(override.aes = list(size = 8)))+
  theme_minimal()




##############
# Get species scores
##############
species_scores <- wascores(umap_data[,1:2], bio_hel)
spscores <- data.frame(
  UMAP1 = species_scores[,1],
  UMAP2 = species_scores[,2],
  species=colnames(bio_hel)
)


# Plot space with species labelled
library(ggrepel) # install.packages('ggrepel') if necessary
# this package keeps the species names from overlapping


# Plot with species scores
ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 2) +
  guides(color = guide_legend(override.aes = list(size = 8))) +
  geom_text_repel(data = spscores, aes(x = UMAP1, y = UMAP2, label=species),
                  color = "black", size = 5)+
  theme_bw() +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Clusters"   # legend title
  ) 



# Save species and plot scores with cluster assignments
write.csv(spscores,'Data/umap_species_scores.csv')

write.csv(umap_data,'Data/umap_plot_scores.csv')







#****************************
#*Q5: Which method did you prefer? Why?


