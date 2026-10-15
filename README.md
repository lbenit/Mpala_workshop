This repo contains scripts for the biodiversity and floristics day of the SEOSAW Mpala WOrkshop (2026). 



There are four scripts and associated data:

01_Calculate_metrics
This script can be used to calculate basic bidiversity metrics like species richness, rarified richness, and shannon diversity. 

02_Rarefaction
This script shows how to create rarefaction curves using two packages: vegan and iNext.

03_Ordination_clustering
Ordination is used to reduce the complexity of ecological data down to fewer dimensions. In this case, we are reducing the species abundances to identify ecological communities. Once we have reduced the data down to a few dimensions, we use clustering to identify the community (aka clusters). This script compares three ordination techniques and two clustering methods.

04_Indicator_analysis
This script is used to identify the indicator species within each cluster.
It allows for the clusters to be identified as ecological communities with indicator statistics for each species that describes how influential it is in each community. 


The first three scripts can be run independently, but do not run script 4 withour first running script 3.
