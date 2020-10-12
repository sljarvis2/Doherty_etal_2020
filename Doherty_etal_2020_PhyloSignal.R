##### Assessing Phylogenetic Signal #####
### By Stacey J. Doherty ###
### last modified 8 October 2020 ###

## This script was used to assess phylogenetic signal for various abiotic parameters and 16S rRNA data.
## This script is intended to accompany the publication Doherty et al. 2020. "The transition from stochastic to deterministic bacterial community assembly during permafrost thaw succession"
## Some steps are computationally intensive, I recommend running on a super computer.

## Load packages
library(analogue)
library(vegan)
library(openxlsx) 
library(phyloseq)
library(adephylo)
library(dplyr)

## Set working directory
setwd("/path_to_files")

## Load data
#Phylogenetic tree - .nwk format
my_tree <- ape::read.tree("tree.nwk")

# ASV table - .csv format, samples as rows, ASVs as columns (had to transpose the dataframe to fit these requirements)
my_ASVs <- read.csv("feature-table.csv", header=T, row.names=1)
my_ASVs <- t(my_ASVs)

# Environmental data - .xls format, columns are the abiotic parameter, rows are the samples, content will be the measured values for each abiotic parameter.
# This example is assessing pH, gravimetric water content, % nitrogen, % carbon, and C:N
env_data <- openxlsx::read.xlsx("all_env_parameters.xlsx", colNames = TRUE, rowNames = TRUE)

env_data2 <- as.matrix(sapply(env_data[,1:5], as.numeric)) # need to change values to fit how many columns of data you have 

## Calculate the abundance-weighted mean environmental optima for each ASV

# Calculate the weighted average optima for each species for each environmental variable
opt_pH <- optima(my_ASVs, env_data2[,1]) 
opt_pH <- as.data.frame(opt_pH)

# Repeat for the other environmental variables
opt_gwc <- optima(my_ASVs, env_data2[,2]) 
opt_gwc <- as.data.frame(opt_gwc)

opt_nitrogen <- optima(my_ASVs, env_data2[,3]) 
opt_nitrogen <- as.data.frame(opt_nitrogen)

opt_carbon <- optima(my_ASVs, env_data2[,4]) 
opt_carbon <- as.data.frame(opt_carbon)

opt_cn <- optima(my_ASVs, env_data2[,5]) 
opt_cn <- as.data.frame(opt_cn)

## For each environmental variable, create a distance matrix of between ASV differences
dist_pH <- dist(opt_pH, method = "manhattan")

# Repeat for the other environmental variables
dist_gwc <- dist(opt_gwc, method = "manhattan")

dist_nitrogen <- dist(opt_nitrogen, method = "manhattan")

dist_carbon <- dist(opt_carbon, method = "manhattan")

dist_cn <- dist(opt_cn, method = "manhattan")

## Generate euclidean distance matrix that incorporates all physiochemical axes (Dini-Andreote et al., 2015)

# First convert row names to a column in all of the opt dataframes created in lines 40-50

# combine all environmental optima into one dataframe
all_optima <- merge(opt_pH, opt_gwc, by="col1") 
all_optima <- merge(all_optima, next_opt_dataframe) # keep repeating until all are included in one dataframe

# convert column 1 back to row names

## Calculate the between ASV phylogenetic distance matrix - this may take awhile
ASV_distance <- distTips(my_tree, tips = "all", method = "patristic")

# Normalize distances so that they range from 0-1 
ASV_distance_mat <- as.matrix(ASV_distance)
xmax = max(ASV_distance_mat)
xmin = min(ASV_distance_mat)
ASV_distance_normalized <- apply(ASV_distance_mat, MARGIN = c(1,2), FUN = function(x) {((x-xmin)/(xmax-xmin))})

## Generate mantel correlogram for each environmental variable - this may take awhile
mant_pH <- mantel.correlog(dist_pH, ASV_distance_normalized, n.class = 50, r.type="pearson", cutoff = FALSE)

mant_gwc <- mantel.correlog(dist_gwc, ASV_distance_normalized, n.class = 50, r.type="pearson", cutoff = FALSE)

mant_nitrogen <- mantel.correlog(dist_nitrogen, ASV_distance_normalized, n.class = 50, r.type="pearson", cutoff = FALSE)

mant_carbon <- mantel.correlog(dist_carbon, ASV_distance_normalized, n.class = 50, r.type="pearson", cutoff = FALSE)

mant_cn <- mantel.correlog(dist_cn, ASV_distance_normalized, n.class = 50, r.type="pearson", cutoff = FALSE)

# Plot the correlogram
mant_pH_plot <- plot(mant_pH)

mant_gwc_plot <- plot(mant_gwc)

mant_nitrogen_plot <- plot(mant_nitrogen)

mant_carbon_plot <- plot(mant_carbon)

mant_cn_plot <- plot(mant_cn)

 