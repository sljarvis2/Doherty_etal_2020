##### Community composition - ordination and statistical analysis #####
### By Stacey J. Doherty ###
### last modified 8 October 2020 ###

## This script was used to assess beta diversity of microbial community data. Script includes NMDS ordination, PerMANOVA, and dispersion test.
## This script is intended to accompany the publication Doherty et al. 2020. "The transition from stochastic to deterministic bacterial community assembly during permafrost thaw succession"

## Load packages

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("remotes")
#remotes::install_github("jbisanz/qiime2R")

library(vegan)
library(dplyr)
library(ggplot2)
library(readr)
library(rlist)
library(tibble)
library(RColorBrewer)
library(qiime2R)
library(ggrepel)

## Import data

# Set the working directory containing files
setwd("/path_to_folder")

# Read in feature table, note I imported .qza files from qiime. 
asv <- read_qza("feature_table.qza")
asv_table <- asv[["data"]]
asv_table <- t(asv_table)
dim(asv_table)

# Read in metadata, note I imported .tsv file from qiime.
metadata <- read_q2metadata("metadata.tsv")
metadata <- tibble::column_to_rownames(metadata, "SampleID")

# Need to remove rows from metadata containing samples not found in the table!
metadata <- rownames_to_column(metadata, 'ID')
metadata <- filter(metadata, metadata$ID %in% rownames(asv_table))
metadata <- column_to_rownames(metadata, 'ID')

##### Ordination with NMDS #####

## Run NMDS with 2 and 3 dimensions to see which is best based on stress.

NMDS_depth_2 <- metaMDS(asv_table, distance = "bray", k = 2, try = 20, autotransform = FALSE)

NMDS_depth_3 <- metaMDS(asv_table, distance = "bray", k = 3, try = 20, autotransform = FALSE)

# check the stress of each solution and choose which is the best solution (stress  <0.20)
NMDS_depth_2$stress
NMDS_depth_3$stress

# Stress evalution:
# < 0.05 - excellent
# 0.05-0.10 - good
# 0.10-0.20 - usable, potential to result in misinterpretation at upper end

# Let's move forward with the 2D solution...

## Site scores
site.scrs <- as.data.frame(scores(NMDS_depth_2, display = "sites"))
site.scrs <- merge(site.scrs, select(metadata, depth_cm), by.x = "row.names", by.y = "row.names")
site.scrs

# Plot!
NMDS_1 <- ggplot(site.scrs, aes(x = NMDS1, y = NMDS2)) + geom_point(aes(NMDS1, NMDS2, colour = factor(site.scrs$depth_cm)), size = 4) + labs(colour = "Depth") + theme_classic() + scale_color_brewer(palette="Spectral")
NMDS_1

## Add environmental variables as vectors

# Relativize the environmental data. Log transformation was used
depth_env_num_stand <- decostand(select_if(metadata, is.numeric), method = "log")

# Calculate environmental vectors
depth_envfit <- envfit(NMDS_depth_2, depth_env_num_stand, choices = 1:2, permutations = 999)
env.scores <- as.data.frame(scores(depth_envfit, display = "vectors"))
env.scores2 <- cbind(env.scores, env.variables = rownames(env.scores), pval = depth_envfit$vectors$pvals)
env.scores2

# Plot!
NMDS_1_with_env <- NMDS_1 + geom_segment(data = env.scores2, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd = 0.3) + ggrepel::geom_text_repel(data = env.scores2, aes(x = NMDS1, y = NMDS2, label = env.variables), cex = 4, direction = "both", segment.size = 0.25)
NMDS_1_with_env

##### Conduct PerMANOVA to test if groups are significant ####
# convert ASV table as needed
asv_table_par <- as.data.frame(asv_table)
asv_table_par <- tibble::rownames_to_column(asv_table_par, "SampleID")

# combine metadata and ASV table into one dataframe
metadata2 <- tibble::rownames_to_column(metadata, "SampleID")
my.final.data <- inner_join(asv_table_par, metadata2, by = "SampleID")
my.final.data <- my.final.data[, c(2842:2845, 1:2841)] #these numbers will change depending on your data

# make matrix with only ASV counts!
my.final.data.ASV <- as.matrix(sapply(my.final.data[,6:2845],as.numeric)) #these numbers will change depending on your data

# PerMANOVA, plot was used as strata in checking differences by depth.
adonis(formula = my.final.data.ASV~depth_cm, data = my.final.data, strata = my.final.data$plot, distance="bray", permutations = 9999)

##### Dispersion test #####

# Confirm order of rows matches between ASV table and metadata. This is needed when we define the groups below
rownames(asv_table) == rownames(metadata)

# Calculate Bray-Curtis dissimilarity between samples
asv_dist <- vegdist(asv_table, method = "bray")


## Look at ASV table and assign each row to its appropriate factor, in this case it was depth (i.e. 1 through 8 are the different levels)
depth_groups <- factor(x = c(1,2,4,8,1,2,3,4,5,6,7,8,1,2,3,5,6,7,8,2,3,4,5,6,7,8), labels = c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80"))
depth_groups
#confirm you did this successfully by visually checking ASV table and group output to the console.

# Calculate multivariate dispersions. See warning for betadisper() on why the type should be median, not centroid.

mod <- betadisper(asv_dist, depth_groups, type = "median")
mod

## Permutation test for F
pmod <- permutest(mod, permutations = 999, pairwise = TRUE)
pmod




