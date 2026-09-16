# (BIO 5150L) Optional PCA practice: varechem soil data
#
# Use this starter script with the PCA tutorial as a guide. Copy, paste, and
# adapt the relevant sections of the tutorial to use this dataset instead.
#
# Data documentation:
# https://vegandevs.github.io/vegan/reference/varechem.html


# Load packages -----------------------------------------------------------

library(tidyverse)
library(vegan)
library(ggrepel)

# The tutorial also uses corrplot::corrplot(). Install corrplot once if needed:
# install.packages("corrplot")


# Load and prepare the data ----------------------------------------------

data("varechem", package = "vegan")

dat_varechem <- varechem |>
  as.data.frame() |>
  rownames_to_column(var = "sample_ID") |>
  as_tibble()

glimpse(dat_varechem)
head(dat_varechem)

# All columns except sample_ID are continuous soil measurements. The data are
# already in wide format: one row per site and one column per PCA variable.
pca_variables <- dat_varechem |>
  select(-sample_ID) |>
  names()

pca_variables


# Practice adapting the PCA tutorial -------------------------------------

# Adapt the tutorial code to complete the following steps:
#
# 1. Summarize missing values, distributions, and ranges for all PCA variables.
#
# 2. Create a clustered correlation plot. Which variables contain similar
#    information, and which appear to represent different soil gradients?
#
# 3. Construct the sample-by-variable data frame (i.e., "comm_varechem" table)
# required by vegan, with sample_ID as row names.
#
# 4. Fit a PCA with vegan::rda(). Include the argument to center and scale the
# variables because their units and numerical ranges differ substantially.
#
# 5. Calculate and plot the percentage of variation explained by each PC.
#
# 6. Extract PC1 and PC2 sample and variable scores using a consistent scaling.
#    Add the percentage explained to both axis labels.
#
# 7. Create the ggplot w/ vectors and interpret the primary soil gradients.
# Which variables are most strongly associated with PC1 and PC2?
#
# 8. Join the original measurements to the sample scores. Try mapping a soil
#    variable such as pH or Humdepth to point color or size. These overlays are
#    descriptive because the same variables helped construct the PCA.
#
# 9. Decide whether PC1 or PC2 represents a coherent combined soil metric that
#    could be retained for a later analysis. Support the decision using the
#    variance explained and correlations between the original variables and PCs.
#
# There is no categorical grouping variable in varechem. You can label sites
# with sample_ID for exploration purposes.
