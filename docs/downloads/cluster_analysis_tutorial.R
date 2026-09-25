# (BIO 5150L) Hierarchical cluster analysis tutorial
# Community composition and continuous biological measurements
# Date modified: 25 Sep 2026
#
# This tutorial contains two independent examples:
#
# 1. Coral reef invertebrate assemblages
#    - Bray-Curtis dissimilarity
#    - Average-linkage clustering
#
# 2. Darlingtonia plant morphology
#    - Euclidean distance among centered and scaled measurements
#    - Ward's minimum-variance clustering (method = "ward.D2")
#
# Run either example or compare the results from both. Cluster analysis is an
# exploratory classification method. It identifies patterns in the supplied
# data but does not test whether the resulting groups are statistically
# significant. It can complement ordination methods such as nMDS or PCA,
# particularly when a heatmap is displayed with the cluster dendrograms.
#
# Required files in the same folder:
#   Coral_Reef_Data_2014_to_2025.xlsx
#   Darlingtonia_GE_Table12.1.csv
#
# Additional background:
# https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/types-of-cluster-analyses/
# https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/hierarchical-cluster-analysis/


# Load packages -----------------------------------------------------------

library(tidyverse)
library(readxl)
library(vegan)
library(factoextra) # for comparing mean silhouette width across values of k
library(ggdendro) # for ggplot dendrograms
library(ggalign)  # for heatmaps with dendrograms

# Install packages once if any are not already installed:
# install.packages(c(
#   "tidyverse", "readxl", "vegan", "factoextra", "ggdendro", "ggalign"
# ))


# Choosing the distance and linkage method --------------------------------

# A hierarchical cluster analysis requires two separate choices:
#
# 1. Distance or dissimilarity describes how different two samples are.
# 2. Linkage describes how samples or groups are progressively joined.
#
# Species-abundance data:
# - Use Bray-Curtis after applying the same defensible transformation and
#   rare-taxon filtering used for the nMDS.
# - Average linkage is a common ecological choice and can be used with
#   Bray-Curtis dissimilarities.
#
# Continuous environmental, morphological, or physiological data:
# - Center and scale variables when their units or numerical ranges differ.
# - Calculate Euclidean distances and use method = "ward.D2". Ward's method
#   seeks compact groups by minimizing increases in within-cluster variation.
#
# Do not use Ward's method as the default with Bray-Curtis. Its minimum-
# variance interpretation depends on Euclidean geometry. Other linkage methods
# can also be considered (see other resources for options). Always report both
# the distance and linkage methods.


# EXAMPLE 1: Coral reef invertebrate assemblages ---------------------------

# Import and prepare the community data ----------------------------------
# NOTE: If you are adding this on to an existing script after making an nMDS
# your code may already do this part.

dat_inv_wide <- read_excel(
  "Coral_Reef_Data_2014_to_2025.xlsx",
  sheet = "Invertebrates"
) |>
  # Group identifies the student group that collected the data and is not an
  # ecological sampling variable.
  select(-Group)

dat_inv_long <- dat_inv_wide |>
  pivot_longer(
    cols = -c(Year, Site, Replicate),
    names_to = "Taxa",
    values_to = "Abundance"
  )

# Average replicate transects so each sample in the cluster analysis is one
# site-year, matching the analytical unit used in the nMDS tutorial.
site_year_inv_long <- dat_inv_long |>
  group_by(Site, Year, Taxa) |>
  summarise(
    Mean_abundance = mean(Abundance),
    .groups = "drop"
  ) |>
  mutate(
    Sample_ID = str_c(Site, Year, sep = "_")
  ) |>
  relocate(Sample_ID, Site, Year, Taxa)

# Apply the same 10% prevalence rule used in the nMDS tutorial.
taxon_prevalence <- site_year_inv_long |>
  group_by(Taxa) |>
  summarise(
    prevalence_pct = 100 * mean(Mean_abundance > 0),
    .groups = "drop"
  )

taxa_retained <- taxon_prevalence |>
  filter(prevalence_pct >= 10) |>
  pull(Taxa)

wide_dat_inv <- site_year_inv_long |>
  filter(Taxa %in% taxa_retained) |>
  select(Sample_ID, Site, Year, Taxa, Mean_abundance) |>
  pivot_wider(
    names_from = Taxa,
    values_from = Mean_abundance,
    values_fill = 0
  ) |>
  arrange(Year, Site)

taxa_columns <- wide_dat_inv |>
  select(-Sample_ID, -Site, -Year) |>
  names()

# Create a sample x taxon table and apply the same fourth-root transformation
# used for the nMDS. Clustering and ordination should use the same analytical
# choices when their results will be compared.
comm_inv <- wide_dat_inv |>
  select(Sample_ID, all_of(taxa_columns)) |>
  mutate(
    across(-Sample_ID, \(x) x^(1 / 4))
  ) |>
  column_to_rownames(var = "Sample_ID")

# Remove any all-zero taxa. Empty samples must be investigated and removed
# before calculating Bray-Curtis dissimilarities.
comm_inv <- comm_inv[, colSums(comm_inv) > 0, drop = FALSE]

all_zero_samples_inv <- rownames(comm_inv)[rowSums(comm_inv) == 0]
all_zero_samples_inv

# The supplied dataset has no empty site-years. If your data do, determine why
# before uncommenting and using these lines:
# comm_inv <- comm_inv[!rownames(comm_inv) %in% all_zero_samples_inv, , drop = FALSE]
# wide_dat_inv <- wide_dat_inv |>
#   filter(!Sample_ID %in% all_zero_samples_inv)


# Calculate Bray-Curtis dissimilarities and the cluster analysis ---------

distance_inv <- vegan::vegdist(
  comm_inv,
  method = "bray"
)
# View the triangular distance matrix: pairwise dissimilarities between
# assemblages at every pair of samples.
distance_inv

# Create the hierarchical cluster object using average linkage.
cluster_inv <- hclust(
  distance_inv,
  method = "average"
)

cluster_inv

# Key point:
# - The tree depends on the sampling unit, filtering, transformation, distance,
#   and linkage method. These choices must be reported.

# First look at the dendrogram --------------------------------------------

# Convert the hclust object into line segments and tip labels that ggplot2
# can draw. The fusion height shows how dissimilar samples or groups were when
# they joined; larger joins indicate more distinct groups.
dendro_inv <- ggdendro::dendro_data(
  cluster_inv,
  type = "rectangle"
)

dendro_segments_inv <- ggdendro::segment(dendro_inv)

dendro_labels_inv <- ggdendro::label(dendro_inv) |>
  rename(Sample_ID = label)

# Draw the tree before deciding where it should be cut into groups.
plot_cluster_simple_inv <- ggplot() +
  geom_segment(
    data = dendro_segments_inv,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.5,
    color = "grey25"
  ) +
  geom_text(
    data = dendro_labels_inv,
    aes(x = x, y = y, label = Sample_ID),
    angle = 90,
    hjust = 1.05,
    size = 2.5
  ) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.25, 0.05))) +
  coord_cartesian(clip = "off") +
  labs(
    x = NULL,
    y = "Bray-Curtis fusion height"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(5.5, 5.5, 70, 5.5)
  )

plot_cluster_simple_inv # EXPLORATORY PLOT


# Evaluate candidate numbers of clusters with silhouette width ------------

# Compare the mean silhouette width across candidate numbers of clusters.
# The average-linkage method must match the hclust() analysis above.
# Larger values indicate more compact groups that are better separated from
# their nearest neighboring group.
plot_silhouette_profile_inv <- factoextra::fviz_nbclust(
  distance_inv,
  FUNcluster = factoextra::hcut,
  method = "silhouette",
  k.max = 8,
  hc_method = "average",
  mark_optimal = FALSE,
  print.summary = FALSE
)

plot_silhouette_profile_inv # EXPLORATORY/CONFIRMATION PLOT

# The number of clusters (k) with the highest mean silhouette width is a useful
# starting point, not an automatic answer. A nearby value of k may be preferable
# when its mean silhouette width is only slightly lower and it provides a more
# useful biological description of the data.
#
# Compare nearby solutions. The main pattern is relatively stable if most
# samples remain together and changing k mainly splits one broad group or
# combines two similar groups. A solution is less convincing if small changes in
# k, filtering, transformation, or linkage produce very different biological
# conclusions.
#
# A biologically interpretable solution should produce groups that can be
# described using the original taxa or continuous variables. For example,
# clusters might represent assemblages dominated by different taxa or plants
# with distinct combinations of morphological traits. Avoid choosing a
# solution that creates several one-sample groups, lacks recognizable
# biological differences, or is selected only because it matches an expected
# site or year pattern.
#
# When several values of k perform similarly, prefer the simpler solution that
# addresses the study objective and leaves enough samples in each group for
# meaningful interpretation. Report the silhouette results and explain why the
# selected biological resolution was chosen. Also recognize that a continuous
# ecological gradient may not contain strongly separated clusters at all.

# Select a value of k after considering the silhouette results and biology.
# Change this object to compare another grouping throughout the code below.
n_clusters_inv <- 2

# Place the cut line between the two fusion heights that define the selected k.
# Changing n_clusters_inv above automatically moves this line.
cut_height_inv <- cluster_inv$height |>
  tail(n_clusters_inv) |>
  head(2) |>
  mean()

plot_cluster_cut_inv <- plot_cluster_simple_inv +
  geom_hline(
    yintercept = cut_height_inv,
    linetype = "dashed",
    color = "firebrick",
    linewidth = 0.7
  )

plot_cluster_cut_inv # EXPLORATORY/CONFIRMATION PLOT

# Assign each sample to a cluster at the selected value of k.
cluster_membership_inv <- cutree(
  cluster_inv,
  k = n_clusters_inv
)


# Save and summarize the selected cluster membership ---------------------
# Record which samples are in each cluster and add metadata for interpretation.

cluster_table_inv <- tibble(
  Sample_ID = names(cluster_membership_inv),
  Cluster = factor(unname(cluster_membership_inv))
) |>
  left_join(
    wide_dat_inv |>
      select(Sample_ID, Site, Year),
    by = "Sample_ID"
  ) |>
  arrange(Cluster, Year, Site)

cluster_table_inv

cluster_table_inv |>
  count(Cluster, Site)

cluster_table_inv |>
  count(Cluster, Year)

# Compare clusters with known metadata after clustering. Site and year were
# not used to construct the groups, so any correspondence is descriptive.

# Add sample metadata to the dendrogram labels ---------------------------
# These point mappings can match the colors and shapes used in an nMDS plot.
dendro_labels_metadata_inv <- dendro_labels_inv |>
  left_join(cluster_table_inv, by = "Sample_ID")

# Place the metadata points and rotated sample labels below the tree.
point_height_inv <- -0.03 * max(cluster_inv$height)
label_height_inv <- -0.07 * max(cluster_inv$height)

plot_cluster_metadata_inv <- ggplot() +
  # Draw the dendrogram branches.
  geom_segment(
    data = dendro_segments_inv,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.5,
    color = "grey25"
  ) +
  # Show the selected cut height.
  geom_hline(
    yintercept = cut_height_inv,
    linetype = "dashed",
    color = "firebrick",
    linewidth = 0.7
  ) +
  # Map the same metadata used to distinguish samples in the nMDS.
  geom_point(
    data = dendro_labels_metadata_inv,
    aes(x = x, color = Year, shape = Site),
    y = point_height_inv,
    size = 2.7
  ) +
  # Add the sample IDs directly below their corresponding tips.
  geom_text(
    data = dendro_labels_metadata_inv,
    aes(x = x, label = Sample_ID),
    y = label_height_inv,
    angle = 90,
    hjust = 1,
    size = 2.5
  ) +
  scale_color_viridis_c(option = "C", end = 0.9) +
  scale_x_continuous(breaks = NULL) +
  labs(
    x = NULL,
    y = "Bray-Curtis fusion height",
    color = "Year",
    shape = "Site"
  ) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(5.5, 5.5, 90, 5.5),
    legend.position = "right"
  )

plot_cluster_metadata_inv # POTENTIAL RESULTS PLOT

# OPTIONAL: Heatmaps with sample and variable dendrograms =================

# Heatmaps can help identify which taxa or continuous variables characterize
# the main sample groups. Keep the preprocessing, distance, and linkage choices
# consistent with the corresponding cluster analysis.
#
# Documentation and examples:
# https://yunuuuu.github.io/ggalign-gallery/basics.html

# Community heatmap: Bray-Curtis and average linkage ----------------------

# The heatmap displays fourth-root abundance. Rows are site-years and columns
# are taxa. Both dendrograms use Bray-Curtis and average linkage, although the
# selected number of clusters applies only to the site-year dendrogram.

plot_heatmap_inv <- ggalign::ggheatmap(as.matrix(comm_inv)) +
  # Format the taxon names (x) and site-year IDs (y).
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7)
  ) +
  # Reserve space on the left, then cluster the rows (site-years).
  ggalign::anno_left(size = 0.20) +
  ggalign::align_dendro(
    aes(color = branch),
    distance = \(x) vegan::vegdist(x, method = "bray"),
    method = "average",
    k = n_clusters_inv # color branches using the selected number of clusters
  ) +
  # Reserve space above the heatmap, then independently cluster the columns
  # (taxa) with the same dissimilarity and linkage choices.
  ggalign::anno_top(size = 0.20) +
  ggalign::align_dendro(
    distance = \(x) vegan::vegdist(x, method = "bray"),
    method = "average"
  ) +
  # Return to the heatmap panel before setting the abundance color scale.
  ggalign::quad_active() +
  scale_fill_viridis_c(
    option = "C",
    name = "Fourth-root\nabundance"
  ) &
  # Apply a discrete palette to the colored row-dendrogram branches.
  scale_color_brewer(
    palette = "Dark2",
    guide = "none"
  )

plot_heatmap_inv # EXPLORATORY/INTERPRETIVE PLOT

# EXAMPLE 2: Darlingtonia plant morphology --------------------------------

# Import and prepare the continuous measurements -------------------------

dat_darl <- read_csv(
  "Darlingtonia_GE_Table12.1.csv",
  show_col_types = FALSE
) |>
  rename(
    plant_height_mm = height,
    mouth_diameter_mm = mouth.diam,
    tube_diameter_mm = tube.diam,
    keel_diameter_mm = keel.diam,
    wing1_length_mm = wing1.length,
    wing2_length_mm = wing2.length,
    wing_spread_mm = wingsprea,
    hood_mass_g = hoodmass.g,
    tube_mass_g = tubemass.g,
    wing_mass_g = wingmass.g
  ) |>
  mutate(
    site = factor(site),
    sample_ID = str_c(site, plant, sep = "_")
  ) |>
  relocate(sample_ID, site, plant)

pca_variables <- c(
  "plant_height_mm",
  "mouth_diameter_mm",
  "tube_diameter_mm",
  "keel_diameter_mm",
  "wing1_length_mm",
  "wing2_length_mm",
  "wing_spread_mm",
  "hood_mass_g",
  "tube_mass_g",
  "wing_mass_g"
)

dat_darl_complete <- dat_darl |>
  drop_na(all_of(pca_variables))

comm_darl <- dat_darl_complete |>
  select(sample_ID, all_of(pca_variables)) |>
  column_to_rownames(var = "sample_ID")

# Center and scale every variable before calculating Euclidean distances.
# Distances are therefore measured in standard-deviation units, and variables
# with large numerical values or variance do not dominate the analysis.
comm_darl_scaled <- comm_darl |>
  scale() |>
  as.data.frame()


# Calculate Euclidean distances and Ward clustering ----------------------

distance_darl <- dist(
  comm_darl_scaled,
  method = "euclidean"
)

# Show the triangular distance matrix: pairwise distances among plants based
# on their standardized morphological and biomass measurements.
distance_darl

# Create the hierarchical cluster object using Ward's linkage.
cluster_darl <- hclust(
  distance_darl,
  method = "ward.D2"
)

cluster_darl

# Ward.D2 is appropriate here because the input distances are Euclidean. The
# resulting clusters seek compact groups of plants with similar standardized
# combinations of morphological and biomass measurements.

# First look at the dendrogram --------------------------------------------

# Convert the hclust object into line segments and tip labels for ggplot2.
dendro_darl <- ggdendro::dendro_data(
  cluster_darl,
  type = "rectangle"
)

dendro_segments_darl <- ggdendro::segment(dendro_darl)

dendro_labels_darl <- ggdendro::label(dendro_darl) |>
  rename(sample_ID = label)

# Draw the tree before deciding where it should be cut into groups.
plot_cluster_simple_darl <- ggplot() +
  geom_segment(
    data = dendro_segments_darl,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.5,
    color = "grey25"
  ) +
  geom_text(
    data = dendro_labels_darl,
    aes(x = x, y = y, label = sample_ID),
    angle = 90,
    hjust = 1.05,
    size = 2.5
  ) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.25, 0.05))) +
  coord_cartesian(clip = "off") +
  labs(
    x = NULL,
    y = "Ward fusion height"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(5.5, 5.5, 70, 5.5)
  )

plot_cluster_simple_darl # EXPLORATORY PLOT


# Evaluate candidate numbers of clusters with silhouette width ------------

# Compare the mean silhouette width across candidate numbers of clusters.
# Ward.D2 must match the hclust() analysis above. Use the general selection
# guidance in Example 1 rather than choosing k automatically.
plot_silhouette_profile_darl <- factoextra::fviz_nbclust(
  distance_darl,
  FUNcluster = factoextra::hcut,
  method = "silhouette",
  k.max = 8,
  hc_method = "ward.D2",
  mark_optimal = FALSE,
  print.summary = FALSE
)

plot_silhouette_profile_darl # EXPLORATORY/CONFIRMATION PLOT

# Select k after considering both silhouette widths and the biological meaning
# of the groups. Change this object to compare another solution below.
n_clusters_darl <- 7

# Place the cut line between the fusion heights that define the selected k.
# Changing n_clusters_darl above automatically moves this line.
cut_height_darl <- cluster_darl$height |>
  tail(n_clusters_darl) |>
  head(2) |>
  mean()

plot_cluster_cut_darl <- plot_cluster_simple_darl +
  geom_hline(
    yintercept = cut_height_darl,
    linetype = "dashed",
    color = "firebrick",
    linewidth = 0.7
  )

plot_cluster_cut_darl # EXPLORATORY/CONFIRMATION PLOT

# Assign each plant to a cluster at the selected value of k.
cluster_membership_darl <- cutree(
  cluster_darl,
  k = n_clusters_darl
)


# Save and summarize the selected cluster membership ---------------------
# Record which plants are in each cluster and add site metadata.

cluster_table_darl <- tibble(
  sample_ID = names(cluster_membership_darl),
  Cluster = factor(unname(cluster_membership_darl))
) |>
  left_join(
    dat_darl_complete |>
      select(sample_ID, site, plant),
    by = "sample_ID"
  ) |>
  arrange(Cluster, site, plant)

cluster_table_darl

cluster_table_darl |>
  count(Cluster, site)

# Add site metadata to the dendrogram labels -----------------------------
# Point colors and shapes can match those used in the PCA plot.
dendro_labels_metadata_darl <- dendro_labels_darl |>
  left_join(cluster_table_darl, by = "sample_ID")

# Place the metadata points and rotated sample labels below the tree.
point_height_darl <- -0.03 * max(cluster_darl$height)
label_height_darl <- -0.07 * max(cluster_darl$height)

plot_cluster_metadata_darl <- ggplot() +
  # Draw the dendrogram branches.
  geom_segment(
    data = dendro_segments_darl,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.5,
    color = "grey25"
  ) +
  # Show the selected cut height.
  geom_hline(
    yintercept = cut_height_darl,
    linetype = "dashed",
    color = "firebrick",
    linewidth = 0.7
  ) +
  # Show the site of each plant using both color and shape.
  geom_point(
    data = dendro_labels_metadata_darl,
    aes(x = x, color = site, shape = site),
    y = point_height_darl,
    size = 2.7
  ) +
  # Add the sample IDs directly below their corresponding tips.
  geom_text(
    data = dendro_labels_metadata_darl,
    aes(x = x, label = sample_ID),
    y = label_height_darl,
    angle = 90,
    hjust = 1,
    size = 2.5
  ) +
  scale_x_continuous(breaks = NULL) +
  labs(
    x = NULL,
    y = "Ward fusion height",
    color = "Site",
    shape = "Site"
  ) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(5.5, 5.5, 90, 5.5),
    legend.position = "right"
  )

plot_cluster_metadata_darl # POTENTIAL RESULTS PLOT


# Continuous-variable heatmap: Euclidean and Ward.D2 ---------------------

# The heatmap displays standardized measurements. A value of zero is the
# variable mean, and positive or negative values indicate measurements above
# or below that mean in standard-deviation units.
# 
# Examples: https://yunuuuu.github.io/ggalign-gallery/basics.html

plot_heatmap_darl <- ggalign::ggheatmap(as.matrix(comm_darl_scaled)) +
  # Format the variable names (x) and plant IDs (y).
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7)
  ) +
  # Reserve space on the left, then cluster the rows (plants).
  ggalign::anno_left(size = 0.20) +
  ggalign::align_dendro(
    aes(color = branch),
    distance = "euclidean",
    method = "ward.D2",
    k = n_clusters_darl
  ) +
  # Reserve space above the heatmap, then independently cluster the columns
  # (morphological and biomass variables).
  ggalign::anno_top(size = 0.20) +
  ggalign::align_dendro(
    distance = "euclidean",
    method = "ward.D2"
  ) +
  # Return to the heatmap panel before setting the standardized-value scale.
  ggalign::quad_active() +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Standardized\nvalue"
  ) &
  # Apply a discrete palette to the colored row-dendrogram branches.
  scale_color_brewer(
    palette = "Dark2",
    guide = "none"
  )

plot_heatmap_darl # EXPLORATORY/INTERPRETIVE PLOT


# Reporting guidance -----------------------------------------------------

# A Methods description should report:
# - the analytical unit and any averaging of replicate observations;
# - variables or taxa included and any exclusions or rare-taxon filtering;
# - transformations and standardization;
# - the distance or dissimilarity measure;
# - the hierarchical linkage method used with hclust();
# - the number of clusters chosen and how the final cut was evaluated; and
# - the R packages and functions used.
#
# A Results description should report:
# - the selected number of clusters and how the silhouette pattern supported it;
# - samples assigned to each cluster (summarize and cite plots or tables);
# - associations between clusters and known site, year, or other metadata;
# - taxa or continuous variables that characterize the groups.
#
# Describe the main biological pattern(s) before listing individual memberships.
# Treat the clusters as an exploratory summary unless a separate analysis that
# matches the study design evaluates a specific hypothesis about group structure.

