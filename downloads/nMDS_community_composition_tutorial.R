# (BIO 5150L) nMDS of community composition tutorial
# Coral reef invertebrate assemblages among sites and years
# Date modified: 11 Sep 2026
#
# Guiding question:
# How do invertebrate assemblages vary across sites and years?
#
# Required file in the same folder:
#   Coral_Reef_Data_2014_to_2025.xlsx
#
# This script follows the data-summary and diversity tutorial. See that script
# for more examples of sampling-effort summaries and community plots.


# Load packages -----------------------------------------------------------

library(tidyverse)
library(readxl)
library(vegan)
library(ggrepel)

# Used only for the stacked-bar color palette near the end of the tutorial.
# Install it once if it is not already installed:
# install.packages("Polychrome")


# Import and briefly inspect the data -------------------------------------

dat_inv_wide <- read_excel(
  "Coral_Reef_Data_2014_to_2025.xlsx",
  sheet = "Invertebrates"
) |>
  # Group identifies the student group that collected the data. It is not an
  # ecological sampling variable and will not be used in this analysis.
  select(-Group)

glimpse(dat_inv_wide)
head(dat_inv_wide)

# Key points to check:
# - One row represents one transect.
# - Year, Site, and Replicate identify the sampling unit.
# - All remaining columns contain taxon abundances.
# - Missing abundances must be resolved before the analysis.

data_quality_check <- dat_inv_wide |>
  pivot_longer(
    cols = -c(Year, Site, Replicate),
    names_to = "Taxa",
    values_to = "Abundance"
  ) |>
  summarise(
    n_transects = n_distinct(Year, Site, Replicate),
    n_taxa_columns = n_distinct(Taxa),
    n_missing = sum(is.na(Abundance))
  )

data_quality_check

# If missing abundance values are present (i.e., there are NAs), they need to 
# be resolved before moving forward. Need to figure out what an NA means, and 
# do not automatically replace missing observations with zero unless zero 
# truly means that a taxon was searched for and not observed.


# Convert to long format and summarize replicate transects ----------------

dat_inv_long <- dat_inv_wide |>
  pivot_longer(
    cols = -c(Year, Site, Replicate),
    names_to = "Taxa",
    values_to = "Abundance"
  )



# Summarize sampling effort by site and year ------------------------------
sampling_effort_table <- dat_inv_long |>
  distinct(Site, Year, Replicate) |>
  count(
    Site,
    Year,
    name = "n_transects"
  ) |>
  pivot_wider(
    names_from = Year,
    values_from = n_transects,
    values_fill = 0,
    names_sort = TRUE
  ) |>
  arrange(Site)

sampling_effort_table


## Main analytical choice: Sampling unit (i.e., what will the points be in your
# nMDS plot)?
# 
# Before running an nMDS ordination, decide what each point will represent. A point
# could represent an individual transect or plot, or replicate samples could
# first be averaged to a higher level, such as a site, site-year, habitat, or
# treatment.

# This decision should match the research question and study design. Retaining
# individual replicates shows within-group variation, while averaging focuses the
# analysis on differences among the higher-level sampling units. Averaging
# replicates also prevents groups with greater sampling effort from receiving
# more weight simply because they contain more samples. However, it reduces the
# number of observations and removes information about variation among
# replicates. It also depends on if you think a single transect, plot, or
# experimental unit, does a good enough job of representing the assemblage, and
# if not you are often best off pooling or averaging across replicates to get a
# more robust estimate of the assemblage at that site or treatment.
#
# (This tutorial) Each nMDS point in this tutorial will represent one site-year.
#  Average the replicate transects within each Site x Year x Taxa combination 
#  first. This gives every transect sampled within a site-year equal weight.

# BUT NOTE - typically when working with your own data you probably will try
# running multiple nMDS ordinations with different sampling units to see how the
# results related to each other, but eventually you will make a decision on what
# is most applicable for your research question and study design and go with that.

# average across transects to each site + year -----------------------------

site_year_inv_long <- dat_inv_long |>
  group_by(Site, Year, Taxa) |>
  summarise(
    Mean_abundance = mean(Abundance),
    n_transects = n_distinct(Replicate),
    .groups = "drop"
  ) |>
  mutate(
    Sample_ID = str_c(Site, Year, sep = "_")
  ) |>
  relocate(Sample_ID, Site, Year, Taxa)

site_year_inv_long

# Confirm sampling effort before interpreting site-year means.
site_year_effort <- site_year_inv_long |>
  distinct(Sample_ID, Site, Year, n_transects) |>
  arrange(Year, Site)

site_year_effort

# Key points to interpret:
# - The site-year is the sampling unit used in this nMDS.
# - A mean describes the typical transect rather than the total number counted.
# - Site-year means prevent years with more transects from contributing larger
#   totals simply because sampling effort was greater.
# - Unequal precision can remain because means based on fewer transects are
#   less precisely estimated. Record and report differences in effort.


# Explore transformation options -----------------------------------------

# Community data commonly contain a few very abundant taxa and many less
# abundant taxa. A transformation controls how strongly dominant taxa affect
# Bray-Curtis dissimilarities and the resulting nMDS.
#
# Possible options (in most cases pick between square root and fourth root):
# - No transformation: preserves abundance differences most strongly, but is
# rarely used for species-abundance data with strongly dominant taxa.
# - Square root: reduces dominance while retaining substantial information
# about relative abundance. This is often preferred when differences in
# abundance and which taxa dominate are central to the question.
# - Fourth root: downweights dominant taxa more strongly and gives uncommon
# taxa more influence. This is often preferred when overall biodiversity and the
# contributions of less abundant taxa are important.
# - log(x + 1): also strongly downweights high values. For these data it is
# generally very similar to the fourth root, so the fourth root is simpler to
# explain and is preferred here.
# - Presence/absence: the most extreme transformation. It discards abundance
# information and compares only which taxa were present. Jaccard dissimilarity
# is commonly used for presence/absence community data.
#
# The choice must be made from the study objectives and documented. Do not
# choose a transformation only because it creates the most visually distinct
# groups.

# Compare a high-abundance taxon with a lower-abundance taxon.
focal_taxa <- c(
  "Pale Rock-boring Urchin",
  "Sea Cucumber"
)

transformation_comparison <- site_year_inv_long |>
  filter(Taxa %in% focal_taxa) |>
  transmute(
    Taxa,
    Raw = Mean_abundance,
    `Square root` = sqrt(Mean_abundance),
    `Fourth root` = Mean_abundance^(1 / 4),
    `log(x + 1)` = log1p(Mean_abundance)
  ) |>
  pivot_longer(
    cols = -Taxa,
    names_to = "Transformation",
    values_to = "Transformed_abundance"
  ) |>
  mutate(
    Transformation = factor(
      Transformation,
      levels = c("Raw", "Square root", "Fourth root", "log(x + 1)")
    )
  )

plot_transformation_comparison <- transformation_comparison |>
  ggplot(aes(x = Transformed_abundance)) +
  geom_histogram(bins = 16, color = "white", fill = "steelblue") +
  facet_grid(
    rows = vars(Taxa),
    cols = vars(Transformation),
    scales = "free_x"
  ) +
  labs(
    x = "Abundance on the displayed scale",
    y = "Number of site-years"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.text.x = element_text(size = 9)
  )

plot_transformation_comparison # EXPLORATORY/CONFIRMATION PLOT

# This histogram comparison documents what each transformation does. It is
# useful during analysis to think about how the transformation you pick will
# affect the influence of abundant versus uncommon species, but it would not
# normally be included as a thesis or paper figure.
# Look for how the high values are compressed and how similar they become
# between an abundant versus less abundant species and what makes more sense
# for your analysis.
# Then pick square root if you want to emphasize the dominant taxa, or fourth
# root if you want to emphasize overall biodiversity and the contributions of
# less abundant taxa.

# For the main analysis below, we'll use square-root transformation, to let
# abundant species have a greater impact on the nMDS because the most abundant
# urchin species (rock-boring urchins) cause issues with bioerosion and can
# break down the coral reefs so their high abundance has ecological consequences
# we want to highlight. However, if you were more interested in overall
# biodiversity and the contributions of less abundant taxa, you would choose
# fourth root transformation instead.


# Optionally filter uncommon taxa using prevalence ------------------------

# Very uncommon taxa can add many joint absences and can be strongly affected
# by chance detections or inconsistent identification. Filtering them may make
# broad assemblage patterns more stable and interpretable. However, uncommon
# taxa may be biologically important, so filtering depends on the study goal.
# Always calculate prevalence on the untransformed data and report the rule.
# A filter threshold you commonly see used is occurring in at least 5% of 
# samples (plots, transects, etc). But this could be adjusted up or down 
# depending on the nature of your data set and research question. For example,
# if you have a very large data set with many samples, you may want to use a
# higher threshold (e.g., 10% or 20%) to focus on the more common taxa.
# Conversely, if your data set is small or you are interested in rare taxa, you
# may choose a lower threshold (e.g., 1% or 2%).
# 
# Here, prevalence is the percentage of transects in which a taxon occurred.

taxon_prevalence <- dat_inv_long |>
  group_by(Taxa) |>
  summarise(
    n_transects = n(),
    n_transects_present = sum(Abundance > 0),
    prevalence_pct = 100 * n_transects_present / n_transects,
    .groups = "drop"
  ) |>
  arrange(desc(prevalence_pct))

taxon_prevalence

# Retain taxa observed in at least 5% of transects.
rare_taxon_cutoff <- 5

taxa_retained <- taxon_prevalence |>
  filter(prevalence_pct >= rare_taxon_cutoff) |>
  pull(Taxa)

taxa_removed <- taxon_prevalence |>
  filter(prevalence_pct < rare_taxon_cutoff)

# Always check and see what is being removed vs. retained
taxa_retained
taxa_removed

# filter for retained taxa
site_year_inv_filtered <- site_year_inv_long |>
  filter(Taxa %in% taxa_retained)

# Key points to interpret:
# - Other defensible rules include a minimum number of samples, a minimum total
#   abundance, or no filtering when rare taxa are the focus of the study.
# - Try a small number of defensible thresholds as a sensitivity check. Major
#   conclusions should not depend entirely on one arbitrary cutoff.
# - Filtering and transformation answer different problems. Filtering removes
#   taxa; transformation changes the influence of the taxa that remain.


# Construct the sample x taxa tables --------------------------------------

# Keep a wide table with sample metadata and untransformed mean abundances.
# We will join this complete table to the nMDS scores so any retained taxon can
# later be mapped to point size or used to help interpret the ordination.

wide_dat_inv <- site_year_inv_filtered |>
  select(Sample_ID, Site, Year, Taxa, Mean_abundance) |>
  pivot_wider(
    names_from = Taxa,
    values_from = Mean_abundance,
    values_fill = 0
  ) |>
  relocate(Sample_ID, Site, Year) |>
  arrange(Year, Site)

glimpse(wide_dat_inv)

taxa_columns <- setdiff(
  names(wide_dat_inv),
  c("Sample_ID", "Site", "Year")
)

taxa_columns

# Create the community table required by vegan. Rows are site-years, columns
# are taxa, and cells are transformed mean abundances.

comm_inv <- wide_dat_inv |>
  select(Sample_ID, all_of(taxa_columns)) |>
  # transform abundances in the community table
  mutate(
    across(-Sample_ID, ~ .x^(1 / 2)) 
    #x^(1/2) for square root transformation; x^(1/4) for fourth root transformation 
  ) |>
  column_to_rownames(var = "Sample_ID")

dim(comm_inv)
View(comm_inv)
# note in comm_inv, rows are site-years (Sample_ID) and columns are taxa, and
# the values are the transformed mean abundances. This is the format required
# for vegan functions like metaMDS() and adonis2().

# Key points to interpret:
# - comm_inv is used to calculate Bray-Curtis dissimilarities and the nMDS below.
# - wide_dat_inv retains metadata and raw site-year mean abundances.
# - The transformation is applied once here. metaMDS() will be told not to
#   apply another automatic transformation.

# Check for all-zero taxa and samples -------------------------------------

# Identify if any taxa have zero total abundance across all samples
all_zero_taxa <- comm_inv |>
  # sum abundance for each taxon across all samples
  summarise(
    across(
      where(is.numeric),
      ~ sum(.x, na.rm = TRUE)
    )
  ) |>
  pivot_longer(
    everything(),
    names_to = "Taxa",
    values_to = "Total_abundance"
  ) |>
  filter(Total_abundance == 0) |>
  pull(Taxa)

all_zero_taxa
# if character(0) then there are no taxa with 0 abundance in the data table

# If needed, remove all-zero taxa
# comm_inv <- comm_inv |>
#   select(-any_of(all_zero_taxa))


# Identify samples with zero total abundance across all taxa
all_zero_samples <- comm_inv |>
  rownames_to_column(var = "Sample_ID") |>
  mutate(
    Total_abundance = rowSums(
      pick(where(is.numeric)),
      na.rm = TRUE
    )
  ) |>
  filter(Total_abundance == 0) |>
  pull(Sample_ID)

all_zero_samples
# if character(0) then there are no samples with 0 abundance for all taxa in the data table

# If needed, remove all-zero samples from both tables
# comm_inv <- comm_inv |>
#   rownames_to_column(var = "Sample_ID") |>
#   filter(!Sample_ID %in% all_zero_samples) |>
#   column_to_rownames(var = "Sample_ID")

# wide_dat_inv <- wide_dat_inv |>
#   filter(!Sample_ID %in% all_zero_samples)

# Run a two-dimensional nMDS ----------------------------------------------

# Bray-Curtis is widely used for species-abundance data. It compares the
# composition and abundance of taxa and is not inflated by joint absences.
#
# nMDS begins from multiple random configurations and may find local optima.
# set.seed() makes the random starts and final result reproducible.

set.seed(38)

nmds_inv <- metaMDS(
  comm_inv,
  distance = "bray",
  k = 2,
  try = 50,
  trymax = 200,
  autotransform = FALSE,
  trace = 1
)

nmds_inv

# stress metric from the nMDS
nmds_inv$stress

# Key points to interpret:
# - The console output compares repeated random starts. A repeated best solution
# provides evidence that a stable configuration was found. If the best solution
# was not repeated, increase trymax and rerun with the same seed.
# - Stress measures mismatch between the ranked Bray-Curtis dissimilarities and
# distances on the ordination plot (i.e., how hard was it to get the
# multidimensional space to be represented in 2 dimensions). Lower stress is
# better.
# - Traditional guidance treated stress > 0.20 as potentially poor. Current
# ecological guidance is often more flexible: 0.20 to 0.30 may be usable for
# complex data if interpreted cautiously, while > 0.30 is generally poor.
# - Stress is guidance, not a universal pass/fail threshold. Consider sample
# size, complexity, stability across starts, and whether conclusions are robust
# to reasonable analytical choices. 
# Very low stress values (near zero) should also be checked because they can
# result from too few samples, duplicate samples, or limited variation (e.g.,
# all/most samples are very similar).

# Display the Shepard stress plot -----------------------------------------

stressplot( # EXPLORATORY/CONFIRMATION PLOT
  nmds_inv,
  p.col = "gray40",
  l.col = "firebrick",
  lwd = 2
)

# Each point compares an original Bray-Curtis dissimilarity with its fitted
# distance in the ordination. Points close to the monotonic fitted line are
# represented well. Large or systematic departures indicate distortion.
# This diagnostic is worth checking, but it usually receives less emphasis
# than stress, convergence, and biological interpretation. Just a diagnostic
# plot, not a primary result. This would not be something you'd include in a
# thesis.


# Create table of nMDS result to use with ggplot --------------------------
# Extract scores and add metadata and taxon abundances

# In vegan, "site scores" are the x- and y-axis coordinates of the 
# samples (points) in your ordination plot, even when the sampling units are transects, plots, or other units.

nmds_scores_inv <- scores(
  nmds_inv,
  display = "sites"
) |>
  as.data.frame() |>
  rownames_to_column(var = "Sample_ID") |>
  as_tibble() |>
  left_join(
    wide_dat_inv,
    by = "Sample_ID"
  )

glimpse(nmds_scores_inv)

# Key points to interpret:
# - Each row represents one sample (point) included in the nMDS.
#
# - NMDS1 and NMDS2 are the coordinates used to position each sample in the
#   two-dimensional ordination (ggplot).
#   
# - Sample_ID was stored as the row name of the community matrix. Converting it
#   back to a column allows the ordination scores to be joined to the original
#   sample information.
#
# - Adding metadata columns, such as Site and Year, allows these variables to be
#   mapped to point color, shape, labels, or facets when plotting the nMDS.
#
# - Retaining the taxon-abundance columns also allows the abundance of selected
#   taxa to be added to the plot or used to help interpret assemblage patterns.

# Format used for the nMDS result plots ------------------------------

# nMDS plot axes do not have independent biological meanings (so we remove axes
# ticks and numbering). The configuration can be rotated, reflected, or shifted
# without changing the result. Interpret the relative distances among points:
# nearby points have more similar assemblages, and distant points have more
# different assemblages.

nmds_plot_theme <- theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "right"
  )

# extract and format stress metric to add to the plot
stress_label <- str_c(
  "2D stress = ",
  format(round(nmds_inv$stress, 2), nsmall = 2) #first round, then show 2 decimal places, This ensures a stress value such as 0.1 is displayed as "Stress = 0.10"
)
stress_label

# Main nMDS plot: year gradient and site shapes ---------------------------

plot_nmds_year_site <- nmds_scores_inv |>
  ggplot(
    aes(
      x = NMDS1,
      y = NMDS2,
      color = Year,
      shape = Site
    )
  ) +
  geom_point(size = 3.5, alpha = 0.9) +
  scale_color_viridis_c(
    option = "C",
    end = 0.9,
    breaks = sort(unique(nmds_scores_inv$Year))
  ) +
  scale_shape_manual(values = c(16, 17, 15, 18, 8)) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = stress_label,
    hjust = 1.1,
    vjust = -0.7,
    size = 3.5
  ) +
  coord_equal() +
  labs(color = "Year", shape = "Site") +
  nmds_plot_theme

plot_nmds_year_site

# PRIME RESULTS PLOT
# This is a strong starting point for a thesis or paper figure because it shows
# both spatial and temporal structure without assigning meaning to the axes.
# Look for clustering by point shape and directional change along the color
# gradient. Overlap indicates similar assemblages, not identical samples.


# Option: add label each point with its year ----------------------------------

plot_nmds_year_labels <- plot_nmds_year_site +
  geom_text_repel(
    aes(label = Year),
    size = 3,
    show.legend = FALSE
  )

plot_nmds_year_labels

# EXPLORATORY OR ALTERNATE RESULTS PLOT
# Labels make individual years explicit but can become crowded. They are useful
# for checking individual points, often just used in the exploratory stage, then
# removed for the final plot (or sometimes can be good on a final plot too).


# Alternative mapping: emphasize sites and label every point -------------

plot_nmds_site_labels <- nmds_scores_inv |>
  ggplot(aes(x = NMDS1, y = NMDS2, color = Site)) +
  geom_point(size = 3.2, alpha = 0.85) +
  geom_text_repel(
    aes(label = Site),
    size = 3,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  scale_color_brewer(palette = "Dark2") +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = stress_label,
    hjust = 1.1,
    vjust = -0.7,
    size = 3.5
  ) +
  coord_equal() +
  labs(color = "Site") +
  nmds_plot_theme

plot_nmds_site_labels

# ALTERNATE RESULTS PLOT
# This version makes site membership immediately visible. Repeated site labels
# can be redundant, so choose this format only when it improves interpretation.


# Add site centroids ------------------------------------------------------

site_centroids <- nmds_scores_inv |>
  group_by(Site) |>
  summarise(
    centroid_NMDS1 = mean(NMDS1),
    centroid_NMDS2 = mean(NMDS2),
    .groups = "drop"
  )

nmds_scores_centroids <- nmds_scores_inv |>
  left_join(site_centroids, by = "Site")

plot_nmds_centroids <- nmds_scores_centroids |>
  ggplot(aes(x = NMDS1, y = NMDS2, color = Site)) +
  geom_segment(
    aes(
      xend = centroid_NMDS1,
      yend = centroid_NMDS2
    ),
    linewidth = 0.4,
    alpha = 0.45,
    show.legend = FALSE
  ) +
  geom_point(size = 2.8, alpha = 0.85) +
  geom_point(
    data = site_centroids,
    aes(x = centroid_NMDS1, y = centroid_NMDS2, fill = Site),
    shape = 23,
    color = "black",
    size = 5,
    inherit.aes = FALSE
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2", guide = "none") +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = stress_label,
    hjust = 1.1,
    vjust = -0.7,
    size = 3.5
  ) +
  coord_equal() +
  labs(color = "Site") +
  nmds_plot_theme

plot_nmds_centroids

# EXPLORATORY OR SUPPORTING RESULTS PLOT
# Centroids summarize the average ordination location for each site. Segments
# show within-site variation through time. Centroids are descriptive and are
# not a statistical test of group differences.


# Add minimum convex polygons around years --------------------------------

# Identify the outermost points for each year.
# A polygon requires at least three samples within a year.

year_hulls <- nmds_scores_inv |>
  group_by(Year) |>
  filter(n() >= 3) |>
  slice(
    chull(NMDS1, NMDS2)
  ) |>
  ungroup()


plot_nmds_year_polygons <- plot_nmds_year_site +
  geom_polygon(
    data = year_hulls,
    aes(
      x = NMDS1,
      y = NMDS2,
      group = Year,
      fill = Year
    ),
    inherit.aes = FALSE,
    alpha = 0.15,
    color = NA,
    show.legend = FALSE
  ) +
  # Redraw points so they remain visible over the polygon fills
  geom_point(
    data = nmds_scores_inv,
    aes(
      x = NMDS1,
      y = NMDS2,
      color = Year,
      shape = Site
    ),
    size = 3.5,
    alpha = 0.9
  ) +
  scale_fill_viridis_c(
    option = "C",
    end = 0.9
  )

plot_nmds_year_polygons

# Add site ellipses -------------------------------------------------------

# possibly useful if you have lots of points per group to summarize
# not suggested if you have a more limited number of points per group in your nMDS plot (maybe less than 10 or even 20 points per group).

plot_nmds_ellipses <- nmds_scores_inv |>
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  # Draw transparent fills first.
  stat_ellipse(
    aes(group = Site, fill = Site),
    geom = "polygon",
    type = "t",
    level = 0.95,
    alpha = 0.22,
    color = NA
  ) +
  # White outlines on top make boundaries visible where ellipses overlap.
  stat_ellipse(
    aes(group = Site),
    type = "t",
    level = 0.95,
    color = "white",
    linewidth = 1
  ) +
  geom_point(
    aes(fill = Site),
    shape = 21,
    color = "black",
    size = 3,
    alpha = 0.9
  ) +
  scale_fill_brewer(palette = "Dark2") +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = stress_label,
    hjust = 1.1,
    vjust = -0.7,
    size = 3.5
  ) +
  coord_equal() +
  labs(fill = "Site") +
  nmds_plot_theme

plot_nmds_ellipses

# PRIME OR SUPPORTING RESULTS PLOT
# Ellipses summarize the location and spread of each site's points. They require
# enough observations per group and should not be interpreted as a hypothesis
# test. Large overlap suggests weak separation; different ellipse sizes may
# indicate different within-site dispersion, which also matters for PERMANOVA.
# FOR THIS EXAMPLE PLOT, THERE ARE PROBABLY NOT ENOUGH POINTS TO USE ELLIPSES,
# BUT IT IS INCLUDED HERE AS AN EXAMPLE OF HOW TO ADD THEM.


# Scale point size to the abundance of one taxon --------------------------

selected_taxon <- "Blue-black Urchin"

plot_nmds_selected_taxon <- nmds_scores_inv |>
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(
    aes(
      size = .data[[selected_taxon]],
      color = Year
    ),
    alpha = 0.8
  ) +
  scale_size_area(max_size = 9) +
  scale_color_viridis_c(option = "C", end = 0.9) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = stress_label,
    hjust = 1.1,
    vjust = -0.7,
    size = 3.5
  ) +
  coord_equal() +
  labs(
    size = str_c("Mean abundance\n", selected_taxon),
    color = "Year"
  ) +
  nmds_plot_theme

plot_nmds_selected_taxon

# EXPLORATORY OR SUPPORTING RESULTS PLOT
# Larger points show greater untransformed mean abundance of the selected taxon.
# This helps connect assemblage positions to a focal species.
# One way to use this, is to make a panel of 4-6 nMDS plots, one for each taxa
# that are ecologically important or abundant, and then compare the patterns
# across those taxa to help interpret the nMDS plot.


# Fit species-abundance vectors with envfit() -----------------------------

# envfit() finds the direction of the strongest linear increase for each taxon
# across the ordination and evaluates that association using permutations.
# These vectors describe correlations with the configuration. They do not show
# causation, and the p-values should not be treated as independent confirmatory
# tests for every species.

set.seed(38)

nmds_species_fit <- envfit(
  nmds_inv,
  comm_inv,
  permutations = 9999
)

nmds_species_fit

envfit_vectors <- scores(
  nmds_species_fit,
  display = "vectors"
) |>
  as.data.frame() |>
  rownames_to_column(var = "Taxa") |>
  as_tibble() |>
  mutate(
    r2 = unname(nmds_species_fit$vectors$r[Taxa]),
    p_value = unname(nmds_species_fit$vectors$pvals[Taxa])
  ) |>
  arrange(p_value)

envfit_vectors

# Retain vectors with p <= 0.05 for the figure. If too many overlap, use a more
# restrictive prespecified cutoff such as p <= 0.01 and/or an r2 threshold.
envfit_vectors_plot <- envfit_vectors |>
  filter(p_value <= 0.05)

# # Can scale arrows to better fit the plotting region (but keep relative
# # proportions of # arrows to each other). This multiplier affects only the
# # display.
# arrow_multiplier <- 0.8 * min(
#   diff(range(nmds_scores_inv$NMDS1)),
#   diff(range(nmds_scores_inv$NMDS2))
# ) / max(
#   sqrt(envfit_vectors_plot$NMDS1^2 + envfit_vectors_plot$NMDS2^2)
# )
# 
# envfit_vectors_plot <- envfit_vectors_plot |>
#   mutate(
#     NMDS1 = NMDS1 * arrow_multiplier,
#     NMDS2 = NMDS2 * arrow_multiplier
#   )

plot_nmds_envfit <- plot_nmds_year_site +
  geom_segment(
    data = envfit_vectors_plot,
    aes(
      x = 0,
      y = 0,
      xend = NMDS1,
      yend = NMDS2
    ),
    inherit.aes = FALSE,
    arrow = grid::arrow(length = grid::unit(0.18, "cm")),
    color = "black",
    linewidth = 0.6
  ) +
  geom_text_repel(
    data = envfit_vectors_plot,
    aes(
      x = NMDS1,
      y = NMDS2,
      label = Taxa
    ),
    inherit.aes = FALSE,
    size = 3,
    color = "black",
    max.overlaps = Inf
  )

plot_nmds_envfit

# PRIME OR SUPPORTING RESULTS PLOT
# Arrow direction shows increasing transformed abundance. Longer arrows have
# stronger associations with the displayed configuration, but arrow lengths
# are only comparable within this plot. Projection of a site-year point onto an
# arrow approximates its position along that taxon's abundance gradient.


# Plot species scores -----------------------------------------------------

species_scores_inv <- scores(
  nmds_inv,
  display = "species",
  shrink = TRUE
) |>
  as.data.frame() |>
  rownames_to_column(var = "Taxa") |>
  as_tibble()


# Add species names at their exact nMDS locations
plot_nmds_species_scores <- plot_nmds_year_site +
  geom_text(
    data = species_scores_inv,
    aes(
      x = NMDS1,
      y = NMDS2,
      label = Taxa
    ),
    inherit.aes = FALSE,
    size = 3,
    fontface = "bold",
    color = "black"
  )

plot_nmds_species_scores

# EXPLORATORY OR SUPPORTING RESULTS PLOT
# Species scores are weighted-average locations. A taxon is placed near the
# site-years where it tends to be more abundant. In contrast, envfit arrows show
# a direction of increasing abundance plus r2 and a permutation p-value. Species
# scores are locations, not directional tests. Avoid displaying so many labels
# that the community pattern becomes unreadable.


# Supporting stacked bar plot for biological interpretation --------------

# Plot the same site-year means and retained taxa used for the nMDS. The bars
# show abundance, not proportions, so both dominant taxa and changes in total
# mean abundance remain visible.

taxa_present <- site_year_inv_filtered |>
  group_by(Taxa) |>
  summarise(total_mean_abundance = sum(Mean_abundance), .groups = "drop") |>
  arrange(desc(total_mean_abundance)) |>
  pull(Taxa)

# Generate more candidate colors than needed, remove white, and assign a stable
# named color to each taxon. Naming prevents colors from changing across plots.
taxa_color_candidates <- Polychrome::glasbey.colors(
  length(taxa_present) + 3
)

taxa_color_candidates <- taxa_color_candidates[
  toupper(taxa_color_candidates) != "#FFFFFF"
]

taxa_colors <- taxa_color_candidates[seq_along(taxa_present)] |>
  setNames(taxa_present)

plot_stacked_abundance <- site_year_inv_filtered |>
  mutate(
    Taxa = factor(Taxa, levels = rev(taxa_present)),
    Year = factor(Year)
  ) |>
  ggplot(
    aes(
      x = Year,
      y = Mean_abundance,
      fill = Taxa
    )
  ) +
  geom_col(color = "gray25", linewidth = 0.15) +
  facet_wrap(vars(Site), ncol = 2) +
  scale_fill_manual(
    values = taxa_colors,
    breaks = taxa_present,
    drop = TRUE
  ) +
  labs(
    x = "Year",
    y = "Mean abundance per transect",
    fill = "Taxon"
  ) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.key.height = grid::unit(0.35, "cm"),
    legend.key.width = grid::unit(0.35, "cm")
  ) +
  guides(
    fill = guide_legend(ncol = 1, byrow = TRUE)
  )

plot_stacked_abundance

# SUPPORTING OR PRIME RESULTS PLOT
# Use this plot to identify which taxa contribute to site and year patterns in
# the nMDS. A stacked bar plot may be a main result when assemblage composition
# is central to the study. See the data-summary and diversity tutorial for
# additional abundance, proportional-composition, and effort plots.


# PERMANOVA with adonis2() ------------------------------------------------

# nMDS visualizes differences in community composition. PERMANOVA tests
# whether community composition differs among groups or along continuous
# gradients using the same transformed community data and Bray-Curtis
# dissimilarities as the nMDS.

# Arrange metadata in the same order as the community table (it probably is in
# the same order already, but best to make sure).

permanova_metadata <- wide_dat_inv |>
  slice(match(rownames(comm_inv), Sample_ID)) |>
  transmute(
    Sample_ID,
    Site = factor(Site),
    Year = as.numeric(Year)
  )




# Test for differences among sites ----------------------------------------

# IMPORTANT: Because the same sites were sampled repeatedly, the unrestricted
# permutation p-values below are simplified teaching examples. The R2 values
# describe effect magnitude, but a formal hypothesis test would require a
# permutation structure appropriate for the repeated sampling design.

# Site is categorical, so this tests whether assemblage composition differs
# among the sampled sites.

set.seed(38)

permanova_site <- adonis2(
  comm_inv ~ Site,
  data = permanova_metadata,
  permutations = 9999,
  method = "bray",
  by = "margin"
)

permanova_site


# Test for a linear trend through time ------------------------------------

# Year is numerical, so this tests whether assemblage composition changes
# along a linear temporal trend. This does not test for any possible
# difference among individual years.

set.seed(38)

permanova_year <- adonis2(
  comm_inv ~ Year,
  data = permanova_metadata,
  permutations = 9999,
  method = "bray",
  by = "margin"
)

permanova_year

## WARNING:
# Before implementing PERMANOVA with your own data, review that section of this
# course:
# https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/permanova/ 
# And carefully consider your sampling design, model structure, permutation
# restrictions, and how to interpret the results in more detail.


# Key points and considerations -------------------------------------------

# - R2 is the proportion of total multivariate variation associated with the
#   predictor. Consider its magnitude and biological relevance in addition to
#   the permutation p-value.

# - The permutation p-value evaluates whether the observed association is
#   stronger than expected when sample labels are permuted.

# - by = "margin" tests each predictor after accounting for all other terms in
#   the model. With only one predictor, it gives the same result as the overall
#   test. It is included here so the code can be adapted to models containing
#   multiple predictors. by = "margin" is particularly important if you have 2
#   or more predictors that are correlated or unbalanced, because it tests the
#   marginal effect of each predictor while controlling for the others.

# - Marginal tests are especially useful when predictors share explanatory
#   variation, such as in unbalanced designs, datasets with missing sampling
#   combinations, or models containing correlated predictors.

# - In a balanced design with orthogonal predictors, adding another predictor
#   may not change a term's R2. The combined model can still explain more total
#   variation and have less residual variation.

# - A numerical Year predictor tests a linear temporal trend and uses one degree
#   of freedom. Using factor(Year) instead would test for any difference among
#   years without assuming a linear pattern.

# - PERMANOVA assumes that samples are exchangeable under the null hypothesis.
#   Unrestricted permutations are not appropriate when observations are nested,
#   blocked, spatially paired, or repeatedly sampled.

# - Averaging transects means that each observation in this analysis is a
#   site-year mean. Variation among transects within a site-year is no longer
#   represented in the community table or included in the PERMANOVA residual
#   variation.

# - PERMANOVA can respond to differences in group centroids, differences in
#   within-group dispersion, or both. Also use betadisper() when differences in
#   multivariate dispersion are a concern.

# - More complicated models require enough independent sampling units and
#   residual degrees of freedom for all predictors and interactions.


# Methods Section - What to write? ------------------------------------------

# A Methods description should identify the sampling unit and averaging
# procedure; taxon consolidation and filtering; transformation (square or 4th
# root); Bray-Curtis dissimilarity; final stress (typically on the plot); (if
# applicable) species-vector display criteria; PERMANOVA terms and number of
# permutations;
#
# Here’s example methods text from my grad student's paper Calderon et al. 2024
# https://www.sciencedirect.com/science/article/pii/S0022098124000637
#
# To quantify the overall diet composition of Garibaldi, the weight for each
# diet content category was averaged across specimens within a given location,
# reef type, and sex. A similarity matrix was constructed with the Bray-Curtis
# similarity coefficient using square-root transformed weight values. These
# relationships were visualized with a two-dimensional non-metric
# multidimensional scaling plot (nMDS) using the metaMDS function in the ‘vegan’
# package (Oksanen et al., 2022). Ordination scores were then visualized as
# vectors on the nMDS for individual stomach content categories that
# significantly (p < 0.05) correlated with the diet composition patterns in the
# nMDS, with arrow length indicating the relative strength of each correlation.
# The effect of reef type, temperature, and sex on the observed variation in
# diet composition was assessed using a permutational multivariate analysis of
# variance (PERMANOVA) with the adonis2 function in the ‘vegan’ package (Oksanen
# et al., 2022). Because our dataset was unbalanced (i.e., artificial reef types
# were only present at some locations), the by = “margin” argument was included
# in the function to analyze the marginal effects of each variable in the model
# after accounting for all other variables.
