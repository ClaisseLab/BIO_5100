#(BIO 5150L) Multivariate ecological data Intro Tutorial: 
# data summaries, exploration, visualization, and diversity metrics
# Date modified: 9 Sep 2026
# 
# Required files in the same folder:
#   Coral_Reef_Data_2014_to_2025.xlsx
#   invertebrate_taxonomy_lookup.csv

# note - a lot of what we do initially with these data 
# (e.g., summarize sampling effort, and species abundance)
# is generally useful for any ecological data set and something 
# you'd want to do with your own data at the start of any analysis.


# Load packages -----------------------------------------------------------

library(tidyverse)
library(readxl)

library(vegan) #multivariate analysis function package

# will need to install Polychrome package for color palettes if you don't have it already
# install.packages("Polychrome")


# Import and examine community data ---------------------------------------

# Import the invertebrate worksheet
dat_inv_wide <- read_excel(
  "Coral_Reef_Data_2014_to_2025.xlsx",
  sheet = "Invertebrates"
)

# Examine the structure and first several rows
glimpse(dat_inv_wide)
head(dat_inv_wide)


# Questions ----------------------------------------------------------------
# 1. What does one row of dat_inv_wide represent?
#
# 2. Which columns together describe the sampling unit (transect in this case)?
#
# 3. Which columns contain invertebrate abundance data?


# Convert from wide to long format ----------------------------------------

# The vegan package provides functions for analyzing ecological communities,
# including diversity metrics, ordinations, and other multivariate analyses.
#
# Many vegan functions require a wide community matrix in which rows are
# samples and columns are taxa. However, the raw wide table above also
# contains sample metadata and is not yet formatted appropriately for vegan.
# We will construct the required community matrix later in the tutorial.
#
# For now, convert the data to long format because it is easier to summarize,
# modify, and visualize using dplyr and ggplot2.

dat_inv_long <- dat_inv_wide |>
  pivot_longer(
    cols = -c(Year, Site, Replicate, Group),
    names_to = "Taxa",
    values_to = "Abundance"
  )

# Examine the structure and first several rows
glimpse(dat_inv_long)
head(dat_inv_long)


# Question ----------------------------------------------------------------
#
# 1. What does one row of dat_inv_long represent?


# filter out taxa that don't occur in data set ---------------------------

# Identify taxa that occur at least once in the complete data set.
# i.e., if there are taxa columns in the raw data, but end up 
# having 0's in every transects (was never seen)

taxa_present <- dat_inv_long |>
  filter(Abundance > 0) |>
  distinct(Taxa) |>
  pull(Taxa)

# creates a vector of taxa names in the data set
taxa_present

# filter out any taxa that are not present in the data set
dat_inv_long <- dat_inv_long |>
  filter(Taxa %in% taxa_present)

# Check the number of samples, taxa, and problematic values ---------------

data_quality_summary <- dat_inv_long |>
  summarise(
    n_rows = n(),
    n_samples = n_distinct(Year, Site, Replicate),
    n_taxa = n_distinct(Taxa),
    n_missing_abundance = sum(is.na(Abundance))
  )

data_quality_summary


# Identify samples (transects) with no invertebrates (all 0's) --------------

empty_samples <- dat_inv_long |>
  group_by(Year, Site, Replicate) |>
  summarise(
    total_abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(total_abundance == 0)

empty_samples
# if there are any, they should remain in the data for now since not seeing any
# taxa on a transect is a valid observation. and should be included in
#  calculations of mean abundances etc. Although they may need to be 
#  filtered out later in they end up as a weird point in your nMDS etc.

# Summarize sampling effort ------------------------------------------------

site_year_effort <- dat_inv_long |>
  # Each sample occurs once for every taxon in the long table
  distinct(Year, Site, Replicate) |>
  count(Site, Year, name = "n_transects")

site_year_effort

# pivot wider for easier viewing (no longer "tidy" but easier to read)
sampling_effort <- site_year_effort |>
  pivot_wider(
    names_from = Year,
    values_from = n_transects
  )

sampling_effort


# Questions and tasks ------------------------------------------------------

# 1. Was sampling effort equal among all sites and years?
#
# If sampling effort varies (e.g., unequal number of transects among sites or
#  years), then this is an important consideration for later analyses in how you
# calculate and compare diversity metrics across sites, years, etc.

# Consolidate taxonomy ----------------------------------------------------

# Differences in taxonomic names or identification resolution can create
# artificial differences among samples. For example, one survey may identify
# organisms to species while another records only a genus or family. This
# could occur when different people with different taxonomic expertise
# conduct the surveys or process the samples, or when some taxa are difficult
# to identify in the field. In these cases, it may be appropriate to
# consolidate taxa into broader groups for analysis.
#
# Taxa should be consolidated only when there is a biological or sampling
# reason to do so, and careful consideration is needed and depends on your
# research objectives. Preserve the original names, use the finest 
# resolution consistently available across samples, and document 
# every change.
#
# In most studies, this process should use scientific names and an accepted
# taxonomic reference. The broad groups below are simplified examples.
# 
# Most likely you'll want to report the original level of taxa resolution in
# summary/appendix tables in your thesis, and then potentially use a
# consolidated taxa table for your analyses.


# Manually consolidate selected taxa --------------------------------------

# Combine every taxon containing "Urchin" into a single "Sea urchins" group.
# Save the result as a new object so dat_inv_long remains unchanged.

dat_inv_urchin_consolidated <- dat_inv_long |>
  mutate(
    Taxa_Consolidated = case_when(
      # Match any taxon name containing "urchin", regardless of capitalization
      # (i.e.,regex part will match to "Urchin", "URCHIN" as well)
      str_detect(Taxa, regex("urchin", ignore_case = TRUE)) ~ "Sea urchins",
      .default = Taxa
    )
  ) |>
  # Several original taxa may now have the same consolidated name within
  # a sample, so their abundances must be summed.
  group_by(
    Year,
    Site,
    Replicate,
    Group,
    Taxa_Consolidated
  ) |>
  summarise(
    Abundance = sum(Abundance),
    .groups = "drop"
  )

dat_inv_urchin_consolidated

# how many were taxa were changed to "Sea Urchin"
dat_inv_long |> 
  filter(str_detect(Taxa, regex("urchin", ignore_case = TRUE))) |> 
  distinct(Taxa)

# Questions and task -------------------------------------------------------

# 1. How many original taxa were combined into the "Sea urchins" category?
#
# 2. Why was summarise() needed after changing the taxon names?
# 
# 3. Compare the number of obs. between dat_inv_long and
#  dat_inv_urchin_consolidated. Why is the number of rows in the latter
#   smaller?
#
# 4. What information was lost by combining the individual urchin taxa?

# Use a taxonomic lookup table ---------------------------------------------

# A lookup table is preferable when many taxa need to be standardized or
# consolidated. It also creates a permanent record of the decisions.
#
# The lookup table contains:
#   Taxa             = original name in dat_inv_long
#   Taxonomic_Group  = broad group used for this example

taxon_lookup <- read_csv(
  "invertebrate_taxonomy_lookup.csv"
)

glimpse(taxon_lookup)
taxon_lookup


# Check the lookup table before joining -----------------------------------

# Every original taxon should occur exactly once in the lookup table.
# An empty result indicates that there are no duplicated taxon names.

taxon_lookup |>
  count(Taxa) |>
  filter(n > 1)

# Check whether any taxa in the community data are missing from the lookup.
# An empty result indicates that every taxon has a match.

dat_inv_long |>
  distinct(Taxa) |>
  anti_join(taxon_lookup, by = "Taxa")


# Add broad taxonomic groups to the observations --------------------------

dat_inv_with_groups <- dat_inv_long |>
  left_join(
    taxon_lookup,
    by = "Taxa",
  )

glimpse(dat_inv_with_groups)

# Examine the original and consolidated names together
dat_inv_with_groups |>
  distinct(Taxa, Taxonomic_Group) |>
  arrange(Taxonomic_Group, Taxa)

# Questions and task -------------------------------------------------------

# 1. View dat_inv_with_groups and taxon_lookup in the data viewer. Why are 
# there NA's in some rows for the various taxa columns?


# Sum abundance within each broad taxonomic group -------------------------

dat_inv_taxonomic_groups <- dat_inv_with_groups |>
  group_by(
    Year,
    Site,
    Replicate,
    Group,
    Taxonomic_Group
  ) |>
  summarise(
    Abundance = sum(Abundance),
    .groups = "drop"
  )

glimpse(dat_inv_taxonomic_groups)


# Confirm that consolidation did not change total abundance ---------------

dat_inv_long |>
    summarise(total_abundance = sum(Abundance))

dat_inv_taxonomic_groups |>
    summarise(total_abundance = sum(Abundance)) 



# Questions ----------------------------------------------------------------
#
# 1. Why is it useful to retain both Taxa and Taxonomic_Group in
#    dat_inv_with_groups?
#
# 2. Did consolidation change total abundance? Why should it remain
#    unchanged?
#
# We will return to the original taxa in dat_inv_long for the remaining
# tutorial sections.


# Summarize and visualize community data ---------------------------------

# Summary tables are useful for exploring  data and checking that
# later analyses are based on the intended sampling units. Well-organized
# summaries may also be useful as thesis tables (could go in main thesis 
# text or an appendix)

# Summarize sampling effort by site ---------------------------------------

site_sampling_summary <- site_year_effort |>
  group_by(Site) |>
  summarise(
    n_years = n_distinct(Year),
    total_transects = sum(n_transects),
    min_transects_per_year = min(n_transects),
    max_transects_per_year = max(n_transects),
    .groups = "drop"
  ) |>
  arrange(Site)

site_sampling_summary

# Sampling effort differs among some site-year combinations. Therefore,
# comparisons among sites and years should generally use mean abundance per
# transect rather than abundance summed across all transects.

# Summarize total abundance (across all taxa) within each transect --------------------------

transect_totals <- dat_inv_long |>
  group_by(Year, Site, Replicate) |>
  summarise(
    total_abun_trans = sum(Abundance),
    .groups = "drop"
  )

# Calculate mean total abundance per transect for each site and year
site_year_total_summary <- transect_totals |>
  group_by(Site, Year) |>
  summarise(
    n_transects = n(),
    mean_total_abun_trans = mean(total_abun_trans),
    sd_total_abun_trans = sd(total_abun_trans),
    min_total_abundance = min(total_abun_trans),
    max_total_abundance = max(total_abun_trans),
    .groups = "drop"
  ) |>
  arrange(Site, Year)

site_year_total_summary

# Summarize abundance of each taxon ---------------------------------------

# Overall summaries using individual transects as the sampling units.
taxa_overall_summary <- dat_inv_long |>
  group_by(Taxa) |>
  summarise(
    n_transects = n(),
    n_transects_present = sum(Abundance > 0),
    total_abundance = sum(Abundance),
    mean_abun_trans = round(mean(Abundance), 2),
    sd_abundance = round(sd(Abundance), 2 ),
    .groups = "drop"
  ) |>
  mutate(
    proportion_transects_present = (round(n_transects_present / n_transects, 2))
  ) |>
  arrange(desc(mean_abun_trans)) |>
  mutate(
    abundance_rank = row_number(),
    relative_abund_perc = round(total_abundance / sum(total_abundance) * 100, 2)
  )

taxa_overall_summary

# Tasks ----------------------------------------------------------------
#
# 1. Interpret what each column in taxa_overall_summary represents.


# 
# Calculate mean taxon abundance across transects within each site and year.
# this will then be the basis for subsequent averaging to overall site or year means
# so that variation in sampling effort (some sites or years were sampled with more
# transects than others) do not influence the overall means.

taxa_site_year_summary <- dat_inv_long |>
  group_by(Site, Year, Taxa) |>
  summarise(
    n_transects = n_distinct(Replicate),
    mean_abun_trans = mean(Abundance),
    sd_abun_trans = sd(Abundance),
    min_abun_trans = min(Abundance),
    max_abun_trans = max(Abundance),
    .groups = "drop"
  )

taxa_site_year_summary


# Mean abundance by taxon -------------------------------------------------

plot_taxa_mean <- taxa_site_year_summary |>
  group_by(Taxa) |>
  #calculate mean abundance across all site-year combinations
  summarise(
    mean_abun_trans = mean(mean_abun_trans),
    .groups = "drop") |> 
  ggplot(
    aes(
      x = mean_abun_trans,
      y = fct_reorder(Taxa, mean_abun_trans)
    )
  ) +
  geom_col(fill = "steelblue") +
  labs(
    x = "Mean abundance per transect",
    y = NULL
  ) +
  theme_classic()

plot_taxa_mean


#We'll use taxa_overall_summary to then make Taxa a factor with levels 
#in their order of abundance (so they will show up this way in the plots, etc.)
dat_inv_long <- dat_inv_long |> 
  mutate(
    Taxa = factor(Taxa, levels = taxa_overall_summary |> pull(Taxa))
  )

glimpse(dat_inv_long)

# also create factor levels in this table
taxa_site_year_summary <- taxa_site_year_summary |> 
  mutate(
    Taxa = factor(Taxa, levels = taxa_overall_summary |> pull(Taxa))
  )

# Plotting assemblages ---------------------------------------------------

# Create a named vector of distinct colors.
# Naming the colors ensures that each taxon receives the same color in
# every plot, even when some taxa are absent from a particular subset.

taxa_colors <- Polychrome::glasbey.colors(
  length(taxa_present)
) |>
  setNames(taxa_present)

# one color in this palette is white, which is difficult to see 
# on a white background. Replace it with a dark gray.
taxa_colors[toupper(taxa_colors) == "#FFFFFF"] <- "#666666"

taxa_colors

# Can use the same compact taxon legend for all stacked bar plots below. Keeping the
# legend on the right in one column makes it easier to match each color to a
# taxon without using too much plotting space if you need plots to be tall.

compact_taxa_legend <- theme(
  legend.position = "right",
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 7),
  legend.key.height = grid::unit(0.30, "cm"),
  legend.key.width = grid::unit(0.30, "cm"),
  legend.key.spacing.y = grid::unit(0, "cm"),
  legend.margin = margin(0, 0, 0, 0)
)
# NOTE - this is something you can do more generally for ggplot if you want to 
# have standardized formating across multiple plots you are going to make.


# Stacked abundance by individual transect -------------------------------

# Display only the most recent year to keep individual transects readable.
most_recent_year <- max(dat_inv_long$Year)

plot_transect_abundance <- dat_inv_long |>
  filter(Year == most_recent_year) |>
  ggplot(
    aes(
      x = factor(Replicate),
      y = Abundance,
      fill = Taxa
    )
  ) +
  geom_col() +
  scale_fill_manual(
    values = taxa_colors,
    drop = TRUE
  ) +
  facet_wrap(~ Site, nrow = 1, scales = "free_x") +
  labs(
    title = paste("Invertebrate abundance by transect in", most_recent_year),
    x = "Transect",
    y = "Abundance",
    fill = "Taxon"
  ) +
  theme_classic() +
  compact_taxa_legend +
  guides(
    fill = guide_legend(ncol = 1)
  )

plot_transect_abundance

# TO-DO - uncomment the scale_fill_manual() line above to use the taxa_colors palette
# then re-run the plot to see how it changes


# Stacked mean abundance by site and year --------------------------------

plot_site_year_abundance <- taxa_site_year_summary |>
  ggplot(
    aes(
      x = factor(Year),
      y = mean_abun_trans,
      fill = Taxa
    )
  ) +
  geom_col() +
  scale_fill_manual(
    values = taxa_colors,
    drop = TRUE
  ) +
  facet_wrap(~ Site, nrow = 1) +
  labs(
    x = "Year",
    y = "Mean abundance per transect",
    fill = "Taxon"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  compact_taxa_legend +
  guides(
    fill = guide_legend(ncol = 1)
  )

plot_site_year_abundance


# Proportional assemblage composition ------------------------------------

plot_site_year_proportions <- taxa_site_year_summary |>
  ggplot(
    aes(
      x = factor(Year),
      y = mean_abun_trans,
      fill = Taxa
    )
  ) +
  geom_col(position = "fill") +
  scale_fill_manual(
    values = taxa_colors,
    drop = TRUE
  ) +
  facet_wrap(~ Site, nrow = 1) +
  scale_y_continuous(
    labels = scales::label_percent()
  ) +
  labs(
    x = "Year",
    y = "Proportion of mean abundance",
    fill = "Taxon"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  compact_taxa_legend +
  guides(
    fill = guide_legend(ncol = 1)
  )

plot_site_year_proportions


# Overall mean abundance and composition by site -------------------------

# taxa_site_year_summary already contains one mean for each taxon within
# each site and year, calculated across the transects sampled that year.
# Averaging those site-year means across years gives every sampled year at a
# site equal weight, even when the number of transects differs among years.
# Sites with different numbers or sets of sampled years should still be
# compared cautiously.

taxa_site_summary <- taxa_site_year_summary |>
  group_by(Site, Taxa) |>
  summarise(
    n_years = n_distinct(Year),
    mean_abun_trans = mean(mean_abun_trans),
    .groups = "drop"
  ) |>
  arrange(Site, desc(mean_abun_trans))

taxa_site_summary


# Overall mean abundance by site -----------------------------------------

plot_site_abundance <- taxa_site_summary |>
  ggplot(
    aes(
      x = Site,
      y = mean_abun_trans,
      fill = Taxa
    )
  ) +
  geom_col() +
  scale_fill_manual(
    values = taxa_colors,
    drop = TRUE
  ) +
  labs(
    x = "Site",
    y = "Mean abundance per transect",
    fill = "Taxon"
  ) +
  theme_classic() +
  compact_taxa_legend +
  guides(
    fill = guide_legend(ncol = 1)
  )

plot_site_abundance


# Overall proportional assemblage composition by site -------------------

plot_site_proportions <- taxa_site_summary |>
  ggplot(
    aes(
      x = Site,
      y = mean_abun_trans,
      fill = Taxa
    )
  ) +
  geom_col(position = "fill") +
  scale_fill_manual(
    values = taxa_colors,
    drop = TRUE
  ) +
  scale_y_continuous(
    labels = scales::label_percent()
  ) +
  labs(
    x = "Site",
    y = "Proportion of mean abundance",
    fill = "Taxon"
  ) +
  theme_classic() +
  compact_taxa_legend +
  guides(
    fill = guide_legend(ncol = 1)
  )

plot_site_proportions


# Questions and task ------------------------------------------------------

# 1. Which sites have the greatest and lowest overall mean abundance?
#    Which taxa account for most of that difference?
#
# 2. Do the abundance and proportional plots lead to the same conclusions
#    about differences among sites?
#
# TASK: Copy and paste the code used to create taxa_site_summary and its two
# plots. Revise the copies so that each bar represents a year, with taxon
# abundance averaged across sites. Begin with taxa_site_year_summary, group
# by Year and Taxa, and record the number of sites contributing to each
# annual mean. Use new object and plot names so you do not overwrite the
# site-level versions. Why is it important to average the site-year means
# rather than pooling all transects within each year?


# Abundance patterns for common taxa -------------------------------------

# Select the six taxa with the greatest overall mean abundance.
top_taxa <- taxa_overall_summary |>
  slice_max(
    order_by = mean_abun_trans,
    n = 6,
    with_ties = FALSE
  ) |>
  pull(Taxa)

plot_common_taxa <- taxa_site_year_summary |>
  filter(Taxa %in% top_taxa) |>
  ggplot(
    aes(
      x = Year,
      y = mean_abun_trans,
      color = Site,
      group = Site
    )
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  facet_wrap(
    ~ Taxa,
    scales = "free_y"
  ) +
  scale_x_continuous(
    breaks = sort(unique(dat_inv_long$Year))
  ) +
  labs(
    x = "Year",
    y = "Mean abundance per transect",
    color = "Site"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot_common_taxa


# Questions and tasks ------------------------------------------------------

# 1. What do the mean-abundance plots reveal about
#    dominant and uncommon taxa in this assemblage?
#
# 2. Compare the absolute and proportional stacked bar plots. Identify a
#    site-year combination for which the two plots suggest different
#    impressions of assemblage change. Why?
#
# TASK: Revise top_taxa to display the four taxa that interest you most,
# rather than automatically selecting the six most abundant taxa.
#
# TASK: Modify plot_common_taxa so that each site is displayed in a separate
# facet. Decide whether Taxa or Site should define the rows and columns.


# Create a community input data table for vegan -------------------------------

# Most vegan functions require a community data table in which:
#   - each row is a unique sampling unit
#   - each column is a taxon
#   - each cell contains the abundance of that taxon

# Convert the long data back to wide format
wide_dat_inv <- dat_inv_long |>
  pivot_wider(
    id_cols = c(Year, Site, Replicate, Group),
    names_from = Taxa,
    values_from = Abundance,
    values_fill = 0
  ) |>
  # create a unique sample identifier for each transect used later as row names
  mutate(
    Sample_ID = str_c(Year, Site, Replicate, sep = "_"),
    .before = 1
  ) |>
  arrange(Year, Site, Replicate)

glimpse(wide_dat_inv)


# Create the community data table ----------------------------------------

# vegan functions require the community table to contain only numeric taxon
# abundances. Unique sample identifiers (in this case transects) need to be 
# stored as row names so the results can later be matched back to the sample metadata.
#
# column_to_rownames() converts the tibble to a data frame because tibbles
# do not support row names.

comm_dat_inv <- wide_dat_inv |>
  # remove columns that are not taxa or the unique sample identifier
  select(
    -Year,
    -Site,
    -Replicate,
    -Group
  ) |>
  column_to_rownames(var = "Sample_ID")

# Examine the dimensions and a small portion of the community table
dim(comm_dat_inv)
view(comm_dat_inv)


# Questions ----------------------------------------------------------------

# 1. What does each row of comm_dat_inv represent?
#
# 2. What does each column of comm_dat_inv represent?
#
# 3. Why were Year, Site, Replicate, and Group removed from comm_dat_inv?


# Calculate univariate diversity metrics ---------------------------------

# vegan includes functions for calculating taxon richness and several
# abundance-based measures of diversity.
#
# Additional references:
# https://cran.r-project.org/web/packages/vegan/vignettes/diversity-vegan.pdf
# https://www.webpages.uidaho.edu/veg_measure/modules/lessons/module%209(composition&diversity)/9_3_Estimating%20Biodiversity.htm
# https://en.wikipedia.org/wiki/Diversity_index
#
# Calculate these metrics using the complete community table before removing
# rare taxa. Rare taxa are important components of richness and diversity.
#
# Because some categories in this dataset are broader than individual
# species, we refer to this metric as taxon richness rather than species
# richness.

diversity_metrics_inv <- tibble(
  Sample_ID = rownames(comm_dat_inv),

  # Total abundance across all taxa
  total_abundance = rowSums(comm_dat_inv),

  # vegan package function specnumber: Number of taxa observed in 
  # each transect (row in data table)
  taxa_richness = specnumber(comm_dat_inv),

  # vegan package function diversity: Shannon diversity incorporates 
  # richness and evenness and is relatively sensitive to uncommon and 
  # moderately abundant taxa.
  shannon_diversity = diversity(
    comm_dat_inv,
    index = "shannon"
  ),

  # vegan's Simpson index is 1 - D. It places greater weight on common taxa
  # and ranges from 0 toward 1 as diversity increases.
  simpson_diversity = diversity(
    comm_dat_inv,
    index = "simpson"
  ),

  # Inverse Simpson can be interpreted as the effective number of common or
  # dominant taxa and is often easier to compare among samples.
  inverse_simpson = diversity(
    comm_dat_inv,
    index = "invsimpson"
  )
) |>
  mutate(
    # Pielou's evenness describes how evenly abundance is distributed among
    # the taxa present. It is undefined when richness is 0 or 1.
    pielou_evenness = if_else(
      taxa_richness > 1,
      shannon_diversity / log(taxa_richness),
      NA_real_
    )
  )

glimpse(diversity_metrics_inv)


# Add sample metadata columns back to the diversity results table -----
dat_inv_diversity <- wide_dat_inv |>
  select(Sample_ID, Year, Site, Replicate, Group) |> 
  left_join(diversity_metrics_inv,
    by = "Sample_ID",
    relationship = "one-to-one"
  )

glimpse(dat_inv_diversity)


# Summarize diversity by site and year ------------------------------------

# Diversity was calculated separately for each transect and is then
# summarized across replicate transects. This estimates mean transect-level
# diversity rather than the diversity of all observations pooled together.
# This is one common approach to account for differences in sampling effort 
# (differences among sites and years in this data set)

diversity_site_year_summary <- dat_inv_diversity |>
  group_by(Site, Year) |>
  summarise(
    n_transects = n(),
    across(
      c(
        total_abundance,
        taxa_richness,
        shannon_diversity,
        simpson_diversity,
        inverse_simpson,
        pielou_evenness
      ),
      list(
        mean = \(x) mean(x, na.rm = TRUE),
        sd = \(x) sd(x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),

    .groups = "drop"
  ) |>
  mutate(
    shannon_diversity_se =
      shannon_diversity_sd / sqrt(n_transects)
  ) |>
  arrange(Site, Year)

glimpse(diversity_site_year_summary)


# Calculate an overall summary for each site ------------------------------

diversity_site_summary <- dat_inv_diversity |>
  group_by(Site) |>
  summarise(
    n_years = n_distinct(Year),
    n_transects = n(),

    across(
      c(
        total_abundance,
        taxa_richness,
        shannon_diversity,
        simpson_diversity,
        inverse_simpson,
        pielou_evenness
      ),
      \(x) mean(x, na.rm = TRUE),
      .names = "mean_{.col}"
    ),

    .groups = "drop"
  ) |>
  arrange(desc(mean_shannon_diversity))

glimpse(diversity_site_summary)


# Questions ----------------------------------------------------------------

# 1. Which site-year combination had the greatest mean taxon richness?
#
# 2. Which site had the greatest overall mean Shannon diversity?
#
# 3. Do richness, Shannon diversity, and Simpson diversity rank the sites in
#    the same order? Why might these metrics produce different rankings?
#
# 4. Why did we calculate diversity for individual transects before
#    summarizing the results by site and year?


# Taxon richness among sites and years ------------------------------------

# Small points show individual transects. Larger white points show mean
# richness for each site-year combination.

plot_taxa_richness <- dat_inv_diversity |>
  ggplot(
    aes(
      x = Site,
      y = taxa_richness
    )
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    alpha = 0.6,
    size = 2,
    color = "steelblue"
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 3,
    fill = "white",
    color = "black"
  ) +
  facet_wrap(~ Year) +
  labs(
    x = "Site",
    y = "Taxon richness per transect"
  ) +
  theme_classic()

plot_taxa_richness


# Mean Shannon diversity through time ------------------------------------

# Points and lines show site-year means. Error bars show plus or minus one
# standard error across replicate transects. These are descriptive summaries
# and do not account for the repeated sampling of sites through time.

plot_shannon_diversity <- diversity_site_year_summary |>
  ggplot(
    aes(
      x = Year,
      y = shannon_diversity_mean,
      color = Site,
      group = Site
    )
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  scale_x_continuous(
    breaks = sort(unique(dat_inv_diversity$Year))
  ) +
  labs(
    x = "Year",
    y = "Mean Shannon diversity",
    color = "Site"
  ) +
  theme_classic()

plot_shannon_diversity

