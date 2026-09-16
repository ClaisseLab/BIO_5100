# (BIO 5150L) Principal Component Analysis (PCA) tutorial
# Darlingtonia plant morphology and biomass among sites
# Date modified: 15 Sep 2026
#
# PCA is commonly used to describe correlation structure (relationships) among
# multiple different variables, e.g. environmental variables measured for each
# sample, site, plot, transect or individual, or species characteristics
# (traits) measured for individual species. Variables are typically centered
# and scaled when they have different units or ranges; otherwise, variables
# with greater variation can dominate the analysis.
#
# PCA concepts and scaling:
# https://www.davidzeleny.net/anadat-r/doku.php/en:pca
# 
# Online course including detailed PCA background (good background info, but
# note it uses a different R package and function to do the PCA):
# https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/pca/
#
# Data source:
# Required file in the same folder:
#   Darlingtonia_GE_Table12.1.csv
#   
# Original source: https://github.com/jon-bakker/appliedmultivariatestatistics/blob/main/Darlingtonia_GE_Table12.1.csv
#
# PCA is used here to summarize multiple continuous morphological variables.
# The same workflow is commonly used for environmental, habitat, and
# physiological variables. In contrast, the nMDS tutorial uses species-
# abundance data to describe differences in community composition.


# Load packages -----------------------------------------------------------

library(tidyverse)
library(vegan)
library(ggrepel)

# Used only for the clustered correlation plot below. Install it once if it
# is not already installed:
# install.packages("corrplot")


# Import and inspect the data ---------------------------------------------

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

glimpse(dat_darl)
head(dat_darl)

# Key points (Review with dat_darl before starting the analysis):
# - One row represents one individual Darlingtonia plant.
# - sample_ID uniquely identifies each plant.
# - site is categorical metadata. It can be used to group or color samples in
# plots, but it is not included as a variable in the PCA.
# - The remaining ten columns are continuous morphology and biomass variables
#   measured on each plant.
# - Measurements use different units and numerical scales, which will affect
# the decision about standardization.

sample_ID_check <- dat_darl |>
  summarise(
    n_rows = n(),
    n_unique_sample_IDs = n_distinct(sample_ID)
  )

sample_ID_check


# Identify the variables to include in the PCA ----------------------------

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

# Before conducting a PCA, clearly distinguish among:
# - identifiers, such as sample_ID and plant;
# - grouping variables, such as site; and
# - continuous variables (measurements on each plant) included in the PCA.
#
# Do not include categorical variables simply by converting their factor
# levels to numbers. Also avoid including a derived total together with all
# of the component variables used to calculate it. That would give one aspect
# of the data extra weight.


# Check missing values and variable distributions -------------------------

pca_variable_summary <- dat_darl |>
  select(all_of(pca_variables)) |>
  # convert data from wide to long format for analysis w/ tidyverse functions
  pivot_longer(
    cols = everything(),
    names_to = "Variable",
    values_to = "Value"
  ) |>
  group_by(Variable) |>
  summarise(
    n_observed = sum(!is.na(Value)),
    n_missing = sum(is.na(Value)),
    mean = mean(Value, na.rm = TRUE),
    median = median(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    minimum = min(Value, na.rm = TRUE),
    maximum = max(Value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(Variable)

pca_variable_summary

# Save the identity of any samples containing missing PCA measurements.
samples_with_missing_pca_data <- dat_darl |>
  filter(if_any(all_of(pca_variables), is.na)) |>
  select(sample_ID, site, plant)

samples_with_missing_pca_data

# rda() cannot fit this PCA with missing measurements. In this example there is
# no missing data, but the code to do so is included here anyway. With your own
# data, first determine why values are missing and whether exclusion or a
# justified imputation method (i.e., fill in missing values) is appropriate.
dat_darl_complete <- dat_darl |>
  drop_na(all_of(pca_variables))

# View the distribution of each variable.
plot_pca_variable_distributions <- dat_darl_complete |>
  select(all_of(pca_variables)) |>
  pivot_longer(
    cols = everything(),
    names_to = "Variable",
    values_to = "Value"
  ) |>
  mutate(
    Variable = Variable |>
      str_replace_all("_", " ") |>
      str_to_sentence()
  ) |>
  ggplot(aes(x = Value)) +
  geom_histogram(
    bins = 15,
    color = "white",
    fill = "steelblue"
  ) +
  facet_wrap(
    vars(Variable),
    scales = "free",
    ncol = 3
  ) +
  labs(
    x = "Observed value",
    y = "Number of plants"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 8)
  )

plot_pca_variable_distributions # EXPLORATORY/CONFIRMATION PLOT

# Key points:
# - Look for strongly skewed variables, extreme observations, impossible
# values, and very different numerical ranges.
# - The centering and scaling procedure (done later) will place variables
#   on comparable scales (mean = 0 and sd = 1), but it does not correct skew or
#   reduce the influence of extreme observations.
# - A scientifically defensible transformation may be appropriate before
#   PCA when variables are strongly skewed, although it is not always necessary.
#   Do not transform a variable only to make the PCA produce a preferred pattern.
# - This distribution plot is a data-quality/exploration figure. It potentially
# could be a primary thesis or paper result/figure if describing the 
# distributions of each variable is an important aspect of the study, however in
# some cases it may not be included in your results.


# Construct the sample-by-variable data frame  -----------------------------
# (i.e., "comm" table format vegan package functions require as inputs)
# 
# Keep sample_ID as row names so sample identities remain attached to the
# PCA matrix and can be recovered when scores are extracted.
comm_darl <- dat_darl_complete |>
  select(sample_ID, all_of(pca_variables)) |>
  column_to_rownames(var = "sample_ID")

glimpse(comm_darl)
view(comm_darl) # notice sample IDs are now row names

# Examine correlations among variables -----------------------------------

cor_darl <- comm_darl |>
  cor()

round(cor_darl, digits = 2)

# Arrange correlated variables together using hierarchical clustering and
# draw boxes around three visually useful groups. The number of rectangles
# is selected for display and is not a formal test of the number of variable
# groups.


# Format variable names for the correlation plot.
cor_variable_labels <- colnames(cor_darl) |>
  str_replace("_mm$", " (mm)") |>
  str_replace("_g$", " (g)") |>
  str_replace_all("_", " ") |>
  str_to_sentence()

# Add the formatted names to a copy of the correlation matrix.
cor_darl_plot <- cor_darl

dimnames(cor_darl_plot) <- list(
  cor_variable_labels,
  cor_variable_labels
)

n_correlation_groups <- 3 # pick by trial and error (see what appears to create the most cohesive groups of correlated variables)

# Plot the correlation matrix.
corrplot::corrplot(
  cor_darl_plot,
  method = "color",
  order = "hclust", # runs a cluster analysis to order the variables by strength of correlations with each other
  hclust.method = "ward.D2",
  addrect = n_correlation_groups,
  rect.col = "black",
  rect.lwd = 1.5,
  addCoef.col = "black",
  number.digits = 2,
  number.cex = 0.5,
  tl.col = "black",
  tl.srt = 45,
  tl.cex = 0.7,
  diag = TRUE
)

# Key points to interpret:
# - Variables with strong positive correlations increase together and may
#   contribute to the same biological gradient.
# - Strong negative correlations indicate opposing patterns.
# - Weakly correlated variables may contribute to different principal
#   components.
# - PCA is most effective as a dimension-reduction method when some variables
#   contain overlapping information.
# - A clustered correlation plot is useful to identify and describe patterns
#  (relationships among the variables), but the boxes are descriptive rather
#  than inferential (i.e., they are for visualization purposes only - no
#  statistical test is run).


# Brief connection to Euclidean distance ---------------------------------

# PCA represents each sample in a Euclidean ("straight line distance") variable
# space. Samples with similar combinations of measurements occur near one
# another, while samples with different measurements occur farther apart.
#
# When all principal components are retained, PCA preserves the Euclidean
# geometry of the centered data. A PC1 versus PC2 plot is an approximation
# because it displays only the portion of total variation represented by
# those two axes.
#
# When variables are standardized, distances reflect differences measured in
# standard deviations rather than differences in the original units.


# Main analytical choice: Center and scale the variables -----------------

# rda() always centers variables for a PCA. Setting scale = TRUE also divides
# each variable by its standard deviation. This gives every variable equal
# variance before the PCA and is equivalent to conducting the PCA from a
# correlation matrix.
#
# Centering and scaling are usually appropriate when variables have different
# units or very different numerical ranges. Otherwise, a variable measured on
# a larger scale or with greater variance can dominate the PCA.
#
# Using scale = FALSE produces a covariance-based PCA. This can be appropriate
# when all variables use comparable units and scales and differences in their
# variances are biologically meaningful. That should be an intentional choice,
# not a default.
#
# The morphology variables include measurements in millimeters and masses in
# grams, so the main analysis below uses scale = TRUE.
# 
# Note: You would also typically center and scale variables if you are going to
# use them as explanatory variables in a multiple-regression analysis. In that
# case, you would do it manually and use the centered and scaled variables in
# the regression. In contrast, rda() performs the centering and scaling for this
# PCA when scale = TRUE.



# Fit an unconstrained PCA with vegan::rda() ------------------------------

pca_darl <- rda(
  comm_darl,
  scale = TRUE #center and scale all variables to mean = 0 and sd = 1
)

pca_darl
summary(pca_darl, scaling = 2)

# Key points to interpret:
# - Supplying only the continuous data frame fits an "unconstrained" PCA.
# - scaling = 2 changes the displayed scores, not the fitted PCA or eigenvalues;
# it emphasizes relationships among variables in the biplot.
# - PCA is deterministic (calculated) and does not use random starting
# configurations, so set.seed() is not needed.
# - The grouping variable site was not used to calculate the PCA. It will be
# added later only to help interpret patterns among plants.
# - Other R functions and packages can fit PCA. We use vegan::rda() to
# maintain a consistent ordination framework across this multivariate module (it
# is main package we use for nMDS), but there may be benefits in some cases to
# using a different function to fit your PCA. More options described here:
# https://www.davidzeleny.net/anadat-r/doku.php/en:pca_r

# Calculate variance explained by each principal component ---------------

pca_eigenvalues <- eigenvals(pca_darl)

pca_variance_table <- tibble(
  PC = names(pca_eigenvalues),
  Eigenvalue = as.numeric(pca_eigenvalues)
) |>
  mutate(
    Proportion_explained = Eigenvalue / sum(Eigenvalue),
    Percent_explained = 100 * Proportion_explained,
    Cumulative_percent = cumsum(Percent_explained)
  ) |>
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, digits = 2)
    )
  )

pca_variance_table

# Eigenvalues measure how much variance is represented by each principal
# component. Each eigenvector contains the coefficients used to construct a
# PC from the original variables. See the linked course chapter for a more
# detailed explanation of the underlying calculations.


# Create a scree plot -----------------------------------------------------

plot_pca_scree <- pca_variance_table |>
  mutate(
    PC = factor(PC, levels = PC)
  ) |>
  ggplot(aes(x = PC, y = Percent_explained, group = 1)) +
  geom_col(fill = "steelblue", width = 0.75) +
  geom_line(color = "black", linewidth = 0.6) +
  geom_point(color = "black", size = 2) +
  labs(
    x = "Principal component",
    y = "Variance explained (%)"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

plot_pca_scree # EXPLORATORY/CONFIRMATION PLOT

# Key points to interpret:
# - The scree plot helps evaluate dimensionality but is usually an analytical
# check rather than a primary results figure.
# - PC1 explains the greatest possible amount of variation, and each later PC
# (principal component) explains the greatest possible amount of the remaining
# variation that is orthogonal (perpendicular) to the previous PCs. 
# - Examine the percentage and cumulative percentage explained by PC1 and PC2.
# - The first two PCs are emphasized because they provide the most informative
# two-dimensional view, but they may not capture every important pattern.
# - Other PC combinations can be examined when justified, but they are not
# plotted in this introductory tutorial.
# (for example sometimes PC1 and PC3 are also plotted together and presented
# alongside the typical PC1 vs. PC2 plot, if PC3 explains a relatively high
# amount of variation and interpretation reveals it is also important to the
# biological interpretation).



# Prepare axis labels containing the percentage explained by each PC ---------

pc1_percent <- pca_variance_table |>
  filter(PC == "PC1") |>
  pull(Percent_explained) 

pc2_percent <- pca_variance_table |>
  filter(PC == "PC2") |>
  pull(Percent_explained)

pc1_axis_label <- sprintf("PC1 (%.0f%%)", pc1_percent)
pc1_axis_label
pc2_axis_label <- sprintf("PC2 (%.0f%%)", pc2_percent)
pc2_axis_label

# Extract sample scores and join metadata --------------------------------

# scaling = 2 emphasizes relationships among the original variables and is
# commonly used for a correlation-focused PCA biplot. Use the same scaling
# when extracting both sample scores and variable scores.
pca_scores_darl <- vegan::scores(
  pca_darl,
  display = "sites",
  choices = c(1, 2),
  scaling = 2
) |>
  as.data.frame() |>
  rownames_to_column(var = "sample_ID") |>
  as_tibble() |>
  left_join(
    dat_darl_complete,
    by = "sample_ID"
  )

glimpse(pca_scores_darl)

# Key points to interpret:
# - PC1 and PC2 are sample scores (point locations).
# - Scores are centered around zero. Positive and negative values indicate
#   opposite ends of a gradient, not good versus bad conditions.
# - The sign of an entire PC can reverse without changing the PCA. A reflected
#   plot represents the same relationships among samples and variables.
# - Joining the original data and metadata allows site, plant identity, or an
#   original measurement to be mapped onto a plot. These joined columns were
#   not necessarily used to calculate the PCA.


# Define a consistent PCA plot theme -------------------------------------

pca_plot_theme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.box = "vertical"
  )


# Plot sample scores ------------------------------------------------------

plot_pca_samples <- pca_scores_darl |>
  ggplot(aes(x = PC1, y = PC2)) +
  geom_hline(
    yintercept = 0,
    color = "grey75",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = 0,
    color = "grey75",
    linewidth = 0.5
  ) +
  geom_point(
    aes(color = site, shape = site),
    size = 3,
    alpha = 0.8
  ) +
  coord_equal() +
  labs(
    x = pc1_axis_label,
    y = pc2_axis_label,
    color = "Site",
    shape = "Site"
  ) +
  pca_plot_theme

plot_pca_samples # POTENTIAL PRIMARY RESULTS PLOT

# Key points to interpret:
# - With scaling = 2, the plot prioritizes relationships among variable vectors.
# - Relative sample positions can reveal broad patterns, but distances among
#   plotted samples do not directly approximate their original Euclidean distances.
# - Use scaling = 1 when the primary goal is to approximate Euclidean distances
#   among samples in the ordination plot.
# - Look for separation, overlap, and within-site variation rather than only
#   whether groups form completely distinct clusters.
# - Because site was added after fitting the PCA, apparent site separation is
#   descriptive and is not a statistical test of site differences.
# - The axis labels report how much of the total standardized variation is
#   represented in this two-dimensional plot.


# Option: Label individual samples ---------------------------------------

plot_pca_sample_labels <- plot_pca_samples +
  geom_text_repel(
    aes(label = sample_ID, color = site),
    size = 2.5,
    show.legend = FALSE,
    max.overlaps = Inf
  )

plot_pca_sample_labels # EXPLORATORY PLOT

# Labels are useful for identifying unusual observations and checking the
# original data. They may be too crowded for a final figure.



# Extract variable scores -------------------------------------------------

# vegan retains community-ecology terminology and calls these "species"
# scores. In this PCA they represent the ten morphology and biomass variables,
# not species.
pca_variable_scores_darl <- vegan::scores(
  pca_darl,
  display = "species",
  choices = c(1, 2),
  scaling = 2
) |>
  as.data.frame() |>
  rownames_to_column(var = "Variable") |>
  as_tibble() |>
  mutate(
    Variable_label = Variable |>
      str_replace_all("_", " ") |>
      str_to_sentence()
  )

pca_variable_scores_darl


# Scale variable arrows to fit the ggplot --------------------------------

# The sample and variable scores use the same PCA scaling, but an additional
# constant is needed to make the arrows readable within the plotting region.
# Multiplying every arrow by the same constant changes only their displayed
# size, not their directions, relative lengths, or interpretations.
# You can adjust the multiplier up or down to make the arrow lengths fit better
# in your plot
arrow_multiplier <- 0.70 * min(
  diff(range(pca_scores_darl$PC1)) /
    diff(range(pca_variable_scores_darl$PC1)),
  diff(range(pca_scores_darl$PC2)) /
    diff(range(pca_variable_scores_darl$PC2))
)

pca_variable_scores_plot <- pca_variable_scores_darl |>
  mutate(
    PC1_plot = PC1 * arrow_multiplier,
    PC2_plot = PC2 * arrow_multiplier
  )


# Add variable vectors to create a PCA biplot ----------------------------

plot_pca_biplot <- plot_pca_samples +
  geom_segment(
    data = pca_variable_scores_plot,
    aes(
      x = 0,
      y = 0,
      xend = PC1_plot,
      yend = PC2_plot
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.6,
    arrow = grid::arrow(
      length = grid::unit(0.18, "cm")
    )
  ) +
  geom_text_repel(
    data = pca_variable_scores_plot,
    aes(
      x = PC1_plot,
      y = PC2_plot,
      label = Variable_label
    ),
    inherit.aes = FALSE,
    color = "black",
    size = 3,
    fontface = "bold",
    min.segment.length = 0,
    max.overlaps = Inf
  )

plot_pca_biplot # PRIMARY RESULTS PLOT

# Key points to interpret:
# - Each arrow points in the direction where that variable increases.
# - Variables with arrows pointing in similar directions are positively
#   correlated. Variables pointing in opposite directions are negatively
#   correlated. Variables near right angles are weakly correlated.
# - Plants in the direction of an arrow tend to have higher values of that
#   variable relative to other plants.
# - Longer arrows are more strongly represented by PC1 and PC2. Because all
#   arrows were multiplied by the same display constant, compare their lengths
#   relative to one another rather than reading them as raw correlations.
# - Interpret the variable arrows together with the sample positions. First
#   describe the main biological gradient, then describe how sites or samples
#   are distributed along it.


# Calculate correlations between variables and PCs -----------------------

# These correlations provide a compact numerical interpretation of the same
# relationships displayed by the biplot vectors. Correlation is unchanged by
# multiplying a PC score by a positive scaling constant.
pca_variable_pc_correlations <- cor(
  comm_darl,
  vegan::scores(
    pca_darl,
    display = "sites",
    choices = c(1, 2),
    scaling = 0
  )
) |>
  as.data.frame() |>
  rownames_to_column(var = "Variable") |>
  as_tibble() |>
  mutate(
    Absolute_PC1_correlation = abs(PC1),
    Absolute_PC2_correlation = abs(PC2)
  ) |>
  arrange(desc(Absolute_PC1_correlation)) |>
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, digits = 2)
    )
  )

pca_variable_pc_correlations

# Key points to interpret:
# - Large positive correlations identify variables that increase toward the
#   positive end of a PC; large negative correlations increase toward its
#   negative end.
# - Variables near zero are weakly associated with that PC.
# - Several variables with large correlations in the same direction can
#   define a shared biological gradient, such as overall plant size.
# - Variables with large correlations in opposite directions identify a
#   contrast between the two ends of a PC.
# - There is no universal loading or correlation threshold that defines an
#   important variable. Interpret magnitude, consistency, biology, and the
#   percentage of variation explained together.
# - This table is useful for interpreting the PCA and may be reported in the
#   main text, a table, or an appendix.

# Option: Add site ellipses ----------------------------------------------

# Ellipses may be useful when there are many observations per group
# (e.g., >10–20) and you want to summarize each group's distribution.

plot_pca_site_ellipses <- plot_pca_samples +
  stat_ellipse(
    aes(
      fill = site,
      group = site
    ),
    geom = "polygon",
    alpha = 0.15,
    color = "white",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  # Replot points over the filled ellipses
  geom_point(
    aes(
      color = site,
      shape = site
    ),
    size = 3,
    alpha = 0.8,
    show.legend = FALSE
  )

plot_pca_site_ellipses # POTENTIAL RESULTS PLOT

# Ellipses summarize the location and spread of each site in the displayed
# PC space. They help show overlap and within-site variation, but they do not
# by themselves test for a difference among sites.




# Summarize PC scores by site ---------------------------------------------

pca_site_summary <- pca_scores_darl |>
  group_by(site) |>
  summarise(
    n_plants = n(),
    mean_PC1 = mean(PC1),
    sd_PC1 = sd(PC1),
    mean_PC2 = mean(PC2),
    sd_PC2 = sd(PC2),
    .groups = "drop"
  ) |>
  arrange(site)

pca_site_summary

# Site means summarize locations along the PCA gradients, while the standard
# deviations describe variation among individual plants. Differences among
# site means should be interpreted alongside the overlap visible in the plot.


# Retain PC scores as combined biological metrics ------------------------

dat_darl_with_pca_scores <- pca_scores_darl |>
  relocate(sample_ID, site, plant, PC1, PC2)

glimpse(dat_darl_with_pca_scores)

# PC1 or PC2 can sometimes serve as a combined metric in a later analysis.
# This is appropriate only when the variables associated with that PC form a
# coherent biological gradient and the PC explains a meaningful portion of
# their variation.
#
# Important considerations:
# - PCA is unsupervised. It summarizes variation in the PCA variables without
#   considering a later response variable.
# - Report which variables define each PC and how much variation that
#   PC explains.
# - Do not include both a PC and the original variables used to construct it
#   as predictors in the same (future) model.
# - If an entire PC is multiplied by -1 to make its direction more intuitive,
#   document the change. Reversing the sign does not change the PCA.


# Methods and results reporting guidance ---------------------------------

# A methods description should report:
# - the analytical unit and number of samples;
# - the continuous variables included and any variables excluded;
# - transformations if used before the PCA on certain variables and biological
#  justification for doing so;
# - whether variables were centered and scaled;
# - the scaling used to extract and display PCA scores;
# - use of vegan::rda() to create the PCA (w/ citation for the vegan package);
# - whether any retained PC scores were used in later analyses (e.g., PC1 was
#   used as a predictor in a subsequent model).
# 
# A results description should report:
# - the percentage of variation explained by PC1 and PC2;
# - the main variables associated with each displayed PC;
# - a biological interpretation of the main gradient;
# - separation, overlap, and within-group variation among samples; and
# - important exceptions or unusual samples.
#
# Describe the broad biological pattern first, followed by the quantitative
# support from variance explained and variable-PC correlations. Avoid treating
# visually separated groups as a formal statistical result unless they are
# evaluated with an analysis that matches the study design.
