# Practice nMDS analysis: Oribatid mite assemblages -----------------------

# Dataset description:
# https://www.davidzeleny.net/anadat-r/doku.php/en:data:mite
#
# This terrestrial dataset contains abundances of 35 oribatid mite taxa from
# 70 soil cores, along with environmental information for each sample.


# Load packages -----------------------------------------------------------

library(tidyverse)
library(vegan)


# Load the example datasets ----------------------------------------------

# mite contains the taxon abundances.
# mite.env contains the environmental variables for the same samples.
data("mite", package = "vegan")
data("mite.env", package = "vegan")

glimpse(mite)
glimpse(mite.env)


# Prepare the data --------------------------------------------------------

# Convert the row names into an explicit sample identifier. The rows of the
# abundance and environmental datasets refer to the same soil-core samples.
mite_abundance <- mite |>
  rownames_to_column(var = "Sample_ID") |>
  as_tibble()

mite_metadata <- mite.env |>
  rownames_to_column(var = "Sample_ID") |>
  as_tibble()

# Join the metadata and taxon abundances into one wide dataset. Metadata
# columns occur first, followed by one abundance column for each mite taxon.
mite_wide <- mite_metadata |>
  left_join(
    mite_abundance,
    by = "Sample_ID"
  )

glimpse(mite_wide)


# Adapt the nMDS tutorial -------------------------------------------------

# Use the completed coral invertebrate nMDS tutorial as your starting point.
# Copy, paste, and revise its code to analyze the mite assemblage data.
#
# Sampling unit:
# - Each row represents one soil core, and each nMDS point should represent
#   one soil-core sample.
# - These data do not contain replicate transects to average (i.e., omit the
#   site-year averaging step used for the coral invertebrate data).
#
# Community data:
# - Sample_ID, SubsDens, WatrCont, Substrate, Shrub, and Topo are metadata.
#   DO NOT include these columns in the sample-by-taxa community matrix.
# - The remaining columns contain mite abundances.
# - Work through the same decisions and checks used in the tutorial,
#   including choice of transformation, optional rare-taxon filtering,
#   nMDS stress, convergence, and the stress plot.
#
# Categorical explanatory variables:
# - Substrate identifies seven substrate categories.
# - Shrub is an ordered categorical variable with three levels.
# - Topo identifies two microtopographic categories: Blanket and Hummock.
# 
# - Try different combinations of these variables for point color, shape,
#   labels, ellipses, or minimum convex polygons. Choose mappings that remain
#   readable and that address a clear ecological comparison.
#
# Numeric explanatory variables:
# - SubsDens measures substrate density.
# - WatrCont measures substrate water content.
# 
# - Use envfit() to fit these two continuous environmental variables to the
#   nMDS configuration. Interpret an arrow as the direction in which that
#   variable increases and its fit as the strength of its association with
#   the ordination configuration.
#   
# - Also try mapping either SubsDens or WatrCont to point size. Use point size
# to show one variable at a time so the figure remains interpretable.
#
# Interpretation:
# - Describe the main assemblage patterns shown by distances among samples.
# - Evaluate whether samples separate or overlap among substrate, shrub, or
#   microtopography categories.
# - Evaluate whether substrate density or water content is associated with
#   the observed assemblage pattern. An envfit association does not by itself
#   demonstrate that an environmental variable caused the pattern.
# - Identify which figures are exploratory and which could communicate the
#   primary biological result in a report, thesis, or paper.



