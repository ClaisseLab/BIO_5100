# REMEMBER TO FIRST RENAME THE .R File name then delete this comment

# Author: add your name and maybe email here

# Date Created: DD MMM YYYY


# Purpose: 
# Briefly describe what this script does.

# TO-DO/Questions/issues:
# Note next steps, unresolved questions, problems, or decisions here.


# Load packages ---------------------------------------------------------------

library(tidyverse)

# Load data files -------------------------------------------------------------

# Example:
# survey_raw <- read_csv("survey_data.csv")


# Inspect data ---------------------------------------------------------------

# glimpse(survey_raw)




# Filter and tidy data --------------------------------------------------------

# Example:
# survey_clean <- survey_raw |>
#   filter(...) |>
#   mutate(...)


# Analysis --------------------------------------------------------------------

# Explore the data, create figures, and conduct statistical analyses here.


# Save outputs ---------------------------------------------------------------

# Revise and insert this code throughout the script as needed to save tables,
# figures, cleaned data, or other outputs.

# Save a table:
# results_table |>
#  write_csv("results_table.csv")

# Or if you using file structure for GitHub or best practice and want to write
# the file to a "tables" subfolder that is in your main project folder:
# results_table |>
#   write_csv("tables/results_table.csv")

# Save a ggplot as a .png file:
# count_temp_plot |>
#   ggsave(
#     filename = "count_temp_plot.png", 
#     # or if writing to figures folder: "figures/count_temp_plot.png",
#     width = 6, #in inches, think about the final size in a word doc page etc.
#     height = 5,
#     dpi = 600 #or 300 if you are worried about file size
# )