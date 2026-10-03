# GLM count-data tutorial: Poisson and negative binomial regression
#
# This tutorial uses data on male satellite counts associated with female
# horseshoe crabs. It uses female width as a continuous predictor and female
# color as a categorical predictor.
#
# Watch these background videos before working through the tutorial.
# Rewatch them as needed while working through the tutorial.
# GLM Part 1: A New Perspective (4 min)
# https://youtu.be/6Nfv0Xr44y8
# GLM Part 2: Count Regression (7 min)
# https://youtu.be/i62gffPrZYA
# Poisson Regression in R (25 min)
# https://youtu.be/FfBnX5dfxXw
#
# Here is a book you can check out for more information about GLMs.
# (The CPP Library has a copy.)
# Zuur, Hilbe, and Ieno (2013), A Beginner's Guide to GLM and GLMM with R
# https://csu-cpp.primo.exlibrisgroup.com/permalink/01CALS_PUP/1fn8cvg/alma991008774460302915
# Book website:
# https://www.highstat.com/index.php/books2?catid=18&id=21&view=article
#
# This tutorial focuses on Poisson and negative binomial models for count data.

# Load packages ------------------------------------------------------------

library(tidyverse)
library(glm2) # for the crabs data set
library(DHARMa) # for GLM diagnostics
library(ggeffects) # for response-scale predictions and plots
library(emmeans) # for estimated marginal means and contrasts
library(AICcmodavg) # for AICc model comparison
library(broom) # for tidying model outputs
library(skimr) # for data summaries and checks
library(cowplot) # for combining plots

# Load and format the data -------------------------------------------------

# The data set contains 173 female horseshoe crabs.
# Satellites = number of male partners
# Width = female width in centimeters
# Dark = whether the female has dark coloring ("yes" or "no")
#
# The data are derived from Agresti (2007, Table 3.2).
# Documentation:
# https://www.rdocumentation.org/packages/glm2/versions/1.2.1/topics/crabs

dat_crab <- glm2::crabs |>
  select(Satellites, Width, Dark) |>
  transmute(
    count_males = Satellites,
    fem_color = recode(
      as.character(Dark),
      "no" = "light",
      "yes" = "dark"
    ) |>
      factor(levels = c("light", "dark")),
    fem_width_cm = Width
  )

glimpse(dat_crab)

# The light-color category is the reference level. In a model containing
# fem_color, the coefficient for fem_colordark compares dark females with
# light females.

# Basic variable summaries -------------------------------------------------
dat_crab |>
  skim()

# Compare simple summaries among female-color categories.
dat_crab |>
  group_by(fem_color) |>
  summarise(
    n = n(),
    mean = mean(count_males),
    median = median(count_males),
    variance = var(count_males),
    zeros = sum(count_males == 0),
    .groups = "drop"
  )

# Explore the response and explanatory variables ----------------------------

# Count responses are non-negative integers. Integer-aligned bins in your
# histogram make the response distribution easy to see.
dat_crab |>
  ggplot(aes(x = count_males)) +
  geom_histogram(
    binwidth = 1,
    boundary = -0.5, # align bins with integers
    colour = "white" # make bar outlines white
  ) +
  facet_grid(rows = vars(fem_color)) +
  labs(
    x = "Number of male crabs"
  )

# Note whether the response is positively-skewed, whether there are many zeros,
# and whether the distribution appears different between color categories.

# Plot individual observations among the color categories.
dat_crab |>
  ggplot(aes(x = fem_color, y = count_males, colour = fem_color)) +
  geom_jitter(width = 0.08, height = 0, alpha = 0.6) +
  scale_colour_manual(values = c(light = "grey60", dark = "black")) +
  labs(
    x = "Female color",
    y = "Number of male satellite crabs",
    colour = "Female color"
  ) +
  theme_classic()

# Explore the relationship between female width and the number of male
# satellites. Color is used to see whether the relationship may differ between
# light and dark females.
dat_crab |>
  ggplot(aes(
    x = fem_width_cm,
    y = count_males,
    colour = fem_color
  )) +
  geom_point(alpha = 0.6) +
  scale_colour_manual(values = c(light = "grey60", dark = "black")) +
  labs(
    x = "Female width (cm)",
    y = "Number of male satellite crabs",
    colour = "Female color"
  ) +
  theme_classic()

# View summarized distributions as boxplots.
dat_crab |>
  ggplot(aes(
    x = fem_color,
    y = count_males,
    colour = fem_color
  )) +
  geom_boxplot() +
  labs(
    y = "Number of male satellite crabs per female"
  ) +
  theme_classic()

# Questions to ask during exploration:
# - Are there many zeros?
# - Is the response strongly positively-skewed?
# - Does the response appear more variable for wider females?
# - Does the response differ between light and dark females?
# - Does the width relationship appear different between color categories?


# What if you fit an ordinary linear regression to the count response? The
# plot below shows a linear fit to the count data. This is not a recommended
# approach for count data, but it is shown here for comparison with the GLM
# approach.
dat_crab |>
  ggplot(aes(
    x = fem_width_cm,
    y = count_males
  )) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(
    x = "Female width (cm)",
    y = "Number of male satellite crabs"
  ) +
  theme_classic()
# The data are not equally distributed around the line. The line extends into
# negative values for small widths, which is not biologically possible. The
# linear model does not account for the non-negative integer nature of the
# response; you cannot have a negative number of male crabs.

# Connect the GLM (Generalized Linear Model) components to this analysis -----

# Every GLM combines:
# 1. A response distribution describing variation around the expected mean.
# 2. A predictor structure describing how explanatory variables affect that
#    expected mean.
# 3. A link function connecting the predictors to the expected response.
#
# For this tutorial:
# - the response distribution is Poisson or negative binomial; negative binomial
# is used if the amount of variability is too much for the Poisson distribution
# to handle ("overdispersion");
# - the predictors can include female width, female color, and their
#   interaction (in this context, separate curves for light and dark females);
# - the link function is the log link because it is used with Poisson and
#   negative binomial count models.
#
# The log link keeps predicted counts positive. Predictors combine linearly on
# the log expected-count scale, while plots and predictions are later returned
# ("back-transformed") on the original count scale.
#
# Poisson and negative binomial models are for integer counts.

# Fit a Poisson model with female width -------------------------------------

# Research question: Does female width predict the number of male crabs?

# Start with one continuous predictor. The Poisson model assumes that the
# conditional variance is approximately equal to the conditional mean.
# (i.e., as the predicted count of males increases, the spread also increases in
# a predictable way)

# Here is the base plot for this question (we are ignoring female color for now).
dat_crab |>
  ggplot(aes(
    x = fem_width_cm,
    y = count_males
  )) +
  geom_point(alpha = 0.6) +
  labs(
    x = "Female width (cm)",
    y = "Number of male satellite crabs"
  ) +
  theme_classic()

# Fit the Poisson model with female width as the only predictor.
mod_pois_count_width <- glm(
  count_males ~ fem_width_cm,
  family = poisson(link = "log"),
  data = dat_crab
)

summary(mod_pois_count_width)

# Review the model coefficients and their uncertainty.
broom::tidy(
  mod_pois_count_width,
  conf.int = TRUE,
  conf.level = 0.95,
  exponentiate = TRUE # converts log-scale coefficients to the response scale
)

# An exponentiated width coefficient is the multiplicative change in expected
# male-satellite count associated with a 1-cm increase in female width.
# For example, 1.18 would mean 1.18 times the expected count, or an 18% increase,
# for each additional centimeter of width.
# Here the 95% confidence interval is also on the response scale. So we can say
# we are 95% confident the true rate of increase is 13% to 23% per additional
# centimeter of width.

# Model diagnostics for the Poisson width model with DHARMa ------------------

# Watch these videos before or during the diagnostic sections:
# GLM Part 4: Overdispersion
# https://youtu.be/xJm6eN5ZDzk
# GLM Part 5: Diagnostics
# https://youtu.be/SnZDysPBaW4
# DHARMa vignettes:
# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html

# Note: A significant p-value in a DHARMa test does not automatically mean
# that your model is "wrong" or "bad." It indicates that the model does not
# adequately describe the data in a particular way. In some cases, you should
# try a different model; for example, if a Poisson model is overdispersed, you
# might try a negative binomial model. If one or two tests are marginally
# significant, consider whether the model is biologically sensible, whether it
# has otherwise reasonable diagnostics, and whether it answers your research
# question. You can also describe your data and study, then copy and paste the
# model summary and DHARMa results into an AI tool to discuss options for
# investigating the diagnostic issue further.

set.seed(38)
sim_pois_count_width <- DHARMa::simulateResiduals(
  fittedModel = mod_pois_count_width,
  n = 3000
)

plot(sim_pois_count_width)

# This test looks for overall distribution mismatches by checking whether the
# observed residuals fit the overall pattern expected by the model.
DHARMa::testUniformity(sim_pois_count_width)
# A significant p-value indicates that the model may not adequately describe
# the data in some way.

# This test looks for incorrect variance by checking if the data spreads out
# significantly more (overdispersion) or less (underdispersion) than what your
# model predicts.
DHARMa::testDispersion(sim_pois_count_width)

# This test looks for an excess of zeros by checking if your data contains more
# zero values than would be expected based on the model you are checking
DHARMa::testZeroInflation(sim_pois_count_width)


# A quick Pearson-dispersion summary:
pearson_dispersion_pois_width <- sum(
  residuals(mod_pois_count_width, type = "pearson")^2
) / df.residual(mod_pois_count_width)

pearson_dispersion_pois_width

# A value substantially greater than 1 can indicate overdispersion, but it
# is not a final decision rule. Also consider whether the mean structure is
# missing a nonlinear pattern, an important predictor, or dependence among
# observations.

# Fit a negative binomial model with female width -----------------------------

# If the Poisson model is overdispersed, fit a negative binomial model. The
# negative binomial allows more variation (spread in the response) than the
# Poisson distribution.
mod_nb_count_width <- MASS::glm.nb(
  count_males ~ fem_width_cm,
  link = "log",
  data = dat_crab
)

summary(mod_nb_count_width)

broom::tidy(
  mod_nb_count_width,
  conf.int = TRUE,
  conf.level = 0.95,
  exponentiate = TRUE # converts log-scale coefficients to the response scale
)

# An exponentiated width coefficient is the multiplicative change in expected
# male-satellite count associated with a 1-cm increase in female width.
# For example, 1.21 would mean 1.21 times the expected count, or a 21% increase,
# for each additional centimeter of width.
# Here the 95% confidence interval is also on the response scale. So we can say
# we are 95% confident the true rate of increase is 11% to 33% per additional
# centimeter of width.



# Check diagnostics for the negative binomial width model.
set.seed(38)
sim_nb_count_width <- DHARMa::simulateResiduals(
  fittedModel = mod_nb_count_width,
  n = 2000
)

plot(sim_nb_count_width)

# This test looks for overall distribution mismatches by checking whether the
# observed residuals fit the overall pattern expected by the model.
DHARMa::testUniformity(sim_nb_count_width)
# A significant p-value indicates that the model may not adequately describe
# the data in some way.

# This test looks for incorrect variance by checking if the data spreads out
# significantly more (overdispersion) or less (underdispersion) than what your
# model predicts.
DHARMa::testDispersion(sim_nb_count_width)

# This test looks for an excess of zeros by checking if your data contains more
# zero values than would be expected based on the model you are checking
DHARMa::testZeroInflation(sim_nb_count_width)

# Compare Poisson and negative binomial width models -------------------------

AICc(mod_pois_count_width)
AICc(mod_nb_count_width)

# AICc is a relative comparison of the models; lower AICc indicates more support.
# The top model should also have reasonable DHARMa diagnostics and make biological sense.

# Choose the supported width model for plotting.
# Change this to mod_pois_count_width if the Poisson model is supported.
mod_final_count_width <- mod_nb_count_width

# Plot the fitted relationship for the width model ---------------------------

# For plotting the model, use ggeffects to estimate the expected counts
# (regression line or curve) and 95% confidence intervals.
pred_count_width <- ggeffects::predict_response(
  mod_final_count_width,
  terms = "fem_width_cm [all]"
)

pred_count_width_df <- as_tibble(pred_count_width)

ggplot(pred_count_width_df, aes(x = x, y = predicted)) +
  # Add a ribbon for the 95% confidence interval
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.2
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(
    data = dat_crab,
    aes(x = fem_width_cm, y = count_males),
    inherit.aes = FALSE,
    alpha = 0.55
  ) +
  # Expand the y-axis to include zero
  expand_limits(y = 0) +
  labs(
    x = "Female width (cm)",
    y = "Count of male satellite crabs"
  ) +
  theme_minimal()

# Note: Predictions (the curved line) and 95% confidence intervals are
# back-transformed to the original count scale.

# The confidence interval is often asymmetric, i.e., the upper bound is further
# from the predicted mean than the lower bound. This is a common feature of
# log-linked (Poisson and negative binomial) models.


# Fit a Poisson model for female color --------------------------------------

# Research question: Does female color predict the number of male crabs?
# (Ignore width for this model.)

# Now use female color as the only predictor. Because fem_color is a factor,
# the model compares the dark category with the light reference category.

# Here is the base plot for this question (we are ignoring female width for now).
dat_crab |>
  ggplot(aes(
    x = fem_color,
    y = count_males,
    colour = fem_color
  )) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.5, size = 2) +
  scale_colour_manual(values = c(light = "grey60", dark = "black")) +
  labs(
    x = "Female color",
    y = "Number of male satellite crabs",
    colour = "Female color"
  ) +
  theme_classic()

# Fit the Poisson model with female color as the only predictor
mod_pois_count_color <- glm(
  count_males ~ fem_color,
  family = poisson(link = "log"),
  data = dat_crab
)

summary(mod_pois_count_color)

broom::tidy(
  mod_pois_count_color,
  conf.int = TRUE,
  conf.level = 0.95,
  exponentiate = TRUE
)
# An exponentiated fem_colordark coefficient is the multiplicative change in
# expected male-satellite count for dark females compared with light females. In
# this case 0.64 means that dark females are expected to have 0.64 times the
# number of male satellites as light females, or a 36% decrease in expected
# count. The 95% confidence interval is also on the response scale, so we can
# say we are 95% confident that the true ratio of expected counts for dark
# versus light females is between 0.53 and 0.78.


# Model diagnostics for the Poisson color model.
set.seed(38)
sim_pois_count_color <- DHARMa::simulateResiduals(
  fittedModel = mod_pois_count_color,
  n = 2000
)

plot(sim_pois_count_color)

# This test looks for overall distribution mismatches by checking whether the
# model's residuals are uniformly distributed.
DHARMa::testUniformity(sim_pois_count_color)

# This test looks for incorrect variance by checking whether the data spreads
# out too much or too little.
DHARMa::testDispersion(sim_pois_count_color)

# This test looks for more zeros than the model would predict.
DHARMa::testZeroInflation(sim_pois_count_color)

# Fit a negative binomial model with female color ----------------------------

# Because the Poisson model is overdispersed, fit a negative binomial model
# instead. The negative binomial allows more variation (spread in y) than the
# Poisson distribution.

mod_nb_count_color <- MASS::glm.nb(
  count_males ~ fem_color,
  link = "log",
  data = dat_crab
)

summary(mod_nb_count_color)

broom::tidy(
  mod_nb_count_color,
  conf.int = TRUE,
  conf.level = 0.95,
  exponentiate = TRUE
)

# An exponentiated fem_colordark coefficient is the multiplicative change in
# expected male-satellite count for dark females compared with light females. In
# this case, the coefficient is 0.64, meaning that dark females are expected to
# have 0.64 times the number of male satellites as light females, or a 36%
# decrease in expected count. The 95% confidence interval is also on the
# response scale, so we can say we are 95% confident that the true ratio of
# expected counts for dark versus light females is between 0.43 and 0.95

# Diagnose the negative binomial color model.
set.seed(38)
sim_nb_count_color <- DHARMa::simulateResiduals(
  fittedModel = mod_nb_count_color,
  n = 2000
)

plot(sim_nb_count_color)

# This test looks for overall distribution mismatches by checking whether the
# model's residuals are uniformly distributed.
DHARMa::testUniformity(sim_nb_count_color)

# This test looks for incorrect variance by checking whether the data spreads
# out too much or too little.
DHARMa::testDispersion(sim_nb_count_color)

# This test looks for more zeros than the model would predict.
DHARMa::testZeroInflation(sim_nb_count_color)

# Interpreting the DHARMa uniformity test above -------------------------------

# The significant uniformity test suggests that the model residuals are not
# perfectly consistent with the residual pattern expected under the fitted
# negative-binomial model.
#
# This is an omnibus test: it indicates that some aspect of model fit may be
# imperfect, but it does not identify the source or indicate how biologically
# important the problem is.
#
# The nonsignificant dispersion and zero-inflation tests provide no evidence
# for those specific problems. We should also inspect the DHARMa plots and
# residuals by female color. If there are no strong residual patterns or
# influential observations, we can retain this model as a useful candidate,
# while acknowledging that the color-only model shows some lack of fit.
#
# For this color-only model, the significant uniformity test could indicate that
# female color alone does not capture all of the variation in male satellite
# counts. For example, female width may explain important variation that is
# omitted from this model. The width-only or width-by-color model may therefore
# provide better overall fit.
#
# We will also compare this model with models that include female width and
# the width-by-color interaction. A better-fitting model may indicate that
# female width explains variation not captured by female color alone.


# Compare Poisson and negative binomial color models using AICc ---------------

AICc(mod_pois_count_color)
AICc(mod_nb_count_color)
# The negative binomial model has a lower AICc and is therefore preferred over
# the Poisson model. The negative binomial model also has more acceptable 
# DHARMa diagnostics.


# Choose the supported color model for the estimated marginal means.
# Change this to mod_pois_count_color if the Poisson model is supported.
mod_final_count_color <- mod_nb_count_color

# Estimate and plot expected counts for each color category -------------------

# emmeans estimates the expected count for each color category on the response
# scale. The confidence intervals cannot extend below zero and are often
# asymmetric.
emm_count_color <- emmeans::emmeans(
  mod_final_count_color,
  ~ fem_color,
  type = "response",
  CIs = TRUE
)

# The response-scale estimate is the expected count for each color category.

# Plot the estimated means and 95% CIs ----------------------------------------

plot_count_color_means <- dat_crab |>
  ggplot(
    aes(
      x = fem_color,
      y = count_males,
      color = fem_color
    )
  ) +
  geom_jitter(
    width = 0.08,
    height = 0,
    alpha = 0.45,
    size = 2
  ) +
  scale_colour_manual(values = c(light = "grey60", dark = "black")) +  
  geom_crossbar(
    data = as_tibble(emm_count_color),
    aes(
      x = fem_color,
      y = response,
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.7
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.10))
  ) +
  labs(
    x = "Female color",
    y = "Number of male satellite crabs"
  ) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

plot_count_color_means

# The confidence intervals are often asymmetric, i.e., the upper bound is
# further from the predicted mean than the lower bound. This is a common feature
# of log-linked (Poisson and negative binomial) models.

# Pairwise contrast between dark and light females ----------------------------
color_contrast_response <- emmeans::contrast(
  emm_count_color,
  method = list("Dark / Light" = c(-1, 1)),
  ratios = TRUE
) |>
  summary(
    infer = c(TRUE, TRUE),
    type = "response"
  )

color_contrast_response

# Plot the dark-to-light response ratio and its 95% CI ------------------------

plot_count_color_contrast <- as_tibble(color_contrast_response) |>
  ggplot(
    aes(
      y = contrast,
      x = ratio,
      label = format.pval(
        p.value,
        digits = 2,
        eps = 0.001
      )
    )
  ) +
  geom_pointrange(
    aes(
      xmin = asymp.LCL,
      xmax = asymp.UCL
    )
  ) +
  geom_vline(
    xintercept = 1,
    linetype = 2
  ) +
  scale_x_continuous(
    position = "top",
    expand = expansion(mult = c(0.10, 0.15))
  ) +
  labs(
    x = "Response ratio",
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 11),
    plot.subtitle = element_text(face = "italic")
  ) +
  # Plot the p-value above the point estimate and confidence interval.
  geom_text(
    vjust = 0,
    nudge_y = 0.2,
    size = 3
  )

plot_count_color_contrast

# The pairwise contrast is commonly presented as a count ratio because of the
# log link function.

# Combine the contrast and estimated-means plots ------------------------------

final_count_color_plot <- cowplot::plot_grid(
  plot_count_color_contrast,
  plot_count_color_means,
  nrow = 2,
  align = "v",
  axis = "rl",
  rel_heights = c(0.35, 1)
)

final_count_color_plot



# Fit the "full" model with both width and color -----------------------------

# Research question: Do female width and color together predict the number of
# male crabs?

# We will fit the fem_width_cm * fem_color interaction model, which creates
# separate lines (curves) with different intercepts and slopes for light and
# dark females.

mod_pois_count_width_color <- glm(
  count_males ~ fem_width_cm * fem_color,
  family = poisson(link = "log"),
  data = dat_crab
)

summary(mod_pois_count_width_color)

broom::tidy(
  mod_pois_count_width_color,
  conf.int = TRUE,
  conf.level = 0.95,
  exponentiate = TRUE
)
# With this full model where they have separate slopes and intercepts, it
# becomes harder to interpret the coefficients by themselves. Plotting and
# interpreting the separate curves from the plots is typically more useful for
# interpretation than trying to interpret the intercepts and slope coefficients.


# Diagnose the Poisson interaction model.
set.seed(38)
sim_pois_count_width_color <- DHARMa::simulateResiduals(
  fittedModel = mod_pois_count_width_color,
  n = 2000
)

plot(sim_pois_count_width_color)
DHARMa::testUniformity(sim_pois_count_width_color)
DHARMa::testDispersion(sim_pois_count_width_color)
DHARMa::testZeroInflation(sim_pois_count_width_color)

# These DHARMa tests show that the model does not capture the variability in
# the data well. The Poisson model is overdispersed, and the uniformity test is
# significant. The negative binomial model may be a better choice for this
# interaction model.


# Fit a negative binomial full model ------------------------------------------

mod_nb_count_width_color <- MASS::glm.nb(
  count_males ~ fem_width_cm * fem_color,
  link = "log",
  data = dat_crab
)

summary(mod_nb_count_width_color)

broom::tidy(
  mod_nb_count_width_color,
  conf.int = TRUE,
  conf.level = 0.95,
  exponentiate = TRUE
)
# Interpreting the coefficients in an interaction model is more complex than
# interpreting coefficients in a simple model. Plotting and interpreting the
# separate curves for light and dark females are typically more useful than
# trying to interpret the intercept and slope coefficients individually.


# Diagnose the negative binomial interaction model.
set.seed(38)
sim_nb_count_width_color <- DHARMa::simulateResiduals(
  fittedModel = mod_nb_count_width_color,
  n = 2000
)

plot(sim_nb_count_width_color)

# This test looks for overall distribution mismatches by checking whether the
# model's residuals are uniformly distributed.
DHARMa::testUniformity(sim_nb_count_width_color)

# This test looks for whether the model captures the expected variance by
# checking whether the data spreads out too much or too little.
DHARMa::testDispersion(sim_nb_count_width_color)

# This test looks for more zeros than the model would predict.
DHARMa::testZeroInflation(sim_nb_count_width_color)

# Compare the full Poisson and negative binomial models with AICc -------------

AICc(mod_pois_count_width_color)
AICc(mod_nb_count_width_color)
# The negative binomial model has a lower AICc and is therefore preferred over
# the Poisson model. The negative binomial model also has more acceptable DHARMa
# diagnostics.


# Choose the supported interaction model for plotting.
# Change this to mod_pois_count_width_color if the Poisson model is supported.
mod_final_count_width_color <- mod_nb_count_width_color

# Plot separate response-scale curves for light and dark females -------------

# For plotting the model, use ggeffects to estimate the expected counts
# (curves) and 95% confidence intervals for light and dark females separately.
pred_count_width_color <- ggeffects::predict_response(
  mod_final_count_width_color,
  terms = c("fem_width_cm [all]", "fem_color")
)

pred_count_width_color_df <- as_tibble(pred_count_width_color)

ggplot(
  pred_count_width_color_df,
  aes(
    x = x,
    y = predicted,
    colour = group,
    fill = group,
    group = group
  )
) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(
    data = dat_crab,
    aes(
      x = fem_width_cm,
      y = count_males,
      colour = fem_color
    ),
    inherit.aes = FALSE,
    alpha = 0.5
  ) +
  scale_colour_manual(values = c(light = "grey60", dark = "black")) +
  scale_fill_manual(values = c(light = "grey60", dark = "black")) +
  expand_limits(y = 0) +
  labs(
    x = "Female width (cm)",
    y = "Count of male satellite crabs",
    colour = "Female color"
  ) +
  theme_minimal()
# The confidence interval is often asymmetric, i.e., the upper bound is further
# from the predicted mean than the lower bound. This is a common feature of
# log-linked (Poisson and negative binomial) models.


# Interpret the separate curves rather than the interaction coefficient alone:
# - Do the curves have different slopes?
# - Are the predicted counts meaningfully different across the observed width
#   range?
# - Do the confidence bands provide evidence of a biologically important
#   difference?
# - Are the predictions within the observed range of female widths?

# Compare the three predictor structures with AICc --------------------------

# After reviewing the Poisson-versus-negative-binomial comparisons above,
# choose the response distribution to use for comparing predictor structures.
#
# The table below uses negative binomial models as an example. If the Poisson
# model is adequately dispersed and otherwise supported, replace these three
# objects with mod_pois_count_width, mod_pois_count_color, and
# mod_pois_count_width_color.
models_predictor_structure <- list(
  Width = mod_nb_count_width,
  Color = mod_nb_count_color,
  Width_by_color = mod_nb_count_width_color
)

aictab_predictor_structure <- AICcmodavg::aictab(
  cand.set = models_predictor_structure,
  modnames = names(models_predictor_structure),
  second.ord = TRUE,
  sort = TRUE
)

aictab_predictor_structure
# The model including female width alone had the lowest AICc and received the
# most support, but the width-by-color interaction model was also plausible
# because it differed by only 1.27 AICc units. These models produced very
# similar predictions across the observed range of female widths, and the
# interaction model suggested only a small difference between light and dark
# females. Thus, the data provide limited evidence that the width-count
# relationship differs by female color. The simpler width-only model would be a
# reasonable final choice, while recognizing that a relatively small
# color-specific effect cannot be ruled out without more data. The color-only
# model (predicting a mean for each color) had substantially less support,
# indicating that female color alone did not explain the variation in male
# satellite counts as well as models including female width.


# AICc compares relative support among these three candidate explanations
# (biological hypotheses). The selected model should also have acceptable DHARMa
# diagnostics and a biologically sensible interpretation. A more complex
# interaction model is not automatically preferable just because it fits the
# observed data better.
