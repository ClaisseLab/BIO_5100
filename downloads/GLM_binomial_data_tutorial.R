# GLM binomial-data tutorial: logistic regression for presence/absence data
#
# This tutorial uses presence/absence data for the mite PHTH from the
# vegan mite data set. It uses topography as a categorical predictor and water
# content as a continuous predictor.
#
# The main analysis fits and compares three binomial GLMs:
# 1. mite presence/absence as a function of topography only
# 2. mite presence/absence as a function of water content only
# 3. mite presence/absence as a function of topography, water content, and their
# interaction (i.e., separate water-content slopes and intercepts for each
# topography level)
#
# A separate example at the end of the script introduces grouped binomial data
# using seed germination observations recorded as successes out of trials.
#
# Watch these background videos before working through the tutorial ----------
# Rewatch them as needed while working through the tutorial.
# 1) StatQuest: Logistic Regression (another name for Binomial GLM)
# https://youtu.be/yIYKR4sgzI8
# 2) GLM Part 3: Logistic Regression
# https://youtu.be/Gvf815Mz-mU
# 3) Binomial GLM (Logistic Regression) using tidyverse in R:
# https://youtu.be/E7J3M1oYVlc
#
# This tutorial focuses on binomial generalized linear models with the logit
# link, also called logistic regression. For individual observations, the
# response is coded as 0 or 1 and is sometimes called Bernoulli data. For
# grouped observations, the response can be entered as successes and failures
# using cbind(successes, failures).

# Load packages --------------------------------------------------------------

library(tidyverse)
library(vegan) # for the mite and mite.env data sets
library(DHARMa) # for simulation-based GLM diagnostics
library(ggeffects) # for response-scale predictions and plots
library(emmeans) # for estimated marginal means and contrasts
library(broom) # for tidying model outputs
library(AICcmodavg) # for AICc model-comparison tables
library(skimr) # for compact data summaries
library(cowplot) # for combining plots
library(agridat) # for the grouped-binomial germination example

# Load and format the mite presence/absence data -----------------------------

# The vegan package contains two linked data objects:
# - mite: species observations for 70 soil cores
# - mite.env: environmental variables for the same soil cores
#
# We select only the PHTH column and immediately convert it to a binary
# presence/absence response. After this step, mite_present is the only PHTH
# variable used in the rest of the analysis.
data("mite", package = "vegan")
data("mite.env", package = "vegan")

dat_mite_response <- mite |>
  as_tibble() |>
  select(PHTH) |>
  transmute(
    mite_present = as.integer(PHTH > 0)
  )

# Select only the two predictors used in this tutorial from the environmental
# data. Topography is categorical, while water content is continuous and
# measured in g/L (see the vegan mite data documentation).
# https://vegandevs.github.io/vegan/reference/mite.html
dat_mite <- dat_mite_response |>
  bind_cols(
    mite.env |>
      as_tibble() |>
      select(
        topography = Topo,
        water_content = WatrCont
      )
  ) |>
  mutate(
    topography = factor(topography)
  )

# Check the final analysis data frame.
dat_mite
glimpse(dat_mite)

# Confirm that the response contains only 0 and 1 values.
sort(unique(dat_mite$mite_present))

# General data summaries -----------------------------------------------------
# Summarize the number and proportion of presence and absence observations.
# For a 0/1 response, the mean equals the observed proportion present.
# Here we summarize the response without considering the predictors.
dat_mite |>
  count(mite_present) |>
  mutate(
    response_label = if_else(
      mite_present == 1,
      "Present",
      "Absent"
    ),
    proportion = n / sum(n)
  )

# Summarize each variable separately.
skimr::skim(dat_mite)

# Check whether each topography level includes both presence and absence
# observations. A level containing only one response value can result in
# unstable estimates or complete separation in a binomial GLM.
dat_mite |>
  count(topography, mite_present) |>
  complete(
    topography,
    mite_present = c(0, 1),
    fill = list(n = 0)
  ) |>
  arrange(topography, mite_present)

# Explore the response and explanatory variables -----------------------------

# Explore mite occurrence by topography. The jittered observations show the
# raw binary response, while the larger points show the observed proportion
# present within each topography level.
topography_summary <- dat_mite |>
  group_by(topography) |>
  summarise(
    n_cores = n(),
    n_present = sum(mite_present),
    proportion_present = mean(mite_present),
    .groups = "drop"
  )

topography_summary

plot_mite_topo <- ggplot() +
  geom_jitter(
    data = dat_mite,
    aes(
      x = topography,
      y = mite_present
    ),
    width = 0.10, # to spread overlapping points horizontally for visibility
    height = 0,
    alpha = 0.50,
    colour = "grey35",
    size = 2
  ) +
  geom_point(
    data = topography_summary,
    aes(
      x = topography,
      y = proportion_present
    ),
    size = 3,
    colour = "darkorange3"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1),
    expand = c(0, 0.05)
  ) +
  labs(
    x = "Topography",
    y = "Mite occurrence or observed proportion present"
  ) +
  theme_classic()

plot_mite_topo

# Explore mite occurrence against the continuous predictor. Because the
# response is binary, points occur at 0 and 1 rather than forming a
# continuous response cloud. Transparency helps reveal overlapping points.
plot_mite_water <- ggplot(
  dat_mite,
  aes(
    x = water_content,
    y = mite_present
  )
) +
  geom_point(
    alpha = 0.50, # to improve visibility of overlapping points
    colour = "steelblue4",
    size = 2
  ) +
  scale_y_continuous(
    breaks = c(0, 1),
    limits = c(-0.08, 1.08),
    labels = c("Absent", "Present")
  ) +
  labs(
    x = "Water content of substrate (g/L)",
    y = "Mite occurrence"
  ) +
  theme_classic()

plot_mite_water

# Explore whether water content differs among topography levels. This is
# useful context because predictors can be associated with one another in
# observational ecological data.
plot_water_by_topo <- ggplot(
  dat_mite,
  aes(
    x = topography,
    y = water_content
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    fill = "grey85"
  ) +
  geom_jitter(
    width = 0.10,
    alpha = 0.70,
    colour = "steelblue4"
  ) +
  labs(
    x = "Topography",
    y = "Water content of substrate (g/L)"
  ) +
  theme_classic()

plot_water_by_topo

# Explore mite occurrence against the continuous predictor by topography.
plot_mite_water_topo <- ggplot(
  dat_mite,
  aes(
    x = water_content,
    y = mite_present
  )
) +
  geom_point(
    alpha = 0.50,
    colour = "steelblue4",
    size = 2
  ) +
  scale_y_continuous(
    breaks = c(0, 1),
    limits = c(-0.08, 1.08),
    labels = c("Absent", "Present")
  ) +
  labs(
    x = "Water content of substrate (g/L)",
    y = "Mite occurrence"
  ) +
  facet_grid(rows = vars(topography)) +
  theme_classic()

plot_mite_water_topo

# Questions to consider during exploration:
# - Are mite presence and absence observations both sufficiently represented?
# - Does the observed proportion present differ between topography levels?
# - Does mite occurrence appear to change across water content?
# - Does the water-content relationship appear different among topography
#   levels?

# Connect the GLM components to this analysis --------------------------------

# Every GLM combines:
# 1. A response distribution describing variation around the expected mean.
# 2. A predictor structure describing how explanatory variables affect the
#    expected mean.
# 3. A link function connecting the predictors to the expected response.
#
# For this tutorial:
# - the response distribution is binomial; for individual 0/1 observations,
#   this is also called a Bernoulli response;
# - the predictors are topography, water content, or both with an interaction
#   that allows a separate water-content relationship for each topography;
# - the link is the logit, a commonly used link for binomial data.
#
# The binomial distribution is appropriate when each observation represents a
# success or failure, such as presence/absence or yes/no. The model estimates
# the probability of presence (or success), p.

# A binomial GLM with a categorical predictor --------------------------------

# Research question:
# Does the probability of mite presence differ among topography levels?

mod_topo <- glm(
  mite_present ~ topography,
  family = binomial(link = "logit"),
  data = dat_mite
)

summary(mod_topo)

# Probabilities of mite presence by topography with 95% CIs ------------------

# Use emmeans to summarize the model on the response (probability) scale.
# Estimated marginal means provide model-estimated probabilities of PHTH
# presence for each topography level. type = "response" back-transforms the
# estimates from the logit scale to the probability scale.
emm_topo <- emmeans::emmeans(
  mod_topo,
  ~ topography,
  type = "response"
)

emm_topo

# The estimated mean probability of PHTH presence is 0.18 (18%) for
# Blanket topography and 0.62 (62%) for Hummock topography.
#
# The 95% confidence intervals describe uncertainty in these estimated
# probabilities. Based on this model, the probability of PHTH presence is
# plausibly between approximately 9% and 32% in Blanket topography, and
# between approximately 42% and 78% in Hummock topography.
#
# These are confidence intervals for the estimated average probability, not
# intervals containing the outcomes of individual soil cores. Each individual
# core still has a binary outcome: PHTH is either present or absent.

# Pairwise comparisons among topography levels -------------------------------

# We have already estimated the probability of PHTH presence for each
# topography level. We can now directly compare those two estimated
# probabilities:
#
# - Blanket: 18% (95% CI 9% to 32%)
# - Hummock: 62% (95% CI 42% to 78%)
#
# The separated confidence intervals suggest a difference. Use a pairwise
# contrast for the formal comparison, rather than judging significance from
# whether the individual confidence intervals overlap.
#
# For a binomial GLM with a logit link, the default response-scale contrast
# is expressed as a ratio of odds. The contrast provides an odds ratio,
# its 95% confidence interval, and a p-value for testing equal odds.
# With a logit link, equal odds also imply
# equal probabilities.

# Odds describe how likely an event is to occur relative to not occurring:
#
# odds = probability of presence / probability of absence
# odds = p / (1 - p)
#
# For example:
# - A probability of 0.50 gives odds of 0.50 / 0.50 = 1.
# - A probability of 0.75 gives odds of 0.75 / 0.25 = 3.
# - A probability of 0.20 gives odds of 0.20 / 0.80 = 0.25.
#
# Odds of 1 mean that presence and absence are equally likely. Odds greater
# than 1 indicate that presence is more likely than absence, while odds less
# than 1 indicate that presence is less likely than absence.
#

pairs_topo <- emmeans::contrast(
  emm_topo,
  method = "revpairwise", # see note below
  adjust = "none",
  type = "response",
  infer = c(TRUE, TRUE) # request the 95% CI and significance test
)
# NOTE: We use method = "revpairwise" to reverse the usual factor-level order.
# Because the levels are ordered Blanket, then Hummock, this produces the
# contrast Hummock / Blanket. The method = "pairwise" argument would produce
# Blanket / Hummock instead.
#
# revpairwise reverses the factor-level order; it does not automatically sort
# levels according to their estimated probabilities. Always check the
# contrast label before interpreting the result.

pairs_topo

# The odds ratio is 7.2 for Hummock relative to Blanket. This means that the
# odds of PHTH presence in Hummock are estimated to be 7.2 times those in
# Blanket topography.
#
# This is not the same as saying that the probability is 7.2 times higher.
# To see the calculation, use the estimated probabilities from emm_topo:
#
# Hummock probability = 0.615
# Blanket probability = 0.182
#
# Convert each probability to odds:
#
# Hummock odds = 0.615 / (1 - 0.615) = 1.60
# Blanket odds = 0.182 / (1 - 0.182) = 0.22
#
# Then calculate the odds ratio:
#
# Using the unrounded odds, Hummock odds / Blanket odds = approximately 7.2.
#
# The 95% confidence interval for the odds ratio is approximately 2.4 to
# 21.6. Because this interval does not include 1, the model provides evidence
# that the odds of PHTH presence differ between the two topography levels.
# The wide interval indicates uncertainty in the precise size of the
# difference.
#
# For an odds ratio, the null value is 1 (equal odds in the two topographies).
# The p-value tests this same null hypothesis. The odds ratio and its
# confidence interval are especially useful because they describe both the
# estimated effect size and its uncertainty.

# Example results statement:
# The estimated probability that a core sample contained PHTH was higher in
# Hummock topography (62%, 95% CI: 42–78%) than in Blanket topography (18%,
# 95% CI: 9–32%; see figure below). The odds of a core sample containing PHTH
# in Hummock were 7.2 times those in Blanket topography (95% CI: 2.4–21.6;
# p = 0.0004).

# Plot the topography model results ------------------------------------------
# Probability of presence for each topography level with 95% CI error bars.
# The individual 0/1 observations for each soil core are jittered horizontally
# for visibility.

plot_topo <- dat_mite |>
  ggplot() +
  geom_jitter(
    aes(
      x = topography,
      y = mite_present
    ),
    width = 0.10,
    height = 0,
    alpha = 0.50,
    colour = "grey40",
    size = 2
  ) +
  geom_errorbar(
    data = as_tibble(emm_topo),
    aes(
      x = topography,
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    width = 0.12,
    linewidth = 0.8,
    colour = "dodgerblue4"
  ) +
  geom_point(
    data = as_tibble(emm_topo),
    aes(
      x = topography,
      y = prob
    ),
    size = 3,
    colour = "dodgerblue4"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1),
    expand = c(0, 0.05)
  ) +
  labs(
    x = "Topography",
    y = "Estimated probability of mite presence"
  ) +
  theme_classic()

plot_topo

# Diagnostic checks for the categorical-only model. DHARMa simulates
# residuals from the fitted model and compares the observed residual pattern
# with the simulated reference distribution.
set.seed(21)
sim_topo <- DHARMa::simulateResiduals(
  fittedModel = mod_topo,
  n = 2000
)

# The diagnostic plot provides an overall visual check of the simulated
# residuals, including their distribution, dispersion, and possible outliers.
plot(sim_topo)

# Test whether the residuals follow the expected uniform distribution.
DHARMa::testUniformity(sim_topo)

# Interpretation: The test was not significant (p = 0.5675), so there is no
# evidence that the residuals deviate from the expected uniform distribution.

# Test for more or less residual variability than expected under the model.
DHARMa::testDispersion(sim_topo)

# Interpretation: The test was not significant (p = 0.916), so there is no
# evidence of overdispersion or underdispersion relative to the model.

# Test whether the number of simulation outliers differs from expectations.
DHARMa::testOutliers(sim_topo)

# Interpretation: No outliers were detected (0 out of 70 observations;
# p = 1.000), so there is no evidence of an unusual number of extreme
# observations.

# With individual binary observations, diagnostics can have limited power to
# detect model problems. Use these results together with the raw data, the
# predictor structure, and the response-scale plots.

# A binomial GLM with a continuous predictor ---------------------------------

# Research question:
# Does the probability of mite presence change with water content?

mod_water <- glm(
  mite_present ~ water_content,
  family = binomial(link = "logit"),
  data = dat_mite
)

summary(mod_water)

# The water-content coefficient describes the change in log-odds associated
# with a 1 g/L increase in water content. Its exponentiated value is the
# multiplicative change in the odds of presence for the same increase.
coef_water <- broom::tidy(
  mod_water,
  conf.int = TRUE,
  conf.level = 0.95,
  exponentiate = TRUE
)

coef_water

# For water_content, the odds ratio is 0.987. This means that the odds of
# PHTH presence are multiplied by 0.987 for each 1 g/L increase in water
# content. Equivalently, the odds decrease by approximately 1.3% per g/L:
#
# 1 - 0.987 = 0.013, or 1.3%
#
# The 95% confidence interval for the odds ratio is 0.979 to 0.993. Thus, the
# estimated decrease in the odds is between approximately 0.7% and 2.1% per
# 1 g/L increase in water content:
#
# 1 - 0.993 = 0.007, or 0.7%
# 1 - 0.979 = 0.021, or 2.1%
#
# Because the confidence interval is entirely below 1, the model provides
# evidence that higher water content is associated with lower odds of PHTH
# presence. The p-value is also small (p = 0.0001), but the odds ratio and its
# confidence interval more directly describe the size and uncertainty of the
# association.
#
# The exponentiated intercept is the estimated odds of PHTH presence when
# water content is zero. This may not be biologically meaningful if zero is
# outside, or near the edge of, the observed range, so interpretation usually
# focuses on the water_content coefficient.

# Example results sentence:
#
# The odds of a core sample containing PHTH decreased by approximately 1.3%
# for each 1 g/L increase in water content (odds ratio = 0.987, 95% CI:
# 0.979–0.993; p = 0.0001).
#
# This is a percentage decrease in odds, rather than a percentage-point
# decrease in the probability of PHTH presence.

# A 1 g/L change may not be the most biologically useful scale for
# interpreting water content. Predicted probabilities let us describe the
# estimated response across the observed range of the predictor.
pred_water <- ggeffects::predict_response(
  mod_water,
  terms = "water_content [all]"
)

pred_water_df <- as_tibble(pred_water)

pred_water_df

# Plot the observed binary data and predicted probabilities with 95%
# confidence intervals on the response scale.
plot_water <- ggplot(
  pred_water_df,
  aes(
    x = x,
    y = predicted
  )
) +
  geom_ribbon(
    aes(
      ymin = conf.low,
      ymax = conf.high
    ),
    alpha = 0.20,
    fill = "steelblue"
  ) +
  geom_line(
    linewidth = 0.9,
    colour = "steelblue4"
  ) +
  geom_point(
    data = dat_mite,
    aes(
      x = water_content,
      y = mite_present
    ),
    inherit.aes = FALSE,
    alpha = 0.55,
    colour = "grey35",
    size = 2
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1),
    expand = c(0, 0.05)
  ) +
  labs(
    x = "Water content of substrate (g/L)",
    y = "Predicted probability of mite presence"
  ) +
  theme_classic()

plot_water

# Diagnostic checks for the continuous-only model.
set.seed(22)
sim_water <- DHARMa::simulateResiduals(
  fittedModel = mod_water,
  n = 2000
)
# The diagnostic plot provides an overall visual check of the simulated
# residuals, including their distribution, dispersion, and possible outliers.
plot(sim_water)

# Test whether the residuals follow the expected uniform distribution.
DHARMa::testUniformity(sim_water)

# Interpretation: The test was not significant (p = 0.652), so there is no
# evidence that the residuals deviate from the expected uniform distribution.

# Test for more or less residual variability than expected under the model.
DHARMa::testDispersion(sim_water)

# Interpretation: The test was not significant (p = 0.733), so there is no
# evidence of overdispersion or underdispersion relative to the model.

# Test whether the number of simulation outliers differs from expectations.
DHARMa::testOutliers(sim_water)

# Interpretation: No outliers were detected (0 out of 70 observations;
# p = 1.000), so there is no evidence of an unusual number of extreme
# observations.

# A binomial GLM with both predictors and their interaction ------------------
# The interaction allows each topography level to have its own water-content
# slope on the log-odds scale, producing separate probability curves.

# Research question:
# Does the relationship between water content and mite presence differ among
# topography levels?

mod_topo_water_full <- glm(
  mite_present ~ topography * water_content,
  family = binomial(link = "logit"),
  data = dat_mite
)

summary(mod_topo_water_full)

# Estimate the water-content trend separately for each topography. The trends
# are initially estimated on the log-odds scale. We exponentiate the trends
# and their confidence limits to express them as odds ratios.
#
# Each odds ratio describes the multiplicative change in the odds of PHTH
# presence associated with a 1 g/L increase in water content within that
# topography. infer = c(TRUE, TRUE) requests both confidence intervals and
# significance tests for the topography-specific slopes.
slopes_odds_ratio <- emmeans::emtrends(
  mod_topo_water_full,
  ~ topography,
  var = "water_content"
) |>
  summary(infer = c(TRUE, TRUE)) |>
  as_tibble() |>
  dplyr::mutate(
    odds_ratio = exp(water_content.trend),
    odds_ratio.LCL = exp(asymp.LCL),
    odds_ratio.UCL = exp(asymp.UCL)
  ) |>
  dplyr::select(
    topography,
    odds_ratio,
    odds_ratio.LCL,
    odds_ratio.UCL,
    p.value
  )

slopes_odds_ratio

# In Blanket samples, the odds ratio was 0.981 (95% CI: 0.968–0.995).
# Therefore, the odds of PHTH presence decreased by approximately 1.9% for
# each 1 g/L increase in water content:
#
# 1 - 0.981 = 0.019, or 1.9%
#
# The confidence interval corresponds to an estimated decrease of approximately
# 0.5% to 3.2% in the odds, and the slope differed from zero based on the
# significance test (p = 0.007).
#
# In Hummock samples, the odds ratio was 0.992 (95% CI: 0.984–1.001).
# This corresponds to an estimated decrease of approximately 0.8% in the odds
# for each 1 g/L increase in water content. However, the confidence interval
# includes 1, and the slope was not clearly different from zero (p = 0.081).
#
# These p-values test whether the water-content slope differs from zero within
# each topography. They do not test whether the slopes differ from one another.
# The topography-by-water-content interaction tests whether the two slopes
# differ.
#
# Example results sentence:
#
# The odds of a core sample containing PHTH decreased by approximately 1.9% for
# each 1 g/L increase in water content in Blanket topography (odds ratio =
# 0.981, 95% CI: 0.968–0.995; p = 0.007). The estimated decline was smaller
# in Hummock topography (odds ratio = 0.992, 95% CI: 0.984–1.001; p = 0.081),
# with the confidence interval including no change in the odds.

# Generate predictions for each topography level across the observed range of
# water content. The predictions and confidence intervals are on the response
# scale because predict_response() back-transforms the logit predictions.
pred_topo_water_full <- ggeffects::predict_response(
  mod_topo_water_full,
  terms = c("water_content [all]", "topography")
)

plot_topo_water_full <- ggplot(
  as_tibble(pred_topo_water_full), # convert to a tibble for plotting
  aes(
    x = x,
    y = predicted,
    colour = group,
    fill = group,
    group = group
  )
) +
  geom_ribbon(
    aes(
      ymin = conf.low,
      ymax = conf.high
    ),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(
    data = dat_mite,
    aes(
      x = water_content,
      y = mite_present,
      colour = topography
    ),
    inherit.aes = FALSE,
    alpha = 0.50,
    size = 2
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1),
    expand = c(0, 0.05)
  ) +
  labs(
    x = "Water content of substrate (g/L)",
    y = "Predicted probability of mite presence",
    colour = "Topography",
    fill = "Topography"
  ) +
  theme_classic()

plot_topo_water_full

# Estimated probabilities at representative water-content values make it
# easier to compare topography levels on the response scale. Here we use the
# 25th, 50th, and 75th percentiles of the observed water-content values.
representative_water_content <- quantile(
  dat_mite$water_content,
  probs = c(0.25, 0.50, 0.75),
  na.rm = TRUE,
  names = FALSE
)

emm_topo_water_full <- emmeans::emmeans(
  mod_topo_water_full,
  ~ topography | water_content,
  at = list(water_content = representative_water_content),
  type = "response"
)

emm_topo_water_full

emm_topo_water_full_df <- as_tibble(emm_topo_water_full)

plot_topo_water_full_emmeans <- ggplot(
  emm_topo_water_full_df,
  aes(
    x = topography,
    y = prob,
    colour = topography
  )
) +
  geom_errorbar(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    width = 0.12,
    linewidth = 0.8
  ) +
  geom_point(size = 3) +
  facet_wrap(
    ~ water_content,
    labeller = label_bquote(
      "Water content = " ~ .(water_content)
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1),
    expand = c(0, 0.05)
  ) +
  labs(
    x = "Topography",
    y = "Estimated probability of mite presence",
    colour = "Topography"
  ) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

plot_topo_water_full_emmeans

# Pairwise comparisons among topography levels at representative values of
# water content. For a logit model, type = "response" gives odds-ratio
# contrasts, while emm_topo_water_full_df contains the estimated probabilities.
pairs_topo_water_full <- emmeans::contrast(
  emm_topo_water_full,
  method = "pairwise",
  adjust = "none"
)

summary(
  pairs_topo_water_full,
  infer = c(TRUE, TRUE),
  type = "response"
)

# Diagnostic checks for the interaction model.
set.seed(23)
sim_topo_water_full <- DHARMa::simulateResiduals(
  fittedModel = mod_topo_water_full,
  n = 2000
)

# Visually check uniformity, residual patterns, and simulation outliers.
plot(sim_topo_water_full)

# Test whether the residuals follow the expected uniform distribution.
DHARMa::testUniformity(sim_topo_water_full)

# Test for more or less residual variability than expected under the model.
DHARMa::testDispersion(sim_topo_water_full)

# Test whether the number of simulation outliers differs from expectations.
DHARMa::testOutliers(sim_topo_water_full)

# Check for residual patterns across water content.
DHARMa::plotResiduals(
  sim_topo_water_full,
  form = dat_mite$water_content
)

# Compare the three binomial models with AICc --------------------------------
# Compare the topography-only, water-content-only, and interaction models.

# AICc compares the relative support for models fitted to the same response
# data. Like AIC, it balances model fit against model complexity, but AICc
# applies an additional correction for smaller sample sizes.
#
# Smaller AICc values indicate stronger relative support among the models being
# compared. The AICc weight gives the relative support for each model compared
# with all models in this candidate set.
#
# AICc is not a test that proves one model is correct. Model comparison should
# also consider the biological question, model diagnostics, estimated effect
# sizes, confidence intervals, and patterns in the raw data.

models_binomial <- list(
  Topography_only = mod_topo,
  Water_content_only = mod_water,
  Topography_by_water_content = mod_topo_water_full
)

# These models must use the same response observations for AICc comparison.
# Create an AICc model-selection table.
aicc_binomial <- AICcmodavg::aictab(
  cand.set = models_binomial,
  modnames = names(models_binomial),
  sort = TRUE
)

aicc_binomial

# The full topo by water model had the strongest relative support among the
# three candidate models (AICc = 65.59; AICc weight = 0.73). The
# water-content-only model also received some support (ΔAICc = 2.01; AICc weight
# = 0.27), so the AICc results do not provide decisive evidence that the
# interaction model is clearly better than the water-content-only model.
#
# The topography-only model had little relative support (ΔAICc = 14.96;
# AICc weight approximately 0), indicating that models including water content
# were much better supported than a model including topography alone.
#
# The AICc weights describe relative support among only these three candidate
# models. They are not probabilities that any model is true. The biological
# interpretation should therefore also consider the estimated effects,
# confidence intervals, model diagnostics, and predicted probabilities.

# Grouped-binomial example: seed germination ---------------------------------

# The first part of the tutorial used one 0/1 response for each soil core.
# Grouped binomial data have a different structure: each row records the
# number of successes and failures among a known number of trials.
#
# The mead.germination data set records the number of seeds germinating out of
# 50 seeds tested under different temperature regimens and concentrations.
# This example shows how to use this data structure in glm().
data("mead.germination", package = "agridat")

dat_germination <- agridat::mead.germination |>
  as_tibble() |>
  mutate(
    temp = factor(temp),
    failures = seeds - germ,
    germination_proportion = germ / seeds,
    log_conc = log10(conc + 0.01) # add a small constant to avoid log10(0)
  )

glimpse(dat_germination)

# Check the grouped-response variables. The number of failures must be
# nonnegative, and the number of successes cannot exceed the number tested.
dat_germination |>
  summarise(
    n_rows = n(),
    minimum_successes = min(germ),
    maximum_successes = max(germ),
    minimum_failures = min(failures),
    maximum_failures = max(failures),
    unique_trial_sizes = n_distinct(seeds),
    minimum_proportion = min(germination_proportion),
    maximum_proportion = max(germination_proportion)
  )

levels(dat_germination$temp)

# Explore the observed germination proportions. The small constant added
# before log10() allows observations at concentration 0 to remain in the
# analysis. The transformation should be explained clearly in a report.
plot_germination_observed <- ggplot(
  dat_germination,
  aes(
    x = log_conc,
    y = germination_proportion,
    colour = temp
  )
) +
  geom_point(
    alpha = 0.50,
    size = 2
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1),
    expand = c(0, 0.05)
  ) +
  labs(
    x = "Log10(chemical concentration + 0.01)",
    y = "Observed germination proportion",
    colour = "Temperature regimen"
  ) +
  theme_classic()

plot_germination_observed

# Fit a full grouped-binomial model (concentration * temp allows each
# temperature regimen to have its own relationship with concentration). The
# cbind() response tells glm() how many successes and failures occurred in each
# row. This preserves the denominator of the proportion and is not equivalent to
# fitting an ordinary linear model to germination_proportion.
# The interaction allows separate concentration slopes on the log-odds scale.
mod_germination <- glm(
  cbind(germ, failures) ~ log_conc * temp,
  family = binomial(link = "logit"),
  data = dat_germination
)

summary(mod_germination)

coef_germination <- broom::tidy(
  mod_germination,
  conf.int = TRUE,
  conf.level = 0.95,
  exponentiate = TRUE
)

coef_germination

# The exponentiated intercept gives the germination odds for the reference
# temperature regimen when log_conc = 0. The exponentiated log_conc coefficient
# gives the odds ratio for a one-unit increase in log_conc in that regimen.
# Interaction terms describe how this odds ratio changes in other regimens.
# Predicted germination percentages below are easier to interpret biologically.

# Estimate germination probabilities and 95% confidence intervals for each
# concentration-temperature combination. Explicitly specify the observed
# predictor values so that numeric predictors are not held at their means.
# The * here requests every combination; the interaction was fitted in glm().
emm_germination <- emmeans::emmeans(
  mod_germination,
  ~ log_conc * temp,
  at = list(
    log_conc = sort(unique(dat_germination$log_conc)),
    temp = sort(unique(dat_germination$temp))
  ),
  type = "response",
  level = 0.95
)

emm_germination

# Convert probabilities and confidence limits to percentages for reporting.
emm_germination_df <- emm_germination |>
  tibble::as_tibble() |>
  dplyr::transmute(
    log_conc,
    temp,
    percent_germinated = 100 * prob,
    lower_95_CI = 100 * asymp.LCL,
    upper_95_CI = 100 * asymp.UCL
  )

emm_germination_df

# These estimates describe the expected percentage of seeds germinating at
# each combination. The confidence intervals describe uncertainty in that
# expected percentage, rather than the range expected in a new batch of seeds.

# Check diagnostics for the grouped-binomial model.
set.seed(24)
sim_germination <- DHARMa::simulateResiduals(
  fittedModel = mod_germination,
  n = 2000
)

# Visually check uniformity, residual patterns, and simulation outliers.
plot(sim_germination)

# Test whether the residuals follow the expected uniform distribution.
DHARMa::testUniformity(sim_germination)

# Test for more or less residual variability than expected under the model.
DHARMa::testDispersion(sim_germination)

# Test whether the number of simulation outliers differs from expectations.
DHARMa::testOutliers(sim_germination)

# Plot predicted germination probabilities with 95% confidence intervals.
pred_germination <- ggeffects::predict_response(
  mod_germination,
  terms = c("log_conc [all]", "temp")
)

pred_germination_df <- as_tibble(pred_germination)

plot_germination_predictions <- ggplot(
  pred_germination_df,
  aes(
    x = x,
    y = predicted,
    colour = group,
    fill = group,
    group = group
  )
) +
  geom_ribbon(
    aes(
      ymin = conf.low,
      ymax = conf.high
    ),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(
    data = dat_germination,
    aes(
      x = log_conc,
      y = germination_proportion,
      colour = temp
    ),
    inherit.aes = FALSE,
    alpha = 0.50,
    size = 2
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1),
    expand = c(0, 0.05)
  ) +
  labs(
    x = "Log10(chemical concentration + 0.01)",
    y = "Predicted probability of germination",
    colour = "Temperature regimen",
    fill = "Temperature regimen",
    caption = paste(
      "Lines show predicted probabilities;",
      "ribbons show 95% confidence intervals."
    )
  ) +
  theme_classic()

plot_germination_predictions
