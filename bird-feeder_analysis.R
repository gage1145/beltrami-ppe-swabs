library(tidyverse)
library(janitor)
library(lme4)
library(emmeans)
library(broom)
library(ggpubr)



df_ <- read.csv("data/bird-feeder/data.csv", check.names = FALSE) %>%
  clean_names(replace = c("Sample IDs" = "sample_ids", "TtT" = "ttt"))

# SS field feeder swabs with Top/Bottom sides — main analysis dataset.
# Plastic swabs have no side assignment and are analyzed separately below.
df_ss <- df_ %>%
  filter(!control, side %in% c("Top", "Bottom")) %>%
  mutate(
    feeder      = as.factor(feeder),
    side        = factor(side, levels = c("Top", "Bottom")),
    site        = factor(site, levels = c("perimeter", "dump site")),
    reextracted = as.logical(reextracted),
    time_in_env = as.numeric(time_in_env)
  )



# Helper: tidy pairwise contrasts for a categorical factor
get_pairwise_comparisons <- function(factor_name, model) {
  em <- emmeans(model, as.formula(paste("pairwise ~", factor_name)))
  tidy(em$contrasts) %>%
    mutate(factor = factor_name) %>%
    {if ("adj.p.value" %in% colnames(.)) rename(., p.value = adj.p.value) else .}
}



# LME models (SS swabs) -------------------------------------------------------
# surface is excluded: all SS swabs, and Plastic lacks side info (see below).
# Fixed effects: time_in_env (continuous), side, reextracted (methodological covariate).
# Random effect: feeder (repeated swabs of the same sentinel).
# lme4 does not provide p-values for fixed effects without lmerTest;
# use the t statistic (|t| > 2 is a common rule of thumb for significance).



metrics <- c("mpr", "ms", "auc")

model_list <- lapply(metrics, function(m) {
  lmer(
    as.formula(paste(m, "~ time_in_env + side + reextracted + site + (1 | feeder)")),
    data = df_ss
  )
})
names(model_list) <- metrics

# Fixed effects table
fixed_effects <- map_dfr(metrics, function(m) {
  coefs <- as.data.frame(coef(summary(model_list[[m]])))
  coefs$term <- rownames(coefs)
  rownames(coefs) <- NULL
  coefs$metric <- m
  coefs
}) %>%
  rename(estimate = Estimate, std.error = `Std. Error`, statistic = `t value`) %>%
  mutate(significant = abs(statistic) > 2) 

# Pairwise contrasts for side
side_results <- map_dfr(metrics, function(m) {
  get_pairwise_comparisons("side", model_list[[m]]) %>%
    mutate(metric = m)
})

# Pairwise contrasts for reextracted
extraction_results <- map_dfr(metrics, function(m) {
  get_pairwise_comparisons("reextracted", model_list[[m]]) %>%
    mutate(metric = m)
})

# Pairwise contrasts for site (perimeter vs. dump site)
site_results <- map_dfr(metrics, function(m) {
  get_pairwise_comparisons("site", model_list[[m]]) %>%
    mutate(metric = m)
})

write.csv(fixed_effects, "data/bird-feeder/fixed_effects.csv", row.names = FALSE)
write.csv(side_results, "data/bird-feeder/side_results.csv", row.names = FALSE)
write.csv(extraction_results, "data/bird-feeder/extraction_results.csv", row.names = FALSE)
write.csv(site_results, "data/bird-feeder/site_results.csv", row.names = FALSE)



# Surface comparison (SS vs. Plastic) -----------------------------------------
# Plastic swabs have no side data; compare as a separate Wilcoxon test.



df_surface <- df_ %>%
  filter(!control, surface %in% c("SS", "Plastic"))

surface_results <- map_dfr(metrics, function(m) {
  ss_vals  <- df_surface[[m]][df_surface$surface == "SS"]
  pl_vals  <- df_surface[[m]][df_surface$surface == "Plastic"]
  test <- wilcox.test(ss_vals, pl_vals, exact = FALSE)
  data.frame(
    metric    = m,
    statistic = test$statistic,
    p.value   = test$p.value,
    n_ss      = length(ss_vals),
    n_plastic = length(pl_vals)
  )
})

write.csv(surface_results, "data/bird-feeder/surface_comparison.csv", row.names = FALSE)



# Control vs. field comparison ------------------------------------------------



df_ctrl_compare <- df_ %>%
  filter(side %in% c("Top", "Bottom") | control)

ctrl_results <- map_dfr(metrics, function(m) {
  field_vals <- df_ctrl_compare[[m]][!df_ctrl_compare$control]
  ctrl_vals  <- df_ctrl_compare[[m]][df_ctrl_compare$control]
  test <- wilcox.test(field_vals, ctrl_vals, exact = FALSE)
  data.frame(
    metric    = m,
    statistic = test$statistic,
    p.value   = test$p.value,
    n_field   = length(field_vals),
    n_ctrl    = length(ctrl_vals)
  )
})

write.csv(ctrl_results, "data/bird-feeder/ctrl_comparison.csv", row.names = FALSE)

df_ctrl_compare %>%
  pivot_longer(c("mpr", "ms", "auc"), names_to = "metric") %>%
  mutate(control = ifelse(control, "Negative Site", "Positive Site")) %>%
  ggplot(aes(control, value, fill=control)) +
  geom_violin() +
  stat_compare_means() +
  facet_wrap(vars(metric), ncol=3, scales="free")
