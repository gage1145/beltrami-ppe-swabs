library(tidyverse)
library(janitor)
library(skimr)
library(broom)


df_ <- read.csv("data/bird-feeder/data.csv", check.names = FALSE) %>%
  clean_names(replace=c("Sample IDs" = "sample_ids", "TtT" = "ttt")) %>%
  filter(!(sample_ids %in% c("N", "P"))) %>%
  mutate(
    across(c(sample_ids: tech, feeder: collector, site), as.factor),
    across(start_date: end_date, as.Date),
    across(c(time_in_env, time_since_start), ~ as.difftime(., units = "days"))
  )

skim(df_)

df_sum <- df_ %>%
  summarize(
    across(c("mpr", "ms", "auc"), list(mean=mean, stdev=sd)),
    n = n(),
    .by = c(feeder, site, side, surface, control, reextracted, time_in_env, time_since_start)
  )

df_neg <- df_sum %>%
  filter(control, !reextracted) %>%
  pivot_longer(c("mpr_mean", "ms_mean", "auc_mean"), names_to = "metric") %>%
  summarize(neg_mean = mean(value), neg_sd = sd(value), .by = metric)

df_site <- df_sum %>%
  mutate(time_since_start = as.numeric(time_since_start)) %>%
  filter(!control, !reextracted, time_in_env > 0, surface != "Plastic")

df_site %>%
  pivot_longer(c("mpr_mean", "ms_mean", "auc_mean"), names_to = "metric") %>%
  ggplot(aes(time_since_start, value, color=site)) +
  geom_point() +
  stat_smooth(method="lm") +
  # geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=neg_mean-neg_sd, ymax=neg_mean+neg_sd),
  #           data=df_neg, inherit.aes=FALSE, alpha=0.2) +
  geom_hline(aes(yintercept=neg_mean), data=df_neg, linetype="dashed") +
  facet_grid(vars(metric), vars(side), scales="free") 

metrics <- c("mpr_mean", "ms_mean", "auc_mean")

get_mods <- function(x) {
  form <- as.formula(paste(x, " ~ site * time_since_start * side"))
  mod <- glm(form, data = df_site)
  tidy(mod) %>%
    mutate(metric = x)
}

site_mods <- map_dfr(metrics, get_mods)
# Looks like as time goes on, MPR and AUC detect a significant negative trend in activity.
# There are no significant differences between the two sites as a whole.
# Slopes do not diverge between sites. They have parallel slopes.

# Compare against negative site.
