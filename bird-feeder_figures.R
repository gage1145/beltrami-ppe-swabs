library(tidyverse)
library(janitor)
library(skimr)
library(broom)
library(ggpubr)
library(ggrepel)


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
    across(c("mpr", "ms", "auc"), mean),
    n = n(),
    .by = c(feeder, site, side, surface, control, reextracted, time_in_env, time_since_start)
  ) %>%
  filter(!control, time_in_env > 0, surface != "Plastic") %>%
  pivot_longer(c("mpr", "ms", "auc"), names_to = "metric")

df_neg <- df_ %>%
  filter(control) %>%
  pivot_longer(c("mpr", "ms", "auc"), names_to = "metric")

df_neg_sum <- df_neg %>%
  summarize(neg_mean = mean(value), neg_sd = sd(value), .by = metric)

df_site <- df_ %>%
  mutate(time_since_start = as.numeric(time_since_start)) %>%
  filter(!control, time_in_env > 0, surface != "Plastic") %>%
  pivot_longer(c("mpr", "ms", "auc"), names_to = "metric")

groupings <- expand.grid(
    m = unique(df_site$metric),
    st = unique(df_site$site),
    sd = unique(df_site$side),
    tss = unique(df_site$time_since_start)
)


get_comps <- function(m, st, sd, tss) {
  temp <- filter(df_sum, metric == m, site == st, side == sd, time_since_start == tss)
  if (nrow(temp) < 1) return(data.frame())
  temp_neg <- filter(df_neg, metric == m)
  value_map <- list(mpr = 25, ms = 7.7, auc = 800)
#   scalar <- (st == "dump site") * 0.12 + 1

  tidy(wilcox.test(temp$value, temp_neg$value, alternative = "greater", exact = FALSE)) %>%
    mutate(
        metric = m, 
        site = st, 
        side = sd,
        time_since_start = tss,
        value=value_map[[m]], 
        significant = p.value <= 0.05,
        p.value=paste("p =", signif(p.value, 3))
    )
}

comps <- pmap_dfr(groupings, get_comps)

df_sum %>%
  ggplot(aes(x=time_since_start, y=value, color=site)) +
  geom_point() +
  stat_smooth(method="lm", linewidth=1.5) +
  geom_hline(aes(yintercept=neg_mean), data=df_neg_sum, linetype="dashed", linewidth=1.2) +
  geom_label_repel(aes(time_since_start, value, label=p.value, color=site), size = 6, data=comps, inherit.aes = FALSE) +
#   geom_label(aes(y=value, label=p.value), data=comps, size=6, position = position_jitter(0.1, 2, seed = 82313498)) +
  facet_grid(vars(metric), vars(side), scales="free") +
  theme(
    axis.text = element_text(size=16),
    axis.title = element_text(size=20),
    strip.text = element_text(size=20),
    legend.title = element_blank(),
    legend.text = element_text(size=16),
  )

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
