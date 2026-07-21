# bird-feeder_figures.R
#
# Generates the figures for the sentinel ("bird feeder") analysis and saves
# them to figures/bird-feeder/. Reads the curated data produced by
# bird-feeder_curate.R (data/bird-feeder/data.csv).
#
# Sentinels are stainless-steel plates mounted on bird-feeder stands that were
# left in the environment at dump-site/perimeter locations (and negative
# control sites) and periodically swabbed. RT-QuIC metrics (MPR, MS, AUC) are
# collapsed to principal components; PC1 is used as the contamination score and
# compared against the negative controls and tracked over time.
#
# Scope: stainless-steel surfaces only; the per-plate well controls (sample
# IDs "N"/"P") are dropped.
#
# Outputs:
#   figures/bird-feeder/pca.png            - PCA score scatter (PC pairs)
#   figures/bird-feeder/pca_vs_tss.png     - PC1 vs. year, by site/side
#   figures/bird-feeder/plate_cols.png     - mean PC1 per plate/site, by year
#   data/bird-feeder/tss_results.csv       - per-group lm(PC1 ~ time) fit stats


library(tidyverse)
library(janitor)   # clean_names
library(skimr)     # skim
library(broom)     # tidy, glance
library(ggpubr)    # ggarrange, stat_cor, stat_compare_means
library(magrittr)  # %T>% (tee pipe)
library(Rdimtools)  # supervised pca


main_theme <- theme(
  plot.title = element_text(size = 30, hjust = 0.5),
  axis.title = element_text(size = 20),
  axis.text = element_text(size = 16),
  strip.text = element_text(size = 20),
  legend.title = element_text(size = 20, hjust = 0.5),
  legend.text = element_text(size = 16),
  panel.grid = element_blank()
)


# Load data --------------------------------------------------------------
#
# Read the curated sentinel data, drop plate controls and non-stainless
# surfaces, and coerce columns to their proper types (factors for the
# categorical fields, dates, and day-valued time differences).


df_ <- read.csv("data/bird-feeder/data.csv", check.names = FALSE) %>%
  clean_names(replace = c("Sample IDs" = "sample_ids", "TtT" = "ttt")) %>%
  filter(
    !(sample_ids %in% c("N", "P")),
    surface == "SS"
  ) %>%
  select(-c("dilutions", "tech", "collector", "surface")) %>%
  mutate(
    across(c(mpr, ms, auc, raf, ttt), as.numeric),
    across(c(sample_ids:reader, feeder:side, site), as.factor),
    across(start_date:end_date, as.Date),
    across(c(time_in_env, time_since_start), ~ as.difftime(., units = "days")),
    year = as.factor(year(end_date)),
  )

df_ %<>%
  mutate(
    site_pca = do.spc(select(df_, mpr, ms, auc),
      pull(df_, site),
      1
    )$Y
  )

df_long <- df_ %>%
  pivot_longer(
    cols = c(mpr, ms, auc, site_pca),
    names_to = "metric",
    values_to = "value"
  )


# Summary / group comparisons --------------------------------------------
#
# Build per-group means and compare each site/side/timepoint against the
# pooled negative-control distribution with a one-sided Wilcoxon test
# (greater). comps holds the per-group test results; these feed the
# significance annotations used downstream. The figures focus on PC1.


df_sum <- df_long %>%
  summarize(
    across(c(value), mean),
    n = n(),
    .by = c(metric, feeder, site, side, control, reextracted, time_in_env, time_since_start)
  ) %>%
  filter(!control, time_in_env > 0)

df_neg <- df_long %>%
  filter(control)

df_neg_sum <- df_neg %>%
  summarize(neg_mean = mean(value), neg_sd = sd(value), .by = metric)

df_site <- df_long %>%
  filter(!control, time_in_env > 0)

groupings <- expand.grid(
  m = unique(df_site$metric),
  st = unique(df_site$site),
  sd = unique(df_site$side),
  tss = unique(df_site$time_since_start)
)

# One group (metric x site x side x timepoint) vs. the negative controls.
get_comps <- function(m, st, sd, tss) {
  temp <- filter(df_sum, metric == m, site == st, side == sd, time_since_start == tss)
  if (nrow(temp) < 1) {
    return(data.frame())
  }
  temp_neg <- filter(df_neg, metric == m)
  value_map <- list(
    mpr = df_neg_sum$neg_mean[df_neg_sum$metric == "mpr"],
    ms  = df_neg_sum$neg_mean[df_neg_sum$metric == "ms"],
    auc = df_neg_sum$neg_mean[df_neg_sum$metric == "auc"],
    pca = df_neg_sum$neg_mean[df_neg_sum$metric == "site_pca"]
  )

  tidy(wilcox.test(temp$value, temp_neg$value, alternative = "greater", exact = FALSE)) %>%
    mutate(
      metric = m,
      site = st,
      side = sd,
      time_since_start = tss,
      value = value_map[[m]],
      significant = p.value <= 0.05,
      p.value = paste("p =", signif(p.value, 3))
    )
}

comps <- pmap_dfr(groupings, get_comps)


# PC1 vs. time-since-start -----------------------------------------------
#
# Track PC1 over sampling years for each site/side. Per-group linear fits
# (PC1 ~ time_since_start) are written to tss_results.csv, and the trend is
# plotted with a smooth and correlation labels. pca_vs_tss.png.

# The time-trend figures use PC1, relative to the mean negative-control PC1.
fig_metric <- "site_pca"
neg_value <- df_neg_sum$neg_mean[df_neg_sum$metric == fig_metric]

# Established sentinels only: long-deployed (>700 days since start),
# non-control, PC1, excluding re-extracted samples.
filter_for_fig <- function(df) {
  df %>%
    filter(
      time_since_start > 700,
      !control,
      metric == fig_metric,
      !reextracted
    )
}

mods <- df_long %>%
  filter_for_fig() %>%
  nest(.by = c(site)) %>%
  mutate(
    model = map(data, ~ lm(value ~ time_since_start, data = .x)),
    glanced = map(model, broom::glance)
  ) %>%
  unnest(c(glanced))
write.csv(select(mods, -c(data, model)), "data/bird-feeder/tss_results.csv")

# Build the plotting frame once so the same rows feed both the plot and the
# per-time-point comparisons below.
fig_data <- df_long %>%
  filter_for_fig() %>%
  mutate(
    site = factor(site, levels = c("dump site", "perimeter"), labels = c("Dump Site", "Perimeter")),
    year = as.integer(as.character(year))
  )

# Negative-site signal as plotted: negated PC1, so larger = more signal. The
# negative site has a single visit, so its points are pooled and reused as the
# shared baseline that every sentinel time-point is tested against.
neg_signal <- -df_neg$value[df_neg$metric == fig_metric]
neg_sd     <- sd(neg_signal)

# Two-sided Wilcoxon per site x year: does the sentinel signal differ from the
# negative-site baseline? The label sits just above each group's highest point.
year_comps <- fig_data %>%
  summarize(
    p.value = wilcox.test(neg_signal, value, exact = FALSE, alternative="greater")$p.value,
    .by = c(site, year)
  ) %>%
  mutate(
    label = sprintf("p = %.3g", p.value),
    y = case_when(site == "Dump Site" ~ 5.6, TRUE ~ 5.1),
  )

fig_data %>%
  ggplot(aes(x = year, y = neg_value - value, color=site, fill=site)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -neg_sd, ymax = neg_sd, 
           fill = "blue", alpha = 0.05) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2, color="blue") +
  geom_point(size=3, position=position_jitterdodge(0.04, 0, 0.1, seed=56543)) +
  geom_smooth(method="lm",linewidth=1.5, alpha=0.5) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    size = 5, label.x.npc=c(0.25, 0.75), label.y=6, show.legend = FALSE, hjust = 0.5
  ) +
  geom_text(
    aes(x = year, y = y, label = label, group = site),
    data = year_comps,
    size = 5, fontface = "bold", show.legend = FALSE
  ) +
  # scale_x_continuous(breaks = seq(2023, 2025)) +
  # scale_y_continuous(breaks = seq(-2, 7), limits=c(NA, 6)) +
  scale_fill_manual(values = c("darkorange", "darkcyan")) +
  scale_color_manual(values = c("darkorange", "darkcyan")) +
  labs(
    x = "Year of Sampling",
    y = "PC1 - mean Neg PC1"
  ) +
  main_theme +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = "top"
  )
ggsave("pca_vs_tss.png", path = "figures/bird-feeder", width = 12, height = 8)


# PC1 by plate -----------------------------------------------------------
#
# Mean PC1 per individual feeder plate, stacked within each site and faceted by
# year, with site-level means/error bars and a Wilcoxon site comparison.
# PC1 is negated so larger bars mean more signal. plate_cols.png.


pca_plates <- df_ %>%
  filter(
    str_detect(feeder, "BF"),
    time_in_env > 0
  )

# Per-plate and per-site mean PC1 (negated) with standard deviations.
plate_sum <- pca_plates %>%
  summarize(
    across(site_pca, c(mean=~-mean(.), stdev=sd)),
    .by = c(feeder, year, site)
  )

site_sum <- pca_plates %>%
  summarize(
    across(site_pca, c(mean=~-mean(.), stdev=sd)),
    .by = c(year, site)
  )

plate_sum %>%
  ggplot(aes(site, site_pca_mean)) +
  geom_col(aes(fill=feeder, group=site)) +
  geom_point(
    aes(site, site_pca_mean, group=site),
    data=site_sum, inherit.aes = FALSE, position = position_dodge(0.9), shape="square", size=3
  ) +
  geom_errorbar(
    aes(site, site_pca_mean, ymin=site_pca_mean - site_pca_stdev, ymax=site_pca_mean + site_pca_stdev, group=site),
    data=site_sum, inherit.aes = FALSE, position = position_dodge(0.9), width=0.2
  ) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=1) +
  stat_compare_means(
    aes(site, site_pca),
    data = pca_plates,
    method = "wilcox.test",
    label.y.npc = 0.5, label.x.npc=0.5, hjust=0.5,
    label = "p.format", size = 6, show.legend = FALSE
  ) +
  facet_grid(cols=vars(year)) +
  scale_y_continuous(breaks = seq(-400, 400, 100)) +
  labs(
    y="PC1"
  ) +
  scale_fill_manual(values=c("cyan4", "cyan3", "cyan2", "orange4", "orange3", "orange2")) +
  main_theme +
  theme(
    axis.title.x = element_blank(),
    legend.title = element_blank()
  )
ggsave("plate_cols.png", path="figures/bird-feeder", width=12, height=10)
