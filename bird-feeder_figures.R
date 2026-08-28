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
library(janitor)   
library(skimr)     
library(broom)     
library(ggpubr)   
library(ggrepel)
library(ggnewscale)
library(magrittr)  
library(plotly)


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
    feeder = str_replace_all(feeder, "BF", "S"),
    time_since_start = ifelse(time_since_start == 699, 760, time_since_start)
  )


# Principal Component Analysis -------------------------------------------
#
# Collapse the three log-metrics into principal components. PC1 is the shared
# "amount of signal" axis used as the contamination score below. Scores are
# appended to the data, along with the sampling year derived from end_date.


pca_res <- df_ %>%
  select(mpr, ms, auc) %>%
  mutate(across(everything(), log)) %>%
  prcomp(scale. = TRUE)
pca_sum <- summary(pca_res)

pc_lookup <- c(
  paste0("PC1", "(", round(100 * pca_sum$importance[2, "PC1"], 1), "%)") = PC1,
  paste0("PC2", "(", round(100 * pca_sum$importance[2, "PC2"], 1), "%)") = PC2,
  paste0("PC3", "(", round(100 * pca_sum$importance[2, "PC3"], 1), "%)") = PC3
)

get_pc_lab <- function(pc) paste0(paste0(pc, "(", round(100 * pca_sum$importance[2, pc], 1), "%)"))

pca_df <- df_ %>%
  bind_cols(pca_res$x) %>%
  mutate(
    year = as.factor(year(end_date)),
  ) %>%
  rename_with(
    ~ paste0(.x, "(", round(100 * pca_sum$importance[2, .x], 1), "%)"),
    starts_with("PC")
  )

pca_long <- pca_df %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(c("mpr", "ms", "auc", "PC1", "PC2", "PC3"), names_to = "metric", values_to = "value")


# PCA score scatter ------------------------------------------------------
#
# Scatter of each PC pair, points coloured by whether the sample came from a
# negative control site, with a 99% k-means cluster ellipse. Shows whether the
# contaminated sites separate from the controls in PC space. pca.png.


pca_combos <- as.data.frame(t(combn(paste0("PC", seq(length(pca_res$scale))), 2))) %>%
  mutate(
    across(everything(), ~ paste0(., "(", round(100 * pca_sum$importance[2, as.integer(str_remove(., "PC"))], 1), "%)"))
  )
colnames(pca_combos) <- c("x", "y")

make_pca_plot <- function(x, y, alpha=0.5) {
  cls <- pca_df %>%
    select(!!sym(x), !!sym(y)) %>%
    kmeans(centers=2, trace=TRUE) %T>%
    summary()
  pca_df %>%
    mutate(
      cluster = as.character(cls$cluster),
      site = factor(site, levels = c("dump site", "perimeter"), labels = c("Core", "Perimeter")),
    ) %>%
    filter(!control) %>%
    # arrange(control) %>%
    ggplot(aes(x = !!sym(x), y = !!sym(y))) +
    geom_point(aes(color=site), alpha = alpha, size=3) +
    stat_ellipse(aes(fill=cluster), level=0.99, geom="polygon", alpha=0, type="t", color="black", show.legend = FALSE) +
    scale_color_manual(values=c("darkorange", "darkcyan")) +
    scale_fill_manual(values=c("darkorange", "darkcyan")) +
    facet_grid(cols=vars(time_since_start)) +
    scale_x_continuous(n.breaks=3) +
    # coord_fixed() +
    labs(
      x = toupper(x),
      y = toupper(y)
    ) +
    main_theme +
    theme(
      legend.position="none"
    )
}

pca_plots <- pmap(pca_combos, make_pca_plot)
pca_plot <- ggarrange(plotlist = pca_plots, align = "hv", ncol=1, common.legend = FALSE) 
# ggsave("pca.png", path = "figures/bird-feeder", bg="white", width = 12, height = 20)
pca_plot

# Summary / group comparisons --------------------------------------------
#
# Build per-group means and compare each site/side/timepoint against the
# pooled negative-control distribution with a one-sided Wilcoxon test
# (greater). comps holds the per-group test results; these feed the
# significance annotations used downstream. The figures focus on PC1.


df_sum <- pca_df %>%
  filter(!control) %>%
  mutate(
    time_since_start = as.factor(time_since_start),
    site = factor(site, levels = c("dump site", "perimeter"), labels = c("Core", "Perimeter")),
  ) %>%  summarize(n=n(), .by=c(site, time_since_start)) %>%
  mutate(y = ifelse(site=="Core", 2.2, 1.7))

pc_box <- pca_df %>%
  filter(!control) %>%
  select(PC1, site, time_since_start) %>%
  mutate(
    time_since_start = as.factor(time_since_start),
    site = factor(site, levels = c("dump site", "perimeter"), labels = c("Core", "Perimeter")),
  ) %>%
  ggplot(aes(time_since_start, PC1, fill=site, color=site)) +
  geom_hline(yintercept = 0, linetype="dashed", color="gray50", linewidth=1) +
  geom_boxplot(position=position_dodge(0.7), outlier.shape=NA,  color="black", width=0.4) +
  geom_text(
    aes(label = paste("n =", n), group = time_since_start, y=y), fontface = "bold",
    data = df_sum, position = position_dodge(0.7), size=5, vjust=0.5, show.legend=FALSE
  ) +
  stat_compare_means(
    aes(group=site, label = paste0("p = ", after_stat(p.format))),
    size=4, vjust=0.5, fontface="bold", show.legend=FALSE
  ) +
  scale_y_continuous(limits=c(-5, 2.7)) +
  scale_color_manual(values=c("darkorange", "darkcyan")) +
  scale_fill_manual(values=c("darkorange", "darkcyan")) +
  labs(
    x = "Days Since Start",
    y = "PC1"
  ) +
  main_theme +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size=24),
    legend.key.size = unit(1, "cm"),
    # axis.title.x = element_blank()
  )
pc_box

ggarrange(
  pca_plot, pc_box, ncol = 2, common.legend = TRUE, legend="bottom", labels = c("A", "B"),
  font.label = list(size = 24), widths = c(3,2)
)

ggsave("pca_plot.png", path="figures/bird-feeder", width = 12, height = 10, bg="white")


df_sum <- pca_long %>%
  summarize(
    across(c(value), mean),
    n = n(),
    .by = c(metric, feeder, site, side, control, reextracted, time_in_env, time_since_start)
  ) %>%
  filter(!control, time_in_env > 0)

df_neg <- pca_long %>%
  filter(control)

df_neg_sum <- df_neg %>%
  summarize(neg_mean = median(value), neg_sd = sd(value), .by = metric)

df_site <- pca_long %>%
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
    PC1 = df_neg_sum$neg_mean[df_neg_sum$metric == "PC1"],
    PC2 = df_neg_sum$neg_mean[df_neg_sum$metric == "PC2"],
    PC3 = df_neg_sum$neg_mean[df_neg_sum$metric == "PC3"]
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
fig_metric <- "PC1"
neg_value <- df_neg_sum$neg_mean[df_neg_sum$metric == fig_metric]

# Established sentinels only: long-deployed (>700 days since start),
# non-control, PC1, excluding re-extracted samples.
filter_for_fig <- function(df) {
  df %>%
    filter(
      # time_since_start > 700,
      time_since_start > 0,
      !control,
      metric == fig_metric,
      !reextracted
    ) %>%
    mutate(
      time_since_start = ifelse(time_since_start == 699, 760, time_since_start),
      site = factor(site, levels = c("dump site", "perimeter"), labels = c("Core", "Perimeter")),
      year = as.integer(as.character(year)),
      time_since_start = as.integer(time_since_start),
      time_in_env = ifelse(time_in_env < 400, 365, time_in_env)
    )
}


# PC1 vs TSS Plot --------------------------------------------------------


# Build the plotting frame once so the same rows feed both the plot and the
# per-time-point comparisons below.
fig_data <- pca_long %>%
  filter_for_fig()

tss_mods <- fig_data %>%
  nest(.by = c(site)) %>%
  mutate(
    model = map(data, ~ lm(value ~ time_since_start, data = .x)),
    glanced = map(model, broom::glance)
  ) %>%
  unnest(c(glanced))
write.csv(select(tss_mods, -c(data, model)), "data/bird-feeder/tss_results.csv")

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
    .by = c(site, time_since_start)
  ) %>%
  mutate(
    label = sprintf("p = %.3g", p.value),
    y = case_when(site == "Core" ~ 6.6, TRUE ~ 6.1),
  )

fig_data %>%
  ggplot(aes(x = time_since_start, y = neg_value - value, color=site, fill=site)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -neg_sd, ymax = neg_sd, 
           fill = "blue", alpha = 0.05) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2, color="blue") +
  geom_point(size=3, position=position_jitterdodge(10, 0, 20, seed=56543)) +
  geom_smooth(method="lm",linewidth=1.5, alpha=0.5) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    size = 5, label.x.npc=c(0.25, 0.75), label.y=7, show.legend = FALSE, hjust = 0.5
  ) +
  geom_text(
    aes(x = time_since_start, y = y, label = label, group = site),
    data = year_comps,
    size = 5, fontface = "bold", show.legend = FALSE
  ) +
  scale_x_continuous(breaks = seq(750, 1500, 125)) +
  scale_y_continuous(breaks = seq(-2, 8), limits=c(NA, 7)) +
  scale_fill_manual(values = c("darkorange", "darkcyan")) +
  scale_color_manual(values = c("darkorange", "darkcyan")) +
  labs(
    x = "Days Since Project Start",
    y = "PC1 - mean Neg PC1"
  ) +
  main_theme +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 20),
    # axis.title.x = element_blank(),
    strip.text = element_text(size = 20),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.background = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.75, 0.75)
  )
ggsave("pca_vs_tss.png", path = "figures/bird-feeder", width = 12, height = 8)


# PC1 vs TIE Plot --------------------------------------------------------


tie_mods <- fig_data %>%
  nest(.by = c(site)) %>%
  mutate(
    model = map(data, ~ lm(value ~ time_in_env, data = .x)),
    glanced = map(model, broom::glance)
  ) %>%
  unnest(c(glanced))
write.csv(select(tie_mods, -c(data, model)), "data/bird-feeder/tie_results.csv")

time_in_env_comps <- fig_data %>%
  summarize(
    p.value = wilcox.test(neg_signal, value, exact = FALSE, alternative="greater")$p.value,
    .by = c(time_in_env)
  ) %>%
  mutate(
    point = factor(time_in_env, levels = unique(time_in_env), labels = paste(unique(time_in_env), "days")),
    label = sprintf("p = %.3g", p.value),
    y = 6.6,
  )

fig_data %>%
  arrange(time_in_env) %>%
  mutate(
    time_in_env = as.integer(time_in_env),
    point = factor(time_in_env, levels = unique(time_in_env), labels = paste(unique(time_in_env), "days"))
  ) %>%
  ggplot(aes(x = time_in_env, y = neg_value - value)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -neg_sd, ymax = neg_sd, 
           fill = "blue", alpha = 0.05) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2, color="blue") +
  geom_point(aes(color = point),size=3, position=position_jitterdodge(0.04, 0, 0.1, seed=56543)) +
  geom_smooth(method="lm", linewidth=1.5, alpha=0.5) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    size = 5, color="blue", label.x.npc=0.5, label.y=7, show.legend = FALSE, hjust = 0.5
  ) +
  geom_text_repel(
    aes(x = time_in_env, y = y, label = label, color = point),
    data = time_in_env_comps, direction = "x", force = 0.1, hjust = 0.5,
    size = 5, fontface = "bold", show.legend = FALSE
  ) +
  scale_x_continuous(breaks = seq(300, 1500, 300), limits = c(300, 1500)) +
  scale_y_continuous(breaks = seq(-2, 8), limits=c(NA, 7)) +
  labs(
    x = "Days in Environment",
    y = "PC1 - mean Neg PC1"
  ) +
  main_theme +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = "inside",
    legend.position.inside = c(0.75, 0.75)
  )
ggsave("pca_vs_tie.png", path = "figures/bird-feeder", width = 12, height = 8)


# Plots by plate ---------------------------------------------------------


pca_plates <- pca_df %>%
  filter(
    str_detect(feeder, "^S"),
    time_in_env > 0
  )

# Per-plate and per-site mean PC1 (negated) with standard deviations.
plate_sum <- pca_plates %>%
  summarize(
    across(PC1, c(mean=~-mean(.), median=~-median(.), stdev=sd)),
    .by = c(feeder, year, site)
  )

site_sum <- pca_plates %>%
  summarize(
    p.value = wilcox.test(PC1 ~ site)$p.value,
    .by = c(year)
  ) %>%
  mutate(
    label = paste("p =", signif(p.value, 2)),
    y     = 6,
    color = ifelse(p.value < 0.05, "darkred", "black")
  )

# year_sum <- pca_plates %>%
#   summarize(
#     p.value = pairwise.wilcox.test(PC1, year, p.adjust.method = "bonferroni", exact = FALSE)$p.value,
#     .by = c(site)
#   ) 

df_pairwise <- pca_plates %>%
  select(PC1, year, site) %>%
  na.omit()
  
mod <- aov(PC1 ~ site * year, data = df_pairwise)
tukey_results <- TukeyHSD(mod)
summary(tukey_results)
print(tukey_results)

dodge_w <- 0.5

comps <- tukey_results$`site:year` %>%
  as.data.frame() %>%
  rownames_to_column(var = "comparison") %>%
  separate(comparison, c("var1", "var2"), sep = "-") %>%
  separate(var1, c("site1", "year1"), sep = ":") %>%
  separate(var2, c("site2", "year2"), sep = ":") %>%
  filter(site1 == site2) %>%
  mutate(
    label      = paste("p =", signif(`p adj`, 2)),
    color      = ifelse(`p adj` < 0.05, "darkred", "black"),
    y          = c(4.1, 4.9, 4.5, 5.3, 4.1, 4.5),
    site_offset = ifelse(site1 == "dump site", -dodge_w ** 2 / 2, dodge_w ** 2 / 2),
    x1         = as.integer(as.character(year2)) - 2022 + site_offset + 0.01,
    x2         = as.integer(as.character(year1)) - 2022 + site_offset - 0.01,
    x_mid      = (x1 + x2) / 2
  )

pca_plates %>%
  summarize(across(PC1, mean), .by = c(sample_ids, year, site)) %>%
  ggplot(aes(year, -PC1, fill=site)) +
  geom_boxplot(width=0.4, box.color = NA, position = position_dodge(dodge_w)) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=1) +
  scale_fill_manual(values = c("darkorange", "darkcyan")) +
  scale_color_manual(values = c("darkorange", "darkcyan")) +
  new_scale_color() +
  geom_text(# Site Comparisons
    aes(x = year,y = y, label = label, color = color),
    data = site_sum, inherit.aes = FALSE, size = 5, fontface = "bold"
  ) +
  geom_text(# Year Comparisons
    aes(x = x_mid, y = y + 0.2, label = label, color = color),
    data = comps, inherit.aes = FALSE, size = 5, fontface = "bold"
  ) +
  geom_segment(# Horizontal bracket
    aes(x = x1, xend = x2, y = y, yend = y, color = color),
    data = comps, inherit.aes = FALSE, linewidth = 0.5
  ) +
  geom_segment(# Left tick
    aes(x = x1, xend = x1, y = y - 0.15, yend = y, color = color),
    data = comps, inherit.aes = FALSE, linewidth = 0.5
  ) +
  geom_segment(# Right tick
    aes(x = x2, xend = x2, y = y - 0.15, yend = y, color = color),
    data = comps, inherit.aes = FALSE, linewidth = 0.5
  ) +
  scale_color_identity() +
  scale_y_continuous(limits=c(NA, 6.4)) +
  labs(
    y="PC1"
  ) +
  main_theme +
  theme(
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom"
  )
# ggsave("plate_cols.png", path="figures/bird-feeder", width=12, height=10)
