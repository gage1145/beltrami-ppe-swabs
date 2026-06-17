# ppe_figures.R
#
# Generates all figures for the PPE swab analysis and saves them to
# figures/ppe/. Reads the curated, metadata-joined data produced by
# ppe_curate.R (data/ppe/updated_metadata.csv).
#
# The analysis collapses the RT-QuIC metrics (MPR, MS, AUC) into principal
# components and compares contamination signal across location (Core/dump site
# vs. negative site), timing (pre- vs. post-sampling), and surface side
# (top vs. bottom of the boot).
#
# Scope: only blank and boot swabs from the negative and dump sites are kept;
# perimeter samples and post-decon timing are excluded.
#
# Exported figures (figures/ppe/):
#   corplot.png                  - pairwise metric correlations (hex + lm fit)
#   distributions.png            - metric histograms split by a manual cutoff
#   pca_plot.png                 - PCA score scatter (PC pairs)
#   comparisons.png              - scaled metrics by location, faceted by timing
#   roc.png                      - ROC curves per metric/timing
#   timing_box.png               - PC1 by timing, faceted by location
#   location_box.png             - PC1 by location, faceted by timing
#   boxplots.png                 - timing + location boxplots stacked (A/B)
#   surface-location_box.png     - PC1 by boot side (top/bottom)
#   location_timing_density.png  - PC1 density ridges by location/timing
#
# Several intermediate plots (metric ridges, normality diagnostics, logistic
# smooths) are drawn for inspection only and are deliberately not exported.

library(tidyverse)
library(magrittr)   # %T>% (tee pipe) used in the PCA normality diagnostic
library(ggpubr)     # ggarrange, stat_compare_means
library(ggridges)   # geom_density_ridges
library(pROC)       # roc, coords, ci.auc
library(scales)     # rescale
library(modelr)     # data_grid, add_predictions


set.seed(321)

# Shared theme applied to every exported figure.
main_theme <- theme(
  plot.title = element_text(size = 30, hjust = 0.5),
  axis.title = element_text(size = 20),
  axis.text = element_text(size = 16),
  strip.text = element_text(size = 20),
  legend.title = element_text(size = 20, hjust = 0.5),
  legend.text = element_text(size = 16),
)


# Load data --------------------------------------------------------------
#
# Read the curated metadata and restrict to the comparison of interest:
# blank/boot swabs from the negative and dump (Core) sites, pre/post sampling.
# Factor levels are set so the negative site and pre-sampling act as the
# reference (baseline) levels in plots and models.


data <- read.csv("data/ppe/updated_metadata.csv", check.names = FALSE) %>%
  mutate(
    across(c("location", "surface", "timing"), tolower),
    timing = str_replace(timing, " ", "-")
  ) %>%
  filter(
    if_any(c(mpr, ms, auc, qr), ~ !is.na(.)),
    location != "perimeter",
    timing != "post-decon",
    surface %in% c("blank", "boot")
  ) %>%
  mutate(
    timing = factor(
      timing,
      levels = c("pre-sampling", "post-sampling"),
      labels = c("Pre-sampling", "Post-sampling")
    ),
    location = factor(
      location,
      levels = c("negative site", "dump site"),
      labels = c("Negative Site", "Core Site")
    ),
    surface = as.factor(surface),
    surface_location = factor(
        tolower(surface_location),
        levels = c("top", "bottom"),
        labels = c("Top", "Bottom")
    )
  )


# Metric Correlations ----------------------------------------------------
#
# For each pair of metrics, plot the joint (log) distribution as a hex
# density with a linear fit overlaid, to show how strongly the metrics
# co-vary. corplot.png.


metrics <- c("ms", "mpr", "auc")
metric_combos <- as.data.frame(t(combn(metrics, 2)))
colnames(metric_combos) <- c("x", "y")

make_cor_plot <- function(xval, yval) {
  form <- as.formula(sprintf("%s ~ %s", yval, xval))
  mod <- lm(form, data)

  grid <- data %>%
    data_grid(!!sym(xval)) %>%
    add_predictions(mod)

  data %>%
    ggplot(aes(x = log(!!sym(xval)), y = log(!!sym(yval)))) +
    geom_hex(bins=50) +
    geom_line(aes(y=log(pred)), data=grid, color="red", linewidth=1.5) +
    coord_fixed() +
    scale_color_manual(values=c("darkcyan", "darkorange")) +
    scale_fill_gradient(low="darkcyan", high="darkorange") +
    main_theme +
    theme(
      legend.key.width = unit(7, "cm")
    )
}

cor_plots <- pmap(metric_combos, make_cor_plot)
ggarrange(plotlist = cor_plots, align = "hv", ncol=3, common.legend = TRUE, legend="bottom")
ggsave("corplot.png", path = "figures/ppe", width = 16, height = 5.6)


# Metric Distributions (diagnostic, not exported) ------------------------
#
# Density ridges of each metric by timing and location. Used to eyeball
# whether the metrics separate by group before modelling.


data %>%
  pivot_longer(c(mpr, ms, auc, qr), names_to = "metric") %>%
  ggplot(aes(value, timing, fill=location)) +
  geom_density_ridges(alpha=0.5, scale=1.25) +
  facet_wrap(vars(metric), scales="free_x") +
  scale_color_manual(values=c("darkcyan", "darkorange")) +
  scale_fill_manual(values=c("darkcyan", "darkorange")) +
  scale_x_log10(expand=c(0.01, 0)) +
  scale_y_discrete(expand=expansion(c(0.05, 1.5))) +
  main_theme +
  theme(
    axis.title = element_blank(),
    legend.position = "top"
  )


# Metric normality -------------------------------------------------------
#
# Each metric is roughly bimodal (a negative cluster and a positive
# contamination cluster). Split each metric at a manual log cutoff and run
# Shapiro-Wilk on each side to check whether the two modes are individually
# normal. df_norm holds the per-mode p-values; the histogram is exported as
# distributions.png.


df_norm <- data %>%
  mutate(
    across(c(mpr, ms, auc), list(log=log)),
  ) %>%
  pivot_longer(c(mpr, ms, auc, mpr_log, ms_log, auc_log), names_to = "metric") %>%
  mutate(
    pos = case_when(
      metric == "mpr_log" ~ value > 0.9,
      metric == "mpr" ~ log(value) > 0.9,
      metric == "ms_log" ~ value > -2,
      metric == "ms" ~ log(value) > -2,
      metric == "auc_log" ~ value > 4.5,
      metric == "auc" ~ log(value) > 4.5
    )
  ) %>%
  group_by(metric, pos) %>%
  summarize(
    norm = shapiro.test(value)$p.value,
    std = sd(value)
  )
df_norm

# Shapiro p-values per metric/mode (diagnostic, not exported).
df_norm %>%
  ggplot(aes(norm, metric, color=pos)) +
  geom_point(size=3) +
  scale_x_log10() +
  main_theme +
  theme(
    axis.title = element_blank(),
    legend.position = "top"
  )

# Histograms of each metric coloured by the cutoff (the two modes). distributions.png.
data %>%
  mutate(
    across(c(mpr, ms, auc), list(log=log)),
  ) %>%
  pivot_longer(c(mpr, ms, auc, mpr_log, ms_log, auc_log), names_to = "metric") %>%
  mutate(
    pos = case_when(
      metric == "mpr_log" ~ value > 0.9,
      metric == "mpr" ~ log(value) > 0.9,
      metric == "ms_log" ~ value > -2,
      metric == "ms" ~ log(value) > -2,
      metric == "auc_log" ~ value > 4.5,
      metric == "auc" ~ log(value) > 4.5
    )
  ) %>%
  ggplot(aes(value, color=pos, fill=pos)) +
  geom_histogram(position = "stack", alpha=0.5, bins=50) +
  facet_grid(cols=vars(metric), scales="free") +
  scale_color_manual(values=c("darkcyan", "darkorange")) +
  scale_fill_manual(values=c("darkcyan", "darkorange")) +
  main_theme +
  theme(
    axis.title = element_blank(),
    legend.position = "top"
  )
ggsave("distributions.png", path="figures/ppe", width=16, height=12)


# Principal Component Analysis -------------------------------------------
#
# Collapse the three log-metrics into principal components. PC1 captures the
# shared "amount of signal" axis and is used as the summary contamination
# score in the boxplots and density plots below. Scores are appended to the
# data as PC1-PC3.


pca_data <- data %>%
  mutate(across(c(mpr, ms, auc), log))

pca_res <- prcomp(~ mpr + ms + auc, data = pca_data, scale. = TRUE)
summary(pca_res)

pca_df <- pca_data %>%
  cbind(pca_res$x)


# PCA score normality (diagnostic, not exported) -------------------------
#
# Histogram of each PC split at zero, plus Shapiro-Wilk per side. The %T>%
# tee draws the plot mid-pipeline while the summary p-values flow through.


pca_df %>%
  pivot_longer(c(PC1, PC2, PC3), names_to = "metric") %>%
  mutate(pos = value > 0) %T>%
  {
    print(
      ggplot(., aes(value, color=pos, fill=pos)) +
        geom_histogram(position = "identity", alpha=0.5, bins=50) +
        facet_grid(cols=vars(metric), scales="free") +
        scale_color_manual(values=c("darkcyan", "darkorange")) +
        scale_fill_manual(values=c("darkcyan", "darkorange")) +
        main_theme +
        theme(
          axis.title = element_blank(),
          legend.position = "top"
        )
      )
  } %>%
  summarise(
    norm = shapiro.test(value)$p.value,
    std = sd(value),
    .by = c(metric, pos)
  )


# PCA score scatter ------------------------------------------------------
#
# Scatter of each PC pair, points coloured by location and overlaid with a
# 2-cluster k-means ellipse, to show whether locations separate in PC space.
# pca_plot.png.


make_pca_plot <- function(x, y, alpha=0.5) {
  cls <- pca_df %>%
    select(!!sym(x), !!sym(y)) %>%
    kmeans(centers=2, trace=TRUE)
  pca_df %>%
    mutate(cluster = as.character(cls$cluster)) %>%
    arrange(desc(location)) %>%
    ggplot(aes(x = !!sym(x), y = !!sym(y))) +
    geom_point(aes(color=location), alpha = alpha, size=3) +
    stat_ellipse(aes(fill=cluster), geom="polygon", alpha=0, type="t", color="black", show.legend = FALSE) +
    scale_color_manual(values=c("darkcyan", "darkorange")) +
    scale_fill_manual(values=c("darkcyan", "darkorange")) +
    coord_fixed() +
    labs(
      x = toupper(x),
      y = toupper(y)
    ) +
    main_theme
}

pca_combos <- as.data.frame(t(combn(paste0("PC", seq(length(pca_res$scale))), 2)))
colnames(pca_combos) <- c("x", "y")

pca_plots <- pmap(pca_combos, make_pca_plot)
ggarrange(plotlist = pca_plots, align = "hv", ncol=1, common.legend = TRUE, legend="bottom")
ggsave("pca_plot.png", path="figures/ppe", width = 12, height = 5)


# Metric comparisons by location -----------------------------------------
#
# Long-format metrics + PCs, each rescaled to 0-1 for a common axis, compared
# between locations within each timing via Wilcoxon tests (significance stars).
# comparisons.png.


pca_long <- pca_df %>%
  pivot_longer(c(mpr, ms, auc, PC1, PC2, PC3), names_to="metric") %>%
  mutate(metric = factor(metric, levels=c("mpr", "ms", "auc", "PC1", "PC2", "PC3")))

pca_df %>%
  mutate(across(c(mpr, ms, auc, PC1, PC2, PC3), ~ rescale(., c(0, 1)))) %>%
  pivot_longer(c(mpr, ms, auc, PC1, PC2, PC3), names_to="metric") %>%
  ggplot(aes(metric, value, fill=location, color=location)) +
  geom_point(position=position_jitterdodge(0.4, 0, 0.7), alpha=0.3) +
  geom_boxplot(outlier.color=NA, color="black", width=0.6, position=position_dodge(0.7)) +
  stat_compare_means(
    aes(label = after_stat(p.signif)),
    method="wilcox.test", size=8, vjust=0.5, show.legend=FALSE
  ) +
  facet_grid(~ timing) +
  scale_fill_manual(values = c("darkcyan", "darkorange")) +
  scale_color_manual(values = c("darkcyan", "darkorange")) +
  main_theme +
  labs(
    y = "Scaled Value (0-1)"
  ) +
  theme(
    legend.position = "top",
  )
ggsave("comparisons.png", path="figures/ppe", width=12, height=10)

# Logistic smooth of each metric against location (diagnostic, not exported):
# how well each metric separates Core vs. negative site.
pca_df %>%
  mutate(across(c(mpr, ms, auc, PC1, PC2, PC3), ~ rescale(., c(-1, 1)))) %>%
  pivot_longer(c(mpr, ms, auc, PC1, PC2, PC3), names_to="metric") %>%
  mutate(location = ifelse(location == "Core Site", 1, 0)) %>%
  ggplot(aes(value, location, color=metric, fill=metric)) +
  geom_point(position=position_jitter(0, 0.01), alpha=0.1) +
  stat_smooth(method="glm", se=FALSE, method.args=list(family="binomial")) +
  facet_wrap(vars(metric), ncol=3, scales="free_x") +
  main_theme


# ROC Analysis -----------------------------------------------------------
#
# For each metric x timing, build an ROC curve for discriminating Core vs.
# negative site and record the AUC (with CI). The curves are overlaid per
# timing to compare which metric discriminates best. roc.png.


get_rocs <- function(m, t, df) {
    temp_df <- df %>%
      filter(
        timing == t,
        metric == m
      ) %>%
      roc(location, value)
    coords(temp_df) %>%
      mutate(
        timing = t,
        metric = m,
        auc = ci.auc(temp_df)[[2]]
      )
}

var_list <- expand.grid(m = unique(pca_long$metric), t = unique(pca_long$timing))
roc_df <- pmap_dfr(var_list, get_rocs, df=pca_long)

roc_df %>%
  arrange(auc, sensitivity) %>%
  ggplot(aes(specificity, sensitivity, color=metric)) +
  geom_step(linewidth=1.2) +
  geom_abline(
    slope = 1, intercept = 1, linetype = "dashed",
    color = "gray50", inherit.aes = FALSE
  ) +
  scale_x_reverse(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  facet_grid(~ timing) +
  labs(
    x = "Specificity", y = "Sensitivity"
  ) +
  main_theme +
  theme(
    legend.title = element_blank(),
  )
ggsave("roc.png", path="figures/ppe", width=12, height=12)


# PC1 Box Plots ----------------------------------------------------------
#
# PC1 (the contamination score) compared across the experimental factors with
# Wilcoxon tests. timing_box and location_box are also stacked into a single
# labelled A/B panel (boxplots.png), and PC1 is compared across boot side.


# Wilcoxon comparison pairs for each factor.
smpl_comps <- list(c("Pre-sampling", "Post-sampling"))
site_comps <- list(c("Negative Site", "Core Site"))
sl_comps   <- list(c("Bottom", "Top"))

# PC1 by timing, faceted by location. timing_box.png.
tim_box <- pca_df %>%
  ggplot(aes(timing, PC1, fill = timing)) +
  geom_boxplot(show.legend=FALSE, outliers = FALSE) +
  geom_dotplot(
    show.legend=FALSE, binaxis="y", stackdir="center",
    dotsize=0.2, stackratio=0.3, color="#333", fill="#333",
    position=position_jitter(0.01, 0.1), alpha=0.25
  ) +
  stat_compare_means(
    aes(group = timing),
    comparisons = smpl_comps,
    method = "wilcox.test",
    label = "p.format", size = 6, show.legend = FALSE
  ) +
  facet_grid(cols=vars(location)) +
  scale_fill_manual(values = c("darkcyan", "darkorange", "maroon")) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  main_theme +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    axis.title.x = element_blank(),
  )
ggsave("timing_box.png", tim_box, path = "figures/ppe", width = 12, height = 8)

# PC1 by location, faceted by timing (one-sided: Core > negative). location_box.png.
loc_box <- pca_df %>%
  ggplot(aes(location, PC1, fill = timing, group = location)) +
  geom_boxplot(show.legend=FALSE, outliers = FALSE) +
  geom_dotplot(
    show.legend=FALSE, binaxis="y", stackdir="center",
    dotsize=0.2, stackratio=0.3, color="#333", fill="#333",
    position=position_jitter(0.01, 0.1), alpha=0.25
  ) +
  stat_compare_means(
    method = "wilcox.test", method.args = list(alternative = "less"),
    comparisons = site_comps, label = "p.format", ref.group = "negative site",
    bracket.size = 0.5, size = 6, show.legend = FALSE,
  ) +
  scale_fill_manual(values = c("darkcyan", "darkorange", "maroon")) +
  facet_grid(cols = vars(timing)) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  main_theme +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    axis.title.x = element_blank(),
  )
ggsave("location_box.png", loc_box, path = "figures/ppe", width = 12, height = 8)

# Stacked location (A) + timing (B) panel. boxplots.png.
ggarrange(
  loc_box, tim_box, align = "hv", ncol=1, common.legend = TRUE, legend="bottom",
  labels = c("A", "B"), font.label = list(size = 24)
)
ggsave("boxplots.png", path = "figures/ppe", width = 12, height = 12)

# PC1 by boot side (top vs. bottom), faceted by timing. surface-location_box.png.
pca_df %>%
  filter(!is.na(surface_location)) %>%
  ggplot(aes(surface_location, PC1, fill = timing, group=surface_location)) +
  geom_boxplot(show.legend=FALSE, outliers = FALSE) +
  geom_dotplot(
    show.legend=FALSE, binaxis="y", stackdir="center",
    dotsize=0.2, stackratio=0.3, color="#333", fill="#333",
    position=position_jitter(0.01, 0.1), alpha=0.25
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = sl_comps, label = "p.format",
    bracket.size = 0.5, size = 6, show.legend = FALSE,
  ) +
  scale_fill_manual(values = c("darkcyan", "darkorange")) +
  facet_grid(cols = vars(timing)) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  main_theme +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    axis.title.x = element_blank(),
  )
ggsave("surface-location_box.png", path = "figures/ppe", width = 12, height = 8)


# Density Plots ----------------------------------------------------------
#
# PC1 density ridges by timing, filled by location, summarising the overall
# shift in contamination signal. location_timing_density.png.


pca_df %>%
  ggplot(aes(PC1, timing, fill = location)) +
  geom_density_ridges(alpha = 0.5, scale = 1.2) +
  scale_fill_manual(values = c("darkcyan", "darkorange", "red")) +
  scale_y_discrete(expand = expansion(mult = c(0.025, 0.5))) +
  main_theme +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.title.y = element_blank()
  )
ggsave("location_timing_density.png", path = "figures/ppe", width = 14, height = 8)
