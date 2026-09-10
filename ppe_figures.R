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
#   corplot.jpg                  - pairwise metric correlations (hex + lm fit)
#   distributions.jpg            - metric histograms split by a manual cutoff
#   pca_plot.jpg                 - PCA score scatter (PC pairs)
#   comparisons.jpg              - scaled metrics by location, faceted by timing
#   roc.jpg                      - ROC curves per metric/timing
#   timing_box.jpg               - PC1 by timing, faceted by location
#   location_box.jpg             - PC1 by location, faceted by timing
#   boxplots.jpg                 - timing + location boxplots stacked (A/B)
#   surface-location_box.jpg     - PC1 by boot side (top/bottom)
#   location_timing_density.jpg  - PC1 density ridges by location/timing
#
# Several intermediate plots (metric ridges, normality diagnostics, logistic
# smooths) are drawn for inspection only and are deliberately not exported.

library(tidyverse)
library(magrittr)   # %T>% (tee pipe) used in the PCA normality diagnostic
library(ggpubr)     # ggarrange, stat_compare_means
library(ggridges)   # geom_density_ridges
library(ggnewscale)
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
      labels = c("Negative Site", "Study Site")
    ),
    surface = as.factor(surface),
    surface_location = factor(
        tolower(surface_location),
        levels = c("top", "bottom"),
        labels = c("Top", "Bottom")
    )
  )
write.csv(data, "data/ppe/filtered_data.csv", row.names=FALSE)


# Principal Component Analysis -------------------------------------------
#
# Collapse the three log-metrics into principal components. PC1 captures the
# shared "amount of signal" axis and is used as the summary contamination
# score in the boxplots and density plots below. Scores are appended to the
# data as PC1-PC3.


pca_data <- data %>%
  summarize(across(c(mpr, ms, auc), ~ median(log(.))), .by = c(sample_id, location, timing, surface_location))

pca_res <- prcomp(~ mpr + ms + auc, data = pca_data, scale. = TRUE)
pca_sum <- summary(pca_res)
pca_sum

pca_df <- pca_data %>%
  cbind(pca_res$x)


# PCA score normality (diagnostic, not exported) -------------------------
#
# Histogram of each PC split at zero, plus Shapiro-Wilk per side. The %T>%
# tee draws the plot mid-pipeline while the summary p-values flow through.


pca_df %>%
  pivot_longer(c(PC1, PC2, PC3), names_to = "metric") %>%
  mutate(
    pc = as.integer(str_remove(metric, "PC")),
    metric = paste0(metric, "(", round(100 * pca_sum$importance[2, pc], 1), "%)")
  ) %>%
  mutate(pos = value > 0) %T>%
  {
    print(
      ggplot(., aes(value, color=pos, fill=pos)) +
        geom_histogram(position = "stack", alpha=0.5, bins=50, show.legend = FALSE) +
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
# pca_plot.jpg.


pca_combos <- as.data.frame(t(combn(paste0("PC", seq(length(pca_res$scale))), 2))) %>%
  mutate(
    across(everything(), ~ paste0(., " (", round(100 * pca_sum$importance[2, as.integer(str_remove(., "PC"))], 1), "%)"))
  )
colnames(pca_combos) <- c("x", "y")

rename_pc <- function(cols) {
  for (i in 1:length(cols)) {
    if (str_detect(cols[i], "PC")) {
      pc <- as.integer(str_remove(cols[i], "PC"))
      # pc_format <- round(pc * 100, 2)
      label <- paste0("PC", pc, " (", round(100 * pca_sum$importance[2, pc], 1), "%)")
      cols[i] <- label
    }
  }
  return(cols)
}

make_pca_plot <- function(x, y, alpha=0.5) {
  centers <- 2

  cls <- pca_df %>%
    rename_with(rename_pc) %>%
    select(!!sym(x), !!sym(y)) %>%
    kmeans(centers=centers, trace=TRUE)
  
  df_centers <- cls$centers %>%
    as.data.frame() %>%
    rename(x = 1, y = 2) %>%
    mutate(cluster = paste0("Cluster ", seq(centers)))

  temp_df <- pca_df %>%
    mutate(cluster = as.character(cls$cluster))

  temp_df %>%
    rename_with(rename_pc) %>%
    arrange(desc(location)) %>%
    ggplot(aes(x = !!sym(x), y = !!sym(y))) +
    geom_point(aes(color=location), alpha = alpha, size=3) +
    stat_ellipse(aes(fill=cluster), level=0.95, geom="polygon", alpha=0, type="t", color="black", show.legend = FALSE) +
    # geom_text(aes(x = x, y = y, label = cluster), data = df_centers) +
    scale_color_manual(values=c("darkcyan", "darkorange")) +
    scale_fill_manual(values=c("darkcyan", "darkorange")) +
    facet_wrap(~ timing, ncol=2) +
    labs(
      x = toupper(x),
      y = toupper(y)
    ) +
    main_theme +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom"
    )
}

pca_plots <- pmap(pca_combos, make_pca_plot)
pca_plot <- ggarrange(plotlist = pca_plots, ncol=1, legend="none")
pca_plot


# Boxplot Comparisons ----------------------------------------------------


timing_pvals <- pca_df %>%
  summarize(
    p.value = wilcox.test(PC1 ~ timing)$p.value,
    .by = location
  ) %>%
  mutate(
    label = paste("p =", signif(p.value, 2)),
    xmin  = ifelse(location == "Negative Site", 0.875, 1.125),
    xmax  = ifelse(location == "Negative Site", 1.875, 2.125),
    y     = ifelse(location == "Negative Site", 5.5, 6.5),
    color = ifelse(p.value < 0.05, "darkred", "black")
  )

location_pvals <- pca_df %>%
  summarize(
    p.value = wilcox.test(PC1 ~ location)$p.value,
    .by = timing
  ) %>%
  mutate(
    label = paste("p =", signif(p.value, 2)),
    y     = 9,
    color = ifelse(p.value < 0.05, "darkred", "black")
  )

df_sum <- pca_df %>%
  summarize(
    n = n(), .by = c(location, timing)
  )

dodge_amt <- 0.7

pc_box <- pca_df %>% 
  select(PC1, location, timing) %>% 
  ggplot(aes(timing, PC1, fill=location, color=location)) +
  geom_hline(yintercept = 0, linetype="dashed", color="gray50", linewidth=1) +
  geom_point(position=position_jitterdodge(0.3, 0, dodge_amt)) +
  geom_boxplot(position = position_dodge(dodge_amt), outlier.shape=NA,  color="black", width=0.4) +
  geom_text(# Number of samples
    aes(label = paste("n =", n), group = location, y = 9.5), fontface = "bold",
    data = df_sum, position = position_dodge(dodge_amt), size=5, vjust=0.5
  ) +
  scale_color_manual(values=c("darkcyan", "darkorange")) +
  scale_fill_manual(values=c("darkcyan", "darkorange")) +
  new_scale_color() +
  geom_text(# Location Comparisons
    aes(x = timing, y = y, label = label, color = color),
    data = location_pvals, inherit.aes = FALSE, size = 5, fontface = "bold"
  ) +
  geom_text(# Timing Comparisons
    aes(x = (xmin + xmax) / 2, y = y + 0.3, label = label, color = color),
    data = timing_pvals, inherit.aes = FALSE, size = 5, fontface = "bold"
  ) +
  geom_segment(
    aes(x = xmin, xend = xmax, y = y, yend = y, color = color),
    data = timing_pvals, inherit.aes = FALSE, linewidth = 0.5
  ) +
  geom_segment(
    aes(x = xmin, xend = xmin, y = y - 0.2, yend = y, color = color),
    data = timing_pvals, inherit.aes = FALSE, linewidth = 0.5
  ) +
  geom_segment(
    aes(x = xmax, xend = xmax, y = y - 0.2, yend = y, color = color),
    data = timing_pvals, inherit.aes = FALSE, linewidth = 0.5
  ) +
  scale_color_identity() +
  scale_y_continuous(breaks = seq(-10, 20, 2)) +
  coord_cartesian(ylim=c(-1.5, 10)) +
  labs(
    x = "Timing",
    y = "PC1"
  ) +
  main_theme +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size=24),
    legend.key.size = unit(1, "cm"),
    axis.title.x = element_blank()
  )
pc_box

ggarrange(
  pca_plot, pc_box, ncol = 2, common.legend = TRUE, legend="bottom", labels = c("A", "B"),
  font.label = list(size = 24), widths = c(3,2)
)

ggsave("pca_plot.jpg", path="figures/ppe", width = 12, height = 10, bg="white", dpi=600)


# Metric comparisons by location -----------------------------------------
#
# Logistic smooth of each metric against location (diagnostic, not exported):
# how well each metric separates Core vs. negative site.
pca_df %>%
  mutate(across(c(mpr, ms, auc, PC1, PC2, PC3), ~ rescale(., c(0, 1)))) %>%
  # pivot_longer(c(mpr, ms, auc, PC1, PC2, PC3), names_to="metric") %>%
  mutate(location = ifelse(location == "Study Site", 1, 0)) %>%
  ggplot(aes(PC1, location, color=timing, fill=timing)) +
  geom_point(position=position_jitter(0, 0.01), alpha=0.5) +
  stat_smooth(method="glm", se=TRUE, method.args=list(family="binomial"), fullrange = T) +
  # facet_grid(vars(timing), scales="free_x") +
  main_theme


# ROC Analysis -----------------------------------------------------------
#
# For each metric x timing, build an ROC curve for discriminating Core vs.
# negative site and record the AUC (with CI). The curves are overlaid per
# timing to compare which metric discriminates best. roc.jpg.
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

# Long-format metrics + PCs, each rescaled to 0-1 for a common axis, compared
# between locations within each timing via Wilcoxon tests (significance stars).
# comparisons.jpg.
pca_long <- pca_df %>%
  pivot_longer(c(mpr, ms, auc, PC1, PC2, PC3), names_to="metric") %>%
  mutate(metric = factor(metric, levels=c("mpr", "ms", "auc", "PC1", "PC2", "PC3")))

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
ggsave("roc.jpg", path="figures/ppe", width=12, height=12, dpi=600)


# PC1 Box Plots ----------------------------------------------------------
#
# PC1 (the contamination score) compared across the experimental factors with
# Wilcoxon tests. timing_box and location_box are also stacked into a single
# labelled A/B panel (boxplots.jpg), and PC1 is compared across boot side.

sl_pvals <- pca_df %>%
  summarize(
    p.value = wilcox.test(PC1 ~ surface_location)$p.value,
    .by = c(location, timing)
  ) %>%
  mutate(
    label = paste("p =", signif(p.value, 2)),
    xmin  = ifelse(location == "Negative Site", 0.875, 1.125),
    xmax  = ifelse(location == "Negative Site", 1.875, 2.125),
    y     = ifelse(location == "Negative Site", 5.5, 6.5),
    color = ifelse(p.value < 0.05, "darkred", "black")
  )

location_pvals2 <- pca_df %>%
  na.omit() %>%
  summarize(
    p.value = wilcox.test(PC1 ~ location)$p.value,
    .by = c(timing, surface_location)
  ) %>%
  mutate(
    label = paste("p =", signif(p.value, 2)),
    y     = 9,
    color = ifelse(p.value < 0.05, "darkred", "black")
  )
location_pvals2

# PC1 by boot side (top vs. bottom), faceted by timing. surface-location_box.jpg.
sl_box <- pca_df %>%
  filter(!is.na(surface_location)) %>%
  ggplot(aes(surface_location, PC1, fill=location, color = location)) +
  geom_point(position=position_jitterdodge(0.3, 0, dodge_amt)) +
  geom_boxplot(outliers = FALSE, color="black") +
  scale_fill_manual(values = c("darkcyan", "darkorange")) +
  scale_color_manual(values = c("darkcyan", "darkorange")) +
  new_scale_color() +
  geom_text(# Location Comparisons
    aes(x = surface_location, y = y, label = label, color = color),
    data = location_pvals2, inherit.aes = FALSE, size = 5, fontface = "bold"
  ) +
  geom_text(# Timing Comparisons
    aes(x = (xmin + xmax) / 2, y = y + 0.3, label = label, color = color),
    data = sl_pvals, inherit.aes = FALSE, size = 5, fontface = "bold"
  ) +
  geom_segment(
    aes(x = xmin, xend = xmax, y = y, yend = y, color = color),
    data = sl_pvals, inherit.aes = FALSE, linewidth = 0.5
  ) +
  geom_segment(
    aes(x = xmin, xend = xmin, y = y - 0.2, yend = y, color = color),
    data = sl_pvals, inherit.aes = FALSE, linewidth = 0.5
  ) +
  geom_segment(
    aes(x = xmax, xend = xmax, y = y - 0.2, yend = y, color = color),
    data = sl_pvals, inherit.aes = FALSE, linewidth = 0.5
  ) +
  scale_color_identity() +
  facet_grid(cols = vars(timing)) +
  labs(
    x = "Surface Location"
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  main_theme +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    # axis.title.x = element_blank(),
  )
sl_box
ggsave("surface-location_box.jpg", path = "figures/ppe", width = 8, height = 8, dpi=600)

# Surface Location Histogram ---------------------------------------------

# Determine if values are normally distributed
# Lets us know if we can use t.test or wilcox.test
pc_bottom <- pca_df %>%
  filter(surface_location == "Bottom") %>%
  pull(PC1)
pc_top <- pca_df %>%
  filter(surface_location == "Top") %>%
  pull(PC1)

k_test <- ks.test(pc_bottom, pc_top, simulate.p.value = TRUE)
wilcox_test <- wilcox.test(pc_bottom, pc_top)

sl_hist <- pca_df %>%
  filter(!is.na(surface_location)) %>%
  ggplot(aes(PC1, after_stat(ndensity), fill=surface_location)) +
  stat_bin(position="identity", alpha=0.5, color="black") +
  annotate(
    "text",
    label = paste("Kolmogorov\tp =", signif(k_test$p.value, 2)), 
    x=2.5, y=1, size = 6, vjust = 0.5, hjust=0
  ) +
  annotate(
    "text",
    label = paste("Mann-Whitney\tp =", signif(wilcox_test$p.value, 2)), 
    x=2.5, y=0.95, size = 6, vjust = 0.5, hjust=0
  ) +
  scale_fill_manual(values = c("darkcyan", "darkorange")) +
  scale_color_manual(values = c("darkcyan", "darkorange")) +
  labs(
    x = "PC1",
    y = "Density",
    fill = "Surface Location"
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  main_theme +
  theme(
    # legend.title = element_blank(),
    legend.position = "bottom",
    legend.position.inside = c(0.9, 0.95),
    legend.justification = c(0.5, 0.63),
    legend.background = element_blank()
  )
sl_hist

ggsave("surface-location_histogram.jpg", path = "figures/ppe", width = 12, height = 8, dpi=600)


ggarrange(sl_box, sl_hist, ncol = 2, widths=c(1, 1.5), labels = c("A", "B"), font.label = list(size = 24))
ggsave("surface-location_combined.jpg", path = "figures/ppe", width = 16, height = 12, bg="white", dpi=600)
