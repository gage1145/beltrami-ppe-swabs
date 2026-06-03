library(tidyverse)
library(ggpubr)
library(ggfortify)
library(quicR)
library(ggridges)
library(ggforce)
library(ggpointdensity)
library(readr)
library(emmeans)
library(pROC)


main_theme <- theme(
  plot.title = element_text(size = 30, hjust = 0.5),
  axis.title = element_text(size = 20),
  axis.text = element_text(size = 16),
  strip.text = element_text(size = 20),
  legend.title = element_text(size = 20, hjust = 0.5),
  legend.text = element_text(size = 16),
  # panel.grid = element_blank()
)

filter_outliers <- function(x) {
  filter(
    x,
    value >= (quantile(value, 0.25, na.rm=TRUE) - 1.5 * IQR(value, na.rm=TRUE)),
    value <= (quantile(value, 0.75, na.rm=TRUE) + 1.5 * IQR(value, na.rm=TRUE)),
    .by = c(metric, Location, Timing)
  )
}


# Estimate Plots ----------------------------------------------------------


# Total Results
# results <- read.csv("data/ppe/results.csv", check.names = FALSE)
# results %>%
#   arrange(contrast) %>%
#   mutate(p.value = signif(p.value, 4)) %>%
#   ggplot(
#     aes(
#       fct_inorder(contrast),
#       estimate,
#       ymin = estimate - std.error,
#       ymax = estimate + std.error,
#       color = p.value < 0.05
#     )
#   ) +
#   geom_pointrange(
#     position = position_dodge(width = 0.4), size = 1, linewidth = 1
#   ) +
#   geom_hline(yintercept = 0, linewidth = 1.2, linetype = "dashed") +
#   geom_label(aes(label = p.value), nudge_x = 0.3, size = 3) +
#   facet_grid(term ~ metric, scale = "free", space = "free_y") +
#   scale_colour_manual(values = c("navy", "darkred")) +
#   coord_flip() +
#   theme(
#     axis.title.y = element_blank(),
#     axis.text = element_text(size = 12),
#     axis.title.x = element_text(size = 14, face = "bold"),
#     strip.text = element_text(size = 14, face = "bold"),
#     legend.title = element_text(size = 14),
#     legend.text = element_text(size = 12),
#     legend.position = "bottom"
#   )
# ggsave("results.png", path = "figures/ppe", width = 8, height = 8)

# # Boot Results
# boot_results <- read.csv("data/ppe/boot_results.csv", check.names = FALSE)
# boot_results %>%
#   filter(term != "Surface_Location") %>%
#   arrange(contrast) %>%
#   mutate(p.value = signif(p.value, 4)) %>%
#   ggplot(
#     aes(
#       fct_inorder(contrast),
#       estimate,
#       ymin = estimate - std.error,
#       ymax = estimate + std.error,
#       color = p.value
#     )
#   ) +
#   geom_pointrange(
#     position = position_dodge(width = 0.4), size = 1, linewidth = 1
#   ) +
#   geom_hline(yintercept = 0, linewidth = 1.2, linetype = "dashed") +
#   geom_label(aes(label = p.value), nudge_x = 0.3, size = 3) +
#   facet_grid(term ~ metric, scale = "free", space = "free_y") +
#   scale_colour_gradient(high = "navy", low = "darkred") +
#   coord_flip() +
#   theme(
#     axis.title.y = element_blank(),
#     axis.text = element_text(size = 12),
#     axis.title.x = element_text(size = 14, face = "bold"),
#     strip.text = element_text(size = 14, face = "bold"),
#     legend.title = element_text(size = 14),
#     legend.text = element_text(size = 12),
#     legend.position = "bottom"
#   )
# ggsave("boot_results.png", path = "figures/ppe", width = 6, height = 8)


# Metric Correlations ----------------------------------------------------


data <- read.csv("data/ppe/updated_metadata.csv", check.names = FALSE) %>%
  mutate(
    across(c("Location", "Surface", "Timing"), tolower),
    Timing = str_replace(Timing, " ", "-")
  ) %>%
  filter(
    Location != "perimeter", 
    Timing != "post-decon", 
    Surface %in% c("blank", "boot")
  ) %>%
  mutate(
    Timing = factor(
      Timing,
      levels = c("pre-sampling", "post-sampling", "post-decon"),
      labels = c("Pre-sampling", "Post-sampling", "Post-decon")
    ),
    Location = factor(
      Location,
      levels = c("negative site", "dump site"),
      labels = c("Negative Site", "Core Site")
    ),
    Surface = as.factor(Surface),
    `Surface Location` = factor(
        tolower(`Surface Location`),
        levels = c("top", "bottom"),
        labels = c("Top", "Bottom")
    )
  )

metric_combos <- as.data.frame(t(combn(c("MPR", "MS", "AUC", "QR"), 2)))
colnames(metric_combos) <- c("x", "y")

make_cor_plot <- function(x, y, alpha=0.2) {
  data %>%
    ggplot(aes(x = .data[[x]], y = .data[[y]], color=Location, fill=Location)) +
    # stat_smooth(
    #   alpha = 0.2,
    #   method.args=list(degree=2, span=200, control=loess.control(surface="direct"))
    # ) +
    geom_point(alpha = alpha) +
    # stat_cor(
    #   aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~','~")),
    #   label.y.npc = 1, size = 6
    # ) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values=c("darkcyan", "darkorange")) +
    scale_fill_manual(values=c("darkcyan", "darkorange")) +
    labs(
      x = toupper(x),
      y = toupper(y)
    ) +
    main_theme
}

cor_plots <- pmap(metric_combos, make_cor_plot, alpha=0.75)
ggarrange(plotlist = cor_plots, align = "hv", common.legend = TRUE)
ggsave("corplot.png", path = "figures/ppe", width = 12, height = 8)

data %>%
  pivot_longer(c(MPR, MS, AUC, QR), names_to = "metric") %>%
  # filter_outliers() %>%
  ggplot(aes(value, Location, fill=Timing)) +
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
ggsave("distributions.png", path="figures/ppe", width=16, height=12)


# Principal Component Analysis -------------------------------------------


pca_data <- data %>%
  filter(!is.na(MPR))

pca_res <- prcomp(~ MPR + MS + AUC + QR, data = pca_data, scale. = TRUE)
summary(pca_res)

make_pca_plot <- function(x, y, alpha=0.5) {
  autoplot(
    pca_res, x = x, y = y, data = pca_data, 
    color = "Location", alpha = alpha,
    loadings = TRUE, loadings.color="blue", loadings.label = TRUE
  ) +
    facet_grid(cols=vars(Timing)) +
    scale_color_manual(values=c("darkcyan", "darkorange")) +
    main_theme
}

pca_combos <- as.data.frame(t(combn(seq(length(pca_res$scale)), 2)))
colnames(pca_combos) <- c("x", "y")

pca_plots <- pmap(pca_combos, make_pca_plot)
ggarrange(plotlist = pca_plots, align = "hv", common.legend = TRUE)

pca_df <- pca_data %>%
  cbind(pca_res$x)

pca_long <- pca_df %>%
  pivot_longer(c(MPR, MS, AUC, QR, PC1, PC2, PC3, PC4), names_to="metric") %>%
  mutate(metric = factor(metric, levels=c("MPR", "MS", "AUC", "QR", "PC1", "PC2", "PC3", "PC4")))

xval <- "Location"
colval <- "Timing"

pca_long %>%
  # filter_outliers() %>%
  ggplot(aes(!!sym(xval), value, color=!!sym(colval), fill=!!sym(colval))) +
  geom_point(position=position_jitterdodge(0.4, 0, 0.7), alpha=0.3) +
  geom_boxplot(outlier.color=NA, color="black", width=0.6, position=position_dodge(0.7)) +
  facet_wrap(vars(metric), scales="free_y", ncol=4) +
  stat_compare_means(
    aes(label = after_stat(p.signif)),
    method="wilcox.test", size=8, vjust=0.5, show.legend=FALSE
  ) +
  scale_fill_manual(values = c("darkcyan", "darkorange")) +
  scale_color_manual(values = c("darkcyan", "darkorange")) +
  main_theme +
  theme(
    legend.position = "top",
    axis.title = element_blank()
  )
ggsave("comparisons.png", path="figures/ppe", width=16, height=10)

pca_long %>%
  filter_outliers() %>%
  mutate(Location = ifelse(Location == "Core Site", 1, 0)) %>%
  ggplot(aes(value, Location, color=Timing, fill=Timing)) +
  geom_point(position=position_jitter(0, 0.01), alpha=0.1) + 
  stat_smooth(method="glm", se=FALSE, method.args=list(family="binomial")) +
  facet_wrap(vars(metric), ncol=4, scales="free_x") +
  scale_fill_manual(values = c("darkcyan", "darkorange")) +
  scale_color_manual(values = c("darkcyan", "darkorange")) +
  main_theme



# ROC Analysis -----------------------------------------------------------


get_rocs <- function(m, df) {
  temp_df <- df %>%
    filter(
      Timing == "Post-sampling",
      metric == m
    ) %>%
    roc(Location, value)
  coords(temp_df) %>%
    mutate(
      metric = m,
      auc = ci.auc(temp_df)[[2]]
    )
}

metrics <- unique(pca_long$metric)
roc_df <- map_dfr(metrics, get_rocs, df=pca_long)

auc_df <- roc_df %>%
  summarize(auc = signif(mean(auc), 3), .by = metric) %>%
  mutate(
    y = c(0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15),
    label = paste(metric, "=", auc)
  )

roc_df %>%
  arrange(auc, sensitivity) %>%
  ggplot(aes(specificity, sensitivity, color=metric)) +
  geom_step(linewidth=1.2) +
  geom_abline(
    slope = 1, intercept = 1, linetype = "dashed",
    color = "gray50", inherit.aes = FALSE
  ) +
  geom_text(aes(0.2, y, label=label, color=metric), data=auc_df, size=8, hjust=0, inherit.aes=FALSE, show.legend=FALSE) +
  scale_x_reverse(
    limits = c(0, 1), breaks = seq(0, 1, 0.2),
    # sec.axis = sec_axis(~., name = "Log Dilution Factor", breaks = NULL, labels = NULL)
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    x = "Specificity", y = "Sensitivity"
  ) +
  main_theme +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
  )

ggsave("roc.png", path="figures/ppe", width=12, height=12)


# Box Plots --------------------------------------------------------------


box_theme <- theme(
  plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
  axis.title.x = element_blank(),
  axis.title.y = element_text(size = 14),
  axis.text = element_text(size = 12),
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 12),
  legend.position = "bottom",
  strip.text = element_text(size = 14)
)

density_theme <- theme(
  plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
  axis.title.x = element_text(size = 14),
  axis.title.y = element_blank(),
  axis.text = element_text(size = 12),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 12),
  legend.position = "bottom"
)


# Box Plots ---------------------------------------------------------------


# comp <- wilcox.test()
smpl_comps <- list(c("Pre-sampling", "Post-sampling"))
site_comps <- list(c("Negative Site", "Core Site"))
sl_comps   <- list(c("Bottom", "Top"))
pc1_scalar <- abs(min(pca_df$PC1)) + 1

pca_prepped <- pca_df %>%
  mutate(PC1 = log(PC1 + pc1_scalar))

# Timing Boxplot
pca_prepped %>%
  ggplot(aes(Timing, PC1, fill = Timing)) +
  geom_boxplot(show.legend=FALSE, outliers = FALSE) +
  geom_dotplot(
    show.legend=FALSE, binaxis="y", stackdir="center", 
    dotsize=0.2, stackratio=0.3, color="#333", fill="#333",
    position=position_jitter(0.01, 0.1), alpha=0.25
  ) +
  stat_compare_means(
    aes(group = Timing),
    comparisons = smpl_comps,
    method = "wilcox.test",
    label = "p.signif", size = 6, show.legend = FALSE
  ) +
  facet_grid(cols=vars(Location)) +
  scale_fill_manual(values = c("darkcyan", "darkorange", "maroon")) +
  scale_y_continuous(limits = c(0, 3.1)) +
  labs(
    y = sprintf("Log(PC1 + %s)", round(pc1_scalar, 2)),
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  main_theme +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    axis.title.x = element_blank(),
  )
ggsave("timing_box.png", path = "figures/ppe", width = 12, height = 8)

# Location Boxplot
pca_prepped %>%
  ggplot(aes(Location, PC1, fill = Timing, group = Location)) +
  geom_boxplot(show.legend=FALSE, outliers = FALSE) +
  geom_dotplot(
    show.legend=FALSE, binaxis="y", stackdir="center", 
    dotsize=0.2, stackratio=0.3, color="#333", fill="#333",
    position=position_jitter(0.01, 0.1), alpha=0.25
  ) +
  stat_compare_means(
    label.y = c(3, 3.4, 3.2),
    method = "wilcox.test", method.args = list(alternative = "less"),
    comparisons = site_comps, label = "p.signif", ref.group = "negative site",
    bracket.size = 0.5, size = 6, show.legend = FALSE,
  ) +
  scale_y_continuous(limits = c(0, 3.7)) +
  scale_fill_manual(values = c("darkcyan", "darkorange", "maroon")) +
  facet_grid(cols = vars(Timing)) +
  labs(
    y = sprintf("Log(PC1 + %s)", round(pc1_scalar, 2)),
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  main_theme +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    axis.title.x = element_blank(),
  )
ggsave("location_box.png", path = "figures/ppe", width = 12, height = 8)

# Side Boxplot
pca_prepped %>%
  filter(!is.na(`Surface Location`)) %>%
  ggplot(aes(`Surface Location`, PC1, fill = Timing, group=`Surface Location`)) +
  geom_boxplot(show.legend=FALSE, outliers = FALSE) +
  geom_dotplot(
    show.legend=FALSE, binaxis="y", stackdir="center", 
    dotsize=0.2, stackratio=0.3, color="#333", fill="#333",
    position=position_jitter(0.01, 0.1), alpha=0.25
  ) +
  stat_compare_means(
    label.y = c(3, 3.4, 3.2),
    method = "wilcox.test",
    comparisons = sl_comps, label = "p.signif",
    bracket.size = 0.5, size = 6, show.legend = FALSE,
  ) +
  scale_y_continuous(limits = c(0, 3.7)) +
  scale_fill_manual(values = c("darkcyan", "darkorange")) +
  facet_grid(cols = vars(Timing)) +
  labs(
    y = sprintf("Log(PC1 + %s)", round(pc1_scalar, 2)),
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  main_theme +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    axis.title.x = element_blank(),
  )
ggsave("surface-location_box.png", path = "figures/ppe", width = 12, height = 8)


# Density Plots -----------------------------------------------------------


pca_prepped %>%
  ggplot(aes(PC1, Location, fill = Timing)) +
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


# Boot figures ------------------------------------------------------------


# boot_data <- read.csv("data/boot_analyzed.csv", check.names = FALSE)
results <- read.csv("data/ppe/boot_results.csv", check.names = FALSE)


results %>%
  arrange(contrast) %>%
  mutate(
    term = ifelse(term == "Surface_Location", "Surface", term),
    p.value = signif(p.value, 4)
  ) %>%
  ggplot(
    aes(
      fct_inorder(contrast),
      estimate,
      ymin = estimate - std.error,
      ymax = estimate + std.error,
      color = p.value < 0.05
    )
  ) +
  geom_pointrange(
    position = position_dodge(width = 0.4), size = 1, linewidth = 1
  ) +
  geom_hline(yintercept = 0, linewidth = 1.2, linetype = "dashed") +
  geom_label(aes(label = p.value), nudge_x = 0.3, size = 3) +
  facet_grid(term ~ metric, scale = "free", space = "free_y") +
  scale_colour_manual(values = c("navy", "darkred")) +
  coord_flip() +
  theme(
    axis.title.y = element_blank(),
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  )
ggsave("boot_results.png", path = "figures/ppe", width = 8, height = 8)


# Plate Views -------------------------------------------------------------


for (file in list.files(c("raw/ppe", "raw/spike_raw"), recursive = TRUE, full.names = TRUE)) {
  file_list <- str_split(file, "/")[[1]]
  depth <- length(file_list)
  file_name <- file_list[depth] %>%
    str_remove(".xlsx")

  save_name <- paste0(file_name, "_plate_view.png")

  if (file.exists(paste0("figures/plate_views/", save_name))) {
    print(paste("Plate view for", file_name, "already exists."))
    next
  }
  print(paste("Working on reaction:", file_name))

  plate_view(get_quic(file), plot_deriv = FALSE) +
    ggtitle(file_name)

  ggsave(save_name, path = "figures/plate_views", width = 12, height = 8)
}
  