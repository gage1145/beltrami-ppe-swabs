library(dplyr)
library(ggplot2)
library(ggpubr)
library(quicR)
library(ggridges)
library(readr)
library(forcats)
library(stringr)
library(ggh4x)


# Plot themes ------------------------------------------------------------

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

# Load in data -----------------------------------------------------------



results <- read.csv("data/results.csv", check.names = FALSE)
boot_results <- read.csv("data/boot_results.csv", check.names = FALSE)

data <- read.csv("data/analyzed.csv", check.names = FALSE) 



# Estimate Plots ----------------------------------------------------------



# Total Results
results %>%
  arrange(contrast) %>%
  mutate(p.value = signif(p.value, 4)) %>%
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
ggsave("results.png", path = "figures", width = 8, height = 8)

# Boot Results
boot_results %>%
  filter(term != "surface_location") %>%
  arrange(contrast) %>%
  mutate(p.value = signif(p.value, 4)) %>%
  ggplot(
    aes(
      fct_inorder(contrast), 
      estimate, 
      ymin = estimate - std.error, 
      ymax = estimate + std.error,
      color = p.value
    )
  ) +
  geom_pointrange(
    position = position_dodge(width = 0.4), size = 1, linewidth = 1
  ) +
  geom_hline(yintercept = 0, linewidth = 1.2, linetype = "dashed") +
  geom_label(aes(label = paste("p =", p.value)), nudge_x = 0.3, size = 3) +
  facet_grid(term ~ metric, scale = "free", space = "free_y") +
  scale_colour_gradient(high="navy", low="darkred") +
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
ggsave("boot_results.png", path = "figures", width = 6, height = 8)



# Box Plots ---------------------------------------------------------------



# data %>%
#   ggplot(aes(surface, raf, color = timing)) +
#   geom_boxplot(fill = NA, size = 1) +
#   facet_grid(cols = vars(location)) +
#   # scale_color_manual(values = c("darkcyan", "darkorange", "maroon")) +
#   ggtitle("surface rafs by location") +
#   scale_y_continuous(limits=c(0.013,0.04)) +
#   box_theme
# ggsave("STL_raf.png", path = "figures", width = 14, height = 8)

data %>%
  pivot_longer(c(mpr, ms, auc), names_to="metric") %>%
  filter(surface %in% c("blank", "boot")) %>%
  filter(timing == "post-sampling") %>%
  ggplot(aes(location, value, color = surface)) +
  geom_boxplot(fill = NA, size = 1, outliers = FALSE) +
  stat_compare_means(
    method = "wilcox.test",
    method.args = list(alternative = "greater"),
    label = "p.signif",
    label.y.npc = 0.9,
    tip.length = 0.01,
    size=6,
    show.legend = FALSE
  ) +
  # stat_compare_means(
  #   method = "t.test",
  #   comparisons = list(
  #     c("negative site", "dump site"), 
  #     c("negative site", "perimeter"), 
  #     c("dump site", "perimeter")
  #   ),
  #   label = "p.signif",
  #   label.y.npc = 0.8,
  #   tip.length = 0.01,
  #   size=7,
  #   show.legend = FALSE,
  #   alpha = 0
  # ) +
  facet_wrap(vars(metric), scales="free") +
  facetted_pos_scales(
    y = list(
      metric == "mpr" ~ scale_y_continuous(limits = c(0, 12), expand = expansion(0)),
      metric == "ms"  ~ scale_y_continuous(limits = c(0, 2.5), expand = expansion(0)),
      metric == "auc" ~ scale_y_continuous(limits = c(50, 300), expand = expansion(0)),
      metric == "raf" ~ scale_y_continuous(limits = c(0.013, 0.05), expand = expansion(0))
    )
  ) +
  scale_color_manual(values = c("darkcyan", "darkorange", "maroon")) +
  box_theme +
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle=0, hjust=0.5, vjust=1),
    plot.title = element_blank()
  )
ggsave("LTS.png", path = "figures", width = 12, height = 8)

data %>%
  filter(surface %in% c("blank", "boot")) %>%
  # filter(timing != "post-decon") %>%
  filter(timing == "post-sampling") %>%
  # group_by(location, surface) %>%
  # summarize(raf = mean(raf)) 
  # ggplot(aes(location, mpr, color = timing)) +
  ggplot(aes(location, mpr, color = surface)) +
  geom_boxplot(fill = NA, size = 1) +
  # geom_col() +
  stat_compare_means(
    method = "anova",
    label = "p",
    label.y = 16,
    tip.length = 0.01,
    size=5,
    show.legend = FALSE
  ) +
  # stat_compare_means(
  #   method = "t.test",
  #   comparisons = list(
  #     c("negative site", "dump site"), 
  #     c("negative site", "perimeter"), 
  #     c("dump site", "perimeter")
  #   ),
  #   label = "p.signif",
  #   label.y = c(10,9,9.5),
  #   tip.length = 0.01,
  #   size=7,
  #   show.legend = FALSE
  # ) +
  # facet_grid(cols = vars(surface)) +
  scale_color_manual(values = c("darkcyan", "darkorange", "maroon")) +
  ggtitle("surface rafs by location") +
  # scale_y_continuous(breaks=seq(0.01, 0.05, 0.01)) +
  # coord_cartesian(ylim=c(min(data$raf), 0.05), xlim=c(0.5, 3.5), expand=FALSE) +
  box_theme +
  theme(
    axis.text.x = element_text(angle=0, hjust=0.5, vjust=1),
    plot.title = element_blank()
  )
ggsave("LTS_mpr.png", path = "figures", width = 6, height = 8)



# Density Plots -----------------------------------------------------------



# surface Density Plots
data %>%
  arrange(location) %>%
  filter(timing != "post-decon") %>%
  filter(surface == "boot") %>%
  select(mpr, location, timing) %>%
  gather("key", "value", -c(location, timing)) %>%
  ggplot(aes(value, location, fill = timing)) +
  geom_density_ridges(alpha = 0.5, scale = 1.2) +
  facet_wrap(vars(key), nrow=2, scales = "free", strip.position = "right") +
  scale_fill_manual(values = c("darkcyan", "darkred", "darkorange")) +
  scale_y_discrete(expand = expansion(mult = c(0.025,0.5))) +
  scale_x_continuous(limits=c(0,10)) +
  ggtitle("Density Plots of Boot mprs by Sample timing") +
  density_theme
ggsave("location_timing_density.png", path = "figures", width = 14, height = 8)


# surface Density Plots
data %>%
  filter(timing != "post decon") %>%
  arrange(surface) %>%
  select(mpr, surface, location, timing) %>%
  gather("key", "value", -c(surface, location, timing)) %>%
  ggplot(aes(value, location, fill = surface)) +
  geom_density_ridges(alpha = 0.5, linewidth = 0.75, scale = 1.2) +
  facet_wrap(vars(key), nrow=2, scales = "free", strip.position = "right") +
  scale_fill_manual(values = c("darkcyan", "darkred", "darkorange", "cyan")) +
  scale_y_discrete(expand = expansion(mult = c(0.025,0.35))) +
  scale_x_continuous(limits=c(0,10)) +
  # scale_x_continuous(limits = c(0, 10)) +
  ggtitle("Density Plots of surface mprs by location & timing") +
  density_theme
ggsave(
  "surface_location_timing_density.png",
  path = "figures",
  width = 14,
  height = 8
)



# Boot figures ------------------------------------------------------------



# boot_data <- read.csv("data/boot_analyzed.csv", check.names = FALSE)
results <- read.csv("data/boot_results.csv", check.names = FALSE)


results %>%
  arrange(contrast) %>%
  mutate(
    term = ifelse(term == "surface_location", "surface", term),
    p.value = signif(p.value, 4)) %>%
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
ggsave("boot_results.png", path = "figures", width = 8, height = 8)





# Plate Views -------------------------------------------------------------



for (file in list.files(c("raw", "spike_raw"), recursive=TRUE, full.names=TRUE)) {

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

  plate_view(get_quic(file), plot_deriv=FALSE) +
    ggtitle(file_name)

  ggsave(save_name, path = "figures/plate_views", width = 12, height = 8)

}

