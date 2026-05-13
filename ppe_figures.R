library(dplyr)
library(ggplot2)
library(ggpubr)
library(quicR)
library(ggridges)
library(ggpointdensity)
library(readr)
library(forcats)
library(stringr)


main_theme <- theme(
  plot.title = element_text(size = 30, hjust = 0.5),
  axis.title = element_text(size = 20),
  axis.text = element_text(size = 16),
  strip.text = element_text(size = 20),
  legend.title = element_text(size = 20, hjust = 0.5),
  legend.text = element_text(size = 16),
  panel.grid = element_blank()
)


# Estimate Plots ----------------------------------------------------------



# Total Results
results <- read.csv("data/ppe/results.csv", check.names = FALSE)
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
ggsave("results.png", path = "figures/ppe", width = 8, height = 8)

# Boot Results
boot_results <- read.csv("data/ppe/boot_results.csv", check.names = FALSE)
boot_results %>%
  filter(term != "Surface_Location") %>%
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
  geom_label(aes(label = p.value), nudge_x = 0.3, size = 3) +
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
ggsave("boot_results.png", path = "figures/ppe", width = 6, height = 8)



data <- read.csv("data/ppe/updated_metadata.csv", check.names = FALSE) %>%
  mutate_at(c("Location", "Surface", "Timing"), tolower) %>%
  mutate(Timing = str_replace(Timing, " ", "-")) %>%
  filter(!is.na(Technician), `Sample Type` == "PPE Swab") %>%
  mutate(
    Timing = factor(
      Timing, levels = c("pre-sampling", "post-sampling", "post-decon")
    ),
    Location = factor(
      Location, levels = c("negative site", "perimeter", "dump site")
    ),
    Surface = factor(
      Surface, levels = unique(Surface)
    )
  )
# results <- read.csv("data/ppe/results.csv", check.names = FALSE)



# Principal Component Analysis -------------------------------------------


pca_data <- data %>%
  filter(!is.na(MPR)) %>%
  filter(Surface == "boot")

pca_res <- pca_data %>%
  select(MPR, MS, AUC) %>%
#   mutate(time_since_start = as.numeric(time_since_start)) %>%
  prcomp(scale. = TRUE)
summary(pca_res)

pca_df <- pca_res$x %>%
  as.data.frame() %>%
  mutate(
    Zone = pca_data$Zone,
    Location = pca_data$Location,
    Timing = pca_data$Timing,
    Surface = pca_data$Surface,
    `Surface Location` = pca_data$`Surface Location`
  )

pca_long <- pca_df %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(starts_with("PC"), names_to = "PC", values_to = "value")

pca_comp_df <- inner_join(
  pca_long %>% rename(x_pc = PC, x_val = value),
  pca_long %>% rename(y_pc = PC, y_val = value),
  relationship = "many-to-many"
) %>%
  filter(x_pc != y_pc)

pca_comp_df %>%
  ggplot(aes(x_val, y_val, color = `Surface Location`)) +
  geom_point() +
#   stat_ellipse(geom="polygon", level=0.95, alpha=0.2) +
  facet_grid(rows = vars(y_pc), cols = vars(x_pc)) +
  # scale_color_manual(values=c("cyan", "darkorange", "darkgreen")) +
  coord_fixed() +
  labs(
    x = NULL, y = NULL,
    title = "Principal Component Pairs"
  ) +
  main_theme
ggsave("pca.png", path="figures/ppe", width=12, height=12)


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

comps <- list(c("negative site", "perimeter"), c("negative site", "dump site"), c("perimeter", "dump site"))
pc1_scalar <- abs(min(pca_df$PC1)) + 1

pca_prepped <- pca_df %>%
  mutate(PC1 = log(PC1 + pc1_scalar)) %>%
  filter(Surface %in% c("blank", "boot")) %>%
  filter(Timing != "post-decon") 

# Timing Boxplot
pca_prepped %>%
  ggplot(aes(Location, PC1, fill=Timing)) +
  geom_boxplot() +
  stat_compare_means(
    color = c("red", "red", "black"),
    method = "wilcox.test", method.args = list(alternative = "greater"),
    label.y.npc=0.95,
    label = "p", size=6, show.legend = FALSE
  ) +
  scale_fill_manual(values = c("darkcyan", "darkorange", "maroon")) +
  scale_y_continuous(limits=c(0, 3.7)) +
  labs(
    y = sprintf("Log(PC1 + %s)", round(min_pc1, 2)),
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  main_theme +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    axis.title.x = element_blank()
  )
ggsave("timing_box.png", path = "figures/ppe", width = 6, height = 8)

# Location Boxplot
pca_prepped %>%
  ggplot(aes(Location, PC1, fill=Timing, group=Location)) +
  geom_boxplot() +
  stat_compare_means(
    aes(color = ifelse(after_stat(p.value) < 0.05, "red", "black")),
    label.y = c(3, 3.4, 3.2),
    method = "wilcox.test", method.args = list(alternative = "less"),
    comparisons = comps, label = "p", ref.group = "negative site",
    bracket.size = 0.5, size=6, show.legend = FALSE,
  ) +
  scale_y_continuous(limits=c(0, 3.7)) +
  scale_fill_manual(values = c("darkcyan", "darkorange", "maroon")) +
  facet_grid(cols = vars(Timing)) +
  labs(
    y = sprintf("Log(PC1 + %s)", round(min_pc1, 2)),
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  main_theme +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    axis.title.x = element_blank(),
    strip.background = element_blank(), 
    strip.text = element_blank()
  )
ggsave("location_box.png", path = "figures/ppe", width = 10, height = 8)


data %>%
  filter(Surface %in% c("blank", "boot")) %>%
  # filter(Timing != "post-decon") %>%
  filter(Timing == "post-sampling") %>%
  # group_by(Location, Surface) %>%
  # summarize(RAF = mean(RAF)) 
  # ggplot(aes(Location, MPR, color = Timing)) +
  ggplot(aes(Location, RAF, color = Surface)) +
  geom_boxplot(fill = NA, size = 1) +
  # geom_col() +
  # stat_compare_means(
  #   method = "anova",
  #   label = "p.signif",
  #   # label.y = 12.5,
  #   tip.length = 0.01,
  #   size=7,
  #   show.legend = FALSE
  # ) +
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
  # facet_grid(cols = vars(Surface)) +
  scale_color_manual(values = c("darkcyan", "darkorange", "maroon")) +
  ggtitle("Surface RAFs by Location") +
  scale_y_continuous(breaks=seq(0.01, 0.05, 0.01)) +
  coord_cartesian(ylim=c(min(data$RAF), 0.05), xlim=c(0.5, 3.5), expand=FALSE) +
  box_theme +
  theme(
    axis.text.x = element_text(angle=0, hjust=0.5, vjust=1),
    plot.title = element_blank()
  )
ggsave("LTS_raf.png", path = "figures/ppe", width = 6, height = 8)

data %>%
  filter(Surface %in% c("blank", "boot")) %>%
  # filter(Timing != "post-decon") %>%
  filter(Timing == "post-sampling") %>%
  # group_by(Location, Surface) %>%
  # summarize(RAF = mean(RAF)) 
  # ggplot(aes(Location, MPR, color = Timing)) +
  ggplot(aes(Location, MPR, color = Surface)) +
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
  # facet_grid(cols = vars(Surface)) +
  scale_color_manual(values = c("darkcyan", "darkorange", "maroon")) +
  ggtitle("Surface RAFs by Location") +
  # scale_y_continuous(breaks=seq(0.01, 0.05, 0.01)) +
  # coord_cartesian(ylim=c(min(data$RAF), 0.05), xlim=c(0.5, 3.5), expand=FALSE) +
  box_theme +
  theme(
    axis.text.x = element_text(angle=0, hjust=0.5, vjust=1),
    plot.title = element_blank()
  )
ggsave("LTS_mpr.png", path = "figures/ppe", width = 6, height = 8)



# Density Plots -----------------------------------------------------------



# Surface Density Plots
data %>%
  arrange(Location) %>%
  filter(Timing != "post-decon") %>%
  filter(Surface == "boot") %>%
  select(MPR, Location, Timing) %>%
  gather("key", "value", -c(Location, Timing)) %>%
  ggplot(aes(value, Location, fill = Timing)) +
  geom_density_ridges(alpha = 0.5, scale = 1.2) +
  facet_wrap(vars(key), nrow=2, scales = "free", strip.position = "right") +
  scale_fill_manual(values = c("darkcyan", "darkred", "darkorange")) +
  scale_y_discrete(expand = expansion(mult = c(0.025,0.5))) +
  scale_x_continuous(limits=c(0,10)) +
  ggtitle("Density Plots of Boot MPRs by Sample Timing") +
  density_theme
ggsave("location_timing_density.png", path = "figures/ppe", width = 14, height = 8)



# Boot figures ------------------------------------------------------------



# boot_data <- read.csv("data/boot_analyzed.csv", check.names = FALSE)
results <- read.csv("data/ppe/boot_results.csv", check.names = FALSE)


results %>%
  arrange(contrast) %>%
  mutate(
    term = ifelse(term == "Surface_Location", "Surface", term),
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
ggsave("boot_results.png", path = "figures/ppe", width = 8, height = 8)





# Plate Views -------------------------------------------------------------



for (file in list.files(c("raw/ppe", "raw/spike_raw"), recursive=TRUE, full.names=TRUE)) {

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

