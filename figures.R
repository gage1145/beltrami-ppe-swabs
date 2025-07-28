library(dplyr)
library(ggplot2)
# library(quicR)
library(ggridges)
library(readr)
library(forcats)




# Estimate Plots ----------------------------------------------------------



# Total Results
results <- read.csv("data/results.csv", check.names = FALSE)
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
boot_results <- read.csv("data/boot_results.csv", check.names = FALSE)
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
  # scale_y_continuous(limits = c(-1, 1)) +
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



data <- read.csv("data/updated_metadata.csv", check.names = FALSE) %>%
  filter(!is.na(Technician), `Sample Type` == "PPE Swab") %>%
  mutate(
    Timing = factor(
      Timing, levels = c("pre sampling", "post sampling", "post decon")
    ),
    Location = factor(
      Location, levels = c("negative site", "perimeter", "dump site")
    ),
    Surface = factor(
      Surface, levels = unique(Surface)
    )
  )
results <- read.csv("data/results.csv", check.names = FALSE)

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
    # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  )



# Box Plots ---------------------------------------------------------------



data %>%
  filter(Surface %in% c("blank", "boot")) %>%
  ggplot(aes(Timing, RAF, color = Location)) +
  geom_boxplot(fill = NA, size = 1) +
  facet_grid(cols = vars(Surface)) +
  scale_color_manual(values = c("darkcyan", "darkorange", "maroon")) +
  scale_y_continuous(limits = c(0.013, 0.04)) +
  ggtitle("Surface RAFs by Location") +
  box_theme
ggsave("TLS_raf.png", path = "figures", width = 6, height = 8)

data %>%
  ggplot(aes(Timing, RAF, color = Surface)) +
  geom_boxplot(fill = NA, size = 1) +
  facet_grid(cols = vars(Location)) +
  scale_color_manual(values = c("darkcyan", "cyan3", "forestgreen", "darkorange", "maroon")) +
  ggtitle("Surface RAFs by Location") +
  scale_y_continuous(limits=c(0.013,0.04)) +
  box_theme
ggsave("TSL_raf.png", path = "figures", width = 14, height = 8)

data %>%
  ggplot(aes(Surface, RAF, color = Timing)) +
  geom_boxplot(fill = NA, size = 1) +
  facet_grid(cols = vars(Location)) +
  scale_color_manual(values = c("darkcyan", "darkorange", "maroon")) +
  ggtitle("Surface RAFs by Location") +
  scale_y_continuous(limits=c(0.013,0.04)) +
  box_theme
ggsave("STL_raf.png", path = "figures", width = 14, height = 8)

data %>%
  filter(Surface %in% c("blank", "boot")) %>%
  filter(Timing != "post decon") %>%
  ggplot(aes(Location, RAF, color = Timing)) +
  geom_boxplot(fill = NA, size = 1) +
  facet_grid(cols = vars(Surface)) +
  scale_color_manual(values = c("darkcyan", "darkorange", "maroon")) +
  ggtitle("Surface RAFs by Location") +
  scale_y_continuous(limits=c(0.013,0.04)) +
  box_theme
ggsave("LTS_raf.png", path = "figures", width = 6, height = 8)



# Density Plots -----------------------------------------------------------



# Surface Density Plots
data %>%
  arrange(Location) %>%
  filter(Timing != "post decon") %>%
  select(MPR, MS, Location, Timing) %>%
  gather("key", "value", -c(Location, Timing)) %>%
  ggplot(aes(value, Location, fill = Timing)) +
  geom_density_ridges(alpha = 0.5, scale = 1.2) +
  facet_wrap(vars(key), nrow=2, scales = "free", strip.position = "right") +
  scale_fill_manual(values = c("darkcyan", "darkred", "darkorange")) +
  scale_y_discrete(expand = expansion(mult = c(0.025,0.5))) +
  ggtitle("Density Plots of Surface MPRs by Sample Timing") +
  density_theme
ggsave("location_timing_density.png", path = "figures", width = 14, height = 8)


# Surface Density Plots
data %>%
  filter(Timing != "post decon") %>%
  arrange(Surface) %>%
  select(MPR, MS, Surface, Location, Timing) %>%
  gather("key", "value", -c(Surface, Location, Timing)) %>%
  ggplot(aes(value, fct_rev(Surface), fill = Location, linetype = Timing)) +
  geom_density_ridges(alpha = 0.5, linewidth = 0.75, scale = 1.2) +
  facet_wrap(vars(key), nrow=2, scales = "free", strip.position = "right") +
  scale_fill_manual(values = c("darkcyan", "darkred", "darkorange")) +
  scale_y_discrete(expand = expansion(mult = c(0.025,0.35))) +
  # scale_x_continuous(limits = c(0, 10)) +
  ggtitle("Density Plots of Surface MPRs by Location & Timing") +
  density_theme
ggsave(
  "surface_location_timing_density.png",
  path = "figures",
  width = 14,
  height = 8
)





# Plate Views -------------------------------------------------------------



for (file in list.files("raw", "spike_raw")) {
  
  file_name <- paste0(str_remove(file, ".xlsx"), "_plate_view.png")
  
  if (file.exists(paste0("figures/plate_views/", file_name))) {
    print(paste(file_name, "already exists."))
    next
  }
  
  plate_view(
    get_real(paste0("raw/", file)),
    get_sample_locations(paste0("raw/", file))
  ) +
    ggtitle(file)
  
  ggsave(file_name, path = "figures/plate_views", width = 12, height = 8)
  
}

