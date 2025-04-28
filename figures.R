library(dplyr)
library(ggplot2)
# library(quicR)
library(ggridges)
library(readr)




data <- read.csv("data/updated_metadata.csv", check.names = FALSE) %>%
  filter(!is.na(Technician)) %>%
  mutate(
    Timing = factor(
      Timing, levels = c("pre sampling", "post sampling", "post decon")
    ),
    Location = factor(
      Location, levels = c("negative site", "dump site", "perimeter")
    ),
    Surface = factor(
      Surface, levels = c("blank sample", "tyvek leg", "plastic boot", "boot side", "boot bottom")
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
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  )



# Box Plots ---------------------------------------------------------------



data %>%
  ggplot(aes(Timing, MPR, color = Location)) +
  geom_boxplot(fill = NA, size = 1) +
  facet_grid(cols = vars(Surface)) +
  scale_color_manual(values = c("darkcyan", "darkorange", "maroon")) +
  # ggtitle("Surface RAFs by Location") +
  # scale_y_continuous(limits=c(1,6)) +
  box_theme

data %>%
  ggplot(aes(Timing, MPR, color = Surface)) +
  geom_boxplot(fill = NA, size = 1) +
  facet_grid(cols = vars(Location)) +
  scale_color_manual(values = c("darkcyan", "cyan", "forestgreen", "darkorange", "maroon")) +
  ggtitle("Surface RAFs by Location") +
  scale_y_continuous(limits=c(1,6)) +
  box_theme

data %>%
  ggplot(aes(Surface, MPR, color = Timing)) +
  geom_boxplot(fill = NA, size = 1) +
  facet_grid(cols = vars(Location)) +
  scale_color_manual(values = c("darkcyan", "darkorange", "maroon")) +
  ggtitle("Surface RAFs by Location") +
  scale_y_continuous(limits=c(1,6)) +
  box_theme

data %>%
  ggplot(aes(Location, MPR, color = Timing)) +
  geom_boxplot(fill = NA, size = 1) +
  facet_grid(cols = vars(Surface)) +
  scale_color_manual(values = c("darkcyan", "darkorange", "maroon")) +
  ggtitle("Surface RAFs by Location") +
  scale_y_continuous(limits=c(0,7)) +
  box_theme



# Density Plots -----------------------------------------------------------



# Surface Density Plots
data %>%
  ggplot(aes(RAF, Timing, fill = Location)) +
  geom_density_ridges(alpha = 0.5) +
  scale_fill_manual(values = c("forestgreen", "darkorange", "maroon")) +
  coord_cartesian(expand = FALSE) +
  # scale_x_continuous(limits = c(1, 6)) +
  # facet_grid(cols = vars(Location), scales = "free_x", space = "free_x") +
  ggtitle("Density Plots of Surface RAFs by Sample Timing") +
  density_theme
# ggsave("surface_timing_density.png", path = "figures", width = 12, height = 8)


# Surface Density Plots
data %>%
  # mutate(
  #   Location = as.character(Location),
  #   Location = ifelse(
  #     Location == "perimeter", "negative site", Location
  #   ) %>%
  #     factor(levels = c("negative site", "dump site"))
  # ) %>%
  ggplot(
    aes(RAF, forcats::fct_rev(Surface), fill = Location, linetype = Timing)
  ) +
  geom_density_ridges(alpha = 0.5, linewidth = 0.75) +
  scale_fill_manual(values = c("cyan", "darkorange", "maroon")) +
  ggtitle("Density Plots of Surface RAFs by Location & Timing") +
  density_theme
# ggsave(
#   "surface_location_timing_density.png", 
#   path = "figures", 
#   width = 12, 
#   height = 8
# )





# Plate Views -------------------------------------------------------------



for (file in list.files("raw")) {
  
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

