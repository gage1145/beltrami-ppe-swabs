library(dplyr)
library(ggplot2)
library(quicR)
library(ggridges)
library(readr)




data <- read_csv("data/updated_metadata.csv", name_repair = "minimal") %>%
  mutate_at(c("Sample IDs", "Technician"), as.factor) %>%
  mutate(
    Timing = factor(Timing, levels = c("pre sampling", "post sampling")),
    Surface = factor(
      Surface,
      levels = c(
        "blank sample", "tyvek leg", "plastic boot", "boot side", "boot bottom"
      )
    ),
    Location = factor(
      Location, 
      levels = c("negative site", "perimeter", "dump site")
    )
  )

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
  filter(Timing != "post boot decon") %>%
  ggplot(aes(Surface, RAF, color = Location)) +
  geom_boxplot(fill = NA, size = 1) +
  scale_color_manual(values = c("darkcyan", "darkorange", "maroon")) +
  ggtitle("Surface RAFs by Location") +
  facet_grid(
    cols = vars(Timing), 
    rows = vars(Technician), 
    scales = "free_x", 
    space = "free_x"
  ) +
  box_theme
# ggsave("RAF.png", path = "figures", width = 16, height = 6)


data %>%
  filter(Timing != "post boot decon") %>%
  ggplot(aes(Surface, RAF, color = Location)) +
  geom_boxplot(fill = NA, size = 1) +
  scale_color_manual(values = c("darkcyan", "darkorange", "maroon")) +
  facet_grid(cols = vars(Timing), scales = "free_x", space = "free_x") +
  box_theme


data %>%
  filter(Timing != "post boot decon") %>%
  ggplot(aes(Location, RAF, fill = Timing)) +
  geom_violin(size = 1) +
  # scale_color_manual(values = c("darkcyan", "darkorange", "maroon")) +
  facet_grid(cols = vars(Surface), scales = "free_x", space = "free_x") +
  box_theme





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

