library(dplyr)
library(ggplot2)
library(quicR)




data <- readr::read_csv("data/updated_metadata.csv", name_repair = "minimal")


ggplot(metadata, aes(`Sample IDs`, RAF, fill = Timing)) +
  geom_boxplot() +
  facet_grid(rows = ~Zone, scales = "free_x", space = "free_x") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "bottom"
  )


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

