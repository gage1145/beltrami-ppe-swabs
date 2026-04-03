library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
library(rnaturalearthdata)
library(quicR)



df_ <- read.csv("data/bird-feeder/data.csv", check.names = FALSE)

# Square Format
df_ %>%
  filter(!is.na(Feeder) & !is.na(Side)) %>%
  ggplot(aes(time_since_start, AUC, color = Feeder)) +
  geom_point() +
  facet_grid(vars(Side))

local_map <- ne_download("large", "states", category="raster") %>%
  filter(name == "Minnesota")

df_map <- df_ %>%
  filter(!is.na(Lat) & str_detect(Feeder, "ctrl", negate = TRUE))

local_map %>%
  ggplot() +
  geom_sf(fill="lightgrey", color="black") +
  geom_point(aes(Long, Lat, size=AUC), data = df_map)

ggsave("boxplot.png", path="figures/bird-feeder", width=16, height=8)



# Plate Views -------------------------------------------------------------



multi_plate <- function(file) {
  plate_view(
    get_quic(file)
  ) +
    ggtitle(file)
  ggsave(
    paste0(
      "figures/", 
      str_split_i(file, "/", 2) %>% str_remove(".xlsx"), 
      "_plate_view.png"
    ), 
    width = 12, 
    height = 8
  )
}

sapply(list.files("raw/bird-feeder", full.names = TRUE), multi_plate)


