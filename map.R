library(tidyverse)
library(janitor)
library(ggpubr)
library(png)
library(scales)
library(ggnewscale)
library(ggforce)


main_theme <- theme(
  plot.title = element_text(size = 30, hjust = 0.5),
  axis.title = element_text(size = 20),
  axis.text = element_text(size = 16),
  strip.text = element_text(size = 20),
  legend.title = element_text(size = 20, hjust = 0.5),
  legend.text = element_text(size = 12),
  panel.grid = element_blank()
)


site_img <- readPNG("figures/dump_site.png")
img_dims <- dim(site_img)
dim_ratio <- img_dims[1] / img_dims[2]

df_ <- read.csv("data/bird-feeder/data.csv", check.names = FALSE) %>%
  clean_names(replace = c("Sample IDs" = "sample_ids", "TtT" = "ttt")) %>%
  filter(
    !(sample_ids %in% c("N", "P")),
    surface == "SS",
    !control,
    str_detect(feeder, "BF"),
    # time_since_start != 0
  ) %>%
  select(-c("dilutions", "tech", "collector", "surface")) %>%
  arrange(feeder) %>%
  mutate(
    across(c(mpr, ms, auc, raf, ttt), as.numeric),
    across(c(sample_ids:reader, feeder:side, site), as.factor),
    across(start_date:end_date, as.Date),
    across(c(time_in_env, time_since_start), ~ as.difftime(., units = "days")),
    lat = rescale(lat, c(0, 1), c(47.69795, 47.70075)),
    long = rescale(long, c(0, 1), c(-94.65457, -94.64921)),
    time_since_start = ifelse(time_since_start == 699, 760, time_since_start),
    time_since_start = factor(time_since_start, levels = sort(unique(time_since_start)), labels = paste(sort(unique(time_since_start)), "days"))
  )

pc <- prcomp(df_[c("mpr", "ms", "auc")])
summary(pc)
df_$pc1 <- pc$x[,1]
df_$pc1 <- rescale(df_$pc1, c(1, 0))

df_sum <- df_ %>%
  summarize(
    across(pc1, mean),
    .by = c(feeder, lat, long)
  ) %>%
  arrange(feeder) %>%
  mutate(
    label_x = long,
    label_y = lat + 0.04
  )

df_ %>%
  arrange(pc1) %>%
  ggplot(aes(long, lat, color=pc1)) +
  background_image(site_img) +
  geom_point(size=3, position=position_jitter(0.02, 0.02)) +
  scale_color_gradient(low="lightblue", high="darkred") +
  new_scale_color() +
  # geom_circle(aes(x0=long, y0=lat, r=0.04), data=df_sum, color="white") +
  geom_text(aes(label_x, label_y, label=feeder), data=df_sum, color = "white", size=6, fontface=2) +
  facet_wrap(vars(time_since_start)) +
  coord_fixed(dim_ratio) +
  main_theme +
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(5, "cm"),
    legend.title = element_text(vjust = 1)
  )

ggsave("map.png", path="figures/bird-feeder", width=12, height=13.5)
