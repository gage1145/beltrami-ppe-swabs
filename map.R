library(tidyverse)
library(magrittr)
library(janitor)
library(ggpubr)
library(png)
library(scales)
library(ggnewscale)
library(ggforce)
library(skimr)


set.seed(987654321)

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
  # left_join(coord_ref) %>%
  mutate(
    across(c(mpr, ms, auc, raf, ttt), as.numeric),
    across(c(sample_ids:reader, feeder:side, site), as.factor),
    across(start_date:end_date, as.Date),
    across(c(time_in_env, time_since_start), ~ as.difftime(., units = "days")),
    lat = rescale(lat, c(0, dim_ratio), c(47.69795, 47.70075)),
    long = rescale(long, c(0, 1), c(-94.65457, -94.64921)),
    time_since_start = ifelse(time_since_start == 699, 760, time_since_start),
    time_since_start = factor(time_since_start, levels = sort(unique(time_since_start)), labels = paste(sort(unique(time_since_start)), "days"))
  )

pc <- prcomp(df_[c("mpr", "ms", "auc")])
# summary(pc)

df_ %<>% mutate(pc1 = rescale(pc$x[,1], c(1, 0)))

df_sum <- df_ %>%
  summarize(
    across(pc1, mean),
    .by = c(feeder, lat, long)
  ) %>%
  arrange(feeder) %>%
  mutate(
    label_x = long,
    label_y = lat
  )

jitter_amt <- 0.03

# skim(df_)

# df_ %>%
#   # mutate(pc1=log(pc1)) %>%
#   ggplot(aes(pc1)) +
  # geom_histogram()

df_ %>%
  arrange(pc1) %>%
  ggplot(aes(long, lat, color=pc1)) +
  background_image(site_img) +
  geom_circle(aes(x0=long, y0=lat, r=0.05), color="white", linewidth=0.3) +
  geom_point(size=6, alpha=0.7, position=position_jitter(jitter_amt, jitter_amt, 10548237)) +
  geom_text_repel(
    aes(label_x, label_y, label=feeder), 
    data=df_sum, color = "white", size=6, fontface=2,
    nudge_y=0.08, direction="x", min.segment.length = 10
  ) +
  facet_wrap(vars(time_since_start)) +
  scale_color_gradient2(low="lightblue", mid="darkorange", high="darkred", midpoint=0.6) +
  scale_x_continuous(limits=c(0, 1)) +
  scale_y_continuous(limits=c(0, dim_ratio)) +
  coord_fixed() +
  main_theme +
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(5, "cm"),
    legend.title = element_text(vjust = 1)
  )

ggsave("map.png", path="figures/bird-feeder", width=12, height=11)
