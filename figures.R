library(tidyverse)
library(quicR)



meta <- read.csv("data/meta.csv", check.names = FALSE) %>%
  mutate_at("Sample IDs", as.character) %>%
  na.omit()

df_ <- read.csv("data/data.csv", check.names = FALSE) %>%
  mutate_at("Sample IDs", ~str_remove(., " ")) %>%
  filter(!(`Sample IDs` %in% c("N", "P"))) %>%
  left_join(meta, "Sample IDs")

mean_neg_raf <- {df_ %>%
  filter(is.na(ID)) %>%
  summarize(mean_raf = mean(RAF))}[[1,1]]

dodge_amount = 0.04

# Square Format
df_ %>%
  separate_wider_delim(
    ID, 
    " ", 
    names = c("Feeder", "Plate"), 
    too_few = "align_start", 
    too_many = "merge"
  ) %>%
  mutate(RAF = RAF - mean_neg_raf) %>%
  mutate(RAF = ifelse(RAF <= 0, 0, RAF)) %>%
  filter(Date != unique(Date)[2]) %>%
  filter(`Period in environment` != 'Aug 2021 - July 2024') %>%
  filter(Plate != 'Plastic') %>%
  group_by(Feeder, Year, Locale) %>%
  summarize(
    mean_RAF = mean(RAF),
    sd_RAF = sd(RAF),
    ymin=ifelse(mean_RAF-sd_RAF < 0, 0, mean_RAF-sd_RAF), 
    ymax=mean_RAF+sd_RAF, 
  ) %>%
  ggplot(aes(
    Year, 
    mean_RAF,
    ymin=ymin,
    ymax=ymax,
    color=Feeder,
    fill=Feeder,
    linetype=Locale
  )) +
  geom_line(
    linewidth=1.5, 
    position = position_dodge(width=dodge_amount)
  ) +
  geom_point(
    size=4.5, 
    position=position_dodge(width=dodge_amount)
  ) +
  geom_errorbar(
    linewidth=0.6, 
    linetype='solid', 
    position = position_dodge(width=dodge_amount)
  ) +
  scale_x_continuous(minor_breaks = seq(0, 3, 1/12)) +
  labs(
    y = "Background Corrected RAF",
    x = "Years"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size=20),
    axis.text = element_text(size=16),
    legend.title = element_text(size=20),
    legend.text = element_text(size=16)
  )
ggsave("raf_over_time.png", path="figures", width=12, height=12)



# Plate Views -------------------------------------------------------------



multi_plate <- function(file) {
  plate_view(
    get_real(file)[[1]], 
    get_sample_locations(file)
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

sapply(list.files("raw", full.names = TRUE), multi_plate)


