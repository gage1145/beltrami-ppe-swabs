library(tidyverse)
library(janitor)
library(ggpubr)
library(modelr)
library(magrittr)
library(broom)
library(scales)
library(plotly)

main_theme <- theme(
  plot.title = element_text(size = 30, hjust = 0.5),
  axis.title = element_text(size = 20),
  axis.text = element_text(size = 16),
  strip.text = element_text(size = 20),
  legend.title = element_text(size = 20, hjust = 0.5),
  legend.text = element_text(size = 16),
)


df_ppe <- read.csv("data/ppe/updated_metadata.csv", check.names = FALSE) %>%
  filter(
    surface == "boot",
    location != "perimeter",
    timing != "post-decon"
  ) %>%
  select(sample_id, mpr, ms, auc, location) %>%
  mutate(
    type = "ppe",
    sample_id = as.character(sample_id)
  ) %>%
  na.omit()

df_bf <- read.csv("data/bird-feeder/data.csv", check.names = FALSE) %>%
  filter(surface == "SS") %>%
  clean_names(replace = c("Sample IDs" = "sample_id", "TtT" = "ttt")) %>% 
  select(sample_id, mpr, ms, auc, site) %>%
  mutate(
    type = "bf",
    sample_id = as.character(sample_id)
  ) %>%
  rename(location = site) %>%
  mutate(location = case_when(
    location == "perimeter" ~ "outer",
    location == "dump site" ~ "inner",
    TRUE ~ location
  )) %>%
  na.omit()

df_ <- bind_rows(df_ppe, df_bf)

df_log <- mutate(df_, across(mpr:auc, log))

mprauc_plot <-  df_log %>%
    ggplot(aes(mpr, auc)) +
    geom_point(alpha=0.1) +
    # stat_smooth(method="loess", color="blue", linetype="dashed") +
    stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")), label.y.npc=1) +
    labs(
      x = "log(MPR)", y = "log(AUC)"
    ) +
    main_theme
mprms_plot <- df_log %>%
    ggplot(aes(mpr, ms)) +
    geom_point(alpha=0.1) +
    # stat_smooth(method="loess", color="blue", linetype="dashed") +
    stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")), label.y.npc=1) +
    labs(
      x = "log(MPR)", y = "log(MS)"
    ) +
    main_theme
msauc_plot <- df_log %>%
    ggplot(aes(auc, ms)) +
    geom_point(alpha=0.2) +
    # stat_smooth(method="loess", color="blue", linetype="dashed") +
    stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")), label.y.npc=1) +
    labs(
      x = "log(AUC)", y = "log(MS)"
    ) +
    main_theme

ggarrange(mprauc_plot, mprms_plot, msauc_plot, ncol=3)
ggsave("cor_log.jpg", path="figures/pca", width=12, height=6, dpi=600)


# 3D Plot ----------------------------------------------------------------


fig <- df_log %>%
  mutate(location = ifelse(location %in% c("dump site", "inner"), "Core", "Non-Core")) %>%
  plot_ly() %>%
  add_trace(
    x = ~ms,
    y = ~auc,
    z = ~mpr,
    color = ~location,
    colors = c("darkorange", "darkcyan"),
    opacity = 0.5,
    type = "scatter3d",
    mode = "markers",
    size = 1
  )
fig


# Principal Component Analysis -------------------------------------------


pca_res <- df_log %>%
  select(mpr, ms, auc) %>%
  prcomp(scale. = TRUE) %>%
  summary()

df_pca <- df_ %>%
  bind_cols(pca_res$x)

bf_pca_plot <- df_pca %>%
  filter(type == "bf") %>%
  mutate(
    # across(c(PC1, PC2), rescale),
    location = factor(location, levels = c("control", "inner", "outer"))
) %>%
  arrange(location) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color=location), size=5, alpha=0.5) +
  scale_color_manual(values=c("darkcyan", "red", "orange2")) +
  main_theme +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom"
  ) 
bf_pca_plot

ppe_pca_plot <- df_pca %>%
  filter(type == "ppe") %>%
  mutate(
    # across(c(PC1, PC2), rescale),
    location = factor(location, levels = c("negative site", "dump site", "perimeter"))
) %>%
  arrange(desc(location)) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color=location), size=5, alpha=0.5) +
  scale_color_manual(values=c("darkcyan", "red", "orange2")) +
  main_theme +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom"
  ) 
ppe_pca_plot

ggarrange(bf_pca_plot, ppe_pca_plot, ncol = 2)
