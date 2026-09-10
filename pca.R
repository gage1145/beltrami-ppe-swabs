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
    stat_smooth(method="lm", color="blue", linetype="dashed") +
    stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")), label.y=6.8) +
    labs(
      fill = "Count"
    ) +
    main_theme
mprms_plot <- df_log %>%
    ggplot(aes(mpr, ms)) +
    geom_point(alpha=0.1) +
    stat_smooth(method="lm", color="blue", linetype="dashed") +
    stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")), label.y=3.4) +
    labs(
      fill = "Count"
    ) +
    main_theme
msauc_plot <- df_log %>%
    ggplot(aes(ms, auc)) +
    geom_point(alpha=0.1) +
    stat_smooth(method="lm", color="blue", linetype="dashed") +
    stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")), label.y=6.8) +
    labs(
      fill = "Count"
    ) +
    main_theme

ggarrange(mprauc_plot, mprms_plot, msauc_plot, ncol=3)
ggsave("cor_log.png", path="figures/pca", width=12, height=6)


df_ <- df_ %>%
  mutate(across(mpr:auc, list(log=log)))

write.csv(df_, "data/pca/combined.csv", row.names=FALSE)


# Linear regression
met_mod <- lm(mpr ~ ms + auc + transformation, data = df_long)

rr <- modelr::rsquare(met_mod, df_long %>% filter(transformation == "none"))
rr_log <- modelr::rsquare(met_mod, df_long %>% filter(transformation == "log"))

df_long %<>%
  mutate(pred = predict(met_mod, type = "response"))

df_long_sum <- df_long %>%
  summarize(
    across(c(mpr, ms, auc), list(min = min, max = max)),
    .by = transformation
  )

df_pred <- tibble(
  ms = seq(df_long_sum$ms_min[1], df_long_sum$ms_max[1], length.out = 10),
  auc = seq(df_long_sum$auc_min[1], df_long_sum$auc_max[1], length.out = 10),
  transformation = "none"
) %>%
  add_predictions(met_mod)

df_pred_log <- tibble(
  ms = seq(df_long_sum$ms_min[2], df_long_sum$ms_max[2], length.out = 10),
  auc = seq(df_long_sum$auc_min[2], df_long_sum$auc_max[2], length.out = 10),
  transformation = "log"
) %>%
  add_predictions(met_mod)


df_long %>%
  filter(transformation == "none") %>%
  arrange(pred) %>%
  plot_ly() %>%
  add_trace(
    x = ~ms,
    y = ~auc,
    z = ~mpr,
    type = "scatter3d",
    size = 1
  ) %>%
  add_trace(
    data = df_pred,
    x = ~ms,
    y = ~auc,
    z = ~pred,
    type = "scatter3d",
    mode = "lines"
  )

df_long %>%
  filter(transformation == "log") %>%
  arrange(pred) %>%
  plot_ly() %>%
  add_trace(
    x = ~ms,
    y = ~auc,
    z = ~mpr,
    type = "scatter3d",
    size = 1
  ) %>%
  add_trace(
    data = df_pred_log,
    x = ~ms,
    y = ~auc,
    z = ~pred,
    type = "scatter3d",
    mode = "lines"
  )


metrics <- c("mpr", "mpr_log", "ms", "ms_log", "auc", "auc_log")
metric_combos <- as.data.frame(t(combn(metrics, 2)))
colnames(metric_combos) <- c("xval", "yval")
metric_combos %<>%
  filter(
    str_detect(xval, yval, negate=TRUE),
    str_detect(yval, xval, negate=TRUE)
  )

# Get correlation values for each pair of metrics.
cor_mat <- df_ %>%
  select(mpr, ms, auc) %>%
  cor()
write.csv(cor_mat, "data/pca/cor_mat.csv")




# Principal Component Analysis -------------------------------------------


pca_res <- df_ %>%
  select(mpr, ms, auc) %>%
  prcomp(scale. = TRUE) %>%
  summary()

pca_res_log <- df_ %>%
  select(mpr, ms, auc) %>%
  mutate(across(everything(), log)) %>%
  prcomp(scale. = TRUE) %>%
  summary()

sprintf("Untransformed PC1 = %.2f", pca_res$importance[2,1])
sprintf("Log-transformed PC1 = %.2f", pca_res_log$importance[2,1])

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
