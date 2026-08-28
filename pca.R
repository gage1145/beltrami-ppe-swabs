library(tidyverse)
library(janitor)
library(ggpubr)
library(modelr)
library(magrittr)
library(broom)
library(scales)

main_theme <- theme(
  plot.title = element_text(size = 30, hjust = 0.5),
  axis.title = element_text(size = 20),
  axis.text = element_text(size = 16),
  strip.text = element_text(size = 20),
  legend.title = element_text(size = 20, hjust = 0.5),
  legend.text = element_text(size = 16),
)


df_ppe <- read.csv("data/ppe/updated_metadata.csv", check.names = FALSE) %>%
  filter(surface == "boot") %>%
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

df_ %>%
  pivot_longer(c(mpr, ms, auc, mpr_log, ms_log, auc_log), names_to = "metric") %>%
  ggplot(aes(value)) +
  geom_histogram(position = "stack", bins=50) +
  facet_grid(cols=vars(metric), scales="free") 



metrics <- c("mpr", "mpr_log", "ms", "ms_log", "auc", "auc_log")
metric_combos <- as.data.frame(t(combn(metrics, 2)))
colnames(metric_combos) <- c("xval", "yval")
metric_combos %<>%
  filter(
    str_detect(xval, yval, negate=TRUE),
    str_detect(yval, xval, negate=TRUE)
  )

# Get correlation values for each pair of metrics.
cor_df <- pmap_dfr(
  metric_combos,
  function(xval, yval) {
    if (str_detect(xval, yval) | str_detect(yval, xval)) {
      return()
    }
    cor_val <- df_ %>%
      select(!!sym(xval), !!sym(yval)) %>%
      cor(use="pairwise.complete.obs")
    cor_val <- cor_val[2,1]
    tibble(
      xval = xval, 
      yval = yval, 
      cor_val = cor_val, 
      log_x = str_detect(xval, "log"), 
      log_y = str_detect(yval, "log"))
  } %>%
    mutate(across(c(xval, yval), ~ str_remove(., "_log")))
) %>%
  arrange(desc(cor_val)) %>%
  filter(log_x == log_y)

# df_cor_sum <- cor_df %>%
#   filter(
#     cor_val == max(cor_val), 
#     .by = c(xval, yval)
#   ) %>%
#   mutate(
#     xval = ifelse(log_x, paste0(xval, "_log"), xval),
#     yval = ifelse(log_y, paste0(yval, "_log"), yval)
  # )

metric_combos <- cor_df %>%
  arrange(log_x, xval, yval) %>%
  mutate(
    xval = ifelse(log_x, paste0(xval, "_log"), xval),
    yval = ifelse(log_y, paste0(yval, "_log"), yval)
  ) %>%
  select(xval, yval)

mc_no_log <- metric_combos %>%
  filter(str_detect(xval, "_log") == FALSE & str_detect(yval, "_log") == FALSE)

make_cor_plot <- function(xval, yval) {
  df_log <- df_ %>%
    mutate(across(c(!!sym(xval), !!sym(yval)), log))

  form <- as.formula(sprintf("%s ~ %s", yval, xval))

  mod <- lm(form, df_)
  mod_log <- lm(form, df_log)

  grid <- df_ %>%
    data_grid(!!sym(xval)) %>%
    add_predictions(mod, var = "untransformed")

  grid_log <- df_log %>%
    data_grid(!!sym(xval)) %>%
    add_predictions(mod_log, var = "log_transformed") %>%
    mutate(log_transformed = exp(log_transformed))

  grid <- grid %>%
    mutate(log_transformed = grid_log$log_transformed) %>%
    pivot_longer(c(untransformed, log_transformed), names_to = "type", values_to = "pred")

  rr <- glance(mod)$r.squared
  rr_log <- glance(mod_log)$r.squared

  max_y <- max(df_[[yval]], na.rm = TRUE)
  min_y <- min(df_[[yval]], na.rm = TRUE)
  max_x <- max(df_[[xval]], na.rm = TRUE)
  min_x <- min(df_[[xval]], na.rm = TRUE)

  rr_label_df <- tibble(
    x = c(min_x, min_x),
    y = c(0.95 * max_y, 0.65 * max_y),
    hjust = c(0, 0),
    label = c(paste("Untransformed R^2 =", round(rr, 3)), paste("Log-Transformed R^2 =", round(rr_log, 3))),
    type = c("untransformed", "log_transformed")
  )

  df_ %>%
    ggplot(aes(x = !!sym(xval), y = !!sym(yval))) +
    geom_hex(bins=50) +
    geom_line(aes(y=pred, color=type), data=grid, linewidth=2) +
    geom_label(
      aes(x = x, y = y, label = label, color=type, hjust=hjust), data=rr_label_df,
      size = 5, vjust = 0, show.legend = FALSE
    ) +
    scale_fill_gradient(low="cyan4", high="cyan2") +
    scale_color_manual(values=c("#d70500", "darkred")) +
    scale_x_log10(limits=c(min(df_[[xval]], na.rm=TRUE), max(df_[[xval]], na.rm=TRUE))) +
    scale_y_log10(limits=c(min(df_[[yval]], na.rm=TRUE), max(df_[[yval]], na.rm=TRUE))) +
    labs(
      fill = "Count"
    ) +
    main_theme +
    theme(
      legend.title = element_text(vjust = 1),
      legend.key.width = unit(7, "cm")
    ) 
}

cor_plots <- pmap(mc_no_log, make_cor_plot)
ggarrange(plotlist = cor_plots, align = "hv", ncol = 3, common.legend = TRUE, legend="none")
ggsave("corplot.png", path = "figures/pca", bg="white", width = 12, height = 6)



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
