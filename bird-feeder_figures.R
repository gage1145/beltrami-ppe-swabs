library(tidyverse)
library(janitor)
library(skimr)
library(broom)
library(ggpubr)
library(ggrepel)


main_theme <- theme(
  plot.title = element_text(size = 30, hjust = 0.5),
  axis.title = element_text(size = 20),
  axis.text = element_text(size = 16),
  strip.text = element_text(size = 20),
  legend.title = element_text(size = 20, hjust = 0.5),
  legend.text = element_text(size = 16),
  panel.grid = element_blank()
)


# Load data --------------------------------------------------------------


df_ <- read.csv("data/bird-feeder/data.csv", check.names = FALSE) %>%
  clean_names(replace = c("Sample IDs" = "sample_ids", "TtT" = "ttt")) %>%
  filter(
    !(sample_ids %in% c("N", "P")),
    surface == "SS"
  ) %>%
  select(-c("dilutions", "tech", "collector", "surface")) %>%
  mutate(
    across(c(mpr, ms, auc, raf, ttt), as.numeric),
    across(c(sample_ids:reader, feeder:side, site), as.factor),
    across(start_date:end_date, as.Date),
    across(c(time_in_env, time_since_start), ~ as.difftime(., units = "days"))
  )

skim(df_)


# Principal Component Analysis -------------------------------------------


pca_res <- df_ %>%
  select(mpr, ms, auc) %>%
  prcomp(scale. = TRUE)
summary(pca_res)

pca_df <- bind_cols(df_, pca_res$x) %>%
  mutate(across(matches("PC"), ~ -. + max(.) + 1))

pca_long <- pca_df %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(c("mpr", "ms", "auc", "PC1", "PC2", "PC3"), names_to = "metric", values_to = "value")

# pc1_scalar <- max(pca_df$PC1) + 1

pca_comp_df <- inner_join(
  pca_long %>% rename(x_pc = metric, x_val = value),
  pca_long %>% rename(y_pc = metric, y_val = value),
  relationship = "many-to-many"
) %>%
  filter(x_pc != y_pc) %>%
  arrange(x_pc, desc(y_pc))

# metric_combos <- combn(unique(pca_comp_df$x_pc), 2) %>%
#   t() %>%
#   as.data.frame() %>%
#   rename(x_pc = V1, y_pc = V2)

# pca_comp_df <- right_join(pca_comp_df, metric_combos)

# pca_comp_df <- pca_comp_df %>%
#   mutate(
#     x_val = ifelse(x_pc %in% metric_combos$x_pc & y_pc %in% metric_combos$y_pc, x_val, NA),
#     y_val = ifelse(x_pc %in% metric_combos$x_pc & y_pc %in% metric_combos$y_pc, y_val, NA)
#   )

pca_comp_df %>%
  filter(str_detect(x_pc, "PC") & str_detect(y_pc, "PC")) %>%
  mutate(control = ifelse(control, "Negative Site", "Dump Site/Perimeter")) %>%
  ggplot(aes(x_val, y_val, color = site)) +
  geom_point(size = 3, alpha = 0.25) +
  # stat_ellipse(geom="polygon", level=0.95, alpha=0.2) +
  facet_grid(rows = vars(y_pc), cols = vars(x_pc), scales = "free", switch = ("y")) +
  scale_y_continuous(position = "right") +
  scale_color_manual(values = c("darkcyan", "red", "orange")) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  stat_ellipse() +
  # coord_fixed() +
  labs(
    x = NULL, y = NULL,
  ) +
  main_theme +
  theme(
    legend.title = element_blank(),
    legend.position = "top"
  )
ggsave("pca.png", path = "figures/bird-feeder", width = 12, height = 12.5)


# Summary ----------------------------------------------------------------


df_sum <- pca_long %>%
  summarize(
    across(c(value), mean),
    n = n(),
    .by = c(metric, feeder, site, side, control, reextracted, time_in_env, time_since_start)
  ) %>%
  filter(!control, time_in_env > 0)

df_neg <- pca_long %>%
  filter(control)

df_neg_sum <- df_neg %>%
  summarize(neg_mean = mean(value), neg_sd = sd(value), .by = metric)

df_site <- pca_long %>%
  filter(!control, time_in_env > 0)

groupings <- expand.grid(
  m = unique(df_site$metric),
  st = unique(df_site$site),
  sd = unique(df_site$side),
  tss = unique(df_site$time_since_start)
)


get_comps <- function(m, st, sd, tss) {
  temp <- filter(df_sum, metric == m, site == st, side == sd, time_since_start == tss)
  if (nrow(temp) < 1) {
    return(data.frame())
  }
  temp_neg <- filter(df_neg, metric == m)
  value_map <- list(
    mpr = df_neg_sum$neg_mean[df_neg_sum$metric == "mpr"],
    ms  = df_neg_sum$neg_mean[df_neg_sum$metric == "ms"],
    auc = df_neg_sum$neg_mean[df_neg_sum$metric == "auc"],
    PC1 = df_neg_sum$neg_mean[df_neg_sum$metric == "PC1"],
    PC2 = df_neg_sum$neg_mean[df_neg_sum$metric == "PC2"],
    PC3 = df_neg_sum$neg_mean[df_neg_sum$metric == "PC3"]
  )

  tidy(wilcox.test(temp$value, temp_neg$value, alternative = "greater", exact = FALSE)) %>%
    mutate(
      metric = m,
      site = st,
      side = sd,
      time_since_start = tss,
      value = value_map[[m]],
      significant = p.value <= 0.05,
      p.value = paste("p =", signif(p.value, 3))
    )
}

comps <- pmap_dfr(groupings, get_comps)
fig_metric <- "PC1"
neg_value <- df_neg_sum$neg_mean[df_neg_sum$metric == fig_metric]

filter_for_fig <- function(df) {
  df %>%
    filter(
      time_since_start > 700,
      !control,
      metric == fig_metric,
      !reextracted
    )
}


# PC1 vs TSS Plot --------------------------------------------------------


mods <- pca_long %>%
  filter_for_fig() %>%
  nest(.by = c(side, site)) %>%
  mutate(
    model = map(data, ~ lm(value ~ time_since_start, data = .x)),
    # tidied = map(model, broom::tidy),
    glanced = map(model, broom::glance)
  ) %>%
  unnest(c(glanced))
write.csv(select(mods, -c(data, model)), "data/bird-feeder/tss_results.csv")

pca_long %>%
  filter_for_fig() %>%
  mutate(
    side = factor(side, levels = c("Top", "Bottom")),
    site = factor(site, levels = c("dump site", "perimeter"), labels = c("Dump Site", "Perimeter"))
  ) %>%
  ggplot(aes(x = time_since_start, y = value - neg_value)) +
  geom_point() +
  stat_smooth(method = "lm", fill = "darkcyan", color = "darkcyan", linewidth = 1.5) +
  stat_cor(size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2) +
  facet_grid(vars(side), vars(site)) +
  scale_x_continuous(breaks = seq(800, 1500, 200)) +
  labs(
    x = "Days Since Start",
    y = "PC1 - mean Neg PC1"
  ) +
  main_theme +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
  )
ggsave("pca_vs_tss.png", path = "figures/bird-feeder", width = 12, height = 10)
