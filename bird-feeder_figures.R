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
  mutate(across(everything(), log)) %>%
  prcomp(scale. = TRUE)
summary(pca_res)

pca_df <- df_ %>%
  bind_cols(pca_res$x) %>%
  mutate(
    # across(matches("PC"), ~ -. + max(.) + 1),
    year = as.factor(year(end_date)),
  )

pca_long <- pca_df %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(c("mpr", "ms", "auc", "PC1", "PC2", "PC3"), names_to = "metric", values_to = "value")

# pc1_scalar <- max(pca_df$PC1) + 1

# pca_comp_df <- inner_join(
#   pca_long %>% rename(x_pc = metric, x_val = value),
#   pca_long %>% rename(y_pc = metric, y_val = value),
#   relationship = "many-to-many"
# ) %>%
#   filter(x_pc != y_pc) %>%
#   arrange(x_pc, desc(y_pc))

pca_combos <- as.data.frame(t(combn(paste0("PC", seq(length(pca_res$scale))), 2)))
colnames(pca_combos) <- c("x", "y")

make_pca_plot <- function(x, y, alpha=0.5) {
#   pcs <- paste0("PC", c(x, y))
  cls <- pca_df %>%
    select(!!sym(x), !!sym(y)) %>%
    kmeans(centers=2, trace=TRUE) %T>%
    summary()
  pca_df %>%
    mutate(
      cluster = as.character(cls$cluster),
      control = ifelse(control, "Negative Site", "Dump Site/Perimeter")
    ) %>%
    arrange(control) %>%
    ggplot(aes(x = !!sym(x), y = !!sym(y))) +
    geom_point(aes(color=control), alpha = alpha, size=3) +
    stat_ellipse(aes(fill=cluster), level=0.99, geom="polygon", alpha=0, type="t", color="black", show.legend = FALSE) +
    scale_color_manual(values=c("darkorange", "darkcyan")) +
    scale_fill_manual(values=c("darkorange", "darkcyan")) +
    # facet_grid(~year) +
    # scale_y_continuous(limits=c(-2, 6)) +
    # scale_x_continuous(limits=c(-2, 6)) +
    coord_fixed() +
    labs(
      x = toupper(x),
      y = toupper(y)
    ) +
    main_theme
}

pca_plots <- pmap(pca_combos, make_pca_plot)
ggarrange(plotlist = pca_plots, align = "hv", ncol=3, common.legend = TRUE, legend="bottom") 

ggsave("pca.png", path = "figures/bird-feeder", bg="white", width = 16, height = 4)

# pca_comp_df %>%
#   filter(str_detect(x_pc, "PC") & str_detect(y_pc, "PC")) %>%
#   arrange(control) %>%
#   mutate(control = ifelse(control, "Negative Site", "Dump Site/Perimeter")) %>%
#   ggplot(aes(x_val, y_val, color = control)) +
#   geom_point(size = 3, alpha = 0.75) +
#   # stat_ellipse(geom="polygon", level=0.95, alpha=0.2) +
#   facet_grid(rows = vars(y_pc), cols = vars(x_pc), scales = "free", switch = ("y")) +
#   scale_y_continuous(position = "right") +
#   scale_color_manual(values = c("red", "darkcyan")) +
#   guides(color = guide_legend(override.aes = list(alpha = 1))) +
#   # stat_ellipse() +
#   # coord_fixed() +
#   labs(
#     x = NULL, y = NULL,
#   ) +
#   main_theme +
#   theme(
#     legend.title = element_blank(),
#     legend.position = "top"
#   )


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
    site = factor(site, levels = c("dump site", "perimeter"), labels = c("Dump Site", "Perimeter")),
    year = as.integer(as.character(year))
  ) %>%
  ggplot(aes(x = year, y = -value - neg_value, color=site, fill=site)) +
  geom_point(size=3, position=position_jitterdodge(0.04, 0, 0.1, seed=56543)) +
  # stat_smooth(method = "lm", linewidth = 1.5, position=position_dodge(0.1)) +
  geom_smooth(linewidth=1.5, position=position_dodge(0.1)) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
    size = 6, label.x.npc=0.7, label.y.npc=1, show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2, color="blue") +
  scale_x_continuous(breaks = seq(2023, 2025)) +
  # scale_fill_manual(values=c("darkred", "darkorange")) +
  # scale_color_manual(values=c("darkred", "darkorange")) +
  labs(
    x = "Year of Sampling",
    y = "PC1 - mean Neg PC1"
  ) +
  main_theme +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = "top"
  )
ggsave("pca_vs_tss.png", path = "figures/bird-feeder", width = 12, height = 10)


# Plots by plate ---------------------------------------------------------


pca_plates <- pca_df %>%
  filter(
    str_detect(feeder, "BF"),
    time_in_env > 0
  ) %>%
  mutate(
    # proximity = ifelse(feeder %in% paste0("BF", seq(1,3)), "dump site", "perimeter")
  ) 

plate_sum <- pca_plates %>%
  summarize(
    across(PC1, c(mean=~-mean(.), stdev=sd)),
    .by = c(feeder, year, site)
  )

prox_sum <- plate_sum %>%
  summarize(
    stack_top = sum(PC1_mean),
    .by = c(year, site)
  ) %>%
  left_join(
    pca_plates %>%
      summarize(PC1_stdev = sd(PC1), .by = c(year, site)),
    by = c("year", "site")
  )

plate_sum %>%
  ggplot(aes(site, PC1_mean, fill=feeder, group=site)) +
  geom_col() +
  geom_point(
    aes(site, stack_top, group=site),
    data=prox_sum, inherit.aes = FALSE, position = position_dodge(0.9), shape="square", size=3
  ) +
  geom_errorbar(
    aes(site, stack_top, ymin=stack_top - PC1_stdev, ymax=stack_top + PC1_stdev, group=site),
    data=prox_sum, inherit.aes = FALSE, position = position_dodge(0.9), width=0.2
  ) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=1) +
  stat_compare_means(
    aes(site, PC1),
    data = pca_plates,
    method = "wilcox.test",
    label.y.npc = 1, label.x.npc=0.5, hjust=0.5,
    label = "p.format", size = 6, show.legend = FALSE
  ) +
  facet_grid(cols=vars(year)) +
  scale_y_continuous(limits=c(NA, 6.4)) +
  labs(
    y="PC1"
  ) +
  scale_fill_manual(values=c("cyan4", "cyan3", "cyan2", "orange4", "orange3", "orange2")) +
  main_theme +
  theme(
    axis.title.x = element_blank(),
    legend.title = element_blank()
  )
ggsave("plate_cols.png", path="figures/bird-feeder", width=12, height=10)
