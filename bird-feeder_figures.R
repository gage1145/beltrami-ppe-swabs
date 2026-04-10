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
  clean_names(replace=c("Sample IDs" = "sample_ids", "TtT" = "ttt")) %>%
  filter(!(sample_ids %in% c("N", "P"))) %>%
  mutate(
    across(c(sample_ids: tech, feeder: collector, site), as.factor),
    across(start_date: end_date, as.Date),
    across(c(time_in_env, time_since_start), ~ as.difftime(., units = "days"))
  )

skim(df_)


# Principal Component Analysis -------------------------------------------


pca_res <- df_ %>%
#   mutate(time_since_start = as.numeric(time_since_start)) %>%
  prcomp(select(df_, mpr, ms, auc), scale. = TRUE)
# summary(pca_res)

pca_df <- pca_res$x %>%
  as.data.frame() %>%
  mutate(
    time_since_start = df_$time_since_start,
    site = df_$site,
    control = df_$control,
    side = df_$side,
    reextracted = df_$reextracted
  )

pca_long <- pca_df %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(starts_with("PC"), names_to = "PC", values_to = "value")

pca_comp_df <- inner_join(
  pca_long %>% rename(x_pc = PC, x_val = value),
  pca_long %>% rename(y_pc = PC, y_val = value),
  relationship = "many-to-many"
) %>%
  filter(x_pc != y_pc)

pca_comp_df %>%
  ggplot(aes(x_val, y_val, color = site)) +
  geom_point(size=3, alpha=0.5) +
#   stat_ellipse(geom="polygon", level=0.95, alpha=0.2) +
  facet_grid(rows = vars(y_pc), cols = vars(x_pc)) +
  scale_color_manual(values=c("cyan", "darkorange", "darkgreen")) +
  coord_fixed() +
  labs(
    x = NULL, y = NULL,
    title = "Principal Component Pairs"
  ) +
  main_theme
ggsave("pca.png", path="figures", width=12, height=12)


# Summary ----------------------------------------------------------------


df_sum <- df_ %>%
  summarize(
    across(c("mpr", "ms", "auc"), mean),
    n = n(),
    .by = c(feeder, site, side, surface, control, reextracted, time_in_env, time_since_start)
  ) %>%
  filter(!control, time_in_env > 0, surface != "Plastic") %>%
  pivot_longer(c("mpr", "ms", "auc"), names_to = "metric")

df_neg <- df_ %>%
  filter(control) %>%
  pivot_longer(c("mpr", "ms", "auc"), names_to = "metric")

df_neg_sum <- df_neg %>%
  summarize(neg_mean = mean(value), neg_sd = sd(value), .by = metric)

df_site <- df_ %>%
  mutate(time_since_start = as.numeric(time_since_start)) %>%
  filter(!control, time_in_env > 0, surface != "Plastic") %>%
  pivot_longer(c("mpr", "ms", "auc"), names_to = "metric")

groupings <- expand.grid(
    m = unique(df_site$metric),
    st = unique(df_site$site),
    sd = unique(df_site$side),
    tss = unique(df_site$time_since_start)
)


get_comps <- function(m, st, sd, tss) {
  temp <- filter(df_sum, metric == m, site == st, side == sd, time_since_start == tss)
  if (nrow(temp) < 1) return(data.frame())
  temp_neg <- filter(df_neg, metric == m)
  value_map <- list(
    mpr = df_neg_sum$neg_mean[df_neg_sum$metric == "mpr"], 
    ms = df_neg_sum$neg_mean[df_neg_sum$metric == "ms"], 
    auc = df_neg_sum$neg_mean[df_neg_sum$metric == "auc"]
  )
#   scalar <- (st == "dump site") * 0.12 + 1

  tidy(wilcox.test(temp$value, temp_neg$value, alternative = "greater", exact = FALSE)) %>%
    mutate(
        metric = m, 
        site = st, 
        side = sd,
        time_since_start = tss,
        value=value_map[[m]], 
        significant = p.value <= 0.05,
        p.value=paste("p =", signif(p.value, 3))
    )
}

comps <- pmap_dfr(groupings, get_comps)

df_sum %>%
  ggplot(aes(x=time_since_start, y=value, color=site, label=feeder, group=feeder)) +
  geom_text() +
#   geom_point() +
  stat_smooth(method="lm", linewidth=1.5) +
  geom_hline(aes(yintercept=neg_mean), data=df_neg_sum, linetype="dashed", linewidth=1.2) +
  geom_label_repel(aes(time_since_start, value, label=p.value, color=site), size = 6, data=comps, inherit.aes = FALSE) +
  facet_grid(vars(metric), vars(side), scales="free") +
  theme(
    axis.text = element_text(size=16),
    axis.title = element_text(size=20),
    strip.text = element_text(size=20),
    legend.title = element_blank(),
    legend.text = element_text(size=16),
  )

metrics <- c("mpr", "ms", "auc")

get_mods <- function(x) {
  form <- as.formula(paste(x, " ~ site * time_since_start * side"))
  df_ <- pivot_wider(df_site, names_from = "metric", values_from = "value")
  mod <- glm(form, data = df_)
  tidy(mod) %>%
    mutate(metric = x)
}

site_mods <- map_dfr(metrics, get_mods)
# Looks like as time goes on, MPR and AUC detect a significant negative trend in activity.
# There are no significant differences between the two sites as a whole.
# Slopes do not diverge between sites. They have parallel slopes.

# Compare against negative site.
