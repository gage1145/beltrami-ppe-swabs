library(tidyverse)


df_ <- read.csv("data/bird-feeder/data.csv", check.names = FALSE) %>%
  clean_names(replace=c("Sample IDs" = "sample_ids", "TtT" = "ttt"))

df_sum <- df_ %>%
  summarize(
    across(c("mpr", "ms", "auc", "raf", "ttt"), c(mean = mean, sd = sd)),
    n = n(),
    .by = c(sample_ids, control, lat, long, rxn_date, reextracted, feeder, side, surface)
  )

bkg_mpr <- mean(df_$mpr[which(df_$control)])
bkg_ms <- mean(df_$ms[which(df_$control)])
bkg_auc <- mean(df_$auc[which(df_$control)])
bkg_raf <- mean(df_$raf[which(df_$control)])

df_norm <- df_ %>%
  mutate(
    mpr = mpr - bkg_mpr,
    ms = ms - bkg_ms,
    auc = auc - bkg_auc,
    raf = raf - bkg_raf
  )

df_norm %>%
  pivot_longer(c(mpr, ms, auc, raf), names_to = "metric") %>%
  ggplot(aes(feeder, value, fill=side)) +
  geom_boxplot() +
  facet_wrap(vars(metric), scales="free")

df_sum %>%
  filter(!control, lat > 47.695) %>%
  group_by(feeder) %>%
  ggplot(aes(long, lat, label=feeder, size=auc_mean)) +
  # geom_point() + 
  geom_text() +
  theme(
    axis.title.x = element_blank(),
    plot.title = element_text(size=20, face="bold", hjust=0.5),
    strip.text = element_text(size=16),
    axis.title = element_text(size=16),
    axis.text = element_text(size=12),
    legend.title = element_text(size=16, hjust=0.5),
    legend.text = element_text(size=12)
    
  )
ggsave("feeder_box.png", path="figures/bird-feeder", width=12, height=8)
 