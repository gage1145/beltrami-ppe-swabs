library(tidyverse)
library(quicR)

main_theme <- theme(
  plot.title = element_text(size = 30, hjust = 0.5),
  axis.title = element_text(size = 20),
  axis.text = element_text(size = 16),
  strip.text = element_text(size = 20),
  legend.title = element_text(size = 20, hjust = 0.5),
  legend.text = element_text(size = 16),
  panel.grid = element_blank()
)

df_raw <- get_quic("raw/surface_recovery/20220929_gage_salivary_swabs.xlsx")
df_calc <- calculate_metrics(df_raw, threshold=3) %>%
  mutate(Dilutions = -log10(as.numeric(Dilutions))) %>%
  separate(`Sample IDs`, c("Sample IDs", "Surface"), " ") %>%
  filter(
    Surface %in% c("Ctrl", "Glass", "Steel"),
    `Sample IDs` %in% c("P"),
    Dilutions > -6 & Dilutions < -1
  )

df_ctrl <- filter(df_calc, Surface == "Ctrl")
df_var <- filter(df_calc, Surface != "Ctrl")
df_best <- df_var %>%
  summarize(
    across(RAF, c(mean = mean, sd = sd)),
    .by = c("Dilutions", "Surface")
  ) %>%
  filter(
    RAF_mean == max(RAF_mean), .by = "Surface"
  )

form <- RAF ~ Dilutions
model_data <- df_calc %>%
  filter(Surface == "Ctrl")
mod <- lm(form, model_data)
coefs <- coef(mod)
m <- coefs[2]
b <- coefs[1]

preds <- df_best %>%
  mutate(
    pred_raf = predict.lm(mod, list(Dilutions=Dilutions)),
    pred_dil = (RAF_mean - b) / m,
    recovery = RAF_mean / pred_raf,
    label = paste("Recovery =", round(recovery * 100, 1), "%"),
    label_x = c(-3, -3),
    label_y = c(0.0985, 0.095)
  )

df_ctrl %>%
  ggplot(aes(Dilutions, RAF, fill=Surface, color=Surface)) +
  geom_point(size=3) +
  geom_point(aes(y=pred_raf), data = preds, size=3, position = position_jitter(0.05, 0, 14654)) +
  geom_errorbar(
    aes(y=pred_raf, ymin=pred_raf - RAF_sd, ymax = pred_raf + RAF_sd), 
    data = preds, width=0.05, position = position_jitter(0.05, 0, 14654)
  ) +
  geom_segment(aes(pred_dil, pred_raf, xend=Dilutions), data=preds, linewidth = 1.2, linetype="dashed") +
  geom_segment(aes(pred_dil, RAF_mean, yend=pred_raf), data=preds, linewidth = 1.2, linetype="dashed") +
  geom_smooth(method="lm", fullrange=TRUE) +
  geom_text(aes(label_x, label_y, label=label), data=preds, size=6, show.legend = FALSE) +
  scale_x_continuous(limits=c(-5.1, -1.9), breaks=seq(-11, 0)) +
  coord_cartesian(ylim=c(0.019, 0.105), expand = FALSE) +
  main_theme

ggsave("bird-feeder_recovery.png", path="figures/bird-feeder", width = 12, height=10)
