library(quicR)
library(tidyverse)


# Main figure theme ------------------------------------------------------


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


file <- "raw/spike_raw/20250702_LW_BootRecovery2.xlsx"

df_ <- get_quic(file) %>%
  calculate_metrics("Sample IDs", "Dilutions", "Wells") %>%
  separate_wider_delim(cols=c(`Sample IDs`), delim = "_", names = c("Sample IDs", "control"), too_few = "align_start") %>%
  mutate(
    Dilutions = -log10(as.numeric(Dilutions)),
    control = !is.na(control)
  )

df_sum <- summarize(
  df_,
  sd_AUC = sd(AUC),
  AUC = mean(AUC),
  .by = c(`Sample IDs`, Dilutions, control)
)


# Generate standard curve ------------------------------------------------


df_standard <- filter(df_sum, `Sample IDs` == "P", Dilutions != -3)

model <- lm(AUC ~ Dilutions, data = df_standard)

unknowns <- df_sum %>%
  filter(!(`Sample IDs` %in% c("P", "N", "B")), !control) %>%
  summarize(
    sd_AUC = sd(AUC),
    AUC = mean(AUC),
    .by = Dilutions
  )


# Generate predictions ---------------------------------------------------


preds <- predict(model, unknowns)

std_curve <- data.frame(
  Dilutions = seq(-8, -2, 0.01),
  AUC = predict(model, list(Dilutions = seq(-8, -2, 0.01)))
)

predictions <- unknowns %>%
  mutate(
    # predicted = preds,
    real_auc = AUC,
    AUC = preds,
    loss = real_auc - AUC,
    percent_recovery = (real_auc / AUC) * 100,
    dil_factors = as.factor(Dilutions)
  )

get_pred_dil <- function (x) return (min(std_curve$Dilutions[std_curve$AUC >= x]))
dil_preds <- sapply(predictions$real_auc, get_pred_dil)

predictions$real_dilutions <- predictions$Dilutions
predictions$Dilutions <- dil_preds


# Titration Curve ---------------------------------------------------------


dotsize <- 6

df_standard %>%
  ggplot(aes(Dilutions, AUC, ymin = AUC - sd_AUC, ymax = AUC + sd_AUC)) +
  stat_smooth(aes(Dilutions, AUC), method="lm", fullrange=TRUE, se=FALSE, linetype="dashed") +
  geom_point(size = dotsize) +
  geom_errorbar(width = 0.1) +
  geom_segment(# hline
    aes(y=AUC, x=Dilutions, xend=real_dilutions, color=dil_factors), 
    data = predictions, show.legend = FALSE, linewidth=1.2
  ) +
  geom_segment(# vline
    aes(y=AUC, yend=real_auc + 60, x=Dilutions, color=dil_factors), 
    data = predictions,
    arrow = arrow(type="closed", angle=20), show.legend = FALSE, linewidth=1.2
  ) +
  geom_point(
    aes(x = real_dilutions, y = AUC, color=dil_factors), 
    predictions, size = 3, show.legend = FALSE
  ) +
  geom_errorbar(
    aes(x=Dilutions, ymin = real_auc - sd_AUC, ymax = real_auc + sd_AUC, color=dil_factors),
    predictions, width = 0.1, show.legend = FALSE
  ) +
  geom_point(
    aes(x = Dilutions, y = real_auc, color=dil_factors),
    predictions, size = dotsize, show.legend = FALSE
  ) +
  geom_text(
    aes(
      label=sprintf("Recovery = %s%%", signif(percent_recovery, 2)), 
      x = (Dilutions + real_dilutions) / 2, y = AUC + 60, color=dil_factors
    ), 
    data = predictions, size=6, show.legend = FALSE, fontface=2
  ) +
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")),size=6, fontface=2) +
  scale_color_manual(values=c("darkcyan", "darkorange")) +
  xlab("Log Dilution Factor") +
  ylab("Area Under the Curve") +
  main_theme +
  theme(
    panel.grid = element_line(colour = "lightgrey"),
  )

ggsave("figures/ppe/titration_curve.png", width = 12, height = 8)  


# Save output ------------------------------------------------------------


write.csv(df_sum, "data/spike_data/summary.csv", row.names = FALSE)
write.csv(predictions, "data/spike_data/predictions.csv", row.names = FALSE)