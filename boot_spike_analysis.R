library(quicR)
library(tidyverse)


main_theme <- theme(
  plot.title = element_text(size = 30, hjust = 0.5),
  axis.title = element_text(size = 20),
  axis.text = element_text(size = 16),
  strip.text = element_text(size = 20),
  legend.title = element_text(size = 20, hjust = 0.5),
  legend.text = element_text(size = 16),
  panel.grid = element_blank()
)


# file <- list.files("raw/spike_raw", full.name=TRUE)[[1]]
file <- "raw/spike_raw/20250702_LW_BootRecovery2.xlsx"

metrics <- get_quic(file) %>%
  calculate_metrics("Sample IDs", "Dilutions", "Wells") %>%
  mutate(Dilutions = -log10(as.numeric(Dilutions)))

metrics %>%
  filter(`Sample IDs` != "B",) %>%
  mutate_at("Dilutions", as.factor) %>%
  ggplot(aes(`Sample IDs`, AUC, fill=Dilutions)) +
  geom_boxplot() +
  scale_fill_manual(values=c("green", "orange", "cyan", "red", "blue")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
ggsave("figures/recovery_AUC.png", width=8, height=8)
write.csv(metrics, "data/spike_data/metrics.csv", row.names=FALSE)




# Titration Curve ---------------------------------------------------------



df_ <- metrics %>%
  separate_wider_delim(cols=c(`Sample IDs`), delim = "_", names = c("Sample IDs", "control"), too_few = "align_start") %>%
  mutate(control = !is.na(control))

df_sum <- df_ %>%
  summarize(
    sd_AUC = sd(AUC),
    AUC = mean(AUC),
    .by = c(`Sample IDs`, Dilutions, control)
  )

df_standard <- df_sum %>%
  filter(`Sample IDs` == "P", Dilutions != -3)

model <- lm(
  Dilutions ~ AUC,
  data = df_standard
) 

unknowns <- df_sum %>%
  filter(!(`Sample IDs` %in% c("P", "N", "B")), !control) %>%
  summarize(
    sd_AUC = sd(AUC),
    AUC = mean(AUC),
    .by = Dilutions
  )

# sample_names <- unique(unknowns$`Sample IDs`)

preds <- predict(model, unknowns)

predictions <- unknowns %>%
  mutate(
    # predicted = preds,
    real_dilutions = Dilutions,
    Dilutions = preds,
    loss = Dilutions - real_dilutions,
    percent_recovery = 10^(loss) * 100,
    dil_factors = as.factor(real_dilutions)
  )


# Make the figure ---------------------------------------------------------

std_curve <- data.frame(
  AUC = seq(0, 1500, 1),
  Dilutions = predict(model, list(AUC = seq(0, 1500, 1)))
)
dotsize <- 4

df_standard %>%
  ggplot(aes(Dilutions, AUC, ymin = AUC - sd_AUC, ymax = AUC + sd_AUC)) +
  geom_line(aes(Dilutions, AUC), data=std_curve, linetype="dashed", inherit.aes=FALSE) +
  stat_smooth(aes(AUC, Dilutions), method="lm", fullrange=TRUE) +
  geom_point(size = dotsize) +
  geom_errorbar(width = 0.05) +
  geom_point(
    aes(x = Dilutions, y = AUC, color=dil_factors), 
    predictions, size = dotsize, show.legend = FALSE
  ) +
  geom_errorbar(
    aes(ymin = AUC - sd_AUC, ymax = AUC + sd_AUC, color=dil_factors),
    predictions, width = 0.05, show.legend = FALSE
  ) +
  geom_point(
    aes(x = real_dilutions, y = AUC, color=dil_factors),
    predictions, size = dotsize, show.legend = FALSE
  ) +
  geom_text(
    aes(
      label=sprintf("Recovery = %s%%\nFold-reduction = %s", signif(percent_recovery, 2), signif(10**-loss, 2)), 
      x = Dilutions - loss / 2, y = AUC - 250, color=dil_factors
    ),
    data = predictions, size=6, show.legend = FALSE, fontface=2
  ) +
  geom_segment(
    aes(y=AUC, x=Dilutions, xend=real_dilutions, color=dil_factors), 
    data = predictions, show.legend = FALSE
  ) +
  geom_segment(
    aes(y=0, yend=AUC, x=Dilutions, color=dil_factors), 
    data = predictions,
    arrow = arrow(ends="first", type="closed", angle=20), show.legend = FALSE
  ) +
  scale_color_manual(values=c("darkcyan", "darkorange")) +
  scale_x_continuous(limits=c(-7, -3), expand=expansion(add=c(0.1, 0.1))) +
  scale_y_continuous(limits=c(0, 1200), expand=expansion(c(0, 0.01))) +
  xlab("Log Dilution Factor") +
  ylab("Area Under the Curve") +
  main_theme +
  theme(
    panel.grid = element_line(colour = "lightgrey"),
  )

ggsave("figures/titration_curve.png", width = 12, height = 6.25)  



# Save output ------------------------------------------------------------


write.csv(df_sum, "data/spike_data/summary.csv", row.names = FALSE)
write.csv(predictions, "data/spike_data/predictions.csv", row.names = FALSE)