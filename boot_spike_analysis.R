library(quicR)
library(tidyverse)



file <- list.files("spike_raw")[[1]]

data <- get_real(file)[[0]] %>%
  normalize_RFU()

meta <- organize_tables(file) %>%
  convert_tables()

metrics <- data.frame(
  `Sample IDs` = meta$`Sample IDs`,
  Dilutions = -log10(as.numeric(meta$Dilutions)),
  check.names = FALSE
) %>%
  mutate(
    MPR = calculate_MPR(data),
    MS  = calculate_MS(data),
    TtT = calculate_TtT(data, 2),
    RAF = 1/TtT,
    # crossed = TtT != 48
  )

metrics %>%
  gather("key", "value", -c("Sample IDs", "Dilutions")) %>%
  mutate_at("Dilutions", as.factor) %>%
  ggplot(aes(`Sample IDs`, value, fill=Dilutions)) +
  geom_boxplot() +
  facet_wrap(vars(key), scales="free")
write.csv(metrics, "data/metrics.csv", row.names=FALSE)




# Titration Curve ---------------------------------------------------------



df_sum <- metrics %>%
  # filter(Dilutions != -3) %>%
  group_by_at(1:2) %>%
  summarize(
    mean_RAF = mean(RAF),
    median_RAF = median(RAF),
    sd_RAF = sd(RAF)
  )
write.csv(df_sum, "data/summary.csv", row.names = FALSE)

model <- glm(
  RAF ~ Dilutions,
  data = filter(metrics, `Sample IDs` == "P")
) %>%
  summary()

r_square <- with(model, 1 - deviance / null.deviance)

# predict.glm()



# Make the figure ---------------------------------------------------------



df_sum %>%
  filter(`Sample IDs` != "N") %>%
  ggplot(aes(
    Dilutions, median_RAF,
    ymin = median_RAF - sd_RAF, ymax = median_RAF + sd_RAF
  )) +
  geom_point(size = 3) +
  geom_errorbar(width = 0.05) +
  geom_smooth(
    method = glm,
    se = FALSE,
    linetype = "dashed",
    colour = "black"
  ) +
  geom_text(label = expression(R^2 == 0.8026), x = -6.5, y = 0.1) +
  ggtitle("Positive Titration Curve") +
  xlab("Log Dilution Factor") +
  ylab("Rate of Amyloid Formation (1/s)") +
  # scale_y_continuous(breaks = seq(0, 3e-5, 5e-6)) +
  # coord_cartesian(xlim = c(-7.2, -2.8), ylim = c(0, 3.3e-5), expand = FALSE) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.background = element_blank(),
    panel.grid = element_line(colour = "lightgrey"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    strip.background = element_blank(),
    strip.placement = "outside"
  )

ggsave("figures/titration_curve.png", width = 12, height = 6.25)  