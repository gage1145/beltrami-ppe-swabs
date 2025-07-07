library(quicR)
library(tidyverse)



# file <- list.files("spike_raw", full.name=TRUE)[[1]]
file <- "spike_raw/20250702_LW_BootRecovery2.xlsx"

data <- get_real(file)[[1]] %>%
  normalize_RFU()

meta <- organize_tables(file) %>%
  convert_tables()

metrics <- data.frame(
  `Sample IDs` = meta$`Sample IDs`,
  Dilutions = -log10(as.integer(meta$Dilutions)),
  check.names = FALSE
) %>%
  mutate(
    MPR = calculate_MPR(data),
    MS  = calculate_MS(data),
    TtT = calculate_TtT(data, 2),
    RAF = 1/TtT,
  )

metrics %>%
  filter(
    `Sample IDs` != "B",
    # Dilutions == -3
  ) %>%
  # gather("key", "value", -c("Sample IDs", "Dilutions")) %>%
  mutate_at("Dilutions", as.factor) %>%
  ggplot(aes(`Sample IDs`, RAF, fill=Dilutions)) +
  geom_boxplot() +
  scale_fill_manual(values=c("green", "orange", "cyan", "red", "blue")) +
  # facet_wrap(vars(key), scales="free") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
ggsave("figures/recovery_RAF.png", width=8, height=8)
write.csv(metrics, "spike_data/metrics.csv", row.names=FALSE)




# Titration Curve ---------------------------------------------------------



df_sum <- metrics %>%
  # filter(Dilutions != -3) %>%
  group_by_at(1:2) %>%
  summarize(
    mean_RAF = mean(RAF),
    median_RAF = median(RAF),
    sd_RAF = sd(RAF)
  )
write.csv(df_sum, "spike_data/summary.csv", row.names = FALSE)

df_ <- metrics %>%
  filter(
    `Sample IDs` == "P",
    Dilutions != -3 # This dilution throws off the standard curve.
  )

model <- glm(
  Dilutions ~ RAF,
  data = df_
) 
  
model_sum <- summary(model)

r_square <- with(model, 1 - deviance / null.deviance)

unknowns <- metrics %>%
  filter(!(`Sample IDs` %in% c("P", "N", "B")))

sample_names <- unique(unknowns$`Sample IDs`)

predict_unknowns <- function(sample_name, data) {
  data <- filter(data, `Sample IDs` == sample_name) %>%
    group_by(`Sample IDs`, Dilutions) %>%
    summarize(
      median_RAF = median(RAF),
      sd_RAF = sd(RAF),
      RAF = median(RAF),
      Dilutions = mean(Dilutions)
    )
  split_name <- str_split(sample_name, "_")
  sample_id <- split_name[[1]][1]
  timing <- split_name[[1]][2]
  # sample_id <- sample_name
  data.frame(
    `Sample IDs` = c(sample_id),
    Timing = timing,
    real_dilutions = c(data$Dilutions),
    Dilutions = c(predict(model, data)),
    mean_RAF = data$RAF,
    median_RAF = c(data$median_RAF),
    sd_RAF = c(data$sd_RAF),
    check.names = FALSE
  ) %>%
    mutate(
      loss = Dilutions - real_dilutions,
      percent_recovery = 10^(loss) * 100,
    )
}

predictions <- lapply(
  sample_names, 
  predict_unknowns, 
  data=unknowns
) %>%
  bind_rows()

write.csv(predictions, "spike_data/predictions.csv", row.names = FALSE)



# Make the figure ---------------------------------------------------------



df_sum %>%
  filter(
    `Sample IDs` == "P",
    Dilutions != -3
  ) %>%
  ggplot(aes(
    Dilutions, mean_RAF,
    ymin = median_RAF - sd_RAF, ymax = median_RAF + sd_RAF
  )) +
  geom_smooth(
    method = glm,
    se = TRUE,
    linetype = "dashed",
    colour = "black"
  ) +
  geom_point(size = 3) +
  geom_errorbar(width = 0.05) +
  geom_point(
    aes(
      x = Dilutions, 
      y = mean_RAF,
      color = `Sample IDs`
    ), 
    predictions %>% filter(is.na(Timing) & real_dilutions == -3),
    size = 3
  ) +
  geom_errorbar(
    aes(
      ymin = median_RAF - sd_RAF, ymax = median_RAF + sd_RAF,
      color = `Sample IDs`
    ),
    predictions %>% filter(is.na(Timing) & real_dilutions == -3),
    width = 0.05
  ) +
  geom_point(
    aes(
      x = real_dilutions,
      y = mean_RAF,
      color = `Sample IDs`,
    ),
    predictions %>% filter(is.na(Timing) & real_dilutions == -3),
    shape = "triangle",
    size = 3
  ) +
  geom_hline(
    aes(
      yintercept=mean_RAF, color = `Sample IDs`
    ), 
    data = predictions %>% filter(is.na(Timing) & real_dilutions == -3),
    linetype = "dashed"
  ) +
  geom_text(label = expression(R^2 == 0.5834), x = -6.5, y = 0.15) +
  ggtitle("Positive Standard Curve") +
  xlab("Log Dilution Factor") +
  ylab("Rate of Amyloid Formation (1/s)") +
  theme_classic() +
  theme(
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
