library(lme4)
library(emmeans)
library(broom)
library(tidyverse)



# Load in the data --------------------------------------------------------



data <- read_csv("data/updated_metadata.csv", name_repair = "minimal") %>%
  mutate_at(c("Sample IDs", "Technician"), as.factor) %>%
  select(c(
    `Sample IDs`, Surface, Location, Timing, Technician, Side, MPR, MS, RAF, TtT
  )) %>%
  mutate_at(c("Surface", "Location", "Timing", "Technician"), as.factor) %>%
  filter(!is.na(Technician))

data_sum <- data %>%
  group_by(Surface, Location, Timing) %>%
  summarise(
    count = n(),
    mean_MPR = mean(MPR),
    mean_MS  = mean(MS),
    mean_RAF = mean(RAF),
    mean_TtT = mean(TtT)
  )



# Generate the models -----------------------------------------------------



model_MPR <- lmer(# Model for MPR
  MPR ~ Surface + Location + Timing + (1 | Technician), 
  # MPR ~ Surface + Location + Timing + (1 | `Sample IDs`) + (1 | Technician) + (1 | Side), 
  data = data
)
summary(model_MPR)

data.frame(
  observed = data$MPR,
  predicted = predict(model_MPR)
) %>%
  ggplot(aes(observed, predicted)) +
  geom_point() +
  geom_smooth(method="lm")
  # scale_x_continuous(limits = c(0,10)) +
  # scale_y_continuous(limits = c(0,10))

model_MS <- lmer(# Model for MS
  MS ~ Surface + Location + Timing + (1 | Technician) + (1 | Side),
  # MS ~ Surface + Location + Timing + (1 | `Sample IDs`) + (1 | Technician) + (1 | Side), 
  data = data
)

# Function to get tidy pairwise comparisons for a given factor
get_pairwise_comparisons <- function(factor_name, model) {
  em <- emmeans(model, as.formula(paste("pairwise ~", factor_name)))
  tidy(em$contrasts) %>%
    mutate(Factor = factor_name) %>%
    {if ("adj.p.value" %in% colnames(.)) rename(., p.value = adj.p.value) else .}
}

# Apply to each fixed effect
factors <- c("Surface", "Location", "Timing")

results <- bind_rows(
  # MPR results
  lapply(factors, get_pairwise_comparisons, model = model_MPR) %>%
    bind_rows() %>%
    mutate(metric = "MPR"),
  # MS results
  lapply(factors, get_pairwise_comparisons, model = model_MS) %>%
    bind_rows() %>%
    mutate(metric = "MS")
) %>%
  arrange(p.value)

write.csv(data, "data/analyzed.csv", row.names = FALSE)
write.csv(results, "data/results.csv", row.names = FALSE)


results %>%
  arrange(desc(p.value)) %>%
  mutate(p.value = signif(p.value, 4)) %>%
  ggplot(
    aes(
      fct_inorder(contrast), 
      estimate, 
      ymin = estimate - std.error, 
      ymax = estimate + std.error,
      color = p.value < 0.05
    )
  ) +
  geom_pointrange(
    position = position_dodge(width = 0.4), size = 1, linewidth = 1
  ) +
  geom_hline(yintercept = 0, linewidth = 1.2, linetype = "dashed") +
  geom_label(aes(label = p.value), nudge_x = 0.3) +
  facet_grid(term ~ metric, scale = "free", space = "free_y") +
  scale_colour_manual(values = c("navy", "darkred")) +
  # scale_y_continuous(limits = c(-1, 1)) +
  coord_flip() +
  theme(
    axis.title.y = element_blank(),
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  )
ggsave("results.png", path = "figures", width = 8, height = 12)
  


# Random Effects ----------------------------------------------------------


tech_effects <- bind_rows(
  tech_effects_MPR <- ranef(model_MPR)$Technician %>%
    rename(intercept_shift = `(Intercept)`) %>%
    tibble::rownames_to_column("technician") %>%
    mutate(metric = "MPR"),
  
  tech_effects_MS <- ranef(model_MS)$Technician %>%
    rename(intercept_shift = `(Intercept)`) %>%
    tibble::rownames_to_column("technician") %>%
    mutate(metric = "MS")
)

# Plot the intercept shifts
tech_effects %>%
  ggplot(aes(x = reorder(technician, intercept_shift), y = intercept_shift, fill = metric)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_y_continuous(limits = c(-2, 2)) +
  coord_flip(expand = FALSE) +
  labs(
    title = "Random Intercepts by Technician",
    x = "Technician",
    y = "Intercept Shift (from overall mean)"
  ) +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  )
ggsave("tech_effects.png", path = "figures", width = 14, height = 7)
