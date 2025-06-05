library(lme4)
library(emmeans)
library(broom)
library(tidyverse)



# Function to get tidy pairwise comparisons for a given factor
get_pairwise_comparisons <- function(factor_name, model) {
  em <- emmeans(model, as.formula(paste("pairwise ~", factor_name)))
  tidy(em$contrasts) %>%
    mutate(Factor = factor_name) %>%
    {if ("adj.p.value" %in% colnames(.)) rename(., p.value = adj.p.value) else .}
}



# Load in the data --------------------------------------------------------



data <- read_csv("data/updated_metadata.csv", name_repair = "minimal") %>%
  mutate_at(c("Sample IDs", "Technician"), as.factor) %>%
  filter(
    !is.na(Technician), 
    `Sample Type` == "PPE Swab",
    Timing != "post decon" # Exclude post decon since it only occurs with dump site boots?
  ) %>%
  select(c(
    `Sample IDs`, 
    Surface, 
    `Surface Location`, 
    Location, 
    Timing, 
    Technician, 
    Side, 
    MPR, 
    MS, 
    RAF, 
    TtT
  )) %>%
  mutate_at(c("Surface", "Location", "Timing", "Technician"), as.factor) %>%
  mutate(
    Surface = relevel(Surface, "tyvek"),
    Location = relevel(Location, "negative site"),
    Timing = relevel(Timing, "pre sampling")
  )

data_sum <- data %>%
  group_by(`Sample IDs`, Surface, `Surface Location`, Location, Timing) %>%
  summarise(
    count = n(),
    mean_MPR = mean(MPR),
    mean_MS  = mean(MS),
    mean_RAF = mean(RAF),
    mean_TtT = mean(TtT)
  )



# Generate the models -----------------------------------------------------



# model_MPR <- lmer(# Model for MPR
#   MPR ~ Surface + Location + Timing + (1 | Technician),
#   # MPR ~ Surface + Location + Timing + (1 | `Sample IDs`) + (1 | Technician) + (1 | Side),
#   data = data
# )

model_MPR <- glm(
  MPR ~ Surface + Location + Timing,
  data = data
)
# summary(model_MPR)

# data.frame(
#   observed = data$MPR,
#   predicted = predict(model_MPR)
# ) %>%
#   ggplot(aes(observed, predicted)) +
#   geom_point() +
#   geom_smooth(method="glm")

model_MS <- glm(
  MS ~ Surface + Location + Timing,
  data = data
)
# model_MS <- lmer(# Model for MS
#   # MS ~ Surface + Location + Timing + (1 | `Surface Location`),
#   MS ~ Surface + Location + Timing + (1 | Technician),
#   data = data
# )

# contrast(emmeans(model_MPR, pairwise ~ Location), interaction = TRUE)

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
)
  # arrange(p.value)

write.csv(data, "data/analyzed.csv", row.names = FALSE)
write.csv(results, "data/results.csv", row.names = FALSE)



# Boot Analysis -----------------------------------------------------------



boot_data <- read_csv("data/updated_metadata.csv", name_repair = "minimal") %>%
  mutate_at(c("Sample IDs", "Technician"), as.factor) %>%
  filter(Surface == "boot", !is.na(`Surface Location`)) %>%
  select(c(
    `Sample IDs`, 
    Surface, 
    `Surface Location`, 
    Location, 
    Timing, 
    Technician, 
    Side, 
    MPR, 
    MS, 
    RAF, 
    TtT
  )) %>%
  mutate_at(c("Surface", "Location", "Timing", "Technician"), as.factor) %>%
  rename_with(~str_replace(., " ", "_")) %>%
  mutate(
    Location = relevel(Location, "negative site"),
    Timing = relevel(Timing, "pre sampling")
  )

model_MPR <- glm(MPR ~ Surface_Location + Location + Timing, data = boot_data)
model_MS <- glm(MS ~ Surface_Location + Location + Timing, data = boot_data)

factors <- c("Surface_Location", "Location", "Timing")
boot_results <- bind_rows(
  # MPR results
  lapply(factors, get_pairwise_comparisons, model = model_MPR) %>%
    bind_rows() %>%
    mutate(metric = "MPR"),
  # MS results
  lapply(factors, get_pairwise_comparisons, model = model_MS) %>%
    bind_rows() %>%
    mutate(metric = "MS")
)

write.csv(boot_data, "data/boot_analyzed.csv", row.names = FALSE)
write.csv(boot_results, "data/boot_results.csv", row.names = FALSE)



