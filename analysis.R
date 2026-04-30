library(lme4)
library(emmeans)
library(broom)
library(tidyverse)
library(skimr)
library(lubridate)
library(janitor)



# Load in the data --------------------------------------------------------



data <- read_csv("data/updated_metadata.csv", name_repair = "minimal")
skim(data)

df_ <- data %>%
  rename_with(tolower) %>%
  clean_names() %>%
  rename_with(~ map_vec(., ~ ifelse(. != "ms", str_remove(., "s$"), .))) %>%
  mutate(
    across(where(is.character), tolower),
    across(c(sample_id, sample_type, side, technician, surface, location, timing, well), as.factor),
    across(c(mpr, ms, auc, raf, ttt), as.numeric),
    across(c(collected_date, process_date, run_date), ~ parse_date(., "%m%.%d%.%Y")),
    timing = str_replace(timing, " ", "-"),
    surface = relevel(surface, "tyvek"),
    location = relevel(location, "negative site")
  )



# Define model variables -------------------------------------------------



metrics <- c("mpr", "ms", "auc")
factors <- c("surface", "location", "timing")
boot_factors <- c("surface_location", "location", "timing")


# Functions for generating models ----------------------------------------



get_models <- function(m, f, d) {
  terms <- paste(f, collapse = " + ")
  frmla <- as.formula(paste(m, "~", terms))
  mod <- glm(frmla, data = d)
  mod$metric <- m
  mod
}

get_pairwise_comparisons <- function(mod, f) {
  em <- emmeans(mod, as.formula(paste("pairwise ~", f)))
  tidy(em$contrasts) %>%
    mutate(Factor = f, metric = mod$metric) %>%
    {if ("adj.p.value" %in% colnames(.)) rename(., p.value = adj.p.value) else .}
}



# Generate the swab models -----------------------------------------------



df_swab <- df_ %>%
  filter(
    str_detect(sample_type, "swab"),
    # Exclude post decon since it only occurs with dump site boots.
    timing != "post-decon" 
  )

models <- map(metrics, get_models, f = factors, d = df_swab)
comps <- expand.grid(mod = models, f = factors)
results <- pmap_dfr(comps, get_pairwise_comparisons)

write.csv(df_swab, "data/analyzed.csv", row.names = FALSE)
write.csv(results, "data/results.csv", row.names = FALSE)



# Boot Analysis -----------------------------------------------------------



df_boot <- df_ %>%
  filter(surface == "boot", !is.na(surface_location))

boot_models <- map(metrics, get_models, f = boot_factors, d = df_boot)
boot_comps <- expand.grid(mod = boot_models, f = boot_factors)
boot_results <- pmap_dfr(boot_comps, get_pairwise_comparisons)

write.csv(df_boot, "data/boot_analyzed.csv", row.names = FALSE)
write.csv(boot_results, "data/boot_results.csv", row.names = FALSE)
