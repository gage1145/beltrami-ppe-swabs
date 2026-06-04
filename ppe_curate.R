library(quicR)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(readr)
library(janitor)


# Extract the raw and calculated data. -----------------------------------


files <- list.files("raw/ppe", ".xlsx", recursive = TRUE, full.names = TRUE)

get_data <- function(file) {
  print(file)
  data <- get_quic(file)
  calcs <- calculate_metrics(data, "Sample IDs", "Wells", threshold = 4, flip_ratio=TRUE)
  return(c(list(data), list(calcs)))
}

df_list <- lapply(files, get_data)
raw <- map_dfr(df_list, ~ .[[1]]) %>%
  clean_names(replace=c("Sample IDs"="sample_id"))
analyzed <- map_dfr(df_list, ~ .[[2]]) %>%
  clean_names(replace=c("Sample IDs"="sample_id", "TtT"="ttt"))

write.csv(analyzed, "data/ppe/analyzed.csv")


# Metadata ----------------------------------------------------------------


new_metadata <- read_csv("data/ppe/new_metadata.csv", name_repair = "minimal") %>%
  clean_names(replace=c("Sample IDs"="sample_id")) %>%
  mutate(
    across(where(is.character), trimws),
    across(c(location, timing, surface, surface_location, side, sample_type), tolower),
    sample_id = as.character(sample_id),
    side = na_if(side, ""),
    location = ifelse(is.na(location), "dump site", location)
  ) 
  # relocate(c(sample_id, technician), .before = 1)

metadata <- read_csv("data/ppe/metadata.csv", name_repair = "minimal") %>%
  clean_names(replace=c("Sample IDs"="sample_id")) %>%
  mutate(
    across(where(is.character), trimws),
    sample_id = str_split_i(process_id, "-", 2),
    technician = str_sub(field_id, 1, 2)
  ) %>%
  separate_wider_delim(
    comment,
    "; ",
    names = c("surface", "condition"),
    too_few = "align_start",
    too_many = "merge"
  ) %>%
  separate_wider_regex(
    condition,
    c(zone = "[^ ]+ [^ ]+", " ", timing = ".*")
  ) %>%
  separate_wider_delim(
    timing,
    " of ",
    names = c("timing", "location"),
    too_few = "align_start"
  ) %>%
  separate_wider_delim(
    location,
    "; ",
    names = c("location", "extra_note"),
    too_few = "align_start"
  ) %>%
  separate_wider_regex(
    surface, c(side = "(?:left |right )?", surface = ".*"),
    too_few = "align_start"
  ) %>%
  separate_wider_delim(
    surface, " ",
    names = c("surface", "surface_location"), too_few = "align_start", too_many = "merge"
  ) %>%
  mutate(
    side = na_if(side, ""),
    location = ifelse(zone == "negative site", "negative site", location),
    location = ifelse(is.na(location), "dump site", location)
  ) %>%
  relocate(c(sample_id, technician), .before = 1) %>%
  bind_rows(new_metadata) %>%
  left_join(analyzed, by = "sample_id", relationship = "many-to-many") %>%
  mutate(
    across(c(location, timing, surface, surface_location, side, sample_type), tolower),
    surface_location = ifelse(surface == "boot" & surface_location != "bottom", "top", surface_location)
  )

write.csv(metadata, "data/ppe/updated_metadata.csv", row.names = FALSE)
