library(quicR)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(readr)


files <- list.files("raw/ppe", ".xlsx", recursive = TRUE, full.names = TRUE)
lists <- list(list(), list())

for (file in files) {
  print(file)
  data <- get_quic(file)
  calcs <- calculate_metrics(data, "Sample IDs", "Wells", threshold = 3, flip_ratio=TRUE)

  lists[[1]] <- append(lists[[1]], list(data))
  lists[[2]] <- append(lists[[2]], list(calcs))
}

raw <- bind_rows(lists[[1]])
analyzed <- bind_rows(lists[[2]])

write.csv(analyzed, "data/ppe/analyzed.csv")


# Metadata ----------------------------------------------------------------


new_metadata <- read_csv("data/ppe/new_metadata.csv", name_repair = "minimal") %>%
  mutate(
    across(where(is.character), ~ tolower(trimws(.))),
    `Sample IDs` = as.character(`Sample IDs`),
    Side = na_if(Side, ""),
    Location = ifelse(is.na(Location), "dump site", Location)
  ) %>%
  relocate(c(`Sample IDs`, Technician), .before = 1)

metadata <- read_csv("data/ppe/metadata.csv", name_repair = "minimal") %>%
  mutate(
    across(where(is.character), ~ tolower(trimws(.))),
    `Sample IDs` = str_split_i(`Process ID`, "-", 2),
    Technician = str_sub(`Field ID`, 1, 2)
  ) %>%
  separate_wider_delim(
    Comment,
    "; ",
    names = c("Surface", "Condition"),
    too_few = "align_start",
    too_many = "merge"
  ) %>%
  separate_wider_regex(
    Condition,
    c(Zone = "[^ ]+ [^ ]+", " ", Timing = ".*")
  ) %>%
  separate_wider_delim(
    Timing,
    " of ",
    names = c("Timing", "Location"),
    too_few = "align_start"
  ) %>%
  separate_wider_delim(
    Location,
    "; ",
    names = c("Location", "extra_note"),
    too_few = "align_start"
  ) %>%
  separate_wider_regex(
    Surface, c(Side = "(?:left |right )?", Surface = ".*"),
    too_few = "align_start"
  ) %>%
  separate_wider_delim(
    Surface, " ",
    names = c("Surface", "Surface Location"), too_few = "align_start", too_many = "merge"
  ) %>%
  mutate(
    Side = na_if(Side, ""),
    Location = ifelse(Zone == "negative site", "negative site", Location),
    Location = ifelse(is.na(Location), "dump site", Location)
  ) %>%
  relocate(c(`Sample IDs`, Technician), .before = 1) %>%
  bind_rows(new_metadata) %>%
  left_join(analyzed, by = "Sample IDs", relationship = "many-to-many") %>%
  mutate(
    `Surface Location` = ifelse(Surface == "boot" & `Surface Location` != "bottom", "top", `Surface Location`)
  )

write.csv(metadata, "data/ppe/updated_metadata.csv", row.names = FALSE)
