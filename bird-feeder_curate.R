library(tidyverse)
library(quicR)
library(skimr)



# Metadata ---------------------------------------------------------------



df_meta <- read.csv("data/bird-feeder/meta.csv", check.names=F) %>%
  filter(!is.na(`Sample IDs`)) %>%
  separate(ID, c("feeder", "side"), sep=" ", extra="merge", fill="right") %>%
  separate(side, c("side", "surface"), sep=" ", extra="merge", fill="left") %>%
  separate(`Period in environment`, c("start_date", "end_date"), sep=" - ", remove=F) %>%
  mutate(
    across(everything(), str_squish),
    control = str_detect(tolower(feeder), "ctrl"),
    site = ifelse(control, "control", ifelse(feeder %in% c("BF1", "BF2", "BF3"), "perimeter", "dump site")),
    across(c("start_date", "end_date"), ~ parse_date_time(., "my")),
    across(c("Sample IDs"), as.character),
    time_in_env = end_date - start_date,
    time_in_env = replace_na(time_in_env, difftime(0, 0, units="days")),
    time_since_start = difftime(end_date, min(start_date, na.rm=TRUE), units="days")
  )

files <- list.files("raw/bird-feeder", full.names = TRUE)

# Curation function for vectorized operation.
get_raw <- function(file) {
  reader <- str_extract(file, "r\\d+")
  rxn_date <- as.Date(str_extract(file, "\\d{8}"), "%Y%m%d")
  tech <- str_remove_all(str_extract(file, "_[[:upper:]]{2,3}_"), "_")
  get_quic(file) %>%
    mutate(
      across(c("Sample IDs", "Wells", "Dilutions"), ~ str_squish(as.character(.))),
      reader = reader,
      rxn_date = rxn_date,
      tech = tech,
      reextracted = rxn_date == date("2025-01-03")
    )
}

df_raw <- map_dfr(files, get_raw)

df_calc <- calculate_metrics(
  df_raw,
  "Sample IDs", "Wells", "Dilutions", "reader", "rxn_date", "tech", "reextracted",
  threshold = 3
) %>%
  left_join(df_meta)

skim(df_calc)

write.csv(df_raw, "data/bird-feeder/raw.csv", row.names = FALSE)
write.csv(df_calc, "data/bird-feeder/data.csv", row.names = FALSE)
