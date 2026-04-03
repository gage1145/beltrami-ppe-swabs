library(tidyverse)
library(quicR)
library(skimr)



# Metadata ---------------------------------------------------------------



df_meta <- read.csv("data/meta.csv", check.names=F) %>%
  separate(ID, c("Feeder", "Side"), sep=" ", extra="merge", fill="right") %>%
  separate(Side, c("Side", "Surface"), sep=" ", extra="merge", fill="left") %>%
  separate(`Period in environment`, c("start_date", "end_date"), sep=" - ", remove=F) %>%
  mutate(
    across(everything(), str_squish),
    across(c("start_date", "end_date"), ~ parse_date_time(., "my")),
    across(c("Sample IDs", "Field IDs"), as.character),
    time_in_env = end_date - start_date,
    time_in_env = replace_na(time_in_env, difftime(0,0)),
    time_since_start = end_date - min(start_date, na.rm=TRUE)
  )

files <- list.files("raw", full.names = TRUE)

# Curation function for vectorized operation.
get_raw <- function(file) {
  reader <- str_extract(file, "r\\d+")
  date <- as.Date(str_extract(file, "\\d{8}"), "%Y%m%d")
  tech <- str_remove_all(str_extract(file, "_[[:upper:]]{2,3}_"), "_")
  get_quic(file) %>%
    mutate(
      across(c("Sample IDs", "Wells", "Dilutions"), ~ str_squish(as.character(.))),
      reader = reader,
      date = as.Date(date, "%Y%m%d"),
      tech = tech
    )
}

df_raw <- map_dfr(files, get_raw)

df_calc <- calculate_metrics(
  df_raw,
  "Sample IDs", "Wells", "Dilutions", "reader", "date", "tech",
  threshold = 3
) %>%
  left_join(df_meta, by = "Sample IDs")

skim(df_calc)

write.csv(df_raw, "data/bird-feeder_raw.csv", row.names = FALSE)
write.csv(df_calc, "data/bird-feeder_data.csv", row.names = FALSE)
