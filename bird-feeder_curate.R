# Import Functions -------------------------------------------------------


box::use(
  ./R/birdfeeder,
  dplyr[mutate, filter, across, rename, left_join, select, everything],
  lubridate[date, parse_date_time],
  magrittr[`%>%`],
  purrr[map_dfr],
  quicR[calculate_metrics],
  skimr[skim],
  stringr[str_detect,str_extract, str_remove_all, str_squish],
  tidyr[separate, pivot_longer],
)


# Metadata ---------------------------------------------------------------


height_key <- c(
  BF1 = 5, BF2 = 5, BF3 = 5, BF4 = 1, BF5 = 5, BF6 = 2, "PL18-1" = 0, "PL18-2" = 0
)

site_key <- c(
  Ctrl1 = "control", Ctrl2 = "control", Ctrl3 = "control",
  BF1 = "perimeter", BF2 = "perimeter", BF3 = "perimeter",
  BF4 = "dump site", BF5 = "dump site", BF6 = "dump site",
  "PL18-1" = "dump site", "PL18-2" = "dump site"
)

df_meta <- read.csv("data/bird-feeder/meta.csv", check.names = F) %>%
  filter(!is.na(`Sample IDs`)) %>%
  separate(ID, c("feeder", "side"), sep = " ", extra = "merge", fill = "right") %>%
  separate(side, c("side", "surface"), sep = " ", extra = "merge", fill = "left") %>%
  separate(`Period in environment`, c("start_date", "end_date"), sep = " - ", remove = F) %>%
  mutate(
    across(everything(), str_squish),
    control = str_detect(tolower(feeder), "ctrl"),
    site = site_key[feeder],
    across(c("start_date", "end_date"), ~ parse_date_time(., "my")),
    across(c("Sample IDs"), as.character),
    time_in_env = difftime(end_date, start_date, units = "days"),
    time_since_start = difftime(end_date, min(start_date, na.rm = TRUE), units = "days"),
    height = height_key[feeder]
  )


# Import and analyze raw data --------------------------------------------


files <- list.files("raw/bird-feeder", full.names = TRUE)
df_raw <- map_dfr(files, birdfeeder$get_raw)

df_calc <- calculate_metrics(
  df_raw,
  "Sample IDs", "Wells", "Dilutions", "reader", "rxn_date", "tech", "reextracted",
  threshold = 3
) %>%
  left_join(df_meta)

skim(df_calc)

write.csv(df_raw, "data/bird-feeder/raw.csv", row.names = FALSE)
write.csv(df_calc, "data/bird-feeder/data.csv", row.names = FALSE)
