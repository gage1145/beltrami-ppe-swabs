box::use(
  dplyr[mutate, across],
  lubridate[date],
  magrittr[`%>%`],
  quicR[get_quic],
  stringr[str_extract, str_remove_all, str_squish]
)


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