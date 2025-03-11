library(quicR)
library(dplyr)
library(stringr)
library(ggplot2)



files <- list.files("raw", ".xlsx", full.names = TRUE)

df_ <- list()

for (file in files) {
  data <- get_real(file)[[1]] %>%
    normalize_RFU()
  
  meta <- organize_tables(file) %>%
    convert_tables()
  
  reader <- str_split_i(file, "_", 2)
  date <- str_split_i(file, "_", 1) %>%
    str_split_i("/", 2)
  
  analyzed <- calculate_metrics(data, meta) %>%
    mutate(
      Dilutions = -log10(as.numeric(Dilutions)),
      Date = date,
      Reader = reader
    )
  
  df_ <- append(df_, list(analyzed))
}

df_ <- bind_rows(df_) %>%
  relocate(c(Date, Reader), .before = MPR)

write.csv(df_, "data/data.csv", row.names = FALSE)









