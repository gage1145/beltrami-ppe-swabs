library(quicR)
library(dplyr)
library(stringr)
library(ggplot2)



files <- list.files("raw", ".xlsx", full.names = TRUE)
lists <- list(list(), list(), list())

for (file in files) {
  
  # Retrieve the metadata.
  meta <- organize_tables(file) %>%
    convert_tables() %>%
    mutate(
      Dilutions = -log10(as.numeric(Dilutions)),
      Date = str_split_i(file, "_", 1) %>% str_split_i("/", 2),
      Reader = str_split_i(file, "_", 2)
    )
  
  # Retrieve and normalize the raw data.
  real_list <- list()
  
  # Raw data
  real_list[[1]] <- get_real(file)[[1]] %>% transpose_real()
  
  # Normalized data
  real_list[[2]] <- normalize_RFU(real_list[[1]], transposed = TRUE)
  
  # Analyzed data
  real_list[[3]] <- calculate_metrics(real_list[[2]], meta)
  
  # Add the sample IDs, date, and reader ID to each data frame.
  for (i in 1:length(real_list)) {
    real_list[[i]] <- real_list[[i]] %>%
      mutate(
        `Sample IDs` = meta$`Sample IDs`,
        Date         = meta$Date,
        Reader       = meta$Reader,
        Dilutions    = meta$Dilutions
      ) %>%
      relocate(c(Dilutions, Date, Reader), .after = `Sample IDs`)
    lists[[i]] <- append(lists[[i]], list(real_list[[i]]))
  }
}

# Bind the data frames together and write the CSV files.
file_names <- c("data/raw.csv", "data/normalized.csv", "data/analyzed.csv")
for (i in 1:length(lists)) {
  lists[[i]] <- bind_rows(lists[[i]])
  write.csv(lists[[i]], file_names[i], row.names = FALSE)
}



# Metadata ----------------------------------------------------------------



metadata <- readr::read_csv("data/metadata.csv", name_repair = "minimal") %>%
  mutate(
    `Sample IDs` = str_split_i(`Process ID`, "-", 2)
  ) %>%
  tidyr::separate_wider_delim(
    Comment, "; ", names = c("Surface", "Condition"), too_few = "align_start", too_many = "merge"
  ) %>%
  tidyr::separate_wider_regex(
    Condition, 
    c(Zone = "[^ ]+ [^ ]+", " ", Timing = ".*") 
  ) %>%
  mutate(
    Timing = str_replace(Timing, " of dump site", "") %>% str_split_i(";", 1)
  ) %>%
  relocate(c(`Sample IDs`), .before = 1) %>%
  filter(!is.na(`Run date`)) %>%
  left_join(lists[[3]], by = "Sample IDs")

write.csv(metadata, "data/updated_metadata.csv", row.names = FALSE)

ggplot(metadata, aes(`Sample IDs`, MPR, fill = Zone)) +
  geom_boxplot() +
  facet_grid(rows = ~Timing, scales = "free_x")


