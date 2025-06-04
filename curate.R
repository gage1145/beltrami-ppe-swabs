library(quicR)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(readr)



files <- list.files("raw", ".xlsx", full.names = TRUE)
lists <- list(list(), list(), list())

for (file in files) {
  
  print(file)
  
  # Retrieve the metadata.
  meta <- organize_tables(file) %>%
    convert_tables() %>%
    mutate(
      `Sample IDs` = trimws(`Sample IDs`, "left"),
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
        `Rxn Date`   = meta$Date,
        Reader       = meta$Reader,
        Dilutions    = meta$Dilutions
      ) %>%
      relocate(c(Dilutions, `Rxn Date`, Reader), .after = `Sample IDs`)
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



new_metadata <- read_csv("data/new_metadata.csv", name_repair = "minimal") %>%
  mutate(`Sample IDs` = as.character(`Sample IDs`))

metadata <- read_csv("data/metadata.csv", name_repair = "minimal") %>%
  mutate(
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
    Surface, c(Side = "(?:left |right )?", Surface = ".*"), too_few = "align_start" 
  ) %>%
  separate_wider_delim(
    Surface, " ", names = c("Surface", "Surface Location"), too_few = "align_start", too_many = "merge"
  ) %>%
  mutate(
    Side = na_if(Side, ""),
    Location = ifelse(Zone == "negative site", "negative site", Location),
    Location = ifelse(is.na(Location), "dump site", Location)
  ) %>%
  relocate(c(`Sample IDs`, Technician), .before = 1) %>%
  bind_rows(new_metadata) %>%
  # filter(!is.na(`Run date`)) 
  left_join(lists[[3]], by = "Sample IDs") %>%
  mutate_all(trimws)


write.csv(metadata, "data/updated_metadata.csv", row.names = FALSE)




