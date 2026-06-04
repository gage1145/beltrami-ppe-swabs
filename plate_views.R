library(quicR)


# Plate Views -------------------------------------------------------------


files <- list.files("raw", recursive = TRUE, full.names = TRUE)

make_plate_view <- function(file) {
  file_list <- str_split(file, "/")[[1]]
  depth <- length(file_list)
  file_name <- file_list[depth] %>%
    str_remove(".xlsx")

  save_name <- paste0(file_name, "_plate_view.png")

  if (file.exists(paste0("figures/plate_views/", save_name))) {
    print(paste("Plate view for", file_name, "already exists."))
    next
  }
  print(paste("Working on reaction:", file_name))

  pv <- plate_view(get_quic(file), plot_deriv = FALSE) +
    ggtitle(file_name)

  ggsave(save_name, path = "figures/plate_views", width = 12, height = 8)
}

lapply(files, make_plate_view)