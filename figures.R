library(dplyr)
library(ggplot2)
library(quicR)




data <- readr::read_csv("data/analyzed.csv", name_repair = "minimal")


data %>%
  filter(!(`Sample IDs` %in% c("N", "P"))) %>%
  select(-c(Dilutions, Date)) %>%
  ggplot(aes(`Sample IDs`, MPR, fill = Reader)) +
  geom_boxplot() +
  scale_y_continuous(breaks = seq(2, 16, 2), limits = c(0.9, 12.1)) +
  coord_cartesian(expand = FALSE) +
  # facet_grid(cols = vars(Reader), scales = "free_x") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)
  )



# Plate Views -------------------------------------------------------------




