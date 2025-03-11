library(dplyr)
library(ggplot2)




data <- read.csv("data/data.csv")


data %>%
  ggplot(aes(`Sample IDs`, MPR)) +
  geom_boxplot() +
  facet_grid(cols = vars(Reader), scales = "free_x")