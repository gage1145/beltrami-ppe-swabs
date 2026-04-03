library(tidyverse)



meta <- read.csv("data/bird-feeder/meta.csv", check.names = FALSE) %>%
  mutate_all(as.character)
df_ <- read.csv("data/bird-feeder/data.csv", check.names = FALSE) %>%
  filter(Date == "11/18/2024") %>%
  mutate(`Sample IDs` = str_remove(`Sample IDs`, " ")) %>%
  left_join(meta, "Sample IDs") %>%
  separate_wider_delim(ID, " ", names=c("ID", "Side", "Type"), too_few="align_start") %>%
  rename("Months" = "Length of time in env (months)") %>%
  mutate(Months = ordered(Months)) %>%
  na.omit()



summary <- df_ %>%
  group_by(`Sample IDs`, ID, Side, Type, Date, Reader, Months) %>%
  summarize_at(
    vars(MPR:RAF),
    c(mean, sd)
  ) 
  # mutate(Months = as.numeric(Months))


df_ %>%
  mutate_at("Months", as.factor) %>%
  ggplot(aes(ID, RAF, fill=Months)) +
  geom_boxplot() +
  facet_grid(rows=vars(fct_rev(Side))) +
  labs(
    y="Rate of Amyloid Formation (1/s)",
    title="Longitudinal Bird Feeder Swab RAFs"
  ) +
  theme(
    axis.title.x = element_blank(),
    plot.title = element_text(size=20, face="bold", hjust=0.5),
    strip.text = element_text(size=16),
    axis.title = element_text(size=16),
    axis.text = element_text(size=12),
    legend.title = element_text(size=16, hjust=0.5),
    legend.text = element_text(size=12)
    
  )
ggsave("feeder_box.png", path="figures/bird-feeder", width=12, height=8)

summary %>%
  mutate_at("Months", ~ as.numeric(as.character(.))) %>%
  arrange(Months) %>%
  ggplot(aes(Months, RAF_fn1, ymin=RAF_fn1 - RAF_fn2, ymax=RAF_fn1 + RAF_fn2, color=ID, group=ID)) +
  geom_point(size=3) +
  geom_errorbar(size=0.5) +
  geom_path() +
  facet_grid(cols=vars(Side)) +
  scale_x_continuous(breaks=seq(0, 15, 3), minor_breaks = seq(0,15)) +
  labs(
    x="Months in Environment",
    y="Rate of Amyloid Formation (1/s)",
    title="Longitudinal Bird Feeder Swab RAFs"
  ) +
  theme(
    # axis.title.x = element_blank(),
    plot.title = element_text(size=20, face="bold", hjust=0.5),
    strip.text = element_text(size=16),
    axis.title = element_text(size=16),
    axis.text = element_text(size=12),
    legend.title = element_text(size=16, hjust=0.5),
    legend.text = element_text(size=12)
    
  )
ggsave("feeder_line.png", path="figures/bird-feeder", width=12, height=8)
  
  
  