library(agricolae)
library(tidyverse)






data <- read_csv("data/updated_metadata.csv", name_repair = "minimal") %>%
  mutate_at(c("Sample IDs", "Technician"), as.factor) %>%
  mutate(
    Timing = factor(
      Timing, 
      levels = c(
        "pre sampling", 
        "post sampling"
      )
    ),
    Surface = factor(
      Surface,
      levels = c(
        "blank sample", 
        "tyvek leg", 
        "plastic boot", 
        "boot side", 
        "boot bottom"
      )
    ),
    # Location = ifelse(Location == "perimeter", "negative site", Location),
    Location = factor(
      Location, 
      levels = c(
        "negative site", 
        "perimeter",
        "dump site"
      )
    )
  )

df_1 <- filter(data, !is.na(`Surface`))
model1 <- aov(
  MPR ~ Surface + Location + Timing,
  data = df_1
)


post <- TukeyHSD(model1) %>%
  lapply(as.data.frame) %>%
  bind_rows()

post$trt <- rownames(post)

post %>%
  arrange(desc(`diff`)) %>%
  filter(`p adj` <= 0.05) %>%
  ggplot(aes(fct_inorder(trt), `p adj`)) +
  geom_dotplot(binaxis = "y", stackdir = "center") +
  # geom_pointrange(aes(ymin = lwr, ymax = upr)) +
  coord_flip()

# post <- LSD.test(
#   model,
#   c("Surface", "Location"),
#   group = FALSE
# )










