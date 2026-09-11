library(tidyverse)
library(slider)
library(ggrepel)
library(latex2exp)

main_theme <- theme(
  plot.title = element_text(size = 30, hjust = 0.5),
  axis.title = element_text(size = 20),
  axis.text = element_text(size = 16),
  strip.text = element_text(size = 20),
  legend.title = element_text(size = 20, hjust = 0.5),
  legend.text = element_text(size = 16),
  panel.grid = element_blank()
)

generate_curve <- function(x, coefs) {
  (coefs[1] / (1 + exp(coefs[2] * (coefs[3] - x)))) + 
  (coefs[4] / (1 + exp(coefs[5] * (coefs[6] - x))))
}

generate_deriv <- function(x, coefs) {
  ratio1 <- coefs[1] * coefs[2] * exp(coefs[2] * (coefs[3] - x)) / (exp(coefs[2] * (coefs[3] - x)) + 1)^2
  ratio2 <- coefs[4] * coefs[5] * exp(coefs[5] * (coefs[6] - x)) / (exp(coefs[5] * (coefs[6] - x)) + 1)^2
  ratio1 + ratio2
}

coefs <- c(20, 2, 15, -20/3, 0.3, 23.5)
threshold <- 2
mpr_color <- "red"
ms_color  <- "darkgreen"
auc_color <- "darkcyan"
ttt_color <- "black"
seg_lw <- 0.7
start <- 0
end <- 48

{
df <- tibble(
  x = seq(start, end, by = 0.1),
  y = generate_curve(x, coefs),
  d = generate_deriv(x, coefs)
)

df_sum <- tibble(
  max_y = max(df$y),
  time_to_max = df[[which(df$y == max(df$y)), "x"]],
  max_d = max(df$d),
  time_to_max_d = df[[which(df$d == max(df$d)), "x"]],
  threshold = threshold,
  time_to_threshold = min(df[which(df$y > threshold), "x"]),
  auc = integrate(generate_curve, start, end, coefs)$value
)
}

df %>%
  ggplot(aes(x)) +
  geom_area(aes(y=y), color=mpr_color, fill=auc_color, alpha=0.5, linewidth=1) +
  geom_line(aes(y=d), color=ms_color, linewidth=1) +
  geom_hline(yintercept = 0, color="black", linewidth=1) +
  geom_segment(
    x = 0,
    xend = df_sum$time_to_max,
    y = df_sum$max_y,
    yend = df_sum$max_y,
    color = mpr_color,
    linewidth = seg_lw,
    linetype = "dashed"
  ) +
  geom_segment(
    x = df_sum$time_to_max,
    xend = df_sum$time_to_max,
    y = 0,
    yend = df_sum$max_y,
    color = mpr_color,
    linewidth = seg_lw,
    linetype = "dashed"
  ) +
  geom_segment(
    x = 0,
    xend = df_sum$time_to_max_d,
    y = df_sum$max_d,
    yend = df_sum$max_d,
    color = ms_color,
    linewidth = seg_lw,
    linetype = "dashed"
  ) +
  geom_segment(
    x = df_sum$time_to_max_d,
    xend = df_sum$time_to_max_d,
    y = 0,
    yend = df_sum$max_d,
    color = ms_color,
    linewidth = seg_lw,
    linetype = "dashed"
  ) +
  geom_label_repel(
    aes(label = TeX(paste0("MPR = ", signif(max_y, 3), "\\textit{f}"), output="character"), x = time_to_max, y = max_y),
    data=df_sum, 
    hjust = 0, 
    nudge_x = df_sum$time_to_max_d - df_sum$time_to_max - 12,
    nudge_y = -2, 
    color = mpr_color, 
    size = 6, 
    segment.size = seg_lw,
    label.padding = 0.5, 
    label.r = 0.15,
    parse = TRUE
  ) +
  geom_label_repel(
    aes(
      label = TeX(paste0("MS = ", signif(max_d, 3), "$\\textit{f\\cdot{}h^{-1}}$"), output = "character"), 
      x = time_to_max_d, y = max_d
    ),
    data=df_sum, 
    hjust = 0, 
    nudge_x = -12, 
    nudge_y = 2,
    color = ms_color, 
    size = 6, 
    segment.size = seg_lw, 
    label.padding = 0.5, 
    label.r = 0.15, 
    parse = TRUE
  ) +
  geom_label(
    aes(label = TeX(paste0("AUC = ", signif(auc, 3), "\\textit{f\\cdot{}h}"), output = "character")),
    data = df_sum, 
    x = 28, 
    y = 8, 
    color = auc_color, 
    size = 6,
    parse = TRUE
  ) +
  scale_x_continuous(breaks = seq(0, 48, 4), expand = expansion()) +
  scale_y_continuous(breaks = seq(-4, 30, 2)) +
  labs(
    x = "Time (h)",
    y = "Normalized Fluorescence (f)"
  ) +
  main_theme +
  theme(
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
  )
  
ggsave("figures/metric_example.png", width = 12, height = 8)
