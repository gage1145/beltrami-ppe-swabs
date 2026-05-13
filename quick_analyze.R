scripts <- pmap(
  list(pattern = c("curate.R", "analysis.R", "figures.R")),
  list.files
)

map(
  1:3,
  function(i) {
    sprintf("Now running %s.\n", scripts[[i]])
    sapply(scripts[[i]], source)
    rm(list = ls())
  }
)
