library(here)
library(rmarkdown)

setwd(here::here("R"))

render("SiMBA_Demonstration.Rmd", output_format = "all")

file.rename("SiMBA_Demonstration.md", "README.md")

setwd(here::here())
