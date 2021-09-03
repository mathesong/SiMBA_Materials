library(rmarkdown)

render("R/SiMBA_Demonstration.Rmd")

file.rename("R/SiMBA_Demonstration.md", "R/README.md")