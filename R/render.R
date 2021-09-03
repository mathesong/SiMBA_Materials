library(rmarkdown)

render("R/SiMBA_Demonstration.Rmd", output_format = "all")

file.rename("R/SiMBA_Demonstration.md", "R/README.md")