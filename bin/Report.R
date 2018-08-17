# Usage: Rscript Report.R <report.rmd> <out.html>
args = commandArgs(T)

library(magrittr)
library(kableExtra)
rmarkdown::render(args[1],output_file=args[2])
