# Run manually and locally to rebuild examples.Rmd from examples.Rmd.orig.
# Requires the CILS4EU feathers at local directory (outside the repo, not committed).
old <- setwd("vignettes")
knitr::knit("examples.Rmd.orig", "examples.Rmd")
setwd(old)