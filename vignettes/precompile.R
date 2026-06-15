# Rebuild examples.Rmd from examples.Rmd.orig.
#
# The models are NOT fit here: they are fit once by fit_models.R and cached in
# vignettes/.fits/ (git-ignored). examples.Rmd.orig loads those cached fits, so
# this step runs in seconds and needs neither JAGS nor the CILS4EU data. Re-run
# fit_models.R only when a model's specification, data, or seed changes.
#
# Run from the repository root:  Rscript -e 'source("vignettes/precompile.R")'
old <- setwd("vignettes")
knitr::knit("examples.Rmd.orig", "examples.Rmd")
setwd(old)
