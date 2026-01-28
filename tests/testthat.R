# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(bml)

test_check("bml")

# devtools::document()
devtools::load_all()
library(tidyverse)
data(coalgov)

bml(
  sim.y ~
    1 +
    majority +
    hm(id = id(cid), type = "RE", ar = F) +
    mm(
      id = id(pid, gid),
      vars = vars(fdep),
      fn = fn(w ~ 1 / n),
      RE = T,
      ar = F
    ),
  family = "Gaussian",
  monitor = F,
  data = coalgov
) |>
  summary()
