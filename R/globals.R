# ================================================================================================ #
# Global variable declarations for NSE (non-standard evaluation)
# ================================================================================================ #

# Declare global variables used in tidyverse NSE contexts to avoid R CMD check NOTEs
# These variables are used in dplyr and tidyr functions with unquoted column names

utils::globalVariables(c(
  # mcmcDiag.R
  "Parameter",
  "Chain",
  ".",

  # summary.bml.R
  "where",

  # monetPlot.R
  "value",

  # createData.R
  "mmid",
  "mainid",
  "hmid",
  "val",
  "rn",
  "mmn",
  "mmi1",
  "mmi2",

  # createJagsVars.R
  "n",
  "ev",

  # formatJags.R
  "name",
  "sd",
  "lb",
  "ub",
  "j",
  "hmname",
  "hmlabel",
  "component",

  # helpers.R (mm() default weight helper)
  "importance",

  # tidy.bml.R / bmlCompare.R / coefPlot.R
  "term",
  "estimate",
  "conf.low",
  "conf.high",
  "model",
  "N",
  "DIC"
))
