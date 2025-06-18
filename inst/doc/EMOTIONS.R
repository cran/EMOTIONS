## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(EMOTIONS)

## ----eval=FALSE---------------------------------------------------------------
# devtools::install_github("https://github.com/pablobio/EMOTIONS")

## -----------------------------------------------------------------------------
library(EMOTIONS)

## -----------------------------------------------------------------------------
# Load the dummy dataset
data("LacData")

# Display the first rows
head(LacData)

## ----warning=FALSE------------------------------------------------------------
# Running model fitting and ensemble modeling
out.ensemble <- LacCurveFit(
  data = LacData, ID = "ID", trait = "DMY",
  dim = "DIM", alpha = 0.1,
  models = "All", param_list = NULL, silent=TRUE
)

## -----------------------------------------------------------------------------
head(out.ensemble$converged_models$ID2)

## -----------------------------------------------------------------------------
head(out.ensemble$models_weight$ID2)

## -----------------------------------------------------------------------------
head(out.ensemble$production$ID2)

## ----fig.width=7, fig.height=9------------------------------------------------
RidgeModels(out.ensemble, metric = "AIC_rank")

## ----fig.width=9, fig.height=9------------------------------------------------
ModelRankRange(out.ensemble, metric = "AIC_rank")

## ----fig.width=7, fig.height=8------------------------------------------------
PlotWeightLac(
  data = out.ensemble, ID = "ID2",
  trait = "DMY", metric = "weight_AIC",
  dim = "DIM", col = c("red", "blue")
)

## -----------------------------------------------------------------------------
data("models_EMOTIONS")
head(models_EMOTIONS)

## ----warning=FALSE------------------------------------------------------------
out.ensemble.sub <- LacCurveFit(
  data = LacData,
  ID = "ID",
  trait = "DMY",
  dim = "DIM",
  alpha = 0.1,
  models = c("wil", "wilk", "wilycsml", "DiG", "DiGpw", "legpol3", 
             "legpol4", "legpolWil", "cubsplin3", "cubsplin4", "cubsplin5", 
             "cubsplindef", "wilminkPop", "qntReg"),
  param_list = NULL
)

## ----fig.width=7, fig.height=8------------------------------------------------
RidgeModels(out.ensemble.sub, metric = "AIC_rank")

## -----------------------------------------------------------------------------
data(model_pars)
head(model_pars)

## ----warning=FALSE------------------------------------------------------------
edited_list <- list(
  MM = c(a = 20, b = -2),
  wil = c(a = 35, b = -5, c = -0.01, k = 0.2)
)

out.ensemble.edited <- LacCurveFit(
  data = LacData,
  ID = "ID",
  trait = "DMY",
  dim = "DIM",
  alpha = 0.1,
  models = "All",
  param_list = edited_list
)

## -----------------------------------------------------------------------------
out.ensemble.edited$converged_models$ID2[["MM"]]
out.ensemble.edited$converged_models$ID2[["wil"]]

## -----------------------------------------------------------------------------
out.res <- ResInd(
  out.ensemble$production,
  dim_filter_range = c(1, 7, 203, 210),
  outlier_sd_threshold = 4,
  weight = "weight_AIC",
  trait = "DMY",
  DIM = "DIM",
  ID_col = "ID"
)

## -----------------------------------------------------------------------------
head(out.res$ri_filtered)

## -----------------------------------------------------------------------------
out.res$ri_stats

