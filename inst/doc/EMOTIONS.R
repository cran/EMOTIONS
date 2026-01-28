## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(EMOTIONS)

## ----eval=FALSE---------------------------------------------------------------
# install.packages("EMOTIONS")

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

## ----warning=FALSE------------------------------------------------------------

#Saving the original values
original.MY<-LacData[which(LacData$ID=="ID2" & LacData$DIM
                           %in%c(32,34,36,37,40)),
                     "DMY"]

#Masking the original values
LacData.imp<-LacData

LacData.imp[which(LacData.imp$ID=="ID2" & LacData.imp$DIM
                           %in%c(32,34,36,37,40)),
                     "DMY"]<-NA

LacData.imp<-LacData.imp[!is.na(LacData.imp$DMY),]

# Running model fitting and ensemble modeling
out.ensemble.imp <- LacCurveFit(
  data = LacData.imp, ID = "ID", trait = "DMY",
  dim = "DIM", alpha = 0.1,
  models = "All", param_list = NULL, silent=TRUE
)

#Imputing missing masked MY
out.imp<-imp_my(out.ind = out.ensemble.imp, dim = 1:305)

#Checking output
out.imp[which(out.imp$ID=="ID2" & out.imp$DIM
                           %in%c(32,34,36,37,40)),]

#Comparing with the original values
original.MY

## -----------------------------------------------------------------------------
#Creating a input file based on data frame with all ensemble predictions
out.imp.ensemb <- do.call(rbind, out.ensemble$production)

#Estimating the milk loss events
res.ensem <- milkloss_detect(
  data = out.imp.ensemb,
  id_col = "ID",
  dim_col = "DIM",
  MY_col = "DMY",
  MY_pred = "weight_AIC",
  dim_start = 11, dim_end = 294,
  drop_pct = 0.10, min_len = 1,
  rec_mode = "pctbase", rec = 0.99, stick = 2
)

#Checking the output
str(res.ensem)

## ----fig.width=17, fig.height=8-----------------------------------------------
PlotMilkLoss(data=out.imp.ensemb, 
                       ID="ID2", 
                       res.milkloss=res.ensem,
                       MY_col="DMY", 
                       MY_pred="weight_AIC",
                       col = c("red", "blue", "darkgreen"),
                       id_col = "ID")

