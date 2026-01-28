#' A wrap function to the ModelsLac function that allows the fit of lactation curve models based on daily production and days in milk records simultaneously for a list of animals
#'
#' The function uses a data frame containing the daily milking records as input
#' @param data A data frame containing the daily milking records
#' @param ID The name of the column containing the unique IDs of the individuals
#' @param trait The name of the column containing daily milking records
#' @param dim  The name of the column containing days in milk records
#' @param alpha A penalization factor, ranging from 0 to 1, for the estimation of the model`s weight
#' @param models A vector describing the models to be included in the analysis. In total, 47 models are included in EMOTIONS. The default option is "All", which results in the inclusion of the 47 models. Alternatively, a vector containing any subset of the following models can be provided: "MMR","MME","brody23","brody24", "SCH","SCHL","PBE","wood","DHA", "CB","QP","CLD","PapBo1","PapBo2", "PapBo3", "PapBo4", "PapBo6", "GS1",  "GS2","LQ", "wil", "wilk", "wilycsml", "BC", "DJK","MG2", "MG4", "MG", "KHN", "AS", "FRP","PTmult","PTmod", "MonoG", "MonoGpw", "DiG", "DiGpw","legpol3", "legpol4", "legpolWil", "cubsplin3", "cubsplin4", "cubsplin5", "cubsplindef", "wilminkPop", "qntReg".
#' @param param_list A list composed by the models, named as in the models parameter, and the repective parameters included in the models.
#' @param silent A logical string defining if warning should be printed or not during the model fitting. The default value is TRUE (not printing warnings).
#'@importFrom stats na.exclude
#' @return A list containing the fitted models, the model`s weigths and ranks for each weighting strategy, and the predicted daily production obtained through the model ensemble for each weighting strategy
#' @export

LacCurveFit<-function(data,ID,trait,dim,alpha=0.5, models="All",param_list=NULL,silent=T){

  param_list<-param_list

  ids<-unique(data[,ID])

  out_main_list <- list(
    converged_models = list(),
    models_weight = list(),
    production = list()
  )

  k<-0
  for(i in ids){

    k<-k+1

    x<-data[which(data[,ID]==i),]

    out.ModelsLac<-ModelsLac(data=x,ID_col=ID,ID=i, trait=trait, dim=dim, alpha=alpha,models=models,param_list=param_list,silent=silent)

    out_main_list$converged_models[[i]]<-out.ModelsLac$converged_models[[i]]

    out_main_list$models_weight[[i]]<-out.ModelsLac$models_weight[[i]]

    out_main_list$production[[i]]<-out.ModelsLac$production[[i]]

  }

  return(out_main_list)
}
