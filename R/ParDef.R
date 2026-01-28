#' Define the parameters for the lactation curve models to be fitted
#'
#' @param x  data frame containing the daily milking records for the individual
#' @param trait The name of the column containing daily milking records
#' @param dim The name of the column containing days in milk records
#' @return A list containing the parameters to be included in the lactation curve models
#' @keywords internal
ParDef<-function(x,trait,dim){

  par_list<-list(MM=c(a = 19.8, b = -1.65),
                 MMR=c(a = -10, b = -1.4),
                 MME=c(a = -0.06608, b = 317.49, c = -328.06, d = -0.027),
                 brody23=c(a = 25.6, b = 0.0015),
                 brody24=c(a = 26.127, b = 0.0017, c = 0.2575),
                 SCH=c(a = 1, b = 16),
                 SCHL=c(a = min(x[,trait]), b = -2.6),
                 PBE=c(a = 25, b = -8e-04, c = 2e-06),
                 wood=c(a = median(x[,trait]), b = log(max(x[,trait]) / median(x[,trait])) / max(x[,dim]), c = 0.05),
                 DHA=c(a = median(x[,trait]), b = log(max(x[,trait]) / median(x[,trait])) / max(x[,dim]), c = 0.05),
                 CB=c(a = min(x[,trait]), b = 0.07, c = 0.002),
                 QP=c(a = 25, b = -0.02, c = -1.5e-05),
                 CLD=c(a = 25, b = 0.03, c = 0.28),
                 PapBo1=c(a = min(x[,trait]), b = -0.0098, c = 0.0033),
                 PapBo2=c(a = 23.896, b = 0.398, c = 0.0034),
                 PapBo3=c(a = 15.27, b = 2.587, c = 0.00346),
                 PapBo4=c(a = 1.763, b = 76690, c = 0.00219),
                 PapBo6=c(a = 17.25, b = 0.493, c = 0.0018),
                 GS1=c(a = 18.28, b = -1.58, c = 4.33),
                 GS2=c(a = 17.75, b = -1.72, c = 4.81, d = 7e-11),
                 LQ=c(a = min(x[,trait]), b = (max(x[,trait], na.rm = TRUE) - min(x[,trait], na.rm = TRUE)) / max(x[,dim], na.rm = TRUE), c = (max(x[,trait], na.rm = TRUE) - min(x[,trait], na.rm = TRUE)) / 2),
                 wil=c(a = 25, b = -7, c = -0.03, k = 0.1),
                 wilk=c(a = 25, b = -7, c = -3, k = 0.1),
                 wilycsml=c(a = 20, b = 30, c = 0.005, k = 0.08),
                 BC=c(a = min(x[,trait], na.rm = TRUE), b =0.0017, c = -7.68, d = 0.08),
                 DJK=c(a = 19, b = 0.027, c = 0.08, d = 0.0017),
                 MG2=c(a = 3, b = -0.18, c = 0.002, d = -1.4),
                 MG4=c(a = 20.5, b = -0.38, c = -1.45, d = -0.005),
                 MG=c(a = 3, b = 0.18, c = 0.002, d = -1.4),
                 KHN=c(a = 15, b = -0.15, c = 0.00043, d = 5e-07, f = 4.05),
                 AS=c(a = 19, b = -5, c = 0.003, d = 5.3, f = -1.1),
                 FRP=c(a = 15, b = -0.08, c = 8, d = -3, f = -4e-05),
                 PTmult=c(a = 25, b = -0.001),
                 PTmod=c(a = min(x[,trait],na.rm = TRUE), b = (max(x[,trait], na.rm = TRUE) - min(x[,trait], na.rm = TRUE)) / max(x[,dim], na.rm = TRUE), c = 0.13, d = -0.001),
                 MonoG=c(a = 7703, b = 0.002, c = 277),
                 MonoGpw=c(a = 24.5707, b = 0.6292, c = 1.9977, k = 0.1978),
                 DiG=c(a = 353.8, b = 0.016, c = 36.3, d = 7371, f = 0.003, g = 113.5),
                 DiGpw=c(a = 15.99, b = 0.3351, c = 4.693, d = 8687, f = 0.00226, g = 80.33, k = 0.435),
                 legpolWil=c(a= -0.8, b = -0.6, c = 0.1, d = 25.7, k = 0.002),
                 cubsplindef=c(knots="c(49, 78, 112, 157, 210)"),
                 wilminkPop=c(start_values = list(
                   b0 = min(x[,trait], na.rm = TRUE),
                   b1 = (max(x[,trait], na.rm = TRUE) - min(x[,trait], na.rm = TRUE)) / max(x[,dim], na.rm = TRUE),
                   b2 = (max(x[,trait], na.rm = TRUE) - min(x[,trait], na.rm = TRUE)) / 2
                 )),
                 Legpol4Poppe=c(b0 = 0.5, b1 = 0.5, b2 = 0.5, b3 = 0.5, b4 = 0.5)
                 )

  return(par_list)
}
