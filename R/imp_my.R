#' Impute missing daily milk yields using the ensemble created
#'
#' @param out.ind The list containing the data frames with the daily production records obtained from the LacCurveFit function
#' @param dim A vector with the days in milk where the milk yield will be imputed. It can contain observed and missing DIM
#' @importFrom orthopolynom polynomial.values
#' @importFrom orthopolynom legendre.polynomials
#' @return A data frame containing the imputed milk yields for the days in milk informed.
#' @export

imp_my<-function(out.ind=NULL,dim=NULL){

list_ids<-names(out.ind$models_weight)

# Apply the function to each data frame in the list
out.imp.ensemb <- lapply(seq_along(list_ids), function(i) {
  id <- list_ids[[i]]

  models.ids<-names(out.ind$converged_models[[id]])

  pred.tab<-data.frame(DIM=dim)
  for(m in models.ids){

    if(m%in%"legpol3"){

      fit<-out.ind$converged_models[[id]][[m]]

      leg3_df <- as.data.frame(polynomial.values(
        legendre.polynomials(3, normalized = TRUE),
        x = scaleX(dim, u = -1, v = 1)
      ))
      names(leg3_df) <- c("leg0","leg1","leg2","leg3")
      leg3_mat <- as.matrix(leg3_df[, c("leg1","leg2","leg3")])

      train_term <- fit$model$leg3
      leg3_fix <- leg3_mat
      class(leg3_fix) <- class(train_term)

      for (nm in setdiff(names(attributes(train_term)), c("class","dim","dimnames"))) {
        attr(leg3_fix, nm) <- attr(train_term, nm)
      }

      nd <- data.frame(leg3 = I(leg3_fix))

      pred.tab[,m] <- predict(fit, newdata =nd)

      pred.tab[,m]<-pred.tab[,m]*out.ind$models_weight[[id]][which(out.ind$models_weight[[id]][,1]==m),"BIC_weight"]

    }

    if(m %in% c("legpol4")){

      fit<-out.ind$converged_models[[id]][[m]]

      leg4_df <- as.data.frame(polynomial.values(
        legendre.polynomials(4, normalized = TRUE),
        x = scaleX(dim, u = -1, v = 1)
      ))
      names(leg4_df) <- c("leg0","leg1","leg2","leg3","leg4")
      leg4_mat <- as.matrix(leg4_df[, c("leg1","leg2","leg3","leg4")])

      train_term <- fit$model$leg4
      leg4_fix <- leg4_mat
      class(leg4_fix) <- class(train_term)

      for (nm in setdiff(names(attributes(train_term)), c("class","dim","dimnames"))) {
        attr(leg4_fix, nm) <- attr(train_term, nm)
      }

      nd <- data.frame(leg4 = I(leg4_fix))

      pred.tab[,m] <- predict(fit, newdata =nd)

      pred.tab[,m]<-pred.tab[,m]*out.ind$models_weight[[id]][which(out.ind$models_weight[[id]][,1]==m),"BIC_weight"]

    }


    if(m %in% c("legpolWil")){

      fit<-out.ind$converged_models[[id]][[m]]

      DIM_new <- dim

      leg4_df <- as.data.frame(polynomial.values(
        legendre.polynomials(4, normalized = TRUE),
        x = scaleX(DIM_new, u = -1, v = 1)
      ))
      names(leg4_df) <- c("leg0","leg1","leg2","leg3","leg4")
      leg4_mat <- as.matrix(leg4_df[, c("leg1","leg2","leg3","leg4")])

      train_term <- fit$model$leg4
      leg4_fix <- leg4_mat
      class(leg4_fix) <- class(train_term)
      for (nm in setdiff(names(attributes(train_term)), c("class","dim","dimnames"))) {
        attr(leg4_fix, nm) <- attr(train_term, nm)
      }

      nd <- data.frame(
        leg4 = I(leg4_fix),
        DIM  = DIM_new
      )

      pred.tab[,m] <- predict(fit, newdata =nd)

      pred.tab[,m]<-pred.tab[,m]*out.ind$models_weight[[id]][which(out.ind$models_weight[[id]][,1]==m),"BIC_weight"]

    }

    if(!m %in% c("legpol3","legpol4","legpolWil")){
      pred.tab[,m]<-predict(out.ind$converged_models[[id]][[m]],newdata=pred.tab)

      pred.tab[,m]<-pred.tab[,m]*out.ind$models_weight[[id]][which(out.ind$models_weight[[id]][,1]==m),"BIC_weight"]
    }

  }

  ensemb.imp<-rowSums(pred.tab[,-1])

  tab.pred.ensemb<-data.frame(ID=id,DIM=dim,MY_ensem=ensemb.imp)

  tab.prod<-out.ind$production[[id]]

  tab.pred.ensemb$MY_real<-tab.prod[match(tab.pred.ensemb$DIM,tab.prod$DIM),"weight_BIC"]

  return(tab.pred.ensemb)
})

out.imp.ensemb <- do.call(rbind, out.imp.ensemb)
out.imp.ensemb <- as.data.frame(out.imp.ensemb)

return(out.imp.ensemb)
}
