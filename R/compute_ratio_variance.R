#' Function to compute the variance bounds
#'
#'
#' @param x a matrix containing the directions to compute the variance bounds on the radial function.
#' @param Xp the observations of the noncommon regressor (possibly conditional on Xc).
#' @param Yp the observations of the outcome variable.
#' @param dimX2  the dimension of the noncommon regressors Xnc.
#' @param weights_xp the sampling or bootstrap weights for the dataset (Xnc,Xc).
#' @param weights_yp  the sampling or bootstrap weights for the dataset (Y,Xc).
#' @param T_xy the apparent sample size taking into account the difference in the two datasets.
#' @param es the numerical bootstrap parameter delta_n.
#' @param boot_par the numerical bootstrap parameter. Default is 0.3 without Xc, 0.35 with Xc.
#' @param Opt the option of the method to compute the critical values, either "boot" (numerical bootstrap) or "subsampling". Default is numerical bootstrap.
#' @param ratio_ref the variances ratio of the point estimate for the computation of the numerical bootstrap replications. Default is NULL.
#'
#' @return
#' the value of the ratio of the variance entering the variance bounds.
#'
#'
#' @export
#'
#' @examples

compute_ratio_variance <- function(x,Xp,Yp,dimX2,weights_xp ,weights_yp, T_xy,es, boot_par, Opt="boot", ratio_ref=NULL){
   nb = x[1]
   x0 = x[-c(1)]
   x0 <- x0/Norm(x0)
   XX = Xp - matrix(rep(1,dim(Xp)[1]))%*%apply(Xp*(matrix(weights_xp,dim(Xp)[1],1)%*%rep(1,dim(Xp)[2])),2,sum)
   XX =as.matrix(XX,dim(XX)[1],dimX2)%*%matrix(x0)

   den = wtd.var(c(XX), c(weights_xp), normwt=TRUE)
   if( den==0){
      lambda=-Inf
   }else{
      lambda= sqrt(wtd.var(c(Yp), weights_yp, normwt=TRUE)/ den)
   }

   if(!is.null( ratio_ref)& Opt=="boot"){
      lambda= ratio_ref[nb] +  es*sqrt(T_xy)*( ratio_ref[nb]-  lambda)
   }
   return(pmax(lambda,0) )
}
