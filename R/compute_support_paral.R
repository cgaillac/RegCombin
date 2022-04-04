#' Function to minimize to compute the function sigma for the projections of the identified set
#'
#' @param dir_nb the reference for the considered direction e in sam0
#' @param sam0 the directions q to compute the radial function.
#' @param Xnc the noncommon regressor on the dataset  (Xnc,Xc). No default
#' @param kp0 the matrix containing the directions q and the selected epsilon(q)
#' @param data_k the number of points for the grid search on epsilon. Default is 30. If NULL, then epsilon is taken fixed equal to kp.
#' @param dimXc the dimension of the common regressors Xc.
#' @param dimXnc the dimension of the noncommon regressors Xnc.
#' @param Xc_xb the possibly bootstraped/subsampled common regressor on the dataset  (Xnc,Xc). Default is NULL.
#' @param Xncb the possibly bootstraped/subsampled noncommon regressor on the dataset (Xnc,Xc). No default.
#' @param Xc_yb the possibly bootstraped/subsampled common regressor on the dataset  (Y,Xc). Default is NULL.
#' @param Yb the possibly bootstraped/subsampled outcome variable  on the dataset  (Y,Xc). No default.
#' @param values the different unique points of support of the common regressor Xc.
#' @param Opt the option of the method to compute the critical values, either "boot" (numerical bootstrap) or "subsampling". Default is numerical bootstrap.
#' @param weights_x the bootstrap or sampling weights for the dataset (Xnc,Xc).
#' @param weights_y  the bootstrap or sampling weights for the dataset (Y,Xc).
#' @param c_sign sign restrictions on the commonly observed regressors: -1 for a minus sign, 1 for a plus sign, 0 otherwise. Default is NULL, i.e. no constraints.
#' @param nc_sign sign restrictions on the non-commonly observed regressors Xnc: -1 for a minus sign, 1 for a plus sign, 0 otherwise. Default is NULL, i.e. no constraints.
#' @param refs0 indicating the positions in the vector values corresponding to the components of betac.
#' @param meth the method for the choice of epsilon, either "adapt", i.e. adapted to the direction or "min" the minimum over the directions. Default is "adapt".
#' @param boot_par the numerical bootstrap parameter. Default is 0.3 without Xc, 0.35 with Xc.
#' @param T_xy the apparent sample size the taking into account the difference in the two datasets.
#' @param bc if TRUE compute also the bounds on betac. Default is FALSE.
#' @param winsor indicates if winsorisation. Default is FALSE.
#' @param version version of the computation of the ratio, "first" is a degraded version but fast; "second" is a correct version but slower. Default is "second".
#' @param weights_xs the sampling weights for the dataset (Xnc,Xc).
#' @param weights_ys the bootstrap or sampling weights for the dataset (Y,Xc).
#'
#' @return
#' the value of the support function in the specifed direction dir_nb.
#'
#' @export
#'
#' @examples
compute_support_paral <- function(dir_nb , sam0,Xnc,  kp0, data_k,dimXc,dimXnc,Xc_xb=NULL ,Xncb,Xc_yb=NULL,Yb,
                                  values,Opt,weights_x,weights_y, c_sign,
                                  nc_sign,refs0,meth, boot_par,  T_xy ,
                                  bc,winsor, version, weights_xs, weights_ys){
  sam11 = sam0[dir_nb,]
  sam11[sam0[dir_nb,]==0] <- 1
  sam11[sam0[dir_nb,]!=0] <- 0

  ### compute starting point.
  v12 = cov(Xnc%*%sam0[dir_nb,],Xnc%*% sam11)
  v1= var(Xnc%*% sam11)
  if(dimXnc==2){
    start= -v12/v1
  }else{
    start= c(-v12/v1, rep(1,dimXnc-2))
  }


  if(dimXc==0){

    if(meth=="min"){
      kp1= kp0
    }else{
      kp1= kp0[dir_nb,1]
    }
    optf <-  optim(par= start ,  fn=objective_support, gr=NULL,   dir_nb=dir_nb,  sam0 = sam0, kp1= kp1,data_k = data_k, Xc_xb=NULL ,Xncb=Xncb,Xc_yb=NULL,Yb=Yb,
                   values=values,Opt,weights_x=weights_x,weights_y=weights_y, c_sign=c_sign,
                   nc_sign=nc_sign,refs0=refs0,meth=meth, boot_par=boot_par,  T_xy = T_xy ,
                   bc=bc,winsor = winsor, version=version,  weights_x_ref= weights_xs, weights_y_ref=weights_ys,
                   method="BFGS")
  }else{

    if(meth=="min"){
      kp1= kp0
    }else{
      kp1= kp0[dir_nb,1]
    }

    optf <- optim(par= start, objective_support, gr=NULL,   dir_nb=dir_nb,  sam0 = sam0, kp1= kp1,data_k = data_k,Xc_xb=NULL ,Xncb=Xncb,Xc_yb=NULL,Yb=Yb,
                  values=values,Opt,weights_x=weights_x,weights_y=weights_y, c_sign=c_sign,
                  nc_sign=nc_sign,refs0=refs0,meth=meth, boot_par=boot_par,  T_xy = T_xy ,
                  bc=bc,winsor = winsor, version=version,  weights_x_ref= weights_xs, weights_y_ref=weights_ys,
                  method="BFGS")

  }

  return(1/optf$value)
}
