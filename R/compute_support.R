#' Compute the support function for the projections of the identified set
#'
#'
#' @param sample1  if NULL compute the point estimate, if a natural number then evaluate a bootstrap or subsampling replication.
#' @param Xc_x the common regressor on the dataset  (Xnc,Xc). Default is NULL.
#' @param Xnc the noncommon regressor on the dataset  (Xnc,Xc). No default.
#' @param Xc_y the common regressor on the dataset  (Y,Xc). Default is NULL.
#' @param Y  the outcome variable. No default.
#' @param values  the different unique points of support of the common regressor Xc.
#' @param dimXc the dimension of the common regressors Xc.
#' @param dimXnc the dimension of the noncommon regressors Xnc.
#' @param nb_pts the constant C in DGM for the epsilon_0, the lower bound on the grid for epsilon, taken equal to nb_pts*ln(n)/n. Default is 1 without regressors Xc, 3 with Xc.
#' @param sam0 the directions q to compute the variance bounds on the radial function.
#' @param kp0 the matrix containing the directions q and the selected epsilon(q).
#' @param data_k the number of points for the grid search on epsilon. Default is 30. If NULL, then epsilon is taken fixed equal to kp.
#' @param Opt the option of the method to compute the critical values, either "boot" (numerical bootstrap) or "subsampling". Default is numerical bootstrap.
#' @param lim the limit number of observations under which we do no compute the conditional variance.
#' @param weights_x the sampling weights for the dataset (Xnc,Xc).
#' @param weights_y  the sampling weights for the dataset (Y,Xc).
#' @param c_sign sign restrictions on the commonly observed regressors: -1 for a minus sign, 1 for a plus sign, 0 otherwise. Default is NULL, i.e. no constraints.
#' @param nc_sign sign restrictions on the non-commonly observed regressors Xnc: -1 for a minus sign, 1 for a plus sign, 0 otherwise. Default is NULL, i.e. no constraints.
#' @param refs0 indicating the positions in the vector values corresponding to the components of betac.
#' @param type Equal to "both".
#' @param meth  the method for the choice of epsilon, either "adapt", i.e. adapted to the direction or "min" the minimum over the directions. Default is "adapt".
#' @param trunc equal to 2, for the definition of epsilon.
#' @param boot_par the numerical bootstrap parameter. Default is 0.3 without Xc, 0.35 with Xc.
#' @param bc  if TRUE compute also the bounds on betac. Default is FALSE.
#' @param winsor indicates if winsorisation. Default is FALSE.
#' @param version version of the computation of the ratio, "first" is a degraded version but fast; "second" is a correct version but slower. Default is "second".
#'
#' @return
#' a matrix containing the considered directions and the computed value of the support function.
#'
#' @export
#'
#' @examples
compute_support <- function(sample1 = NULL,Xc_x,Xnc,Xc_y,Y,
                            values ,dimXc,dimXnc,nb_pts,
                            sam0, kp0,data_k, Opt="subsampling",
                            lim = 30,weights_x = NULL,weights_y = NULL,
                            c_sign = NULL, nc_sign= NULL,
                            refs0=NULL,type="both",meth="adapt",trunc=2,boot_par= 0.45, bc = FALSE,
                            winsor = FALSE, version="first"){
  mat_var_low= NULL
  if(is.null(weights_x)){
    weights_x= rep(1/dim(Xnc)[1],dim(Xnc)[1])
  }
  if(is.null(weights_y)){
    weights_y= rep(1/length(Y),length(Y))
  }

  weights_xs <-  weights_x
  weights_ys <-  weights_y


  if(!is.null(sample1)){
    if(Opt == "subsampling"){
      n_x = dim(Xnc)[1]
      bsx = sampling_rule(n_x)
      bb = sample(1:n_x,bsx, replace=FALSE)
      if(!is.null(Xc_x)){
        Xc_xb = matrix(Xc_x[bb,],bsx,dimXc)
      }
      Xncb = matrix(Xnc[bb,],bsx,dimXnc)
      weights_x =  matrix(weights_x[bb],bsx,1)
      weights_x = weights_x/sum(weights_x)

      n_y = dim(Y)[1]
      bsy = sampling_rule(n_y)
      bby = sample(1:n_y,bsy, replace=FALSE)
      if(!is.null(Xc_y)){
        Xc_yb = matrix(Xc_y[bby,],bsy,dimXc)
      }
      Yb = matrix(Y[bby,],bsy,1)
      weights_y =  matrix(weights_y[bby],bsy,1)
      weights_y = weights_y/sum(weights_y)

    }else{

      n_x = dim(Xnc)[1]
      n_y = dim(Y)[1]
      es = min(n_y,n_x)^(-boot_par)  ## epsilon_n
      drawsx = rmultinom(1, n_x, rep(1/n_x,n_x))
      prob0 =  drawsx
      prob0 = prob0/sum( prob0,na.rm=T)

      if(!is.null(Xc_x)){
        Xc_xb =matrix(Xc_x,n_x,dimXc)
      }
      Xncb = matrix(Xnc,n_x,dimXnc)
      weights_x =  matrix(weights_x*prob0,n_x,1)
      weights_x = weights_x/sum(weights_x)

      drawsy = rmultinom(1, n_y, rep(1/n_y,n_y))
      prob0 =drawsy
      prob0 = prob0/sum( prob0,na.rm=T)
      if(!is.null(Xc_y)){
        Xc_yb = matrix(Xc_y,n_y,dimXnc)
      }
      Yb =matrix(Y,n_y,1)
      weights_y =  matrix(weights_y*prob0,n_y,1)
      weights_y = weights_y/sum(weights_y)


    }
  }else{
    ## point estimate
    weights_xs = NULL
    weights_ys = NULL
    Xc_xb =Xc_x
    Xncb = Xnc
    Xc_yb =   Xc_y
    Yb = Y

  }

  dir_nb=1
  sam1 = matrix(NA,dim(sam0)[1],1)
  varn = function(x){var(x,na.rm=TRUE)}
  bnd =1.5*sqrt( var(Yb,na.rm=TRUE)/min(apply(Xncb,2,varn)))

  n_x = dim(Xnc)[1]
  n_y = dim(Y)[1]
  T_xy=n_x*(n_y/(n_x+n_y))

  nbCores_dir = 1

  if(nbCores_dir==1){
    out1 = lapply(1:dim(sam0)[1],compute_support_paral,sam0,Xnc,  kp0, data_k,dimXc,dimXnc,Xc_xb=NULL ,Xncb,Xc_yb=NULL,Yb,
                  values,Opt,weights_x,weights_y, c_sign,
                  nc_sign,refs0,meth, boot_par,  T_xy  ,
                  bc, winsor, version, weights_xs, weights_ys)
  }else{
    out1 =  sfLapply(1:dim(sam0)[1],compute_support_paral,sam0,Xnc,  kp0, data_k,dimXc,dimXnc,Xc_xb=NULL ,Xncb,Xc_yb=NULL,Yb,
                     values,Opt,weights_x,weights_y, c_sign,
                     nc_sign,refs0,meth, boot_par,  T_xy ,
                     bc,winsor, version, weights_xs, weights_ys)
  }
  out11 <- unlist(out1)

  sam1 = cbind(sam0, out11)

  return(sam1)

}
