#' Internal function to minimize to compute the function sigma for the projections of the identified set
#'
#' @param x value at which the function is evaluated.
#' @param dir_nb the index of the considered direction.
#' @param sam0 the set of directions e where to compute the support function
#' @param kp1 the matrix of directions q, along the canonical axis, and the selected epsilon(q)
#' @param Xc_xb the possibly bootstraped/subsampled common regressor on the dataset  (Xnc,Xc). Default is NULL.
#' @param Xncb the possibly bootstraped/subsampled noncommon regressor on the dataset (Xnc,Xc). No default.
#' @param Xc_yb the possibly bootstraped/subsampled common regressor on the dataset  (Y,Xc). Default is NULL.
#' @param Yb the possibly bootstraped/subsampled outcome variable  on the dataset  (Y,Xc). No default.
#' @param values the different unique points of support of the common regressor Xc.
#' @param data_k the number of points for the grid search on epsilon. Default is 30. If NULL, then epsilon is taken fixed equal to kp.
#' @param Opt the option of the method to compute the critical values, either "boot" (numerical bootstrap) or "subsampling". Default is numerical bootstrap.
#' @param weights_x the bootstrap or sampling weights for the dataset (Xnc,Xc).
#' @param weights_y  the bootstrap or sampling weights for the dataset (Y,Xc).
#' @param c_sign sign restrictions on the commonly observed regressors: -1 for a minus sign, 1 for a plus sign, 0 otherwise. Default is NULL, i.e. no constraints.
#' @param nc_sign sign restrictions on the non-commonly observed regressors Xnc: -1 for a minus sign, 1 for a plus sign, 0 otherwise. Default is NULL, i.e. no constraints.
#' @param refs0 indicating the positions in the vector values corresponding to the components of betac.
#' @param meth  the method for the choice of epsilon, either "adapt", i.e. adapted to the direction or "min" the minimum over the directions. Default is "adapt".
#' @param boot_par the numerical bootstrap parameter. Default is 0.3 without Xc, 0.35 with Xc.
#' @param T_xy the apparent sample size the taking into account the difference in the two datasets.
#' @param bc  if TRUE compute also the bounds on betac. Default is FALSE.
#' @param winsor indicates if winsorisation. Default is FALSE.
#' @param version version of the computation of the ratio, "first" is a degraded version but fast; "second" is a correct version but slower. Default is "second".
#' @param weights_x_ref the sampling weights for the dataset (Xnc,Xc).
#' @param weights_y_ref the sampling weights for the dataset (Y,Xc).
#'
#' @return
#' the value the support function
#'
#'
#' @export
#'
#' @examples
objective_support <- function(x,dir_nb,sam0, kp1,
                              Xc_xb ,Xncb,Xc_yb,Yb,
                              values, data_k,Opt="boot",weights_x,weights_y,
                              c_sign, nc_sign,refs0,meth="adapt", boot_par, T_xy ,bc=FALSE,
                              winsor = FALSE, version="first",
                              weights_x_ref =NULL, weights_y_ref =NULL){
  # enforce the constraint q'e=1
  sam1 = sam0[dir_nb,]
  if(sum(sam1==0)>0){
    sam1[ sam1==0] <- x
  }else{
    dd = (1-sam1[1]*x)/sam1[2]
    sam1[1] <- x
    sam1[2] <- dd
  }
  sam1 <- matrix(sam1,1, dim(Xncb)[2])
  XX = Xncb%*%t(sam1)
  n_x = 1
  if(!is.null( values)){
    n_c=  dim(Xc_xb )[2]
  }else{
    n_c=0
  }
  lim=1
  # sample1 =NULL
  nb_pts=1

  ### compute point estimate
  if(is.null(values)){
    mat_var_out1 <- compute_radial(sample1 =NULL,Xc_x=NULL ,Xnc= XX,Xc_y=NULL,Y=Yb,
                                   values,n_c,n_x,
                                   nb_pts,sam0= matrix(1,1,1), kp0=kp1, data_k,Opt="boot",lim,
                                   weights_x,weights_y,
                                   c_sign, nc_sign,refs0,type="both",meth=meth,trunc=2,boot_par=boot_par,
                                   winsor= winsor, version =version)

    ##### if the bootstrap, compute the "reference" value : simply taking the original weights
    if(!is.null( weights_x_ref)){

      mat_var_out1_ref <- compute_radial(sample1 =NULL,Xc_x=NULL ,Xnc= XX,Xc_y=NULL,Y=Yb,
                                         values,n_c,n_x,
                                         nb_pts,sam0= matrix(1,1,1), kp0=kp1, data_k,Opt="boot",lim,
                                         weights_x_ref,weights_y_ref,
                                         c_sign, nc_sign,refs0,type="both",meth=meth,trunc=2,boot_par=boot_par,
                                         winsor= winsor, version =version)

    }




  }else{

    mat_var_out1 <- compute_radial(sample1 =NULL,Xc_x=Xc_xb ,Xnc= XX,Xc_y=Xc_xb,Y=Yb,
                                   values,n_c,n_x,
                                   nb_pts,sam0= matrix(1,1,1), kp0=kp1, data_k ,Opt="boot",lim,
                                   weights_x,weights_y,
                                   c_sign, nc_sign,refs0,type="both",meth=meth,trunc=2,boot_par=boot_par,
                                   winsor= winsor, version =version)
    if(!is.null( weights_x_ref)){

      mat_var_out1_ref <- compute_radial(sample1 =NULL,Xc_x=Xc_xb ,Xnc= XX,Xc_y=Xc_xb,Y=Yb,
                                         values,n_c,n_x,
                                         nb_pts,sam0= matrix(1,1,1), kp0=kp1, data_k,Opt="boot",lim,
                                         weights_x_ref,weights_y_ref,
                                         c_sign, nc_sign,refs0,type="both",meth=meth,trunc=2,boot_par=boot_par,
                                         winsor= winsor, version =version)


    }

  }

  #### changer en fonction de ce qu'on veut S, Sc, Scon.

  es =   T_xy ^(-boot_par)
  if(bc ==TRUE){

    if(!is.null( weights_x_ref)){
      res= 1/(mat_var_out1_ref$upper + sqrt( T_xy )*es*(mat_var_out1$upper- mat_var_out1_ref$upper ))

    }else{
      res = 1/mat_var_out1$upper
    }


  }else{

    if(!is.null( weights_x_ref)){
      res= 1/(mat_var_out1_ref$unconstr + sqrt( T_xy )*es*(mat_var_out1$unconstr - mat_var_out1_ref$unconstr ))
    }else{
      res = 1/mat_var_out1$unconstr
    }
  }
  if(is.na(res)){
    res=100000000
  }else{
    if(abs(res)==Inf){
      res=100000000
    }
  }
  return(res)

}
#
