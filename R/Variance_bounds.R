#' Function to compute the variance bounds for Xnc
#'
#' @param Ldata dataset containing (Y,Xc) where Y is the outcome, Xc are potential common regressors
#' @param Rdata dataset containing (Xnc,Xc) where Xnc are the non commonly observed regressors, Xc are potential common regressors
#' @param out_var label of the outcome variable Y.
#' @param c_var label of the commonly observed regressors Xc.
#' @param nc_var label of the non commonly observed regressors Xnc.
#' @param c_sign sign restrictions on the commonly observed regressors: -1 for a minus sign, 1 for a plus sign, 0 otherwise. Default is NULL, i.e. no constraints.
#' @param nc_sign sign restrictions on the non-commonly observed regressors Xnc: -1 for a minus sign, 1 for a plus sign, 0 otherwise. Default is NULL, i.e. no constraints.
#' @param projections if FALSE compute the identified set along some directions or the confidence regions. Default is FALSE.
#' @param values the different unique points of support of the common regressor Xc.
#' @param sam0 the directions q to compute the radial function.
#' @param refs0 indicating the positions in the vector values corresponding to the components of betac.
#' @param nb_pts the constant C in DGM for the epsilon_0, the lower bound on the grid for epsilon, taken equal to nb_pts*ln(n)/n. Default is 1 without regressors Xc, 3 with Xc.
#' @param kp If data_k =NULL, then epsilon is taken equal to kp.
#' @param Opt the option of the method to compute the critical values, either "boot" (numerical bootstrap) or "subsampling". Default is numerical bootstrap.
#' @param nbCores number of cores for the parallel computation. Default is 1.
#' @param Bsamp the number of bootstrap/subsampling replications. Default is 1000.
#' @param list_ex the list of the data that should not be used on the different clusters. This reduces the computational cost. Default is empty list.
#' @param weights_x the sampling weights for the dataset (Xnc,Xc).
#' @param weights_y  the sampling weights for the dataset (Y,Xc).
#' @param outside if TRUE indicates that the parallel computing has been launched outside of the function. Default is FALSE.
#' @param boot_par the numerical bootstrap parameter. Default is 0.3 without Xc, 0.35 with Xc.
#' @param alpha for the level of the confidence region. Default is 0.05.
#'
#' @return a list containing, in order:
#'     - ci : a list with all the information on the confidence intervals
#'
#'     - upper: upper bound of the confidence interval on betanc at level alpha, possibly with sign constraints
#'
#'     - lower: lower bound upper bound of the confidence interval on betanc, possibly with sign constraints
#'
#'     - unconstr: confidence interval on betanc, without sign constraints
#'
#'     - betac_ci: confidence intervals on each coefficients related to the common regressor, possibly with sign constraints
#'
#'     - betac_ci_unc: confidence intervals on each coefficients related to the common regressor without sign constraints
#'
#'  - point : a list with all the information on the point estimates
#'
#'      - upper: the upper bounds on betanc, possibly with sign constraints
#'
#'      - lower: the lower bounds on betanc, possibly with sign constraints
#'
#'      -unconstr: bounds on betanc without sign constraints
#'
#'      -betac_pt: bounds on betanc, possibly with sign constraints
#'
#'      -betac_pt_unc: bounds on betanc without sign constraints
#'
#' @export
#'
#' @examples
Variance_bounds <- function(Ldata,Rdata,
                            out_var, c_var, nc_var,
                            c_sign=NULL, nc_sign=NULL,
                            projections=TRUE,
                            values,sam0, refs0,nb_pts,kp, Opt,nbCores,Bsamp=2000,
                            list_ex=NULL,weights_x = NULL,weights_y = NULL, outside=FALSE,boot_par=0.37,alpha=0.05){


  #### get sizes
  dimXc = length(c_var)
  dimXnc = length(nc_var)

  if( dimXc!=0){
    ### dataset 1
    Xc_x = as.matrix(Rdata[,c_var],dim(Rdata)[1],dimXc)
    Xnc = as.matrix(Rdata[,nc_var],dim(Rdata)[1],dimXnc)

    ### dataset 2
    Xc_y = as.matrix(Ldata[,c_var],dim(Ldata)[1],dimXc)
    Y = as.matrix(Ldata[,out_var],dim(Ldata)[1],1)
  }else{

    Xc_x = NULL
    Xnc =  as.matrix(Rdata[,nc_var],dim(Rdata)[1],dimXnc)

    ### dataset 2
    Xc_y = NULL
    Y = as.matrix(Ldata[,out_var],dim(Ldata)[1],1)

  }

  n_x = dim( Xc_x)[1]
  n_y = dim( Xc_y)[1]
  n_xy = min(n_x,n_y)
  T_xy =   n_xy

  output  <- vector("list")
  hull_point=NULL
  hull_sharp=NULL
  #####################################################################
  if(projections==FALSE){


    if(dim(sam0)[1]==1){
      sam1 = rbind(sam0,sam0)
    }else{
      sam1= sam0
    }

    if(outside ==FALSE & nbCores >1){
      ### subsampling (B samples) parallel nbCores=4
      sfInit(parallel = TRUE, cpus = nbCores, type = "SOCK")
      sfExportAll( except = list_ex)
      sfLibrary(R.matlab)
      sfLibrary(pracma)
      sfLibrary(Hmisc)
    }


    mat_var_out1<- compute_stat_variance(sample1 =NULL, Xc_x,Xnc,Xc_y,Y,values, refs0,dimXc,dimXnc,nb_pts,sam1,kp,Opt,lim =1,
                                      weights_x = weights_x,weights_y = weights_y,boot_par=boot_par,c_sign=c_sign,nc_sign=nc_sign)
    hull_point <- mat_var_out1
    if(length(c_sign)>0){
      if(sum(abs(c_sign))>0){
        effective_c_sign = c_sign
        for(j in 1:length(c_sign)){
          tests=matrix(hull_point$tests[,-c(1)],dim(sam0)[1],dim(values)[1]-1)
          if(dimXc >1){
            if(prod(!is.na(tests[,j])) ==0){
              effective_c_sign[j] = 0
            }else{
              if( prod(tests[,j]< 10^(-2))==0 & c_sign[j]==1 ){
                effective_c_sign[j] = 0
              }
            }
          }else{
            if(prod(!is.na(tests[j])) ==0){
              effective_c_sign[j] = 0
            }else{
              if( prod(tests[j]< 10^(-2))==0 & c_sign[j]==1 ){
                effective_c_sign[j] = 0
              }
            }
          }
        }
        hull_point [["effective_c_sign"]] <-  effective_c_sign
        c_sign =  effective_c_sign
      }
    }

    if(Bsamp>0){
      ##### for subsampling
      if(nbCores>1){
        res0 <- sfLapply(1:Bsamp, compute_stat_variance,X1_x=Xc_x,X2=Xnc,X1_y=Xc_y,Y,values, refs0,dimXc,
                         dimXnc,nb_pts,sam1,kp,Opt,lim =1,weights_x = weights_x,weights_y = weights_y,c_sign=c_sign,nc_sign=nc_sign, ratio_ref=hull_point[["ratio_ref"]],
                         Dmat_Yk_ref= hull_point[["DYk"]], Dmat_Xk_ref=hull_point[["DXk"]] )

      }else{
        res0 <- lapply(1:Bsamp, compute_stat_variance,X1_x=Xc_x,X2=Xnc,X1_y=Xc_y,Y,values, refs0,dimXc,
                       dimXnc,nb_pts,sam0=sam1,kp,Opt,lim =1,weights_x = weights_x,weights_y= weights_y,c_sign=c_sign,nc_sign=nc_sign, ratio_ref=hull_point[["ratio_ref"]],
                       Dmat_Yk_ref= hull_point[["DYk"]], Dmat_Xk_ref=hull_point[["DXk"]] )


      }
      mat_varb = matrix(0,Bsamp,dim(sam0)[1])
      mat_varb0 = matrix(0,Bsamp,dim(sam0)[1])
      mat_varb_unc = matrix(0,Bsamp,dim(sam0)[1])


      for(b in 1:Bsamp){
        mat_varb[b,] = res0[[b]][["upper"]]
        mat_varb0[b,] = res0[[b]][["lower"]]
        mat_varb_unc[b,] = res0[[b]][["unconstr"]]
      }


      q95 <-  function(x){quantile(x,1-alpha,na.rm=T)}
      q05 <-  function(x){quantile(x,alpha,na.rm=T)}
      q95s <-  function(x){quantile(x,1-alpha,na.rm=T)}
      q05s <-  function(x){quantile(x,alpha,na.rm=T)}

      n_x = dim(Xnc)[1]
      n_y = dim(Y)[1]

      if(Opt=="subsampling"){
        bsx = sampling_rule(n_x)
        bsy = sampling_rule(n_y)
        bs0 = sqrt(min( bsx, bsy)/  T_xy)
        mat_var_out =   hull_point[["upper"]] - apply((mat_varb -  matrix(1,dim(mat_varb)[1],1)%*% hull_point[["upper"]])*bs0,2,q05s)
        mat_varb_out_unc  =   hull_point[["unconstr"]] - apply((mat_varb_unc -   matrix(1,dim(mat_varb_unc)[1],1)%*% hull_point[["unconstr"]])*bs0  ,2,q05)
      }else{ ## numerical bootstrap
        es =    T_xy^(-boot_par)
        mat_var_out = hull_point[["upper"]] -   apply(mat_varb -  matrix(1,dim(mat_varb)[1],1)%*% hull_point[["upper"]],2,q05s)*es^(-1)/sqrt(  T_xy)
        mat_varb_out_unc  =  hull_point[["unconstr"]]- apply( mat_varb_unc -   matrix(1,dim(mat_varb_unc)[1],1)%*% hull_point[["unconstr"]],2,q05)*es^(-1)/sqrt(  T_xy)
      }
       if(Opt=="subsampling"){
        mat_var_in =   hull_point[["lower"]] - apply(mat_varb0   -  matrix(1,dim(mat_varb0)[1],1)%*% hull_point[["lower"]],2,q95s)*bs0
      }else{
        mat_var_in =    hull_point[["lower"]] - apply(mat_varb0   -  matrix(1,dim(mat_varb0)[1],1)%*% hull_point[["lower"]],2,q95s)*es^(-1)/sqrt(  T_xy)
      }

      ##
      hull_sharp[["upper"]] <-  mat_var_out
      hull_sharp[["lower"]] <-  mat_var_in
      hull_sharp[["unconstr"]] <-  mat_varb_out_unc

      ## stock replications
      hull_sharp[["upper_repli"]] <-  mat_varb
      hull_sharp[["lower_repli"]] <-  mat_varb0
      hull_sharp[["unconstr_repli"]] <-  mat_varb_unc

      #### compute the projections for beta_c
      if(Opt=="subsampling"){

        ### with sign constraints #############################################################
        if(!is.null(values)){
          beta1_pt <- compute_bnds_betac(sample1 =NULL, Opt1=Opt, info0 = mat_var_out1, values,  c_sign0 = c_sign, nc_sign0=nc_sign ,  refs0, c_var, nc_var,sam0)
          mat_beta1_l = matrix(0,Bsamp,length(refs0))
          mat_beta1_u = matrix(0,Bsamp,length(refs0))
          if(nbCores>1){
            res0b <- sfLapply(1:Bsamp,compute_bnds_betac, Opt1=Opt, info0 = res0 ,values,  c_sign0 = c_sign, nc_sign0=nc_sign ,  refs0, c_var, nc_var,sam0)
          }else{
            res0b <- lapply(1:Bsamp,compute_bnds_betac, Opt1=Opt, info0 = res0, values,  c_sign0 = c_sign, nc_sign0=nc_sign ,  refs0, c_var, nc_var,sam0)
          }
          for(b in 1:Bsamp){
            mat_beta1_l[b,] = res0b[[b]][,1]
            mat_beta1_u[b,] = res0b[[b]][,2]
          }
          beta1_l_ci =  beta1_pt[,1] +   apply(mat_beta1_l -  matrix(1,dim(mat_beta1_l)[1],1)%*% beta1_pt[,1]*bs0,2,q05)
          beta1_u_ci   = beta1_pt[,2] +   apply(mat_beta1_u -  matrix(1,dim(mat_beta1_l)[1],1)%*% beta1_pt[,2]*bs0,2,q95)
          beta1_ci <- cbind(beta1_l_ci,  beta1_u_ci )
          hull_sharp[["betac_ci"]] <-  beta1_ci
          hull_point[["betac_pt"]] <-  beta1_pt

          ### without sign constraints #############################################################
          beta1_pt <- compute_bnds_betac(sample1 =NULL, Opt1=Opt, info0 = mat_var_out1, values,  c_sign0 = c_sign, nc_sign0=nc_sign ,  refs0, c_var, nc_var,  sam0 ,constr=FALSE)
          mat_beta1_l = matrix(0,Bsamp,length(refs0))
          mat_beta1_u = matrix(0,Bsamp,length(refs0))
          if(nbCores>1){
            res0b <- sfLapply(1:Bsamp,compute_bnds_betac, Opt1=Opt, info0 = res0, values,  c_sign0 = c_sign, nc_sign0=nc_sign ,  refs0, c_var, nc_var, sam0 , constr=FALSE)
          }else{
            res0b <- lapply(1:Bsamp,compute_bnds_betac, Opt1=Opt, info0 = res0, values,  c_sign0 = c_sign, nc_sign0=nc_sign ,  refs0, c_var, nc_var, sam0 , constr=FALSE)

          }
          for(b in 1:Bsamp){
            mat_beta1_l[b,] = res0b[[b]][,1]
            mat_beta1_u[b,] = res0b[[b]][,2]
          }
          beta1_l_ci =  beta1_pt[,1] +   apply(mat_beta1_l -  matrix(1,dim(mat_beta1_l)[1],1)%*% beta1_pt[,1]*bs0,2,q05)
          beta1_u_ci   = beta1_pt[,2] +   apply(mat_beta1_u -  matrix(1,dim(mat_beta1_l)[1],1)%*% beta1_pt[,2]*bs0,2,q95)
          beta1_ci <- cbind(beta1_l_ci,  beta1_u_ci )
          hull_sharp[["betac_ci_unc"]] <-  beta1_ci
          hull_point[["betac_pt_unc"]] <-  beta1_pt
        }
      }else{
        ### with sign constraints #############################################################
        if(!is.null(values)){
          beta1_pt <- compute_bnds_betac(sample1 =NULL, Opt1=Opt, info0 = mat_var_out1, values,  c_sign0 = c_sign, nc_sign0=nc_sign ,  refs0, c_var, nc_var,sam0)
          mat_beta1_l = matrix(0,Bsamp,length(refs0))
          mat_beta1_u = matrix(0,Bsamp,length(refs0))
          if(nbCores>1){
            res0b <- sfLapply(1:Bsamp,compute_bnds_betac,Opt1=Opt, info0 = res0, values,  c_sign0 = c_sign, nc_sign0=nc_sign ,  refs0, c_var, nc_var,sam0, info1=mat_var_out1, T_xy,es)
          }else{
            res0b <- lapply(1:Bsamp,compute_bnds_betac,Opt1=Opt, info0 = res0, values,  c_sign0 = c_sign, nc_sign0=nc_sign ,  refs0, c_var, nc_var,sam0, info1=mat_var_out1,T_xy,es)
          }
          for(b in 1:Bsamp){
            mat_beta1_l[b,] = res0b[[b]][,1]
            mat_beta1_u[b,] = res0b[[b]][,2]
          }
          es =  T_xy^(-boot_par)
          term = matrix(es^(-1)/sqrt(T_xy),dim(mat_beta1_l)[1],length(refs0))
          beta1_l_ci =  beta1_pt[,1] +   apply( (mat_beta1_l -  matrix(1,dim(mat_beta1_l)[1],1)%*% beta1_pt[,1])* term,2,q05)
          beta1_u_ci   = beta1_pt[,2] +   apply( (mat_beta1_u -  matrix(1,dim(mat_beta1_l)[1],1)%*% beta1_pt[,2])*  term,2,q95)
          beta1_ci <- cbind(beta1_l_ci,  beta1_u_ci )
          hull_sharp[["betac_ci"]] <-  beta1_ci
          hull_point[["betac_pt"]] <-  beta1_pt

          ### without sign constraints #############################################################
          beta1_pt <- compute_bnds_betac(sample1 =NULL,Opt1=Opt, info0 = mat_var_out1, values,  c_sign0 = c_sign, nc_sign0=nc_sign , refs0, c_var, nc_var,  sam0, info1=NULL , T_xy=NULL,es=NULL,constr=FALSE)
          mat_beta1_l = matrix(0,Bsamp,length(refs0))
          mat_beta1_u = matrix(0,Bsamp,length(refs0))

          if(nbCores>1){
            res0b <- sfLapply(1:Bsamp,compute_bnds_betac,Opt1=Opt, info0 = res0, values,  c_sign0 = c_sign, nc_sign0=nc_sign ,  refs0, c_var,
                              nc_var, sam0, info1=mat_var_out1 , T_xy,es, constr=FALSE)
          }else{
            res0b <- lapply(1:Bsamp,compute_bnds_betac,Opt1=Opt, info = res0, values,  c_sign0 = c_sign, nc_sign0=nc_sign ,  refs0, c_var,
                            nc_var, sam0 , info1=mat_var_out1, T_xy,es, constr=FALSE)
          }

          for(b in 1:Bsamp){
            mat_beta1_l[b,] = res0b[[b]][,1]
            mat_beta1_u[b,] = res0b[[b]][,2]
          }
          beta1_l_ci   =  beta1_pt[,1] +   apply((mat_beta1_l -  matrix(1,dim(mat_beta1_l)[1],1)%*% beta1_pt[,1])*  term,2,q05)
          beta1_u_ci   =  beta1_pt[,2] +   apply((mat_beta1_u -  matrix(1,dim(mat_beta1_l)[1],1)%*% beta1_pt[,2])*  term,2,q95)
          beta1_ci <- cbind(beta1_l_ci,  beta1_u_ci )
          hull_sharp[["betac_ci_unc"]] <-  beta1_ci
          hull_point[["betac_pt_unc"]] <-  beta1_pt
        }
      }
    }else{
      hull_sharp <- matrix(NA,1,1)
    }

    if(outside ==FALSE){
      sfStop()
    }

    ########################################### dim(X) >1, compute the convex Hull using sampled or fixed directions ########################################################################""
  }else{


    if(outside ==FALSE){
      sfInit(parallel = TRUE, cpus = nbCores, type = "SOCK")
      sfExportAll()
      sfLibrary(R.matlab)
      sfLibrary(pracma)
      sfLibrary(Hmisc)
    }


    sam = NULL
    for(k in 1:dimXnc){
      sam = cbind(sam,runif(nb_pts,-10,10))
    }
    sam1 <- t(apply(sam,1,Norm))


    ### compute point estimate
    mat_var0<-  compute_stat_variance(sample1 =NULL, Xc_x,Xnc,Xc_y,Y,values, refs0,dimXc,dimXnc,nb_pts,sam1,kp,Opt,weights_x = weights_x,weights_y = weights_y)
    mat_var= mat_var0[["upper"]]*sam1
    pp <-  convhulln(mat_var, options = "Tv", output.options = "FA",return.non.triangulated.facets = FALSE)
    p0 = NULL
    for( j in 1: dimXnc){
      p0  <- cbind(p0,pp$p[pp$hull[,j],j])
    }

    hull_point <-   matrix(0,1,dim(sam0)[1])

    for(kk in 1:dim(sam0)[1]){
      hull_point [1,kk] <- max(sam0[kk,]%*%t(p0))
    }

    ##### bootstrap or subsampling
    res0 <- sfLapply(1:Bsamp, compute_stat_variance,Xc_x,Xnc,Xc_y,Y,values, refs0,dimXc,dimXnc,nb_pts,sam1,kp,Opt,weights_x = weights_x,weights_y = weights_y)
    mat_varb = matrix(0,Bsamp,dim(sam1)[1])
    for(b in 1:Bsamp){
      mat_varb[b,] = res0[[b]][[1]]
    }

    q95 <-  function(x){quantile(x,0.95,na.rm=T)}
    mat_var_out =   apply(mat_varb,2,q95)
    mat_var = sam1* mat_var_out

    pp <-  convhulln(mat_var, options = "Tv", output.options = "FA",return.non.triangulated.facets = FALSE)
    p1 = NULL
    for( j in 1: dimXnc){
      p1  <- cbind(p1,pp$p[pp$hull[,j],j])
    }
    # compute the bounds taking quantiles
    hull_sharp  <-   matrix(0,1,dim(sam0)[1])

    for(kk in 1:dim(sam0)[1]){
      hull_sharp  [1,kk] <- max(sam0[kk,]%*%t(p1))
    }

    if(outside ==FALSE){
      sfStop()
    }


  } #################################################################### end of if dimXnc ==1 or >1

  output[["ci"]] <- hull_sharp
  output[["point"]] <- hull_point


  return(output)
}
