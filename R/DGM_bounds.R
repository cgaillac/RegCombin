#' This function compute the DGM bounds for all the different coefficients.
#'
#' @param Ldata dataset containing (Y,Xc) where Y is the outcome, Xc are potential common regressors.
#' @param Rdata dataset containing (Xnc,Xc) where Xnc are the non commonly observed regressors, Xc are potential common regressors.
#' @param values the different unique points of support of the common regressor Xc.
#' @param sam0 the directions q to compute the radial function.
#' @param refs0 indicating the positions in the vector values corresponding to the components of betac.
#' @param out_var label of the outcome variable Y.
#' @param nc_var label of the non commonly observed regressors Xnc.
#' @param c_var label of the commonly observed regressors Xc.
#' @param nc_sign sign restrictions on the non-commonly observed regressors Xnc: -1 for a minus sign, 1 for a plus sign, 0 otherwise. Default is NULL, i.e. no constraints.
#' @param c_sign sign restrictions on the commonly observed regressors: -1 for a minus sign, 1 for a plus sign, 0 otherwise. Default is NULL, i.e. no constraints.
#' @param nbCores number of cores for the parallel computation. Default is 1.
#' @param kp If data_k =NULL, then epsilon is taken equal to kp.
#' @param nb_pts the constant C in DGM for the epsilon_0, the lower bound on the grid for epsilon, taken equal to nb_pts*ln(n)/n. Default is 1 without regressors Xc, 3 with Xc.
#' @param Opt the option of the method to compute the critical values, either "boot" (numerical bootstrap) or "subsampling". Default is numerical bootstrap.
#' @param Bsamp the number of bootstrap/subsampling replications. Default is 1000.
#' @param data_k the number of points for the grid search on epsilon. Default is 30. If NULL, then epsilon is taken fixed equal to kp.
#' @param C0 the upper bound on the grid for epsilon. Default is 0.5.
#' @param list_ex the list of the data that should not be used on the different clusters. This reduces the computational cost. Default is empty list.
#' @param weights_x the sampling weights for the dataset (Xnc,Xc). Default is NULL.
#' @param weights_y  the sampling weights for the dataset (Y,Xc). Default is NULL.
#' @param outside if TRUE indicates that the parallel computing has been launched outside of the function. Default is FALSE.
#' @param trunc equal to 2, for the definition of epsilon.
#' @param boot_par the numerical bootstrap parameter. Default is 0.3 without Xc, 0.35 with Xc.
#' @param meth the method for the choice of epsilon, either "adapt", i.e. adapted to the direction or "min" the minimum over the directions. Default is "adapt".
#' @param modeNA indicates if NA introduced if the interval is empty. Default is FALSE.
#' @param winsor indicates if winsorisation. Default is FALSE.
#' @param version version of the computation of the ratio, "first" is a degraded version but fast; "second" is a correct version but slower. Default is "second".
#' @param alpha for the level of the confidence region. Default is 0.05.
#' @param projections if FALSE compute the identified set along some directions or the confidence regions. Default is FALSE
#'
#' @return a list containing, in order:
#'  - ci : a list with all the information on the confidence intervals
#'
#'     * upper: upper bound of the confidence interval on the radial function S in the specified direction at level alpha, possibly with sign constraints
#'
#'     * lower: lower bound upper bound of the confidence interval on the radial function S, possibly with sign constraints
#'
#'     * unconstr: confidence interval on the radial function S, without sign constraints
#'
#'     * betac_ci: confidence intervals on each coefficients related to the common regressor, possibly with sign constraints
#'
#'     * betac_ci_unc: confidence intervals on each coefficients related to the common regressor without sign constraints
#'
#'     If projection is TRUE:
#'
#'        * support: confidence bound on the support function in each specified direction
#'
#'  - point : a list with all the information on the point estimates
#'
#'     * upper: the upper bounds on betanc, possibly with sign constraints
#'
#'     * lower: the lower bounds on betanc, possibly with sign constraints
#'
#'     * unconstr: bounds on betanc without sign constraints
#'
#'     * betac_pt: bounds on betanc, possibly with sign constraints
#'
#'     * betac_pt_unc: bounds on betanc without sign constraints
#'      If projection ==TRUE:
#'
#'     * support: point estimate of the support function in each specified direction
#'
#'  - epsilon : the values of the selected epsilon(q)
#'
#' @export
#'
#' @examples

DGM_bounds <- function(Ldata, Rdata,values,
                       sam0,
                       refs0,
                       out_var,  nc_var, c_var =NULL,
                       nc_sign  = NULL, c_sign = NULL,
                       nbCores=1,
                       kp=0.5,nb_pts=1,Opt="boot",Bsamp=1000,data_k=30,C0=0.5,
                       list_ex=c(),  weights_x = NULL, weights_y = NULL,outside = FALSE,trunc=2, boot_par=0.3, meth="adapt",
                       modeNA =FALSE, winsor = FALSE , version = "second", alpha=0.05, projections = FALSE){



  ## nb of obs min to use the conditioning.
  limit = 1
  lim = limit
  #### get sizes
  dimXc = length(c_var)
  dimXnc = length(nc_var)

  if(dimXc >=1 & nb_pts <3){
    nb_pts =3
  }
  if(dimXc >=1 &  boot_par<0.35){
    boot_par=0.35
  }
  if(dimXnc>1 & Bsamp>=1000){
    Bsamp=200
  }
  if(dimXnc>1 & data_k>15){
    data_k=15
  }

  ## data_k is the number of points on the grid of the selection rule of epsilon.

  if( dimXc!=0){
    ### dataset 1
    Xc_x = as.matrix(Rdata[,c_var],dim(Rdata)[1],dimXc)
    Xnc = as.matrix(Rdata[,nc_var],dim(Rdata)[1],dimXnc)

    ### dataset 2
    Xc_y = as.matrix(Ldata[,c_var],dim(Ldata)[1],dimXc)
    Y = as.matrix(Ldata[,out_var],dim(Ldata)[1],1)
  }else{

    Xc_x = NULL
    Xnc = as.matrix(Rdata[,nc_var],dim(Rdata)[1],dimXnc)

    ### dataset 2
    Xc_y = NULL
    Y = as.matrix(Ldata[,out_var],dim(Ldata)[1],1)

  }


  if(outside ==FALSE & nbCores>1){
    ### subsampling (B samples) parallel nbCores=10
    sfInit(parallel = TRUE, cpus = nbCores, type = "SOCK")
    sfExportAll( except = list_ex)
    sfLibrary(R.matlab)
    sfLibrary(pracma)
    sfLibrary(Hmisc)
  }

  output  <- vector("list")
  hull_point=NULL
  hull_sharp=NULL

  #################################################################################################################
  ##########  epsilon selection   #################################################################################
  if(!is.null(data_k)){


    if( projections == FALSE){
      if(dimXnc==1){
        sam1 = sam0
      }else{
        sam00 <- eye(dimXnc)
        sam00 <- rbind(-sam00,sam00)
        sam1 = sam00
      }
    }else{
      sam00 <- eye(dimXnc)
      sam00 <- rbind(-sam00,sam00)
      sam1 = sam00
    }

    set.seed(2131)
    kp0 =  select_epsilon(sam1,kp,  Xc_x,Xnc,Xc_y,Y,
                          values,dimXc,dimXnc, nb_pts,Opt,lim =limit ,
                          weights_x,weights_y, refs0,  trunc, boot_par =boot_par, data_k,C0,c_sign,nc_sign,meth=meth,
                          nbCores=nbCores, winsor = winsor, version =version, alpha=alpha )


    set.seed(NULL)

  }else{

    if(meth=="min"){
      kp0 = kp
    }else{
      kp0 = matrix(kp,dim(sam0)[1],1)
    }
  }

  #################################################################################################################
  #################################################################################################################

  #### compute the radial function
  if(projections == FALSE){
    # limit = 100
    sample1 =NULL
    minsel="normal"
    lim =limit

    ### compute point estimate
    type="both"
    set.seed(2131)
    mat_var_out1 <- compute_radial(sample1 =sample1,Xc_x,Xnc,Xc_y,Y,
                                   values,dimXc,dimXnc,
                                   nb_pts, sam0, kp0,data_k,Opt,lim =limit,
                                   weights_x,weights_y,
                                   c_sign, nc_sign,refs0,type="both",meth=meth,trunc=trunc,
                                   boot_par=boot_par,winsor= winsor, version = version )




    set.seed(NULL)

    hull_point <-    mat_var_out1
    if(modeNA ==TRUE){
      no_inter =  hull_point [["upper"]] <  hull_point[["lower"]]
      hull_point [["upper"]][no_inter] <- NA
      hull_point [["lower"]][no_inter] <- NA

      if(dimXnc>1 & length(c_sign)>0){
        if(sum(abs(c_sign))>0){
          no_inter = apply(matrix((hull_point [["tests"]][,-c(1)]< 10^(-5)), dim(hull_point [["tests"]])[1],  dim(hull_point [["tests"]])[2]-1),1,prod)*TRUE
          hull_point [["upper"]][!no_inter] <- NA
          hull_point [["lower"]][!no_inter] <- NA
        }
      }
    }


    ##### adjust if non relevant sign constraints ################################
    if(length(c_sign)>0){
      if(sum(abs(c_sign))>0){
        effective_c_sign = c_sign
        for(j in 1:length(c_sign)){
          tests=matrix(hull_point$tests[,-c(1)],dim(sam0)[1],length(c_sign))
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


    ##### compute point estimate of betac if common regressors Xc ###################################################################################""
    if(!is.null(values)){
      beta1_pt <- compute_bnds_betac(sample1 =NULL,Opt1 = Opt, info0 = mat_var_out1, values,  c_sign0 = c_sign, nc_sign0=nc_sign ,  refs0, c_var, nc_var,sam0)
      mat_beta1_l = matrix(0,Bsamp,length(refs0))
      mat_beta1_u = matrix(0,Bsamp,length(refs0))
      hull_point[["betac_pt_unc"]] <-  beta1_pt
    }



    #### replications numerical bootstrap or subsampling ############################################################################################
    if(Bsamp >0){
      if(nbCores>1){
        res0 <- sfLapply(1:Bsamp, compute_radial, Xc_x,Xnc,Xc_y,Y,values,dimXc,dimXnc,nb_pts,
                         sam0, kp0,data_k,Opt,lim =limit ,
                         weights_x,   weights_y ,c_sign , nc_sign,refs0,type="both",meth=meth,
                         trunc=trunc,boot_par=boot_par, winsor= winsor , version = version,
                         ratio_ref=  hull_point[["ratio_ref"]], ratio_ref_m=  hull_point[["ratio_ref_m"]],
                         Dmat_Yk_ref= hull_point[["DYk"]], Dmat_Xk_ref=hull_point[["DXk"]], grid_I_ref=hull_point[["grid_ref"]] )
      }else{

        res0 <- lapply(1:Bsamp, compute_radial, Xc_x,Xnc,Xc_y,Y,
                       values,dimXc,dimXnc,nb_pts,sam0, kp0,data_k ,Opt,lim =limit ,
                       weights_x,   weights_y,c_sign , nc_sign,refs0,type="both",meth=meth,
                       trunc=trunc,boot_par=boot_par,winsor= winsor , version = version,
                       ratio_ref=  hull_point[["ratio_ref"]], ratio_ref_m=  hull_point[["ratio_ref_m"]],
                       Dmat_Yk_ref= hull_point[["DYk"]], Dmat_Xk_ref=hull_point[["DXk"]], grid_I_ref=hull_point[["grid_ref"]])
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
      alpha1= alpha
      q95s <-  function(x){quantile(x,1-alpha1,na.rm=T)}
      q05s <-  function(x){quantile(x,alpha1,na.rm=T)}

      n_x = dim(Xnc)[1]
      n_y = dim(Y)[1]

      if(Opt=="subsampling"){
        T_xy = (n_y/(n_x+n_y))*n_x
        bs = sampling_rule(T_xy)
        bs0 = sqrt(bs/T_xy)
        mat_var_out =   hull_point[["upper"]] - apply((mat_varb -  matrix(1,dim(mat_varb)[1],1)%*% hull_point[["upper"]])*bs0,2,q05s)
        mat_varb_out_unc  =   hull_point[["unconstr"]] - apply((mat_varb_unc -   matrix(1,dim(mat_varb_unc)[1],1)%*% hull_point[["unconstr"]])*bs0  ,2,q05)
      }else{ ### numerical bootstrap
        T_xy = (n_y/(n_x+n_y))*n_x
        es =  T_xy^(-boot_par)
        term_unc = matrix(es^(-1)/sqrt(T_xy),dim(mat_varb)[1],dim(sam0)[1])
        term =   term_unc
        mat_var_out = hull_point[["upper"]] -   apply( (mat_varb -  matrix(1,dim(mat_varb)[1],1)%*% hull_point[["upper"]])*term ,2,q05s)
        mat_varb_out_unc  =  hull_point[["unconstr"]]- apply( (mat_varb_unc -   matrix(1,dim(mat_varb_unc)[1],1)%*% hull_point[["unconstr"]])*term_unc ,2,q05)
      }

      if(Opt=="subsampling"){
        mat_var_in =   hull_point[["lower"]] - apply(mat_varb0   -  matrix(1,dim(mat_varb0)[1],1)%*% hull_point[["lower"]],2,q95s)*bs0
      }else{
        term_low = matrix(es^(-1)/sqrt(T_xy),dim(mat_varb)[1],dim(sam0)[1])
        mat_var_in =    hull_point[["lower"]] - apply( (mat_varb0   -  matrix(1,dim(mat_varb0)[1],1)%*% hull_point[["lower"]])*term_low ,2,q95s)
      }

      ##
      hull_sharp[["upper"]] <-  mat_var_out
      hull_sharp[["lower"]] <-  mat_var_in
      hull_sharp[["unconstr"]] <-  mat_varb_out_unc

      ## stock replications
      hull_sharp[["upper_repli"]] <-  mat_varb
      hull_sharp[["lower_repli"]] <-  mat_varb0
      hull_sharp[["unconstr_repli"]] <-  mat_varb_unc

      #### compute the replications for beta_c ###################################################################################
      if(Opt=="subsampling"){

        if(!is.null(values)){

          ### with sign constraints ######################################################################
          if(nbCores>1){
            res0b <- sfLapply(1:Bsamp,compute_bnds_betac,Opt1 = Opt, info0 = res0, values,  c_sign0 = c_sign, nc_sign0=nc_sign ,
                              refs0, c_var, nc_var,sam0, info1=mat_var_out1, T_xy,es)
          }else{
            res0b <- lapply(1:Bsamp,compute_bnds_betac, Opt1 = Opt,info0 = res0, values,  c_sign0 = c_sign, nc_sign0=nc_sign ,
                            refs0, c_var, nc_var,sam0, info1=mat_var_out1,T_xy,es)
          }
          for(b in 1:Bsamp){
            mat_beta1_l[b,] = res0b[[b]][,1]
            mat_beta1_u[b,] = res0b[[b]][,2]
          }
          bs = floor(sampling_rule(n_x*(n_y/(n_x+n_y))))
          term = sqrt(bs/T_xy)

          beta1_l_ci =  beta1_pt[,1] +   apply( (mat_beta1_l -  matrix(1,dim(mat_beta1_l)[1],1)%*% beta1_pt[,1])* term,2,q05)
          beta1_u_ci   = beta1_pt[,2] +   apply( (mat_beta1_u -  matrix(1,dim(mat_beta1_l)[1],1)%*% beta1_pt[,2])*  term,2,q95)
          beta1_ci <- cbind(beta1_l_ci,  beta1_u_ci )
          hull_sharp[["betac_ci"]] <-  beta1_ci
          hull_point[["betac_pt"]] <-  beta1_pt


          ### without sign constraints ######################################################################"
          beta1_pt <- compute_bnds_betac(sample1 =NULL,Opt1 = Opt, info0 = mat_var_out1, values,  c_sign0 = c_sign, nc_sign0=nc_sign , refs0, c_var, nc_var,  sam0, info1=NULL , T_xy=NULL,es=NULL,constr=FALSE)
          mat_beta1_l = matrix(0,Bsamp,length(refs0))
          mat_beta1_u = matrix(0,Bsamp,length(refs0))

          if(nbCores>1){
            res0b <- sfLapply(1:Bsamp,compute_bnds_betac,Opt1 = Opt, info0 = res0, values,  c_sign0 = c_sign, nc_sign0=nc_sign ,  refs0, c_var,
                              nc_var, sam0, info1=mat_var_out1 ,  T_xy,es, constr=FALSE)
          }else{
            res0b <- lapply(1:Bsamp,compute_bnds_betac, Opt1 = Opt,info = res0, values,  c_sign0 = c_sign, nc_sign0=nc_sign ,  refs0, c_var,
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
        }

      }else{ ### numerical bootstrap
        ### with sign constraints ######################################################################
        if(!is.null(values)){
          if(nbCores>1){
            res0b <- sfLapply(1:Bsamp,compute_bnds_betac, Opt1 = Opt,info0 = res0, values,
                              c_sign0 = c_sign, nc_sign0=nc_sign ,  refs0, c_var, nc_var,sam0,
                              info1=mat_var_out1, T_xy,es)
          }else{
            res0b <- lapply(1:Bsamp,compute_bnds_betac, Opt1 = Opt,info0 = res0, values,
                            c_sign0 = c_sign, nc_sign0=nc_sign ,  refs0, c_var, nc_var,sam0,
                            info1=mat_var_out1,T_xy,es)
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


          ### without sign constraints ######################################################################
          beta1_pt <- compute_bnds_betac(sample1 =NULL, Opt1 = Opt,info0 = mat_var_out1, values,  c_sign0 = c_sign, nc_sign0=nc_sign , refs0, c_var, nc_var,  sam0, info1=NULL , T_xy=NULL,es=NULL,constr=FALSE)
          mat_beta1_l = matrix(0,Bsamp,length(refs0))
          mat_beta1_u = matrix(0,Bsamp,length(refs0))

          if(nbCores>1){
            res0b <- sfLapply(1:Bsamp,compute_bnds_betac,Opt1 = Opt, info0 = res0, values,  c_sign0 = c_sign, nc_sign0=nc_sign ,  refs0, c_var,
                              nc_var, sam0, info1=mat_var_out1 ,  T_xy,es, constr=FALSE)
          }else{
            res0b <- lapply(1:Bsamp,compute_bnds_betac, Opt1 = Opt,info = res0, values,  c_sign0 = c_sign, nc_sign0=nc_sign ,  refs0, c_var,
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
        }
      }

      if(modeNA ==TRUE){
        no_inter =   hull_sharp[["upper"]] <  hull_sharp[["lower"]]
        hull_sharp[["upper"]][no_inter] <- NA
        hull_sharp[["lower"]][no_inter] <- NA
      }
    }else{
      hull_sharp <- matrix(NA,1,1)

    }

    if(outside ==FALSE & nbCores>1){
      sfStop()
    }

    ###########################################################################################################################################################
    ########################################### Compute projections (Support function) ########################################################################
    ###########################################################################################################################################################
  }else{


    minsel="normal"
    mat_var_out1 <- compute_support(sample1 =NULL,Xc_x,Xnc,Xc_y,Y,
                                    values,dimXc,dimXnc,
                                    nb_pts, sam0, kp0, data_k, Opt,lim =limit,
                                    weights_x,weights_y,
                                    c_sign, nc_sign,refs0,type="both",meth=meth,trunc=trunc,boot_par=boot_par, bc=TRUE,
                                    winsor= winsor, version = version)

    hull_point[["support"]] <-      mat_var_out1

    if(Bsamp!=0){
      # mat_var_out1
      Bsamp1 = Bsamp
      # start_time <- Sys.time()
      if(nbCores>1){
        res0 <- sfLapply(1:Bsamp1, compute_support,Xc_x,Xnc,Xc_y,Y,values,dimXc,dimXnc,nb_pts,sam0, kp0, data_k,Opt,lim =limit ,
                         weights_x,   weights_y,
                         c_sign , nc_sign,refs0,type="both",meth=meth,
                         trunc=trunc,boot_par=boot_par,bc=TRUE,winsor= winsor, version = version)
      }else{

        res0 <- lapply(1:Bsamp1, compute_support, Xc_x,Xnc,Xc_y,Y,values,dimXc,dimXnc,nb_pts,sam0, kp0, data_k,Opt,lim =limit ,
                       weights_x,   weights_y,
                       c_sign , nc_sign,refs0,type="both",meth=meth,
                       trunc=trunc,boot_par=boot_par,bc=TRUE,winsor= winsor, version = version)
      }
      mat_varb = matrix(0,Bsamp1,dim(sam0)[1])
      for(b in 1:Bsamp1){
        mat_varb[b,] = res0[[b]][,3]
      }

      q95 <-  function(x){quantile(x,1-alpha,na.rm=T)}
      q05 <-  function(x){quantile(x,alpha,na.rm=T)}

      n_x = dim(Xnc)[1]
      n_y = dim(Y)[1]
      if(Opt=="boot"){ ### numerical bootstrap
        T_xy=n_x*(n_y/(n_x+n_y))
        es =  T_xy^(-boot_par)
        mat_var_out = mat_var_out1[,3] -   apply(mat_varb -  matrix(1,dim(mat_varb)[1],1)%*% mat_var_out1[,3],2,q05)*es^(-1)/sqrt(T_xy)
      }else{ #### subsampling
        T_xy=n_x*(n_y/(n_x+n_y))
        bs = sampling_rule(T_xy)
        bs0 = sqrt(bs/T_xy)
        mat_var_out = mat_var_out1[,3] -   apply(mat_varb -  matrix(1,dim(mat_varb)[1],1)%*% mat_var_out1[,3],2,q05)* bs0
      }
      hull_sharp[["support"]] <-  cbind(sam0,mat_var_out)
    }else{
      hull_sharp[["support"]] <-  matrix(NA,1,1)
    }

    if(outside ==FALSE & nbCores>1){
      sfStop()
    }

  } ###########################

  output[["ci"]] <- hull_sharp
  output[["point"]] <- hull_point
  output[["epsilon"]] <- kp0
  return(output)
}
