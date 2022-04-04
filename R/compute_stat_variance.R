#' Function to compute the Variance bounds on the noncommon regressor Xnc
#'
#' @param sample1 if NULL compute the point estimate, if a natural number then evaluate a bootstrap or subsampling replication.
#' @param X1_x the common regressor on the dataset  (Xnc,Xc). Default is NULL.
#' @param X2 the noncommon regressor on the dataset  (Xnc,Xc). No default.
#' @param X1_y the common regressor on the dataset  (Y,Xc). Default is NULL.
#' @param Y the outcome variable. No default.
#' @param values the different unique points of support of the common regressor Xc.
#' @param refs0 indicating the positions in the vector values corresponding to the components of betac.
#' @param dimX1 the dimension of the common regressors Xc.
#' @param dimX2 the dimension of the noncommon regressors Xnc.
#' @param nb_pts the constant C in DGM for the epsilon_0, the lower bound on the grid for epsilon, taken equal to nb_pts*ln(n)/n. Default is 1 without regressors Xc, 3 with Xc.
#' @param sam0 the directions q to compute the variance bounds on the radial function.
#' @param kp If data_k =NULL, then epsilon is taken equal to kp.
#' @param Opt the option of the method to compute the critical values, either "boot" (numerical bootstrap) or "subsampling". Default is numerical bootstrap.
#' @param lim the limit number of observations under which we do no compute the conditional variance.
#' @param weights_x the sampling weights for the dataset (Xnc,Xc).
#' @param weights_y  the sampling weights for the dataset (Y,Xc).
#' @param boot_par the numerical bootstrap parameter. Default is 0.3 without Xc, 0.35 with Xc.
#' @param c_sign sign restrictions on the commonly observed regressors: -1 for a minus sign, 1 for a plus sign, 0 otherwise. Default is NULL, i.e. no constraints.
#' @param nc_sign sign restrictions on the non-commonly observed regressors Xnc: -1 for a minus sign, 1 for a plus sign, 0 otherwise. Default is NULL, i.e. no constraints.
#' @param ratio_ref the variances ratio of the point estimate for the computation of the numerical boostrap replications. Default is NULL.
#' @param Dmat_Yk_ref the vector of means of Y| Xc on the initial dataset.
#' @param Dmat_Xk_ref  the vector of means of Xnc| Xc on the initial dataset.
#'
#' @return a list containing:
#'
#'  - upper: the upper bound in the specified directions, possibly with sign constraints
#'
#'  - lower: the lower bound in the specified directions, possibly with sign constraints
#'
#'  - unconstr: the bounds without sign constraints in the specified directions
#'
#'  - Ykmean: the means of Y|Xc for the considered sample
#'
#'  - Xkmean: the means of Xnc|Xc for the considered sample
#'
#'  - DYk: the difference of means of Y|Xc =k -  Y|Xc =0 for the considered sample
#'
#'  - DXk: the difference of means of Xnc|Xc =k -  Xnc|Xc =0 for the considered sample
#'
#'  - tests: the pvalues of the tests H0 : DXk =0
#'
#'  - ratio_ref: the ratio R in the radial function computed for the initial sample
#'
#'
#' @export
#'
#' @examples
compute_stat_variance <- function(sample1 = NULL,X1_x,X2,X1_y,Y,values,
                            refs0,dimX1,dimX2,
                            nb_pts,sam0,kp,Opt="boot",
                            lim = 1,
                            weights_x = NULL,weights_y = NULL,
                            boot_par=0.35,c_sign= NULL,nc_sign= NULL,
                            ratio_ref=NULL,
                            Dmat_Yk_ref=NULL, Dmat_Xk_ref=NULL){

  if( dimX1==0){
    X1_x = NULL
    # X2 =X2

    ### dataset 2
    X1_y = NULL
    # Y = Y

  }

  n_x = dim(X2)[1]
  n_y = dim(Y)[1]
  n_xy= min(n_x,n_y)
  T_xy= n_xy

  if(is.null(weights_x)){
    weights_x= rep(1/dim(X2)[1],dim(X2)[1])
  }
  if(is.null(weights_y)){
    weights_y= rep(1/length(Y),length(Y))
  }

  if(!is.null(sample1)){
    if(Opt == "subsampling"){
      ## subsampling

      bsx = floor(0.3*dim(X1_x)[1])
      bb = sample(1:dim(X1_x)[1],bsx, replace=FALSE)
      X1_xb = matrix(X1_x[bb,],bsx,dimX1)
      X2b = matrix(X2[bb,],bsx,dimX2)
      weights_x =  matrix(weights_x[bb],bsx,1)
      weights_x = weights_x/sum(weights_x)

      bsy = floor(0.3*dim(Y)[1])
      bby = sample(1:dim(X1_y)[1],bsy, replace=FALSE)
      X1_yb = matrix(X1_y[bb,],bsy,dimX1)
      Yb = matrix(Y[bby,],bsy,1)
      weights_y =  matrix(weights_y[bby],bsy,1)
      weights_y = weights_y/sum(weights_y)

    }else{
      ## bootstrap
      es =  T_xy^(-boot_par)
      drawsx = rmultinom(1, n_x, rep(1/n_x,n_x))
      if(!is.null(X1_x)){
        X1_xb =matrix(X1_x,n_x,dimX1)
      }
      X2b = matrix(X2,n_x,dimX2)
      weights_x =  matrix(weights_x*drawsx,n_x,1)
      weights_x = weights_x/sum(weights_x)
      drawsy = rmultinom(1, n_y, rep(1/n_y,n_y))
      if(!is.null(X1_y)){
        X1_yb = matrix(X1_y,n_y,dimX1)
      }
      Yb =matrix(Y,n_y,1)
      weights_y =  matrix(weights_y* drawsy,n_y,1)
      weights_y = weights_y/sum(weights_y)

    }

  }else{
    ## point estimate
    X1_xb =X1_x
    X2b = X2
    X1_yb =   X1_y
    Yb = Y

  }


  if(!is.null(values)){### if there is a discrete common regressor X1


    n_x_all = n_x
    n_y_all = n_y
    T_xy_all=   min( n_x_all,  n_y_all)
    es =  T_xy_all^(-boot_par)

    ind = NULL
    mat_var= matrix(0,dim(values)[1],dim(sam0)[1])
    mat_var_unc= matrix(0,dim(values)[1],dim(sam0)[1])
    mat_var_low= matrix(0,dim(values)[1],dim(sam0)[1])
    mat_Yk= matrix(NA,dim(values)[1],1)
    mat_Xk= matrix(NA,dim(values)[1],dimX2)
    Dmat_Yk= matrix(NA,dim(values)[1],1)
    Dmat_Xk= vector("list")

    inds=0
    if(dimX1==1){
      val0 = values[1,]
      sel0_x =(X1_xb==val0)
      sel0_y =(X1_yb==val0)
    }else{
      val0 = t(as.matrix(values[1,]))
      sel0_x = matrix(1,dim(X1_xb)[1],1)
      sel0_y = matrix(1,dim(X1_yb)[1],1)
      for(ddd in 1:dimX1){
        sel0_x =  sel0_x & (X1_xb[,ddd]==val0[ddd])
        sel0_y =  sel0_y & (X1_yb[,ddd]==val0[ddd])
      }
      sel0_x = matrix( sel0_x,dim(X1_xb)[1],1)
      sel0_y = matrix( sel0_y,dim(X1_yb)[1],1)
    }

    weights_xp0 =  matrix(weights_x[sel0_x],sum(sel0_x),1)
    weights_xp0 = weights_xp0/sum(weights_xp0)
    weights_yp0 =  matrix(weights_y[sel0_y],sum(sel0_y),1)
    weights_yp0 = weights_yp0/sum(weights_yp0)
    Xp0 = matrix(X2b[sel0_x,],sum(sel0_x),dimX2)
    X_0 =((sam0%*%t(Xp0))*(matrix(1,dim(sam0)[1],1)%*%matrix(weights_xp0,1,sum(sel0_x))))%*%matrix(1,dim(Xp0)[1],1)
    Yp0 = Yb[sel0_y]
    Y0 = sum(Yp0*weights_yp0)
    tests = matrix(NA,dim(sam0)[1],dim(values)[1])
    T_n = matrix(NA,1,dim(values)[1])

    if(is.null(sample1)){
      ratio_ref= vector("list")
    }

    for(k in 1:dim(values)[1]){


      if(dimX1==1){
        val = values[k,]
        sel_x = (X1_xb==val)
        sel_y = (X1_yb==val)
      }else{
        val = t(as.matrix(values[k,]))
        sel_x = matrix(1,dim(X1_xb)[1],1)
        sel_y = matrix(1,dim(X1_yb)[1],1)
        for(ddd in 1:dimX1){
          sel_x =  sel_x & (X1_xb[,ddd]==val[ddd])
          sel_y =  sel_y & (X1_yb[,ddd]==val[ddd])
        }
        sel_x = matrix( sel_x,dim(X1_xb)[1],1)
        sel_y = matrix( sel_y,dim(X1_yb)[1],1)

      }


      if(sum(sel_x)> lim & sum(sel_y) >  lim){
        Xp= matrix(X2b[sel_x,],sum(sel_x),dimX2);
        Yp = Yb[sel_y]
        weights_yp =   weights_y[sel_y]
        if(sum(weights_yp!=0)<=1){
          weights_yp=    weights_yp + 1/length( weights_yp)
          weights_yp = weights_yp /sum( weights_yp )
        }
        weights_yp  =  weights_yp /sum( weights_yp )
        weights_xp =   weights_x[sel_x]

        if(sum(weights_xp!=0)<=1){
          weights_xp=    weights_xp + 1/length( weights_xp)
          weights_xp = weights_xp /sum( weights_xp )
        }
        weights_xp  =  weights_xp /sum( weights_xp )
        n_x = sum(sel_x)
        n_y = sum(sel_y)
        n_xy = min(n_x,n_y)
        T_xy=   n_xy
        T_n[1,k] =  T_xy
        sam1 = cbind(1:dim(sam0)[1],sam0)

        if(is.null(sample1)){
          bsharp_beta2 <- na.omit(t(apply(sam1,1,compute_ratio_variance,Xp=Xp,Yp=Yp,dimX2=dimX2,
                                          weights_xp ,weights_yp,T_xy_all, es,boot_par, Opt )))
          ratio_ref[[k]] <-bsharp_beta2
        }else{
          ratio_ref_k = ratio_ref[[k]]
          bsharp_beta2 <- na.omit(t(apply(sam1,1,compute_ratio_variance,Xp=Xp,Yp=Yp,dimX2=dimX2,
                                          weights_xp ,weights_yp,T_xy_all, es,boot_par,Opt, ratio_ref=ratio_ref_k)))

        }

        mat_Yk[k,1] <- sum(Yp*weights_yp)
        mat_Xk[k,] <- sum(Xp*weights_xp)

        mat_var[k,] <- bsharp_beta2
        mat_var_unc[k,] <-mat_var[k,]
        mat_var_low[k,] <- - rev(mat_var[k,])
        bsharp_beta2_m  = rev(bsharp_beta2)

        EU = FALSE
        EL=  FALSE

        Ybarre = sum(Yp*weights_yp);
        X_k = ((sam0%*%t(Xp))*(matrix(1,dim(sam0)[1],1)%*%weights_xp))%*%matrix(1,dim(Xp)[1],1)
        dXk = X_k- X_0
        dYk = Ybarre- Y0

        for( dir in 1:dim(sam0)[1]){

          tt = t.test((sam0[dir,]%*%t(Xp)), y = (sam0[dir,]%*%t(Xp0)),alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
          tests[dir,k] = tt$p.value #< 0.001

        }


        Dmat_Yk[k,1] <- dYk
        Dmat_Xk[[k]] <- dXk

        mean_ratio =  dYk/dXk

        off = TRUE
        if(!is.null(sample1) & Opt=="boot"){
          mean_ratio_ref =  Dmat_Yk_ref[k,1]/Dmat_Xk_ref[[k]]
          es = T_xy_all^(-boot_par)

          mean_ratio =  mean_ratio_ref + ( mean_ratio -mean_ratio_ref )*es*sqrt(T_xy_all)
          if(sum(is.na(mean_ratio) >0)){
            off=FALSE
          }else{
            for(l in 1:length(mean_ratio)){
              if(sign(mean_ratio[l])!=sign(mean_ratio_ref[l])){off=FALSE}
            }
          }

        }


        if(is.null(sample1)){
          ind = c(ind,k)
        }else{
          if(sum(mat_var_unc[k,]==0)==0){
            ## add signal if not, raise
            ind = c(ind,k)
          }
        }


        indk = grep(1,values[k,])
        if(sum(k==refs0)>0 & length(indk) ==1 & k!=1 & !is.null(c_sign) & off ){
          if(c_sign[indk]!=0){
            EL = (c_sign[indk]>0 & dXk<0  ) |  (c_sign[indk]<0 &dXk>0  )
            EU = (c_sign[indk]>0 & dXk>0  ) |  (c_sign[indk]<0 & dXk<0  )
            for(jj in 1:dim(sam0)[1]){
              if(EU[jj] &  (tests[jj,k] < 10^(-2))){
                mat_var[k,jj] <-min( bsharp_beta2[jj] ,    mean_ratio[jj] )
              }else{
                mat_var[k,jj] <- bsharp_beta2[jj]
              }
               if(EL[jj] &  (tests[jj,k] <  10^(-2))){
                mat_var_low[k,jj] <- max( - bsharp_beta2_m[jj] , mean_ratio[jj])
              }
            }


          }else{
            mat_var[k,] <- bsharp_beta2
          }


        }else{
          mat_var[k,] <- bsharp_beta2

        }
      }else{
        mat_var[k,] <- bsharp_beta2

      }
    }
    # names[138]
    mat_var_unc1<- mat_var_unc[ind ,]
    mat_var_unc1 <- matrix(mat_var_unc1,length(ind),dim(sam0)[1])
    mat_var1 <-   mat_var[ind,]
    mat_var1 <- matrix(mat_var1,length(ind),dim(sam0)[1])
    mat_var_low1 <-   mat_var_low[ind,]
    mat_var_low1 <- matrix(mat_var_low1,length(ind),dim(sam0)[1])

    if(!is.null(values)){
      mat_var <- apply(  mat_var1,2,min)
      mat_var_low <- apply(  mat_var_low1,2,max)
      mat_var_unc <- apply(  mat_var_unc1,2,min)
    }

    if(!is.null(nc_sign)){
      ee = eye(dimX2)
      for(jj in 1:dimX2){
        if(nc_sign[jj]!=0){
          cond = (matrix(ee[jj,]*nc_sign[jj],1,dimX2)%*%t(sam0))<0
          mat_var[cond]<- pmin(mat_var[cond],0)
          mat_var_low[!cond]<- pmax(mat_var_low[!cond],0)
        }
      }
    }

  }else{

    mat_var= matrix(0,1,dim(sam0)[1])
    mat_var_unc = matrix(0,1,dim(sam0)[1])
    Xp= X2b
    Yp = Yb
    weights_yp =   weights_y/sum( weights_y)
    weights_xp =   weights_x/sum( weights_x )

    n_x = dim(  Xp)[1]
    n_y = dim(Yp)[1]
    n_xy = min(n_x,n_y)
    T_xy= n_xy
    es =    T_xy^(-boot_par)

    sam1 = cbind(1:dim(sam0)[1],sam0)


    if(is.null(sample1)){
      bsharp_beta2 <- na.omit(t(apply(sam1,1,compute_ratio_variance,Xp=Xp,Yp=Yp,dimX2=dimX2,
                                      weights_xp ,weights_yp,T_xy, es ,boot_par, Opt)))
      ratio_ref =  bsharp_beta2
    }else{
      bsharp_beta2 <- na.omit(t(apply(sam1,1,compute_ratio_variance,Xp=Xp,Yp=Yp,dimX2=dimX2,
                                      weights_xp,weights_yp, T_xy, es ,boot_par, Opt,ratio_ref = ratio_ref )))
    }

    mat_var<- bsharp_beta2
    mat_var_low <- -rev( mat_var)

  }


  output <- vector("list")
  output[["upper"]] <- mat_var
  output[["unconstr"]] <-mat_var_unc
  output[["lower"]] <- mat_var_low
  if(!is.null(values)){
    output[["Ykmean"]] <- mat_Yk
    output[["Xkmean"]] <- mat_Xk
    output[["DYk"]] <- Dmat_Yk
    output[["DXk"]] <- Dmat_Xk
    output[["tests"]] <- tests
    output[["T_n"]] <- T_n
  }
  if(is.null(sample1)){
    output[["ratio_ref"]] <-ratio_ref
  }

  return(output)

}
