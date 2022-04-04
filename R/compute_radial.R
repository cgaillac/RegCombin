#' Function to compute the DGM bounds on the noncommon regressor Xnc
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
#' @param Opt  the option of the method to compute the critical values, either "boot" (numerical bootstrap) or "subsampling". Default is numerical bootstrap.
#' @param lim the limit number of observations under which we do no compute the conditional variance.
#' @param weights_x the sampling weights for the dataset (Xnc,Xc).
#' @param weights_y  the sampling weights for the dataset (Y,Xc).
#' @param c_sign sign restrictions on the commonly observed regressors: -1 for a minus sign, 1 for a plus sign, 0 otherwise. Default is NULL, i.e. no contraints.
#' @param nc_sign sign restrictions on the non-commonly observed regressors Xnc: -1 for a minus sign, 1 for a plus sign, 0 otherwise. Default is NULL, i.e. no contraints.
#' @param refs0  indicating the positions in the vector values corresponding to the components of betac.
#' @param type Equal to "both".
#' @param meth the method for the choice of epsilon, either "adapt", i.e. adapted to the direction or "min" the minimum over the directions. Default is "adapt".
#' @param trunc equal to 2, for the definition of epsilon.
#' @param boot_par the numerical bootstrap parameter. Default is 0.3 without Xc, 0.35 with Xc.
#' @param winsor indicates if winsorisation. Default is FALSE.
#' @param version version of the computation of the ratio, "first" is a degraded version but fast; "second" is a correct version but slower. Default is "second".
#' @param ratio_ref  the point estimate of the ratio R in the radial function.
#' @param ratio_ref_m  the point estimate of the ratio R in the radial function on the opposite direction.
#' @param Dmat_Yk_ref the vector of means of Y| Xc on the initial dataset.
#' @param Dmat_Xk_ref  the vector of means of Xnc| Xc on the initial dataset.
#' @param es_all the numerical bootstrap parameter delta_n
#' @param T_xy_all the apparent sample size the taking into account the difference in the two datasets.
#' @param grid_I_ref the grid of alpha on which we evaluate the ratio R to compute the point estimate of the radial function.
#'
#' @return  a list contaning:
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
#
compute_radial <- function(sample1 = NULL,Xc_x,Xnc,Xc_y,Y,values,
                      dimXc,dimXnc,nb_pts,
                      sam0,kp0,data_k = NULL, Opt="boot",
                      lim = 1,weights_x = NULL,weights_y = NULL,
                      c_sign = NULL, nc_sign= NULL,
                      refs0=NULL,type="both",meth="adapt",
                      trunc=2,boot_par= 0.35, winsor= winsor , version = "first",
                      ratio_ref=NULL, ratio_ref_m=NULL, Dmat_Yk_ref=NULL,
                      Dmat_Xk_ref=NULL , grid_I_ref=NULL, es_all = NULL,  T_xy_all = NULL){

  if(is.null(kp0)){
    learn_eps = TRUE
  }else{
    learn_eps = FALSE
  }

  mat_var_low= NULL
  if(is.null(weights_x)){
    weights_x= rep(1/dim(Xnc)[1],dim(Xnc)[1])
  }
  if(is.null(weights_y)){
    weights_y= rep(1/length(Y),length(Y))
  }

  # save original weights
  weights_xs <-  weights_x
  weights_ys <-  weights_y


  if(!is.null(sample1)){

    if(Opt == "subsampling"){
      ## subsampling
      n_x = dim(Xnc)[1]
      n_y = dim(Y)[1]
      bs = floor(sampling_rule(n_x*(n_y/(n_x+n_y))))

      n_x = dim(Xnc)[1]
      bb = sample(1:n_x,n_x-bs, replace=FALSE)
      if(!is.null(Xc_x)){
        Xc_xb = matrix(Xc_x,n_x,dimXc)
      }
      Xncb = matrix(Xnc,n_x,dimXnc)
      weights_x =  matrix(weights_x,n_x,1)
      weights_x[bb] <-0
      weights_x = weights_x/sum(weights_x)

      n_y = dim(Y)[1]
      bby = sample(1:n_y,n_y-bs, replace=FALSE)
      if(!is.null(Xc_y)){
        Xc_yb = matrix(Xc_y,n_y,dimXc)
      }
      Yb = matrix(Y,n_y,1)
      weights_y =  matrix(weights_y,n_y,1)
      weights_y[bby] <-0
      weights_y = weights_y/sum(weights_y)

    }else{

      n_x = dim(Xnc)[1]
      n_y = dim(Y)[1]

      drawsx = rmultinom(1, n_x, rep(1/n_x,n_x))
      if(!is.null(Xc_x)){
        Xc_xb =matrix(Xc_x,n_x,dimXc)
      }
      Xncb = matrix(Xnc,n_x,dimXnc)
      weights_x =  matrix(weights_x*drawsx,n_x,1)
      weights_x = weights_x/sum(weights_x)
      drawsy = rmultinom(1, n_y, rep(1/n_y,n_y))
      if(!is.null(Xc_y)){
        Xc_yb = matrix(Xc_y,n_y,dimXc)
      }
      Yb =matrix(Y,n_y,1)
      weights_y =  matrix(weights_y* drawsy,n_y,1)
      weights_y = weights_y/sum(weights_y)
    }
  }else{
    ## point estimate

    Xc_xb =Xc_x
    Xncb = Xnc
    Xc_yb = Xc_y
    Yb = Y

  }

  if(!is.null(values)){### if there is a discrete common regressor Xc
    # mat_var = NULL
    n_x_all = dim(Xc_xb)[1]
    n_y_all = dim(Xc_yb)[1]
    T_xy_all = min(n_x_all,n_y_all)
    mat_Yk= matrix(NA,dim(values)[1],1)
    mat_Xk= matrix(NA,dim(values)[1],dimXnc)
    Dmat_Yk= matrix(NA,dim(values)[1],1)
    Dmat_Xk= vector("list")


    mat_var_unc= matrix(0,dim(values)[1],dim(sam0)[1])
    mat_var= matrix(0,dim(values)[1],dim(sam0)[1])
    mat_var_low= matrix(0,dim(values)[1],dim(sam0)[1])
    ind = NULL
    inds=0
    if(dimXc==1){
      val0 =values[1,]
      sel0_x =  (Xc_xb==val0)
      sel0_y =  (Xc_yb==val0)
    }else{
      val0 = t(as.matrix(values[1,]))
      sel0_x = matrix(1,dim(Xc_xb)[1],1)
      sel0_y = matrix(1,dim(Xc_yb)[1],1)
      for(ddd in 1:dimXc){
        sel0_x =  sel0_x & (Xc_xb[,ddd]==val0[ddd])
        sel0_y =  sel0_y & (Xc_yb[,ddd]==val0[ddd])
      }
      sel0_x = matrix( sel0_x,dim(Xc_xb)[1],1)
      sel0_y = matrix( sel0_y,dim(Xc_yb)[1],1)
    }

    weights_xp0 =  matrix(weights_x[sel0_x],sum(sel0_x),1)
    weights_xp0 = weights_xp0/sum(weights_xp0)
    weights_yp0 =  matrix(weights_y[sel0_y],sum(sel0_y),1)
    weights_yp0 = weights_yp0/sum(weights_yp0)
    Xp0 = matrix(Xncb[sel0_x,],sum(sel0_x),dimXnc)
    if(!is.null(kp0)){
      X_0 =((sam0%*%t(Xp0))*(matrix(1,dim(sam0)[1],1)%*%matrix(weights_xp0,1,sum(sel0_x))))%*%matrix(1,dim(Xp0)[1],1)
    }
    Yp0 = Yb[sel0_y]
    Y0 = sum(Yp0*weights_yp0)
    grid_I = vector("list")
    #### vector of matrices of point estimate ratios.
    if(is.null(sample1)){
      ratio_ref= vector("list")
      ratio_ref_m= vector("list")
    }
    tests = matrix(NA,dim(sam0)[1],dim(values)[1])
    T_n = matrix(NA,1,dim(values)[1])



    for(k in 1:dim(values)[1]){ ## for all the points of support of Xc


      if(dimXc==1){
        val =values[k,]
        sel_x = (Xc_xb==val)
        sel_y = (Xc_yb==val)
      }else{
        val = t(as.matrix(values[k,]))
        sel_x = matrix(1,dim(Xc_xb)[1],1)
        sel_y = matrix(1,dim(Xc_yb)[1],1)
        for(ddd in 1:dimXc){
          sel_x =  sel_x & (Xc_xb[,ddd]==val[ddd])
          sel_y =  sel_y & (Xc_yb[,ddd]==val[ddd])
        }
        sel_x = matrix( sel_x,dim(Xc_xb)[1],1)
        sel_y = matrix( sel_y,dim(Xc_yb)[1],1)
      }
      Xp= matrix(Xncb[sel_x,],sum(sel_x),dimXnc);
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
      T_xy = min(n_x,n_y)
      T_n[1,k] =  T_xy

      Ybarre = sum(Yp*weights_yp);
      mat_Yk[k,1] <- Ybarre
      mat_Xk[k,] <- sum(Xp*weights_xp)
      Xp_s = Xp
      weights_xp_s = weights_xp


      if(version == "first"){
        ### handle the potentially different size
        if(n_x > n_y){
          sit = sample(1:n_x,n_y,replace = F)
          Xp = matrix(Xp[sit,],n_y,dimXnc)
          weights_xp =   weights_xp[sit]
          weights_xp  =  weights_xp /sum( weights_xp )

        }else if( n_y > n_x){
           sit = sample(1:n_y,n_x,replace = F)
          Yp =Yp[sit]
          weights_yp =   weights_yp[sit]
          weights_yp  =  weights_yp /sum( weights_yp )
        }
      }

      if(is.null(data_k) | is.null(kp0)){
        test0 = TRUE
      }else{
        if(meth=="adapt"){
          test0 = sum(is.na(kp0[,k]))==0
        }else{
          test0 = sum(is.na(kp0[k]))==0
        }
      }

      if(sum(sel_x)> lim & sum(sel_y) >  lim & test0){

        if(winsor ==TRUE){
          qu=quantile(Yp,   max( 0.9, 1-2*log(length(Yp))/length(Yp)) )
          Yp[Yp> qu] <- qu
        }


        if(version == "second"){


          Ybarre = sum(Yp*weights_yp);
          Ys0 = cbind(weights_yp,Yp-Ybarre)
          Ys1 = Ys0[order(Ys0[,2], decreasing = T),]
          indexes = Ys1[,1]>0
          Ysort = Ys1[indexes,2]
          weights_yp = Ys1[indexes,1]
          ##
          for_critY = cumsum(weights_yp*Ysort)
          ## to handle ties
          if(is.null(sample1)){
            grid_I_k = sort(unique(c(1,0,1- cumsum(weights_yp),1- cumsum(weights_xp)) ), decreasing = T)
            grid_I[[k]] <- grid_I_k
          }else{
            grid_I_k =  sort(unique(c(grid_I_ref[[k]],1- cumsum(weights_yp),1- cumsum(weights_xp)) ), decreasing = T)
           }
          for_critY = approx(c(1,1- cumsum(weights_yp)) , c(for_critY,0), xout= grid_I_k  )$y
        }else{
          Ys0 = cbind(weights_yp,Yp-Ybarre)
          Ys0 = Ys0[order(Ys0[,2], decreasing = T),]
          weights_yp = Ys0[,1]
          Ysort = Ys0[,2]
          for_critY=cumsum(weights_yp*Ysort);
        }
        if(!is.null(sample1)){
          ### handle the potentially different size
          if(dimXc==1){
            sel_x = (Xc_x==val)
            sel_y = (Xc_y==val)
          }else{
            sel_x = matrix(1,dim(Xc_x)[1],1)
            sel_y = matrix(1,dim(Xc_y)[1],1)
            for(ddd in 1:dimXc){
              sel_x =  sel_x & (Xc_x[,ddd]==val[ddd])
              sel_y =  sel_y & (Xc_y[,ddd]==val[ddd])
            }
            sel_x = matrix( sel_x,dim(Xc_x)[1],1)
            sel_y = matrix( sel_y,dim(Xc_y)[1],1)
          }
          Yp_ref=  Y[sel_y]
          Xp_ref =  matrix(Xnc[sel_x,],sum(sel_x),dimXnc);
          weights_yp_ref =    weights_ys[sel_y]
          weights_yp_ref  =  weights_yp_ref /sum( weights_yp_ref )
          weights_xp_ref  =  weights_xs[sel_x]
          weights_xp_ref  =  weights_xp_ref /sum( weights_xp_ref )
        }else{
          Yp_ref= NULL
          Xp_ref=NULL
          weights_yp_ref=NULL
          weights_xp_ref=NULL
        }
        if(meth=="min"){
          if(is.null(kp0)){
            sam0_kp0=sam0
          }else{
            if(is.null(data_k)){
              sam0_kp0= cbind(rep(kp0,dim(sam0)[1]),sam0)
            }else{
              sam0_kp0= cbind(rep(kp0[k],dim(sam0)[1]),sam0)
            }
          }
        }else{
          #adapt
          if(is.null(data_k)){
            sam0_kp0= cbind(rep(kp0,dim(sam0)[1]),sam0)
          }else{
            sam0_kp0= cbind(kp0[,k],sam0)
          }
        }

        sam0_kp0 = cbind(1:dim(sam0_kp0)[1],sam0_kp0)
       es =  T_xy_all^(-boot_par)

        if(is.null(sample1)){
          bsharp_beta2_j <- na.omit(t(apply(sam0_kp0,1,compute_ratio_point,Xp=Xp,Yp= Ysort,for_critY=for_critY,
                                            dimXnc=dimXnc,weights_xp=weights_xp,weights_yp=weights_yp,T_xy_all,T_xy, es, version,trunc=trunc,
                                            grid_I = grid_I_k )))

          bsharp_beta2 = matrix(NA,1,dim(sam0_kp0)[1])
          ratio_ref_k = vector("list") #matrix(NA,length(Yp),dim(sam0_kp0)[1])
          # kk =1
          for(kk in 1:dim(sam0_kp0)[1]){
            bsharp_beta2[,kk] =bsharp_beta2_j[[kk]][["lambda"]]
            ratio_ref_k[[kk]] =bsharp_beta2_j[[kk]][["ratio"]]
          }
          ratio_ref[[k]] <- ratio_ref_k

        }else{
          ratio_ref_k = ratio_ref[[k]]
          bsharp_beta2 <- na.omit(t(apply(sam0_kp0,1,compute_ratio,Xp=Xp,Yp= Ysort,for_critY=for_critY,
                                          dimXnc=dimXnc,weights_xp=weights_xp,weights_yp=weights_yp,T_xy_all,T_xy,es, version,trunc=trunc,
                                          ratio_ref=ratio_ref_k, Yp_ref= Yp_ref,
                                          Xp_ref=Xp_ref,weights_yp_ref=weights_yp_ref,weights_xp_ref=weights_xp_ref, Opt=Opt, grid_I =  grid_I_k)))
        }



        if(dimXnc==1){
          mat_var_unc[k,] <- bsharp_beta2
          mat_var_low[k,] <- - rev(mat_var_unc[k,])
          bsharp_beta2_m  = rev(bsharp_beta2)
        }else{
          #### compute S(-q)
          #min
          if(meth=="min"){
            if(is.null(kp0)){
              sam0_kp0_m= cbind(sam0_kp0[,1],-sam0_kp0[,2])
            }else{
              if(is.null(data_k)){
                sam0_kp0_m= cbind(rep(kp0,dim(sam0)[1]),-sam0)
              }else{
                sam0_kp0_m= cbind(rep(kp0[k],dim(sam0)[1]),-sam0)
              }
            }
          }else{
            #adapt
            if(is.null(data_k)){
              sam0_kp0_m= cbind(rep(kp0,dim(sam0)[1]),-sam0)
            }else{
              sam0_kp0_m= cbind(kp0[,k],-sam0)
            }
          }
          sam0_kp0_m = cbind(1:dim(sam0_kp0_m)[1],sam0_kp0_m)

          #################
          if(is.null(sample1)){
            bsharp_beta2_m_j <- na.omit(t(apply(sam0_kp0_m,1,compute_ratio_point,Xp=Xp,Yp= Ysort,for_critY=for_critY,
                                                dimXnc=dimXnc,weights_xp=weights_xp,weights_yp=weights_yp,T_xy_all,T_xy, es, version,trunc=trunc,
                                                grid_I =  grid_I_k)))

            bsharp_beta2_m = matrix(NA,1,dim(sam0_kp0)[1])
            ratio_ref_m_k =  vector("list")
            for(kk in 1:dim(sam0_kp0_m)[1]){
              bsharp_beta2_m[,kk] = bsharp_beta2_m_j[[kk]][["lambda"]]
              ratio_ref_m_k[[kk]] =bsharp_beta2_m_j[[kk]][["ratio"]]
            }

            ratio_ref_m[[k]] <- ratio_ref_m_k


            ######################"
          }else{

            ratio_ref_m_k = ratio_ref_m[[k]]
            grid_I_ref_k = grid_I_ref[[k]]
            bsharp_beta2_m <- na.omit(t(apply(sam0_kp0_m,1,compute_ratio,Xp=Xp,Yp= Ysort,for_critY=for_critY,
                                              dimXnc=dimXnc,weights_xp=weights_xp,weights_yp=weights_yp,T_xy_all,T_xy, es, version,trunc=trunc,
                                              ratio_ref=ratio_ref_m_k, Yp_ref= Yp_ref, Xp_ref=Xp_ref,
                                              weights_yp_ref=weights_yp_ref,weights_xp_ref=weights_xp_ref, Opt=Opt, grid_I =  grid_I_ref_k)))
          }
        }

        ### save the unconstrained value of  S(q)
        mat_var_unc[k,] <- bsharp_beta2
        mat_var_low[k,] <- -bsharp_beta2_m
        EU = FALSE
        EL=  FALSE


        if(!is.null(kp0)){
          X_k = ((sam0%*%t(Xp_s))*(matrix(1,dim(sam0)[1],1)%*%weights_xp_s))%*%matrix(1,dim(Xp_s)[1],1)
          dXk = X_k- X_0
          dYk = Ybarre- Y0

          # dir=1
          for( dir in 1:dim(sam0)[1]){

            tt = t.test((sam0[dir,]%*%t(Xp_s)), y = (sam0[dir,]%*%t(Xp0)),alternative = c("two.sided", "less", "greater"),
                        mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
            tests[dir,k] = tt$p.value #< 0.001

          }


          Dmat_Yk[k,1] <- dYk
          Dmat_Xk[[k]] <- dXk

          mean_ratio =  dYk/dXk

          off = TRUE
          if(!is.null(sample1) & Opt=="boot"){
            mean_ratio_ref =  Dmat_Yk_ref[k,1]/Dmat_Xk_ref[[k]]
            es = T_xy_all^(-boot_par)

            mean_ratio =  mean_ratio_ref + (mean_ratio -mean_ratio_ref)*es*sqrt(T_xy_all)
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
              EL = (c_sign[indk]>0 & dXk<0  ) |  (c_sign[indk]<0 & dXk>0  )
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
    }

    mat_var_unc1<- mat_var_unc[ind ,]
    mat_var_unc1 <- matrix(mat_var_unc1,length(ind),dim(sam0_kp0)[1])
    mat_var1 <-   mat_var[ind,]
    mat_var1 <- matrix(mat_var1,length(ind),dim(sam0_kp0)[1])
    mat_var_low1 <-   mat_var_low[ind,]
    mat_var_low1 <- matrix(mat_var_low1,length(ind),dim(sam0_kp0)[1])

    if(!is.null(values)){

      mat_var <- apply(  mat_var1,2,min)
      mat_var_low <- apply(  mat_var_low1,2,max)
      mat_var_unc <- apply(  mat_var_unc1,2,min)
    }

    if(!is.null(nc_sign)){
      ee = eye(dimXnc)
      for(jj in 1:dimXnc){
        if(nc_sign[jj]!=0){
          cond = (matrix(ee[jj,]*nc_sign[jj],1,dimXnc)%*%t(sam0))<0
          mat_var[cond]<- pmin(mat_var[cond],0)
          mat_var_low[!cond]<- pmax(mat_var_low[!cond],0)
        }
      }
    }


  #### without Xc
  }else{

    mat_var_unc= matrix(0,1,dim(sam0)[1])
    mat_var_low= matrix(0,1,dim(sam0)[1])
    if(!is.null(values)){
      mat_Yk= matrix(NA,dim(values)[1],1)
      mat_Xk= matrix(NA,dim(values)[1],dimXnc)
    }

    Xp=matrix(Xncb, dim(Xncb)[1] ,dimXnc )
    Yp = matrix(Yb, dim(Yb)[1],1)
    weights_yp =   weights_y/sum( weights_y)
    weights_xp =   weights_x/sum( weights_x )

    n_x = dim(Xp)[1]
    n_y = dim(Yp)[1]
    if(is.null(T_xy_all)){
      T_xy  = min(n_x,n_y)
    }else{
      T_xy = T_xy_all
    }

    if(version=="first"){

      ### handle the potentially different size
      if(n_x > n_y){
        sit = sample(1:n_x,n_y,replace = F)
        Xp = matrix(Xp[sit,],n_y,dimXnc)
        weights_xp =   weights_x[sit]
        weights_xp  =  weights_xp /sum( weights_xp )
      }else if( n_y > n_x){
        sit = sample(1:n_y,n_x,replace = F)
        Yp =Yp[sit]
        weights_yp =   weights_y[sit]
        weights_yp  =  weights_yp /sum( weights_yp )
      }
    }


    if(winsor ==TRUE){
      qu=quantile(Yp,   max( 0.9, 1-2*log(length(Yp))/length(Yp)) )
      Yp[Yp> qu] <- qu
    }


    if(version == "second"){

      Ybarre = sum(Yp*weights_yp);
      Ys0 = cbind(weights_yp,Yp-Ybarre)
      Ys1 = Ys0[order(Ys0[,2], decreasing = T),]
      indexes0 = Ys1[,1]>0
      Ysort = Ys1[ indexes0 ,2]
      weights_yp = Ys1[ indexes0 ,1]
      weights_yp = weights_yp /sum( weights_yp )
      for_critY = cumsum(weights_yp*Ysort)
      if(is.null(sample1)){
        grid_I =  sort(unique(c(1,0,1- cumsum(weights_yp),1- cumsum(weights_xp)) ), decreasing = T)
      }else{
        grid_I =  sort(unique(c(grid_I_ref,1- cumsum(weights_yp),1- cumsum(weights_xp)) ), decreasing = T)
      }
      for_critY = approx( c(1,1- cumsum(weights_yp)) , c(for_critY,0), xout= grid_I )$y

    }else{

      Ybarre = sum(Yp*weights_yp);
      Ys0 = cbind(weights_yp,Yp-Ybarre)
      Ys0 = Ys0[order(Ys0[,2], decreasing = T),]
      weights_yp = Ys0[,1]
      Ysort = Ys0[,2]
      for_critY=cumsum(weights_yp*Ysort);

      grid_I=NULL
    }



    if(!is.null(values)){
      mat_Yk[k,1] <- Ybarre
      mat_Xk[k,] <- sum(Xp*weights_xp)
    }


    if(!is.null(sample1)){
      ### handle the potentially different size
      Yp_ref= Y
      Xp_ref = Xnc
      weights_yp_ref =   weights_ys
      weights_xp_ref  =  weights_xs
    }else{
      # for_critY_ref=NULL
      Yp_ref= NULL
      Xp_ref=NULL
      weights_yp_ref=NULL
      weights_xp_ref=NULL
    }

    if(meth=="min"){
      # compute S(q)
      if(is.null(kp0)){
        ## when select epsilon
        sam0_kp0= sam0
      }else{
        sam0_kp0= cbind(rep(kp0,dim(sam0)[1]),sam0)
      }
    }else{
      # adapt
      # # compute S(q)
      if(is.null(kp0)){
        ## when select epsilon
        sam0_kp0= sam0
      }else{
        sam0_kp0= cbind(kp0,sam0)
      }
      #
    }
    sam0_kp0 = cbind(1:dim(sam0_kp0)[1],sam0_kp0)
    if(is.null(es_all)){
      es = T_xy^(-boot_par)
    }else{
      es = es_all
    }
    ### if the point estimate, compute and save the ratio
    if(is.null(sample1)){
      bsharp_beta2 <- na.omit(t(apply(sam0_kp0,1,compute_ratio_point,Xp=Xp,Yp= Ysort,for_critY=for_critY,
                                      dimXnc=dimXnc,weights_xp=weights_xp,weights_yp=weights_yp,T_xy,T_xy,es, version,trunc=trunc,
                                      grid_I = grid_I)))

      mat_var = matrix(NA,1,dim(sam0_kp0)[1])
      ratio_ref = vector("list")
      for(k1 in 1:dim(sam0_kp0)[1]){
        mat_var[,k1]  = bsharp_beta2[[k1]][["lambda"]]
        ratio_ref[[k1]] =bsharp_beta2[[k1]][["ratio"]]
      }
      ### otherwise, simply compute, using the point ratio.
    }else{
      bsharp_beta2 <- t(apply(sam0_kp0,1,compute_ratio,Xp=Xp,Yp= Ysort,for_critY=for_critY,
                              dimXnc=dimXnc,weights_xp=weights_xp,weights_yp=weights_yp,T_xy,T_xy,es, version,trunc=trunc,
                              ratio_ref=ratio_ref, Yp_ref= Yp_ref, Xp_ref=Xp_ref,
                              weights_yp_ref=weights_yp_ref,weights_xp_ref=weights_xp_ref, Opt=Opt, grid_I = grid_I))

      mat_var <- bsharp_beta2
    }

    if(type=="both" | type=="low"){
      if(dimXnc==1){
        mat_var_low <-  - rev( mat_var)
        bsharp_beta2_m  = rev(bsharp_beta2)
      }else{
        # compute S(-q)
        if(is.null(kp0)){
          # min
          sam0_kp0_m=  sam0
        }else{
          # min
          if(meth=="min"){
            sam0_kp0_m= cbind(rep(kp0,dim(sam0)[1]),-sam0)
          }else{
            # adapt
            sam0_kp0_m= cbind(kp0,-sam0)
          }
        }
       sam0_kp0_m = cbind(1:dim(sam0_kp0_m)[1],sam0_kp0_m)


        if(is.null(sample1)){
          bsharp_beta2_m <- na.omit(t(apply(sam0_kp0_m,1,compute_ratio_point,Xp=Xp,Yp= Ysort,for_critY=for_critY,
                                            dimXnc=dimXnc,weights_xp=weights_xp,weights_yp=weights_yp,T_xy,T_xy,es, version,trunc=trunc,
                                            grid_I = grid_I )))

          mat_var_low = matrix(NA,1,dim(sam0_kp0)[1])
          ratio_ref_m =  vector("list")
          for(k1 in 1:dim(sam0_kp0_m)[1]){
            mat_var_low[,k1] =- bsharp_beta2_m[[k1]][["lambda"]]
            ratio_ref_m[[k1]] =bsharp_beta2_m[[k1]][["ratio"]]
          }

        }else{
          bsharp_beta2_m <- na.omit(t(apply(sam0_kp0_m,1,compute_ratio,Xp=Xp,Yp= Ysort,for_critY=for_critY,
                                            dimXnc=dimXnc,weights_xp=weights_xp,weights_yp=weights_yp,T_xy,T_xy,es, version ,trunc=trunc,
                                            ratio_ref=ratio_ref_m, Yp_ref= Yp_ref, Xp_ref=Xp_ref,
                                            weights_yp_ref=weights_yp_ref,weights_xp_ref=weights_xp_ref, Opt=Opt, grid_I = grid_I)))
          mat_var_low <-  -bsharp_beta2_m
        }
      }
    }

    ind = 1
  }

  if(!is.null(nc_sign)){
    ee = eye(dimXnc)
    for(jj in 1:dimXnc){
      if(nc_sign[jj]!=0){
        cond = (matrix(ee[jj,]*nc_sign[jj],1,dimXnc)%*%t(sam0))<0
        mat_var[ cond]<- 0
      }
    }
  }

  #### return the values of the stats Tinf/Tsup and S(q) for each q
  if(type=="both"){
    output <- vector("list")
    output[["upper"]] <- mat_var
    output[["lower"]] <- mat_var_low
    output[["unconstr"]] <- mat_var_unc
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
      output[["ratio_ref_m"]] <- ratio_ref_m
      output[["grid_ref"]] <- grid_I
    }

    #### return only the values of the stat Tinf for each q
  }else if(type=="low"){
    output <- mat_var_low
  }else{
    #### return only the values of the stat S(q) for each q
    output <- mat_var
  }

  return(output)

}
