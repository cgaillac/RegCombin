#' Function for the data-driven selection of the epsilon tuning parameter
#'
#' @param sam1 the matrix containing the directions q on which to compute the selected rule for epsilon(q)
#' @param kp If data_k =NULL, then epsilon is taken equal to kp.
#' @param Xc_x the common regressor on the dataset  (Xnc,Xc). Default is NULL.
#' @param Xnc the noncommon regressor on the dataset  (Xnc,Xc). No default.
#' @param Xc_y the common regressor on the dataset  (Y,Xc). Default is NULL.
#' @param Y the outcome variable. No default.
#' @param values  the different unique points of support of the common regressor Xc.
#' @param dimXc the dimension of the common regressors Xc.
#' @param dimXnc the dimension of the noncommon regressors Xnc.
#' @param nb_pts  the constant C in DGM for the epsilon_0, the lower bound on the grid for epsilon, taken equal to nb_pts*ln(n)/n. Default is 1 without regressors Xc, 3 with Xc.
#' @param Opt the option of the method to compute the critical values, either "boot" (numerical bootstrap) or "subsampling". Default is numerical bootstrap.
#' @param lim the limit number of observations under which we do no compute the conditional variance.
#' @param weights_x the sampling weights for the dataset (Xnc,Xc).
#' @param weights_y the sampling weights for the dataset (Y,Xc).
#' @param refs0 indicating the positions in the vector values corresponding to the components of betac.
#' @param trunc equal to 2, for the definition of epsilon.
#' @param boot_par the numerical bootstrap parameter. Default is 0.3 without Xc, 0.35 with Xc.
#' @param data_k the number of points for the grid search on epsilon. Default is 30. If NULL, then epsilon is taken fixed equal to kp.
#' @param C0 the upper bound on the grid for epsilon. Default is 0.5.
#' @param c_sign sign restrictions on the commonly observed regressors: -1 for a minus sign, 1 for a plus sign, 0 otherwise. Default is NULL, i.e. no constraints.
#' @param nc_sign sign restrictions on the non-commonly observed regressors Xnc: -1 for a minus sign, 1 for a plus sign, 0 otherwise. Default is NULL, i.e. no constraints.
#' @param meth the method for the choice of epsilon, either "adapt", i.e. adapted to the direction or "min" the minimum over the directions. Default is "adapt".
#' @param nbCores number of cores for the parallel computation. Default is 1.
#' @param winsor indicates if winsorisation. Default is FALSE.
#' @param version version of the computation of the ratio, "first" is a degraded version but fast; "second" is a correct version but slower. Default is "second".
#' @param alpha the level for the confidence regions. Default is 0.05.
#'
#' @return
#' a matrix containing the values of the selected epsilon(q) for q directions in sam1.
#'
#'
#' @export
#'
#' @examples
select_epsilon <- function(sam1, kp,Xc_x,Xnc,Xc_y,Y,
                           values,dimXc,dimXnc,
                           nb_pts,Opt,lim ,weights_x,weights_y,
                           refs0,
                           trunc, boot_par=0.35, data_k=30,C0=0.5,c_sign=NULL,nc_sign=NULL, meth="adapt",
                           nbCores=1, winsor = FALSE, version = "first", alpha=0.05){



  if(is.null(weights_x)){
    weights_x= rep(1/dim(Xnc)[1],dim(Xnc)[1])
  }
  if(is.null(weights_y)){
    weights_y= rep(1/length(Y),length(Y))
  }


  limit=lim
  alpha1= alpha
  q95 <-  function(x){quantile(x,1-alpha,na.rm=T)}
  q05 <-  function(x){quantile(x,alpha,na.rm=T)}
  version0 = "first"


  vart0 = NULL
  vart = NULL
  vart0_v = NULL
  vart2_v = NULL
  vart3_v = NULL
  vart1_v = NULL
  vart_v = NULL
  Bsamp0=50
  boot_par1= boot_par

  if(is.null(values)){

    n_x = dim(Xnc)[1]
    n_y = dim(Y)[1]
    n_xy = min(n_x,n_y)
    T_xy =   n_xy
    eps0 = nb_pts*log(T_xy)
    kp_grid = seq(eps0/T_xy,C0,length.out=data_k)
    sam0_kp0 =matrix(NA,dim(sam1)[1]*length(kp_grid),dimXnc+1)
    for(j0 in 1:dim(sam1)[1]){
      cur_ind = ((j0-1)*length(kp_grid)+1):(j0*length(kp_grid))
      sam0_kp0[ cur_ind,1] <- kp_grid
      sam0_kp0[ cur_ind,2:(dimXnc+1)] <- matrix(1,length(kp_grid),1)%*%sam1[j0,]
    }
    ### compute point estimate
    mat_var_out <- compute_radial(sample1 =NULL,Xc_x,Xnc,Xc_y,Y,
                             values=NULL,dimXc,dimXnc,
                             nb_pts, sam0= sam0_kp0, kp0=NULL, data_k ,Opt,lim =limit,
                             weights_x,weights_y,
                             c_sign = NULL, nc_sign  =NULL,refs0,type="both",meth=meth,
                             trunc=trunc,boot_par=boot_par,  winsor= winsor, version =version0)
    if(nbCores>1){
      res0 <- sfLapply(1:Bsamp0, compute_radial,Xc_x,Xnc,Xc_y,Y,values=NULL,dimXc,dimXnc,nb_pts,
                       sam0= sam0_kp0, kp0=NULL, data_k ,Opt,lim = limit ,weights_x,weights_y,c_sign = NULL, nc_sign  =NULL,
                       refs0=NULL,type="up",meth=meth,trunc=trunc,boot_par=boot_par1, winsor= winsor, version =version0,
                       ratio_ref=    mat_var_out[["ratio_ref"]], ratio_ref_m=    mat_var_out[["ratio_ref_m"]],
                       Dmat_Yk_ref=   mat_var_out[["DYk"]], Dmat_Xk_ref=  mat_var_out[["DXk"]],grid_I_ref= mat_var_out[["grid_ref"]]  )

    }else{
      res0 <- lapply(1:Bsamp0, compute_radial,Xc_x,Xnc,Xc_y,Y,values=NULL,dimXc,dimXnc,nb_pts,
                     sam0= sam0_kp0, kp0=NULL, data_k ,Opt,lim = limit ,weights_x,weights_y,c_sign = NULL, nc_sign  =NULL,
                     refs0=NULL,type="up",meth=meth,trunc=trunc,boot_par=boot_par1,   winsor= winsor, version =version0,
                     ratio_ref=   mat_var_out[["ratio_ref"]], ratio_ref_m=    mat_var_out[["ratio_ref_m"]],  Dmat_Yk_ref=  mat_var_out[["DYk"]],
                     Dmat_Xk_ref=  mat_var_out[["DXk"]], grid_I_ref= mat_var_out[["grid_ref"]] )
    }




    # compute_radial,Xc_x,Xnc,Xc_y,Y,values=NULL,dimXc,dimXnc,nb_pts,
    # sam0= sam0_kp0, kp0=NULL,Opt,lim = limit ,weights_x,weights_y,c_sign = NULL, nc_sign  =NULL,
    # refs0=NULL,type="up",meth=meth,trunc=trunc,boot_par=boot_par1,   winsor= winsor, version =version0,
    # ratio_ref=   mat_var_out[["ratio_ref"]]
    # ratio_ref_m=    mat_var_out[["ratio_ref_m"]]
    # Dmat_Yk_ref=  mat_var_out[["DYk"]]
    # Dmat_Xk_ref=  mat_var_out[["DXk"]]
    # grid_I_ref= mat_var_out[["grid_ref"]]



    mat_varb= matrix(0,Bsamp0,dim(sam0_kp0)[1])
    for(b in 1:Bsamp0){
      mat_varb[b,] <-  res0[[b]]
    }
    S_e = matrix(1,Bsamp0,1)%*%mat_var_out$upper

    if(Opt!="boot"){
      bsx = sampling_rule(n_x)
      bsy = sampling_rule(n_y)
      bs0 = sqrt(min( bsx, bsy)/ T_xy)
      vart0_v= - apply((mat_varb -  S_e)*  bs0,2,q05)
      scal = (mat_varb -  S_e)/  matrix(1,Bsamp0,1)%*%sqrt(vart0_v)
      vart2_v=mat_var_out$upper - apply((mat_varb -  S_e)*  bs0,2,q05)
    }else{
      es =  T_xy^(-boot_par1)  ## epsilon_n
      term = es^(-1)/sqrt(T_xy)
      vart0_v= - apply((mat_varb -  S_e)*term,2,q05)
      scal = (mat_varb -  S_e)/  matrix(1,Bsamp0,1)%*%sqrt(vart0_v)
      vart2_v=mat_var_out$upper - apply((mat_varb -  S_e)*term ,2,q05)
    }
    out_sim = cbind(sam0_kp0,matrix(mat_var_out[[1]],dim(sam0_kp0)[1],1),matrix(vart0_v,dim(sam0_kp0)[1],1), matrix(vart2_v,dim(sam0_kp0)[1],1))
    out_sim_st = matrix(NA,dim(sam0_kp0)[1],5+dimXnc)
    kp1= matrix(0.001,dim(sam1)[1],1)
    C1 = rep(1,length(kp_grid))
    # j0=1
    save_inter=matrix(0,1,dim(out_sim)[1])
    for(j0 in 1:dim(sam1)[1]){
      cur_ind = ((j0-1)*length(kp_grid)+1):(j0*length(kp_grid))
      out_sim_p = out_sim[cur_ind,]
      out_sim_p = cbind(  out_sim_p,  out_sim_p[,2+dimXnc] + C1* out_sim_p[,3+dimXnc])
      out_sim_st[cur_ind,] <- out_sim_p
      save_inter =  save_inter + out_sim_p[,5+dimXnc]
      kp1[j0,1] = kp_grid[which.min(out_sim_p[,5+dimXnc])]
    }
    if(meth=="min"){
      kp_fin = min(kp1)
    }else{
      # adapt
      kp_fin = kp1
    }
  }else{

    if(meth=="min"){
      # min
      kp_fin = matrix(NA,dim(values)[1],1)
    }else{
      # #adapt
      kp_fin = matrix(NA,dim(sam1)[1],dim(values)[1])
    }


    n_x_all =  dim(Xnc)[1]
    n_y_all = dim(Y)[1]
    n_xy_all= min(n_x_all,n_y_all)
    T_xy_all =  n_xy_all
    for(k in 1:dim(values)[1]){

      if(dimXc==1){
        val = values[k,]
        sel_x =   (Xc_x ==val)
        sel_y =   (Xc_y ==val)
      }else{
        val = t(as.matrix(values[k,]))
        sel_x = matrix(1,dim(Xc_x)[1],1)
        sel_y = matrix(1,dim(Xc_y)[1],1)
        for(ddd in 1:dimXc){
          sel_x =  sel_x & (Xc_x[,ddd]==val[ddd])
          sel_y =  sel_y & (Xc_y[,ddd]==val[ddd])
        }
      }

      if(sum(sel_x)> lim & sum(sel_y) >  lim){
        Xp =  matrix(Xnc[sel_x,],sum(sel_x),dimXnc);
        Yp =  matrix(Y[sel_y],sum(sel_y),1)

        weights_yp =   weights_y[sel_y]
        weights_yp  =  weights_yp /sum( weights_yp )
        weights_xp =   weights_x[sel_x]
        weights_xp  =  weights_xp /sum( weights_xp )

        n_x = sum(sel_x)
        n_y = sum(sel_y)
        n_xy=min(n_x,n_y)
        T_xy = n_xy
        eps0 = nb_pts*log(T_xy)
        kp_grid = seq(eps0/T_xy,C0,length.out=data_k)
        sam0_kp0 =matrix(NA,dim(sam1)[1]*length(kp_grid),dimXnc+1)
        for(j0 in 1:dim(sam1)[1]){
          cur_ind = ((j0-1)*length(kp_grid)+1):(j0*length(kp_grid))
          sam0_kp0[ cur_ind,1] <- kp_grid
          sam0_kp0[ cur_ind,2:(dimXnc+1)] <- matrix(1,length(kp_grid),1)%*%sam1[j0,]
        }
        ### compute point estimate
        df = diff(sort(Xp))
        delta0 = min(abs(df[df>0]),na.rm=T)/n_x
        Xp = matrix(Xp + runif(n_x,-delta0,delta0), n_x, dimXnc)

        df = diff(sort(Yp))
        delta0 = min(abs(df[df>0]),na.rm=T)/n_y
        Yp = matrix(Yp + runif(n_y,-delta0,delta0),n_y ,1)

        # Y= Yp
        # Xnc=Xp

        mat_var_out <- compute_radial(sample1 =NULL,Xc_x=NULL,Xp,Xc_y=NULL,Yp,
                                 values=NULL,dimXc=0,dimXnc,
                                 nb_pts, sam0= sam0_kp0, kp0 = NULL, data_k,Opt,lim =limit,
                                 weights_xp,weights_yp,
                                 c_sign = NULL, nc_sign  =NULL,refs0,type="both",meth=meth,trunc=trunc,boot_par=boot_par,
                                 winsor= winsor, version =version0)



        if(nbCores>1){
          res0 <- sfLapply(1:Bsamp0, compute_radial,Xc_x=NULL,Xnc=Xp,Xc_y=NULL,Y=Yp,values=NULL,dimXc=0,dimXnc,nb_pts,
                           sam0=  sam0_kp0 , kp0 =NULL, data_k, Opt,lim = limit ,weights_x=weights_xp,weights_y=weights_yp,c_sign = NULL, nc_sign  =NULL,
                           refs0=NULL,type="up",meth=meth,trunc=trunc,boot_par=boot_par,
                           winsor= winsor, version =version0,
                           ratio_ref=   mat_var_out[["ratio_ref"]], ratio_ref_m=   mat_var_out[["ratio_ref_m"]],  Dmat_Yk_ref=   mat_var_out[["DYk"]], Dmat_Xk_ref=  mat_var_out[["DXk"]],
                           grid_I_ref= mat_var_out[["grid_ref"]], es_all =T_xy_all^(-boot_par) , T_xy_all = T_xy_all )

        }else{
          res0 <- lapply(1:Bsamp0, compute_radial,Xc_x=NULL,Xnc=Xp,Xc_y=NULL,Y=Yp,values=NULL,dimXc=0,dimXnc,nb_pts,
                         sam0=  sam0_kp0 , kp0 =NULL, data_k ,Opt,lim = limit ,weights_x=weights_xp,weights_y=weights_yp,c_sign = NULL, nc_sign  =NULL,
                         refs0=NULL,type="up",meth=meth,trunc=trunc,boot_par=boot_par,
                          winsor= winsor, version =version0,
                         ratio_ref=   mat_var_out[["ratio_ref"]], ratio_ref_m=   mat_var_out[["ratio_ref_m"]],  Dmat_Yk_ref=   mat_var_out[["DYk"]], Dmat_Xk_ref=  mat_var_out[["DXk"]],
                        grid_I_ref= mat_var_out[["grid_ref"]] ,  es_all =T_xy_all^(-boot_par) , T_xy_all = T_xy_all)
        }


        mat_varb= matrix(0,Bsamp0,dim(sam0_kp0)[1])
        for(b in 1:Bsamp0){
          mat_varb[b,] <-  res0[[b]]
        }

        S_e = matrix(1,Bsamp0,1)%*%mat_var_out[[1]]

        if(Opt!="boot"){
          bsx = sampling_rule(n_x)
          bsy = sampling_rule(n_y)
          bs0 = sqrt(min( bsx, bsy)/T_xy)
          vart0_v= - apply((mat_varb -  S_e)*  bs0,2,q05)
          scal = (mat_varb -  S_e)/  matrix(1,Bsamp0,1)%*%sqrt(vart0_v)
          vart2_v=mat_var_out[[1]] - apply((mat_varb -  S_e)*  bs0,2,q05)
        }else{
          es = T_xy_all^(-boot_par)  ## epsilon_n
          vart0_v= - apply((mat_varb -  S_e)*es^(-1)/sqrt(T_xy_all),2,q05)
          scal = (mat_varb -  S_e)/  matrix(1,Bsamp0,1)%*%sqrt(vart0_v)
          vart2_v=mat_var_out[[1]] - apply((mat_varb -  S_e)*es^(-1)/sqrt(T_xy_all),2,q05)
        }


        out_sim = cbind(sam0_kp0,matrix(mat_var_out[[1]],dim(sam0_kp0)[1],1),matrix(vart0_v,dim(sam0_kp0)[1],1), matrix(vart2_v,dim(sam0_kp0)[1],1))
        out_sim_st = matrix(NA,dim(sam0_kp0)[1],5+dimXnc)
        kp1= matrix(0.001,dim(sam1)[1],1)
        C1 = rep(1,length(kp_grid))

        j0=1
        for(j0 in 1:dim(sam1)[1]){
          cur_ind = ((j0-1)*length(kp_grid)+1):(j0*length(kp_grid))
          out_sim_p = out_sim[cur_ind,]
          out_sim_p = cbind(  out_sim_p,  out_sim_p[,3] + C1* out_sim_p[,4])
          out_sim_st[cur_ind,] <- out_sim_p
          kp1[j0,1] = kp_grid[which.min(out_sim_p[,5+dimXnc])]
        }

        if(meth=="min"){
          # min
          kp_fin[k,1] = min(kp1)
        }else{
          # adapt
          kp_fin[,k] = kp1
        }

      }

    }
  }

  return(kp_fin)
}
