#' Function computing all the different bounds : DGM and/or Variance
#'
#' @param Ldata dataset containing (Y,Xc) where Y is the outcome, Xc are potential common regressors.
#' @param Rdata dataset containing (Xnc,Xc) where Xnc are the non commonly observed regressors, Xc are potential common regressors.
#' @param out_var label of the outcome variable Y.
#' @param nc_var label of the non commonly observed regressors Xnc.
#' @param c_var label of the commonly observed regressors Xc.
#' @param nc_sign sign restrictions on the non-commonly observed regressors Xnc: -1 for a minus sign, 1 for a plus sign, 0 otherwise. Default is NULL, i.e. no constraints.
#' @param c_sign sign restrictions on the commonly observed regressors: -1 for a minus sign, 1 for a plus sign, 0 otherwise. Default is NULL, i.e. no constraints.
#' @param nbCores number of cores for the parallel computation. Default is 1.
#' @param methods method used for the bounds: "DGM" (Default) and/or "Variance".
#' @param kp If data_k =NULL, then epsilon is taken equal to kp.
#' @param nb_pts the constant C in DGM for the epsilon_0, the lower bound on the grid for epsilon, taken equal to nb_pts*ln(n)/n. Default is 1 without regressors Xc, 3 with Xc.
#' @param Opt the option of the method to compute the critical values, either "boot" (numerical bootstrap) or "subsampling". Default is numerical bootstrap.
#' @param Bsamp the number of bootstrap/subsampling replications. Default is 1000.
#' @param data_k the number of points for the grid search on epsilon. Default is 30. If NULL, then epsilon is taken fixed equal to kp.
#' @param C0 the upper bound on the grid for epsilon. Default is 0.5.
#' @param list_ex the list of the data that should not be used on the different clusters. This reduces the computational cost. Default is empty list.
#' @param weights_x the sampling weights for the dataset (Xnc,Xc). Default is NULL.
#' @param weights_y  the sampling weights for the dataset (Y,Xc). Default is NULL.
#' @param trunc equal to 2, for the definition of epsilon.
#' @param boot_par the numerical bootstrap parameter. Default is 0.3 without Xc, 0.35 with Xc.
#' @param meth the method for the choice of epsilon, either "adapt", i.e. adapted to the direction or "min" the minimum over the directions. Default is "adapt".
#' @param modeNA indicates if NA introduced if the interval is empty. Default is FALSE.
#' @param winsor indicates if winsorisation. Default is FALSE.
#' @param version version of the computation of the ratio, "first" is a degraded version but fast; "second" is a correct version but slower. Default is "second".
#' @param projections if FALSE compute the identified set along some directions or the confidence regions. Default is FALSE
#'
#' @return a list containing, in order:
#'
#'  - DGM_complete or Variance_complete : the complete outputs of the functions DGM_bounds or Variance_bounds.
#'
#'  and additional pre-treated outputs, replace below "method" by either "DGM" or "Variance":
#'
#' - methodCI:  the confidence region on the betanc without sign constraints
#'
#' - methodpt:  the bounds point estimates on the betanc without sign constraints
#'
#' - methodCI_sign:  the confidence region on the betanc with sign constraints
#'
#' - methodpt_sign:  the bounds point estimates on the betanc with sign constraints
#'
#' - methodkp: the values of epsilon(q)
#'
#' - methodbeta1: the confidence region on the betac corresponding to the common regressors Xc without sign constraints
#'
#' - methodbeta1_pt: the bounds point estimates on the betac corresponding to the common regressors Xc without sign constraints
#'
#' - methodbeta1_sign: the confidence region on the betac corresponding to the common regressors Xc with sign constraints
#'
#' - methodbeta1_sign_pt: the bounds point estimates on the betac corresponding to the common regressors Xc with sign constraints
#'
#' @export
#'
#' @examples
#'
regCombin <- function(Ldata, Rdata,
                      out_var, nc_var, c_var =NULL,
                      nc_sign = NULL, c_sign = NULL,
                      nbCores=1,
                      methods=c("DGM"),
                      kp=0.5, nb_pts=1, Opt="boot",Bsamp=1000,data_k=30,C0=0.5,
                      list_ex=c(),
                      weights_x = NULL,weights_y = NULL, trunc=2,boot_par=0.3, meth="adapt",
                      modeNA=FALSE,  winsor= FALSE, version = "second",  projections= FALSE){

  set = FALSE
  #### get sizes
  dimXc = length(c_var)
  dimXnc = length(nc_var)
  ####
  projections0 = projections
  if(dimXc!=0){
    if(dimXc!=1){
      values = create_values(dimXc,c_var,Rdata)
      refs0=NULL
      for(j in 1:dim(values)[1]){
        if(sum(values[j,]==0)==(dim(values)[2]-1)){
          refs0 = c( refs0,j )
        }
      }
    }else{
      values = matrix(create_values(dimXc,c_var,Rdata))
      refs0 = (1:length(values))[values>0]
    }

  }else{
    values = NULL
    s= NULL
    refs0 = NULL
  }

  weights_xs =  weights_x
  weights_ys =  weights_y
  output <- vector("list")
  idx_output = 1
  names_output = NULL

  for(  method in methods ){
    if( method=="DGM"){

      sam0 <- eye(dimXnc)
      sam0 <- rbind(-sam0,sam0)
      # if(dimXnc==1){projections0=FALSE}else{projections0=TRUE}
      eps0 = 0


      out <- DGM_bounds(Ldata, Rdata, values,sam0,refs0,
                        out_var,  nc_var, c_var,
                        nc_sign, c_sign,
                        nbCores=nbCores,
                        kp,nb_pts, Opt,Bsamp=Bsamp,data_k,C0,
                        list_ex=list_ex,
                        weights_x =weights_x,weights_y =weights_y,outside=FALSE,trunc= trunc,boot_par=boot_par, meth=meth,
                        modeNA =modeNA, winsor = winsor, version =version , projections=projections)


      output[[paste0(method,"_complete")]] <-  out


      if(dimXnc==1){
        mt_sharpCI_proj_sign <- out$ci$upper
        mt_sharpCI_proj_sign_low <- out$ci$lower
        mt_sharpCI_proj  <-  out$ci$unconstr

        mt_sharp0_proj_sign  <- out$point$upper
        mt_sharp0_proj_sign_low  <- out$point$lower
        mt_sharp0_proj  <- out$point$unconstr

        hull0 <-  mt_sharpCI_proj_sign
        hull0_low <-  mt_sharpCI_proj_sign_low


        out00 = matrix(0,dimXnc,2)

        for(id0 in 1:dimXnc){
          out00[id0,1] =-1* mt_sharpCI_proj_sign[id0]*(mt_sharpCI_proj_sign_low[dimXnc+id0]==0) +
            (mt_sharpCI_proj_sign_low[dimXnc+id0]!=0)*mt_sharpCI_proj_sign_low[dimXnc+id0]
          out00[id0,2] = (1* mt_sharpCI_proj_sign[dimXnc+id0] )*(mt_sharpCI_proj_sign_low[id0]==0) -
            (mt_sharpCI_proj_sign_low[id0]!=0)*mt_sharpCI_proj_sign_low[id0]
        }
        output[[paste0(method,"CI_sign")]] <- out00



        for(id0 in 1:dimXnc){
          out00[id0,1] = -1* mt_sharp0_proj_sign[id0]*(mt_sharp0_proj_sign_low[dimXnc+id0]==0) +
            (mt_sharp0_proj_sign_low[dimXnc+id0]!=0)*mt_sharp0_proj_sign_low[dimXnc+id0]
          out00[id0,2] = 1*mt_sharp0_proj_sign[dimXnc+id0]*(mt_sharp0_proj_sign_low[id0]==0) -
            (mt_sharp0_proj_sign_low[id0]!=0)*mt_sharp0_proj_sign_low[id0]
        }
        output[[paste0(method,"pt_sign")]] <- out00

     }else{
       ### dimXnc >1


        if(projections==TRUE){
          ## support function
          ## sign constraints not implemented yet in dimXnc >1
          output[[paste0(method,"CI_sign")]] <- NULL
          output[[paste0(method,"pt_sign")]] <- NULL

          v0 = NULL
          for(d in 1:dimXnc){
          v0 <- c(v0,paste0("q_",d))
          }
          out_ci <-  out$ci$support
          colnames(out_ci) <- c(v0,"Support")
          out_pt <-  out$point$support
          colnames(out_pt) <- c(v0,"Support")

          res_ci <- matrix(NA,dimXnc,2)
          res_pt <- matrix(NA,dimXnc,2)
          for(d in 1:dimXnc){
            sub <- out_ci[out_ci[,d]!=0,c(d,dimXnc+1)]
            res_ci[d,] <- c(apply(sub,1,prod))
            sub <- out_pt[out_pt[,d]!=0,c(d,dimXnc+1)]
            res_pt[d,] <- c(apply(sub,1,prod))
          }


          output[[paste0(method,"support_CI")]] <- out_ci
          output[[paste0(method,"support_pt")]] <- out_pt
          output[[paste0(method,"CI")]] <- res_ci
          output[[paste0(method,"pt")]] <- res_pt

        }else{
          #### the status is different/ convex hull points.
          mt_sharpCI_sign <- out$ci$upper
          mt_sharp0_sign <- out$point$upper

          output[[paste0(method,"CI_all_sign")]] <-  mt_sharpCI_sign
          output[[paste0(method,"pt_all_sign")]] <-  mt_sharp0_sign
          mt_sharpCI_proj <- out$ci$unconstr
          mt_sharp0_proj  <- out$point$unconstr
        }

      }

      output[[paste0(method,"kp_sign")]] <-   out$epsilon

      sam0 <- eye(dimXnc)
      sam0 <- rbind(-sam0,sam0)

      # bounds on betanc not yet implemented in dimXnc > 1
      if(dimXc>0 & dimXnc ==1){

        beta1K <-  out$ci$betac_ci
        beta1K_pt <-  out$point$betac_pt

      }else{
        beta1K = c(NA,NA)
        beta1K_pt = c(NA,NA)
      }

      output[[ paste0(method,"beta1_sign")]] <- beta1K
      output[[ paste0(method,"beta1_sign_pt")]] <- beta1K_pt
    }



    ##### For the variance bounds ####################################################
    if( method=="Variance"){

      sam0 <- eye(dimXnc)
      sam0 <- rbind(-sam0,sam0)

      if(dimXnc==1){projections0=FALSE}else{projections0=TRUE}


      out <- Variance_bounds(Ldata,Rdata,
                             out_var, c_var, nc_var,
                             c_sign, nc_sign,
                             values,sam0, refs0,
                             nb_pts,kp, Opt="boot",nbCores,Bsamp=1000,
                             list_ex=list_ex,weights_x =weights_x,weights_y =weights_y, projections=projections0)

      output[[paste0(method,"_complete")]] <-  out
      mt_sharpCI_proj_sign <- out$ci$upper
      mt_sharp0_proj_sign  <- out$point$upper
      mt_sharpCI_proj <- out$ci$unconstr
      mt_sharp0_proj  <- out$point$unconstr

      out00 = matrix(0,dimXnc,2)
      output[[paste0(method,"CI_all_sign")]] <- mt_sharpCI_proj_sign
      output[[paste0(method,"pt_all_sign")]] <-   mt_sharp0_proj_sign
    }


    #### 1) bounds on beta_nc without any sign constraints
    if(projections==FALSE | dimXnc==1){

      out00 = matrix(0,dimXnc,2)

      for(id0 in 1:dimXnc){
        if(method!="HP"){
          out00[id0,1] = -mt_sharpCI_proj[id0]
        }else{
          out00[id0,1] = mt_sharpCI_proj[id0]
        }
        out00[id0,2] = mt_sharpCI_proj[dimXnc+id0]
      }
      output[[paste0(method,"CI")]] <- out00

      for(id0 in 1:dimXnc){
        if(method!="HP"){
          out00[id0,1] = - mt_sharp0_proj[id0]
        }else{
          out00[id0,1] = mt_sharp0_proj[id0]
        }
        out00[id0,2] = mt_sharp0_proj[dimXnc+id0]
      }
      output[[paste0(method,"pt")]] <- out00
    }



    ######################### bounds on betac
    # # if(dimXnc==1){
    # hull0 <- mt_sharpCI_proj
    # hull00 <- mt_sharp0_proj
    # # }else{
    # #   hull0 <- hull_CI
    # #   hull00 <- hull_pt
    # # }
    # sam0 <- eye(dimXnc)
    # sam0 <- rbind(-sam0,sam0)


    if(dimXc>0 & dimXnc==1){

      if(method=="Variance"){

        beta1K <- out$ci$betac_ci_unc
        beta1K0 <- out$point$betac_pt_unc
        output[[paste0(method,"beta1_sign")]] <- out$ci$betac_ci
        output[[paste0(method,"beta1_sign_pt")]] <- out$point$betac_pt

      }else{

        beta1K <-  out$ci$betac_ci_unc
        beta1K0 <-  out$point$betac_pt_unc

      }

    }else{
      beta1K = c(NA,NA)
      beta1K0 =c(NA,NA)
    }

    output[[paste0(method,"beta1")]] <- beta1K
    output[[paste0(method,"beta1_pt")]] <- beta1K0
  }

  return(output)

}
