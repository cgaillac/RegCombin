#' Function to compute the main statistic for the point estimate
#'
#' @param x_kp0 a matrix containing the directions to compute the radial function, and the associated choice epsilon(q).
#' @param Xp the observations of the noncommon regressor (possibly conditional on Xc).
#' @param Yp the observations of the outcome variable.
#' @param for_critY the numerator of the ratio R for the point estimate of the radial function, on the grid grid_I;
#' @param dimXnc the dimension of the noncommon regressors
#' @param weights_xp the sampling or bootstrap weights for the dataset (Xnc,Xc).
#' @param weights_yp the sampling or bootstrap weights for the dataset (Y,Xc).
#' @param T_xy1 the apparent sample size taking into account the difference in the two datasets.
#' @param T_xy the apparent sample size, possibly conditional on Xc, the taking into account the difference in the two datasets.
#' @param es the numerical bootstrap parameter delta_n.
#' @param version version of the computation of the ratio, "first" is a degraded version but fast; "second" is a correct version but slower. Default is "second".
#' @param trunc equal to 2, for the definition of epsilon.
#' @param grid_I the grid of alpha on which we evaluate the ratio R to compute the point estimate of the radial function.
#' @param winsor indicates if winsorisation. Default is FALSE.
#'
#' @return
#'  the value of the point estimate of the radial function using the DGM method.
#'
#' @export
#'
#' @examples
compute_ratio_point <- function(x_kp0,Xp,Yp,for_critY,dimXnc,
                                  weights_xp,weights_yp,T_xy1,T_xy,es,
                                  version = "first", trunc=2, grid_I = NULL, winsor=FALSE){
  nb =  x_kp0[c(1)]
  q = x_kp0[-c(1,2)]
  eps = x_kp0[c(2)]

  if(dimXnc==1){
    XX = Xp%*%q
  }else{
    XX =as.matrix(Xp,dim(Xp)[1],dimXnc)%*%matrix(q,dimXnc,1)
  }

  #### to add winsorizing
  if(winsor ==TRUE){
    qu=quantile(XX, max(0.9, 1-3*log(length(XX))/length(XX)))
    XX[XX> qu] <- qu
  }

  if(version=="second"){
    indexbarre=sum(XX*weights_xp);
    Xs0 = cbind(weights_xp,XX-  indexbarre)
    Xs1 = Xs0[order(Xs0[,2], decreasing = T),]
    indexes0 = Xs1[,1]>0
    index = Xs1[  indexes0,2]
    weights_xp0 = Xs1[  indexes0,1]
    weights_xp0 = weights_xp0/sum(weights_xp0)
    for_critind=cumsum(weights_xp0*index);
    denX = approx( c(1,1- cumsum(weights_xp0)), c(for_critind,0) , xout = grid_I  )$y

    ### compute the ratio
    ratio =  for_critY/denX
    ratio[is.na(  ratio) ] = 10^6
    ratio_fun =  approxfun(grid_I,   ratio)
    kp00 = max(eps, trunc/T_xy)
    select0 =  (grid_I>=kp00) & (grid_I<=(1-kp00))
    if(sum(select0)==0){
      select0[floor(length(select0)/2)]=TRUE
    }
  }else{

    indexbarre=sum(XX*weights_xp);
    Xs0 = cbind(weights_xp,XX-  indexbarre)
    Xs0 = Xs0[order(Xs0[,2], decreasing = T),]
    weights_xp0 = Xs0[,1]
    index = Xs0[,2]

    for_critind=cumsum( weights_xp0*index);
    ratio = for_critY/for_critind
    kp00 = max(ceil(eps*length(ratio)), trunc)
    kp11 = length(ratio)- kp00

    ratio[is.na(  ratio)] = Inf
    ratio_fun = ratio

    select0 = kp00:kp11
  }
  ratio = ratio[select0]


  lambda=  pmax(min(ratio),0)

  output = vector("list")
  output[["lambda"]] <- lambda
  output[["ratio"]] <-  ratio_fun

  return(output)
}

