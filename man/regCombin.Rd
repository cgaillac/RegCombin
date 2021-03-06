% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regCombin.R
\name{regCombin}
\alias{regCombin}
\title{Function computing all the different bounds : DGM and/or Variance}
\usage{
regCombin(
  Ldata,
  Rdata,
  out_var,
  nc_var,
  c_var = NULL,
  nc_sign = NULL,
  c_sign = NULL,
  nbCores = 1,
  methods = c("DGM"),
  kp = 0.5,
  nb_pts = 1,
  Opt = "boot",
  Bsamp = 1000,
  data_k = 30,
  C0 = 0.5,
  list_ex = c(),
  weights_x = NULL,
  weights_y = NULL,
  trunc = 2,
  boot_par = 0.3,
  meth = "adapt",
  modeNA = FALSE,
  winsor = FALSE,
  version = "second",
  projections = FALSE
)
}
\arguments{
\item{Ldata}{dataset containing (Y,Xc) where Y is the outcome, Xc are potential common regressors.}

\item{Rdata}{dataset containing (Xnc,Xc) where Xnc are the non commonly observed regressors, Xc are potential common regressors.}

\item{out_var}{label of the outcome variable Y.}

\item{nc_var}{label of the non commonly observed regressors Xnc.}

\item{c_var}{label of the commonly observed regressors Xc.}

\item{nc_sign}{sign restrictions on the non-commonly observed regressors Xnc: -1 for a minus sign, 1 for a plus sign, 0 otherwise. Default is NULL, i.e. no constraints.}

\item{c_sign}{sign restrictions on the commonly observed regressors: -1 for a minus sign, 1 for a plus sign, 0 otherwise. Default is NULL, i.e. no constraints.}

\item{nbCores}{number of cores for the parallel computation. Default is 1.}

\item{methods}{method used for the bounds: "DGM" (Default) and/or "Variance".}

\item{kp}{If data_k =NULL, then epsilon is taken equal to kp.}

\item{nb_pts}{the constant C in DGM for the epsilon_0, the lower bound on the grid for epsilon, taken equal to nb_pts*ln(n)/n. Default is 1 without regressors Xc, 3 with Xc.}

\item{Opt}{the option of the method to compute the critical values, either "boot" (numerical bootstrap) or "subsampling". Default is numerical bootstrap.}

\item{Bsamp}{the number of bootstrap/subsampling replications. Default is 1000.}

\item{data_k}{the number of points for the grid search on epsilon. Default is 30. If NULL, then epsilon is taken fixed equal to kp.}

\item{C0}{the upper bound on the grid for epsilon. Default is 0.5.}

\item{list_ex}{the list of the data that should not be used on the different clusters. This reduces the computational cost. Default is empty list.}

\item{weights_x}{the sampling weights for the dataset (Xnc,Xc). Default is NULL.}

\item{weights_y}{the sampling weights for the dataset (Y,Xc). Default is NULL.}

\item{trunc}{equal to 2, for the definition of epsilon.}

\item{boot_par}{the numerical bootstrap parameter. Default is 0.3 without Xc, 0.35 with Xc.}

\item{meth}{the method for the choice of epsilon, either "adapt", i.e. adapted to the direction or "min" the minimum over the directions. Default is "adapt".}

\item{modeNA}{indicates if NA introduced if the interval is empty. Default is FALSE.}

\item{winsor}{indicates if winsorisation. Default is FALSE.}

\item{version}{version of the computation of the ratio, "first" is a degraded version but fast; "second" is a correct version but slower. Default is "second".}

\item{projections}{if FALSE compute the identified set along some directions or the confidence regions. Default is FALSE}
}
\value{
a list containing, in order:

 - DGM_complete or Variance_complete : the complete outputs of the functions DGM_bounds or Variance_bounds.

 and additional pre-treated outputs, replace below "method" by either "DGM" or "Variance":

- methodCI:  the confidence region on the betanc without sign constraints

- methodpt:  the bounds point estimates on the betanc without sign constraints

- methodCI_sign:  the confidence region on the betanc with sign constraints

- methodpt_sign:  the bounds point estimates on the betanc with sign constraints

- methodkp: the values of epsilon(q)

- methodbeta1: the confidence region on the betac corresponding to the common regressors Xc without sign constraints

- methodbeta1_pt: the bounds point estimates on the betac corresponding to the common regressors Xc without sign constraints

- methodbeta1_sign: the confidence region on the betac corresponding to the common regressors Xc with sign constraints

- methodbeta1_sign_pt: the bounds point estimates on the betac corresponding to the common regressors Xc with sign constraints
}
\description{
Function computing all the different bounds : DGM and/or Variance
}
\examples{

}
