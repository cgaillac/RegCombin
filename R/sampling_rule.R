#' Subsampling rule
#'
#' @param n sample size.
#'
#' @return
#' the subsampling size
#'
#' @export
#'
#' @examples
sampling_rule <- function(n){return(0.25*n - 0.05*max(n-800,0)  -
                                      0.05*max(n-1200,0) - 0.05*max(n-1800,0)  -
                                      0.05*max(n-2000,0) - 0.05*(1-log(2400)/log(n))*max(n-2400,0) )} ## en cours
