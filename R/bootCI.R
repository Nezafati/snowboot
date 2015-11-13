#' Build Bootstrap Confidence Intervals for \eqn{\hat{p}_k^*}
#'
#' The function will build bootstrap confidence intervals for the bootstrap
#' estimate of  and \eqn{\mu} with a lower-bound of
#' \code{0.025} and an upper-bound of \code{0.975}.
#'
#' @inheritParams BparametersEst
#' @return A list of two elements
#'  \itme{p_k_CI}{This a list of length \code{length(outBootdeg$num.sam)}, one
#'    element per LSMI. Each element is contains three sets of bootstrap confidence
#'    intervals for \eqn{\hat{p}_k^*} corresponding to the three estimation
#'    methods. See \code{\link{bootdeg} for more on the three estimation methods.}
#'  \item{mean_CI}{This a list of length \code{length(outBootdeg$num.sam)}, one
#'    element per LSMI. Each element is contains three sets of bootstrap confidence
#'    intervals for \eqn{\hat{\mu}} corresponding to the three estimation
#'    methods. See \code{\link{bootdeg} for more on the three estimation methods.}
bootCI <- function(outBootdeg, bootstrap_mean=T){
  # inception apply
  # take each "list" in empd, apply CI function to each distribution in "list".
  p_k_CI <-   lapply(outBootdeg$empd,
                     FUN <- function(x) {
                       lapply(x,
                              FUN <- function(df)
                                apply(df, 2, quantile, c(0.025, 0.975)))})
  mean_CI <- lapply(outBootdeg$empd,
                    FUN <- function(x) {
                      lapply(x,
                             FUN <- function(df) {
                               quantile(bootmeans_from_bootdegdistrib(df),
                                        c(0.025, 0.975))
                               })
                      })


  if(bootstrap_mean){
    list(p_k_CI=p_k_CI, mean_CI=mean_CI)
  } else list(p_k_CI=p_k_CI)
}
