#' @title Build Bootstrap Confidence Intervals for \eqn{\hat{p_{k}}}
#'
#' The function will build bootstrap confidence intervals for the bootstrap
#' estimate of \eqn{\mu}.
#'
#' @inheritParams BparametersEst
#' @param bootstrap_mean A Boolean option to return the bootstrap confidence
#'  interval for the mean.
#' @param alpha \eqn{\alpha}-level of desired significance for statistical
#' inference. \eqn{\alpha < 0.5}.
#' @references Efron, B. (1992). Bootstrap methods: another look at the jackknife.
#' In Breakthroughs in statistics (pp. 569-593). Springer New York.
#' @references Gel, Y. R., Lyubchich, V., & Ramirez Ramirez, L. L. (2017).
#' Bootstrap quantification of estimation uncertainties in
#' network degree distributions. Scientific Reports, 7, 5807.
#' \url{http://doi.org/10.1038/s41598-017-05885-x}
#' @return A list of two elements
#'  \item{p_k_CI}{list of length \code{length(outBootdeg$num.sam)}, one
#'    element per LSMI. Each element contains bootstrap confidence
#'    intervals for \eqn{\hat{p}_k^*} built from the provided bootstrap  sample of
#'    the degree distribution (supplied by the output of
#'    the \code{\link{bootdeg}} function). See \code{\link{bootdeg}} for more
#'    on the available estimation methods.}
#'  \item{mean_CI}{list of length \code{length(outBootdeg$num.sam)}, one
#'    element per LSMI. Each element contains a bootstrap confidence
#'    intervals for \eqn{\hat{\mu}} built from the provided bootstrap sample of
#'    the degree distribution (supplied by the output of
#'    the \code{\link{bootdeg}} function). See \code{\link{bootdeg}} for more
#'    on the available estimation methods.}
#' @export
#' @examples
#' net <- artificial_networks[[1]]
#' sam.out <- Oempdegreedistrib(net = net, n.seeds = 40, n.neigh = 1, num.sam = 1)
#' outBootdeg <- bootdeg(sam.out = sam.out, n.boot = 50)
#' a <- bootCI(outBootdeg)


bootCI <- function(outBootdeg, bootstrap_mean=TRUE, alpha = 0.05){
  # inception apply
  # take each "list" in empd, apply CI function to each distribution in "list".
  lower_bound <-  alpha/2
  upper_bound <- (1 - (alpha/2))
  p_k_CI <-   lapply(outBootdeg$empd,
                     FUN <- function(x) {
                       lapply(x,
                              FUN <- function(df)
                                apply(df, 2, stats::quantile,
                                      c(lower_bound, upper_bound)))})
  mean_CI <- lapply(outBootdeg$empd,
                    FUN <- function(x) {
                      lapply(x,
                             FUN <- function(df) {
                               stats::quantile(bootmeans_from_bootdegdistrib(df),
                                        c(lower_bound, upper_bound))
                               })
                      })


  if(bootstrap_mean){
    list(p_k_CI=p_k_CI, mean_CI=mean_CI)
  } else list(p_k_CI=p_k_CI)
}
