#' Summary of a Network Degree Sequence
#'
#' This function provides summary statistics of a network degree
#' distribution.
#' @inheritParams Oempdegreedistrib
#' @return a list consisting of:
#'    \item{realdd}{a vector of length \code{net$n} where each element
#'    corresponds to the degree of a node in the network.}
#'    \item{rmean}{a numeric vector of length one that is the atrithmetic mean
#'          of \code{realdd}. see (\code{\link[base]{mean}}).}
#'    \item{rquart}{the lower, median, and upper quartiles of degree sequence
#'          \code{realdd}.}
#'    \item{rfreq}{a numeric vector of length five where each element
#'          corresponds to the proportion of seeds with degree: 0, 1, 2, 3, and
#'          4 (respectively).}
#'    \item{rdeci}{a numeric vector of lenth nine containing the deciles of
#'          of degree sequence \code{realdd} in increasing order.}
#' @export
real.parameters <- function(net) {
      # this function obtains the real parameters in a network
      realdd <- net$degree - net$degree.left
      rmean <- mean(realdd)
      rquart <- quantile(realdd, prob = c(0.25, 0.5, 0.75))
      rfreq <- c(sum(realdd == 0), sum(realdd == 1), sum(realdd == 2), sum(realdd == 3), sum(realdd == 4))/length(net$degree)
      rdeci <- quantile(realdd, prob = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
      list(realdd = realdd, rmean = rmean, rquart = rquart, rfreq = rfreq, rdeci = rdeci)
}
