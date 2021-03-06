#' Confidence Intervals from Bootstrapped Network Degree Distribution
#'
#' The function calculates bootstrap confidence intervals for the parameters
#' of network degree distribution: probabilities of node degrees \eqn{f(k)}
#' and mean degree \eqn{\mu}, where \eqn{k = 0, 1, \ldots} are the degrees.
#'
#' @details Currently, the bootstrap intervals can be calculated with two alternative
#' methods: \code{"percentile"} or \code{"basic"}. The \code{"percentile"}
#' intervals correspond to Efron's \eqn{100\cdot}\code{prob}\% intervals
#' \insertCite{@see @efron_1979, also Equation 5.18 by @davison_hinkley_1997 and Equation 3 by @gel_etal_2017, @chen_etal_2018_snowboot}{snowboot}:
#' \deqn{(\theta^*_{[B\alpha/2]}, \theta^*_{[B(1-\alpha/2)]}),}
#' where \eqn{\theta^*_{[B\alpha/2]}} and \eqn{\theta^*_{[B(1-\alpha/2)]}}
#' are empirical quantiles of the bootstrap distribution with \code{B} bootstrap
#' replications for parameter \eqn{\theta}
#' (\eqn{\theta} can be the \eqn{f(k)} or \eqn{\mu}),
#' and \eqn{\alpha = 1 -} \code{prob}.
#'
#' The \code{"basic"} method produces intervals
#' \insertCite{@see Equation 5.6 by @davison_hinkley_1997}{snowboot}:
#' \deqn{(2\hat{\theta} - \theta^*_{[B(1-\alpha/2)]}, 2\hat{\theta} - \theta^*_{[B\alpha/2]}),}
#' where \eqn{\hat{\theta}} is the sample estimate of the parameter.
#' Note that this method can lead to negative confidence bounds, especially
#' when \eqn{\hat{\theta}} is close to 0.
#'
#' @param x a list with bootstrapped results -- output of \code{\link{boot_dd}}.
#' @param prob confidence level for the intervals. Default is 0.95
#' (i.e., 95\% confidence).
#' @param method method for calculating the bootstrap intervals. Default is
#' \code{"percentile"} (see Details).
#'
#' @return A list object of class "\code{snowboot}" with the following elements:
#' \item{fk_ci}{A matrix of dimensions \eqn{2 \times}\code{length(x$fk)}, where
#' the number of columns corresponds to the number of probabilities \eqn{f(k)}
#' estimated from an LSMI sample. Each column of the matrix is a confidence
#' interval for a corresponding \eqn{f(k)}. I.e., the first row of the matrix
#' gives the lower bounds, while the second row contains all upper bounds.}
#' \item{mu_ci}{A numeric vector of length 2 with lower and upper confidence
#' bounds for the network mean degree \eqn{\mu}.}
#' \item{prob}{Confidence level for the intervals.}
#' \item{method}{Method that was used for calculating the bootstrap intervals.}
#' \item{fk}{A vector with an estimate of the degree distribution, copied
#'    from the input \code{x$fk}.}
#' \item{mu}{An estimate of the mean degree, copied from the input \code{x$mu}.}
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{boot_dd}}
#' @export
#' @examples
#' net <- artificial_networks[[1]]
#' lsmiEstimate <- lsmi_dd(net = net, n.seed = 5, n.wave = 3)
#' bootEstimates <- boot_dd(lsmiEstimate, B = 10)
#' bootIntervals1 <- boot_ci(bootEstimates)
#'
#' #Another version of the intervals:
#' bootIntervals2 <- boot_ci(bootEstimates, method = "basic")
#'
boot_ci <- function(x, prob = 0.95, method = c("percentile", "basic")) {
  method <- match.arg(method)
  alpha <- 1 - prob
  fk_ci <- apply(x$fkb, 1, quantile, probs = c(alpha/2, 1 - alpha/2))
  mu_ci <- quantile(x$mub, probs = c(alpha/2, 1 - alpha/2))
  if (method == "basic") {
    tmp <- fk_ci
    fk_ci[1,] <- 2 * x$fk - tmp[2,]
    fk_ci[2,] <- 2 * x$fk - tmp[1,]
    tmp <- mu_ci
    mu_ci[1] <- 2 * x$mu - tmp[2]
    mu_ci[2] <- 2 * x$mu - tmp[1]
  }
  res  <- list(fk_ci = fk_ci, mu_ci = mu_ci, prob = prob, method = method, fk = x$fk, mu = x$mu)
  class(res) <- "snowboot"
  return(res)
}
