#' Build Bootstrap Confidence Intervals from Estimators of Network Degree
#' Distribution.
#'
#' The function will build bootstrap confidence intervals for the bootstrap
#' estimate of \eqn{\mu}.
#'
#' @inheritParams boot_parameters_est
#' @param bootstrap_mean A Boolean option to return the bootstrap confidence
#'  interval for the mean.
#' @param alpha \eqn{\alpha}-level of desired significance for statistical
#' inference. \eqn{\alpha < 0.5}.
#' @references
#' \insertRef{efron_bootstrap_1979}{snowboot}
#' @references
#' \insertRef{gel_bootstrap_2017}{snowboot}
#' @return A list of two elements
#'  \item{p_k_CI}{list of length \code{length(boot_deg_out$num_lsmi)}, one
#'    element per LSMI. Each element contains bootstrap confidence
#'    intervals for \eqn{\hat{p}_k^*} built from the provided bootstrap
#'    sample of the degree distribution (supplied by the output of
#'    the \code{\link{boot_deg}} function).}
#'  \item{mean_CI}{list of length \code{length(boot_deg_out$num_lsmi)}, one
#'    element per LSMI. Each element contains a bootstrap confidence
#'    intervals for \eqn{\hat{\mu}} built from the provided bootstrap sample of
#'    the degree distribution (supplied by the output of
#'    the \code{\link{boot_deg}} function).}
#' @export
#' @examples
#' net <- artificial_networks[[1]]
#' sam.out <- empd_deg_lsmi(
#'   net = net,
#'   n.seeds = 40,
#'   n.neigh = 1,
#'   num_lsmi = 1)
#' boot_deg_out <- boot_deg(sam.out = sam.out, boot_rep = 50)
#' a <- boot_ci(boot_deg_out)


boot_ci <- function(boot_deg_out, bootstrap_mean = TRUE, alpha = 0.05){
  # inception apply
  # take each "list" in empd, apply CI function to each distribution in "list".
  lower_bound <-  alpha / 2
  upper_bound <- (1 - (alpha / 2))
  p_k_CI <- lapply(
    boot_deg_out$empd,
    FUN <- function(x) {
      lapply(
        x,
        FUN <- function(dat) apply(
          dat, 2, stats::quantile,c(lower_bound, upper_bound)
          )
        )})
  mean_CI <- lapply(
    boot_deg_out$empd,
    FUN <- function(x) {
      lapply(
        x,
        FUN <- function(dat) {
          stats::quantile(
            boot_means(dat),
            c(lower_bound, upper_bound))
                             }
        )
      }
    )

  if (bootstrap_mean){
    list(p_k_CI = p_k_CI, mean_CI = mean_CI)
  } else list(p_k_CI = p_k_CI)
}
