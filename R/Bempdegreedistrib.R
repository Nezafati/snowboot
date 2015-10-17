#' Obtaining an Empirical Network Degree Distribution by Bootstrap Sampling from
#' a Labeled Snowball Sample with Multiple Inclusion (LSMI).
#'
#' The function will boostrap the a LSMI to obtain an empirical degree
#' distribution.
#'
#' @param sam.out a list that is the output of \code{\link{Oempdegreedistrib}}.
#' @param num.sam a vector of integers containing the number of
#'    different samples taken from the same network.
#' @param n.boot \code{integer.} The size of the bootstrap sample.
#' @return A list consisting of:
#'    \item{values}{a list of vectors containing unique
#'          sampled degree values. Each sample has its own vector. The output
#'          can be used to identify which degrees are present in the sample.}
#'    \item{empd}{a list of length three}
#'    \item{num.sam}{the number of LSMI from the \code{sam.out} object}
#'    \item{n.boot}{\code{integer.} The size of the bootstrap sample.}
#'    \item{n.neigh}{the number of waves carried out by the snowball sample.
#'          This is the same value from \code{sam.out$n.neigh}.}
#' @export
Bempdegreedistrib <- function(sam.out, num.sam, n.boot) {
      if (sam.out$n.neigh == 0) {
            # only information from the seeds
            res <- Bempdegreedistrib0(sam.out, num.sam, n.boot)
      } else {
            res <- BempdegreedistribK(sam.out, num.sam, n.boot)
      }
      res
}
