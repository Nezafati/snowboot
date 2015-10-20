#' Bootstraping Empirical Degree Distribution
#'
#' This function delivers a boostrap estimate of network degree distribution
#' based on a LSMI sample. Default is one bootstrap replication.
#'
#' @param sam.out a list that is the output of \code{\link{Oempdegreedistrib}}.
#' @param num.sam a vector of integers containing the numeric ids of the LSMI
#'    samples when \code{sam.out$num.sam} is greater than one. When it is an
#'    integer, N, , all LSMI from 1 to N are taken from the input \code{sam.out}.
#' @param n.boot A positive integer number, the number of bootstrap replications.
#' @return A list consisting of:
#'    \item{values}{a list of length \code{num.sam} where each element is a
#'          vector containing the unique degree values sampled in each LSMI.}
#'    \item{empd}{a list of length three where each element contains a different
#'          estimate of degree distribution: empd.w.p0s - weighted bootstrap
#'          with a proportion of isolated nodes p0 being estimated by simple
#'          random sampling of bootstrapped seeds; empd.nw.p0sEkb - non-weighted
#'          bootstrap with a proportion of isolated nodes p0 being estimated by
#'          simple random sampling of bootstrapped seeds; empd.nw.p0sEks -ignore
#'          (see Thompson et al. for details)}
#'    \item{num.sam}{the same object as input argument \code{num.sam}}
#'    \item{n.boot}{the same object as input argument \code{n.boot}.}
#'    \item{n.neigh}{the number of waves carried out by the snowball sample.
#'          This is the same value from \code{sam.out$n.neigh}.}
#' @export
Bempdegreedistrib <- function(sam.out, num.sam = sam.out$num.sam, n.boot = 1) {
      if (sam.out$n.neigh == 0) {
            # only information from the seeds
            res <- Bempdegreedistrib0(sam.out, num.sam, n.boot)
      } else {
            res <- BempdegreedistribK(sam.out, num.sam, n.boot)
      }
      res
}
