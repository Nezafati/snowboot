#' A function to obtain the empirical distribution of the snowball samples.
#'
#' This function obtains the empirical degree distribution from the snowball
#' samples. The user inputs the network, the number of seeds, waves and snowball
#' samples to take. The return value is a list which includes the empirical
#' distribution of the snowball samples \code{Oempd}. Usually one snowball
#' sample is requested (\code{num.sam = 1}) n.seed is
#' the number of seed to set the neighbourhood sample n.neigh is the neighbouhood size around each seed num.sam is the
#' number of different samples taken from the same network idname is to identify from which nets we are sampling and
#' resampling.
#'
#' @param net the network that the function samples from.
#' @param n.seed the number of seeds in the snowball sample.
#' @param n.neigh the number of waves carried out by the snowball sample (e.g.
#'    n.neigh = 2 means each snowball sample will contain the seeds, the
#'    neighbors of the seeds, and the neighors of the neighbors of the seeds).
#'    Recall this algorithm allow mutiple inclusions.
#' @param num.sam a number value for the snowball samples to take.
#' @param idname a single \code{character} to name the object. The default
#'    value is "Temp".
#' @param seeds an optional paramater that is used whenever we know the seeds
#' we would like to initiate the snowball with.
#'
#' @return a list consisting of
#'    \item{idname}{A \code{character} value for the optional id of the object.}
#'    \item{samples}{A list containing lists (one per sample) each with
#'          three tables: the frequency of degrees sampled from seeds,
#'          non-seeds including duplicated nodes,
#'          and non-seeds without duplications.}
#'    \item{values}{A list containing vectors (one per sample) of the unique
#'          degree values sampled.}
#'    \item{Oemp}{The main output. A list of lists (one per sample) each
#'          containing two tables for estimating the empirical distribution
#'          from the sample.
#'          The first table, called \code{Oempd} is of interest}
#'    \item{num.sam}{a number value for the snowball samples to take}
#'    \item{val.seed}{A list containing vectors (one per snowball sample)
#'          of the unique degree values sampled from the seeds.}
#'    \item{val.nonseed}{A list containing vectors (one per snowball sample)
#'          of the unique degree values sampled from the non-seeds.
#'          *Note: This item does not exist when n.neigh equals zero.}
#'    \item{n.seed}{the number of seeds in the snowball sample.}
#'    \item{n.neigh}{the number of waves carried out by the snowball sample. See
#'          input argument for details.}
#'    \item{p0.real}{the proportion of nodes from the network that are
#'          degree zero.
#'          *Note: This item does not exist when n.neigh equals zero.}
#'    \item{p0.seed}{a list of number values (one per snowball sample) equal to
#'          the proportion of seeds with degree zero.}
#'    \item{ekseed}{a list of numeric values (one per snowball sample) equal
#'          to the sample mean of the seeds.
#'          *Note: This item does not exist when n.neigh equals zero.}
#'    \item{seeds1}{a matrix of the seed ids, where each row corresponds to
#'          one snowball sample. }
#' @export
Oempdegreedistrib <- function(net, n.seed, n.neigh, num.sam, idname = "Temp", seeds = NULL) {
      if (n.neigh == 0) {
            # only information from the seeds
            res <- Oempdegreedistrib0(net, n.seed, n.neigh, num.sam, idname = "Temp", seeds)
      } else {
            res <- OempdegreedistribK(net, n.seed, n.neigh, num.sam, idname = "Temp", seeds)
      }
      res
}
