#' Obtaining Empirical Network Degree Distribution from Labeled Snowball
#' Sampling with Multiple Inclusion (LSMI).
#'
#' Oempdegreedistrib is used to obtain the empirical network degree distribution
#' from labeled snowball sampling with multiple inclusion (LSMI). # move to arg# The user inputs the network, the number of seeds, waves and snowball
#' samples to take. The return value is a list which includes the empirical
#' distribution of the snowball samples \code{Oempd}.  n.seed is
#' the number of seed to set the neighbourhood sample n.neigh is the neighbouhood size around each seed num.sam is the
#' number of different samples taken from the same network idname is to identify from which nets we are sampling and
#' resampling.
#'
#' @param net a list that must contain elements $edges (\code{matrix}. a two column matrix),
#'    $n (\code{num}. network order). The object can be created by
#'    \code{\link{local.network.MR.new5}} or it can be imported.
#' @param n.seed a number of seeds in the snowball sample.
#'    It must be a positive integer.
#' @param n.neigh a number of waves to be sampled around each seed in LSMI. (e.g.
#'    n.neigh = 0 corresponds to seeds only, n.neigh = 1 corresponds to sampling
#'    seeds and their first neighbors).
#'    Recall this algorithm allow for mutiple inclusions #reference to paper.
#' @param num.sam a number for the snowball samples to take.
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
