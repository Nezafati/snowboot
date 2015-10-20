#' Snowball sampling with multiple inclusion.
#'
#' The function will conduct snowball sampling.
#'
#' @param net a list that must contain elements
#'    \code{$n} (\code{integer}. network order),
#'    \code{$edges} (\code{matrix}. a \code{n}x\code{2} matrix),
#'    and \code{$degree} (\code{integer} vector of length n).
#'    The object can be created by \code{\link{local.network.MR.new5}} or
#'    it can be imported.
#' @param n.seeds a number of seeds in the snowball sample.
#'    It must be a positive integer.
#' @param n.neigh a number of waves to be sampled around each seed in LSMI.
#'    For example, n.neigh = 0 corresponds to seeds only, and n.neigh = 1
#'    corresponds to sampling seeds and their first neighbors).
#'    Note that the algorithm allow for mutiple inclusions.
#' @param seeds A vector of length \code{n.seeds} containing the
#'    numeric ids of the seeds to initiate sampling. Note that this is an
#'    optional parameter.
#' @return A list containing the following elements:
#'    \item{seeds}{A \code{numeric} a vector containing the numerics ids of
#'          sampled seeds.}
#'    \item{sampleN}{A \code{numeric} vector containg ids of the nodes from
#'          the snowball sampling and the intial seeds' ids. This vector may have
#'          duplicates, since the algorithm allows for multiple inclusions.}
#'    \item{unodes}{a list of length \code{n.seeds} where each element is a
#'          \code{numeric} vector containing the seed's id and
#'          the unique ids of all nodes that were snowball sampled from
#'          that seed using \code{\link{sampleneighAroundOneSeed}}
#'          (one vector per seed).}
#'    \item{nodes.waves}{A list of length \code{n.seeds} where each element is
#'          a list of length \code{n.neigh} (Note: these lists are the output
#'          object \code{$nodes.waves} from
#'          \code{\link{sampleneighAroundOneSeed}}) that contains vectors of
#'          numeric id's of the nodes reached in each respective wave from the
#'          respective seed.}
#' @export
sampleneighSequential <- function(net, n.seeds = 10, n.neigh = 1, seeds = NULL) {

      unodes <- nodes.waves <- as.list(rep(0, n.seeds))
      # Seed selection: is without replacement and at random
      if (is.null(seeds)) {
            seed0 <- sort(sample(1:length(net$degree), n.seeds, rep = FALSE))
      } else {
            seed0 <- seeds
      }
      sampleN <- NULL
      for (i in 1:n.seeds) {
            res <- sampleneighAroundOneSeed(net, seed0[i], n.neigh)
            sampleN <- c(sampleN, res$sampleN)
            unodes[[i]] <- res$unodes
            nodes.waves[[i]] <- res$nodes.waves
      }
      list(seeds = seed0, sampleN = sort(sampleN), unodes = unodes, nodes.waves = nodes.waves)
}
# Examples #we are not really interested in running this function directly but within the next function called empdegree
# distrib6 net<-local.network.MR.new5(n=100,distrib='pois',param=2)
# a<-sampleneighSequential(net,n.seeds=3,n.neigh=3,seed=NULL) a<-sampleneigh(net,n.seeds=3,n.neigh=3,seeds=NULL)
