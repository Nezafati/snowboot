#' Snowball sampling with multiple inclusion.
#'
#' The function will conduct snowball sampling.
#' @references
#' \insertRef{gel_bootstrap_2017}{snowboot}
#' @param classic Option for neighborhoods, i.e. waves, without multiple inclusions.
#' @inheritParams empd_deg_lsmi
#'
#' @return A list containing the following elements:
#'    \item{seeds}{A \code{numeric} a vector containing the numeric ids of
#'          sampled seeds.}
#'    \item{sampleN}{A \code{numeric} vector containing ids of the nodes from
#'          the snowball sampling and the initial seeds' ids. This vector may have
#'          duplicates, since the algorithm allows for multiple inclusions.}
#'    \item{unodes}{A list of length \code{n.seeds} where each element is a
#'          \code{numeric} vector containing the seed's id and
#'          the unique ids of all nodes that were snowball sampled from
#'          that seed using \code{\link{sample_about_one_seed}}
#'          (one vector per seed).}
#'    \item{nodes.waves}{A list of length \code{n.seeds} where each element is
#'          a list of length \code{n.neigh} (Note: these lists are the output
#'          object \code{$nodes.waves} from
#'          \code{\link{sample_about_one_seed}}) that contains vectors of
#'          numeric id's of the nodes reached in each respective wave from the
#'          respective seed.}
#' @export
#' @examples
#' net <- artificial_networks[[1]]
#' a <- lsmi(net, n.seeds = 20, n.neigh = 2)

lsmi <- function(net, n.seeds = 10, n.neigh = 1, seeds = NULL, classic = FALSE) {
  unodes <- nodes.waves <- as.list(rep(0, n.seeds))
  # Seed selection: is without replacement and at random
  if (is.null(seeds)) {
    seed0 <- sort(sample(1:length(net$degree), n.seeds, replace = FALSE))
  } else {
    seed0 <- seeds
  }
  sampleN <- NULL
  if(classic == TRUE){
    g <- igraph::graph_from_edgelist(net$edges)
    sampleN <- unlist(igraph::ego(g, order = n.neigh, nodes = seed0))
  } else {
    for (i in 1:n.seeds) {
      snowball <- sample_about_one_seed(net, seed0[i], n.neigh)
      sampleN <- c(sampleN, snowball$sampleN)
      unodes[[i]] <- snowball$unodes
      nodes.waves[[i]] <- snowball$nodes.waves
    }
  }
  res <- list(seeds = seed0, sampleN = sort(sampleN))
  res
}
