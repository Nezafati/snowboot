#' A function to snowball sample the neighbors around a single seed.
#'
#' This function returns the vertices samples by snowball up to wave n.neigh
#' around a single seed. The functions main purpose is to be called from
#' \code{\link{sampleneighSequential}.}
#' @param net A network object.
#' @param seed0 A number that is the id of the seed whose neighbors are sampled.
#' @param n.neigh A number for the waves of snowball sampling.
#' @return a list containing:
#'    \item{seed}{A number that is the id of the seed whose neighbors are sampled.}
#'    \item{sampleN}{A numeric vector containg ids of the vertices from
#'          the snowball sampling and the intial seed's id. This vector may have
#'          duplicates, since the algorithm allows for multiple inclusions.}
#'    \item{unode}{A numeric vector containg the unique values in \code{$sampleN}.}
#'    \item{nodes.waves}{A list of numeric vectors (one per wave) each
#'          containing the ids of the nodes sampled in that particular wave.}
#' @export
sampleneighAroundOneSeed <- function(net, seed0, n.neigh = 1) {
      # this function returns the vertices samples by snowball up to wave n.neigh around a single seed net is object network
      # (what is important is the component $edges and the length of $degree) this function randomly sample n.seeds and then
      # select the neighbours up to wave n.neigh seed0 the id of the seed sampleN are the possibly repeated elements in the
      # sample unodes are the no repeated elements in the sample nodes.waves are the vertices added in each wave. Vertices may
      # be present in more that one wave and more that once in a single wave. last.added are the vertices that are the most
      # recently added into the set.

      sampleN <- nodes <- seed0
      nodes.waves <- as.list(rep(0, n.neigh))
      effEdges <- net$edges
      more <- TRUE
      nn <- n.neigh
      new.nodes <- 0
      # if(n.neigh==0) we only keep the seeds

      wave <- 1
      while (wave <= n.neigh & more) {
            a <- is.element(effEdges, nodes)  #'nodes' will be accumulating all included vertices (non repeated)
            if (any(a))
            {
                  eedges <- which(matrix(a, dim(effEdges)[1], 2), arr.ind = TRUE)  #now it is the row number and column where they are in the edges matrix
                  nodes.waves[[wave]] <- arr.nodes <- sort(effEdges[cbind(eedges[, 1], sapply(eedges[, 2], FUN = switch, 2,
                                                                                              1))])  #the vertices we arrived to (duplicity is allowed)
                  # I need this specially to know which vertices were the last added:
                  if (!anyDuplicated(eedges[, 1])) {
                        new.nodes <- arr.nodes  #all the vertices we arrive to weren't already included in 'nodes'
                  } else {
                        new.nodes <- setdiff(arr.nodes, nodes)
                  }  #Then, already included vertices are not considered new because of inclusion of edge connecting them
                  ### subEdges<-effEdges[unique(a),] #the subset of edges. The repeated just have to be included once. maybe we are arriving
                  ### to the nodes more than once (due to small cycles) or we can get again to already included vertices (due to larger
                  ### cycles).  We want to include them as may times as they are neighbours of already included vertices. That is why I
                  ### consider arr.nodes.if a originally seed vertex is included more than once, it is because it was selected also by
                  ### following one edge and then it also has the category of non seed.

                  sampleN <- sort(c(sampleN, arr.nodes))
                  nodes <- unique(sampleN)
                  if (nn > 1)
                        effEdges <- effEdges[-unique(eedges[, 1]), ]  #I remove the 'used edges' to facilitate following searches within while,and assure
                  # we do not 'arrive' to a node more times than edges it has.
                  if (length(effEdges) > 0) {
                        if (is.vector(effEdges))
                              effEdges <- t(effEdges)  #I have to this very often in R to make sure effEdges is a matrix and not a vector
                  } else {
                        more <- FALSE
                  }  #when it reduces to become a matrix with one row.
            }  #end if(any(a))
            wave <- wave + 1
      }  #end while
      # browser() list(seed=seed0,sampleN=sampleN,unodes=nodes,nodes.waves=nodes.waves,last.added=sort(new.nodes))
      list(seed = seed0, sampleN = sampleN, unodes = nodes, nodes.waves = nodes.waves)
}
