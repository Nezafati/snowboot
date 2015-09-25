sampleneighSequential <- function(net, n.seeds = 10, n.neigh = 1, seed = NULL) {
      # this function returns the vertices samples by snowball up to wave n.neigh. net is object network (what is important is
      # the component $edges and the length of $degree) this function randomly sample n.seeds and then select the neighbours up
      # to wave n.neigh also give the index of those nodes last added and for which in the next stages of the sampling we
      # assume we do not have their complete degree information. seed0 are the original seed sampleN are the possibly repeated
      # elements in the sample unodes are the no repeated elements in the sample nodes.waves are the vertices added in each
      # wave. Vertices may be present in more that one wave and more that once in a single wave. last.added are the vertices
      # that are the most recently added into the set.


      unodes <- nodes.waves <- as.list(rep(0, n.seeds))
      # Seed selection: is without replacement and at random
      if (is.null(seed)) {
            seed0 <- sort(sample(1:length(net$degree), n.seeds, rep = FALSE))
      } else {
            seed0 <- seed
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
# a<-sampleneighSequential(net,n.seeds=3,n.neigh=3,seed=NULL) a<-sampleneigh(net,n.seeds=3,n.neigh=3,seed=NULL)
