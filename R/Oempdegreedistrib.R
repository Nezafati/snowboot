Oempdegreedistrib <- function(net, n.seed, n.neigh, num.sam, idname = "Temp", seeds) {
      if (n.neigh == 0) {
            # only information from the seeds
            res <- Oempdegreedistrib0(net, n.seed, n.neigh, num.sam, idname = "Temp", seeds)
      } else {
            res <- OempdegreedistribK(net, n.seed, n.neigh, num.sam, idname = "Temp", seeds)
      }
      res
}
