Bempdegreedistrib <- function(sam.out, num.sam, n.boot, idname = "Temp") {
      if (sam.out$n.neigh == 0) {
            # only information from the seeds
            res <- Bempdegreedistrib0(sam.out, num.sam, n.boot, idname = "Temp")
      } else {
            res <- BempdegreedistribK(sam.out, num.sam, n.boot, idname = "Temp")
      }
      res
}
