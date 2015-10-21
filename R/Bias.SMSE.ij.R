Bias.SMSE.ij <- function(net, real.par, n.seeds, n.neigh, sam.size) {
      All.biasSMSE <- 1
      for (i in n.seeds) {
            for (j in n.neigh) {
                  # cat('i= ',i,'\t','j= ',j,'\n')
                  Obs.distrib <- Oempdegreedistrib(net, n.seeds = i, n.neigh = j, num.sam = sam.size)
                  Oparam <- OparametersEst(Obs.distrib)
                  biasSMSE <- Bias(Oparam, real.par)
                  All.biasSMSE <- cbind(All.biasSMSE, rbind(biasSMSE, c(i, j)))
            }
      }
      All.biasSMSE
}
