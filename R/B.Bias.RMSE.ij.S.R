B.Bias.RMSE.ij.S <- function(net, n.seeds, n.neigh, sam.size, n.boot, otherNetParameters = FALSE) {
      # sam.size is the number of different samples taken from the network for each i and j otherNetParameters is true if
      # intervals and fallins for the rest of the parmeters (other than mean) are required.
      seeds2 <- array(0, dim = c(length(as.vector(n.neigh)), length(as.vector(n.seeds)), sam.size, max(n.seeds)))
      All.biasRMSE <- 1
      Mean.intervals.list <- Mean.fallins <- B.mean <- as.list(rep(NA, length(n.seeds) * length(n.neigh)))
      realparam <- real.parameters(net)
      OtherPar.intervals.list <- NULL
      if (otherNetParameters) {
            OtherPar.intervals.list <- OtherPar.fallins <- as.list(rep(NA, length(n.seeds) * length(n.neigh)))
      }
      b <- 1
      for (i in n.seeds) {
            for (j in n.neigh) {

                  if (j == 0) {
                        n.dist <- 1  #n.dist is the number of different emp distr.
                  } else {
                        n.dist <- 3
                  }

                  # cat('i= ',i,'\t','j= ',j,'\n') browser()
                  if (j == 0) {
                        Obs.distrib <- Oempdegreedistrib(net, n.seeds = i, n.neigh = j, num.sam = sam.size)
                        TMP <- Obs.distrib$seeds1
                  } else {
                        Obs.distrib <- Oempdegreedistrib(net, n.seeds = i, n.neigh = j, num.sam = sam.size, seeds = TMP)
                  }

                  Oparam <- OparametersEst(Obs.distrib)
                  # browser()
                  seeds2[which(n.neigh == j), which(n.seeds == i), , 1:dim(Obs.distrib$seeds1)[2]] <- Obs.distrib$seeds1

                  B.distrib <- Bempdegreedistrib(Obs.distrib, num.sam = sam.size, n.boot = n.boot)
                  Bparam <- BparametersEst(B.distrib)
                  # browser()
                  B.mean[[b]] <- Bparam$mean


                  biasRMSE <- B.Bias(Bparam, Oparam)
                  All.biasRMSE <- cbind(All.biasRMSE, biasRMSE)

                  # mean intervals and fallins

                  Mean.intervals.list[[b]] <- intervals <- BMean.intervals(Opar = Oparam$mean, Bpar = Bparam$mean, n.dist)
                  Mean.fallins[[b]] <- list(ij = c(i, j), fallins = fallInsMean(ints = intervals, rpar = realparam$rmean, n.dist))
                  if (otherNetParameters)
                  {
                        OtherPar.intervals.list[[b]] <- OtherPar.intervals <- list(quartiles = Bintervalsi(Opari = Oparam$quart,
                                                                                                           Bpari = Bparam$quart, n.dist, num.par = 3), rfreq = Bintervalsi(Opari = Oparam$rfreq, Bpari = Bparam$rfreq,
                                                                                                                                                                           n.dist, num.par = 5), deciles = Bintervalsi(Opari = Oparam$deciles, Bpari = Bparam$deciles, n.dist, num.par = 9))
                        Other.fallins[[b]] <- list(ij = c(i, j), fi.quart = NA, fi.rfreq = NA, fi.deciles = NA)
                  }  #if(OtherNetParameters)

                  b <- b + 1
            }
      }
      list(All.biasRMSE = All.biasRMSE, Mean.intervals = Mean.intervals.list, Mean.fallins = Mean.fallins, B.mean = B.mean,
           seeds2 = seeds2)
}


