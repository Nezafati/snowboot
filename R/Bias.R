Bias <- function(opar.out, rpar.out) {
      # opar.out the 'estimated parameters'. List with the names, mean, quetiles, rfreq, deciles rpar.out the 'real
      # parameters'. Listi the the names rmean, rquart, rfreq, rdeci
      meanBias <- mean(opar.out$mean - rpar.out$rmean)
      quarBias <- rowMeans(t(opar.out$quartiles) - rpar.out$rquart)
      rfreqBias <- rowMeans(t(opar.out$rfreq) - rpar.out$rfreq)
      deciBias <- rowMeans(t(opar.out$deciles) - rpar.out$rdeci)
      meanSMES <- sqrt(mean((opar.out$mean - rpar.out$rmean)^2))
      quarSMES <- sqrt(rowMeans((t(opar.out$quartiles) - rpar.out$rquart)^2))
      rfreqSMES <- sqrt(rowMeans((t(opar.out$rfreq) - rpar.out$rfreq)^2))
      deciSMES <- sqrt(rowMeans((t(opar.out$deciles) - rpar.out$rdeci)^2))

      rnam <- c("mean", "quart1", "quart2", "quart3", paste("rfreq", 0:4), paste("deci", 1:9))
      BiasSMES <- data.frame(meanBias = c(meanBias, quarBias, rfreqBias, deciBias), sqrtSMES = c(meanSMES, quarSMES, rfreqSMES,
                                                                                                 deciSMES))
      rownames(BiasSMES) <- rnam
      BiasSMES
}

# ---------------------------------------------------------------------------------------#

Bias.SMSE.ij <- function(net, real.par, n.seeds, n.neigh, sam.size) {
      All.biasSMSE <- 1
      for (i in n.seeds) {
            for (j in n.neigh) {
                  # cat('i= ',i,'\t','j= ',j,'\n')
                  Obs.distrib <- Oempdegreedistrib(net, n.seed = i, n.neigh = j, num.sam = sam.size)
                  Oparam <- OparametersEst(Obs.distrib)
                  biasSMSE <- Bias(Oparam, real.par)
                  All.biasSMSE <- cbind(All.biasSMSE, rbind(biasSMSE, c(i, j)))
            }
      }
      All.biasSMSE
}

# ---------------------------------------------------------------------------------------#

B.Bias <- function(bpar.out, opar.out) {
      # bpar.out the bootstrap estimated parameters. List with the names mean, quartiles, rfreq, deciles. This elements are
      # arrays opar.out the estimated parameters form the original samples. List with the names, mean, quetiles, rfreq, deciles

      if (length(bpar.out$num.sam) == 1) {
            sams <- 1:bpar.out$num.sam
      } else {
            sams <- bpar.out$num.sam
      }  #the subset of samples
      # apply(bpar.out$mean-opar.out$mean[num.sam],c(1,3),FUN=mean) #matrix of mean tn.sam X n.dist
      meanBias <- apply(bpar.out$mean - opar.out$mean[sams], 3, FUN = mean)  #the mean error by emp distribution
      if (bpar.out$n.dist > 1) {
            q1Bias <- apply(bpar.out$quartiles[, 1, , ] - opar.out$quartiles[sams, 1], 3, FUN = mean)
            q2Bias <- apply(bpar.out$quartiles[, 2, , ] - opar.out$quartiles[sams, 2], 3, FUN = mean)
            q3Bias <- apply(bpar.out$quartiles[, 3, , ] - opar.out$quartiles[sams, 3], 3, FUN = mean)
            ## quarBias<-apply(bpar.out$quartiles-opar.out$rquart[num.sam],c(2,4),FUN=mean)
            quarBias <- cbind(q1Bias, q2Bias, q3Bias)

            f0Bias <- apply(bpar.out$rfreq[, 1, , ] - opar.out$rfreq[sams, 1], 3, FUN = mean)
            f1Bias <- apply(bpar.out$rfreq[, 2, , ] - opar.out$rfreq[sams, 2], 3, FUN = mean)
            f2Bias <- apply(bpar.out$rfreq[, 3, , ] - opar.out$rfreq[sams, 3], 3, FUN = mean)
            f3Bias <- apply(bpar.out$rfreq[, 4, , ] - opar.out$rfreq[sams, 4], 3, FUN = mean)
            f4Bias <- apply(bpar.out$rfreq[, 5, , ] - opar.out$rfreq[sams, 5], 3, FUN = mean)
            rfreqBias <- cbind(f0Bias, f1Bias, f2Bias, f3Bias, f4Bias)

            d1Bias <- apply(bpar.out$deciles[, 1, , ] - opar.out$deciles[sams, 1], 3, FUN = mean)
            d2Bias <- apply(bpar.out$deciles[, 2, , ] - opar.out$deciles[sams, 2], 3, FUN = mean)
            d3Bias <- apply(bpar.out$deciles[, 3, , ] - opar.out$deciles[sams, 3], 3, FUN = mean)
            d4Bias <- apply(bpar.out$deciles[, 4, , ] - opar.out$deciles[sams, 4], 3, FUN = mean)
            d5Bias <- apply(bpar.out$deciles[, 5, , ] - opar.out$deciles[sams, 5], 3, FUN = mean)
            d6Bias <- apply(bpar.out$deciles[, 6, , ] - opar.out$deciles[sams, 6], 3, FUN = mean)
            d7Bias <- apply(bpar.out$deciles[, 7, , ] - opar.out$deciles[sams, 7], 3, FUN = mean)
            d8Bias <- apply(bpar.out$deciles[, 8, , ] - opar.out$deciles[sams, 8], 3, FUN = mean)
            d9Bias <- apply(bpar.out$deciles[, 9, , ] - opar.out$deciles[sams, 9], 3, FUN = mean)
            deciBias <- cbind(d1Bias, d2Bias, d3Bias, d4Bias, d5Bias, d6Bias, d7Bias, d8Bias, d9Bias)

            meanSMSE <- sqrt(apply((bpar.out$mean - opar.out$mean[sams])^2, 3, FUN = mean))

            q1SMSE <- apply((bpar.out$quartiles[, 1, , ] - opar.out$quartiles[sams, 1])^2, 3, FUN = mean)
            q2SMSE <- apply((bpar.out$quartiles[, 2, , ] - opar.out$quartiles[sams, 2])^2, 3, FUN = mean)
            q3SMSE <- apply((bpar.out$quartiles[, 3, , ] - opar.out$quartiles[sams, 3])^2, 3, FUN = mean)
            ## quarBias<-apply(bpar.out$quartiles-opar.out$rquart[num.sam],c(2,4),FUN=mean)
            quarSMSE <- sqrt(cbind(q1SMSE, q2SMSE, q3SMSE))

            f0SMSE <- apply((bpar.out$rfreq[, 1, , ] - opar.out$rfreq[sams, 1])^2, 3, FUN = mean)
            f1SMSE <- apply((bpar.out$rfreq[, 2, , ] - opar.out$rfreq[sams, 2])^2, 3, FUN = mean)
            f2SMSE <- apply((bpar.out$rfreq[, 3, , ] - opar.out$rfreq[sams, 3])^2, 3, FUN = mean)
            f3SMSE <- apply((bpar.out$rfreq[, 4, , ] - opar.out$rfreq[sams, 4])^2, 3, FUN = mean)
            f4SMSE <- apply((bpar.out$rfreq[, 5, , ] - opar.out$rfreq[sams, 5])^2, 3, FUN = mean)
            rfreqSMSE <- cbind(f0SMSE, f1SMSE, f2SMSE, f3SMSE, f4SMSE)

            d1SMSE <- apply((bpar.out$deciles[, 1, , ] - opar.out$deciles[sams, 1])^2, 3, FUN = mean)
            d2SMSE <- apply((bpar.out$deciles[, 2, , ] - opar.out$deciles[sams, 2])^2, 3, FUN = mean)
            d3SMSE <- apply((bpar.out$deciles[, 3, , ] - opar.out$deciles[sams, 3])^2, 3, FUN = mean)
            d4SMSE <- apply((bpar.out$deciles[, 4, , ] - opar.out$deciles[sams, 4])^2, 3, FUN = mean)
            d5SMSE <- apply((bpar.out$deciles[, 5, , ] - opar.out$deciles[sams, 5])^2, 3, FUN = mean)
            d6SMSE <- apply((bpar.out$deciles[, 6, , ] - opar.out$deciles[sams, 6])^2, 3, FUN = mean)
            d7SMSE <- apply((bpar.out$deciles[, 7, , ] - opar.out$deciles[sams, 7])^2, 3, FUN = mean)
            d8SMSE <- apply((bpar.out$deciles[, 8, , ] - opar.out$deciles[sams, 8])^2, 3, FUN = mean)
            d9SMSE <- apply((bpar.out$deciles[, 9, , ] - opar.out$deciles[sams, 9])^2, 3, FUN = mean)
            deciSMSE <- sqrt(cbind(d1SMSE, d2SMSE, d3SMSE, d4SMSE, d5SMSE, d6SMSE, d7SMSE, d8SMSE, d9SMSE))
      } else {
            meanBias <- mean(bpar.out$mean - opar.out$mean[sams])

            q1Bias <- mean(bpar.out$quartiles[, 1, , ] - opar.out$quartiles[sams, 1])
            q2Bias <- mean(bpar.out$quartiles[, 2, , ] - opar.out$quartiles[sams, 2])
            q3Bias <- mean(bpar.out$quartiles[, 3, , ] - opar.out$quartiles[sams, 3])
            ## quarBias<-apply(bpar.out$quartiles-opar.out$rquart[num.sam],c(2,4),FUN=mean)
            quarBias <- cbind(q1Bias, q2Bias, q3Bias)

            f0Bias <- mean(bpar.out$rfreq[, 1, , ] - opar.out$rfreq[sams, 1])
            f1Bias <- mean(bpar.out$rfreq[, 2, , ] - opar.out$rfreq[sams, 2])
            f2Bias <- mean(bpar.out$rfreq[, 3, , ] - opar.out$rfreq[sams, 3])
            f3Bias <- mean(bpar.out$rfreq[, 4, , ] - opar.out$rfreq[sams, 4])
            f4Bias <- mean(bpar.out$rfreq[, 5, , ] - opar.out$rfreq[sams, 5])
            rfreqBias <- cbind(f0Bias, f1Bias, f2Bias, f3Bias, f4Bias)

            d1Bias <- mean(bpar.out$deciles[, 1, , ] - opar.out$deciles[sams, 1])
            d2Bias <- mean(bpar.out$deciles[, 2, , ] - opar.out$deciles[sams, 2])
            d3Bias <- mean(bpar.out$deciles[, 3, , ] - opar.out$deciles[sams, 3])
            d4Bias <- mean(bpar.out$deciles[, 4, , ] - opar.out$deciles[sams, 4])
            d5Bias <- mean(bpar.out$deciles[, 5, , ] - opar.out$deciles[sams, 5])
            d6Bias <- mean(bpar.out$deciles[, 6, , ] - opar.out$deciles[sams, 6])
            d7Bias <- mean(bpar.out$deciles[, 7, , ] - opar.out$deciles[sams, 7])
            d8Bias <- mean(bpar.out$deciles[, 8, , ] - opar.out$deciles[sams, 8])
            d9Bias <- mean(bpar.out$deciles[, 9, , ] - opar.out$deciles[sams, 9])
            deciBias <- cbind(d1Bias, d2Bias, d3Bias, d4Bias, d5Bias, d6Bias, d7Bias, d8Bias, d9Bias)

            meanSMSE <- sqrt(mean((bpar.out$mean - opar.out$mean[sams])^2))

            q1SMSE <- mean((bpar.out$quartiles[, 1, , ] - opar.out$quartiles[sams, 1])^2)
            q2SMSE <- mean((bpar.out$quartiles[, 2, , ] - opar.out$quartiles[sams, 2])^2)
            q3SMSE <- mean((bpar.out$quartiles[, 3, , ] - opar.out$quartiles[sams, 3])^2)
            ## quarBias<-apply(bpar.out$quartiles-opar.out$rquart[num.sam],c(2,4),FUN=mean)
            quarSMSE <- sqrt(cbind(q1SMSE, q2SMSE, q3SMSE))

            f0SMSE <- mean((bpar.out$rfreq[, 1, , ] - opar.out$rfreq[sams, 1])^2)
            f1SMSE <- mean((bpar.out$rfreq[, 2, , ] - opar.out$rfreq[sams, 2])^2)
            f2SMSE <- mean((bpar.out$rfreq[, 3, , ] - opar.out$rfreq[sams, 3])^2)
            f3SMSE <- mean((bpar.out$rfreq[, 4, , ] - opar.out$rfreq[sams, 4])^2)
            f4SMSE <- mean((bpar.out$rfreq[, 5, , ] - opar.out$rfreq[sams, 5])^2)
            rfreqSMSE <- cbind(f0SMSE, f1SMSE, f2SMSE, f3SMSE, f4SMSE)

            d1SMSE <- mean((bpar.out$deciles[, 1, , ] - opar.out$deciles[sams, 1])^2)
            d2SMSE <- mean((bpar.out$deciles[, 2, , ] - opar.out$deciles[sams, 2])^2)
            d3SMSE <- mean((bpar.out$deciles[, 3, , ] - opar.out$deciles[sams, 3])^2)
            d4SMSE <- mean((bpar.out$deciles[, 4, , ] - opar.out$deciles[sams, 4])^2)
            d5SMSE <- mean((bpar.out$deciles[, 5, , ] - opar.out$deciles[sams, 5])^2)
            d6SMSE <- mean((bpar.out$deciles[, 6, , ] - opar.out$deciles[sams, 6])^2)
            d7SMSE <- mean((bpar.out$deciles[, 7, , ] - opar.out$deciles[sams, 7])^2)
            d8SMSE <- mean((bpar.out$deciles[, 8, , ] - opar.out$deciles[sams, 8])^2)
            d9SMSE <- mean((bpar.out$deciles[, 9, , ] - opar.out$deciles[sams, 9])^2)
            deciSMSE <- sqrt(cbind(d1SMSE, d2SMSE, d3SMSE, d4SMSE, d5SMSE, d6SMSE, d7SMSE, d8SMSE, d9SMSE))
      }
      # browser()
      BiasSMSE <- cbind(meanBias, quarBias, rfreqBias, deciBias, meanSMSE, quarSMSE, rfreqSMSE, deciSMSE)
      rnam <- colnames(BiasSMSE)
      BiasSMSE <- t(BiasSMSE)
      if (bpar.out$n.dist == 1) {
            colnames(BiasSMSE) <- "only seeds"
      } else {
            colnames(BiasSMSE) <- c("empd.w.p0s", "empd.nw.p0sEkb", "empd.nw.p0sEks")
      }
      BiasSMSE
}

# ---------------------------------------------------------------------------------------#

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
                        Obs.distrib <- Oempdegreedistrib(net, n.seed = i, n.neigh = j, num.sam = sam.size)
                        TMP <- Obs.distrib$seeds1
                  } else {
                        Obs.distrib <- Oempdegreedistrib(net, n.seed = i, n.neigh = j, num.sam = sam.size, seeds = TMP)
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


