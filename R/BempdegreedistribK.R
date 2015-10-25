BempdegreedistribK <- function(sam.out, num.sam, n.boot) {
      #see Bempdegreedistrib for description.
      n.seeds <- sam.out$n.seeds
      n.neigh <- sam.out$n.neigh
      if (length(num.sam) == 1)
            num.sam <- 1:num.sam

      empd <- as.list(rep(NA, length(num.sam)))
      i <- 1
      for (m in num.sam) {
            # if(i%%100==1)#cat('Processing bootstrap samples of sample=',i,'\n') #print every 100
            i <- i + 1
            # Boostrap samples of seed, nonseeds-noWeighted and nonseeds-Weighted:
            val.seed <- sam.out$val.seed[[m]]
            val.nonseed <- sam.out$val.nonseed[[m]]
            freq.deg.seed <- sam.out$samples[[m]]$freq.deg.seed
            freq.deg.nonseed <- sam.out$samples[[m]]$freq.deg.nonseed

            bsam.seed <- myBsample(val.seed, n.seeds, n.boot, prob = freq.deg.seed)  #matrix n.boot x n.seeds
            bsam.nonseed.nw <- myBsample(val.nonseed, sum(freq.deg.nonseed), n.boot, prob = freq.deg.nonseed)  #matrix n.boot x sum(freq.deg.nonseed)
            bsam.nonseed.w <- myBsample(val.nonseed, sum(freq.deg.nonseed), n.boot, prob = freq.deg.nonseed/val.nonseed)  #matrix

            p0.B <- rep(0, n.boot)
            if (any(val.seed == 0)) {
                  # if any seed has degree zero
                  p0.B <- rowSums(bsam.seed == 0)/n.seeds
                  #^the estimation from the bootstrap samples
            }

            values <- sam.out$values[[m]]
            # ^all the possible degree values to resample


            ###################### Frequency ##### (Not the relative frequency)
            Fseed <- t(apply(bsam.seed, 1, table.row, vect = values))
            # browser()
            if (is.null(bsam.nonseed.nw)) {
                  # browser()
                  Fnonseed.nw <- 0
            } else {
                  Fnonseed.nw <- t(apply(as.matrix(bsam.nonseed.nw), 1, table.row, vect = values))
            }
            if (is.null(bsam.nonseed.w)) {
                  Fnonseed.n <- 0
            } else {
                  Fnonseed.w <- t(apply(as.matrix(bsam.nonseed.w), 1, table.row, vect = values))
            }

            ############### combining information from seeds and nonseeds ######

            # mean degree computed from the original sampled seeds:
            ekseed <- sam.out$ekseed[[m]]

            colzero <- NULL
            if (any(values == 0)) {
                  colzero <- which(values == 0)
                  vals <- values[-colzero]
                  f.seed <- Fseed[, -colzero]
                  f.nonseed.nw <- Fnonseed.nw[, -colzero]
                  f.nonseed.w <- Fnonseed.w[, -colzero]
            } else {
                  vals <- values
                  f.seed <- Fseed
                  f.nonseed.nw <- Fnonseed.nw
                  f.nonseed.w <- Fnonseed.w
            }
            empd.w.p0s <- (f.seed + f.nonseed.w * (1 - p0.B))/(n.seeds + sum(freq.deg.nonseed))
            empd.nw.p0sEkb <- (f.seed + t(t(f.nonseed.nw)/vals) * (1 - p0.B) * apply(bsam.seed, 1, FUN = mean))/(n.seeds + rowSums(t(t(f.nonseed.nw)/vals)) *
                                                                                                                       apply(bsam.seed, 1, FUN = mean))
            empd.nw.p0sEks <- (f.seed + ekseed * t(t(f.nonseed.nw)/vals) * (1 - p0.B))/(n.seeds + ekseed * rowSums(t(t(f.nonseed.nw)/vals)))

            if (any(values == 0)) {
                  empd.w.p0s <- cbind(`0` = p0.B, empd.w.p0s)
                  empd.nw.p0sEkb <- cbind(`0` = p0.B, empd.nw.p0sEkb)
                  empd.nw.p0sEks <- cbind(`0` = p0.B, empd.nw.p0sEks)
            }
            empd[[m]] <- list(empd.w.p0s = empd.w.p0s, empd.nw.p0sEkb = empd.nw.p0sEkb, empd.nw.p0sEks = empd.nw.p0sEks)
      }
      list(values = sam.out$values, empd = empd, num.sam = num.sam,
           n.boot = n.boot, n.neigh = n.neigh, seeds1 = sam.out$seeds1[num.sam,],
           nodes_of_LSMI = sam.out$nodes_of_LSMI[num.sam])
}
