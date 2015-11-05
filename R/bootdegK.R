bootdegK <- function(sam.out, num.sam, n.boot) {
      #see bootdeg for description.
      n.seeds <- sam.out$n.seeds
      n.neigh <- sam.out$n.neigh
      if (length(num.sam) == 1)
            num.sam <- 1:num.sam

      empd <- as.list(rep(NA, length(num.sam)))
      i <- 1
      for (m in num.sam) {
            # Boostrap samples of seeds, nonseeds-noWeighted and nonseeds-Weighted:
            val.seeds <- sam.out$val.seeds[[m]]
            val.nonseeds <- sam.out$val.nonseeds[[m]]
            freq.deg.seeds <- sam.out$samples[[m]]$freq.deg.seeds
            freq.deg.nonseeds <- sam.out$samples[[m]]$freq.deg.nonseeds

            bsam.seeds <- myBsample(val.seeds, n.seeds, n.boot, prob = freq.deg.seeds)  #matrix n.boot x n.seeds
            bsam.nonseeds.nw <- myBsample(val.nonseeds, sum(freq.deg.nonseeds), n.boot, prob = freq.deg.nonseeds)  #matrix n.boot x sum(freq.deg.nonseeds)
            bsam.nonseeds.w <- myBsample(val.nonseeds, sum(freq.deg.nonseeds), n.boot, prob = freq.deg.nonseeds/val.nonseeds)  #matrix

            p0.B <- rep(0, n.boot)
            if (any(val.seeds == 0)) {
                  # if any seeds has degree zero
                  p0.B <- rowSums(bsam.seeds == 0)/n.seeds
                  #^the estimation from the bootstrap samples
            }

            values <- sam.out$values[[m]]
            # ^all the possible degree values to resample

            ###################### Frequency ##### (Not the relative frequency)
            Fseeds <- t(apply(bsam.seeds, 1, table.row, vect = values))
            if (is.null(bsam.nonseeds.nw)) {
                  # browser()
                  Fnonseeds.nw <- 0
            } else {
                  Fnonseeds.nw <- t(apply(as.matrix(bsam.nonseeds.nw), 1, table.row, vect = values))
            }
            if (is.null(bsam.nonseeds.w)) {
                  Fnonseeds.w <- 0
            } else {
                  Fnonseeds.w <- t(apply(as.matrix(bsam.nonseeds.w), 1, table.row, vect = values))
            }

            ############### combining information from seeds and nonseeds ######

            # mean degree computed from the original sampled seeds:
            ekseed <- sam.out$ekseed[[m]]

            colzero <- NULL
            if (any(values == 0)) {
                  colzero <- which(values == 0)
                  vals <- values[-colzero]
                  f.seeds <- Fseeds[, -colzero]
                  f.nonseeds.nw <- Fnonseeds.nw[, -colzero]
                  f.nonseeds.w <- Fnonseeds.w[, -colzero]
            } else {
                  vals <- values
                  f.seeds <- Fseeds
                  f.nonseeds.nw <- Fnonseeds.nw
                  f.nonseeds.w <- Fnonseeds.w
            }
            empd.w.p0s <- (f.seeds + f.nonseeds.w * (1 - p0.B))/(n.seeds + sum(freq.deg.nonseeds))
            empd.nw.p0sEkb <- (f.seeds + t(t(f.nonseeds.nw)/vals) * (1 - p0.B) * apply(bsam.seeds, 1, FUN = mean))/(n.seeds + rowSums(t(t(f.nonseeds.nw)/vals)) *
                                                                                                                       apply(bsam.seeds, 1, FUN = mean))
            empd.nw.p0sEks <- (f.seeds + ekseed * t(t(f.nonseeds.nw)/vals) * (1 - p0.B))/(n.seeds + ekseed * rowSums(t(t(f.nonseeds.nw)/vals)))

            if (any(values == 0)) {
                  empd.w.p0s <- cbind(`0` = p0.B, empd.w.p0s)
                  empd.nw.p0sEkb <- cbind(`0` = p0.B, empd.nw.p0sEkb)
                  empd.nw.p0sEks <- cbind(`0` = p0.B, empd.nw.p0sEks)
            }
            empd[[i]] <- list(empd.w.p0s = empd.w.p0s, empd.nw.p0sEkb = empd.nw.p0sEkb, empd.nw.p0sEks = empd.nw.p0sEks)
            i <- i + 1
      }
      list(values = sam.out$values, empd = empd, num.sam = num.sam,
           n.boot = n.boot, n.neigh = n.neigh, seeds1 = sam.out$seeds1[num.sam,],
           nodes_of_LSMI = sam.out$nodes_of_LSMI[num.sam])
}
