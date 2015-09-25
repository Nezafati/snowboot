BempdegreedistribK <- function(sam.out, num.sam, n.boot, idname = "Temp") {
      # This function obtains the bootstrap samples for each sample from a network sam.out is the output of Oempdegreedistrib
      # num.sam is the number of different samples taken from the same network. Scalar o vector. n.boot is the number of
      # bootstrap samples taken from each sample
      n.seed <- sam.out$n.seed
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

            bsam.seed <- myBsample(val.seed, n.seed, n.boot, prob = freq.deg.seed)  #matrix n.boot x n.seed
            bsam.nonseed.nw <- myBsample(val.nonseed, sum(freq.deg.nonseed), n.boot, prob = freq.deg.nonseed)  #matrix n.boot x sum(freq.deg.nonseed)
            bsam.nonseed.w <- myBsample(val.nonseed, sum(freq.deg.nonseed), n.boot, prob = freq.deg.nonseed/val.nonseed)  #matrix

            p0.B <- rep(0, n.boot)
            if (any(val.seed == 0)) {
                  # if any seed has degree zero
                  p0.B <- rowSums(bsam.seed == 0)/n.seed  #the estimation from the bootstrap samples
            }
            p0.real <- sam.out$p0.real
            p0.seed <- sam.out$p0.seed[[m]]

            values <- sam.out$values[[m]]  #all the possible degree values to resample
            # if (all(values==0)) {values<-0.1}

            ###################### Frequency ##### (Not the relative frequency)
            Fseed <- t(apply(bsam.seed, 1, table.row, vect = values))  #frequency (sorted according to values)
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

            # Fnonseed.nw<-t(apply(as.matrix(bsam.nonseed.nw),1,table.row,vect=values)) # '
            # Fnonseed.w<-t(apply(as.matrix(bsam.nonseed.w),1,table.row,vect=values)) # '

            #################################################################### combining information from seeds and nonseeds ######

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
            # browser() WB # seeds and weighted nonseeds ###### consider the p0 fixed from the seed information:
            # empd.w.p0s<-(f.seed+(1-p0.seed)*f.nonseed.w)/(n.seed+sum(freq.deg.nonseed))
            empd.w.p0s <- (f.seed + f.nonseed.w * (1 - p0.B))/(n.seed + sum(freq.deg.nonseed))
            #### consider the p0 fixed from the real information:
            #### empd.w.p0r<-(f.seed+(1-p0.real)*f.nonseed.w)/(n.seed+sum(freq.deg.nonseed)) NWB # seeds and non weighted nonseeds ####
            #### p0 estimated from orginal sampled seeds# E(K) estimated from bootstrap samples from the seeds
            #### empd.nw.p0sEkb<-(f.seed+(1-p0.seed)*apply(bsam.seed,1,FUN=mean)*t(t(f.nonseed.nw)/vals))/(n.seed+
            #### apply(bsam.seed,1,FUN=mean)*rowSums(t(t(f.nonseed.nw)/vals)))
            empd.nw.p0sEkb <- (f.seed + t(t(f.nonseed.nw)/vals) * (1 - p0.B) * apply(bsam.seed, 1, FUN = mean))/(n.seed + rowSums(t(t(f.nonseed.nw)/vals)) *
                                                                                                                       apply(bsam.seed, 1, FUN = mean))
            # E(K) estimated from the original seeds sample
            # empd.nw.p0sEks<-(f.seed+(1-p0.seed)*ekseed*t(t(f.nonseed.nw)/vals))/(n.seed+ ekseed*rowSums(t(t(f.nonseed.nw)/vals)))
            empd.nw.p0sEks <- (f.seed + ekseed * t(t(f.nonseed.nw)/vals) * (1 - p0.B))/(n.seed + ekseed * rowSums(t(t(f.nonseed.nw)/vals)))
            ######################################### p0 taken as known # E(K) estimated from bootstrap samples from the seeds
            ######################################### empd.nw.p0rEkb<-(f.seed+(1-p0.real)*apply(bsam.seed,1,FUN=mean)*t(t(f.nonseed.nw)/vals))/
            ######################################### (n.seed+apply(bsam.seed,1,FUN=mean)*rowSums(t(t(f.nonseed.nw)/vals))) E(K) estimated from the original seeds sample
            ######################################### empd.nw.p0rEks<-(f.seed+(1-p0.real)*ekseed*t(t(f.nonseed.nw)/vals))/(n.seed+ ekseed*rowSums(t(t(f.nonseed.nw)/vals)))

            if (any(values == 0)) {
                  empd.w.p0s <- cbind(`0` = p0.B, empd.w.p0s)  #1
                  # empd.w.p0r<-cbind('0'=p0.real,empd.w.p0r) #2
                  empd.nw.p0sEkb <- cbind(`0` = p0.B, empd.nw.p0sEkb)  #3
                  empd.nw.p0sEks <- cbind(`0` = p0.B, empd.nw.p0sEks)  #4
                  # empd.nw.p0rEkb<-cbind('0'=p0.real,empd.nw.p0rEkb) #6 empd.nw.p0rEks<-cbind('0'=p0.real,empd.nw.p0rEks) #7
            }
            empd[[m]] <- list(empd.w.p0s = empd.w.p0s, empd.nw.p0sEkb = empd.nw.p0sEkb, empd.nw.p0sEks = empd.nw.p0sEks)
      }  # for(m in num.sam)
      list(idname = idname, values = sam.out$values, empd = empd, num.sam = num.sam, n.boot = n.boot, n.neigh = n.neigh)
}
