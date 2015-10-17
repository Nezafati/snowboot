Bempdegreedistrib0 <- function(sam.out, num.sam, n.boot) {
      # This function obtains the bootstrap samples for each sample from a network sam.out is the output of Oempdegreedistrib
      # num.sam is the number of different samples taken from the same network (Scalar or vector) n.boot is the number of
      # bootstrap samples taken from each sample
      n.seeds <- sam.out$n.seeds
      n.neigh <- sam.out$n.neigh
      if (length(num.sam) == 1)
            num.sam <- 1:num.sam
      empd <- as.list(rep(NA, length(num.sam)))
      i <- 1
      for (m in num.sam) {
            # if(i%%100==1)#cat('Processing bootstrap samples of sample=',i,'\n') #print every 100
            i <- i + 1
            val.seed <- sam.out$val.seed[[m]]
            freq.deg.seed <- sam.out$samples[[m]]$freq.deg.seed
            bsam.seed <- myBsample(val.seed, n.seeds, n.boot, prob = freq.deg.seed)
            values <- sam.out$values[[m]]  #all the possible degree values toresample
            #### Frequency ##### (Not the relative frequency)
            Fseed <- t(apply(bsam.seed, 1, table.row, vect = values))  #freq (sorted according to values)
            # browser()
            empd.seed <- Fseed/n.seeds
            empd[[m]] <- list(empd.seed = empd.seed)
      }  # for(m in num.sam)
      list(values = sam.out$values, empd = empd, num.sam = num.sam, n.boot = n.boot, n.neigh = n.neigh)
}
