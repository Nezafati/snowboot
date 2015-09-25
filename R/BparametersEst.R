BparametersEst <- function(outBempd) {
      # outBempd is output of Bempdegreedistrib
      tn.sam <- length(outBempd$num.sam)
      if (outBempd$n.neigh == 0) {
            n.dist <- 1  #n.dist is the number of different emp distr.
      } else {
            n.dist <- 3
      }

      mean <- array(NA, c(tn.sam, outBempd$n.boot, n.dist))
      quartiles <- array(NA, c(tn.sam, 3, outBempd$n.boot, n.dist))
      rfreq <- array(NA, c(tn.sam, 5, outBempd$n.boot, n.dist))  #freq of 0,1,2,3,4
      deciles <- array(NA, c(tn.sam, 9, outBempd$n.boot, n.dist))  #10,20,30,40,50,60,70,80,90

      for (m in 1:tn.sam) {
            w <- 1
            in.while <- TRUE
            while (in.while) {
                  if (outBempd$n.neigh == 0) {
                        empd <- outBempd$empd[[m]]$empd.seed
                        in.while <- FALSE
                  } else {
                        if (w == 1) {
                              empd <- outBempd$empd[[m]]$empd.w.p0s
                        } else if (w == 2) {
                              empd <- outBempd$empd[[m]]$empd.nw.p0sEkb
                        } else if (w == 3) {
                              empd <- outBempd$empd[[m]]$empd.nw.p0sEks
                              in.while <- FALSE
                        }
                  }
                  vals <- outBempd$values[[m]]  #as.numeric(colnames(empd))  ###checar

                  if (dim(empd)[2] != length(as.vector(vals))) {
                        # browser()
                        mean[m, , w] <- t(empd) %*% as.vector(vals)
                  } else mean[m, , w] <- empd %*% as.vector(vals)

                  cempd <- t(apply(empd, 1, FUN = cumsum))
                  quartiles[m, , , w] <- t(sapply(X = c(0.25, 0.5, 0.75), FUN = cempdpercentile, cempd = cempd, vals = vals))
                  rfreq[m, , , w] <- t(sapply(X = 0:4, FUN = distribvalsmat, empd = empd, vals = vals))
                  deciles[m, , , w] <- t(sapply(X = seq(0.1, 0.9, by = 0.1), FUN = cempdpercentile, cempd = cempd, vals = vals))
                  w <- w + 1
            }  #while(w<=3)
      }  #for (m in 1:...)
      # browser()
      list(mean = mean, quartiles = quartiles, rfreq = rfreq, deciles = deciles, n.dist = n.dist, num.sam = outBempd$num.sam)
}  #function
