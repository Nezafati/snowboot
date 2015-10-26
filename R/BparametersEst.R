#' Summary of a Network Degree Sequence
#'
#' This function provides summary statistics of a network degree
#' distribution.
#' @param outBootdeg a list that is the output of \code{\link{bootdeg}}
#' @return a list consisting of:
#'    \item{mean}{an array of dimension \code{length(outBootdeg$num.sam) x
#'          outBootdeg$n.boot x 3}. The last dimension, of 3, is for the three
#'          different methods of obtaining the empirical degree distribution
#'          from \code{outBootdeg$empd}
#'          (see output empd from \code{\link{bootdeg}} for details.).
#'          Each element in the array is an estimation of the mean for the
#'          respective LSMI, bootstrap replication, and empirical distribution
#'          from \code{outBootdeg$empd}.}
#'    \item{quartiles}{an array of dimension \code{length(outBootdeg$num.sam) x 3 x
#'          outBootdeg$n.boot x 3}. The last dimension, of 3, is for the three
#'          different methods of estimation from \code{outBootdeg$empd}
#'          (see output empd from \code{\link{bootdeg}} for details.). The second
#'          dimension, of 3, corresponds to the quartiles (.25, .5, .75).
#'          Each element in the array is a quartile for the
#'          respective LSMI, quartile (.25, .5, .75), for a bootstrap replication,
#'          and empirical distribution from \code{outBootdeg$empd}.}
#'    \item{rfreq}{an array of dimension \code{length(outBootdeg$num.sam) x 5 x
#'          outBootdeg$n.boot x 3}. The last dimension, of 3, is for the three
#'          different methods of estimation from \code{outBootdeg$empd}
#'          (see output empd from \code{\link{bootdeg}} for details.). The second
#'          dimension, of 5, corresponds to degree values: 0, 1, 2, 3, 4.
#'          Each element in the array is the proportion of nodes from the
#'          respective LSMI, with degree value (0-4), for a bootstrap replication,
#'          and empirical distribution from \code{outBootdeg$empd}.}
#'    \item{deciles}{an array of dimension \code{length(outBootdeg$num.sam) x 9 x
#'          outBootdeg$n.boot x 3}. The last dimension, of 3, is for the three
#'          different methods of estimation from \code{outBootdeg$empd}
#'          (see output empd from \code{\link{bootdeg}} for details.). The second
#'          dimension, of 9, corresponds to the deciles (.1, .2, ... , .9).
#'          Each element in the array is a decile for the
#'          respective LSMI, decile (.1, .2, ... , .9), for a bootstrap replication,
#'          and empirical distribution from \code{outBootdeg$empd}.}
#'    \item{num.sam}{same as \code{outBootdeg$num.sam}}
#' @export

BparametersEst <- function(outBootdeg) {
      # outBootdeg is output of bootdeg
      tn.sam <- length(outBootdeg$num.sam)
      if (outBootdeg$n.neigh == 0) {
            n.dist <- 1  #n.dist is the number of different emp distr.
      } else {
            n.dist <- 3
      }

      mean <- array(NA, c(tn.sam, outBootdeg$n.boot, n.dist))
      quartiles <- array(NA, c(tn.sam, 3, outBootdeg$n.boot, n.dist))
      rfreq <- array(NA, c(tn.sam, 5, outBootdeg$n.boot, n.dist))  #freq of 0,1,2,3,4
      deciles <- array(NA, c(tn.sam, 9, outBootdeg$n.boot, n.dist))  #10,20,30,40,50,60,70,80,90

      for (m in 1:tn.sam) {
            w <- 1
            in.while <- TRUE
            while (in.while) {
                  if (outBootdeg$n.neigh == 0) {
                        empd <- outBootdeg$empd[[m]]$empd.seed
                        in.while <- FALSE
                  } else {
                        if (w == 1) {
                              empd <- outBootdeg$empd[[m]]$empd.w.p0s
                        } else if (w == 2) {
                              empd <- outBootdeg$empd[[m]]$empd.nw.p0sEkb
                        } else if (w == 3) {
                              empd <- outBootdeg$empd[[m]]$empd.nw.p0sEks
                              in.while <- FALSE
                        }
                  }
                  vals <- outBootdeg$values[[m]]  #as.numeric(colnames(empd))  ###checar

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
      list(mean = mean, quartiles = quartiles, rfreq = rfreq, deciles = deciles, n.dist = n.dist, num.sam = outBootdeg$num.sam)
}  #function
