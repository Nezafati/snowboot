BMean.intervals <- function(Opar, Bpar, n.dist) {

      # Opar is output of OparametersEst Opar$mean Bpar is output of BparametersEst Bpar$mean n.dist is the number of empirical
      # distributions considered browser()
      if (n.dist > 1) {
            mean.int <- apply(Bpar, c(1, 3), FUN = quantile, prob = c(0.025, 0.975))  #debe ser una matriz 2x tn.sam x n.dist
            mean.int1 <- cbind(2 * Opar - mean.int[2, , ], 2 * Opar - mean.int[1, , ])[, c(1, 4, 2, 5, 3, 6)]  #firts to columns dist1, then dist2 and finally dist3
            mean.int2 <- cbind(Opar + apply(Bpar, c(1, 3), FUN = mean) - mean.int[2, , ], Opar + apply(Bpar, c(1, 3), FUN = mean) -
                                     mean.int[1, , ])[, c(1, 4, 2, 5, 3, 6)]
            mean.int3 <- cbind(Opar^2/mean.int[2, , ], 2 * Opar^2/mean.int[1, , ])[, c(1, 4, 2, 5, 3, 6)]
            mean.int <- cbind(t(mean.int[, , 1]), t(mean.int[, , 2]), t(mean.int[, , 3]))
      } else if (n.dist == 1)
      {
            mean.int <- t(apply(Bpar, c(1), FUN = quantile, prob = c(0.025, 0.975)))  #debe ser una matriz tn.sam x2
            mean.int1 <- cbind(2 * Opar - mean.int[, 2], 2 * Opar - mean.int[, 1])
            mean.int2 <- cbind(Opar + apply(Bpar, c(1), FUN = mean) - mean.int[, 2], Opar + apply(Bpar, c(1), FUN = mean) -
                                     mean.int[, 1])
            mean.int3 <- cbind(Opar^2/mean.int[, 2], 2 * Opar^2/mean.int[, 1])
      }  #if (n.dist>1)
      list(mean.int = mean.int, mean.int1 = mean.int1, mean.int2 = mean.int2, mean.int3 = mean.int3)
}
