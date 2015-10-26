fallInsMean <- function(ints, rpar, n.dist) {
      # ints in output of Bintervals (mean.int, mean.int1,mean.int2,mean.int3) rpar is output of summary.net rpar$rmean
      # n.dist is the number of empirical distributions browser()
      if (n.dist > 1) {
            fi.mean <- fi1.mean <- fi2.mean <- fi3.mean <- rep(0, 3)
            col <- c(1, 3, 5)
            for (k in 1:3) {
                  fi.mean[k] <- sum(ints$mean.int[, col[k]] <= rpar & rpar <= ints$mean.int[, col[k] + 1])/dim(ints$mean.int)[1]
                  fi1.mean[k] <- sum(ints$mean.int1[, col[k]] <= rpar & rpar <= ints$mean.int1[, col[k] + 1])/dim(ints$mean.int)[1]
                  fi2.mean[k] <- sum(ints$mean.int2[, col[k]] <= rpar & rpar <= ints$mean.int2[, col[k] + 1])/dim(ints$mean.int)[1]
                  fi3.mean[k] <- sum(ints$mean.int3[, col[k]] <= rpar & rpar <= ints$mean.int3[, col[k] + 1])/dim(ints$mean.int)[1]
            }
      } else {
            fi.mean <- sum(ints$mean.int[, 1] <= rpar & rpar <= ints$mean.int[, 2])/dim(ints$mean.int)[1]
            fi1.mean <- sum(ints$mean.int1[, 1] <= rpar & rpar <= ints$mean.int1[, 2])/dim(ints$mean.int)[1]
            fi2.mean <- sum(ints$mean.int2[, 1] <= rpar & rpar <= ints$mean.int2[, 2])/dim(ints$mean.int)[1]
            fi3.mean <- sum(ints$mean.int3[, 1] <= rpar & rpar <= ints$mean.int3[, 2])/dim(ints$mean.int)[1]
      }
      fi <- rbind(fi.mean = fi.mean, fi1.mean = fi1.mean, fi2.mean = fi2.mean, fi3.mean = fi3.mean)
      fi
}
