Bintervalsi <- function(Opari, Bpari, n.dist, num.par) {
      # Opari is the parameter(s) of interesout from OparametersEst (either: mean,quartiles,rfreq,deciles) Bpari is the
      # parameter(s) of interes from BparametersEst (either: mean,quartiles,rfreq,deciles) n.dist is the number of empirical
      # distributions considered

      reorder <- as.vector(matrix(rep(1:(num.par), each = 2), 2, num.par) + c(0, num.par))
      # browser()
      if (n.dist > 1) {
            par.int <- apply(Bpari, c(1, 2, 4), FUN = quantile, prob = c(0.025, 0.975))  #debe ser una matriz: 2 x num.sam x num param x n.dist

            par.int1.a <- apply(-par.int[2, , , ], 3, FUN = "+", 2 * Opari)  #matrix (num.sam)(num.par) x ndist #the first col is the vect of a matrix num.sam x num.par con dist1
            par.int1.b <- apply(-par.int[1, , , ], 3, FUN = "+", 2 * Opari)  # '
            par.int1 <- cbind(par.int1.a, par.int1.b)[, reorder]  #the first two columns are the lower and upper limit according to n.dist1...

            ##### test for the quartiles plot(1:60,rep(5,60),t='n',ylim=c(-3,6)) segments(1:60,par.int[,1],1:60,par.int[,2])

            par.int2.a <- apply(apply(Bpari, c(1, 2, 4), FUN = mean) - par.int[2, , , ], 3, FUN = "+", Opari)
            par.int2.b <- apply(apply(Bpari, c(1, 2, 4), FUN = mean) - par.int[1, , , ], 3, FUN = "+", Opari)
            par.int2 <- cbind(par.int2.a, par.int2.b)[, reorder]

            par.int3.a <- apply(1/par.int[2, , , ], 3, FUN = "*", Opari)
            par.int3.b <- apply(1/par.int[1, , , ], 3, FUN = "*", Opari)
            par.int3 <- cbind(par.int3.a, par.int3.b)[, reorder]

            ##### segments(1:60,par.int3[,1],1:60,par.int3[,2])
      } else if (n.dist == 1)
      {
            par.int <- apply(Bpari, c(1, 2), FUN = quantile, prob = c(0.025, 0.975))  #debe ser una matriz 2 x num.sam x num.par
            par.int1 <- cbind(2 * Opari - par.int[2, , ], 2 * Opari - par.int[1, , ])[, reorder]  #debe ser una matriz num.sam x 2(num.par)
            par.int2 <- cbind(Opari + apply(Bpari, c(1, 2), FUN = mean) - par.int[2, , ], Opari + apply(Bpari, c(1, 2), FUN = mean) -
                                    par.int[1, , ])[, reorder]
            par.int3 <- cbind(Opari^2/par.int[2, , ], 2 * Opari^2/par.int[1, , ])[, reorder]
            mpar.int <- NULL
            for (w in 1:num.par) mpar.int <- cbind(mpar.int, t(par.int[, , w]))
      }  #if (n.dist>1)
      list(par.int = mpar.int, par.int1 = par.int1, part.int2 = par.int2, par.int3 = par.int3, num.par = num.par)
}
