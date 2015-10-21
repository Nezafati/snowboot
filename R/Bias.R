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
