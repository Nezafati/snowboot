real.parameters <- function(net) {
      # this function obtains the real parameters in a network
      realdd <- net$degree - net$degree.left
      rmean <- mean(realdd)
      rquart <- quantile(realdd, prob = c(0.25, 0.5, 0.7))
      rfreq <- c(sum(realdd == 0), sum(realdd == 1), sum(realdd == 2), sum(realdd == 3), sum(realdd == 4))/length(net$degree)
      rdeci <- quantile(realdd, prob = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
      list(realdd = realdd, rmean = rmean, rquart = rquart, rfreq = rfreq, rdeci = rdeci)
}
