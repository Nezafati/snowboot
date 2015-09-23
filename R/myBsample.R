myBsample <- function(val, length, n.boot, prob = NULL) {
      # val is a vector to resample the sample size is length*n.boot
      if (length(val) == 0) {
            res <- NULL
      } else if (length(val) == 1) {
            # problems with weighted sample if only one term
            res <- matrix(rep(val, length * n.boot), n.boot, length, byrow = TRUE)
      } else if (length(val) > 1) {
            res <- matrix(sample(val, length * n.boot, replace = TRUE, prob = prob), n.boot, length, byrow = TRUE)
      }
      res
}
