myBsample <- function(val, original_sample_size, boot_rep, prob = NULL) {
      # val is a vector to resample the sample size is original_sample_size*boot_rep
      if (length(val) == 0) {
            res <- NULL
      } else if (length(val) == 1) {
            # problems with weighted sample if only one term
            res <- matrix(rep(val, original_sample_size * boot_rep), boot_rep,
                          original_sample_size, byrow = TRUE)
      } else if (length(val) > 1) {
        res <- matrix(sample(val, original_sample_size * boot_rep,
                             replace = TRUE, prob = prob), boot_rep,
                      original_sample_size, byrow = TRUE)
      }
  res
}
