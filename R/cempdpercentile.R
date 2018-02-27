cempdpercentile <- function(cempd, perc, vals) {
      # cempd is a matrix perc is a scalar
  res <- vals[apply(X = cempd, 1, FUN = min.greater, x = perc, ge = TRUE)]
  res
}

