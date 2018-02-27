#' Summary of the Bootstrap Degree Distribution
#'
#' This function provides summary statistics of a bootstrap degree
#' distribution.
#' @param boot_deg_out A list that is the output of \code{\link{boot_deg}}
#' @references
#' \insertRef{efron_bootstrap_1979}{snowboot}
#' @references
#' \insertRef{gel_bootstrap_2017}{snowboot}
#' @return A list consisting of:
#'    \item{mean}{An array of dimension
#'    \code{c(length(boot_deg_out$num_lsmi), boot_deg_out$boot_rep, 3)}.
#'    The last dimension, of 3, is for the three different methods of obtaining
#'    the empirical degree distribution from \code{boot_deg_out$empd}
#'    (see output empd from \code{\link{boot_deg}} for details).
#'    The \eqn{(i,j,k)}th element in the array is an estimate of mean degree for the
#'    \eqn{i}th LSMI sample, \eqn{j}th bootstrap replication, and \eqn{k}th empirical distribution
#'    from \code{boot_deg_out$empd}.}
#'    \item{quartiles}{An array of dimension
#'    \code{c(length(boot_deg_out$num_lsmi), 3, boot_deg_out$boot_rep, 3)}.
#'    The last dimension, of 3, is for the three different methods of estimation
#'    from \code{boot_deg_out$empd} (see output empd from \code{\link{boot_deg}}
#'    for details). The second dimension, of 3, corresponds to
#'    the quartiles (0.25, 0.5, 0.75). The \eqn{(i,j,k,l)}th element in the array is
#'    an estimate of \eqn{j}th quartile for the \eqn{i}th LSMI sample,
#'    \eqn{k}th bootstrap replication, and \eqn{l}th empirical distribution from
#'    \code{boot_deg_out$empd}.}
#'    \item{rfreq}{An array of dimension
#'    \code{c(length(boot_deg_out$num_lsmi), 5, boot_deg_out$boot_rep, 3)}.
#'    The last dimension, of 3, is for the three different methods of estimation
#'    from \code{boot_deg_out$empd}
#'    (see output empd from \code{\link{boot_deg}} for details.).
#'    The second dimension, of 5, corresponds to degree values: 0, 1, 2, 3, 4.
#'    The \eqn{(i,j,k,l)}th element in the array is the proportion of nodes
#'    with degree \eqn{j} in the \eqn{i}th LSMI sample, \eqn{k}th bootstrap replication,
#'    and \eqn{l}th empirical distribution from \code{boot_deg_out$empd}.}
#'    \item{deciles}{An array of dimension
#'          \code{c(length(boot_deg_out$num_lsmi), 9, boot_deg_out$boot_rep, 3)}.
#'          The last dimension, of 3, is for the three
#'          different methods of estimation from \code{boot_deg_out$empd}
#'          (see output empd from \code{\link{boot_deg}} for details.). The second
#'          dimension, of 9, corresponds to the deciles (0.1, 0.2, ..., 0.9).
#'          The \eqn{(i, j, k, l)}th element in the array is an estimate of \eqn{j}th
#'          decile for the \eqn{i}th LSMI sample, \eqn{k}th bootstrap replication,
#'          and \eqn{l}th empirical distribution from \code{boot_deg_out$empd}.}
#'    \item{num_lsmi}{Numeric indices corresponding to LSMI samples used for bootstrap.
#'          See value \code{num_lsmi} from \code{\link{boot_deg}}.}
#'    \item{seeds1}{A matrix of dimension \code{length(num_lsmi)} \eqn{\times} \code{n.seeds} with
#'          the numeric seed IDs. Each row corresponds to one LSMI. The rows are
#'          present in the same order as the IDs in \code{num_lsmi}.
#'          See value \code{seeds1} from \code{\link{boot_deg}}.}
#'    \item{nodes_of_lsmi}{A list of length \code{length(num_lsmi)} where each
#'          element is vector containing the numeric IDs of the nodes sampled
#'          using the respective LSMI. The elements are present in the same
#'          order as the IDs in \code{num_lsmi}.
#'          Note: \code{nodes_of_lsmi} is unreported when n.neigh equals zero.
#'          See value \code{nodes_of_lsmi} from \code{\link{boot_deg}}.}
#' @export
#' @examples
#' net <- artificial_networks[[1]]
#' sam.out <- empd_deg_lsmi(net = net, n.seeds = 40, n.neigh = 1, num_lsmi = 1)
#' boot_deg_out <- boot_deg(sam.out = sam.out, boot_rep = 50)
#' a <- boot_parameters_est(boot_deg_out)

boot_parameters_est <- function(boot_deg_out) {
  # boot_deg_out is output of boot_deg
  tn.sam <- length(boot_deg_out$num_lsmi)
  if (boot_deg_out$n.neigh == 0) {
    n.dist <- 1  #n.dist is the number of different emp distr.
  } else {
    n.dist <- 3
  }

  mean <- array(NA, c(tn.sam, boot_deg_out$boot_rep, n.dist))
  quartiles <- array(NA, c(tn.sam, 3, boot_deg_out$boot_rep, n.dist))
  rfreq <- array(NA, c(tn.sam, 5, boot_deg_out$boot_rep, n.dist))  #freq of 0,1,2,3,4
  deciles <- array(NA, c(tn.sam, 9, boot_deg_out$boot_rep, n.dist))  #10,20,30,40,50,60,70,80,90

  for (m in 1:tn.sam) {
    w <- 1
    in.while <- TRUE
            while (in.while) {
                  if (boot_deg_out$n.neigh == 0) {
                        empd <- boot_deg_out$empd[[m]]$empd.seeds
                        in.while <- FALSE
                  } else {
                        if (w == 1) {
                          empd <- boot_deg_out$empd[[m]]$empd.w.p0s
                          in.while <- FALSE
                        } else if (w == 2) {
                          empd <- boot_deg_out$empd[[m]]$empd.nw.p0sEkb
                        } else if (w == 3) {
                          empd <- boot_deg_out$empd[[m]]$empd.nw.p0sEks
                          in.while <- FALSE
                        }
                  }
              vals <- boot_deg_out$values[[m]]  #as.numeric(colnames(empd))  ###checar

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
      list(mean = mean, quartiles = quartiles, rfreq = rfreq, deciles = deciles,
           n.dist = n.dist, num_lsmi = boot_deg_out$num_lsmi,
           seeds1 = boot_deg_out$seeds1, nodes_of_lsmi = boot_deg_out$nodes_of_lsmi)
}  #function
