#' 200 Simulated Networks of order 2000 with Polylogarithmic (0.1, 2)
#' Degree Distributions
#'
#' A list called "networks" containing 200 network objects of order 2000. These
#' networks were simulated using the polylogarithmic (aka Gutenberg--Richter law)
#' degree distribution (Newman et al., 2001; Newman, 2002) with parameters
#' \eqn{\delta = 0.1} and \eqn{\lambda = 2} as see in the following equations:
#' \deqn{f(k) = k^{-{\delta}}e^{-{k/{\lambda}}}/Li_{\delta}(e^{-{1/\lambda}})}{f(k)=k^-\delta exp(-k/\lambda )/Li[\delta](exp(-1/\lambda))}
#' \deqn{Li=\sum_{j=1}^{\infty} z^{-j}/{j^{\delta}},}{Li=\sum_{j=1}^{\infty} z^{-j}/{j^{\delta}},}
#' where \eqn{\lambda > 0}. Please see refence below for details (Thompson  et al, 2016).
#' @references Gel, Y. R., Lyubchich, V., & Ramirez Ramirez, L. L. (2017).
#' Bootstrap quantification of estimation uncertainties in
#' network degree distributions. Scientific Reports, 7, 5807.
#' \url{http://doi.org/10.1038/s41598-017-05885-x}
#' @references Newman, M. E., Strogatz, S. H., & Watts, D. J. (2001).
#' Random graphs with arbitrary degree distributions and their applications.
#' Physical review E, 64(2), 026118.
#' @format a list containing 200 network objects. Each network object is a list
#' with three elements:
#' \describe{
#'    \item{edges}{edgelist of the network. A two column
#'      \code{matrix} where each row is an edge.}
#'    \item{degree}{degree sequence of the network, which is
#'      an \code{integer} vector of length $n$.}
#'    \item{n}{network order. The order is 2000.}
#' }



"artificial_networks"
