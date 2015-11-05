#' 200 Simulated Networks from Polylog Degree Distributions. (reference)
#'
#' A list called "networks" containing 200 network objects of order 2000. These
#' networks were simulated using the polylog degree distribution with paramaters
#' "0.1" and "2". add reference and formula
#'
#' @format a list containg 200 network objects. Each network object is a list
#' with three elements:
#' \describe{
#'    \item{edges}{the edgelist of the network. A two column
#'      \code{matrix} where each row is an edge.}
#'    \item{degree}{the degree sequence of the network, which is
#'      an \code{integer} vector of length n.}
#'    \item{n}{the network order. The order for every network is 2000.}
#' }


"artificial_networks"
