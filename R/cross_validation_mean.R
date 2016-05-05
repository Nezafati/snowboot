B.EmpDistrib <-function(net, n.seeds, n.neigh, sam.size=1, n.boot,
                        method = "w", otherNetParameters=FALSE){
  #sam.size (==1 for LSMI1) is the number of different samples taken from the network for each i and j
  #otherNetParameters is true if intervals and fallins for the rest of the parmeters
  #  (other than mean) are required.
  Obs.distrib.out <- empd <- as.list(rep(NA, length(n.seeds)*length(n.neigh)))
  counter <- 1
  for(i in n.seeds){
    for(j in n.neigh){

      if(j==0){
        Obs.distrib<-Oempdegreedistrib(net, n.seeds=i, n.neigh=j, num.sam=sam.size)
        TMP <- Obs.distrib$seeds1
      }else{
        Obs.distrib<-Oempdegreedistrib(net, n.seeds=i,n.neigh=j, num.sam=sam.size, seeds=NULL)
      }
      Oparam<-OparametersEst(Obs.distrib)
      #B.distrib<-bootdeg(Obs.distrib, num.sam=sam.size,n.boot=n.boot)
      #return(B.distrib)
      #browser()
      tmp <- bootdeg(Obs.distrib, num.sam=sam.size,n.boot=n.boot, method = method)$empd[[1]]
      Obs.distrib.out[[counter]] <-Obs.distrib
      empd[[counter]] <- tmp
      counter <- counter+1
    }
    #return(B.distrib)
  }
  return(list(Obs.distrib.out=Obs.distrib.out, empd = empd))
}


combineLSMINodes <- function(bootEmpD){
  #function combines the unodes(or seed1) from the elements of
  #Obs.distrib.out object, which is inside the bootEmpD list
  nodes <- NULL
  for(i in 1:length(bootEmpD$Obs.distrib.out)){
    if("nodes_of_LSMI"%in%names(bootEmpD$Obs.distrib.out[[i]])){
      tmp=bootEmpD$Obs.distrib.out[[i]]$nodes_of_LSMI
    } else tmp=bootEmpD$Obs.distrib.out[[i]]$seeds1
    nodes=c(nodes,tmp)
  }
  unlist(nodes)
}
closestCoverNDX <- function(x, alpha=alpha){
  coverage <- 1-alpha
  which(abs(x-coverage)==min(abs(x-coverage)),arr.ind=T)
}

# A function that sorts a matrix of tied optimal seed-wave combinations so that
# the first is the largest number of seeds among the smallest waves.
#
# The function takes a matrix containing tied optimal seed-wave combinations
# from the training proxy and sorts the rows of this matrix so that the
# largest number of seeds for the lowest number of waves is on top
# @param inMat is a matrix of tied seed-wave indices.
# @return A matrix with rows such that the indices for the largest number of
# seeds (first column) for the lowest number of waves (second column) is on top.
# @rdname sort_tied_opti
sort_tied_opti <- function(inMat){
  if(dim(inMat)[1] == 1)
    return(inMat)
  smallestSeedFirst <- inMat[order(inMat[,1],decreasing = F),]
  outMat <- smallestSeedFirst[order(smallestSeedFirst[,2]),]
  outMat
}
# Alternative is to get the smallest sample size which is estimated with mean(degree)
# ##1.2. Find approximate sample size
# SS <- matrix(NA, length(n.neigh), length(n.seeds))
# SS[1,] <- n.seeds
# for(i in 1:5){
#       SS[(i+1),] <- SS[i,] + SS[1,]*MU*(MU-1)^(i-1)
# }
# #SS

#' A function that uses cross-validation to select seed-wave combination for
#' estimation of a degree's frequency.
#'
#' The function's inputs are a network, a vector of possible seed sample-sizes,
#' a vector of possible waves, and a few tuning parameters. The output will
#' contain the best seed-wave combination for each degree and the width of the
#' 95 percent bootstrap confidence intervals at each degree for
#' the best seed-wave combination.
#' @note Only one LSMI per seed-wave combination is currently supported.
#' @references Efron, B. (1979). Bootstrap methods: another look at the
#'  jackknife. The annals of Statistics, 1-26.
#' @references Thompson, M. E., Ramirez Ramirez, L. L., Lyubchich, V. and
#'  Gel, Y. R. (2015), Using the bootstrap for statistical inference
#'  on random graphs. Can J Statistics. doi: 10.1002/cjs.11271
#' @param network A network object that is list containing:
#'  \describe{
#'    \item{edges}{The edgelist of the network. A two column
#'      \code{matrix} where each row is an edge.}
#'    \item{degree}{The degree sequence of the network, which is
#'      an \code{integer} vector of length n.}
#'    \item{n}{The network order.}
#'  }
#'    The object can be created by \code{\link{local.network.MR.new5}} or
#'    it can be imported.
#' @param n.seeds A numeric vector for the different sample sizes of seed to use
#'  in cross-validation.
#' @param n.neigh A numeric vector for the different waves to use
#'  in cross-validation.
#' @param n.boot The number of bootstrap sample.
#' @param method Can be either "w" for weighted bootstrap or "nw" for
#'    non-weighted bootstrap. "w" is recommended and set as the default method.
#' @param proxyRep The number of time to sample a proxy. Default is 19.
#' @param proxyOrder The size of the proxy sample. Default is 30.
#' @return A list consisting of
#'  \item{selected_seed_wave}{A matrices that provides
#'    the best seed-wave combinations (obtained via cross-validation) for
#'    the respective estimation method.}
#'  \item{selected_seed_wave}{A vector of length 2 that provides
#'    the bootstrap confidence intervals for the estimated mean degree
#'    using the best seed-wave combinations (see above).}
#' @export
#' @examples
#' net <- artificial_networks[[1]]
#' a <- cross_validation_mean(network = net, n.seeds = c(10, 20, 30), n.neigh = c(1, 2),
#'  n.boot = 200, method = "w")

cross_validation_mean <- function(network, n.seeds, n.neigh, n.boot,
                             proxyRep = 19, proxyOrder = 30, method = "w", alpha = .05){
  sam.size = 1
  n.seeds <- sort(n.seeds)
  n.neigh <- sort(n.neigh)
  net_order <- network$n
  #make bootEmpD list for seed-wave combos
  bootEmpD <- B.EmpDistrib(net = network, n.seeds = n.seeds, n.neigh = n.neigh,
                          sam.size = sam.size, n.boot = n.boot, method = method)
  used <- unique(combineLSMINodes(bootEmpD))
  count <- 1
  fallin.proxy <- array(0, c(length(n.seeds), length(n.neigh), proxyRep))
  for(i in 1:length(n.seeds)){
      # i=1
    for(j in 1:length(n.neigh)){
        # j=1
        # build proxy from bootEmpD$Obs.empd.out
        tmp <- bootEmpD$empd[[count]][[1]]
        values <- bootEmpD$Obs.distrib.out[[count]]$values[[1]]
        est_means <- rowSums(tmp*rep(values,each=n.boot))
        bootCI_mean <- quantile(est_means, c((alpha/2), 1-(alpha/2)))
        for(k in 1:proxyRep){
          # k=1
          proxyNodes <- sample(used, proxyOrder, replace = F)
          proxy_mean <- mean(network$degree[proxyNodes])
          fallin.proxy[i, j, k] <- ((bootCI_mean[1]<proxy_mean) & (proxy_mean<bootCI_mean[2]))
        }
        count <- count+1
    }
  }

  coverage.proxy <- apply(fallin.proxy, c(1, 2), mean, na.rm=T)

  #output Matrices
  opti.CI.w.p0s <- opti.CI.nw.p0sEkb <- opti.CI.nw.p0sEks <-
      matrix(nrow = 2, ncol = length(estimable_k_from_boot),
             dimnames = list(c("LB", "UB"), estimable_k_from_boot))

    opti.seed_wave.w.p0s <- opti.seed_wave.nw.p0sEkb <- opti.seed_wave.nw.p0sEks <-
      matrix(nrow = 2, ncol = length(estimable_k_from_boot),
             dimnames = list(c("Seeds", "Waves"), estimable_k_from_boot))
    #test if matrix or list
    if (is.list(opti.cover.w.p0s)){

      for(kDegree in estimable_k_from_boot){
        kDegree_num <- as.numeric(kDegree)
        sortedOpti.cover.w.p0s <- sort_tied_opti(opti.cover.w.p0s[[kDegree_num]])
        #browser()
        optimalSeedNDX <-  sortedOpti.cover.w.p0s[1, ][1]
        optimalNeighNDX <- sortedOpti.cover.w.p0s[1, ][2]
        opti.seed_wave.w.p0s[1, kDegree] <- n.seeds[optimalSeedNDX]
        opti.seed_wave.w.p0s[2, kDegree] <- n.neigh[optimalNeighNDX]
        listLocation <- (optimalSeedNDX-1)*length(n.neigh)+optimalNeighNDX
        opti.CI.w.p0s[,kDegree] <- stats::quantile(bootEmpD$w.p0s[[listLocation]][, kDegree], probs=c(0.025, 0.975))

      }
    }else if(is.matrix(opti.cover.w.p0s)){
      for(kDegree in estimable_k_from_boot){
        kDegree_num <- as.numeric(kDegree)
        optimalSeedNDX = opti.cover.w.p0s[1, kDegree_num]
        optimalNeighNDX <- opti.cover.w.p0s[2, kDegree_num]
        opti.seed_wave.w.p0s[1, kDegree] <- n.seeds[optimalSeedNDX]
        opti.seed_wave.w.p0s[2, kDegree] <- n.neigh[optimalNeighNDX]
        listLocation <- (optimalSeedNDX-1)*length(n.neigh)+optimalNeighNDX
        opti.CI.w.p0s[,kDegree] <- stats::quantile(bootEmpD$w.p0s[[listLocation]][, kDegree], probs=c(0.025, 0.975))
        }
    }else print("unknown data type output from optimal function")

    CI_selected_seed_wave <- list(opti.CI.w.p0s = opti.CI.w.p0s,
                                  opti.CI.nw.p0sEkb = opti.CI.nw.p0sEkb,
                                  opti.CI.nw.p0sEks = opti.CI.nw.p0sEks)
    selected_seed_wave <-
      list(opti.seed_wave.w.p0s = opti.seed_wave.w.p0s,
           opti.seed_wave.nw.p0sEkb = opti.seed_wave.nw.p0sEkb,
             opti.seed_wave.nw.p0sEks = opti.seed_wave.nw.p0sEks)
    res <- list(selected_seed_wave = selected_seed_wave,
                CI_selected_seed_wave = CI_selected_seed_wave)
    res
}
