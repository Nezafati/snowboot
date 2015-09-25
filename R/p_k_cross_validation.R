B.EmpDistrib <-function(net,n.seeds,n.neigh,sam.size=1,n.boot,otherNetParameters=FALSE){
    #sam.size (==1 for LSMI1) is the number of different samples taken from the network for each i and j
    #otherNetParameters is true if intervals and fallins for the rest of the parmeters
    #  (other than mean) are required.
    Obs.distrib.out <- w.p0s <- nw.p0sEkb <- nw.p0sEks <- as.list(rep(NA, length(n.seeds)*length(n.neigh)))
    counter <- 1
    for(i in n.seeds){
        for(j in n.neigh){

            if(j==0){n.dist<-1   #n.dist is the number of different emp distr.
            }else{n.dist<-3}

            if(j==0){
                Obs.distrib<-Oempdegreedistrib(net,n.seed=i,n.neigh=j,num.sam=sam.size)
                TMP <- Obs.distrib$seeds1
            }else{
                Obs.distrib<-Oempdegreedistrib(net,n.seed=i,n.neigh=j,num.sam=sam.size, seed=NULL)
            }
            Oparam<-OparametersEst(Obs.distrib)
            #B.distrib<-Bempdegreedistrib(Obs.distrib, num.sam=sam.size,n.boot=n.boot)
            #return(B.distrib)
            #browser()
            tmp <- Bempdegreedistrib(Obs.distrib, num.sam=sam.size,n.boot=n.boot)$empd[[1]]
            Obs.distrib.out[[counter]] <-Obs.distrib
            w.p0s[[counter]] <- tmp$empd.nw.p0sEkb
            nw.p0sEkb[[counter]] <- tmp$empd.nw.p0sEks
            nw.p0sEks[[counter]] <- tmp$empd.w.p0s
            counter <- counter+1
        }
        #return(B.distrib)
    }
    return(list(Obs.distrib.out=Obs.distrib.out, w.p0s=w.p0s, nw.p0sEkb=nw.p0sEkb, nw.p0sEks=nw.p0sEks))
}

combineLSMINodes <- function(bootEmpD){
    #function combines the unodes(or seed1) from the elements of
    #Obs.distrib.out object, which is inside the bootEmpD list
    if("unodes"%in%names(bootEmpD$Obs.distrib.out[[1]])){
        nodes=bootEmpD$Obs.distrib.out[[1]]$unodes
    } else nodes=bootEmpD$Obs.distrib.out[[1]]$seeds1


    for(i in 2:length(bootEmpD$Obs.distrib.out)){
        if("unodes"%in%names(bootEmpD$Obs.distrib.out[[i]])){
            tmp=bootEmpD$Obs.distrib.out[[i]]$unodes
        } else tmp=bootEmpD$Obs.distrib.out[[i]]$seeds1
        nodes=c(nodes,tmp)
    }
    unlist(nodes)
}
closestCoverNDX <- function(x,coverage=.95){
    which(abs(x-coverage)==min(abs(x-coverage)),arr.ind=T)}
p_k_cross_validation <- function(networks, n.seeds, n.neigh, n.boot, sam.size=1,
                                 kmax){
      # Add more inputs for kmax upperbound and lowerbound for number of times
      # each p_k is sampled after all proxyReps are done. kmax upperbound should
      # be 1/4 the max k in LSMI
      for (mc in 1:MC){#mc=1
            bootEmpD=B.EmpDistrib(networks[[mc]],n.seeds,n.neigh,sam.size,n.boot)
            used <- unique(combineLSMINodes(bootEmpD))
            count <- 1
            for(i in 1:length(n.seeds)){
                  # i=1
                  for(j in 1:length(n.neigh)){
                        # j=1
                        # build proxy from bootEmpD$Obs.empd.out
                        tmp.w.p0s <- apply(bootEmpD$w.p0s[[count]], 2, quantile, probs=c(0.025, 0.975))
                        tmp.nw.p0sEkb <- apply(bootEmpD$nw.p0sEkb[[count]], 2, quantile, probs=c(0.025, 0.975))
                        tmp.nw.p0sEks <- apply(bootEmpD$nw.p0sEks[[count]], 2, quantile, probs=c(0.025, 0.975))
                        confidenceIntervalWidth.w.p0s[i, j, mc, ] <- tmp.w.p0s[2, 1:kmax]-tmp.w.p0s[1, 1:kmax]
                        confidenceIntervalWidth.nw.p0sEkb[i, j, mc, ] <- tmp.nw.p0sEkb[2,1:kmax]-tmp.w.p0s[1, 1:kmax]
                        confidenceIntervalWidth.nw.p0sEks[i, j, mc, ] <- tmp.nw.p0sEks[2,1:kmax]-tmp.w.p0s[1, 1:kmax]
                        for(k in 1:proxyRep){
                              # k=1
                              proxyNodes <- sample(used, proxyOrder, replace = F)
                              proxyDegrees <- networks[[mc]]$degree[proxyNodes]
                              proxyPMF <- table(proxyDegrees)/proxyOrder
                              #trasform into vector with "0" where no obs occur
                              PMFvector <- as.vector(proxyPMF[match(1:kmax, dimnames(proxyPMF)[[1]], nomatch=NA)])
                              PMFvector[is.na(PMFvector)] <- 0
                              #take first 5
                              firstK <- PMFvector[1:kmax]

                              fallin.proxy.w.p0s[i, j, k, ] <- ((tmp.w.p0s[1, 1:kmax]<firstK) & (firstK<tmp.w.p0s[2, 1:kmax]))
                              fallin.proxy.nw.p0sEkb[i, j, k, ] <- ((tmp.nw.p0sEkb[1, 1:kmax]<firstK) & (firstK<tmp.nw.p0sEkb[2, 1:kmax]))
                              fallin.proxy.nw.p0sEks[i, j, k, ] <- ((tmp.nw.p0sEks[1, 1:kmax]<firstK) & (firstK<tmp.nw.p0sEks[2, 1:kmax]))
                        }
                        count <- count+1
                  }
            }
            coverage.proxy.w.p0s <- apply(fallin.proxy.w.p0s, c(1, 2, 4), mean, na.rm=T)
            coverage.proxy.nw.p0sEkb <- apply(fallin.proxy.nw.p0sEkb, c(1, 2, 4), mean, na.rm=T)
            coverage.proxy.nw.p0sEks <- apply(fallin.proxy.nw.p0sEks, c(1, 2, 4), mean, na.rm=T)

            opti.cover.w.p0s <- apply(coverage.proxy.w.p0s, 3, closestCoverNDX)
            opti.cover.nw.p0sEkb <- apply(coverage.proxy.nw.p0sEkb, 3, closestCoverNDX)
            opti.cover.nw.p0sEks <- apply(coverage.proxy.nw.p0sEks, 3, closestCoverNDX)
            #test if matrix or list

            if (is.list(opti.cover.w.p0s)){

                  for(kDegree in 1:kmax){
                        optimalSeedNDX = opti.cover.w.p0s[[kDegree]][1, ][1]
                        optimalNeighNDX <- opti.cover.w.p0s[[kDegree]][1, ][2]
                        listLocation <- (optimalSeedNDX-1)*length(n.neigh)+optimalNeighNDX
                        CI <- quantile(bootEmpD$w.p0s[[listLocation]][, kDegree], probs=c(0.025, 0.975))
                        fallin.trudDist.w.p0s[optimalSeedNDX, optimalNeighNDX,kDegree,2]=fallin.trudDist.w.p0s[optimalSeedNDX, optimalNeighNDX, kDegree, 2]+1
                        if((CI[1]<trueDist[kDegree]) & (trueDist[kDegree]<CI[2])){
                              fallin.trudDist.w.p0s[optimalSeedNDX, optimalNeighNDX, kDegree, 1]=fallin.trudDist.w.p0s[optimalSeedNDX, optimalNeighNDX, kDegree, 1]+1
                        }
                  }
            }else if(is.matrix(opti.cover.w.p0s)){
                  for(kDegree in 1:kmax){
                        optimalSeedNDX = opti.cover.w.p0s[1, kDegree]
                        optimalNeighNDX <- opti.cover.w.p0s[2, kDegree]
                        listLocation <- (optimalSeedNDX-1)*length(n.neigh)+optimalNeighNDX
                        CI <- quantile(bootEmpD$w.p0s[[listLocation]][, kDegree], probs=c(0.025, 0.975))
                        fallin.trudDist.w.p0s[optimalSeedNDX, optimalNeighNDX, kDegree, 2]=fallin.trudDist.w.p0s[optimalSeedNDX, optimalNeighNDX, kDegree, 2]+1
                        if((CI[1]<trueDist[kDegree]) & (trueDist[kDegree]<CI[2])){
                              fallin.trudDist.w.p0s[optimalSeedNDX, optimalNeighNDX, kDegree, 1]=fallin.trudDist.w.p0s[optimalSeedNDX, optimalNeighNDX, kDegree, 1]+1
                        }
                  }
            }else print("unknown data type output from optimal function")




            if (is.list(opti.cover.nw.p0sEkb)){

                  for(kDegree in 1:kmax){
                        optimalSeedNDX = opti.cover.nw.p0sEkb[[kDegree]][1,][1]
                        optimalNeighNDX <- opti.cover.nw.p0sEkb[[kDegree]][1,][2]
                        listLocation <- (optimalSeedNDX-1)*length(n.neigh)+optimalNeighNDX
                        CI <- quantile(bootEmpD$nw.p0sEkb[[listLocation]][,kDegree], probs=c(0.025, 0.975))
                        fallin.trudDist.nw.p0sEkb[optimalSeedNDX, optimalNeighNDX,kDegree,2]=fallin.trudDist.nw.p0sEkb[optimalSeedNDX, optimalNeighNDX,kDegree,2]+1
                        if((CI[1]<trueDist[kDegree]) & (trueDist[kDegree]<CI[2])){
                              fallin.trudDist.nw.p0sEkb[optimalSeedNDX, optimalNeighNDX,kDegree,1]=fallin.trudDist.nw.p0sEkb[optimalSeedNDX, optimalNeighNDX,kDegree,1]+1
                        }
                  }
            }else if(is.matrix(opti.cover.nw.p0sEkb)){
                  for(kDegree in 1:kmax){
                        optimalSeedNDX = opti.cover.nw.p0sEkb[1,kDegree]
                        optimalNeighNDX <- opti.cover.nw.p0sEkb[2,kDegree]
                        listLocation <- (optimalSeedNDX-1)*length(n.neigh)+optimalNeighNDX
                        CI <- quantile(bootEmpD$nw.p0sEkb[[listLocation]][,kDegree], probs=c(0.025, 0.975))
                        fallin.trudDist.nw.p0sEkb[optimalSeedNDX, optimalNeighNDX,kDegree,2]=fallin.trudDist.nw.p0sEkb[optimalSeedNDX, optimalNeighNDX,kDegree,2]+1
                        if((CI[1]<trueDist[kDegree]) & (trueDist[kDegree]<CI[2])){
                              fallin.trudDist.nw.p0sEkb[optimalSeedNDX, optimalNeighNDX,kDegree,1]=fallin.trudDist.nw.p0sEkb[optimalSeedNDX, optimalNeighNDX,kDegree,1]+1
                        }
                  }
            }else print("unknown data type output from optimal function")


            if (is.list(opti.cover.nw.p0sEks)){

                  for(kDegree in 1:kmax){
                        optimalSeedNDX = opti.cover.nw.p0sEks[[kDegree]][1,][1]
                        optimalNeighNDX <- opti.cover.nw.p0sEks[[kDegree]][1,][2]
                        listLocation <- (optimalSeedNDX-1)*length(n.neigh)+optimalNeighNDX
                        CI <- quantile(bootEmpD$nw.p0sEks[[listLocation]][,kDegree], probs=c(0.025, 0.975))
                        fallin.trudDist.nw.p0sEks[optimalSeedNDX, optimalNeighNDX,kDegree,2]=fallin.trudDist.nw.p0sEks[optimalSeedNDX, optimalNeighNDX,kDegree,2]+1
                        if((CI[1]<trueDist[kDegree]) & (trueDist[kDegree]<CI[2])){
                              fallin.trudDist.nw.p0sEks[optimalSeedNDX, optimalNeighNDX,kDegree,1]=fallin.trudDist.nw.p0sEks[optimalSeedNDX, optimalNeighNDX,kDegree,1]+1
                        }
                  }
            }else if(is.matrix(opti.cover.nw.p0sEks)){
                  for(kDegree in 1:kmax){
                        optimalSeedNDX = opti.cover.nw.p0sEks[1,kDegree]
                        optimalNeighNDX <- opti.cover.nw.p0sEks[2,kDegree]
                        listLocation <- (optimalSeedNDX-1)*length(n.neigh)+optimalNeighNDX
                        CI <- quantile(bootEmpD$nw.p0sEks[[listLocation]][,kDegree], probs=c(0.025, 0.975))
                        fallin.trudDist.nw.p0sEks[optimalSeedNDX, optimalNeighNDX,kDegree,2]=fallin.trudDist.nw.p0sEks[optimalSeedNDX, optimalNeighNDX,kDegree,2]+1
                        if((CI[1]<trueDist[kDegree]) & (trueDist[kDegree]<CI[2])){
                              fallin.trudDist.nw.p0sEks[optimalSeedNDX, optimalNeighNDX,kDegree,1]=fallin.trudDist.nw.p0sEks[optimalSeedNDX, optimalNeighNDX,kDegree,1]+1
                        }
                  }
            }else print("unknown data type output from optimal function")



      }



      crossValidationResults <- list(fallin.trudDist.w.p0s,fallin.trudDist.nw.p0sEkb,
                                     fallin.trudDist.nw.p0sEks,confidenceIntervalWidth.w.p0s,
                                     confidenceIntervalWidth.nw.p0sEkb,confidenceIntervalWidth.nw.p0sEks)

      #save(crossValidationResults,file=paste("polylog(0.1,2)1000x2000_kmax30_CVresults-",jobNDX,".RData",sep = ""))
}


