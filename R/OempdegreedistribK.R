OempdegreedistribK <- function(net, n.seed, n.neigh, num.sam, idname = "Temp", seeds) {
      # This function obtains the empirical degree distribution from num.sam samples net is the network (only one) n.seed is
      # the number of seed to set the neighbourhood sample n.neigh is the neighbouhood size around each seed num.sam is the
      # number of different samples taken from the same network idname is to identify from which nets we are sampling and
      # resampling.
      seeds1 <- matrix(0, num.sam, n.seed)
      p0.seed.array <- Oempd <- ekseed.array <- values.array <- val.seed.array <- val.nonseed.array <- samples <- as.list(rep(NA,
                                                                                                                              num.sam))

      ## -------the 'real' parameters in the network:-------##
      real <- real.parameters(net)
      realdd <- real$realdd
      # rmeand<-real$mean(realdd) rquart<-real$rquart
      rfreq <- real$rfreq
      # rperc<-real$rperc -------------------------------------##

      for (m in 1:num.sam) {
            # if(m%%100==1)#cat('Obtaining empd of sample ',m,'\n') browser()
            neigh <- sampleneighSequential(net, n.seeds = n.seed, n.neigh = n.neigh, seed = seeds[m, ])
            seeds1[m, ] <- neigh$seeds
            # nodes<-neigh$sampleN[!is.element(neigh$sampleN,neigh$last.added)] #vertices that are not included last (with their
            # duplicities)
            nodes <- neigh$sampleN  #vertices up to distance n.neigh(with their duplicities)!!!!!!!!!!!!!!
            tab.nodes <- table(nodes)  #now it has the info of seeds and non-seed
            # Now we want to distinguish between seeds and non seeds. Remember that some seeds can also be non seed if they were also
            # selected by following one edge.
            tab.seeds <- table(neigh$seeds)  #id seeds
            a <- is.element(names(tab.nodes), names(tab.seeds))
            if (all(names(tab.nodes[a]) == names(tab.seeds))) {
                  tab.nodes[a] <- tab.nodes[a] - tab.seeds
            } else {
                  cat("no mismo orden")
                  browser()
            }
            tab.nodes <- tab.nodes[tab.nodes > 0]  #now I have only the id of non-seeds
            ###################### degrees #####
            deg.seed <- realdd[rep(as.integer(names(tab.seeds)), tab.seeds)]  #in case seed are present more than once as seeds
            deg.nonseed <- realdd[rep(as.integer(names(tab.nodes)), tab.nodes)]  #to incorporate their duplicity
            deg.nonseedU <- realdd[as.integer(names(tab.nodes))]  #it contains the non seeds only one time (I do not use in the rest of the code)
            # --------------------#
            samples[[m]] <- list(freq.deg.seed = freq.deg.seed <- table(deg.seed), freq.deg.nonseed = freq.deg.nonseed <- table(deg.nonseed),
                                 freq.deg.nonseedU = freq.deg.nonseedU <- table(deg.nonseedU))
            ##### resample and extract the degree of selected vertices######
            val.seed.array[[m]] <- val.seed <- sort(unique(deg.seed))  #as.numeric(names(freq.deg.seed))
            val.nonseed.array[[m]] <- val.nonseed <- sort(unique(deg.nonseed))  #as.numeric(names(freq.deg.nonseed))
            val.nonseedU <- sort(unique(deg.nonseedU))  #as.numeric(names(freq.deg.nonseedU))

            p0.real <- rfreq[1]
            p0.seed <- 0
            if (any(val.seed == 0)) {
                  # if any seed has degree zero
                  p0.seed <- sum(deg.seed == 0)/n.seed
            }
            p0.seed.array[[m]] <- p0.seed
            values <- sort(union(val.seed, val.nonseed))  #all the possible degree values toresample
            values.array[[m]] <- values

            ###################### Frequency ##### (Not the relative frequency)
            OFseed <- table.row(deg.seed, values)
            OFnonseed <- table.row(deg.nonseed, values)
            #################################################################### combining information from seeds and nonseeds ###### mean degree computed from the original sampled seeds:
            ekseed.array[[m]] <- ekseed <- sum(as.numeric(names(freq.deg.seed)) * freq.deg.seed)/n.seed
            colzero <- NULL
            if (any(values == 0)) {
                  colzero <- which(values == 0)
                  vals <- values[-colzero]
                  Of.seed <- OFseed[-colzero]
                  Of.nonseed.nw <- OFnonseed[-colzero]
            } else {
                  vals <- values
                  Of.seed <- OFseed
                  Of.nonseed.nw <- OFnonseed
            }
            ################################################ NWB # seeds and non weighted nonseeds #### p0 estimated from orginal sampled seeds#
            Oempd.nw.p0sEks <- (Of.seed + (1 - p0.seed) * ekseed * Of.nonseed.nw/vals)/(n.seed + ekseed * sum(Of.nonseed.nw/vals))
            ######################################### p0 taken as known #
            Oempd.nw.p0rEks <- (Of.seed + (1 - p0.real) * ekseed * Of.nonseed.nw/vals)/(n.seed + ekseed * sum(Of.nonseed.nw/vals))
            if (any(values == 0)) {
                  Oempd.nw.p0sEks <- c(p0.seed, Oempd.nw.p0sEks)  #5
                  Oempd.nw.p0rEks <- c(p0.real, Oempd.nw.p0rEks)  #8
            }
            Oempd[[m]] <- list(Oempd = Oempd.nw.p0sEks, Oempd.nw.p0rEks = Oempd.nw.p0rEks)
      }  # for(m in 1:num.sam)
      # browser()
      list(idname = idname, samples = samples, values = values.array, Oempd = Oempd, num.sam = num.sam, val.seed = val.seed.array,
           val.nonseed = val.nonseed.array, n.seed = n.seed, n.neigh = n.neigh, p0.real = p0.real, p0.seed = p0.seed.array,
           ekseed = ekseed.array, seeds1 = seeds1)
}
