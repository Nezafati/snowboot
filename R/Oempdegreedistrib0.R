Oempdegreedistrib0 <- function(net, n.seed, n.neigh, num.sam, idname = "Temp", seeds) {
      p0.seed.array <- Oempd <- values.array <- val.seed.array <- samples <- as.list(rep(NA, num.sam))
      seeds1 <- matrix(NA, num.sam, n.seed)
      ## -------the 'real' parameters in the network:-------##
      real <- real.parameters(net)
      realdd <- real$realdd
      ## ---------------------------------------------------##
      for (m in 1:num.sam) {
            # if(m%%100==1)#cat('Obtaining empd of sample ',m,'\n')
            neigh.seeds <- sort(sample(1:length(net$degree), n.seed, rep = FALSE))  #n.neigh=0!!!!!!!
            tab.seeds <- table(neigh.seeds)  #id seeds
            seeds1[m, ] <- neigh.seeds
            ###### degrees #####
            deg.seed <- realdd[rep(as.integer(names(tab.seeds)), tab.seeds)]  #duplicities allowed
            # --------------------#
            samples[[m]] <- list(freq.deg.seed = freq.deg.seed <- table(deg.seed))
            ##### resample and extract the degree of selected vertices######
            values.array[[m]] <- values <- val.seed.array[[m]] <- val.seed <- sort(unique(deg.seed))
            p0.seed.array[[m]] <- sum(deg.seed == 0)/n.seed
            #### Frequency ##### (Not the relative frequency)
            OFseed <- table.row(deg.seed, values)
            #### Empirical degree distribution ########
            Oempd.seed <- OFseed/n.seed
            Oempd[[m]] <- list(Oempd = Oempd.seed)
      }  # for(m in 1:num.sam)
      # browser()
      list(idname = idname, samples = samples, values = values.array, Oempd = Oempd, num.sam = num.sam, val.seed = val.seed.array,
           n.seed = n.seed, n.neigh = n.neigh, p0.seed = p0.seed.array, seeds1 = seeds1)
}
