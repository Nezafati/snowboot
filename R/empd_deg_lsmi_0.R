empd_deg_lsmi_0 <- function(net, n.seeds, n.neigh, num_lsmi, seeds) {
  p0.seeds.array <- Oempd <- values.array <- val.seeds.array <- samples <- as.list(rep(NA, num_lsmi))
  seeds1 <- matrix(NA, num_lsmi, n.seeds)
  ## -------the 'real' parameters in the network:-------##
  real <- summary_net(net)
  realdd <- real$realdd
  ## ---------------------------------------------------##
  for (m in 1:num_lsmi) {
    # if(m%%100==1)#cat('Obtaining empd of sample ',m,'\n')
    if(is.null(seeds)){
      neigh.seeds <- sort(base::sample(1:length(net$degree), n.seeds, replace = FALSE))  #n.neigh=0!!!!!!!
    }
    else{
      neigh.seeds<-sort(as.vector(seeds))
    }
    tab.seeds <- table(neigh.seeds)  #id seeds
    seeds1[m, ] <- neigh.seeds
    ###### degrees #####
    deg.seeds <- realdd[rep(as.integer(names(tab.seeds)), tab.seeds)]  #duplicities allowed
    # --------------------#
    samples[[m]] <- list(freq.deg.seeds = freq.deg.seeds <- table(deg.seeds))
    ##### resample and extract the degree of selected vertices######
    values.array[[m]] <- values <- val.seeds.array[[m]] <- val.seeds <- sort(unique(deg.seeds))
    p0.seeds.array[[m]] <- sum(deg.seeds == 0)/n.seeds
    #### Frequency ##### (Not the relative frequency)
    OFseed <- table.row(deg.seeds, values)
    #### Empirical degree distribution ########
    Oempd.seeds <- OFseed/n.seeds
    Oempd[[m]] <- list(Oempd = Oempd.seeds)
  }  # for(m in 1:num_lsmi)
  # browser()
  list(samples = samples, values = values.array, Oempd = Oempd, num_lsmi = num_lsmi, val.seeds = val.seeds.array,
       n.seeds = n.seeds, n.neigh = n.neigh, p0.seeds = p0.seeds.array, seeds1 = seeds1)
}
