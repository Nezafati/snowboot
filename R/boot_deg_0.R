###################################################################################################
boot_deg_0 <- function(sam.out, num_lsmi, boot_rep) {
  # This function obtains the bootstrap samples for each sample from a network sam.out is the output of empd_deg_lsmi
  # num_lsmi is the number of different samples taken from the same network (Scalar or vector) boot_rep is the number of
  # bootstrap samples taken from each sample
  n.seeds <- sam.out$n.seeds
  n.neigh <- sam.out$n.neigh
  if (length(num_lsmi) == 1)
    num_lsmi <- 1:num_lsmi
  empd <- as.list(rep(NA, length(num_lsmi)))
  i <- 1
  for (m in num_lsmi) {
    val.seeds <- sam.out$val.seeds[[m]]
    freq.deg.seeds <- sam.out$samples[[m]]$freq.deg.seeds
    bsam.seeds <- myBsample(val.seeds, n.seeds, boot_rep, prob = freq.deg.seeds)
    values <- sam.out$values[[m]]  #all the possible degree values to resample
    #### Frequency ##### (Not the relative frequency)
    #check whether all seeds only have one degree
    if(length(val.seeds)==1){
      Fseed <- t(apply(bsam.seeds, 1, table.row, vect = values))  #freq (sorted according to values)
      Fseed <- t(Fseed)}
    else{
      Fseed <- t(apply(bsam.seeds, 1, table.row, vect = values))
    }
    # browser()
    empd.seeds <- Fseed/n.seeds
    empd[[i]] <- list(empd.seeds = empd.seeds)
    i <- i + 1
  }  # for(m in num_lsmi)
  list(values = sam.out$values[num_lsmi], empd = empd, num_lsmi = num_lsmi,
       boot_rep = boot_rep, n.neigh = n.neigh, seeds1 = sam.out$seeds1[num_lsmi,])
}#safe
##########################################################
