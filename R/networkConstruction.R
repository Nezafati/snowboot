# Copyright ? Lilia Leticia Ramirez Ramirez, llramirezramirez@gmail.com


################################################################### NETWORK CONSTRUCTION FUNCTIONS ###############################

library(VGAM)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

li <- function(delta, b, lim = 10000) {
  # b can be a vector
  res <- rep(0, length(b))
  if (length(delta) != 1) 
    stop("delta debe ser escalar") else res[b == 1] <- zeta(delta)
  # else if(length(delta)==length(b))res[b==1]<-zeta(delta[b==1])
  a <- which(b != 1)
  if (length(a)) 
    for (j in a) res[j] <- sum(b[j]^(1:lim)/(1:lim)^delta)
  res
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

resample <- function(x, size, rep = FALSE, prob = NULL) {
  # Suggested function to improve the sample function in R
  if (length(x) <= 1) {
    if (!missing(size) && size == 0) 
      x[FALSE] else x
  } else sample(x, size, replace = rep, prob = prob)
  # browser()
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

NDIM <- function(a, as.row = FALSE) {
  # I created NDIM,to use with either matrices or vectors takes a vector as column matrix if as.row=F and row matrix if
  # as.row=T
  if ((is.vector(a) | is.list(a)) & as.row) 
    n.dim <- c(1, length(a)) else n.dim <- c(NROW(a), NCOL(a))
  n.dim
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

vector.one <- function(n, q) {
  # function that returns a vector of zeros and ones. The ones are in positions given by q when q==numeric(0) r is a vector
  # of zeros
  r <- rep(0, n)
  r[q] <- 1
  r
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

elements.matrix <- function(matrix, esc, cases = 0, rows.mat = 0, as.row = FALSE, cols = FALSE) {
  #'matrix' can be a matrix or a vector representation of a matrix with rows=rows.mat.
  # If cols=TRUE This function returns: a) The (row,column) position of the elements in 'matrix' with value 'esc', or
  # b)transform the element 'case'-th in the vector 'matrix' into (row, column) for 'matrix' as a matrix if cols=FALSE only
  # returns the rows When 'matrix' is a vector, it can be considered column (as.row=F) or row (as.row=T) cases is the
  # position in the matrix of the elements equal to esc (by columns)
  if (cases == 0) 
    cases <- which(matrix == esc)
  if (rows.mat == 0) 
    rows.mat <- NDIM(matrix, as.row)[1]
  rows <- (cases - 1)%%rows.mat + 1
  result <- rows
  if (cols) {
    columns <- (cases - 1)%/%rows.mat + 1
    result <- list(rows = rows, cols = columns)
  }
  result
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

newrexp <- function(n, p) {
  # Function that deals with special cases for the function rexp (exponential disributed random numbers)
  if (is.null(p) | n == 0) 
    a <- Inf else {
    if (n != length(p)) 
      stop("dimensions do not match in newrexp")
    a <- rep(0, n)
    a[p > 0] <- rexp(length(p[p > 0]), p[p > 0])
    a[p == 0] <- Inf
  }
  a
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

order.edges <- function(edges, ord.col = TRUE) {
  # this function is to order the network matrix (or vector). Output is always a matrix order the network.edges so the rows
  # are increasing ord.col=T if we want the program to order columns as well, so in the first column appears the smaller
  # ID's
  if (is.matrix(edges)) {
    if (ord.col) {
      keep <- edges[, 1] < edges[, 2]
      edges[keep == 0, ] <- cbind(edges[keep == 0, 2], edges[keep == 0, 1])
    }
    edges <- edges[order(edges[, 1], edges[, 2]), ]
  } else edges <- t(as.matrix(sort(edges)))
  dimnames(edges) <- NULL
  edges
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

compare.vectors <- function(y, x) {
  # I define this function us use it with 'sapply'
  a <- sum(y == x)
  a
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

concatenate <- function(edges.hist, vec.mat, constant = NULL) {
  # Function to paste vec.mat (a matrix or vector) as new rows (or row) to a matrix edges.hist constant is a vector of
  # constants (usually: time,arc.width,color) edges.hist is data.frame vec.mat is not empty (is.null is not good)
  if (sum(NDIM(vec.mat) > 0) == 2) {
    if (is.vector(vec.mat)) 
      vec.mat <- t(as.matrix(vec.mat))  #vector (scalars are vectors) is now a matrix with one row 
    vec.mat <- cbind(vec.mat, matrix(rep(constant, nrow(vec.mat)), nrow(vec.mat), length(constant), byrow = TRUE))
    vec.mat <- as.data.frame(vec.mat)
    names(vec.mat) <- names(edges.hist)  #to be able to use rbind next
    edges.hist <- rbind(edges.hist, vec.mat)
  }
  edges.hist
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

grouping <- function(vec.mat, groups) {
  # Function that creates a column of numbers telling to which group each element (row) belongs to vec.mat is a matrix and
  # groups is a vector (numeric or character) with length=dim(vec.mat)[1] used to form the groups
  cbind(vec.mat, match(groups, groups))
}
# I define this function us use it with 'sapply' #Number of elements in column col that are equal to the value x.
number.elements <- function(mat, col, x) {
  length(mat[mat[, col] == x, col])
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

rztpois <- function(n, lambda) {
  # Function so simulate random numbers from: Zero-truncated Poisson distribution
  # P(N=k)=(lambda^k*exp^{-lambda})/(k!*(1-e^{-lambda})), for k=1,2,3,..
  maxval <- max(50, lambda + 10 * sqrt(lambda))
  cdf <- (cumsum(dpois(1:maxval, lambda)))/(1 - exp(-lambda))
  cut(runif(n), unique(c(0, cdf, 1)), labels = FALSE, right = FALSE)
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

rlog <- function(n, param) {
  # Function to generate random numbers from the Logarithmic distribution (poly-logarithmic when alpha=1) where param is a
  # positive real number.
  if (length(param) != 1 | param < 0) 
    stop("the Logarithmic parameter must be a positive real number") else {
    pdf.s <- (1:500)^(-1) * exp(-(1:500)/param)/(-log(1 - exp(-1/param)))
    sim <- cut(runif(n), unique(c(0, cumsum(pdf.s))), labels = FALSE, right = FALSE)
  }
  sim
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

rpower.law <- function(n, param) {
  # Function to generate random numbers from the power law k^{-param} where param a real number >1
  if (length(param) != 1 | param <= 1) 
    stop("the power law parameter must be a real number grater than one") else {
    pdf.s <- (1:1000)^(-param) * (zeta(param))^(-1)
    sim <- cut(runif(n), unique(c(0, cumsum(pdf.s))), labels = FALSE, right = FALSE)
  }
  sim
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

dpower.law <- function(x, param) {
  # Function to generate the distribution for the power law k^{-param}, where param is real number
  if (length(param) != 1) 
    stop("the power law parameter dimension is wrong (it must be of length 1)or the parameter value is incorrect") else pdf.s <- x^(-param) * (zeta(param))^(-1)
  pdf.s
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

rpoly.log <- function(n, param) {
  # Function to generate random numbers with the distribution Gutenberg-Richter law ck^{-delta}e^{-k/lambda}, with:
  # c=(li_lambda(e^{-1/lambda})) param=c(delta,lambda)
  if (length(param) != 2) 
    stop("the Gutenber-Richter law parameter is wrong") else sim <- cut(runif(n), unique(c(0, cumsum(dpoly.log(1:1000, param)))), labels = FALSE, right = FALSE)
  sim
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

dpoly.log <- function(x, param, lim = 10000) {
  # Function to generate the distribution for the Gutenberg-Richter law ck^{-delta}e^{-k/lambda}
  # c=(li_lambda(e^{-1/lambda})) param=c(delta,lambda)
  if (length(param) != 2) 
    stop("the Gutenber-Richter law parameter is wrong") else pdf.s <- x^(-param[1]) * exp(-x/param[2])/li(param[1], exp(-1/param[2]), lim = lim)
  pdf.s
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

sdegree <- function(n, distrib, param = 0) {
  # rutine to generate the degrees based on the specified distribution n is the number of indivudual to connect distrib is
  # the name of the distribution param is the unlist distribution parameter
  degree <- 3
  while (sum(degree)%%2 == 1) {
    # first step to realizable graph (Newman, Stogatz and Watts, 2001)
    if (distrib == "fixed") {
      if (length(param) != n | sum(param)%%2 == 1) 
        stop("the length of the degrees must be equal to the number of nodes and the sum of degree must to be even")
      degree <- param
    } else if (distrib == "pois") 
      degree <- rpois(n, param) else if (distrib == "ztpois") 
      degree <- rztpois(n, param) else if (distrib == "geom") 
      degree <- rgeom(n, 1/(param + 1))  #sometimes called exponential graph. param is then the expected value
 else if (distrib == "ztgeom") 
      degree <- rgeom(n, 1/(param + 1)) + 1 else if (distrib == "nbinom" && length(param) == 2) 
      degree <- rnbinom(size = param[1], prob = param[2]) else if (distrib == "poly.log") 
      degree <- rpoly.log(n, param) else if (distrib == "logarithmic") 
      degree <- rlog(n, param) else if (distrib == "power.law") 
      degree <- rpower.law(n, param) else if (distrib == "full") 
      degree <- rep(n - 1, n) else if (distrib == "none") 
      degree <- rep(0, n) else stop("incorrect distribution specification")
  }
  degree
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

ccheck <- function(nodes.dl, edges.dl) {
  # function that verify that if no other connexion is possible to construct a simple random graph and return a new edge if
  # possible. nodes.dl is the sorted nodes for which its degree left is still greater than zero. They are sorted edges.dl
  # is the subset of edges whose end points are elements of nodes.dl. They are sorted as in edges.
  other <- 0
  flag <- 1
  # browser()
  m <- length(nodes.dl)
  indexes <- cbind(unlist(mapply(rep, 1:(m - 1), (m - 1):1)), unlist(mapply(seq, 2:m, m)))
  a <- !is.element(paste(nodes.dl[indexes[, 1]], "-", nodes.dl[indexes[, 2]]), paste(edges.dl[, 1], "-", edges.dl[, 2]))
  if (any(a)) {
    new.edge <- sort(c(nodes.dl[indexes[a, 1]][1], nodes.dl[indexes[a, 2]][1]))
    other <- 1
  } else {
    new.edge <- NULL
    flag <- 0  #no more nodes can be connected
  }
  list(new.edge = new.edge, other = other, flag = flag)
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

edge.check <- function(edges, new.edge) {
  new.edge <- sort(new.edge)
  a <- which(edges == new.edge[1])
  b <- elements.matrix(edges, new.edge[2])
  if (length(a) == 0 | length(b) == 0) 
    comp <- 0 else comp <- sapply(a, compare.vectors, b)
  comp
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#



###################################################################### NETWORK CONSTRUCTION MAIN FUNCTION ######################


local.network.MR.new5 <- function(n, distrib, param = NULL, degree = NULL, take.p = 0.05, cutnet = 1e+05) {
  # param is not a list We develop the network construction by sampling in two stages. We use the frequency tables of the
  # vertices' degree and arch weights (that are the product of the vertices's degree they connect)
  
  # network creates the network based in the edge distribution and the number of nodes n number of individuals (susceptible
  # and infective) distrib is the degree distribution param is the distribution parameter (if the function is 'fixed' it is
  # a vector of degrees) distrib can be 'fixed' or 'pois' or 'ztpois' or 'geom' or 'nbinom','ztgeom' or 'poly.log' or
  # 'logarithmic' or 'power.law' or 'full' (fully connected) or 'none' (no element connected)
  
  # In this interative algorithm we eliminate the archs that cannot be longer be selected to construct a simple network In
  # order to optimize the algorithm we select more than one arch at each step
  
  if (!is.null(distrib) && is.null(degree)) {
    degree <- sdegree(n, distrib, param)
  } else {
    if (length(degree) != n) 
      stop("degree has to have length ", n, "\n")
  }
  id <- 1:n
  edges <- matrix(NA, ceiling(sum(degree)/2), 2)  #I want to reserve the memory for this variable using the maximum number of edges
  # edges1<-edges2<-rep(NA,ceiling(sum(degree)/2)) #I want to reserve the memory for this variable using the maximum number
  # of edges
  flag <- TRUE
  degree.left <- degree
  edge.row <- 1
  count <- 1
  
  while (sum(degree.left > 0) >= 2 & flag) {
    ## *********************************************************************##
    continue = TRUE
    eff.nodes <- (1:n)[degree.left > 0]
    # eff.degree<-degree.left[degree.left>0] it is the same as degree.left[eff.nodes] as noted next Now, if the network has
    # more than 1e4 I want to consider only cutnet nodes that I select by weight (degree)
    # if(sum(degree.left>0)>cutnet){eff.nodes<-sort(sample(eff.nodes,cutnet,prob=degree.left[eff.nodes]))} #I cut the number
    # of nodes here
    m <- length(eff.nodes)
    # pl<-m*(m-1)/2 #number of elements in products *********************************************************************##
    # a<-Sys.time()
    tab.degree <- table(degree.left[degree.left > 0])  #the frequency of degree.left that are greater than zero
    vals <- as.numeric(names(tab.degree))
    tab.prod.deg <- unlist(sapply(X = 1:(length(tab.degree)), FUN = table.mult.degree, vector = as.numeric(tab.degree), 
      m = length(tab.degree), freq = TRUE))
    vals.prod.deg <- unlist(sapply(X = 1:(length(tab.degree)), FUN = table.mult.degree, vector = vals, m = length(tab.degree), 
      freq = FALSE))
    take <- max(ceiling(length(tab.prod.deg) * take.p), 1)
    
    first.sample <- resample(1:length(tab.prod.deg), take, prob = (tab.prod.deg * vals.prod.deg), rep = TRUE)  #Here I sample the arch group
    # I want to allow to sample from the same group, but if in the second step I repeat a vertex, I remove the rest of
    # sampled groups and repeated vertices ahora tengo que poder traducir el sitio muestreado y el grupo de arcos con uno y
    # otro grados
    degree1 <- vals[id1 <- cut(first.sample, breaks = br <- cumsum(c(0, length(vals):1)), label = FALSE, include.lowest = TRUE)]
    degree2 <- vals[id2 <- id1 - 1 + first.sample - br[id1]]
    sam <- rep(NA, 2 * length(degree1))
    dd <- c(degree1, degree2)
    tdd <- table(dd)
    vtdd <- as.numeric(names(tdd))
    for (ss in 1:length(tdd)) {
      sam[dd == vtdd[ss]] <- resample(id[degree.left == vtdd[ss]], tdd[ss], rep = TRUE)
    }  #I had problems with the fuction sample
    tt <- which(duplicated(as.vector(rbind(sam1 <- sam[1:length(degree1)], sam2 <- sam[(length(degree1) + 1):(2 * length(degree1))]))))
    # tt<-which(duplicated(rbind(sam1<-sam[1:length(degree1)],sam2<-sam[(length(degree1)+1):(2*length(degree1))]))%%length(degree1))
    if (length(tt) > 0) {
      tt <- min(tt)
    } else {
      tt <- NULL
    }
    # now I remove all the edges from tt and up
    if (!is.null(tt)) {
      sam1 <- sam1[-c(ceiling(tt/2):length(degree1))]
      sam2 <- sam2[-c(ceiling(tt/2):length(degree1))]
    }
    if (length(sam1) > 0) {
      new.node1 <- new.node2 <- rep(NA, length(sam1))
      new.node1[sam1 <= sam2] <- sam1[sam1 <= sam2]
      new.node1[sam1 > sam2] <- sam2[sam1 > sam2]
      new.node2[sam1 <= sam2] <- sam2[sam1 <= sam2]
      new.node2[sam1 > sam2] <- sam1[sam1 > sam2]
      # Now new.node1 is always smaller than new.node2
    } else {
      continue = FALSE
      cat("entra \n")
    }
    
    if (continue) 
      {
        # if at least one edge candidate
        new.edge <- cbind(new.node1, new.node2)
        ## ****************************************************************## avoid to repeat new edges when whe have more to
        ## compare #edge.row ==1 for the first pair of nodes connected and >1 for the rest
        if (edge.row > 1 & length(new.node1) > 0) 
          {
          a <- which(is.element(edges[1:edge.row, 1], new.edge[, 1]))  #since new.edges sorted, a's who interest us are in the first column
          if (length(a) > 0) {
            b <- is.element(paste(new.edge[, 1], new.edge[, 2]), paste(edges[a, 1], edges[a, 2]))  #always in the second column
            if (any(b)) {
            # at least one repeated at leat one not repeated (we want to keep that one)
            if (any(!b)) {
              new.edge <- new.edge[!b, ]  #I eliminate the repeated
              if (is.vector(new.edge)) 
              new.edge <- t(new.edge)
            } else {
              new.edge <- NULL
            }  #all are repeated
            }
          }
          }  #if(edge.row>1 & dim(new.edge)[1]>0)
        ## ***************************************************************##
        if (!is.null(new.edge)) {
          edges[edge.row:(edge.row + dim(new.edge)[1] - 1), ] <- new.edge
          edge.row <- edge.row + dim(new.edge)[1]
          degree.left[new.edge[, 1]] <- degree.left[new.edge[, 1]] - 1
          degree.left[new.edge[, 2]] <- degree.left[new.edge[, 2]] - 1  #I can decrease only 1 in both cases because no node are repeated
        } else {
          if (m < 10) 
          count <- count + 1
        }  #few nodes left to connect but unable to find edges that are not repeated
        if (count > 10) 
          flag <- FALSE
        ### cat('productos y grados actualizados en ', Sys.time()-a,'\n')
        if (any(degree.left < 0)) {
          cat("raro\n")
          browser()
        }
      }  #continue
  }  #while
  
  edges <- edges[!is.na(edges[, 1]), ]  #remove NA's I put at the beginning
  edges <- order.edges(edges, ord.col = FALSE)
  list(edges = edges, degree = degree, degree.left = degree.left, n = n)
}  #Main function end
# Examples a<-local.network.MR.new5(1e4,'poly.log',c(2,13))

# tt<-local.network.MR.new4(n=1.2e2,distrib='pois', param=2, degree=NULL)
# real<-realized.degree(tt$n,tt$edges);all(real==tt$degree) plot(tt$degree,real)
# a<-which(tt$degree>real);real[a];tt$degree[a];tt$degree.left[a];sum(tt$degree.left)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#

table.mult.degree <- function(vector, place, m, freq = TRUE) {
  if (freq) {
    if (place < m) 
      res <- c(vector[place] * (vector[place] - 1)/2, vector[place] * vector[(place + 1):m]) else res <- vector[place] * (vector[place] - 1)/2
  } else {
    res <- vector[place] * vector[place:m]
  }
  res
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -#


###################################################################### NETWORK SAMPLING FUNCTION ###############################

sampleneigh <- function(net, n.seeds = 10, n.neigh = 1, seed = NULL) {
  # net is object network (what is important is the component $edges and the length of $degree) this function randomly
  # sample n.seeds o sample around the seeds specified by 'seed' and obtain for each, their n.neigh neighbourhood (set of
  # nodes with distance to the seed at most =n.neigh) also give the index of those nodes last added and for which we do not
  # have their complete degree information. seed0 are the original seed nodes are the no repeated elements in the sample
  # sampleN are the possibly repeated elements in the sample last.added are the vertices that are the most recently added
  # into the set.
  
  if (is.null(seed)) 
    {
      seed <- sort(sample(1:length(net$degree), n.seeds, rep = FALSE))
    }  #the seed selection is random
  sampleN <- seed0 <- seed
  # neighEdges<-NULL
  effEdges <- net$edges
  more <- TRUE
  nn <- n.neigh
  new.nodes <- 0
  if (n.neigh == 0) {
    # only keep the subgraph originated from the seeds
    a <- is.element(effEdges, seed)
    if (any(a)) {
      a <- which(matrix(a, dim(effEdges)[1], 2), arr.ind = TRUE)[, 1]  #now it is the row where they are in the edges matrix
      if (any(duplicated(a))) {
        # We keep the edges between nodes that are seeds and remove the rest. neighEdges<-effEdges[a[duplicated(a)],]
      }
    }
  } else {
    while (nn > 0 & more) {
      a <- is.element(effEdges, seed)  #'seed' will be accumulating all included vertices (non repeated)
      if (any(a)) 
        {
          eedges <- which(matrix(a, dim(effEdges)[1], 2), arr.ind = TRUE)
          a <- sort(eedges[, 1])  #now it is the row number where they are in the edges matrix
          arr.nodes <- sort(effEdges[cbind(eedges[, 1], sapply(eedges[, 2], FUN = switch, 2, 1))])  #the vertices we arrived to (duplicity is allowed)
          # I need this specially to know which vertices were the last added:
          if (!anyDuplicated(a)) {
          new.nodes <- arr.nodes  #all the vertices we arrive to weren't already included in 'seed'
          } else {
          new.nodes <- setdiff(arr.nodes, seed)
          }  #Then, already included vertices are not considered new because of inclusion of edge connecting them
          ### subEdges<-effEdges[unique(a),] #the subset of edges. The repeated just have to be included once.
          
          # maybe we are arriving to the nodes more than once (due to small cycles) or we can get again to already included
          # vertices (due to larger cycles).  We want to include them as may times as they are neighbours of already included
          # vertices. That is why I consider arr.nodes. if a originally seed vertex is included more than once, it is because it
          # was selected also by following one edge and then it also has the category of non seed.
          
          sampleN <- sort(c(sampleN, arr.nodes))
          seed <- unique(sampleN)
          
          if (nn > 1) 
          effEdges <- effEdges[-unique(a), ]  #I remove the 'used edges' to facilitate following searches within while,and assure
          # we do not 'arrive' to a node more times than edges it has.
          
          if (length(effEdges) > 0) {
          if (is.vector(effEdges)) 
            effEdges <- t(effEdges)  #I have to this very often in R to make sure effEdges is a matrix and not a vector 
          } else {
          more <- FALSE
          }  #when it reduces to become a matrix with one row.
          
        }  #end if(any(a))
      nn <- nn - 1
    }  #end while
  }
  list(seeds = seed0, nodes = seed, sampleN = sampleN, last.added = sort(new.nodes))  #neighEdges=neighEdges,
  # Examples #we are not really interested in running this function directly but within the next function called empdegree
  # distrib6 net<-local.network.MR.new5(n=100,distrib='pois',param=2) a<-sampleneigh(net,n.seeds=3,n.neigh=1,seed=NULL)
  # a<-sampleneigh(net,n.seeds=3,n.neigh=3,seed=NULL)
}
 
