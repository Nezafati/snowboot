igraph_to_network <- function(in_graph){
      if(igraph::is.directed(in_graph))
            stop('Please import undirected igraph object, e.g. output of
                 graph_from_edgelist(el,directed = F)')
      edges <- igraph::as_edgelist(in_graph)
      #edges <- order.edges(edges, ord.col = FALSE)
      degree <- igraph::degree(in_graph)
      degree.left <- rep(0,igraph::gorder(in_graph))
      n <- igraph::gorder(in_graph)
      res <- list(edges = edges, degree = degree,
                  degree.left = degree.left, n = n)
      res
}
