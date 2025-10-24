library(igraph)

cliques_connected <- function(A,c1, c2){
  # assuming undirected edge
  cnxtns <- A[c1, c2]
  n.cnxtns <- sum(cnxtns > 0)
  if(n.cnxtns == length(A) | n.cnxtns == 0){
    return(F)
  } else {
    return(T)
  }
}

# checks if clique is connected to any other set
check_cliques <- function(G,cur_clique, cliques_found, min_clique_size){
  if(length(cur_clique) < min_clique_size){
    return(F)
  }

  connect.all = F
  if(length(cliques_found) > 0){
    num.connections <- c()
    n.seq <- c()
    for(cliq in cliques_found){
      num.connections <- c(num.connections,sum(G[cur_clique,cliq]))
      n.seq <- c(n.seq,length(G[cur_clique,cliq]))
    }
    if(any(num.connections > 0)){
      connect.all <- T
    }
  }
  return(connect.all)
}

# return the ego-graph of a node
ego_graph <- function(G,cur_node){

  idx <- which(G[,cur_node] > 0)
  idx.sort <- sort(c(idx,cur_node))
  G.sub <- G[idx.sort,idx.sort]
  out.list <- list("G" = G.sub, "idx" = idx.sort)
  return(out.list)
}

rand_ego_graph <- function(G,cur_node,max_subgraph_size){
  d <- colSums(G)
  full.idx <- which(G[,cur_node] > 0)
  d <- d[full.idx]
  idx <- sample(full.idx,max_subgraph_size,  replace = F, prob = d)
  idx.sort <- sort(c(idx,cur_node))
  G.sub <- G[idx.sort,idx.sort]
  out.list <- list("G" = G.sub, "idx" = idx.sort)
  return(out.list)
}



# adding an option for non-overlapping cliques
rand_clique_walk <- function(A,cur_node,
                             min_clique_size = 5,
                             num_iters = 1000,
                             num_cliques_stop = 50,
                             non_overlap_cliques = T,
                             anchored_reset = F){
  cliques_found = list()
  clique_node_set <- c()
  n.cliques <- 1
  nodes_visited = c()

  for(i in seq(num_iters)){
    cat(paste("Num Cliques:", length(cliques_found), "::: Step:", i, "/", num_iters), end = "\r")
    #print(length(cliques_found))


    if(length(cliques_found) >= num_cliques_stop){
      break
    }
    #cliques_found.extend([c for c in cliques_containing_node(G, cur_node) if check_cliques(G, c, cliques_found, min_clique_size)])
    if(!(cur_node %in% nodes_visited)){

      ego_subgraph = ego_graph(A, cur_node)
      A_sub <- ego_subgraph$G
      idx_sub <- ego_subgraph$idx #mapping

      g_sub <- igraph::graph_from_adjacency_matrix(A_sub, mode = "undirected")
      # cliques in the ego-graph
      cliques_sub <- igraph::max_cliques(g_sub, min = min_clique_size)
      cliques_found <- append(cliques_found,cliques_sub)
      for(clq in cliques_sub){
        main.idx <- idx_sub[clq]

        #no.overlap <- !(any(main.idx %in% clique_node_set))

        # connected.clique <- check_cliques(A, main.idx, cliques_found, min_clique_size) | length(cliques_found) == 0
        # if(non_overlap_cliques){
        #   proceed <- no.overlap & connected.clique
        # } else {
        #   proceed <- connected.clique
        # }
        # if(proceed){
        #   cliques_found[[n.cliques]] <- main.idx
        #   n.cliques <- n.cliques + 1
        #   clique_node_set <- c(clique_node_set, main.idx)
        #   break
        # }
      }
    }


    if(anchored_reset & length(clique_node_set) > 0){
      nodes_visited <- c(nodes_visited, cur_node)
      anchor_node <- sample(clique_node_set,1)
      neighbors = which(A[anchor_node,] > 0 )

      no_backtrack_neighbors <- neighbors[!neighbors %in% clique_node_set]
      cur_node = sample(no_backtrack_neighbors,1)
    } else {
      # no backtracking
      neighbors = which(A[cur_node,] > 0 )
      nodes_visited <- c(nodes_visited, cur_node)
      #no_backtrack_neighbors <- neighbors[!neighbors %in% nodes_visited]
      cur_node = sample(neighbors,1)
    }
  }

  return(cliques_found)
}


rand_clique_walk_2 <- function(A,cur_node,
                               min_clique_size = 10,
                               num_iters = 5000,
                               num_cliques_stop = 50,
                               non_overlap_cliques = T,
                               anchored_reset = F){
  cliques_found = list()
  clique_node_set <- c()
  #n.cliques <- 1
  nodes_visited = c()
  n_nodes = nrow(A)
  for(i in seq(num_iters)){
    cat(paste("Num Cliques:", length(cliques_found), "::: Step:", i, "/", num_iters), end = "\r")
    #print(length(cliques_found))


    if(length(cliques_found) >= num_cliques_stop){
      break
    }
    #cliques_found.extend([c for c in cliques_containing_node(G, cur_node) if check_cliques(G, c, cliques_found, min_clique_size)])
    if(!(cur_node %in% nodes_visited)){
      ego_subgraph = ego_graph(A, cur_node)
      A_sub <- ego_subgraph$G
      idx_sub <- ego_subgraph$idx #mapping
      g_sub <- igraph::graph_from_adjacency_matrix(A_sub, mode = "undirected")
      # cliques in the ego-graph

      #cliques_sub <- igraph::max_cliques(g_sub, min = min_clique_size)
      cliques_sub <- largest_cliques(g_sub)

      cliques_sub_re_idx <- list()

      if(length(cliques_sub) > 0){
        for(k in 1:length(cliques_sub)){
          cliques_sub_re_idx[[k]] <- idx_sub[cliques_sub[[k]]]
        }
      }
      for(k in 1:length(cliques_sub_re_idx)){
        if(length(cliques_sub_re_idx[[k]]) < min_clique_size){
          cliques_sub_re_idx[[k]] = NA
        }
      }
      cliques_sub_re_idx <- cliques_sub_re_idx[!is.na(cliques_sub_re_idx)]
      cliques_found <- append(cliques_found,cliques_sub_re_idx)
      if(length(cliques_found) > 0){
        cliques_found <- clique_partition(cliques_found,F)
      }


      clique_members <- unique(unlist(cliques_found))
    }


    # no backtracking
    neighbors = which(A[cur_node,] > 0 )
    idx1 <- !neighbors %in% nodes_visited
    idx2 <- !neighbors %in% clique_members
    neighbors = neighbors[idx1 & idx2]


    nodes_visited <- c(nodes_visited, cur_node)
    #no_backtrack_neighbors <- neighbors[!neighbors %in% nodes_visited]
    if(length(neighbors) > 0){
      cur_node = sample(neighbors,1)

    } else {
      print("reset node")
      if(length(clique_members) > 0){
        cur_node = sample(clique_members,1)
        idx = seq(n_nodes)
        p.idx = colSums(A)
        new.idx = ! idx %in% clique_members

        cur_node = sample(p.idx[new.idx],1, prob = p.idx[new.idx])
      } else {
        n_nodes <- nrow(A)

        cur_node = sample(seq(n_nodes),1, prob = colSums(A))
      }

    }


  }

  return(cliques_found)
}


rand_clique_walk_3 <- function(A,cur_node,
                               min_clique_size = 10,
                               num_iters = 5000,
                               num_cliques_stop = 50,
                               non_overlap_cliques = T,
                               anchored_reset = F,
                               max_subgraph_size = 1000){
  cliques_found = list()
  clique_node_set <- c()
  #n.cliques <- 1
  nodes_visited = c()
  n_nodes = nrow(A)
  for(i in seq(num_iters)){
    cat(paste("Num Cliques:", length(cliques_found), "::: Step:", i, "/", num_iters), end = "\r")
    #print(length(cliques_found))


    if(length(cliques_found) >= num_cliques_stop){
      break
    }
    #cliques_found.extend([c for c in cliques_containing_node(G, cur_node) if check_cliques(G, c, cliques_found, min_clique_size)])
    if(!(cur_node %in% nodes_visited)){
      ego_size <- sum(A[cur_node, ] > 0)
      if(ego_size < max_subgraph_size){
        ego_subgraph = ego_graph(A, cur_node)
        A_sub <- ego_subgraph$G
        idx_sub <- ego_subgraph$idx #mapping
      } else {
        ego_subgraph = rand_ego_graph(A, cur_node, max_subgraph_size)
        A_sub <- ego_subgraph$G
        idx_sub <- ego_subgraph$idx #mapping
      }


      g_sub <- igraph::graph_from_adjacency_matrix(A_sub, mode = "undirected")
      # cliques in the ego-graph

      #cliques_sub <- igraph::max_cliques(g_sub, min = min_clique_size)
      cliques_sub <- largest_cliques(g_sub)

      cliques_sub_re_idx <- list()

      if(length(cliques_sub) > 0){
        for(k in 1:length(cliques_sub)){
          cliques_sub_re_idx[[k]] <- idx_sub[cliques_sub[[k]]]
        }
      }
      for(k in 1:length(cliques_sub_re_idx)){
        if(length(cliques_sub_re_idx[[k]]) < min_clique_size){
          cliques_sub_re_idx[[k]] = NA
        }
      }
      cliques_sub_re_idx <- cliques_sub_re_idx[!is.na(cliques_sub_re_idx)]
      cliques_found <- append(cliques_found,cliques_sub_re_idx)
      if(length(cliques_found) > 0){
        cliques_found <- clique_partition(cliques_found,F)
      }


      clique_members <- unique(unlist(cliques_found))
    }


    # no backtracking
    neighbors = which(A[cur_node,] > 0 )
    idx1 <- !neighbors %in% nodes_visited
    idx2 <- !neighbors %in% clique_members
    neighbors = neighbors[idx1 & idx2]


    nodes_visited <- c(nodes_visited, cur_node)
    #no_backtrack_neighbors <- neighbors[!neighbors %in% nodes_visited]
    if(length(neighbors) > 0){
      cur_node = sample(neighbors,1)

    } else {
      print("reset node")
      if(length(clique_members) > 0){
        cur_node = sample(clique_members,1)
        idx = seq(n_nodes)
        p.idx = colSums(A)
        new.idx = ! idx %in% clique_members

        cur_node = sample(p.idx[new.idx],1, prob = p.idx[new.idx])
      } else {
        n_nodes <- nrow(A)

        cur_node = sample(seq(n_nodes),1, prob = colSums(A))
      }

    }


  }

  return(cliques_found)
}




clique_finder <- function(G, min_clique_size=10,
                          num_cliques_stop=100,
                          num_iters=1000,
                          non_overlap_cliques = T,
                          anchored_reset = T){

  n_nodes <- nrow(G)
  # more likely to start at high degree nodes
  #probs <- 1*(colSums(G) == max(colSums(G)))
  cur_node = sample(seq(n_nodes),1, prob = colSums(G))
  #print(cur_node)
  cliques_found = rand_clique_walk_2(G, cur_node,
                                     min_clique_size,
                                     num_iters, num_cliques_stop,
                                     non_overlap_cliques,anchored_reset)
  return(cliques_found)
}


clique_finder_2 <- function(G, min_clique_size=10,
                            num_cliques_stop=100,
                            num_iters=1000,non_overlap_cliques = T,
                            anchored_reset = T, max_subgraph_size = 2500){


  n_nodes <- nrow(G)
  # more likely to start at high degree nodes
  #probs <- 1*(colSums(G) == max(colSums(G)))
  cur_node = sample(seq(n_nodes),1, prob = colSums(G))
  #print(cur_node)
  cliques_found = rand_clique_walk_3(G, cur_node,
                                     min_clique_size,
                                     num_iters, num_cliques_stop,
                                     non_overlap_cliques,anchored_reset,
                                     max_subgraph_size)


  return(cliques_found)
}


# G = A.sim
# cluster_clique_search
cluster_clique_search <- function(G, res1 = 1, res2 = 0.9, min_clique_size = 8){
  g <- igraph::graph_from_adjacency_matrix(G, mode= "undirected")
  comm <- igraph::cluster_louvain(g, resolution = res1)
  #comm <- cluster_leading_eigen(g)
  print("Community Sizes")
  print(table(comm$membership))
  K <- max(comm$membership)
  clique.set <- list()
  for(k in seq(K)){
    cat(paste("Block:", k, "/", K), end = "\r")
    idx.sub <- which(comm$membership == k)

    G.sub <- G[idx.sub,idx.sub]
    if(length(G.sub) > 600^2){
      clique.sub <- cluster_clique_search_inner(G.sub,res = res2, min_clique_size = min_clique_size)
      # clique.sub <- clique_finder(G.sub, min_clique_size=min_clique_size,
      #                             num_cliques_stop=100,
      #                             num_iters=1000)

    } else {
      g.sub <- igraph::graph_from_adjacency_matrix(G.sub, mode= "undirected")

      clique.sub <- maximal.cliques(g.sub,min = min_clique_size)
    }


    if(length(clique.sub) > 0){
      cliques.found <- clique_partition(clique.sub,F)
      cliques.found.cor.idx <- list()
      for(s in seq(length(cliques.found))){
        cor.idx.labels <- idx.sub[cliques.found[[s]]]
        cliques.found.cor.idx <- append(cliques.found.cor.idx, list(cor.idx.labels))
      }
      clique.set <- append(clique.set, cliques.found.cor.idx)
    }
  }
  return(clique.set)
}


cluster_clique_search_inner <- function(G, res = 1.1, min_clique_size = 8){
  g <- igraph::graph_from_adjacency_matrix(G, mode= "undirected")
  comm <- igraph::cluster_louvain(g, resolution = res)
  #comm <- cluster_leading_eigen(g)
  print("Inner Community Sizes")
  print(as.numeric(table(comm$membership)))
  K <- max(comm$membership)
  clique.set <- list()
  for(k in seq(K)){
    cat(paste("Inner Block:", k, "/", K), end = "\r")
    idx.sub <- which(comm$membership == k)
    G.sub <- G[idx.sub,idx.sub]
    g.sub <- igraph::graph_from_adjacency_matrix(G.sub, mode= "undirected")

    clique.sub <- maximal.cliques(g.sub,min = min_clique_size)


    if(length(clique.sub) > 0){
      cliques.found <- clique_partition(clique.sub,F)
      cliques.found.cor.idx <- list()
      for(s in seq(length(cliques.found))){
        cor.idx.labels <- idx.sub[cliques.found[[s]]]
        cliques.found.cor.idx <- append(cliques.found.cor.idx, list(cor.idx.labels))
      }
      clique.set <- append(clique.set, cliques.found.cor.idx)
    }
  }
  return(clique.set)
}

#clique.set <- cluster_clique_search_2(G, res = 2, min_clique_size = ell)
cluster_clique_search_2 <- function(G, res = 2, min_clique_size = 8,max.iter = 30, max.cliques = 60, res.scale = 0.8){
  res.tmp = res
  n = nrow(G)
  clique.set <- list()
  iter = 1
  while(length(clique.set) < max.cliques & iter <= max.iter){
    clique.ids <- unlist(clique.set)
    idx.non.clique <- which(! 1:n  %in% clique.ids)
    G.sub <- G[idx.non.clique,idx.non.clique]

    g.sub <- igraph::graph_from_adjacency_matrix(G.sub, mode= "undirected")
    comm <- igraph::cluster_louvain(g.sub, resolution = res.tmp)
    rm(g.sub)
    #comm <- cluster_leading_eigen(g)
    print("Community Sizes")
    memb.counts <- as.numeric(table(comm$membership))
    print(memb.counts)
    K <- max(comm$membership)
    new.cliques <- list()
    for(k in seq(K)){
      cat(paste("Block:", k, "/", K), end = "\r")
      idx.block <- which(comm$membership == k)
      G.block <- G.sub[idx.block,idx.block]
      g.block <- igraph::graph_from_adjacency_matrix(G.block, mode= "undirected")

      clique.block <- maximal.cliques(g.block,min = min_clique_size)


      if(length(clique.block) > 0){
        cliques.found <- clique_partition(clique.block,F)
        cliques.found.cor.idx <- list()
        for(s in seq(length(cliques.found))){
          cor.idx.labels <- idx.block[cliques.found[[s]]]
          cliques.found.cor.idx <- append(cliques.found.cor.idx, list(cor.idx.labels))
        }
        new.cliques <- append(new.cliques, cliques.found.cor.idx)
      }
    }

    if(sum(memb.counts > min_clique_size) == 1){
      print("All Cliques Found")
      break
    }
    if(length(new.cliques) == 0 ){
      res.tmp = res.scale*res.tmp
      print("No new cliques found, coarsening clustering")
      print(paste("New Res:", round(res.tmp,3)))
    } else if(length(new.cliques) == 1 ){
      res.tmp = sqrt(res.scale)*res.tmp
      print("1 new clique found, coarsening clustering")
      print(paste("New Res:", round(res.tmp,3)))
    } else {
      if(max(memb.counts) < 1000){
        res.tmp = res.scale*res.tmp
      }

      new.cliques.re.idx <- list()
      for(k in seq(length(length(new.cliques) ))){
        new.cliques.re.idx[[k]] <- idx.non.clique[new.cliques[[k]]]
      }
      clique.set <- append(clique.set, new.cliques.re.idx)
    }

    iter = iter + 1
    cat(paste("Number of cliques found:", length(clique.set)), end = "\r")
    #sort( sapply(ls(),function(x){object.size(get(x))}))
  }
  if(iter >= max.iter){
    print("Max iterations")
  }

  return(clique.set)
}



cluster_clique_search_3 <- function(G, res = 2, min_clique_size = 8,max.iter = 30, max.cliques = 60, res.scale = 0.8){
  res.tmp = res
  n = nrow(G)
  clique.set <- list()
  iter = 1

  g <- igraph::graph_from_adjacency_matrix(G, mode= "undirected")
  comm <- igraph::cluster_louvain(g, resolution = res.tmp)
  comm <- igraph::cluster_leading_eigen(g)
  K <- max(comm$membership)
  clique.set <- list()
  for(k in seq(K)){
    print(paste("Block:", k,"/",K))
    idx.sub <- which(as.numeric(comm$membership) == k )
    G.sub <- G[idx.sub,idx.sub]
    if(length(idx.sub) > 1){
      n.sub <- nrow(G.sub)
      clique.set.sub <- list()
      clique.vec <- c()
      for(j in seq(n.sub)){
        if(j %% 200 == 0){
          cat(paste0(j,"/",n.sub), end = "\r")
        }

        idx.ego.sub <- which(G.sub[j,] > 0 & ! 1:n.sub %in% clique.vec )
        G.sub.ego <- G.sub[idx.ego.sub,idx.ego.sub]
        g.sub.ego <- igraph::graph_from_adjacency_matrix(G.sub.ego, mode= "undirected")
        cliques.set.ego <- maximal.cliques(g.sub.ego, min = min_clique_size)
        if(length(cliques.set.ego) > 0){
          cliques.set.ego <- clique_partition(cliques.set.ego)
          cliques.sub.re.idx <- list()
          R <- length(cliques.set.ego)
          for(r in seq(R)){
            cliques.sub.re.idx[[r]] <- idx.sub[idx.ego.sub[cliques.set.ego[[r]]]]
          }
          clique.set <- append(clique.set, cliques.sub.re.idx)
          clique.vec <- unlist(clique.set.sub)
        }

        clique.set <- clique_partition(clique.set)
      }
    }
  }
  return(clique.set)
}



cluster_clique_search_4 <- function(G, min_clique_size = 8, res = 1, verbose = F){
  res.tmp = res
  n = nrow(G)
  clique.set <- list()
  iter = 1

  g <- igraph::graph_from_adjacency_matrix(G, mode= "undirected")
  #comm <- igraph::cluster_louvain(g, resolution = res.tmp)
  comm <- igraph::cluster_leading_eigen(g)
  table(comm$membership)
  K <- max(comm$membership)
  clique.set <- list()
  for(k in seq(K)){
    #print(paste("Block:", k,"/",K))
    idx.sub <- which(as.numeric(comm$membership) == k )
    G.sub <- G[idx.sub,idx.sub]
    if(length(idx.sub) > 1){
      n.sub <- nrow(G.sub)
      g.sub <- igraph::graph_from_adjacency_matrix(G.sub, mode= "undirected")
      clique.set.sub <- igraph::maximal.cliques(g.sub, min = min_clique_size)
      if(length(clique.set.sub) > 0){
        clique.set.sub <- clique_partition(clique.set.sub)
        cliques.sub.re.idx <- list()
        R <- length(clique.set.sub)
        for(r in seq(R)){
          cliques.sub.re.idx[[r]] <- idx.sub[clique.set.sub[[r]]]
        }
        clique.set <- append(clique.set,cliques.sub.re.idx)
      }

    }
  }
  if(verbose){
    print(paste("Number of Cliques of size,",min_clique_size,":", length(clique.set)))
  }
  return(clique.set)
}




cluster_clique_search_5 <- function(G, min_clique_size = 8, res = 1, verbose = T, max_subgraph_size = 500){
  res.tmp = res
  n = nrow(G)
  clique.set <- list()
  iter = 1

  g <- igraph::graph_from_adjacency_matrix(G, mode= "undirected")
  comm <- igraph::cluster_louvain(g, resolution = res.tmp)
  K <- max(comm$membership)
  clique.set <- list()
  for(k in seq(K)){
    print(paste("Block:", k,"/",K))
    idx.sub <- which(as.numeric(comm$membership) == k )
    G.sub <- G[idx.sub,idx.sub]
    if(length(idx.sub) > 1){
      n.sub <- nrow(G.sub)
      g.sub <- igraph::graph_from_adjacency_matrix(G.sub, mode= "undirected")
      clique.set.sub <- clique_finder_2(G.sub, min_clique_size = min_clique_size,
                                        num_cliques_stop = 60, num_iters = 2000,
                                        max_subgraph_size = max_subgraph_size)
      if(length(clique.set.sub) > 0){
        clique.set.sub <- clique_partition(clique.set.sub)
        cliques.sub.re.idx <- list()
        R <- length(clique.set.sub)
        for(r in seq(R)){
          cliques.sub.re.idx[[r]] <- idx.sub[clique.set.sub[[r]]]
        }
        clique.set <- append(clique.set,cliques.sub.re.idx)
      }

    }
  }
  if(verbose){
    print(paste("Number of Cliques of size,",min_clique_size,":", length(clique.set)))
  }
  return(clique.set)
}


clique_set_connected_subgraph <- function(G, clique.set){
  K <- length(clique.set)
  A.sub <- matrix(0,K,K)

  for(k1 in seq(K)){
    for(k2 in seq(K)){
      A.sub[k1,k2] = 1*(sum(G[as.numeric(clique.set[[k1]]),as.numeric(clique.set[[k2]])]) > 0 )
    }
  }
  g.sub <- igraph::graph_from_adjacency_matrix(A.sub, mode = "undirected")
  connected.subgraph <- components(g.sub, mode = "strong")
  clust.idx <- which.max(table(connected.subgraph$membership))
  connect.idx <- which(connected.subgraph$membership == clust.idx)

  clique.subset <- clique.set[connect.idx]
  return(clique.subset)
}



