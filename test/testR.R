setwd("~/school/PtDAlgorithms/test")
Rcpp::sourceCpp('../src/ptdalgorithms.cpp')

n <- 10
graph <- create_graph(n)
index <- 1

while (index <= vertices_length(graph)) {
  vertex <- vertex_at(graph, index)
  
  if (vertex$vertex == starting_vertex(graph)$vertex) {
    start <- create_vertex(graph, c(n, rep(0, n-1)))
    add_edge(starting_vertex(graph), start, 1)
    index <- index + 1
    next()
  }
  
  for (i in 1:n) {
    for (j in i:n) {
      state <- vertex$state
      
      if (i == j) {
        if (state[i] < 2) {
          next;
        }
        
        rate <- state[i] * (state[i] - 1) / 2
      } else {
        if (state[i] < 1 || state[j] < 1) {
          next;
        }
        
        rate <- state[i] * state[j]
      }
      
      child_state <- state
      child_state[i] <- child_state[i] - 1
      child_state[j] <- child_state[j] - 1
      child_state[i+j] <- child_state[i+j] + 1
      
      child <- find_or_create_vertex(graph, child_state)
      
      add_edge(vertex, child, rate)
    }
  }
  
  index <- index + 1
}

ptd <- graph_as_matrix(graph)

library(expm)

ptd$IPV %*% solve(-ptd$SIM) %*% rep(1, length(ptd$IPV))
stopifnot(as.integer(ptd$IPV %*% solve(-ptd$SIM) %*% diag(ptd$states[,1]) %*% rep(1, length(ptd$IPV))) == 2)
stopifnot(as.integer(ptd$IPV %*% solve(-ptd$SIM) %*% diag(ptd$states[,2]) %*% rep(1, length(ptd$IPV))) == 1)
rw3 <- as.numeric(ptd$IPV %*% solve(-ptd$SIM) %*% diag(ptd$states[,3]) %*% rep(1, length(ptd$IPV)))
rw3_o <- sum(expected_waiting_time(graph)*sapply(vertices(graph), function(v) {v$state[3]}))
stopifnot(abs(rw3-rw3_o) < 0.0001)

reward_transform(graph, sapply(vertices(graph), function(v) {v$state[3]}))
ptd2 <- graph_as_matrix(graph)

stopifnot(abs(as.numeric(ptd2$IPV %*% solve(-ptd2$SIM) %*% rep(1, length(ptd2$IPV)))-rw3) < 0.0001)
stopifnot(graph_is_acyclic(graph))
add_edge(vertex_at(graph, 5), vertex_at(graph, 3), 4)
stopifnot(!graph_is_acyclic(graph))
ptd3 <- graph_as_matrix(graph)
e <- as.numeric(ptd3$IPV %*% solve(-ptd3$SIM) %*% rep(1, length(ptd3$IPV)))
e2 <- sum(expected_waiting_time(graph))
stopifnot(abs(e-e2) < 0.0001)
