setwd("~/school/PtDAlgorithms/test")
Rcpp::sourceCpp('../src/ptdalgorithms.cpp')
devtools::install_github("TobiasRoikjer/PtDalgorithms")
library(ptdalgorithms)

# Build Kingman graph
n <- 4
graph <- create_graph(n)
starting_state <- c(n, rep(0, n-1))
start <- create_vertex(graph, starting_state)
add_edge(starting_vertex(graph), start, 1)
index <- 2

while (index <= vertices_length(graph)) {
  vertex <- vertex_at(graph, index)
  
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

# Convert to matrix based
ptd <- graph_as_matrix(graph)


# Expected singleton and doubleton length
stopifnot(as.integer(ptd$IPV %*% solve(-ptd$SIM) %*% diag(ptd$states[,1]) %*% rep(1, length(ptd$IPV))) == 2)
stopifnot(as.integer(ptd$IPV %*% solve(-ptd$SIM) %*% diag(ptd$states[,2]) %*% rep(1, length(ptd$IPV))) == 1)
rw3 <- as.numeric(ptd$IPV %*% solve(-ptd$SIM) %*% diag(ptd$states[,3]) %*% rep(1, length(ptd$IPV)))
rw3_o <- sum(expected_waiting_time(graph)*sapply(vertices(graph), function(v) {v$state[3]}))
# Graph based vs matrix
stopifnot(abs(rw3-rw3_o) < 0.0001)

# Reward transform graph
reward_transform(graph, sapply(vertices(graph), function(v) {v$state[3]}))
ptd2 <- graph_as_matrix(graph)

stopifnot(abs(as.numeric(ptd2$IPV %*% solve(-ptd2$SIM) %*% rep(1, length(ptd2$IPV)))-rw3) < 0.0001)


# Make cyclic
stopifnot(graph_is_acyclic(graph))
add_edge(vertex_at(graph, 5), vertex_at(graph, 3), 4)
stopifnot(!graph_is_acyclic(graph))

# Must be able to compute expected waiting time
ptd3 <- graph_as_matrix(graph)
e <- as.numeric(ptd3$IPV %*% solve(-ptd3$SIM) %*% rep(1, length(ptd3$IPV)))
e2 <- sum(expected_waiting_time(graph))
stopifnot(abs(e-e2) < 0.0001)

# Compute new rewards that are equivalent to 2nd moment when taken as
# the expected value
v <- as.numeric(ptd3$IPV %*% solve(-ptd3$SIM)%*% solve(-ptd3$SIM) %*% rep(1, length(ptd3$IPV)))
rw <- moment_rewards(graph, rep(1, vertices_length(graph)))
v2 <- sum(expected_waiting_time(graph)*rw)

stopifnot(abs(v-v2) < 0.0001)
reward_transform(graph, rw)

v3 <- sum(expected_waiting_time(graph))
stopifnot(abs(v-v3) < 0.0001)



#test:
state_vector_length <- 2 
graph <- create_graph(state_vector_length)
starting_vertex <- vertex_at(graph, 1)
initial_state <- c(2, 0)
add_edge(
  starting_vertex,
  create_vertex(graph, initial_state),
  1
)
index <- 2

while (index <= vertices_length(graph)) {
  vertex <- vertex_at(graph, index)
  state <- vertex$state
  
  if (state[1] > 0) {
    # Rabbit jump left to right
    child_state <- c(state[1] - 1, state[2] + 1)
    add_edge(
      vertex,
      find_or_create_vertex(graph, child_state),
      1
    )
    
    # Island flooding
    child_state <- c(0, state[2])
    add_edge(
      vertex,
      find_or_create_vertex(graph, child_state),
      2
    )
  }
  
  if (state[2] > 0) {
    # Rabbit jump right to left
    child_state <- c(state[1] + 1, state[2] - 1)
    add_edge(
      vertex,
      find_or_create_vertex(graph, child_state),
      1
    )
    
    # Island flooding
    child_state <- c(state[1], 0)
    add_edge(
      vertex,
      find_or_create_vertex(graph, child_state),
      4
    )
  }
  
  index <- index + 1
}
print(graph_as_matrix(graph))
reward_transform(graph,sapply(vertices(graph), function(v) {v$state[2]}))
print(graph_as_matrix(graph))

# TODO below
M <- graph_as_matrix(graph)
R <- exp_desc_rewards(graph, rep(1, vertices_length(graph)))
RM <- solve(-M$SIM)%*%rep(1, length(M$IPV))
RR <- exp_desc_rewards(graph, R)
RMR <- solve(-M$SIM)%*%solve(-M$SIM)%*%rep(1, length(M$IPV))
stopifnot(abs(M$IPV %*% RM - R[1]) < 0.01)
stopifnot(abs(M$IPV %*% RMR - RR[1]) < 0.01)



state_vector_length <- 2
graph <- create_graph(state_vector_length)
starting_vertex <- vertex_at(graph, 1)
initial_state <- c(100, 0)
add_edge(
  starting_vertex,
  create_vertex(graph, initial_state),
  1
)
index <- 2

while (index <= vertices_length(graph)) {
  vertex <- vertex_at(graph, index)
  state <- vertex$state
  
  if (state[1] > 0) {
    # Rabbit jump left to right
    child_state <- c(state[1] - 1, state[2] + 1)
    add_edge(
      vertex,
      find_or_create_vertex(graph, child_state),
      1
    )
    
    # Island flooding
    child_state <- c(0, state[2])
    add_edge(
      vertex,
      find_or_create_vertex(graph, child_state),
      1
    )
  }
  
  if (state[2] > 0) {
    # Rabbit jump right to left
    child_state <- c(state[1] + 1, state[2] - 1)
    add_edge(
      vertex,
      find_or_create_vertex(graph, child_state),
      1
    )
    
    # Island flooding
    child_state <- c(state[1], 0)
    add_edge(
      vertex,
      find_or_create_vertex(graph, child_state),
      1
    )
  }
  
  index <- index + 1
}

print(vertices_length(graph))

a <- sapply(vertices(graph), function(v) {(v$state)[2] > 100})
reward_transform(graph, a)
print(vertices_length(graph))


M <- graph_as_matrix(graph)
R <- exp_desc_rewards(graph, rep(1, vertices_length(graph)))
RM <- solve(-M$SIM)%*%rep(1, length(M$IPV))
RR <- exp_desc_rewards(graph, R)
RMR <- solve(-M$SIM)%*%solve(-M$SIM)%*%rep(1, length(M$IPV))
stopifnot(abs(M$IPV %*% RM - R[1]) < 0.01)
stopifnot(abs(M$IPV %*% RMR - RR[1]) < 0.01)


df <- data.frame()
library(dplyr)

for (L in c(2, 10, 20, 50, 100, 150, 200, 350, 400, 450, 500)) {
  print(L)
  tim <- proc.time()
state_vector_length <- 2
graph <- create_graph(state_vector_length)
starting_vertex <- vertex_at(graph, 1)
initial_state <- c(L, 0)
add_edge(
  starting_vertex,
  create_vertex(graph, initial_state),
  1
)
index <- 2

while (index <= vertices_length(graph)) {
  vertex <- vertex_at(graph, index)
  state <- vertex$state
  
  if (state[1] > 0) {
    # Rabbit jump left to right
    child_state <- c(state[1] - 1, state[2] + 1)
    add_edge(
      vertex,
      find_or_create_vertex(graph, child_state),
      1
    )
    
    # Island flooding
    child_state <- c(0, state[2])
    add_edge(
      vertex,
      find_or_create_vertex(graph, child_state),
      1
    )
  }
  
  if (state[2] > 0) {
    # Rabbit jump right to left
    child_state <- c(state[1] + 1, state[2] - 1)
    add_edge(
      vertex,
      find_or_create_vertex(graph, child_state),
      1
    )
    
    # Island flooding
    child_state <- c(state[1], 0)
    add_edge(
      vertex,
      find_or_create_vertex(graph, child_state),
      1
    )
  }
  
  index <- index + 1
}

print(vertices_length(graph))
k <- (proc.time() - tim)[3]
expected_waiting_time(graph, rep(1, vertices_length(graph)))[1]
df <- bind_rows(df, data.frame(rabbits=L, vertices=vertices_length(graph), construction_time = k, time_taken = (proc.time() - tim)[3]  - k))
print(df)

}

library(ggplot2)

ggplot(df, aes(y=time_taken)) +
  geom_point(aes(x=vertices, color="Vertices")) +
  geom_line(aes(x=vertices, color="Vertices"))
  #geom_point(aes(x=rabbits, color="Rabbits"))+
  #scale_x_continuous(sec.axis = sec_axis(~.*5, name = "Rabbits"))+
  scale_colour_manual(values = c("blue", "red"))
  
  
  
  state_vector_length <- 2
  graph <- create_graph(state_vector_length)
  starting_vertex <- vertex_at(graph, 1)
  initial_state <- c(400, 0)
  add_edge(
    starting_vertex,
    create_vertex(graph, initial_state),
    1
  )
  index <- 2
  
  while (index <= vertices_length(graph)) {
    vertex <- vertex_at(graph, index)
    state <- vertex$state
    
    if (state[1] > 0) {
      # Rabbit jump left to right
      child_state <- c(state[1] - 1, state[2] + 1)
      add_edge(
        vertex,
        find_or_create_vertex(graph, child_state),
        1
      )
      
      # Island flooding
      child_state <- c(0, state[2])
      add_edge(
        vertex,
        find_or_create_vertex(graph, child_state),
        1
      )
    }
    
    if (state[2] > 0) {
      # Rabbit jump right to left
      child_state <- c(state[1] + 1, state[2] - 1)
      add_edge(
        vertex,
        find_or_create_vertex(graph, child_state),
        1
      )
      
      # Island flooding
      child_state <- c(state[1], 0)
      add_edge(
        vertex,
        find_or_create_vertex(graph, child_state),
        1
      )
    }
    
    index <- index + 1
  }
  
  print(vertices_length(graph))

tim <- proc.time()

expected_waiting_time(graph)[1]

cat((proc.time() - tim)[3])
cat("\n")

# Parameterized

tim <- proc.time()
state_vector_length <- 2
graph <- create_graph(state_vector_length)
starting_vertex <- vertex_at(graph, 1)
initial_state <- c(400, 0)
add_edge(
  starting_vertex,
  create_vertex(graph, initial_state),
  1
)
index <- 2

while (index <= vertices_length(graph)) {
  vertex <- vertex_at(graph, index)
  state <- vertex$state
  
  if (state[1] > 0) {
    # Rabbit jump left to right
    child_state <- c(state[1] - 1, state[2] + 1)
    add_edge(
      vertex,
      find_or_create_vertex(graph, child_state),
      1
    )
    
    # Island flooding
    child_state <- c(0, state[2])
    add_edge(
      vertex,
      find_or_create_vertex(graph, child_state),
      0,
      c(1)
    )
  }
  
  if (state[2] > 0) {
    # Rabbit jump right to left
    child_state <- c(state[1] + 1, state[2] - 1)
    add_edge(
      vertex,
      find_or_create_vertex(graph, child_state),
      1
    )
    
    # Island flooding
    child_state <- c(state[1], 0)
    add_edge(
      vertex,
      find_or_create_vertex(graph, child_state),
      0,
      c(2)
    )
  }
  
  index <- index + 1
}
cat((proc.time() - tim)[3])
cat("\n")

tim <- proc.time()
expected_waiting_time(graph)[1]

cat((proc.time() - tim)[3])
cat("\n")

tim <- proc.time()
graph_update_weights_parameterized(graph, c(20))

expected_waiting_time(graph)[1]

cat((proc.time() - tim)[3])
cat("\n")

tim <- proc.time()

expected_waiting_time(graph)[1]

cat((proc.time() - tim)[3])
cat("\n")

tim <- proc.time()
graph_update_weights_parameterized(graph, c(10))
expected_waiting_time(graph)[1]
NULL
