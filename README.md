# Fast graph-based algorithms for phase-type distributions

This software library provides fast and scalable algorithms for constructing and computing properties of for the statistical
distributions: continuous and discrete phase-type distributions, rewarded phase-type distributions, and the general multivariate phase-type distributions. The library can compute the moments (e.g. expectation, variance, covariance), the distribution function (pdf, cdf, pmf), the stopping time probability in the Markov jump process or Markov chain, and the distribution function of time-inhomogeneous phase-type distribution, as well as the expectation of rewarded time-inhomogeneous phase-type distributions.

This code is the basis for our paper: `Røikjer, Tobias, Asger Hobolth, and Kasper Munch. "Graph-based algorithms for phase-type distributions." Statistics and Computing 32.6 (2022): 103.`

The algorithms in the package are emperically multiple orders of magnitude faster than the traditional matrix-based equations for phase-type distributions of matrix inversion, matrix multiplication, and the matrix exponential. It can handle billions of vertices (states) and edges (transitions).

The code is written in C, but exposes an interface to both C and R through a C++ layer.


The state-space is built as a graph with an easy-to-use API. The state-space is stored efficiently as a graph, taking up very little space compared to using a matrix.The algorithms are fast, graph-based algorithms that operate of this data structure. The state-space can be converted back to a matrix-based format if desired.

The speed and memory improvements are often in the realm of 1000x faster than traditional matrix-based methods, if they are even possible to store as a matrix without running out of memory.

## Installing/compiling
For the R interface, simply install devtools and install from GitHub:

```R
install.packages("devtools",  INSTALL_opts = c('--no-lock'))
library(devtools)
devtools::install_github("TobiasRoikjer/PtDAlgorithms")
```

That is everything that is needed.

If you want to install the C library only, run e.g.:

```
cmake CMakeLists.txt
sudo make install
```

The library can also be compiled against R directly and easily!

## Source
The C and C++ header files can be found in `api`, and the source
code in `src/c` and `src/cpp`.

The R code can be seen in `src/ptdalgorithms.cpp` and the auto-generated code in `R`.

## Examples and documentation
Examples can be found in the `examples` directory. The full R api is documented as an easy-to-read jupyter notebook and can be seen at [https://github.com/TobiasRoikjer/PtDAlgorithms/blob/master/examples/full_api_example.ipynb](https://github.com/TobiasRoikjer/PtDAlgorithms/blob/master/examples/full_api_example.ipynb)

For the R interface, documentation can also be found in `man` which will automatically be installed when PtDAlgorithms is installed. Thereafter you can run e.g. `?ptdalgorithms` or `?ptdalgorithms::create_graph`.

A small example:

```R
library(ptdalgorithms)
# This is an example of constructing a state-space in R.
# The *construction* API (i.e. adding vertices and edges) is *slow* in R.
# If we need e.g. a hundred thousand+ vertices, the R API for construction
# should *not* be used. The C library should be used instead, and the graph
# produced can be used in R easily.

# State-space description:
# Rabbits stay in two islands, and jump between them. At each time, there is a probability
# of one of the islands being flooded, killing the rabbits. The phase-type distributions
# state-space models the time until no rabbits are left.

construct_rabbit_graph <- function(number_of_rabbits, flooding_rate_l, flooding_rate_r) {
    # We represent the vector as two integers, the number of rabbits on the left and right island
    state_vector_length <- 2
    graph <- create_graph(state_vector_length)
    initial_state <- c(number_of_rabbits, 0)
    # The initial state is the only starting state, with 100% starting probability
    add_edge(
      starting_vertex(graph),
      find_or_create_vertex(graph, initial_state),
      1
    )
    index <- 2
    # Iterate over all unvisited vertices
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

        # Left island flooding
        child_state <- c(0, state[2])
        add_edge(
          vertex,
          find_or_create_vertex(graph, child_state),
          flooding_rate_l
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
        # Right island flooding with rate of 4
        child_state <- c(state[1], 0)
        add_edge(
          vertex,
          find_or_create_vertex(graph, child_state),
          flooding_rate_r
        )
      }
      index <- index + 1
    }
    return(graph)
}
```

We can compute the expectation and variance (and all other moments)

```R
# Extremely fast to compute moments
graph <- construct_rabbit_graph(2, 2, 4)
expectation(graph)
variance(graph)
pdf <- dph(seq(0,5,by=0.01), graph)
cdf <- pph(seq(0,5,by=0.01), graph)

```

![https://github.com/TobiasRoikjer/PtDAlgorithms/blob/master/examples/graphic_rabbits.png?raw=true](https://github.com/TobiasRoikjer/PtDAlgorithms/blob/master/examples/graphic_rabbits.png?raw=true)

## Author
This software is written by Tobias Røikjer and released under the MIT license.
I hope this software is helpful,
