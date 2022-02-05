/*
 * MIT License
 *
 * Copyright (c) 2021 Tobias Røikjer
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/* While it seems very strange to have this in a C file, the R code
* has very strange linking behavior, and we therefore sometimes include
* the same C file...
*/

#ifndef PTDALGORITHMS_PTDCPPR_C
#define PTDALGORITHMS_PTDCPPR_C

#define PTD_RCPP 1
#define PTD_INTEGRATE_EXCEPTIONS 1

#ifndef PTD_DEBUG_1_INDEX
#define PTD_DEBUG_1_INDEX 1
#endif

#include <climits>

#include "ptdalgorithms_types.h"
#include "../api/c/ptdalgorithms.h"
#include "../api/cpp/ptdalgorithmscpp.h"

#include "c/ptdalgorithms.c"

#ifndef SOURCE_CPP
#include "cpp/ptdalgorithmscpp.cpp"
#endif

#include <stdlib.h>

#include <Rcpp.h>
using namespace Rcpp;
using namespace ptdalgorithms;

List _graph_as_matrix(SEXP phase_type_graph);
SEXP get_first_list_entry(SEXP e, std::string message);
List vertex_as_list(Vertex *vertex);
SEXP list_as_vertex(SEXP list);

static void set_c_seed(){
  Function f("runif");
  
  NumericVector res = f(1, 0, 1000000);
  
  srand(res[0]);
}

static int fac(int n) {
  if (n==0) {
    return 1;
  }
  
  return n * fac(n - 1);
}

//' Create a graph representing a phase-type distribution
//' 
//' @description
//' `create_graph` creates a graph representing a phase-type distribution.
//' This is the primary entry-point of the library.
//' 
//' @details
//' There will *always* be a starting vertex added to
//' the graph.
//' 
//' Notice that when the library functions are invoked on
//' this object, the object is *mutated*, i.e. changed, which
//' may be surprising considering the normal behavior of R
//' objects.
//' 
//' @return Simple reference to a CPP object.
//' 
//' @param state_length The length of the integer vector used to represent and reference a state.
//' 
//' @examples
//' graph <- create_graph(4)
// [[Rcpp::export]]
SEXP create_graph(size_t state_length) {
  return Rcpp::XPtr<Graph>(
    new Graph(
        state_length
    )
  );
}

//' Find or create a vertex matching `state`
//' 
//' @description
//' Finds a vertex by the `state` parameter. If no such
//' vertex exists, it creates the vertex and adds it to
//' the graph object instead.
//' 
//' @details
//' A faster and simpler version of calling [ptdalgorithms::find_vertex()] and  [ptdalgorithms::create_vertex()]
//' 
//' @return The newly found or inserted vertex in the graph
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param state An integer vector of what vertex to look for. Has length as given by `state_length` in  [ptdalgorithms::create_graph()]
//' 
//' @examples
//' graph <- create_graph(4)
//' find_or_create_vertex(graph, c(1,2,1,0)) # Adds and returns the vertex
//' find_or_create_vertex(graph, c(1,2,1,0)) # Only returns the vertex
//' # `graph` is now changed permanently
// [[Rcpp::export]]
List find_or_create_vertex(SEXP phase_type_graph, IntegerVector state) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  Vertex *found = graph->find_or_create_vertex_p(as<std::vector<int> >(state));
  
  return vertex_as_list(found);
}


//' Adds an edge between two vertices in the graph
//' 
//' @description
//' The graph represents transitions between states as
//' a weighted direction edge between two vertices.
//' 
//' @seealso [ptdalgorithms::expected_waiting_time()]
//' @seealso [ptdalgorithms::moments()]
//' @seealso [ptdalgorithms::variance()]
//' @seealso [ptdalgorithms::covariance()]
//' @seealso [ptdalgorithms::graph_update_weights_parameterized()]
//' 
//' @param phase_type_vertex_from The vertex that transitions from
//' @param phase_type_vertex_to The vertex that transitions to
//' @param weight The weight of the edge, i.e. the transition rate
//' @param parameterized_edge_state Optional. Associate a numeric vector to an edge, for faster computations of moments when weights are changed.
//' 
//' @examples
//' graph <- create_graph(4)
//' vertex_a <- find_or_create_vertex(graph, c(1,2,1,0))
//' vertex_b <- find_or_create_vertex(graph, c(2,0,1,0))
//' add_edge(vertex_a, vertex_b, 1.5)
// [[Rcpp::export]]
void add_edge(SEXP phase_type_vertex_from, SEXP phase_type_vertex_to, double weight, NumericVector parameterized_edge_state = NumericVector::create()) {
  phase_type_vertex_from = list_as_vertex(phase_type_vertex_from);
  phase_type_vertex_to = list_as_vertex(phase_type_vertex_to);
  
  Rcpp::XPtr<Vertex> from(phase_type_vertex_from);
  Rcpp::XPtr<Vertex> to(phase_type_vertex_to);
  
  if (parameterized_edge_state.length() == 0) {
    from->add_edge(*to.get(), weight);
  } else {
    from->add_edge_parameterized(*to.get(), weight, as<std::vector<double> >(parameterized_edge_state));
  }
}


//' Updates all parameterized edges of the graph by given scalars.
//' 
//' @description
//' Given a vector of scalars, computes a new weight of
//' the parameterized edges in the graph by a simple inner
//' product of the edge state vector and the scalar vector.
//' 
//' @details
//' A faster and simpler version to compute new moments, when
//' the user wants to try multiple different weights.
//'
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param scalars A numeric vector of multiplies for the edge states.
//' 
//' @seealso [ptdalgorithms::expected_waiting_time()]
//' @seealso [ptdalgorithms::add_edge()]
//' 
//' @examples
//' graph <- create_graph(4)
//' v1 <- find_or_create_vertex(graph, c(1,2,1,0))
//' v2 <- find_or_create_vertex(graph, c(2,0,1,0))
//' add_edge(starting_vertex(graph), v1, 5)
//' add_edge(v1, v2, 0, c(5,2))
//' edges(starting_vertex(graph))[[1]]$weight # => 5
//' edges(v1)[[1]]$weight # => 0
//' graph_update_weights_parameterized(graph, c(9,7))
//' edges(starting_vertex(graph))[[1]]$weight # => 5
//' edges(v1)[[1]]$weight # => 59
// [[Rcpp::export]]
void graph_update_weights_parameterized(SEXP phase_type_graph, NumericVector scalars) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  graph->update_weights_parameterized(as<std::vector<double> >(scalars));
}

//' Obtain a list of all vertices in the graph
//' 
//' @description
//' Returns all vertices that have been added to the
//' graph from either calling `find_or_create_vertex` or
//' `create_vertex`. The first vertex in the list is
//' *always* the starting vertex [ptdalgorithms::starting_vertex()].
//' Importantly, for speed, use [ptdalgorithms::vertices_length()] to get the number
//' of added vertices, and use [ptdalgorithms::vertex_at()] to
//' get a vertex at a particular index.
//' 
//' @details
//' The list of vertices contains any added vertex, even
//' if it does not have any in-going / out-going edges.
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' 
//' @seealso [ptdalgorithms::starting_vertex()]
//' @seealso [ptdalgorithms::vertices_length()]
//' @seealso [ptdalgorithms::vertex_at()]
//' 
//' @examples
//' graph <- create_graph(4)
//' vertex_a <- find_or_create_vertex(graph, c(1,2,1,0))
//' vertex_b <- find_or_create_vertex(graph, c(2,0,1,0))
//' vertices(graph)[[1]] == starting_vertex(graph)
//' vertices(graph)[[2]] == vertex_at(graph, 2)
//' vertices_length(graph) == 3
// [[Rcpp::export]]
List vertices(SEXP phase_type_graph) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  std::vector<Vertex*> vertices = graph->vertices_p();
  List res(vertices.size());
  
  for (size_t i = 0; i < vertices.size(); i++) {
    res[i] = vertex_as_list(vertices[i]);
  }
  
  return res;
}

//' Returns a matrix where each row is the state of the vertex at that index
//' 
//' @return A matrix of size [ptdalgorithms::vertices_length()] where the rows match the state of the vertex at that index
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' 
//' @examples
//' graph <- ptdalgorithms::create_graph(4)
//' ptdalgorithms::create_vertex(graph, c(1,2,3,4))
//' ptdalgorithms::create_vertex(graph, c(4,3,3,3))
//' ptdalgorithms::states(graph) # => 
//' # 0 0 0 0
//' # 1 2 3 4
//' # 4 3 3 3
// [[Rcpp::export]]
IntegerMatrix states(SEXP phase_type_graph) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  std::vector<ptdalgorithms::Vertex *> ver = graph->vertices_p();
  
  IntegerMatrix states(ver.size(), graph->state_length());
  Graph g = *graph;
  
  for (size_t i = 0; i < ver.size(); i++) {
    for (size_t j = 0; j < graph->state_length(); j++) {
      states(i,j) = ver[i]->state()[j];
    }
  }
  
  return states;
}

//' Returns the number of vertices in the graph
//' 
//' @description
//' This method is much faster than calling `length(ptdalgorithms::vertices())`
//' 
//' @details
//' There will *always* be a starting vertex, so the returned
//' number is at least 1.
//' 
//' @return An integer of the number of vertices
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
// [[Rcpp::export]]
int vertices_length(SEXP phase_type_graph) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  return (int)graph->c_graph()->vertices_length;
}

//' Returns a vertex at a particular index.
//' 
//' @description
//' This method is much faster than calling `ptdalgorithms::vertices()[i]`
//' 
//' @return The vertex at index `index` in the graph
//'
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param index The index of the vertex to find
// [[Rcpp::export]]
List vertex_at(SEXP phase_type_graph, int index) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  if (index <= 0) {
    char message[1024];
    
    snprintf(
      message,
      1024, 
      "Failed: Index must be 1 or above, was %i",
      index
    );
    
    throw std::runtime_error(
        message
    );
  }
  
  if (index - 1 >= (int)graph->c_graph()->vertices_length) {
    char message[1024];
    
    snprintf(
      message,
      1024, 
      "Failed: Wanted to return vertex at %i, but graph has only %i vertices",
      index,
      (int)graph->c_graph()->vertices_length
    );
    
    throw std::runtime_error(
        message
    );
  }
  
  ptd_vertex *vertex = graph->c_graph()->vertices[index - 1];
  
  return vertex_as_list(new Vertex(*graph, vertex));
}

//' Returns the special starting vertex of the graph
//' 
//' @description
//' The starting vertex is always added to the graph
//' after calling [ptdalgorithms::create_graph()], and
//' always has the first index in [ptdalgorithms::vertex_at()]
//' 
//' @return The starting vertex
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' 
// [[Rcpp::export]]
List starting_vertex(SEXP phase_type_graph) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  Vertex *vertex = graph->starting_vertex_p();
  
  return vertex_as_list(vertex);
}

//' Returns the out-going edges of a vertex
//' 
//' @description
//' Returns a list of edges added by [ptdalgorithms::add_edge()]
//' 
//' @return A list of out-going edges
//' 
//' @param phase_type_vertex The vertex to find the edges for
//' 
// [[Rcpp::export]]
List edges(SEXP phase_type_vertex) {
  phase_type_vertex = list_as_vertex(phase_type_vertex);
  Rcpp::XPtr<Vertex> vertex(phase_type_vertex);
  std::vector<Edge> edges = vertex->edges();
  List r_edges(edges.size());
  
  for (size_t i = 0; i < edges.size(); i++) {
    Vertex child = edges[i].to();
    
    r_edges[i] = List::create(
      Named("weight") = edges[i].weight(),
      _["child"] = vertex_as_list(&child)
    );
  }
  
  return r_edges;
}

//' Create a vertex matching `state`
//' 
//' @description
//' Creates the vertex and adds it to
//' the graph object. Warning: the function [ptdalgorithms::find_or_create_vertex()]
//' should be preferred. This function will *not* update the lookup tree,
//' so [ptdalgorithms::find_vertex()] will *not* return it.
//' 
//' @seealso [ptdalgorithms::find_or_create_vertex()]
//' 
//' @return The newly inserted vertex in the graph
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param state An integer vector of the newly added vertex. Has length as given by `state_length` in  [ptdalgorithms::create_graph()]
// [[Rcpp::export]]
SEXP create_vertex(SEXP phase_type_graph, IntegerVector state) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  std::vector<int> v_state = as<std::vector<int> >(state);
  
  Vertex *vertex = graph->create_vertex_p(v_state);
  
  return vertex_as_list(vertex);
}


//' Finds a vertex matching `state`
//' 
//' 
//' @seealso [ptdalgorithms::find_or_create_vertex()]
//' 
//' @return The found vertex in the graph or NA
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param state An integer vector of the newly added vertex. Has length as given by `state_length` in  [ptdalgorithms::create_graph()]
// [[Rcpp::export]]
SEXP find_vertex(SEXP phase_type_graph, IntegerVector state) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  if (!graph->vertex_exists(as<std::vector<int> >(state))) {
    return List::get_na();
  }
  
  Vertex *found = graph->find_vertex_p(as<std::vector<int> >(state));
  
  return vertex_as_list(found);
}


//' Converts the graph-based phase-type distribution into a traditional sub-intensity matrix and initial probability vector
//' 
//' @details
//' Used to convert to the traditional matrix-based formulation.
//' Has three entries: `$SIM` the sub-intensity matrix, `$IPV` the initial
//' probability vector, `$states` the state of each vertex. Does
//' *not* have the same order as [ptdalgorithms::vertices()]
//' 
//' @seealso [ptdalgorithms::graph_as_dph_matrix()]
//'
//' @return A list of the sub-intensity matrix, states, and initial probability vector
//'
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' 
//' @examples
//' graph <- ptdalgorithms::create_graph(4)
//' v2 <- ptdalgorithms::create_vertex(graph, c(4,0,3,3))
//' v1 <- ptdalgorithms::create_vertex(graph, c(1,2,3,4))
//' a <- ptdalgorithms::create_vertex(graph, c(0,0,0,0))
//' # Note that if a defect is desired, the edge has to be explicitly
//' # added to an absorbing vertex
//' ptdalgorithms::add_edge(ptdalgorithms::starting_vertex(graph), v1, 1)
//' ptdalgorithms::add_edge(v1, v2, 4)
//' ptdalgorithms::add_edge(v2, a, 10)
//' ptdalgorithms::graph_as_matrix(graph) #>
//' # $`states`
//' #         [,1] [,2] [,3] [,4]
//' #   [1,]    1    2    3    4
//' #   [2,]    4    0    3    3
//' # $SIM
//' #         [,1]  [,2]
//' #   [1,]   -4     4
//' #   [2,]    0   -10
//' # $IPV
//' #   [1] 1 0
// [[Rcpp::export]]
List graph_as_matrix(SEXP phase_type_graph) {
  return(_graph_as_matrix(phase_type_graph));
}

//' Converts the the matrix-based representation into a phase-type graph
//' 
//' @details
//' Sometimes the user might want to use the fast graph algorithms,
//' but have some state-space given as a matrix. Therefore we can construct
//' a graph from a matrix. If desired, a discrete phase-type distribution
//' should just have no self-loop given. Note that the function
//' `graph_as_matrix` may reorder the vertices to make the graph represented
//' as strongly connected components in an acyclic manner.
//' 
//' @seealso [ptdalgorithms::matrix_as_graph()]
//'
//' @return A graph object
//'
//' @param IPV The initial probability vector (alpha)
//' @param SIM The sub-intensity matrix (S)
//' @param rewards Optional. The state/rewards of each of the vertices.
//' 
//' @examples
//' g <- matrix_as_graph(
//'     c(0.5,0.3, 0),
//'     matrix(c(-3, 0, 0, 2, -4, 1, 0, 1,-3), ncol=3),
//'     matrix(c(1,4,5,9,2,7), ncol=2)
//' )
//' 
//' graph_as_matrix(g)
// [[Rcpp::export]]
SEXP matrix_as_graph(NumericVector IPV, NumericMatrix SIM, Nullable<NumericMatrix> rewards = R_NilValue) {
  NumericMatrix rw;
  bool has_rewards;
  
  if (IPV.length() <= 0 || IPV.length() != SIM.ncol() || SIM.ncol() != SIM.nrow()) {
    char message[1024];
    
    snprintf(
      message,
      1024, 
      "Failed: IPV must have length > 0, was %i, and SIM must have same dimensions and be square, was %i, %i",
      (int)IPV.length(), (int)SIM.nrow(), (int)SIM.ncol()
    );
    
    throw std::runtime_error(
        message
    );
  }
  
  if (rewards.isNotNull()) {
    rw = NumericMatrix(rewards);
    has_rewards = true;
    
    if (rw.nrow() != SIM.nrow()) {
      char message[1024];
      
      snprintf(
        message,
        1024, 
        "Failed: Rewards must have %i rows, had %i",
        (int)SIM.nrow(), (int)rw.nrow()
      );
      
      throw std::runtime_error(
          message
      );
    }
  } else {
    has_rewards = false;
  }
  
  size_t state_space_size = has_rewards ? rw.ncol() : 1;
  
  Graph *cppGraph = new Graph(state_space_size);
  struct ptd_graph *graph = cppGraph->c_graph();
  
  for (int i = 0; i < SIM.nrow(); i++) {
    if (has_rewards) {
      int *state = (int*) calloc((size_t)rw.ncol(), sizeof(*state));
      
      for (int j = 0; j < rw.ncol(); j++) {
        state[j] = rw.at(i, j);
      }
      
      ptd_vertex_create_state(graph, state);
    } else {
      ptd_vertex_create(graph);
    }
  }
  
  struct ptd_vertex *absorbing_vertex = ptd_vertex_create(graph);
  double sum_outgoing = 0;
  
  for (int i = 0; i < SIM.nrow(); i++) {
    if (IPV[i] != 0) {
      ptd_graph_add_edge(graph->starting_vertex, graph->vertices[i+1], IPV[i]);
      sum_outgoing += IPV[i];
    }
  }
  
  if (sum_outgoing < 0.99999) {
    ptd_graph_add_edge(graph->starting_vertex, absorbing_vertex, 1 - sum_outgoing);
  }
  
  for (int i = 0; i < SIM.nrow(); i++) {
    double s = 0;
    
    for (int j = 0; j < SIM.nrow(); j++) {
      if (i == j) {
        continue;
      } 
      
      double weight = SIM.at(i, j);
      
      if (weight != 0) {
        ptd_graph_add_edge(graph->vertices[i+1], graph->vertices[j+1], weight);
        s += weight;
      }
    }
    
    double w = -(SIM.at(i, i) + s);
    
    if (w >= 0.000001) {
      ptd_graph_add_edge(graph->vertices[i+1], absorbing_vertex, w);
    }
  }
  
  return Rcpp::XPtr<Graph>(
    cppGraph
  ); 
}

//' Clones the graph, returning a new identical graph
//' 
//' @details
//' Many of the functions in this library mutate the graph. Use this to clone.
//' 
//' @return A clone of the input graph
//'
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' 
// [[Rcpp::export]]
SEXP clone_graph(SEXP phase_type_graph) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  return Rcpp::XPtr<Graph>(
    graph->clone_p()
  );
}

//' Converts the graph-based *discrete* phase-type distribution into a traditional sub-transition matrix and initial probability vector
//' 
//' @details
//' Used to convert to the traditional matrix-based formulation.
//' Has three entries: `$STM` the sub-transition matrix, `$IPV` the initial
//' probability vector, `$states` the state of each vertex. Does
//' *not* have the same order as [ptdalgorithms::vertices()]
//' It is expected that all out-going edges have weights summing to 1 or
//' less, and the remaining probability is considered as a self-transition.
//'
//' @seealso [ptdalgorithms::graph_as_matrix()]
//' @return A list of the sub-transition matrix, states, and initial probability vector
//'
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' 
//' @examples
//' graph <- ptdalgorithms::create_graph(4)
//' v2 <- ptdalgorithms::create_vertex(graph, c(4,0,3,3))
//' v1 <- ptdalgorithms::create_vertex(graph, c(1,2,3,4))
//' a <- ptdalgorithms::create_vertex(graph, c(0,0,0,0))
//' ptdalgorithms::add_edge(ptdalgorithms::starting_vertex(graph), v1, 0.5)
//' ptdalgorithms::add_edge(v1, v2, 0.8)
//' ptdalgorithms::add_edge(v2, a, 0.5)
//' # Note that if a defect is desired, the edge has to be explicitly
//' # added to an absorbing vertex
//' ptdalgorithms::add_edge(ptdalgorithms::starting_vertex(graph), a, 0.5)
//' ptdalgorithms::graph_as_dph_matrix(graph) #>
//' # $`states`
//' #         [,1] [,2] [,3] [,4]
//' #   [1,]    1    2    3    4
//' #   [2,]    4    0    3    3
//' # $STM
//' #         [,1]  [,2]
//' #   [1,]   0.2   0.8
//' #   [2,]   0.0   0.5
//' # $IPV
//' #   [1] 0.5 0
// [[Rcpp::export]]
List graph_as_dph_matrix(SEXP phase_type_graph) {
  List ph_res = _graph_as_matrix(phase_type_graph);
  NumericMatrix STM = ph_res["SIM"];
  
  
  for (int i = 0; i < STM.nrow(); i++) {
    double self_loop = 1-(-STM.at(i, i));
    STM(i, i) = self_loop;
  }
  
  return List::create(Named("states") = ph_res["states"] , _["STM"] = STM, _["IPV"] = ph_res["IPV"]);
}

//' Performs a reward transformation, returning a phase-type distribution to model the total accumulated reward until abosorption
//' 
//' @description
//' Returns a new `phase_type_graph` to instead model accumulated rewards until absorption,
//' as described in the paper. 
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param rewards to give to each vertex, zero or positive real number. Must have length equal to [ptdalgorithms::vertices_length()]
//' 
//' @return A new phase-type graph with the reward transformation applied
//' 
//' @examples
//' graph <- ptdalgorithms::create_graph(4)
//' v1 <- ptdalgorithms::create_vertex(graph, c(1,2,3,4))
//' v2 <- ptdalgorithms::create_vertex(graph, c(4,0,3,3))
//' a <- ptdalgorithms::create_vertex(graph, c(0,0,0,0))
//' ptdalgorithms::add_edge(ptdalgorithms::starting_vertex(graph), v1, 1)
//' ptdalgorithms::add_edge(v1, v2, 4)
//' ptdalgorithms::add_edge(v2, a, 10)
//' graph <- ptdalgorithms::reward_transform(graph, ptdalgorithms::states(graph)[,2])
//' # Graph now only has starting state and vertex `v1 and the absorbing vertex `a`.
//' # With edge weight 4/2 = 2
//' 
// [[Rcpp::export]]
SEXP reward_transform(SEXP phase_type_graph, NumericVector rewards) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  if ((int)rewards.size() != (int)graph->c_graph()->vertices_length) {
    char message[1024];
    
    snprintf(
      message,
      1024, 
      "Failed: Rewards must match the number of vertices. Expected %i, got %i",
      (int)graph->c_graph()->vertices_length,
      (int)rewards.size()
    );
    
    throw std::runtime_error(
        message
    );
  }
  
  std::vector<double> r = as<std::vector<double> >(rewards);
  
  return Rcpp::XPtr<Graph>(
    graph->reward_transform_p(r)
  );
}

//' Performs a discrete reward transformation, returning a discrete phase-type distribution to model the total accumulated reward until abosorption
//' 
//' @description
//' Changes `phase_type_graph` to instead model discrete rewarded jumps until absorption,
//' as described in the paper.
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param rewards to give to each vertex, zero or positive integer. Must have length equal to [ptdalgorithms::vertices_length()]
//' 
//' @return A new discrete phase-type graph with the reward transformation applied
//' 
//' @examples
//' graph <- ptdalgorithms::create_graph(4)
//' v2 <- ptdalgorithms::create_vertex(graph, c(4,0,3,3))
//' v1 <- ptdalgorithms::create_vertex(graph, c(1,2,3,4))
//' a <- ptdalgorithms::create_vertex(graph, c(0,0,0,0))
//' ptdalgorithms::add_edge(ptdalgorithms::starting_vertex(graph), v1, 0.5)
//' ptdalgorithms::add_edge(v1, v2, 0.8)
//' ptdalgorithms::add_edge(v2, a, 0.5)
//' # Note that if a defect is desired, the edge has to be explicitly
//' # added to an absorbing vertex
//' ptdalgorithms::add_edge(ptdalgorithms::starting_vertex(graph), a, 0.5)
//' ptdalgorithms::graph_as_dph_matrix(graph) #>
//' # $`states`
//' #         [,1] [,2] [,3] [,4]
//' #   [1,]    1    2    3    4
//' #   [2,]    4    0    3    3
//' # $STM
//' #         [,1]  [,2]
//' #   [1,]   0.2   0.8
//' #   [2,]   0.0   0.5
//' # $IPV
//' #   [1] 0.5 0
//' graph <- ptdalgorithms::dph_reward_transform(graph, ptdalgorithms::states(graph)[,2])
//' ptdalgorithms::graph_as_dph_matrix(graph) #>
//' # $`states`
//' # TODO: how will these look?
//' #         [,1] [,2] [,3] [,4]
//' #   [1,]    1    2    3    4
//' #   [2,]    4    0    3    3
//' # $STM
//' #         [,1]  [,2]
//' #   [1,]   0    0.2
//' #   [2,]   1    0.0
//' # $IPV
//' #   [1] 0 0.5
//' 
// [[Rcpp::export]]
SEXP dph_reward_transform(SEXP phase_type_graph, IntegerVector rewards) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  if ((int)rewards.size() != (int)graph->c_graph()->vertices_length) {
    char message[1024];
    
    snprintf(
      message,
      1024, 
      "Failed: Rewards must match the number of vertices. Expected %i, got %i",
      (int)graph->c_graph()->vertices_length,
      (int)rewards.size()
    );
    
    throw std::runtime_error(
        message
    );
  }
  
  std::vector<int> r = as<std::vector<int> >(rewards);
  
  return Rcpp::XPtr<Graph>(
    graph->dph_reward_transform_p(r)
  );
}

//' Normalizes the phase-type distribution, such that all out-going weights sum to 1 for each vertex
//' 
//' @description
//' Changes `phase_type_graph` such that all out-going weight sum to 1 for each vertex
//' as described in the paper, returns the associated inverse rates (new rewards). Mutates the phase-type graph permanently.
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' 
// [[Rcpp::export]]
NumericVector normalize_graph(SEXP phase_type_graph) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  return wrap(graph->normalize());
}

//' Normalizes the discrete phase-type distribution, by making self-loops go to an auxiliary vertex with reward 0 (such that all out-going weights sum to 1 for each vertex)
//' 
//' @description
//' Changes `phase_type_graph` such that all out-going weight sum to 1 for each vertex
//' and the self-loops go to a newly added vertex, with reward 0, that goes directly back,
//' as described in the paper, returns the associated new rewards. Mutates the phase-type graph permanently.
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' 
// [[Rcpp::export]]
NumericVector normalize_dph_graph(SEXP phase_type_graph) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  return wrap(graph->dph_normalize());
}

//' Computes the expected waiting time (or accumulated rewards) until absorption
//' 
//' @description
//' This function can be used to compute the moments of a 
//' phase-type distribution very fast and without much
//' memory usage compared with the traditional matrix-based
//' equations.
//' 
//' @return A numeric vector where entry `i` is the expected waiting time starting at vertex `i`
//' 
//' @seealso [ptdalgorithms::moments()]
//' @seealso [ptdalgorithms::expectation()]
//' @seealso [ptdalgorithms::variance()]
//' @seealso [ptdalgorithms::covariance()]
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param rewards Optional rewards, which should be applied to the phase-type distribution. Must have length equal to [ptdalgorithms::vertices_length()]
//' 
//' @examples
//' graph <- ptdalgorithms::create_graph(4)
//' v1 <- ptdalgorithms::create_vertex(graph, c(1,2,3,4))
//' v2 <- ptdalgorithms::create_vertex(graph, c(4,0,3,3))
//' a <- ptdalgorithms::create_vertex(graph, c(0,0,0,0))
//' ptdalgorithms::add_edge(ptdalgorithms::starting_vertex(graph), v1, 1)
//' ptdalgorithms::add_edge(v1, v2, 4)
//' ptdalgorithms::add_edge(v2, a, 10)
//' ptdalgorithms::expected_waiting_time(graph) # =>
//' # c(0.35, 0.35, 0.10, 0)
//' # Rewards on absorbing and starting vertex has no effect
//' ptdalgorithms::expected_waiting_time(graph, c(0, 2, 0, 0)) # =>
//' # c(0.5, 0.5, 0, 0)
//' ptdalgorithms::expected_waiting_time(graph, c(9999, 2, 0, 9999)) # =>
//' # c(0.5, 0.5, 0, 0)
// [[Rcpp::export]]
NumericVector expected_waiting_time(SEXP phase_type_graph, Nullable<NumericVector> rewards = R_NilValue) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  if (rewards.isNotNull() && (int)NumericVector(rewards).size() != (int)graph->c_graph()->vertices_length) {
    char message[1024];
    
    snprintf(
      message,
      1024, 
      "Failed: Rewards must match the number of vertices. Expected %i, got %i",
      (int)graph->c_graph()->vertices_length,
      (int)NumericVector(rewards).size()
    );
    
    throw std::runtime_error(
        message
    );
  }
  
  if (rewards.isNull()) {
    return wrap(graph->expected_waiting_time());
  } else {
    return wrap(graph->expected_waiting_time(as<std::vector<double> >(rewards)));
  }
}

//' Computes the expected jumps (or accumulated rewards) until absorption
//' 
//' @description
//' This function can be used to compute the moments of a 
//' disrete phase-type distribution very fast and without much
//' memory usage compared with the traditional matrix-based
//' equations.
//' 
//' The function takes in non-integers as rewards, but to be a *strictly* valid
//' rewarded discrete phase-type distribution these should be integers
//' 
//' @return A numeric vector where entry `i` is the expected rewarded jumps starting at vertex `i`
//' 
//' @seealso [ptdalgorithms::moments()]
//' @seealso [ptdalgorithms::expectation()]
//' @seealso [ptdalgorithms::variance()]
//' @seealso [ptdalgorithms::covariance()]
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param rewards Optional rewards, which should be applied to the discrete phase-type distribution. Must have length equal to [ptdalgorithms::vertices_length()]
//' 
// [[Rcpp::export]]
NumericVector dph_expected_waiting_time(SEXP phase_type_graph, Nullable<NumericVector> rewards = R_NilValue) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  if (rewards.isNotNull() && (int)NumericVector(rewards).size() != (int)graph->c_graph()->vertices_length) {
    char message[1024];
    
    snprintf(
      message,
      1024, 
      "Failed: Rewards must match the number of vertices. Expected %i, got %i",
      (int)graph->c_graph()->vertices_length,
      (int)NumericVector(rewards).size()
    );
    
    throw std::runtime_error(
        message
    );
  }
  
  if (rewards.isNull()) {
    return wrap(graph->dph_expected_waiting_time());
  } else {
    return wrap(graph->dph_expected_waiting_time(as<std::vector<double> >(rewards)));
  }
}

//' Computes the first `k` moments of the phase-type distribution
//' 
//' @description
//' This function invokes [ptdalgorithms::expected_waiting_times()] consequtively to find the first moments,
//' given by the `power` argument
//' 
//' @return A numeric vector of the first `k` moments. The first entry is the first moment (mean)
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param power An integer of the first `k` moments.
//' @param rewards Optional rewards, which should be applied to the phase-type distribution. Must have length equal to [ptdalgorithms::vertices_length()]
//' 
//' @seealso [ptdalgorithms::expected_waiting_time()]
//' @seealso [ptdalgorithms::expectation()]
//' @seealso [ptdalgorithms::variance()]
//' @seealso [ptdalgorithms::covariance()]
//' 
//' @examples
//' graph <- ptdalgorithms::create_graph(4)
//' v1 <- ptdalgorithms::create_vertex(graph, c(1,2,3,4))
//' v2 <- ptdalgorithms::create_vertex(graph, c(4,0,3,3))
//' a <- ptdalgorithms::create_vertex(graph, c(0,0,0,0))
//' ptdalgorithms::add_edge(ptdalgorithms::starting_vertex(graph), v1, 1)
//' ptdalgorithms::add_edge(v1, v2, 4)
//' ptdalgorithms::add_edge(v2, a, 10)
//' ptdalgorithms::moments(graph, 3) # =>
//'   (0.350000 0.097500 0.025375)
//' ptdalgorithms::moments(graph, 3, c(0,2,1,0)) # =>
//'   (0.600 0.160 0.041)
// [[Rcpp::export]]
NumericVector moments(SEXP phase_type_graph, int power, Nullable<NumericVector> rewards = R_NilValue) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  if (rewards.isNotNull() && (int)NumericVector(rewards).size() != (int)graph->c_graph()->vertices_length) {
    char message[1024];
    
    snprintf(
      message,
      1024, 
      "Failed: Rewards must match the number of vertices. Expected %i, got %i",
      (int)graph->c_graph()->vertices_length,
      (int)NumericVector(rewards).size()
    );
    
    throw std::runtime_error(
        message
    );
  }
  
  if (power <= 0) {
    char message[1024];
    
    snprintf(
      message,
      1024, 
      "Failed: power must be a strictly positive integer. Got %i",
      power
    );
    
    throw std::runtime_error(
        message
    );
  }
  
  NumericVector res(power);
  NumericVector rewards2 = expected_waiting_time(graph, rewards);
  NumericVector rewards3(rewards2.size());
  res[0] = rewards2[0];
  
  std::vector<double> rw;
  
  if (!rewards.isNull()) {
    rw = as<std::vector<double> >(rewards);
  }
  
  for (int i = 1; i < power; i++) {
    if (!rewards.isNull()) {
      for (int j = 0; j < rewards2.size(); j++) {
        rewards3[j] = rewards2[j] * rw[j];
      }
    } else {
      rewards3 = rewards2;
    }
    
    rewards2 = expected_waiting_time(graph, rewards3);
    res[i] = fac(i+1)*rewards2[0];
  }
  
  return res;
}

//' Computes the expectation (mean) of the phase-type distribution
//' 
//' @description
//' This function invokes [ptdalgorithms::expected_waiting_times()]
//' and takes the first entry (from starting vertex)
//' 
//' @return The expectation of the distribution
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param rewards Optional rewards, which should be applied to the phase-type distribution. Must have length equal to [ptdalgorithms::vertices_length()]
//' 
//' @seealso [ptdalgorithms::expected_waiting_time()]
//' @seealso [ptdalgorithms::moments()]
//' @seealso [ptdalgorithms::variance()]
//' @seealso [ptdalgorithms::covariance()]
//' 
//' @examples
//' graph <- ptdalgorithms::create_graph(4)
//' v1 <- ptdalgorithms::create_vertex(graph, c(1,2,3,4))
//' v2 <- ptdalgorithms::create_vertex(graph, c(4,0,3,3))
//' a <- ptdalgorithms::create_vertex(graph, c(0,0,0,0))
//' ptdalgorithms::add_edge(ptdalgorithms::starting_vertex(graph), v1, 1)
//' ptdalgorithms::add_edge(v1, v2, 4)
//' ptdalgorithms::add_edge(v2, a, 10)
//' ptdalgorithms::expectation(graph) # =>
//'   0.35
//' ptdalgorithms::expectation(graph, c(0,2,1,0)) # =>
//'   0.6
//' ph <- ptdalgorithms::graph_as_matrix(graph)
//' # This is a much faster version of
//' ph$IPV%*%solve(-ph$SIM) %*% rep(1, length(ph$IPV)) # =>
//'   0.35
//' ph$IPV%*%solve(-ph$SIM) %*% diag(c(2,1))%*% rep(1, length(ph$IPV)) # =>
//'   0.35
// [[Rcpp::export]]
double expectation(SEXP phase_type_graph, Nullable<NumericVector> rewards = R_NilValue) {
  return expected_waiting_time(phase_type_graph, rewards)[0];
}

//' Computes the variance of the phase-type distribution
//' 
//' @description
//' This function invokes [ptdalgorithms::expected_waiting_times()]
//' twice to find the first and second moment
//' 
//' @return The variance of the distribution
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param rewards Optional rewards, which should be applied to the phase-type distribution. Must have length equal to [ptdalgorithms::vertices_length()]
//' 
//' @seealso [ptdalgorithms::expected_waiting_time()]
//' @seealso [ptdalgorithms::expectation()]
//' @seealso [ptdalgorithms::moments()]
//' @seealso [ptdalgorithms::covariance()]
//' 
//' @examples
//' graph <- ptdalgorithms::create_graph(4)
//' v1 <- ptdalgorithms::create_vertex(graph, c(1,2,3,4))
//' v2 <- ptdalgorithms::create_vertex(graph, c(4,0,3,3))
//' a <- ptdalgorithms::create_vertex(graph, c(0,0,0,0))
//' ptdalgorithms::add_edge(ptdalgorithms::starting_vertex(graph), v1, 1)
//' ptdalgorithms::add_edge(v1, v2, 4)
//' ptdalgorithms::add_edge(v2, a, 10)
//' ptdalgorithms::variance(graph) # =>
//'   0.0725
//' ptdalgorithms::variance(graph, c(0,2,1,0)) # =>
//'   0.26
//' ph <- ptdalgorithms::graph_as_matrix(graph)
//' # This is a much faster version of
//' 2*ph$IPV%*%solve(-ph$SIM)%*%solve(-ph$SIM) %*% rep(1, length(ph$IPV)) - ph$IPV%*%solve(-ph$SIM) %*% rep(1, length(ph$IPV)) %*% ph$IPV%*%solve(-ph$SIM) %*% rep(1, length(ph$IPV)) # =>
//'   0.0725
//' 2*ph$IPV%*%solve(-ph$SIM)%*%diag(c(2,1))%*%solve(-ph$SIM)%*%diag(c(2,1)) %*% rep(1, length(ph$IPV)) - ph$IPV%*%solve(-ph$SIM)%*%diag(c(2,1)) %*% rep(1, length(ph$IPV)) %*% ph$IPV%*%solve(-ph$SIM)%*%diag(c(2,1)) %*% rep(1, length(ph$IPV)) # =>
//'   0.26
// [[Rcpp::export]]
double variance(SEXP phase_type_graph, Nullable<NumericVector> rewards = R_NilValue) {
  NumericVector exp = expected_waiting_time(phase_type_graph, rewards);
  NumericVector second;
  
  if (rewards.isNull()) {
    second = expected_waiting_time(phase_type_graph, exp);
  } else {
    NumericVector new_rewards(exp.size());
    std::vector<double> rw = as<std::vector<double> >(rewards);
    
    for (int i = 0; i < exp.size(); i++) {
      new_rewards[i] = exp[i] * rw[i];
    }
    
    second = expected_waiting_time(phase_type_graph, new_rewards);
  }
  
  return (2*second[0]-exp[0]*exp[0]);
}

//' Computes the covariance of the phase-type distribution
//' 
//' @description
//' This function invokes [ptdalgorithms::expected_waiting_times()]
//' twice to find the first and second moment for each of the two rewards
//' 
//' @return The covariance of the distribution given the two rewards
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param rewards1 First vector of rewards, which should be applied to the phase-type distribution. Must have length equal to [ptdalgorithms::vertices_length()]
//' @param rewards2 Second vector of rewards, which should be applied to the phase-type distribution. Must have length equal to [ptdalgorithms::vertices_length()]
//'
//' @seealso [ptdalgorithms::expected_waiting_time()]
//' @seealso [ptdalgorithms::expectation()]
//' @seealso [ptdalgorithms::variance()]
//' @seealso [ptdalgorithms::moments()]
//' 
//' @examples
//' graph <- ptdalgorithms::create_graph(4)
//' v1 <- ptdalgorithms::create_vertex(graph, c(1,2,3,4))
//' v2 <- ptdalgorithms::create_vertex(graph, c(4,0,3,3))
//' a <- ptdalgorithms::create_vertex(graph, c(0,0,0,0))
//' ptdalgorithms::add_edge(ptdalgorithms::starting_vertex(graph), v1, 1)
//' ptdalgorithms::add_edge(v1, v2, 4)
//' ptdalgorithms::add_edge(v2, a, 10)
//' covariance(graph, c(0,2,1,0), c(0,2,1,0)) == variance(graph, c(0,2,1,0))
//' covariance(graph, c(0,2,1,0), c(0,5,2,0)) # =>
//'   0.645
//' ph <- ptdalgorithms::graph_as_matrix(graph)
//' # This is a much faster version of
//' ph$IPV%*%solve(-ph$SIM)%*%diag(c(2,1))%*%solve(-ph$SIM)%*%diag(c(5,2)) %*% rep(1, length(ph$IPV)) + ph$IPV%*%solve(-ph$SIM)%*%diag(c(5,2))%*%solve(-ph$SIM)%*%diag(c(2,1)) %*% rep(1, length(ph$IPV)) - ph$IPV%*%solve(-ph$SIM)%*%diag(c(2,1)) %*% rep(1, length(ph$IPV)) %*% ph$IPV%*%solve(-ph$SIM)%*%diag(c(5,2)) %*% rep(1, length(ph$IPV)) # =>
//'   0.645
// [[Rcpp::export]]
double covariance(SEXP phase_type_graph, NumericVector rewards1, NumericVector rewards2) {
  NumericVector exp1 = expected_waiting_time(phase_type_graph, rewards1);
  NumericVector exp2 = expected_waiting_time(phase_type_graph, rewards2);
  
  NumericVector new_rewards(exp1.size());
  
  for (int i = 0; i < exp1.size(); i++) {
    new_rewards[i] = exp1[i] * rewards2[i];
  }
  
  NumericVector second1 = expected_waiting_time(phase_type_graph, new_rewards);
  
  
  for (int i = 0; i < exp1.size(); i++) {
    new_rewards[i] = exp2[i] * rewards1[i];
  }
  
  NumericVector second2 = expected_waiting_time(phase_type_graph, new_rewards);
  
  return (second1[0]+second2[0]-exp1[0]*exp2[0]);
}

//' Computes the expectation (mean) of the discrete phase-type distribution
//' 
//' @description
//' This function invokes [ptdalgorithms::dph_expected_waiting_times()]
//' and takes the first entry (from starting vertex)
//' 
//' @return The expectation of the distribution
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param rewards Optional rewards, which should be applied to the discrete phase-type distribution. Must have length equal to [ptdalgorithms::vertices_length()]
//' 
//' @seealso [ptdalgorithms::dph_expected_waiting_time()]
//' @seealso [ptdalgorithms::expected_waiting_time()]
//' @seealso [ptdalgorithms::dph_variance()]
//' @seealso [ptdalgorithms::dph_covariance()]
//' 
// [[Rcpp::export]]
double dph_expectation(SEXP phase_type_graph, Nullable<NumericVector> rewards = R_NilValue) {
  return dph_expected_waiting_time(phase_type_graph, rewards)[0];
}

//' Computes the variance of the discrete phase-type distribution
//' 
//' @description
//' This function invokes [ptdalgorithms::dph_expected_waiting_times()]
//' twice to find the first and second moment
//' 
//' @return The variance of the distribution
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param rewards Optional rewards, which should be applied to the discrete phase-type distribution. Must have length equal to [ptdalgorithms::vertices_length()]
//' 
//' @seealso [ptdalgorithms::dph_expected_waiting_time()]
//' @seealso [ptdalgorithms::dph_covariance()]
//' @seealso [ptdalgorithms::dph_expectation()]
//' 
// [[Rcpp::export]]
double dph_variance(SEXP phase_type_graph, Nullable<NumericVector> rewards = R_NilValue) {
  NumericVector exp = dph_expected_waiting_time(phase_type_graph, rewards);
  NumericVector second;
  
  if (rewards.isNull()) {
    second = dph_expected_waiting_time(phase_type_graph, exp);
  } else {
    NumericVector new_rewards(exp.size());
    std::vector<double> rw = as<std::vector<double> >(rewards);
    
    for (int i = 0; i < exp.size(); i++) {
      new_rewards[i] = exp[i] * rw[i];
    }
    
    second = dph_expected_waiting_time(phase_type_graph, new_rewards);
  }
  
  NumericVector sq_rewards;
  
  if (rewards.isNull()) {
    sq_rewards = exp;
  } else {
    NumericVector new_rewards(exp.size());
    std::vector<double> rw = as<std::vector<double> >(rewards);
    
    for (int i = 0; i < (int)rw.size(); i++) {
      new_rewards[i] = rw[i] * rw[i];
    }
    
    sq_rewards = dph_expected_waiting_time(phase_type_graph, new_rewards);
  }
  
  return (2*second[0]-sq_rewards[0]-exp[0]*exp[0]);
}

//' Computes the covariance of the discrete phase-type distribution
//' 
//' @description
//' This function invokes [ptdalgorithms::dph_expected_waiting_times()]
//' twice to find the first and second moment for each of the two rewards
//' 
//' @return The covariance of the discrete distribution given the two rewards
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param rewards1 First vector of rewards, which should be applied to the discrete phase-type distribution. Must have length equal to [ptdalgorithms::vertices_length()]
//' @param rewards2 Second vector of rewards, which should be applied to the discrete  phase-type distribution. Must have length equal to [ptdalgorithms::vertices_length()]
//'
//' @seealso [ptdalgorithms::dph_expected_waiting_time()]
//' @seealso [ptdalgorithms::dph_expectation()]
//' @seealso [ptdalgorithms::dph_variance()]
//' 
// [[Rcpp::export]]
double dph_covariance(SEXP phase_type_graph, NumericVector rewards1, NumericVector rewards2) {
  NumericVector exp1 = dph_expected_waiting_time(phase_type_graph, rewards1);
  NumericVector exp2 = dph_expected_waiting_time(phase_type_graph, rewards2);
  
  NumericVector new_rewards(exp1.size());
  
  for (int i = 0; i < exp1.size(); i++) {
    new_rewards[i] = exp1[i] * rewards2[i];
  }
  
  NumericVector second1 = dph_expected_waiting_time(phase_type_graph, new_rewards);
  
  
  for (int i = 0; i < exp1.size(); i++) {
    new_rewards[i] = exp2[i] * rewards1[i];
  }
  
  NumericVector second2 = dph_expected_waiting_time(phase_type_graph, new_rewards);
  
  NumericVector sq_rewards(rewards1.size());
  
  for (int i = 0; i < rewards1.size(); i++) {
    sq_rewards[i] = rewards1[i] * rewards2[i];
  }
  
  NumericVector sq = dph_expected_waiting_time(phase_type_graph, sq_rewards);
  
  
  return (second1[0]+second2[0]-sq[0]-exp1[0]*exp2[0]);
}

// [[Rcpp::export]]
bool is_graph_acyclic(SEXP phase_type_graph) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  return graph->is_acyclic();
}

/* Utility */

SEXP get_first_list_entry(SEXP e, std::string message) {
  if (Rf_isList(e) || Rf_isNewList(e)) {
    List list = Rcpp::as<List>(e);
    
    if (list.size() != 1) {
      char message[1024];
      
      snprintf(
        message,
        1024, 
        "Failed: When finding %s of a list-type of vertices, list can only contain 1 vertex, it contained %i. Did you use [] as lookup instead of [[]]?",
        message,
        (int)list.size()
      );
      
      throw std::runtime_error(
          message
      );
    }
    
    e = list[0];
  }
  
  return e;
}

List vertex_as_list(Vertex *vertex) {
  std::vector<int> state = vertex->state();
  IntegerVector state_vec(state.begin(), state.end());
  
  return List::create(
    Named("state") = state_vec,
    _["rate"] = vertex->rate(),
    _["index"] = vertex->c_vertex()->index + 1,
    _["vertex"] = (size_t)vertex->c_vertex(),
    _["xptr_vertex"] = Rcpp::XPtr<Vertex>(vertex)
  );
}

List vertex_as_list(Graph &graph, struct ptd_vertex *c_vertex) {
  Vertex *vertex = new Vertex(graph, c_vertex);
  
  return vertex_as_list(vertex);
}

SEXP list_as_vertex(SEXP list) {
  if (!Rf_isList(list) && !Rf_isNewList(list)) {
    char message[1024];
    
    snprintf(
      message,
      1024, 
      "Failed: vertex must be given as the returned lists, the datatype was '%i' (R internal type description)",
      (int)TYPEOF(list)
    );
    
    throw std::runtime_error(
        message
    );
  }
  
  Rcpp::List child = Rcpp::as<Rcpp::List>(list);
  
  return child["xptr_vertex"];
}

List _graph_as_matrix(SEXP phase_type_graph) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  PhaseTypeDistribution dist = graph->phase_type_distribution();
  
  NumericMatrix SIM(dist.length, dist.length);
  NumericVector IPV(dist.length);
  
  for (size_t i = 0; i < dist.length; ++i) {
    IPV(i) = dist.initial_probability_vector[i];
    
    for (size_t j = 0; j < dist.length; ++j) {
      SIM(i, j) = dist.sub_intensity_matrix[i][j];
    }
  }
  
  NumericMatrix states(dist.length, graph->state_length());
  Graph g = *graph;
  
  for (size_t i = 0; i < dist.length; i++) {
    for (size_t j = 0; j < graph->state_length(); j++) {
      states(i,j) = dist.vertices[i].state()[j];
    }
  }
  
  return List::create(Named("states") = states , _["SIM"] = SIM, _["IPV"] = IPV);
}

//' Computes the defect, i.e. the probability of immediately transitioning to the absorbing state, of the phase-type distribution
//'
//' @description
//' Returns the defect of the distribution
//'
//' @return A numeric value
//'
//' @seealso [ptdalgorithms::dph()]
//' @seealso [ptdalgorithms::qph()]
//'
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//'
// [[Rcpp::export]]
double defect(SEXP phase_type_graph) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  return graph->defect();
}

//' Samples `n` outcomes of the phase-type distribution
//' 
//' 
//' @return A vector of numeric values
//' 
//' @seealso [ptdalgorithms::rmph()]
//' 
//' @param n The number of random samples to draw
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param rewards Optional rewards, which should be applied to the phase-type distribution. Must have length equal to [ptdalgorithms::vertices_length()]
//' 
// [[Rcpp::export]]
NumericVector rph(int n, SEXP phase_type_graph, Nullable<NumericVector> rewards = R_NilValue) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  if (rewards.isNotNull() && (int)NumericVector(rewards).size() != (int)graph->c_graph()->vertices_length) {
    char message[1024];
    
    snprintf(
      message,
      1024, 
      "Failed: Rewards must match the number of vertices. Expected %i, got %i",
      (int)graph->c_graph()->vertices_length,
      (int)NumericVector(rewards).size()
    );
    
    throw std::runtime_error(
        message
    );
  }
  
  set_c_seed();
  
  NumericVector res(n);
  
  for (int i = 0; i < n; i++) {
    if (rewards.isNull()) {
      res[i] = (double)(graph->random_sample());
    } else {
      res[i] = (double)(graph->random_sample(as<std::vector<double> >(rewards)));
    }
  }
  
  return res;
}

//' Probability density function of the phase-type distribution
//'
//' @description
//' Returns the density (probability density function) at a specific time. Notice that this
//' is an approximation
//'
//' @return A numeric vector of the density
//'
//' @seealso [ptdalgorithms::rph()]
//' @seealso [ptdalgorithms::qph()]
//'
//' @param x Vector of the number of stopping times
//' @param granularity Optional. Number of jumps per time unit (approximation)
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//'
// [[Rcpp::export]]
NumericVector dph(NumericVector x, SEXP phase_type_graph, int granularity=0) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  NumericVector res(x.length());
  
  for (int i = 0; i < x.length(); i++) {
    double t = x[i];
    
    res[i] = graph->pdf(t, granularity);
  }
  
  return res;
}

//' Cumulative distribution function of the phase-type distribution
//'
//' @description
//' Returns the density (cumulative probability density function) at a specific time. Notice that this
//' is an approximation
//'
//' @return A numeric vector of the distribution function
//'
//' @seealso [ptdalgorithms::rph()]
//' @seealso [ptdalgorithms::dph()]
//'
//' @param q Vector of the quantiles (stopping times)
//' @param granularity Optional. Number of jumps per time unit (approximation)
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//'
// [[Rcpp::export]]
NumericVector pph(NumericVector q, SEXP phase_type_graph, int granularity = 0) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  NumericVector res(q.length());
  
  for (int i = 0; i < q.length(); i++) {
    double t = q[i];
    
    res[i] = graph->cdf(t, granularity);
  }
  
  return res;
}

//' Computes the probability of the Markov Jump Process of the phase-type distribution standing at each vertex after time
//'
//' @return A numeric vector of the stop probabilities, for each vertex
//'
//' @seealso [ptdalgorithms::rph()]
//' @seealso [ptdalgorithms::dph()]
//' @seealso [ptdalgorithms::pph()]
//'
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param time The time to wait
//' @param granularity Optional. Number of jumps per time unit (approximation)
//'
// [[Rcpp::export]]
NumericVector stop_probability(SEXP phase_type_graph, double time, int granularity = 0) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  return wrap(graph->stop_probability(time, granularity));
}

//' Computes the accumulated visiting time (integral of stop probability) of the Markov Jump Process of the phase-type distribution standing at each vertex
//'
//' @return A numeric vector of the accumulated visiting time, for each vertex
//'
//' @seealso [ptdalgorithms::rph()]
//' @seealso [ptdalgorithms::dph()]
//' @seealso [ptdalgorithms::pph()]
//'
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param time The time to wait
//' @param granularity Optional. Number of jumps per time unit (approximation)
//'
// [[Rcpp::export]]
NumericVector accumulated_visiting_time(SEXP phase_type_graph, double time, int granularity = 0) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  return wrap(graph->accumulated_visiting_time(time, granularity));
}

//' Samples `n` outcomes of the multivariate phase-type distribution
//' 
//' @return A matrix of numeric values, with the parameter `rewards` number of rows
//' 
//' @seealso [ptdalgorithms::rph()]
//' 
//' @param n The number of random samples to draw
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param rewards rewards which should be applied to the phase-type distribution. Must have rows equal to [ptdalgorithms::vertices_length()] and columns the dimensions of the multivariate phase-type distribution
//' 
// [[Rcpp::export]]
NumericMatrix rmph(int n, SEXP phase_type_graph, NumericMatrix rewards) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  if ((int)rewards.rows() != (int)graph->c_graph()->vertices_length) {
    char message[1024];
    
    snprintf(
      message,
      1024, 
      "Failed: Rewards rows must match the number of vertices. Expected %i, got %i",
      (int)graph->c_graph()->vertices_length,
      (int)rewards.rows()
    );
    
    throw std::runtime_error(
        message
    );
  }
  
  double *vrewards = (double*) calloc(rewards.rows() * rewards.cols(), sizeof(double));
  
  size_t index = 0;
  
  for (int i = 0; i < rewards.rows(); i++) {
    for (int j = 0; j < rewards.cols(); j++) {
      vrewards[index] = rewards.at(i, j);
      index++;
    }
  }
  
  set_c_seed();
  
  NumericMatrix mat_res(rewards.ncol(), n);
  
  for (int i = 0; i < n; i++) {
    long double *res = ptd_mph_random_sample(graph->c_graph(), vrewards, (size_t)rewards.ncol());
    
    for (int j = 0; j < rewards.cols(); j++) {
      mat_res(j, i) = res[j];
    }
    
    free(res);
  }
  
  free(vrewards);
  
  return mat_res;
}


//' Samples `n` stopping times of the phase-type distribution
//'
//' @return A vector of numeric values, the indices of the vertices at the stopping time
//'
//' @seealso [ptdalgorithms::rph()]
//'
//' @param n The number of random samples to draw
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param time The stopping time
//'
// [[Rcpp::export]]
NumericVector random_sample_stop_vertex(int n, SEXP phase_type_graph, double time) {
    Rcpp::XPtr<Graph> graph(phase_type_graph);

    set_c_seed();

    NumericVector res(n);

    for (int i = 0; i < n; i++) {
        res[i] = graph->random_sample_stop_vertex(time) + 1;
    }

    return res;
}

//' Samples `n` outcomes of the discrete phase-type distribution
//' 
//' @description
//' Allows set.seed as any R function, as it seeds with a random number
//' generated by R.
//' 
//' @return A numeric vector of the draws
//' 
//' @seealso [ptdalgorithms::rph()]
//' @seealso [ptdalgorithms::rmdph()]
//' 
//' @param n The number of random samples to draw
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param rewards Optional rewards, which should be applied to the discrete phase-type distribution. Must have length equal to [ptdalgorithms::vertices_length()]
//' 
// [[Rcpp::export]]
NumericVector rdph(int n, SEXP phase_type_graph, Nullable<NumericVector> rewards = R_NilValue) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  if (rewards.isNotNull() && (int)NumericVector(rewards).size() != (int)graph->c_graph()->vertices_length) {
    char message[1024];
    
    snprintf(
      message,
      1024, 
      "Failed: Rewards must match the number of vertices. Expected %i, got %i",
      (int)graph->c_graph()->vertices_length,
      (int)NumericVector(rewards).size()
    );
    
    throw std::runtime_error(
        message
    );
  }
  
  set_c_seed();
  
  NumericVector res(n);
  
  for (int i = 0; i < n; i++) {
    if (rewards.isNull()) {
      res[i] = (double)(graph->dph_random_sample_c(NULL));
    } else {
      std::vector<double> rw = as<std::vector<double> >(rewards);
      res[i] = (double)(graph->dph_random_sample_c(&rw[0]));
    }
  }
  
  return res;
}

//' Probability mass function of the discrete phase-type distribution
//' 
//' @description
//' Returns the density (probability mass function) at a specific number of jumps.
//' 
//' @return A numeric vector of the density
//' 
//' @seealso [ptdalgorithms::rdph()]
//' @seealso [ptdalgorithms::qdph()]
//' 
//' @param x Vector of the number of jumps (discrete time)
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' 
// [[Rcpp::export]]
NumericVector ddph(IntegerVector x, SEXP phase_type_graph) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  NumericVector res(x.length());
  
  for (int i = 0; i < x.length(); i++) {
    int t = x[i];
    
    res[i] = graph->dph_pmf(t);
  }
  
  return res;
}

//' Cumulative distribution function of the discrete phase-type distribution
//' 
//' @description
//' Returns the density (probability mass function) at a specific number of jumps.
//' 
//' @return A numeric vector of the distribution function
//' 
//' @seealso [ptdalgorithms::rdph()]
//' @seealso [ptdalgorithms::ddph()]
//' 
//' @param q Vector of the quantiles (jumps, discrete time)
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' 
// [[Rcpp::export]]
NumericVector pdph(IntegerVector q, SEXP phase_type_graph) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  NumericVector res(q.length());
  
  for (int i = 0; i < q.length(); i++) {
    int t = q[i];
    
    res[i] = graph->dph_cdf(t);
  }
  
  return res;
}

//' Computes the probability of the Markov Chain of the discrete phase-type distribution standing at each vertex after n jumps
//' 
//' @return A numeric vector of the stop probabilities, for each vertex
//' 
//' @seealso [ptdalgorithms::rdph()]
//' @seealso [ptdalgorithms::ddph()]
//' @seealso [ptdalgorithms::pdph()]
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param jumps Number of jumps, discrete time
//' 
// [[Rcpp::export]]
NumericVector dph_stop_probability(SEXP phase_type_graph, int jumps) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  return wrap(graph->dph_stop_probability(jumps));
}

//' Computes the number of visits of the Markov Chain of the discrete phase-type distribution at each vertex after n jumps
//' 
//' @return A numeric vector of the stop probabilities, for each vertex
//' 
//' @seealso [ptdalgorithms::rdph()]
//' @seealso [ptdalgorithms::ddph()]
//' @seealso [ptdalgorithms::pdph()]
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param jumps Number of jumps, discrete time
//' 
// [[Rcpp::export]]
NumericVector dph_accumulated_visits(SEXP phase_type_graph, int jumps) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  return wrap(graph->dph_accumulated_visits(jumps));
}

//' Samples one outcome of the multivariate discrete phase-type distribution
//' 
//' 
//' @return A matrix of numeric values, with `rewards` number of rows and `n` number of columns
//' 
//' @seealso [ptdalgorithms::rmph()]
//' @seealso [ptdalgorithms::rdph()]
//' 
//' @param n The number of random samples to draw
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param rewards rewards which should be applied to the discrete phase-type distribution. Must have rows equal to [ptdalgorithms::vertices_length()] and columns the dimensions of the multivariate discrete phase-type distribution
//' 
// [[Rcpp::export]]
NumericMatrix rmdph(int n, SEXP phase_type_graph, NumericMatrix rewards) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  if ((int)rewards.rows() != (int)graph->c_graph()->vertices_length) {
    char message[1024];
    
    snprintf(
      message,
      1024, 
      "Failed: Rewards rows must match the number of vertices. Expected %i, got %i",
      (int)graph->c_graph()->vertices_length,
      (int)rewards.rows()
    );
    
    throw std::runtime_error(
        message
    );
  }
  
  
  double *vrewards = (double*) calloc(rewards.rows() * rewards.cols(), sizeof(double));
  
  size_t index = 0;
  
  for (int i = 0; i < rewards.rows(); i++) {
    for (int j = 0; j < rewards.cols(); j++) {
      vrewards[index] = rewards.at(i, j);
      index++;
    }
  }
  
  set_c_seed();
  
  NumericMatrix mat_res(rewards.ncol(), n);
  
  for (int i = 0; i < n; i++) {
    long double *res = ptd_mdph_random_sample(graph->c_graph(), vrewards, (size_t)rewards.ncol());
    
    for (int j = 0; j < rewards.cols(); j++) {
      mat_res(j, i) = res[j];
    }
    
    free(res);
  }
  
  free(vrewards);
  
  return mat_res;
}


//' Samples `n` vertices at the stopping times of the discrete phase-type distribution
//'
//' @return A vector of numeric values, the indices of the vertices at the stopping time (jumps)
//'
//' @seealso [ptdalgorithms::rdph()]
//'
//' @param n The number of random samples to draw
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param jumps The number of jumps before stopping
//'
// [[Rcpp::export]]
NumericVector dph_random_sample_stop_vertex(int n, SEXP phase_type_graph, int jumps) {
    Rcpp::XPtr<Graph> graph(phase_type_graph);

    set_c_seed();

    NumericVector res(n);

    for (int i = 0; i < n; i++) {
        res[i] = graph->dph_random_sample_stop_vertex(jumps) + 1;
    }

    return res;
}

//' Builds a probability distribution context for the phase-type distribution.
//' 
//' @description
//' This allows the user to step through the distribution, computing e.g. the
//' time-inhomogeneous distribution function or the expectation of a multivariate phase-type distribution
//' 
//' 
//' @return A probability distribution context
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' @param granularity Optional. Number of jumps per time unit (approximation)
//' 
// [[Rcpp::export]]
SEXP distribution_context(SEXP phase_type_graph, int granularity = 0) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  return Rcpp::XPtr<ProbabilityDistributionContext>(
    new ProbabilityDistributionContext(*graph, granularity)
  );
}

//' Performs one jump in the probability distribution context for the discrete phase-type distribution.
//' 
//' @description
//' This allows the user to step through the distribution, computing e.g. the
//' time-inhomogeneous distribution function or the expectation of a multivariate phase-type distribution.
//' *mutates* the context
//' 
//' @seealso [ptdalgorithms::distribution_context()]
//' 
//' @param probability_distribution_context The context created by [ptdalgorithms::distribution_context()]
//' 
// [[Rcpp::export]]
void distribution_context_step(SEXP probability_distribution_context) {
  Rcpp::XPtr<ProbabilityDistributionContext> context(probability_distribution_context);
  
  context->step();
}

//' Returns the PDF/PMF/CDF/time for the current probability distribution context for the phase-type distribution.
//' 
//' @description
//' This allows the user to step through the distribution, computing e.g. the
//' time-inhomogeneous distribution function or the expectation of a multivariate phase-type distribution.
//' *mutates* the context
//' 
//' @seealso [ptdalgorithms::distribution_context()]
//' 
//' @param probability_distribution_context The context created by [ptdalgorithms::distribution_context()]
//' 
// [[Rcpp::export]]
List distribution_context_state(SEXP probability_distribution_context) {
  Rcpp::XPtr<ProbabilityDistributionContext> context(probability_distribution_context);
  
  return List::create(Named("pdf") = context->pdf() , _["cdf"] = context->cdf(), _["time"] = context->time());
}


//' Returns the stop probability for the current probability distribution context for the phase-type distribution.
//' 
//' @description
//' This allows the user to step through the distribution, computing e.g. the
//' time-inhomogeneous distribution function or the expectation of a multivariate phase-type distribution.
//' *mutates* the context
//' 
//' @seealso [ptdalgorithms::distribution_context()]
//' 
//' @param probability_distribution_context The context created by [ptdalgorithms::distribution_context()]
//' 
// [[Rcpp::export]]
NumericVector distribution_context_stop_probability(SEXP probability_distribution_context) {
  Rcpp::XPtr<ProbabilityDistributionContext> context(probability_distribution_context);
  
  return wrap(context->stop_probability());
}

//' Returns the accumulated visiting time (integral of stop probability) for the current probability distribution context for the phase-type distribution.
//' 
//' @description
//' This allows the user to step through the distribution, computing e.g. the
//' time-inhomogeneous distribution function or the expectation of a multivariate phase-type distribution.
//' *mutates* the context
//' 
//' @seealso [ptdalgorithms::distribution_context()]
//' 
//' @param probability_distribution_context The context created by [ptdalgorithms::distribution_context()]
//' 
// [[Rcpp::export]]
NumericVector distribution_context_accumulated_visiting_time(SEXP probability_distribution_context) {
  Rcpp::XPtr<ProbabilityDistributionContext> context(probability_distribution_context);
  
  return wrap(context->accumulated_visiting_time());
}

//' Builds a probability distribution context for the discrete phase-type distribution.
//' 
//' @description
//' This allows the user to step through the distribution, computing e.g. the
//' time-inhomogeneous distribution function or the expectation of a multivariate discrete discrete phase-type distribution
//' 
//' @return A probability distribution context
//' 
//' @param phase_type_graph A reference to the graph created by [ptdalgorithms::create_graph()]
//' 
// [[Rcpp::export]]
SEXP dph_distribution_context(SEXP phase_type_graph) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  return Rcpp::XPtr<DPHProbabilityDistributionContext>(
    new DPHProbabilityDistributionContext(*graph)
  );
}

//' Performs one jump in a probability distribution context for the discrete phase-type distribution.
//' 
//' @description
//' This allows the user to step through the distribution, computing e.g. the
//' time-inhomogeneous distribution function or the expectation of a multivariate discrete phase-type distribution.
//' *mutates* the context
//' 
//' @seealso [ptdalgorithms::dph_distribution_context()]
//' 
//' @param probability_distribution_context The context created by [ptdalgorithms::dph_distribution_context()]
//' 
// [[Rcpp::export]]
void dph_distribution_context_step(SEXP probability_distribution_context) {
  Rcpp::XPtr<DPHProbabilityDistributionContext> context(probability_distribution_context);
  
  context->step();
}

//' Returns the PMF/CDF/time for the current probability distribution context for the discrete phase-type distribution.
//' 
//' @description
//' This allows the user to step through the distribution, computing e.g. the
//' time-inhomogeneous distribution function or the expectation of a  discrete multivariate phase-type distribution.
//' *mutates* the context
//' 
//' @seealso [ptdalgorithms::dph_distribution_context()]
//' 
//' @param probability_distribution_context The context created by [ptdalgorithms::dph_distribution_context()]
//' 
// [[Rcpp::export]]
List dph_distribution_context_state(SEXP probability_distribution_context) {
  Rcpp::XPtr<DPHProbabilityDistributionContext> context(probability_distribution_context);
  
  return List::create(Named("pmf") = context->pmf() , _["cdf"] = context->cdf(), _["jumps"] = context->jumps());
}


//' Returns the stop probability for the current probability distribution context for the discrete phase-type distribution.
//' 
//' @description
//' This allows the user to step through the distribution, computing e.g. the
//' time-inhomogeneous distribution function or the expectation of a multivariate discrete phase-type distribution.
//' *mutates* the context
//' 
//' @seealso [ptdalgorithms::dph_distribution_context()]
//' 
//' @param probability_distribution_context The context created by [ptdalgorithms::dph_distribution_context()]
//' 
// [[Rcpp::export]]
NumericVector dph_distribution_context_stop_probability(SEXP probability_distribution_context) {
  Rcpp::XPtr<DPHProbabilityDistributionContext> context(probability_distribution_context);
  
  return wrap(context->stop_probability());
}

//' Returns the accumulated visits for the current probability distribution context for the discrete phase-type distribution.
//' 
//' @description
//' This allows the user to step through the distribution, computing e.g. the
//' time-inhomogeneous distribution function or the expectation of a multivariate discrete phase-type distribution.
//' *mutates* the context
//' 
//' @seealso [ptdalgorithms::dph_distribution_context()]
//' 
//' @param probability_distribution_context The context created by [ptdalgorithms::dph_distribution_context()]
//' 
// [[Rcpp::export]]
NumericVector dph_distribution_context_accumulated_visits(SEXP probability_distribution_context) {
  Rcpp::XPtr<DPHProbabilityDistributionContext> context(probability_distribution_context);
  
  return wrap(context->accumulated_visits());
}


#endif