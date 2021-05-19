#define PTD_RCPP 1
#include <Rcpp.h>

#include "c/ptdalgorithms.h"
#include "cpp/ptdalgorithmscpp.h"

using namespace Rcpp;
using namespace ptdalgorithms;

// TODO: Make all functions exists as Cpp method calls

/*** R
n <- 5
graph <- create_graph(n)
vertices_list <- graph$vertices_list()

while (vertices_list$has_next()) {
  vertex <- vertices_list$get_next()
  
    if (vertex$equals(graph$start_vertex())) {
      start <- create_vertex(graph, c(n, rep(0, n-1)))
      add_edge(graph$start_vertex(), start, 1)
      next()
    }
    
    for (i in 1:n) {
      for (j in i:n) {
        rate <- 0
        state <- vertex$state()

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
}

print(graph_as_matrix(graph))
*/

void *create_matrix(long double **mat, size_t length) {
  throw std::runtime_error(
      "Not implemented"
  );
  
  return NULL;
}
void *matrix_invert(void *matrix, size_t size) {
  throw std::runtime_error(
      "Not implemented"
  );
  
  return NULL;
}

double matrix_get(void *matrix, size_t i, size_t j) {
  throw std::runtime_error(
      "Not implemented"
  );
  
  return 0;
}

// [[Rcpp::export]]
ptdalgorithms::Graph create_graph(size_t state_length) {
  return Graph(
      state_length
  );
}

RCPP_EXPOSED_CLASS(ptdalgorithms::VertexLinkedList)
RCPP_EXPOSED_CLASS(ptdalgorithms::Vertex)
RCPP_EXPOSED_CLASS(ptdalgorithms::Graph)
  
RCPP_MODULE(ptdalgorithms) {
  class_<VertexLinkedList>("VertexLinkedList")
    .method("has_next", &ptdalgorithms::VertexLinkedList::has_next, "Does the list have a next value?")
    .method("get_next", &ptdalgorithms::VertexLinkedList::next, "Obtain the next vertex in list. Updates the list iterator.")
  ;
  
  class_<Vertex>("Vertex")
    .method("state", &ptdalgorithms::Vertex::state, "Obtain vertex state")
    .method("equals", &ptdalgorithms::Vertex::operator==, "Does this vertex equal another vertex?")
  ;
  
  class_<Graph>("Graph")
    .method("start_vertex", &ptdalgorithms::Graph::start_vertex, "Obtain graph start vertex")
    .method("vertices_list", &ptdalgorithms::Graph::vertices_list, "Obtain graph vertices list")
  ;
}

// TODO: add BEGIN_RCPP inside all functions

// [[Rcpp::export]]
void add_edge(Vertex phase_type_vertex_from, Vertex phase_type_vertex_to, double weight) {
  phase_type_vertex_from.add_edge(phase_type_vertex_to, weight);
}

// [[Rcpp::export]]
Vertex create_vertex(Graph phase_type_graph, IntegerVector state) {
  Vertex vertex = phase_type_graph.create_vertex(as<std::vector<size_t> >(state));
  
  return vertex;
}

// [[Rcpp::export]]
bool vertex_exists(Graph phase_type_graph, IntegerVector state) {
  return phase_type_graph.vertex_exists(as<std::vector<size_t> >(state));
}

// [[Rcpp::export]]
Vertex find_vertex(Graph phase_type_graph, IntegerVector state) {
  return phase_type_graph.find_vertex(as<std::vector<size_t> >(state));
}

// [[Rcpp::export]]
Vertex find_or_create_vertex(Graph phase_type_graph, IntegerVector state) {
  return phase_type_graph.find_or_create_vertex(as<std::vector<size_t> >(state));
}

List _graph_as_matrix(Graph graph) {
  PhaseTypeDistribution dist = graph.phase_type_distribution();
  
  NumericMatrix SIM(dist.length, dist.length);
  NumericVector IPV(dist.length);
  
  for (size_t i = 0; i < dist.length; ++i) {
    IPV(i) = dist.initial_probability_vector[i];
    
    for (size_t j = 0; j < dist.length; ++j) {
      SIM(i, j) = dist.sub_intensity_matrix[i][j];
    }
  }
  
  return List::create(Named("vertices") = dist.vertices , _["SIM"] = SIM, _["IPV"] = IPV);
}

// [[Rcpp::export]]
List graph_as_matrix(Graph phase_type_graph) {
  return(_graph_as_matrix(phase_type_graph));
}

