#define PTD_RCPP 1

#include "ptdalgorithms_types.h"
#include "../api/c/ptdalgorithms.h"
#include "../api/cpp/ptdalgorithmscpp.h"
#include "c/ptdalgorithms.cpp"
#include "cpp/ptdalgorithmscpp.cpp"

#include <Rcpp.h>
using namespace Rcpp;
using namespace ptdalgorithms;


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


SEXP get_first_list_entry(SEXP e, std::string message) {
  if (Rf_isList(e) || Rf_isNewList(e)) {
    List list = Rcpp::as<List>(e);
    
    if (list.size() != 1) {
      char message[1024];
      
      snprintf(
        message,
        1024, 
        "Failed: When finding %s of a list-type of vertices, list can only contain 1 vertex, it contained %zu. Did you use [] as lookup instead of [[]]?",
        message,
        (size_t)list.size()
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
  vector<size_t> state = vertex->state();
  IntegerVector state_vec(state.begin(), state.end());
  
  return List::create(
    Named("state") = state_vec,
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
      "Failed: child entries in children list must be a list, the datatype was '%i' (R internal type description)",
      (int)TYPEOF(list)
    );
    
    throw std::runtime_error(
        message
    );
  }
  
  Rcpp::List child = Rcpp::as<Rcpp::List>(list);
  
  return child["xptr_vertex"];
}

// [[Rcpp::export]]
List edges(SEXP phase_type_vertex) {
  phase_type_vertex = get_first_list_entry(phase_type_vertex, (char*)"edges");
  
  Rcpp::XPtr<Vertex> vertex(phase_type_vertex);
  vector<Edge> edges = vertex->edges();
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

// [[Rcpp::export]]
SEXP create_graph(size_t state_length) {
  return Rcpp::XPtr<Graph>(
    new Graph(
        state_length
    )
  );
}

// [[Rcpp::export]]
void add_edge(SEXP phase_type_vertex_from, SEXP phase_type_vertex_to, double weight) {
  phase_type_vertex_from = get_first_list_entry(phase_type_vertex_from, (char*)"edge from");
  phase_type_vertex_to = get_first_list_entry(phase_type_vertex_to, (char*)"edge to");
  
  Rcpp::XPtr<Vertex> from(phase_type_vertex_from);
  Rcpp::XPtr<Vertex> to(phase_type_vertex_to);
  
  from->add_edge(*to.get(), weight);
}

// [[Rcpp::export]]
SEXP create_vertex(SEXP phase_type_graph, IntegerVector state) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  Vertex *vertex = graph->create_vertex_p(as<std::vector<size_t> >(state));
  
  return vertex_as_list(vertex);
}

// [[Rcpp::export]]
void index_topological(SEXP phase_type_graph) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  graph->index_topological();
}

// [[Rcpp::export]]
void index_invert(SEXP phase_type_graph) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  graph->index_invert();
}

// [[Rcpp::export]]
SEXP find_vertex(SEXP phase_type_graph, IntegerVector state) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  if (!graph->vertex_exists(as<std::vector<size_t> >(state))) {
    return List::get_na();
  }
  
  Vertex *found = graph->find_vertex_p(as<std::vector<size_t> >(state));
  
  return vertex_as_list(found);
}


// [[Rcpp::export]]
List find_or_create_vertex(SEXP phase_type_graph, IntegerVector state) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  Vertex *found = graph->find_or_create_vertex_p(as<std::vector<size_t> >(state));
  
  return vertex_as_list(found);
}

// [[Rcpp::export]]
List start_vertex(SEXP phase_type_graph) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  Vertex *vertex = graph->start_vertex_p();
  
  return vertex_as_list(vertex);
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
  
  List vertices(dist.length);
  Graph g = *graph;
  for (size_t i = 0; i < dist.length; i++) {
    vertices[i] = vertex_as_list(&dist.vertices[i]);
  }
  
  return List::create(Named("vertices") = vertices , _["SIM"] = SIM, _["IPV"] = IPV);
}

// [[Rcpp::export]]
List graph_as_matrix(SEXP phase_type_graph) {
  return(_graph_as_matrix(phase_type_graph));
}


// [[Rcpp::export]]
SEXP vertices_list(SEXP phase_type_graph) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  return Rcpp::XPtr<VertexLinkedList>(
      graph->vertices_list_p()
  );
}

// [[Rcpp::export]]
bool list_has_next(SEXP vertex_list) {
  Rcpp::XPtr<VertexLinkedList> list(vertex_list);
  
  return list->has_next();
}

// [[Rcpp::export]]
List list_next(SEXP vertex_list) {
  Rcpp::XPtr<VertexLinkedList> list(vertex_list);
  
  return vertex_as_list(list->next_p());
}