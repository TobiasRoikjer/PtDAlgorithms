#define PTD_RCPP 1


#include <Rcpp.h>

using namespace Rcpp;

void *create_matrix(double **mat, size_t length) {
  double **res = (double**)calloc(length, sizeof(*res));
  
  for (size_t k = 0; k < length; ++k) {
    res[k] = (double*) calloc(length, sizeof(**res));
    
    for (size_t j = 0; j < length; ++j) {
      res[k][j] = mat[k][j];
    }
  }
  
  return (void*)res;
}

void *matrix_invert(void *matrix, size_t size) {
  double **mat = (double**)matrix;
  NumericMatrix r_mat(size, size);
  
  for (size_t k = 0; k < size; ++k) {
    for (size_t j = 0; j < size; ++j) {
      r_mat(k, j) = mat[k][j];
    }
  }
  
  Function f("solve");
  NumericMatrix inverted = f(r_mat);
  
  double **res = (double**)calloc(size, sizeof(*res));
  
  for (size_t k = 0; k < size; ++k) {
    res[k] = (double*) calloc(size, sizeof(**res));
    
    for (size_t j = 0; j < size; ++j) {
      res[k][j] = inverted(k,j);
    }
  }
  
  return (void*)res;
}

double matrix_get(void *matrix, size_t i, size_t j) {
  double **mat = (double**)matrix;
  
  return mat[i][j];
}

void matrix_set(void *matrix, size_t i, size_t j, double x) {
  double **mat = (double**)matrix;
  mat[i][j] = x;
}

void *matrix_init(size_t size) {
  double **res = (double**)calloc(size, sizeof(*res));
  
  for (size_t k = 0; k < size; ++k) {
    res[k] = (double*) calloc(size, sizeof(**res));
  }
  
  return (void*)res;
}

void matrix_destroy(void *matrix, size_t size) {
  double **mat = (double**)matrix;
  
  for (size_t k = 0; k < size; ++k) {
    free(mat[k]);
  }
  
  free(mat);
}

#include "ptdalgorithms_types.h"
#include "../api/c/ptdalgorithms.h"
#include "../api/cpp/ptdalgorithmscpp.h"
#include "c/ptdalgorithms.cpp"
#include "cpp/ptdalgorithmscpp.cpp"

#include <Rcpp.h>
using namespace Rcpp;
using namespace ptdalgorithms;


// TODO: Make all functions exists as Cpp method calls

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
  vector<int> state = vertex->state();
  IntegerVector state_vec(state.begin(), state.end());
  
  return List::create(
    Named("state") = state_vec,
    _["rate"] = vertex->rate(),
    _["vertex"] = (size_t)vertex->c_vertex(),
    _["xptr_vertex"] = Rcpp::XPtr<Vertex>(vertex)
  );
}

List vertex_as_list(Graph &graph, struct ptd_ph_vertex *c_vertex) {
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

// [[Rcpp::export]]
SEXP create_graph(size_t state_length) {
  return Rcpp::XPtr<Graph>(
    new Graph(
        state_length
    )
  );
}

// [[Rcpp::export]]
List edges(SEXP phase_type_vertex) {
  phase_type_vertex = list_as_vertex(phase_type_vertex);
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
List vertices(SEXP phase_type_graph) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  std::vector<Vertex*> vertices = graph->vertices_p();
  List res(vertices.size());
  
  for (size_t i = 0; i < vertices.size(); i++) {
    res[i] = vertex_as_list(vertices[i]);
  }
  
  return res;
}

// [[Rcpp::export]]
int vertices_length(SEXP phase_type_graph) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  return (int)graph->c_graph()->vertices_length;
}

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
  
  ptd_ph_vertex *vertex = graph->c_graph()->vertices[index - 1];
  
  return vertex_as_list(new Vertex(*graph, vertex));
}



// [[Rcpp::export]]
void add_edge(SEXP phase_type_vertex_from, SEXP phase_type_vertex_to, double weight) {
  phase_type_vertex_from = list_as_vertex(phase_type_vertex_from);
  phase_type_vertex_to = list_as_vertex(phase_type_vertex_to);
  
  Rcpp::XPtr<Vertex> from(phase_type_vertex_from);
  Rcpp::XPtr<Vertex> to(phase_type_vertex_to);
  
  from->add_edge(*to.get(), weight);
}

// [[Rcpp::export]]
List starting_vertex(SEXP phase_type_graph) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  Vertex *vertex = graph->starting_vertex_p();
  
  return vertex_as_list(vertex);
}

// [[Rcpp::export]]
SEXP create_vertex(SEXP phase_type_graph, IntegerVector state) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  Vertex *vertex = graph->create_vertex_p(as<std::vector<int> >(state));
  
  return vertex_as_list(vertex);
}

// [[Rcpp::export]]
SEXP find_vertex(SEXP phase_type_graph, IntegerVector state) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  if (!graph->vertex_exists(as<std::vector<int> >(state))) {
    return List::get_na();
  }
  
  Vertex *found = graph->find_vertex_p(as<std::vector<int> >(state));
  
  return vertex_as_list(found);
}


// [[Rcpp::export]]
List find_or_create_vertex(SEXP phase_type_graph, IntegerVector state) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  Vertex *found = graph->find_or_create_vertex_p(as<std::vector<int> >(state));
  
  return vertex_as_list(found);
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

// [[Rcpp::export]]
List graph_as_matrix(SEXP phase_type_graph) {
  return(_graph_as_matrix(phase_type_graph));
}

// [[Rcpp::export]]
void reward_transform(SEXP phase_type_graph, NumericVector rewards) {
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
  
  graph->reward_transform(as<std::vector<double> >(rewards));
}

// [[Rcpp::export]]
NumericVector expected_visits(SEXP phase_type_graph) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  return wrap(graph->expected_visits());
}

// [[Rcpp::export]]
NumericVector expected_waiting_time(SEXP phase_type_graph) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  return wrap(graph->expected_waiting_time());
}

// [[Rcpp::export]]
NumericVector moment_rewards(SEXP phase_type_graph, NumericVector rewards) {
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
  
  return wrap(graph->moment_rewards(as<std::vector<double> >(rewards)));
}


// [[Rcpp::export]]
bool graph_is_acyclic(SEXP phase_type_graph) {
  Rcpp::XPtr<Graph> graph(phase_type_graph);
  
  return graph->is_acyclic();
}
