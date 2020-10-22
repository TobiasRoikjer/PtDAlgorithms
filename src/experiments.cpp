#include <Rcpp.h>
#include "phase.h"
#include "io.h"
#include "generation.h"

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

class Graph {
public:
  Graph(vertex_t *vertex) {
    this->vertex = vertex;
  }
  
  ~Graph() {
    graph_free(this->vertex);
  }
  
  vertex_t *vertex;
};

size_t reward_index;

double reward_by_index(vertex_t *vertex) {
  return vertex->rewards[reward_index];
}

// [[Rcpp::export]]
List gen_mat(unsigned int seed, int n_states, int n_edges, int n_zero_rewards) {
  vertex_t *graph =
    generate_graph(seed, n_states, n_edges, n_zero_rewards);
  
  double **mat;
  size_t size;
  graph_as_mat(&mat, &size, graph);
  
  NumericMatrix SIM(size - 2, size - 2);
  NumericVector IPV(size - 2);
  NumericMatrix RW(size - 2, 1);
  
  for (size_t i = 2; i < size; ++i) {
    IPV(i - 2) = mat[1][i];
    
    for (size_t j = 2; j < size; ++j) {
      SIM(i - 2, j - 2) = mat[i][j];
    }

    if (i < n_zero_rewards) {
        RW(i - 2, 0) = 0;
    } else {
        RW(i - 2, 0) = 2;
    }
  }
  
  for (size_t i = 0; i < size; ++i) {
    free(mat[i]);
  }
  
  free(mat);
  graph_free(graph);
  
  return List::create(Named("IPV") = IPV , _["SIM"] = SIM, _["RW"] = RW);
}

// [[Rcpp::export]]
SEXP gen_c_mat(unsigned int seed, int n_states, int n_edges, int n_zero_rewards) {
  vertex_t *graph =
    generate_graph(seed, n_states, n_edges, n_zero_rewards);
  Graph *g = new Graph(graph);
  return Rcpp::XPtr<Graph>(g);
}

// [[Rcpp::export]]
SEXP reward_c_mat(SEXP pointer) {
  Rcpp::XPtr<Graph> graph(pointer);
  
  reward_index = 0;
  reward_transform(graph->vertex, reward_by_index);
  
  return pointer;
}

// [[Rcpp::export]]
List reward_test(List res) {
  NumericMatrix sim = res["SIM"];
  NumericVector ipv = res["IPV"];
  NumericMatrix rw = res["RW"];
  NumericVector rewards = rw.column(0);

  for (ssize_t u = sim.nrow() - 1; u >= 0; --u) {
    double reward = rewards[u];
    
    if (reward != 0) {
      sim(u,_) = sim(u,_) / reward;
    } else {
      NumericVector adding = sim.row(u) * ipv.at(u) / (-sim.at(u,u));
  
      ipv = ipv + adding;
      
      for (size_t z = 0; z < sim.cols(); ++z) {
        if (z == u) {
          continue;
        }
        
        sim(_,z) = sim(_,z) + (sim.at(u,z)/(-sim.at(u,u)))*sim.column(u);
      }
      
      sim(u,_) = sim(u,_) * 0;
      ipv.at(u) = 0;
      //SIM <- SIM[-u, -u]
      //IPV <- IPV[-u]
    }
  }
  
  return List::create(Named("IPV") = ipv , _["SIM"] = sim);
}