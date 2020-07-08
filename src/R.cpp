#include <Rcpp.h>
#include "phase.h"

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
List kingman_gen_mat(int n, int m) {
  if (n <= 0) {
    throw std::invalid_argument("'n' must be strictly positive");
  }
  if (m <= 0 || m > n) {
    throw std::invalid_argument("'m' must be strictly positive and no larger than 'n'");
  }
  
  vertex_t *graph;
  gen_kingman_graph(&graph, n, m);
  
  double **mat;
  vertex_t **vertices;
  size_t size;
  graph_as_mat(&mat, &vertices, &size, graph);
  
  NumericMatrix SIM(size - 2, size - 2);
  NumericVector IPV(size - 2);
  NumericMatrix RW(size - 2, m);
  
  for (size_t i = 2; i < size; ++i) {
    IPV(i - 2) = mat[1][i];
    
    for (size_t j = 2; j < size; ++j) {
      SIM(i - 2, j - 2) = mat[i][j];
      
      if (j - 2 < (size_t)m) {
        RW(i - 2, j - 2) = vertices[i]->rewards[j - 2];
      }
    }
    
    fprintf(stdout, "\n");
  }
  
  for (size_t i = 0; i < size; ++i) {
    free(mat[i]);
  }
  
  free(mat);
  free(vertices);
  
  graph_free(graph);
  
  return List::create(Named("IPV") = IPV , _["SIM"] = SIM, _["RW"] = RW);
}

// [[Rcpp::export]]
List kingman_gen_reward(int n, int r) {
  if (n <= 0) {
    throw std::invalid_argument("'n' must be strictly positive");
  }
  if (r <= 0 || r > n) {
    throw std::invalid_argument("'r' must be strictly positive and no larger than 'n'");
  }
  
  vertex_t *graph;
  gen_kingman_graph(&graph, n, r+1);
  reward_transform(graph, r-1);
  
  double **mat;
  vertex_t **vertices;
  size_t size;
  graph_as_mat(&mat, &vertices, &size, graph);
  
  NumericMatrix SIM(size - 2, size - 2);
  NumericVector IPV(size - 2);
  
  for (size_t i = 2; i < size; ++i) {
    IPV(i - 2) = mat[1][i];
    
    for (size_t j = 2; j < size; ++j) {
      SIM(i - 2, j - 2) = mat[i][j];
    }
    
    fprintf(stdout, "\n");
  }
  
  for (size_t i = 0; i < size; ++i) {
    free(mat[i]);
  }
  
  free(mat);
  free(vertices);
  
  graph_free(graph);
  
  return List::create(Named("IPV") = IPV , _["SIM"] = SIM);
}

// [[Rcpp::export]]
List kingman_exp_cov(int n, int m) {
  if (n <= 0) {
    throw std::invalid_argument("'n' must be strictly positive");
  }
  if (m <= 0 || m > n) {
    throw std::invalid_argument("'m' must be strictly positive and no larger than 'n'");
  }
  
  vertex_t *graph;
  gen_kingman_graph(&graph, n, m);
  cov_exp_return ret = mph_cov_exp_all(graph, m);
  
  NumericMatrix out_cov(m, m);
  NumericVector out_exp(m);
  
  for (size_t i = 0; i < (size_t)m; ++i) {
    out_exp(i) = (*ret.exp)[i];
    
    for (size_t j = 0; j < (size_t)m; ++j) {
      out_cov(i,j) = (*ret.cov)[max(i,j)][min(i,j)];
    }
  }
  
  graph_free(graph);
  
  return List::create(Named("exp") = out_exp , _["cov"] = out_cov);
}
