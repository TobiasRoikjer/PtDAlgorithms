#include <Rcpp.h>
#include "phase.h"

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

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
