#include <Rcpp.h>
#include "phase.h"
#include "io.h"

using namespace Rcpp;

class PhaseTypeGraph {
public:
    PhaseTypeGraph(vertex_t *start, size_t rewards_length) {
        this->start = start;
        this->rewards_length = rewards_length;
    }

    ~PhaseTypeGraph() {
        graph_free(this->start);
    }

    vertex_t *start;
    size_t rewards_length;
};

// [[Rcpp::export]]
List graph_as_matrix(SEXP phase_type_graph) {
    Rcpp::XPtr<PhaseTypeGraph> graph(phase_type_graph);
    double **mat;
    size_t size;
    vertex_t **vertices;
    graph_as_mat(&mat, &size, &vertices, graph->start);

    NumericMatrix SIM(size - 2, size - 2);
    NumericVector IPV(size - 2);
    NumericMatrix RW(size - 2, graph->rewards_length);

    for (size_t i = 2; i < size; ++i) {
        IPV(i - 2) = mat[1][i];

        for (size_t j = 2; j < size; ++j) {
            SIM(i - 2, j - 2) = mat[i][j];

            if (j - 2 < graph->rewards_length) {
                RW(i - 2, j - 2) = vertices[i]->rewards[j - 2];
            }
        }
    }

    for (size_t i = 2; i < size; ++i) {
        vertex_t *vertex = vertices[i];
      
        for (size_t j = 0; j < graph->rewards_length; j++) {
            RW(i - 2, j) = vertex->rewards[j];
        }
    }
    
    for (size_t i = 0; i < size; ++i) {
        free(mat[i]);
    }

    free(mat);
    free(vertices);

    return List::create(Named("IPV") = IPV , _["SIM"] = SIM, _["RW"] = RW);
}

// [[Rcpp::export]]
SEXP matrix_as_graph(NumericVector initial_probability_vector,NumericMatrix subintensity_matrix,NumericMatrix rewards) {
  size_t entries = (size_t)initial_probability_vector.length();
  
  if (subintensity_matrix.rows() != subintensity_matrix.cols()) {
    throw std::invalid_argument("'subintensity_matrix' must be a square matrix.");
  }
  
  if (subintensity_matrix.rows() != initial_probability_vector.length()) {
    throw std::invalid_argument("'subintensity_matrix' must have the same size as 'initial_probability_vector'.");
  }
  
  if (rewards.rows() != initial_probability_vector.length()) {
    if (!Rcpp::NumericVector::is_na(rewards.at(0,0))) {
      throw std::invalid_argument("'rewards' must have the same number of rows as the length of 'initial_probability_vector'.");
    }
  }
  
  size_t reward_size;
  
  if (Rcpp::NumericVector::is_na(rewards.at(0,0))) {
    reward_size = 0;
  } else {
    reward_size = (size_t)rewards.cols();
  }
  
  vector<double> zero_reward(reward_size);
  
  vertex_t *absorbing = vertex_init(NULL, zero_reward, reward_size);
  vertex_t *start = vertex_init(NULL, zero_reward, reward_size);
  vertex_t **vertices = (vertex_t**)calloc(entries, sizeof(vertex_t*));
  
  for (size_t i = 0; i < entries; ++i) {
    vector<double> reward;
    
    if (reward_size != 0) {
      reward = Rcpp::as<std::vector<double> >(NumericVector(rewards.row(i)));
    }
    
    vertices[i] = vertex_init(NULL, reward, reward_size);
  }
  
  for (size_t i = 0; i < entries; ++i) {
    double total_rate = -subintensity_matrix(i, i);
    double non_absorbing_rate = 0;
    
    for (size_t j = 0; j < entries; ++j) {
      if (j == i) {
        continue;
      }
      
      if (subintensity_matrix(i, j) >= 0.0000001) {
        non_absorbing_rate += subintensity_matrix(i, j);
        vertex_add_edge(vertices[i], vertices[j], subintensity_matrix(i, j));
      }
    }
    
    if (total_rate - non_absorbing_rate >= 0.0001) {
      vertex_add_edge(vertices[i], absorbing, total_rate - non_absorbing_rate);
    }
  }
  
  for (size_t i = 0; i < entries; ++i) {
    if (initial_probability_vector.at(i) >= 0.0000001) {
      vertex_add_edge(start, vertices[i], initial_probability_vector.at(i));
    }
  }
  
  free(vertices);
  
  return Rcpp::XPtr<PhaseTypeGraph>(new PhaseTypeGraph(start, reward_size));
}


// [[Rcpp::export]]
SEXP kingman_gen_graph(int n, int m) {
    if (n <= 0) {
        throw std::invalid_argument("'n' must be strictly positive");
    }
    if (m <= 0 || m > n) {
        throw std::invalid_argument("'m' must be strictly positive and no larger than 'n'");
    }

    vertex_t *graph;
    gen_kingman_graph(&graph, n, m);

    return Rcpp::XPtr<PhaseTypeGraph>(new PhaseTypeGraph(graph, m));
}

size_t reward_index;

double reward_by_index(vertex_t *vertex) {
    return vertex->rewards[reward_index];
}

Rcpp::Function *custom_reward_function;

double custom_reward(vertex_t *vertex) {
    NumericVector rewards(vertex->rewards.size());

    for (ssize_t i = 0; i < rewards.size(); ++i) {
        rewards(i) = vertex->rewards[i];
    }

    return Rcpp::as<double>((*custom_reward_function)(rewards));
}

// [[Rcpp::export]]
SEXP reward_transform(SEXP phase_type_graph, Rcpp::Function reward_function) {
    Rcpp::XPtr<PhaseTypeGraph> graph(phase_type_graph);
    custom_reward_function = &reward_function;
    reward_transform(graph->start, custom_reward);

    return Rcpp::XPtr<PhaseTypeGraph>(graph);
}

// [[Rcpp::export]]
SEXP reward_transform_by_reward(SEXP phase_type_graph, int reward) {
  Rcpp::XPtr<PhaseTypeGraph> graph(phase_type_graph);
  reward_index = reward - 1;
  reward_transform(graph->start, &reward_by_index);
  
  return Rcpp::XPtr<PhaseTypeGraph>(graph);
}

// [[Rcpp::export]]
List graph_exp_cov(SEXP phase_type_graph) {
    Rcpp::XPtr<PhaseTypeGraph> graph(phase_type_graph);
    cov_exp_return ret = mph_cov_exp_all(graph->start, graph->rewards_length);

    NumericMatrix out_cov(graph->rewards_length, graph->rewards_length);
    NumericVector out_exp(graph->rewards_length);

    for (size_t i = 0; i < graph->rewards_length; ++i) {
        out_exp(i) = (*ret.exp)[i];

        for (size_t j = 0; j < graph->rewards_length; ++j) {
            out_cov(i,j) = (*ret.cov)[max(i,j)][min(i,j)];
        }
    }

    return List::create(Named("exp") = out_exp , _["cov"] = out_cov);
}

// [[Rcpp::export]]
List graph_info(SEXP phase_type_graph) {
  Rcpp::XPtr<PhaseTypeGraph> graph(phase_type_graph);
  struct graph_info graph_info = get_graph_info(graph->start);
  
  
  return List::create(Named("vertices") = graph_info.vertices , _["edges"] = graph_info.edges);
}

Rcpp::Function *custom_set_reward_function;

vector<double> custom_set_reward(vector<double> rewards) {
  NumericVector numeric_rewards = wrap(rewards);
    
  NumericVector new_rewards = Rcpp::as<NumericVector >((*custom_set_reward_function)(numeric_rewards));
  
  vector<double> vec_double = Rcpp::as<vector<double> >(new_rewards);
  
  return vec_double;
}

// [[Rcpp::export]]
SEXP set_rewards(SEXP phase_type_graph, Rcpp::Function set_reward_function) {
  Rcpp::XPtr<PhaseTypeGraph> graph(phase_type_graph);

  custom_set_reward_function = &set_reward_function;
  set_graph_rewards(graph->start, custom_set_reward);
  
  graph->rewards_length = graph->start->rewards.size();
  
  return Rcpp::XPtr<PhaseTypeGraph>(graph);
}

// TODO:
// Clone graph when e.g. reward transforming
