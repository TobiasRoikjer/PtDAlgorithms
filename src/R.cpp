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
void label_graph(SEXP phase_type_graph) {
  Rcpp::XPtr<PhaseTypeGraph> graph(phase_type_graph);
  label_vertex_index(NULL, graph->start);
}

// [[Rcpp::export]]
NumericVector probs(SEXP phase_type_graph) {
  Rcpp::XPtr<PhaseTypeGraph> graph(phase_type_graph);
  
  double *probs;
  size_t size;
  
  calculate_prob(graph->start, &size, &probs);
  NumericVector res(size - 2);
  
  for (size_t i = 2; i < size; ++i) {
    res[i - 2] = probs[i];
  }
  
  free(probs);
  
  return res;
}

// [[Rcpp::export]]
double graph_var(SEXP phase_type_graph) {
  Rcpp::XPtr<PhaseTypeGraph> graph(phase_type_graph);
  
  double *vars;
  size_t size;
  
  calculate_var(graph->start, &size, &vars);
  
  double var = vars[1];
  
  free(vars);
  
  return var;
}

List _graph_as_matrix(SEXP phase_type_graph, bool transition_probs) {
    Rcpp::XPtr<PhaseTypeGraph> graph(phase_type_graph);
    double **mat;
    size_t size;
    vertex_t **vertices;
    
    if (!transition_probs) {
      graph_as_mat(&mat, &size, &vertices, graph->start);
    } else {
      graph_as_t_mat(&mat, &size, &vertices, graph->start);
    }

    NumericMatrix SIM(size - 2, size - 2);
    NumericVector IPV(size - 2);
    NumericMatrix RW(size - 2, graph->rewards_length);

    for (size_t i = 2; i < size; ++i) {
      IPV(i - 2) = mat[1][i];

      for (size_t j = 2; j < size; ++j) {
        SIM(i - 2, j - 2) = mat[i][j];
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
List graph_as_matrix(SEXP phase_type_graph) {
  return(_graph_as_matrix(phase_type_graph, false));
}

// [[Rcpp::export]]
List graph_as_t_matrix(SEXP phase_type_graph) {
  return(_graph_as_matrix(phase_type_graph, true));
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

// [[Rcpp::export]]
SEXP kingman_gen_graph_rw(int n, int rw) {
  if (n <= 0) {
    throw std::invalid_argument("'n' must be strictly positive");
  }
  
  vertex_t *graph;
  gen_kingman_graph_rw(&graph, n, rw - 1);
  
  return Rcpp::XPtr<PhaseTypeGraph>(new PhaseTypeGraph(graph, rw - 1 + 2));
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
    NumericVector out_defect(graph->rewards_length);

    for (size_t i = 0; i < graph->rewards_length; ++i) {
        out_exp(i) = (*ret.exp)[i];
        out_defect(i) = (*ret.defect)[i];

        for (size_t j = 0; j < graph->rewards_length; ++j) {
            out_cov(i,j) = (*ret.cov)[max(i,j)][min(i,j)];
        }
    }

    return List::create(Named("exp") = out_exp , _["cov"] = out_cov, _["defect"] = out_defect);
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

size_t custom_state_length;

inline pair<double, vector<size_t> > parse_child(SEXP childS) {
  if (!Rf_isList(childS) && !Rf_isNewList(childS)) {
    char message[1024];
    
    snprintf(
      message,
      1024, 
      "Failed: child entries in children list must be a list, the datatype was '%i' (R internal type description)",
      (int)TYPEOF(childS)
    );
    
    throw std::runtime_error(
        message
    );
  }
  
  Rcpp::List child = Rcpp::as<Rcpp::List>(childS);
  
  if (!child.containsElementNamed("weight")) {
    throw std::runtime_error(
        "Failed: the list of children must have 2 named entries (weight, state), 'weight' not found"
    );
  }
  
  if (!child.containsElementNamed("state")) {
    throw std::runtime_error(
        "Failed: the list of children must have 2 named entries (weight, state), 'state' not found"
    );
  }
  
  double weight = Rcpp::as<double>(child["weight"]);
  
  if (weight <= 0) {
    char message[1024];
    
    snprintf(
      message,
      1024, 
      "Failed: weight is not a strictly positive real (was %f)",
      weight
    );
    
    throw std::runtime_error(
        message
    );
  }
  
  SEXP stateS = child["state"];
  
  if (!Rf_isNumber(stateS)) {
    throw std::runtime_error(
        "Failed: given state is not a number (vector)"
    );
  }
  
  IntegerVector state = Rcpp::as<Rcpp::IntegerVector>(stateS);
  
  if ((int)state.size() != (int)custom_state_length) {
    char message[1024];
    
    snprintf(
      message,
      1024, 
      "Failed: state must have length '%i', state had length '%i'",
      (int)custom_state_length,
      (int)state.size()
    );
    
    throw std::runtime_error(
        message
    );
  }
  
  vector<size_t> child_state = Rcpp::as<vector<size_t> >(state);
  
  return pair<double, vector<size_t> >(weight, child_state);
}

Rcpp::Function *custom_initial_function_R;

vector<pair<double, vector<size_t> > > custom_initial_function() {
  SEXP childrenS = (*custom_initial_function_R)();
  
  if (!Rf_isList(childrenS) && !Rf_isNewList(childrenS)) {
    char message[1024];
    
    snprintf(
      message,
      1024, 
      "Failed: children must be returned as a list, the datatype was '%i' (R internal type description)",
      (int)TYPEOF(childrenS)
    );
    
    throw std::runtime_error(
        message
    );
  }
  
  List children = Rcpp::as<Rcpp::List>(childrenS);
  vector<pair<double, vector<size_t> > > res;
  
  for (size_t i = 0; i < (size_t)children.size(); i++) {
    res.push_back(parse_child(children[i]));
  }
  
  return res;
}

Rcpp::Function *custom_children_function_R;

vector<pair<double, vector<size_t> > > custom_children_function(vector<size_t> state) {
  SEXP childrenS = ((*custom_children_function_R)(state));

  if (!Rf_isList(childrenS) && !Rf_isNewList(childrenS)) {
    char message[1024];
    
    snprintf(
      message,
      1024, 
      "Failed: children must be returned as a list, the datatype was '%i' (R internal type description)",
      (int)TYPEOF(childrenS)
    );
    
    throw std::runtime_error(
        message
    );
  }
  
  List children = Rcpp::as<Rcpp::List>(childrenS);
  vector<pair<double, vector<size_t> > > res;
  
  for (size_t i = 0; i < (size_t)children.size(); i++) {
    res.push_back(parse_child(children[i]));
  }
  
  return res;
}

Rcpp::Function *custom_rewards_from_state_function_R;

vector<double> custom_rewards_from_state_function(vector<size_t> state) {
  NumericVector rewards = Rcpp::as<Rcpp::NumericVector>((*custom_rewards_from_state_function_R)(state));
  
  vector<double> res = Rcpp::as<vector<double> >(rewards);
  
  return res;
}

// [[Rcpp::export]]
SEXP generate_graph(int state_length, Rcpp::Function initial, Rcpp::Function children, Rcpp::Function rewards_from_state) {
  custom_state_length = state_length;
  custom_initial_function_R = &initial;
  custom_children_function_R = &children;
  custom_rewards_from_state_function_R = &rewards_from_state;
  
  vertex_t *graph = generate_state_space(
    (size_t)state_length,
    custom_children_function,
    custom_initial_function,
    custom_rewards_from_state_function
  );
  
  return Rcpp::XPtr<PhaseTypeGraph>(new PhaseTypeGraph(graph, graph->rewards.size()));
}


// [[Rcpp::export]]
List graph_apply(SEXP phase_type_graph, Rcpp::Function func) {
  Rcpp::XPtr<PhaseTypeGraph> graph(phase_type_graph);
  queue<vertex_t *> q = enqueue_vertices(graph->start);

  List list(q.size() - 2);
  
  while (!q.empty()) {
    vertex_t *vertex = q.front();
    q.pop();
    
    if (vertex->vertex_index <= 1) {
      continue;
    }
    
    NumericVector rewards(vertex->rewards.size());
    
    for (ssize_t i = 0; i < rewards.size(); ++i) {
      rewards(i) = vertex->rewards[i];
    }
    
    list[vertex->vertex_index - 2] = func(
      vertex->vertex_index - 1,
      vertex->rate,
      rewards
    );
  }
  
  return list;
}

// TODO:
// Clone graph when e.g. reward transforming


size_t n1, n2, nrow, ncol, t1, t2, matrix_size, state_size;
double m;

void mig_set_state_properties(size_t num_pop1, size_t num_pop2,
                              size_t types_pop1, size_t types_pop2, double mig) {
  n1 = num_pop1;
  n2 = num_pop2;
  nrow = types_pop1 + 1;
  ncol = types_pop2 + 1;
  matrix_size = nrow * ncol;
  state_size = matrix_size * 2;
  t1 = types_pop1;
  t2 = types_pop2;
  m = mig;
}

size_t mig_state_index(size_t p1, size_t p2, size_t population) {
  if (p1 > t1) {
    p1 = t1;
  }
  
  if (p2 > t2) {
    p2 = t2;
  }
  
  return nrow * (p2) + p1 + (population == 1 ? 0 : matrix_size);
}

vector<pair<double, vector<size_t> > > mig_visit_function(vector<size_t> state) {
  vector<pair<double, vector<size_t> > > children;
  
  if (accumulate(state.begin(), state.end(), 0) == 1) {
    return children;
  }
  
  for (size_t p = 1; p <= 2; p++) {
    for (size_t i1 = 0; i1 <= t1; i1++) {
      for (size_t j1 = 0; j1 <= t2; j1++) {
        size_t idx1 = mig_state_index(i1, j1, p);
        
        {
          size_t idx2 = mig_state_index(i1, j1, 2 - p + 1);
          // Migration from the islands
          if (state[idx1] >= 1) {
            double rate = m * state[idx1];
            
            vector<size_t> child_state(state);
            
            child_state[idx1]--;
            child_state[idx2]++;
            
            children.push_back(
              pair<double, vector<size_t> >(rate, child_state)
            );
          }
        }
        
        if (state[idx1] != 0) {
          // Coalescence
          for (size_t i2 = 0; i2 <= t1; i2++) {
            for (size_t j2 = 0; j2 <= t2; j2++) {
              double rate;
              size_t idx2 = mig_state_index(i2, j2, p);
              
              if (i1 == i2 && j1 == j2) {
                if (state[idx1] <= 1) {
                  continue;
                }
                
                rate = (double) (state[idx1] * (state[idx1] - 1)) / 2;
              } else {
                if (state[idx1] == 0 || state[idx2] == 0) {
                  continue;
                }
                
                rate = (double) state[idx1] * state[idx2];
              }
              
              vector<size_t> child_state(state);
              size_t idx3 = mig_state_index(i1 + i2, j1 + j2, p);
              
              child_state[idx1]--;
              child_state[idx2]--;
              child_state[idx3]++;
              
              children.push_back(
                pair<double, vector<size_t> >(rate, child_state)
              );
            }
          }
        }
      }
    }
  }
  
  return (children);
}

vector<pair<double, vector<size_t> > > mig_initial_states(void) {
  vector<pair<double, vector<size_t> > > children;
  vector<size_t> state = vector<size_t>(state_size, 0);
  
  state[mig_state_index(1, 0, 1)] = n1;
  state[mig_state_index(0, 1, 2)] = n2;
  
  children.push_back(pair<double, vector<size_t> >(1.0f, state));
  
  return children;
}

vector<double> mig_rewards(vector<size_t> input) {
  vector<double> rew(input.begin(), input.end());
  
  return rew;
}

// [[Rcpp::export]]
SEXP generate_graphM() {
  mig_set_state_properties(3, 3, 3, 1, 0.1);
  
  vertex_t *graph = generate_state_space(
    state_size,
    mig_visit_function,
    mig_initial_states,
    mig_rewards
  );
  
  return Rcpp::XPtr<PhaseTypeGraph>(new PhaseTypeGraph(graph, graph->rewards.size()));
}


// [[Rcpp::export]]
SEXP generate_graphM2(int n1, int n2, int t1, int t2, double m) {
  mig_set_state_properties((size_t)n1, (size_t)n2, (size_t)t1, size_t(t2), m);
  
  vertex_t *graph = generate_state_space(
    state_size,
    mig_visit_function,
    mig_initial_states,
    mig_rewards
  );
  
  return Rcpp::XPtr<PhaseTypeGraph>(new PhaseTypeGraph(graph, graph->rewards.size()));
}
