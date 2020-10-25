#include <Rcpp.h>
#include "phase.h"
#include "io.h"

using namespace Rcpp;

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
    graph_as_mat(&mat, &size, &vertices, graph);

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
    }

    for (size_t i = 0; i < size; ++i) {
        free(mat[i]);
    }

    free(mat);
    free(vertices);

    graph_free(graph);

    return List::create(Named("IPV") = IPV , _["SIM"] = SIM, _["RW"] = RW);
}

size_t reward_index;

double reward_by_index(vertex_t *vertex) {
    return vertex->rewards[reward_index];
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
    gen_kingman_graph(&graph, n, r);
    reward_index = r-1;
    reward_transform(graph, reward_by_index);

    double **mat;
    vertex_t **vertices;
    size_t size;
    graph_as_mat(&mat, &size, &vertices, graph);

    NumericMatrix SIM(size - 2, size - 2);
    NumericVector IPV(size - 2);

    for (size_t i = 2; i < size; ++i) {
        IPV(i - 2) = mat[1][i];

        for (size_t j = 2; j < size; ++j) {
            SIM(i - 2, j - 2) = mat[i][j];
        }
    }

    for (size_t i = 0; i < size; ++i) {
        free(mat[i]);
    }

    free(mat);
    free(vertices);

    graph_free(graph);

    return List::create(Named("IPV") = IPV , _["SIM"] = SIM);
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
List kingman_gen_reward_by(int n, int m, Rcpp::Function reward_function) {
    if (n <= 0) {
        throw std::invalid_argument("'n' must be strictly positive");
    }
    if (m <= 0 || m > n) {
        throw std::invalid_argument("'m' must be strictly positive and no larger than 'n'");
    }

    vertex_t *graph;
    gen_kingman_graph(&graph, n, m);
    custom_reward_function = &reward_function;
    reward_transform(graph, custom_reward);

    double **mat;
    vertex_t **vertices;
    size_t size;
    graph_as_mat(&mat, &size, &vertices, graph);

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