/*
 * Clone or download the code, and include these files in the repository!
 * Make SURE that the version of the downloaded code is the same as the
 * installed R library!! Otherwise it may crash randomly.
 *
 * The path is currently ../ as we are in the same repository. This path
 * should be something like [full or relative path to cloned code]/api...
 */
#include "../api/c/ptdalgorithms.h"

/*
 * Including a .c file is very strange usually!
 * But the way Rcpp::sourceCpp links is different from what
 * you would usually expect. Therefore this is by far
 * the easiest way of importing the code.
 */
#include "../src/c/ptdalgorithms.c"
#include "../api/cpp/ptdalgorithmscpp.h"

/* This is the binding layer such that R can invoke this function */
#include <Rcpp.h>

/* Basic C libraries */
#include "stdint.h"
#include "stdlib.h"

using namespace std;
using namespace ptdalgorithms;
using namespace Rcpp;

static inline int get_matrix_index(int i, int j, int population, int matrix_size) {
    return (matrix_size * matrix_size * population) + i + j * matrix_size;
}


// [[Rcpp::export]]
int matrix_index(int i, int j, int population, int n1, int n2) {
    const int matrix_size = (n1 > n2 ? n1 : n2) + 1;
    return (matrix_size * matrix_size * population) + i + j * matrix_size;
}

// [[Rcpp::export]]
IntegerVector rewards_at(SEXP phase_type_graph, int i, int j, int n1, int n2) {
    Rcpp::XPtr<Graph> graph(phase_type_graph);
    const int matrix_size = (n1 > n2 ? n1 : n2) + 1;
    IntegerVector res(graph->c_graph()->vertices_length);

    for (int k = 0; k < (int)graph->c_graph()->vertices_length; ++k) {
        int pop1 = graph->c_graph()->vertices[k]->state[get_matrix_index(i, j, 0, matrix_size)];
        int pop2 = graph->c_graph()->vertices[k]->state[get_matrix_index(i, j, 1, matrix_size)];
        res[k] = pop1 + pop2;
    }

    return res;
}

// [[Rcpp::export]]
NumericVector start_prob_from_im(SEXP a_phase_type_graph, SEXP im_phase_type_graph, NumericVector im_stopping) {
    Rcpp::XPtr <Graph> a_graph(a_phase_type_graph);
    Rcpp::XPtr <Graph> im_graph(im_phase_type_graph);

    NumericVector result(a_graph->c_graph()->vertices_length);
    size_t state_length = im_graph->state_length();
    struct ptd_avl_tree *avl_tree = a_graph->c_avl_tree();

    for (int k = 1; k < (int) im_graph->c_graph()->vertices_length; ++k) {
        int *a_state = (int *) calloc(state_length, sizeof(int));
        int *im_state = im_graph->c_graph()->vertices[k]->state;

        for (int i = 0; i < (int) state_length; ++i) {
            if (i >= (int) state_length / 2) {
                a_state[i - state_length / 2] += im_state[i];
            } else {
                a_state[i] += im_state[i];
            }
        }

        vector<int> s2(a_state, a_state + state_length / 2);
        struct ptd_vertex *vertex = (struct ptd_vertex *) ptd_avl_tree_find(avl_tree, a_state)->entry;
        result[vertex->index] += im_stopping[k];
        free(a_state);
    }

    return result;
}

// [[Rcpp::export]]
SEXP construct_ancestral_graph(int n1, int n2) {
    const int matrix_size = (n1 > n2 ? n1 : n2) + 1;
    const int state_length = matrix_size * matrix_size * 2;
    const size_t state_size = sizeof(int) * state_length;

    struct ptd_graph *ancestral_graph = ptd_graph_create((size_t) state_length);
    struct ptd_avl_tree *avl_tree = ptd_avl_tree_create((size_t) state_length);

    int *initial_state = (int *) calloc((size_t) state_length, sizeof(int));
    initial_state[get_matrix_index(1, 0, 0, matrix_size)] = n1;
    initial_state[get_matrix_index(0, 1, 0, matrix_size)] = n2;
    vector<int> s2(initial_state, initial_state + state_length / 2);
    ptd_graph_add_edge(
            ancestral_graph->starting_vertex,
            ptd_find_or_create_vertex(ancestral_graph, avl_tree, initial_state),
            1
    );
    free(initial_state);
    int *child_state = (int *) malloc(state_size);

    for (size_t index = 1; index < ancestral_graph->vertices_length; ++index) {
        struct ptd_vertex *vertex = ancestral_graph->vertices[index];
        int *state = vertex->state;

        int lineages_left = 0;
        for (int i = 0; i < state_length; ++i) {
            lineages_left += state[i];
        }

        if (lineages_left == 0 || lineages_left == 1) {
            // Only one lineage left, absorb
            continue;
        }

        for (int i = 0; i < matrix_size; ++i) {
            for (int j = 0; j < matrix_size; ++j) {
                int index1 = get_matrix_index(i, j, 0, matrix_size);
                int entry1 = state[index1];
                if (entry1 == 0) {
                    continue;
                }
                for (int i2 = 0; i2 < matrix_size; ++i2) {
                    for (int j2 = 0; j2 < matrix_size; ++j2) {
                        if (i + j * matrix_size < i2 + j2 * matrix_size) {
                            continue;
                        }
                        int index2 = get_matrix_index(i2, j2, 0, matrix_size);
                        int entry2 = state[index2];
                        if (entry2 == 0) {
                            continue;
                        }
                        double weight;
                        if (i == i2 && j == j2) {
                            if (entry1 == 1) {
                                continue;
                            }
                            weight = entry1 * (entry1 - 1) / 2;
                        } else {
                            weight = entry1 * entry2;
                        }

                        memcpy(child_state, state, state_size);
                        child_state[get_matrix_index(i, j, 0, matrix_size)]--;
                        child_state[get_matrix_index(i2, j2, 0, matrix_size)]--;
                        child_state[get_matrix_index(i + i2, j + j2, 0, matrix_size)]++;
                        struct ptd_vertex *child = ptd_find_or_create_vertex(
                                ancestral_graph, avl_tree, child_state
                        );
                        ptd_graph_add_edge(vertex, child, weight);
                    }
                }
            }
        }
    }
    free(child_state);

    return Rcpp::XPtr<Graph>(
            new Graph(ancestral_graph, avl_tree)
    );
}

// [[Rcpp::export]]
SEXP construct_im_graph(int n1, int n2, float m12, float m21) {
    const int matrix_size = (n1 > n2 ? n1 : n2) + 1;
    const int state_length = matrix_size * matrix_size * 2;
    const size_t state_size = sizeof(int) * state_length;

    struct ptd_graph *im_graph = ptd_graph_create((size_t) state_length);
    struct ptd_avl_tree *avl_tree = ptd_avl_tree_create((size_t) state_length);

    int *initial_state = (int *) calloc((size_t) state_length, sizeof(int));
    initial_state[get_matrix_index(1, 0, 0, matrix_size)] = n1;
    initial_state[get_matrix_index(0, 1, 1, matrix_size)] = n2;

    ptd_graph_add_edge(
            im_graph->starting_vertex,
            ptd_find_or_create_vertex(im_graph, avl_tree, initial_state),
            1
    );
    free(initial_state);
    int *child_state = (int *) malloc(state_size);

    for (size_t index = 1; index < im_graph->vertices_length; ++index) {
        struct ptd_vertex *vertex = im_graph->vertices[index];
        int *state = vertex->state;
        int lineages_left = 0;
        for (int i = 0; i < state_length; ++i) {
            lineages_left += state[i];
        }
        if (lineages_left == 0 || lineages_left == 1) {
            // Only one lineage left, absorb
            continue;
        }
        // Migration
        for (int i = 0; i < matrix_size; ++i) {
            for (int j = 0; j < matrix_size; ++j) {
                for (int population = 0; population <= 1; population++) {
                    int index1 = get_matrix_index(i, j, population, matrix_size);
                    int index2 = get_matrix_index(i, j, 1 - population, matrix_size);
                    if (state[index1] == 0) {
                        continue;
                    }
                    memcpy(child_state, state, state_size);
                    child_state[index1]--;
                    child_state[index2]++;
                    struct ptd_vertex *child = ptd_find_or_create_vertex(im_graph, avl_tree, child_state);
                    float migration_rate = population == 0 ? m21 : m12;
                    double *edge_state = (double *) calloc(3, sizeof(*edge_state));
                    edge_state[0] = 0;
                    edge_state[1] = 0;
                    edge_state[2] = 0;

                    if (population == 0) {
                        edge_state[2] = state[index1];
                    } else {
                        edge_state[1] = state[index1];
                    }

                    ptd_graph_add_edge_parameterized(
                            vertex, child, migration_rate * state[index1],
                            edge_state
                    );
                }
            }
        }

        for (int population = 0; population <= 1; population++) {
            for (int i = 0; i < matrix_size; ++i) {
                for (int j = 0; j < matrix_size; ++j) {
                    int index1 = get_matrix_index(i, j, population, matrix_size);
                    int entry1 = state[index1];
                    if (entry1 == 0) {
                        continue;
                    }
                    for (int i2 = 0; i2 < matrix_size; ++i2) {
                        for (int j2 = 0; j2 < matrix_size; ++j2) {
                            if (i + j * matrix_size < i2 + j2 * matrix_size) {
                                continue;
                            }
                            int index2 = get_matrix_index(i2, j2, population, matrix_size);
                            int entry2 = state[index2];
                            if (entry2 == 0) {
                                continue;
                            }
                            double weight;
                            if (i == i2 && j == j2) {
                                if (entry1 == 1) {
                                    continue;
                                }
                                weight = entry1 * (entry1 - 1) / 2;
                            } else {
                                weight = entry1 * entry2;
                            }
                            memcpy(child_state, state, state_size);
                            child_state[get_matrix_index(i, j, population, matrix_size)]--;
                            child_state[get_matrix_index(i2, j2, population, matrix_size)]--;
                            child_state[get_matrix_index(i + i2, j + j2, population, matrix_size)]++;
                            struct ptd_vertex *child = ptd_find_or_create_vertex(
                                    im_graph, avl_tree, child_state
                            );
                            double *edge_state = (double *) calloc(3, sizeof(*edge_state));
                            edge_state[0] = weight;
                            edge_state[1] = 0;
                            edge_state[2] = 0;
                            ptd_graph_add_edge_parameterized(vertex, child, weight, edge_state);
                        }
                    }
                }
            }
        }
    }
    free(child_state);
    return Rcpp::XPtr<Graph>(
            new Graph(im_graph, avl_tree)
    );
}
