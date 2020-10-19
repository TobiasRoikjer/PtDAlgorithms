#include <cstddef>
#include <cstdio>
#include <algorithm>
#include "generation.h"

size_t r;

double reward_by_index(vertex_t *vertex) {
    return vertex->rewards[r];
}

int main5(int argv, char **argc) {
    size_t n_states = 2046;//atoi(argc[1]);
    size_t n_edges = (size_t) 50 * 50;//atoi(argc[2]);
    size_t n_zero_rewards = (size_t) (2022 * 0.5);//atoi(argc[3]);

    vertex_t *graph = generate_graph(80002020, n_states, n_edges, n_zero_rewards);
    r = 0;
    reward_transform(graph, reward_by_index);

    double **mat;
    vertex_t **vertices;
    size_t size;

/*    graph_as_mat(&mat, &vertices, &size, graph);

    for (size_t i = 2; i < size; ++i) {
        fprintf(stdout, "%f ", mat[1][i]);
    }

    fprintf(stdout, "\n\n");

    for (size_t i = 2; i < size; ++i) {
        for (size_t j = 2; j < size; ++j) {
            fprintf(stdout, "%f ", mat[i][j]);
        }

        fprintf(stdout, "\n");
    }

    for (size_t i = 0; i < size; ++i) {
        free(mat[i]);
    }

    free(mat);
    free(vertices);
    */

    graph_free(graph);

    return 0;
}

int main(int argc, char **argv) {
    size_t n = 50;
    size_t m = 10;

    vertex_t *graph;
    //gen_kingman_graph(&graph, n, m);
    graph = generate_graph(1234,
                           10,
                           5,
                           0);
    fprintf(stdout, "Done generating\n");
    r = 0;
    reward_transform(graph, reward_by_index);

    double **mat;
    vertex_t **vertices;
    size_t size;

    graph_as_mat(&mat, &vertices, &size, graph);

    for (size_t i = 2; i < size; ++i) {
        fprintf(stdout, "%f ", mat[1][i]);
    }

    fprintf(stdout, "\n\n");

    for (size_t i = 2; i < size; ++i) {
        for (size_t j = 2; j < size; ++j) {
            fprintf(stdout, "%f ", mat[i][j]);
        }

        fprintf(stdout, "\n");
    }

    for (size_t i = 0; i < size; ++i) {
        free(mat[i]);
    }

    free(mat);
    free(vertices);

    graph_free(graph);

    return 0;
}


int main3(int argv, char **argc) {
    size_t n = atoi(argc[1]);
    size_t m = n;

    vertex_t *graph;
    gen_kingman_graph(&graph, n, m);
    r = 0;
    reward_transform(graph, reward_by_index);

    double **mat;
    vertex_t **vertices;
    size_t size;

    graph_as_mat(&mat, &vertices, &size, graph);

    for (size_t i = 2; i < size; ++i) {
        fprintf(stdout, "%f ", mat[1][i]);
    }

    fprintf(stdout, "\n\n");

    for (size_t i = 2; i < size; ++i) {
        for (size_t j = 2; j < size; ++j) {
            fprintf(stdout, "%f ", mat[i][j]);
        }

        fprintf(stdout, "\n");
    }

    for (size_t i = 0; i < size; ++i) {
        free(mat[i]);
    }

    free(mat);
    free(vertices);

    graph_free(graph);

    return 0;
}

int main2(int argv, char **argc) {
    size_t n = atoi(argc[1]);
    size_t m = atoi(argc[2]);

    vertex_t *graph;
    gen_kingman_graph(&graph, n, m);
    cov_exp_return ret = mph_cov_exp_all(graph, m);

    for (size_t i = 0; i < (size_t) m; ++i) {
        for (size_t j = 0; j < (size_t) m; ++j) {
            fprintf(stdout, "Cov %zu %zu %f \n",
                    i, j, (*ret.cov)[max(i, j)][min(i, j)]);
        }
    }

    graph_free(graph);

    return 0;
}