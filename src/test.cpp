#include <cstddef>
#include <cstdio>
#include <algorithm>
#include "phase.h"


thread_local size_t r;

double reward_by_index(vertex_t *vertex) {
    return vertex->rewards[r];
}

int main(int argv, char **argc) {
    size_t n = 7;
    size_t m = n;

    vertex_t *graph;
    gen_kingman_graph(&graph, n, m);
    r = 2;
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
    size_t n = 25;
    size_t m = 5;

    vertex_t *graph;
    gen_kingman_graph(&graph, n, m);
    cov_exp_return ret = mph_cov_exp_all(graph, m);

    for (size_t i = 0; i < (size_t)m; ++i) {
        for (size_t j = 0; j < (size_t)m; ++j) {
            fprintf(stdout, "Cov %zu %zu %f \n",
                    i, j, (*ret.cov)[max(i,j)][min(i,j)]);
        }
    }

    graph_free(graph);

    return 0;
}