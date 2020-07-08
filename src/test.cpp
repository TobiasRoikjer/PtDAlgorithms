#include <cstddef>
#include <cstdio>
#include <algorithm>
#include "phase.h"

int main(int argv, char **argc) {
    size_t n = 5;
    size_t m = 5;

    vertex_t *graph;
    gen_kingman_graph(&graph, n, m);

    double **mat;
    vertex_t **vertices;
    size_t size;
    graph_as_mat(&mat, &vertices, &size, graph);

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
    size_t n = 10;
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