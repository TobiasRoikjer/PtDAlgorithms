#include "io.h"

/*
 * Assumes the graph has been labelled, as they should all be.
 */
int graph_as_mat(double ***weights, size_t *out_size, vertex_t ***vertices, const vertex_t *graph) {
    queue<vertex_t*> q = enqueue_vertices((vertex_t*)graph);

    size_t size = (size_t)q.size();
    *out_size = size;

    *vertices = (vertex_t**) calloc(size, sizeof(vertex_t*));

    *weights = (double **) malloc(sizeof(double *) * size);

    for (size_t i = 0; i < size; ++i) {
        (*weights)[i] = (double *) calloc(size, sizeof(double));
    }

    while (!q.empty()) {
        vertex_t *vertex = q.front();
        q.pop();

        DEBUG_PRINT("Matrix assigning %zu to %p\n", vertex->vertex_index, (void*)vertex);
        (*vertices)[vertex->vertex_index] = vertex;

        for (size_t i = 0; i < vertex->nedges; ++i) {
            llc_t child = vertex->edges[i];

            (*weights)[vertex->vertex_index][child.child->vertex_index] = child.weight;
            (*weights)[vertex->vertex_index][vertex->vertex_index] -= child.weight;
        }
    }

    return 0;
}

/*
 * Assumes the graph has been labelled, as they should all be.
 */
int graph_as_t_mat(double ***weights, size_t *out_size, vertex_t ***vertices, const vertex_t *graph) {
    queue<vertex_t*> q = enqueue_vertices((vertex_t*)graph);

    size_t size = (size_t)q.size();
    *out_size = size;

    *vertices = (vertex_t**) calloc(size, sizeof(vertex_t*));

    *weights = (double **) malloc(sizeof(double *) * size);

    for (size_t i = 0; i < size; ++i) {
        (*weights)[i] = (double *) calloc(size, sizeof(double));
    }

    while (!q.empty()) {
        vertex_t *vertex = q.front();
        q.pop();

        DEBUG_PRINT("Matrix assigning %zu to %p\n", vertex->vertex_index, (void*)vertex);
        (*vertices)[vertex->vertex_index] = vertex;

        for (size_t i = 0; i < vertex->nedges; ++i) {
            llc_t child = vertex->edges[i];

            (*weights)[vertex->vertex_index][child.child->vertex_index] = child.weight / vertex->rate;
        }

        (*weights)[vertex->vertex_index][vertex->vertex_index] = -1;
    }

    return 0;
}
