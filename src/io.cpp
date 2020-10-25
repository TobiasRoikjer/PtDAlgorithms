#include "io.h"
#include "phase.h"

int graph_as_mat(double ***weights, size_t *out_size, vertex_t ***vertices, const vertex_t *graph) {
    queue<vertex_t*> queue = enqueue_vertices((vertex_t*)graph);
    label_vertex_index(NULL, (vertex_t*)graph);

    size_t size = (size_t)queue.size();
    *out_size = size;

    *vertices = (vertex_t**) calloc(size, sizeof(vertex_t*));

    *weights = (double **) malloc(sizeof(double *) * size);

    for (size_t i = 0; i < size; ++i) {
        (*weights)[i] = (double *) calloc(size, sizeof(double));
    }

    while (!queue.empty()) {
        vertex_t *vertex = queue.front();
        queue.pop();

        (*vertices)[vertex->vertex_index] = vertex;

        for (size_t i = 0; i < vertex->nedges; ++i) {
            llc_t child = vertex->edges[i];

            (*weights)[vertex->vertex_index][child.child->vertex_index] = child.weight;
            (*weights)[vertex->vertex_index][vertex->vertex_index] -= child.weight;
        }
    }

    return 0;
}
