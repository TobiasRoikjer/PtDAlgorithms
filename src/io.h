#ifndef PTDALGORITHMS_IO_H
#define PTDALGORITHMS_IO_H

#include "stdlib.h"
#include "phase.h"

int graph_as_mat(double ***weights, size_t *out_size, vertex_t ***vertices, const vertex_t *graph);
int graph_as_t_mat(double ***weights, size_t *out_size, vertex_t ***vertices, const vertex_t *graph);


#endif //PTDALGORITHMS_IO_H
