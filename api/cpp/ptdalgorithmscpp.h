#ifndef PTDALGORITHMS_PTDCPP_H
#define PTDALGORITHMS_PTDCPP_H

#include "../c/ptdalgorithms.h"

namespace ptdalgorithms {
    class Graph {
    public:
        Graph(size_t state_length, size_t reward_length) {
            this->graph = ptd_graph_create(state_length, reward_length);
            this->tree = ptd_avl_tree_create(graph->state_length);
        }

        ~Graph() {
            fprintf(stderr, "Destroying graph\n");
            ptd_avl_tree_vertex_destroy(this->tree);
            ptd_graph_destroy(this->graph);
        }

        void create_vertex();

    public:
        ptd_graph_t *graph;
        ptd_avl_tree_t *tree;
    };

    class Vertex {
    public:
        Vertex(const Graph &graph) : graph(graph) {
            this->vertex = ptd_vertex_create(graph.graph);
        }

        ~Vertex() {
            fprintf(stderr, "Destroying vertex\n");
            //ptd_vertex_destroy(this->vertex);
        }

    private:

        ptd_vertex_t *vertex;
        const Graph& graph;
    };
}

#endif //PTDALGORITHMS_PTDCPP_H