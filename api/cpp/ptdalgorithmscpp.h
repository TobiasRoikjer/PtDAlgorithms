#ifndef PTDALGORITHMS_PTDCPP_H
#define PTDALGORITHMS_PTDCPP_H

#include "../c/ptdalgorithms.h"

namespace ptdalgorithms {
    class Vertex;
    class Edge;

    class Graph {
    public:
        Graph(size_t state_length) {
            fprintf(stderr, "Creating graph\n");
            this->graph = ptd_graph_create(state_length);

            if (this->graph == NULL) {
                throw std::runtime_error("Failed to create ptd_graph\n");
            }

            this->tree = ptd_avl_tree_create(graph->state_length);

            if (this->tree == NULL) {
                throw std::runtime_error("Failed to create ptd_avl_tree\n");
            }
        }

        ~Graph() {
            fprintf(stderr, "Destroying graph\n");
            ptd_avl_tree_vertex_destroy_free(this->tree);
            ptd_graph_destroy(this->graph);
        }

        Vertex create_vertex(std::vector<size_t> state);
        Vertex find_vertex(std::vector<size_t> state) const;
        bool vertex_exists(std::vector<size_t> state);
        Vertex find_or_create_vertex(std::vector<size_t> state);
        Vertex start_vertex() const;
        std::vector<Vertex> vertices() const;

    public:
        ptd_graph_t *graph;
        ptd_avl_tree_t *tree;
    };

    class Vertex {
    private:
        Vertex(const Graph &graph, size_t *state) : graph(graph) {
            this->vertex = ptd_vertex_create_state(graph.graph, state);

            if (this->vertex == NULL) {
                throw std::runtime_error("Failed to create ptd_vertex\n");
            }
        }

        Vertex(const Graph &graph, ptd_vertex_t *vertex) : graph(graph) {
            this->vertex = vertex;
        }

    public:
        ~Vertex() {
            fprintf(stderr, "Destroying vertex\n");
            //ptd_vertex_destroy(this->vertex);
        }

        void add_edge(const Vertex &to, long double weight);
        std::vector<size_t> state();
        std::vector<Edge> edges();

        bool operator==(const Vertex& other) const {
            return vertex == other.vertex;
        }

    private:
        ptd_vertex_t *vertex;
        const Graph& graph;

        friend class Graph;
    };

    class Edge {
    public:
        Edge(const Vertex &to, long double weight) : to(to) {
            this->weight = weight;
        }

        const Vertex& to;
        long double weight;
    };
}

#endif //PTDALGORITHMS_PTDCPP_H