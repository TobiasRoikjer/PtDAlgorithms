#ifndef PTDALGORITHMS_PTDCPP_H
#define PTDALGORITHMS_PTDCPP_H

#include "../c/ptdalgorithms.h"

struct rf_graph {
    ptd_avl_tree_t *tree;
    ptd_graph_t *graph;
    size_t *references;
};

namespace ptdalgorithms {
    class Vertex;

    struct Edge;

    class Graph {
    public:
        Graph(const Graph &o) {
            this->rf_graph = (struct rf_graph *) malloc(sizeof(*this->rf_graph));
            this->rf_graph->references = o.rf_graph->references;
            this->rf_graph->graph = o.rf_graph->graph;
            this->rf_graph->tree = o.rf_graph->tree;
            *(this->rf_graph->references) += 1;
        }

        Graph(size_t state_length) {
            this->rf_graph = (struct rf_graph *) malloc(sizeof(*this->rf_graph));
            this->rf_graph->references = (size_t *) malloc(sizeof(*this->rf_graph->references));
            *this->rf_graph->references = 1;
            fprintf(stderr, "Creating graph\n");
            this->rf_graph->graph = ptd_graph_create(state_length);

            if (this->rf_graph->graph == NULL) {
                throw std::runtime_error("Failed to create ptd_graph\n");
            }

            this->rf_graph->tree = ptd_avl_tree_create(this->rf_graph->graph->state_length);

            if (this->rf_graph->tree == NULL) {
                throw std::runtime_error("Failed to create ptd_avl_tree\n");
            }

        }

        ~Graph() {
            fprintf(stderr, "Destroying graph %zu \n", *(this->rf_graph->references));
            *(this->rf_graph->references) -= 1;


            if (*this->rf_graph->references == 0) {
                ptd_avl_tree_vertex_destroy_free(this->rf_graph->tree);
                ptd_graph_destroy(this->rf_graph->graph);
                free(this->rf_graph->references);
            }

            free(this->rf_graph);
        }

        Vertex create_vertex(std::vector<size_t> state);

        Vertex *create_vertex_p(std::vector<size_t> state);

        Vertex find_vertex(std::vector<size_t> state);

        Vertex *find_vertex_p(std::vector<size_t> state);

        bool vertex_exists(std::vector<size_t> state);

        Vertex find_or_create_vertex(std::vector<size_t> state);

        Vertex *find_or_create_vertex_p(std::vector<size_t> state);

        Vertex start_vertex();

        Vertex *start_vertex_p();

        std::vector<Vertex> vertices();

    public:
        Graph &operator=(const Graph &o) {
            *this->rf_graph->references -= 1;

            if (*this->rf_graph->references == 0) {
                ptd_avl_tree_vertex_destroy_free(this->rf_graph->tree);
                ptd_graph_destroy(this->rf_graph->graph);
                free(this->rf_graph->references);
            }

            free(this->rf_graph);

            //this->rf_graph = o.rf_graph;
            this->rf_graph = (struct rf_graph *) malloc(sizeof(*this->rf_graph));
            this->rf_graph->references = o.rf_graph->references;
            this->rf_graph->graph = o.rf_graph->graph;
            this->rf_graph->tree = o.rf_graph->tree;
            *(this->rf_graph->references) += 1;

            return *this;
        }

    private:
        struct rf_graph *rf_graph;

        friend class Vertex;
    };

    class Vertex {
    private:
        Vertex(Graph &graph, size_t *state) : graph(graph) {
            this->vertex = ptd_vertex_create_state(graph.rf_graph->graph, state);

            if (this->vertex == NULL) {
                throw std::runtime_error("Failed to create ptd_vertex\n");
            }

            fprintf(stderr, "Vertex created\n");
        }

        Vertex(Graph &graph, ptd_vertex_t *vertex) : graph(graph) {
            this->vertex = vertex;

            fprintf(stderr, "Vertex created\n");
        }

    public:
        ~Vertex() {
            fprintf(stderr, "Destroying vertex  \n");
        }

        void add_edge(Vertex &to, long double weight);

        std::vector<size_t> state();

        std::vector<Edge> edges();

        bool operator==(const Vertex &other) const {
            return vertex == other.vertex;
        }

        Vertex &operator=(const Vertex &o) {
            vertex = o.vertex;
            graph = o.graph;

            return *this;
        }

    private:
        ptd_vertex_t *vertex;
        Graph &graph;

        friend class Graph;
    };

    struct Edge {
    private:
        Edge(const Vertex &vertex, long double weight) : to(vertex) {
            this->weight = weight;
        }

    public:
        const Vertex &to;
        long double weight;

        Edge &operator=(const Edge &o) {
            weight = o.weight;
            return *this;
        }

        friend class Vertex;
    };
}

#endif //PTDALGORITHMS_PTDCPP_H