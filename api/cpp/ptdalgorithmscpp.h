#ifndef PTDALGORITHMS_PTDCPP_H
#define PTDALGORITHMS_PTDCPP_H

#include <cstring>
#include <errno.h>
#include "../c/ptdalgorithms.h"

struct rf_graph {
    ptd_avl_tree_t *tree;
    struct ptd_ph_graph *graph;
    size_t *references;
};

namespace ptdalgorithms {
    class Vertex;

    struct Edge;

    class PhaseTypeDistribution;

    class Graph;

    class Graph {
    public:
        Graph(struct ptd_ph_graph *graph) {
            this->rf_graph = (struct rf_graph *) malloc(sizeof(*this->rf_graph));
            this->rf_graph->references = (size_t *) malloc(sizeof(*this->rf_graph->references));
            *this->rf_graph->references = 1;
            this->rf_graph->graph = graph;
            this->rf_graph->tree = ptd_avl_tree_create(this->rf_graph->graph->state_length);

            if (this->rf_graph->tree == NULL) {
                throw std::runtime_error("Failed to create ptd_avl_tree\n");
            }
        }

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
            this->rf_graph->graph = ptd_ph_graph_create(state_length);

            if (this->rf_graph->graph == NULL) {
                throw std::runtime_error("Failed to create ptd_graph\n");
            }

            this->rf_graph->tree = ptd_avl_tree_create(this->rf_graph->graph->state_length);

            if (this->rf_graph->tree == NULL) {
                throw std::runtime_error("Failed to create ptd_avl_tree\n");
            }

        }

        ~Graph() {
            *(this->rf_graph->references) -= 1;


            if (*this->rf_graph->references == 0) {
                fprintf(stderr, "Destroying graph %i \n", (int) (*(this->rf_graph->references)));
                ptd_avl_tree_vertex_destroy(this->rf_graph->tree);
                ptd_ph_graph_destroy(this->rf_graph->graph);
                free(this->rf_graph->references);
            }

            free(this->rf_graph);
        }

        Vertex create_vertex(std::vector<int> state);

        Vertex *create_vertex_p(std::vector<int> state);

        Vertex find_vertex(std::vector<int> state);

        Vertex *find_vertex_p(std::vector<int> state);

        bool vertex_exists(std::vector<int> state);

        Vertex find_or_create_vertex(std::vector<int> state);

        Vertex *find_or_create_vertex_p(std::vector<int> state);

        Vertex starting_vertex();

        Vertex *starting_vertex_p();

        std::vector<Vertex> vertices();

        size_t state_length() {
            return c_graph()->state_length;
        }

        void visit_vertices(int (*visit_func)(Graph &graph, Vertex &vertex),
                            bool include_start = false);


        PhaseTypeDistribution phase_type_distribution();

    public:
        Graph &operator=(const Graph &o) {
            if (this == &o) {
                return *this;
            }

            (*this).rf_graph = (struct rf_graph *) malloc(sizeof(*this->rf_graph));

            *this->rf_graph->references -= 1;

            if (*this->rf_graph->references == 0) {
                ptd_avl_tree_vertex_destroy_free(this->rf_graph->tree);
                ptd_ph_graph_destroy(this->rf_graph->graph);
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

        struct ptd_ph_graph *c_graph() {
            return rf_graph->graph;
        }

    private:
        Graph(struct rf_graph *rf_graph) {
            this->rf_graph = rf_graph;
            rf_graph->references++;
        }

        struct rf_graph *rf_graph;

        friend class VertexLinkedList;

        friend class Vertex;
    };

    class Vertex {
    private:
        Vertex(Graph &graph, int *state) : graph(graph) {
            this->vertex = ptd_ph_vertex_create_state(graph.rf_graph->graph, state);

            if (this->vertex == NULL) {
                throw std::runtime_error("Failed to create ptd_ph_vertex\n");
            }
        }

    public:
        Vertex(Graph &graph, struct ptd_ph_vertex *vertex) : graph(graph) {
            this->vertex = vertex;
        }

        Vertex(const Vertex &o) : graph(o.graph) {
            this->vertex = o.vertex;
        }

        ~Vertex() {
        }

        void add_edge(Vertex &to, double weight);

        std::vector<int> state();

        std::vector<Edge> edges();

        bool operator==(const Vertex &other) const {
            return vertex == other.vertex;
        }

        Vertex &operator=(const Vertex &o) {
            vertex = o.vertex;
            graph = Graph(o.graph.rf_graph);

            return *this;
        }

        struct ptd_ph_vertex *c_vertex() {
            return vertex;
        }

        double rate() {
            return ptd_ph_vertex_rate(vertex);
        }

    private:
        Graph &graph;

        struct ptd_ph_vertex *vertex;

        friend class Graph;
    };

    struct Edge {
    private:
        Edge(struct ptd_ph_vertex *vertex, Graph &graph, double weight) : graph(graph) {
            this->_weight = weight;
            this->_vertex = vertex;
        }

    private:
        Graph &graph;
        struct ptd_ph_vertex *_vertex;
        double _weight;

    public:
        Vertex to() {
            return Vertex(graph, _vertex);
        }

        double weight() {
            return _weight;
        }

        Edge &operator=(const Edge &o) {
            _weight = o._weight;
            _vertex = o._vertex;
            graph = o.graph;

            return *this;
        }

        friend class Vertex;
    };

    class PhaseTypeDistribution {
    private:
        PhaseTypeDistribution(Graph &graph, ptd_phase_type_distribution_t *matrix) {
            this->length = matrix->length;
            this->sub_intensity_matrix = matrix->sub_intensity_matrix;
            this->initial_probability_vector = matrix->initial_probability_vector;

            for (size_t i = 0; i < matrix->length; ++i) {
                this->vertices.push_back(Vertex(graph, matrix->vertices[i]));
            }

            this->distribution = matrix;
        }

        ptd_phase_type_distribution_t *distribution;

    public:
        ~PhaseTypeDistribution() {
            ptd_phase_type_distribution_destroy(distribution);
        }

        size_t length;
        double **sub_intensity_matrix;
        double *initial_probability_vector;
        vector<Vertex> vertices;

        friend class Graph;
    };

    /*class Models {
    public:
        static Graph kingman(size_t n) {
            char msg[1024];

            struct ptd_graph *kingman;

            if ((kingman = ptd_model_kingman(n)) == NULL) {
                snprintf(
                        msg, 1024,
                        "Failed to create Kingman graph: %s \n", std::strerror(errno)
                );

                throw new std::runtime_error(
                        msg
                );
            }

            return Graph(kingman);
        }
    };*/
}

#endif //PTDALGORITHMS_PTDCPP_H