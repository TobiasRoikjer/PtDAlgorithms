/*
 * MIT License
 *
 * Copyright (c) 2021 Tobias Røikjer
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef PTDALGORITHMS_PTDCPP_H
#define PTDALGORITHMS_PTDCPP_H

#include <cstring>
#include <errno.h>
#include <vector>
#include <stdexcept>
#include <cmath>
#include "../c/ptdalgorithms.h"

namespace ptdalgorithms {
    struct rf_graph {
        struct ptd_avl_tree *tree;
        struct ptd_graph *graph;
        size_t *references;
        struct ptd_dph_probability_distribution_context *dph_context;
        struct ptd_probability_distribution_context *ph_context;
        struct ptd_dph_probability_distribution_context *dph_context_markov;
        struct ptd_probability_distribution_context *ph_context_markov;
        int granularity;
        int granularity_markov;
    };

    class Vertex;

    struct Edge;

    struct ParameterizedEdge;

    class PhaseTypeDistribution;

    class Graph;

    class SccGraph;

    class Graph {
    public:
        Graph(struct ptd_graph *graph) {
            this->rf_graph = (struct rf_graph *) malloc(sizeof(*this->rf_graph));
            this->rf_graph->references = (size_t *) malloc(sizeof(*this->rf_graph->references));
            *this->rf_graph->references = 1;
            this->rf_graph->graph = graph;
            this->rf_graph->tree = ptd_avl_tree_create(this->rf_graph->graph->state_length);
            this->rf_graph->dph_context = NULL;
            this->rf_graph->ph_context = NULL;
            this->rf_graph->dph_context_markov = NULL;
            this->rf_graph->ph_context_markov = NULL;

            if (this->rf_graph->tree == NULL) {
                throw std::runtime_error("Failed to create ptd_avl_tree\n");
            }
        }

        Graph(struct ptd_graph *graph, struct ptd_avl_tree *avl_tree) {
            this->rf_graph = (struct rf_graph *) malloc(sizeof(*this->rf_graph));
            this->rf_graph->references = (size_t *) malloc(sizeof(*this->rf_graph->references));
            *this->rf_graph->references = 1;
            this->rf_graph->graph = graph;
            this->rf_graph->tree = avl_tree;
            this->rf_graph->dph_context = NULL;
            this->rf_graph->ph_context = NULL;
            this->rf_graph->dph_context_markov = NULL;
            this->rf_graph->ph_context_markov = NULL;
        }

        Graph(const Graph &o) {
            this->rf_graph = (struct rf_graph *) malloc(sizeof(*this->rf_graph));
            this->rf_graph->references = o.rf_graph->references;
            this->rf_graph->graph = o.rf_graph->graph;
            this->rf_graph->tree = o.rf_graph->tree;
            this->rf_graph->dph_context = o.rf_graph->dph_context;
            this->rf_graph->ph_context = o.rf_graph->ph_context;
            this->rf_graph->dph_context_markov = o.rf_graph->dph_context_markov;
            this->rf_graph->ph_context_markov = o.rf_graph->ph_context_markov;
            this->rf_graph->granularity = o.rf_graph->granularity;
            this->rf_graph->granularity_markov = o.rf_graph->granularity_markov;
            *(this->rf_graph->references) += 1;
        }

        Graph(size_t state_length) {
            this->rf_graph = (struct rf_graph *) malloc(sizeof(*this->rf_graph));
            this->rf_graph->references = (size_t *) malloc(sizeof(*this->rf_graph->references));
            *this->rf_graph->references = 1;
            this->rf_graph->graph = ptd_graph_create(state_length);

            if (this->rf_graph->graph == NULL) {
                throw std::runtime_error("Failed to create ptd_graph\n");
            }

            this->rf_graph->tree = ptd_avl_tree_create(this->rf_graph->graph->state_length);

            if (this->rf_graph->tree == NULL) {
                throw std::runtime_error("Failed to create ptd_avl_tree\n");
            }

            this->rf_graph->dph_context = NULL;
            this->rf_graph->ph_context = NULL;
            this->rf_graph->dph_context_markov = NULL;
            this->rf_graph->ph_context_markov = NULL;
        }

        ~Graph() {
            *(this->rf_graph->references) -= 1;

            if (*this->rf_graph->references == 0) {
                ptd_dph_probability_distribution_context_destroy(this->rf_graph->dph_context);
                ptd_probability_distribution_context_destroy(this->rf_graph->ph_context);
                ptd_dph_probability_distribution_context_destroy(this->rf_graph->dph_context_markov);
                ptd_probability_distribution_context_destroy(this->rf_graph->ph_context_markov);
                ptd_avl_tree_destroy(this->rf_graph->tree);
                ptd_graph_destroy(this->rf_graph->graph);
                this->rf_graph->dph_context = NULL;
                this->rf_graph->ph_context = NULL;
                this->rf_graph->dph_context_markov = NULL;
                this->rf_graph->ph_context_markov = NULL;
                free(this->rf_graph->references);
            }

            free(this->rf_graph);
        }

        void update_weights_parameterized(std::vector<double> scalars);

        std::vector<double> expected_waiting_time(std::vector<double> rewards = std::vector<double>()) {
            double *ptr = ptd_expected_waiting_time(
                    this->c_graph(),
                    rewards.empty() ? NULL : &rewards[0]
            );

            if (ptr == NULL) {
                throw std::runtime_error((char *) ptd_err);
            }

            std::vector<double> res;
            res.assign(ptr, ptr + this->c_graph()->vertices_length);
            free(ptr);

            return res;
        }

        Vertex create_vertex(std::vector<int> state = std::vector<int>());

        Vertex create_vertex(const int *state);

        Vertex *create_vertex_p(std::vector<int> state = std::vector<int>());

        Vertex *create_vertex_p(const int *state);

        Vertex find_vertex(std::vector<int> state);

        Vertex find_vertex(const int *state);

        Vertex *find_vertex_p(std::vector<int> state);

        Vertex *find_vertex_p(const int *state);

        bool vertex_exists(std::vector<int> state);

        bool vertex_exists(const int *state);

        Vertex find_or_create_vertex(std::vector<int> state);

        Vertex find_or_create_vertex(const int *state);

        Vertex *find_or_create_vertex_p(std::vector<int> state);

        Vertex *find_or_create_vertex_p(const int *state);

        Vertex starting_vertex();

        Vertex *starting_vertex_p();

        std::vector<Vertex> vertices();

        std::vector<Vertex *> vertices_p();

        Vertex vertex_at(size_t index);

        Vertex *vertex_at_p(size_t index);

        size_t vertices_length();

        long double random_sample(std::vector<double> rewards = std::vector<double>()) {
            return ptd_random_sample(c_graph(), &rewards[0]);
        }

        std::vector<long double> mph_random_sample(std::vector<double> rewards, size_t vertex_rewards_length) {
            std::vector<long double> res(vertex_rewards_length);
            long double *c_res = ptd_mph_random_sample(c_graph(), &rewards[0], vertex_rewards_length);

            for (size_t i = 0; i < vertex_rewards_length; ++i) {
                res[i] = c_res[i];
            }

            free(c_res);

            return res;
        }

        long double dph_random_sample_c(double *rewards) {
            long double res = ptd_dph_random_sample(c_graph(), rewards);

            if (std::isnan(res)) {
                throw std::runtime_error((char *) ptd_err);
            }

            return res;
        }

        long double dph_random_sample(std::vector<double> rewards = std::vector<double>()) {
            if (rewards.empty()) {
                return this->dph_random_sample_c(NULL);
            } else {
                return this->dph_random_sample_c(&rewards[0]);
            }
        }

        long double *mdph_random_sample_c(double *rewards, size_t vertex_rewards_length) {
            long double *res = ptd_mdph_random_sample(c_graph(), rewards, vertex_rewards_length);

            if (res == NULL) {
                throw std::runtime_error((char *) ptd_err);
            }

            return res;
        }

        std::vector<long double> mdph_random_sample(std::vector<double> rewards, size_t vertex_rewards_length) {
            std::vector<long double> res(vertex_rewards_length);
            long double *c_res;

            if (rewards.empty()) {
                c_res = mdph_random_sample_c(NULL, vertex_rewards_length);
            } else {
                c_res = mdph_random_sample_c(&rewards[0], vertex_rewards_length);
            }

            for (size_t i = 0; i < vertex_rewards_length; ++i) {
                res[i] = c_res[i];
            }

            free(c_res);

            return res;
        }

        std::vector<long double> mdph_random_sample_c(std::vector<double> rewards, size_t vertex_rewards_length) {
            std::vector<long double> res(vertex_rewards_length);
            long double *c_res = ptd_mdph_random_sample(c_graph(), &rewards[0], vertex_rewards_length);

            for (size_t i = 0; i < vertex_rewards_length; ++i) {
                res[i] = c_res[i];
            }

            free(c_res);

            return res;
        }

        size_t random_sample_stop_vertex(double time) {
            struct ptd_vertex *res = ptd_random_sample_stop_vertex(c_graph(), time);

            return res->index;
        }

        size_t dph_random_sample_stop_vertex(int jumps) {
            struct ptd_vertex *res = ptd_dph_random_sample_stop_vertex(c_graph(), jumps);

            if (res == NULL) {
                throw std::runtime_error((char *) ptd_err);
            }

            return res->index;
        }

        size_t state_length() {
            return c_graph()->state_length;
        }

        PhaseTypeDistribution phase_type_distribution();

        bool is_acyclic() {
            return ptd_graph_is_acyclic(c_graph());
        }

        void validate() {
            if (ptd_validate_graph(c_graph())) {
                throw std::runtime_error((char *) ptd_err);
            }
        }

        Graph reward_transform(std::vector<double> rewards);

        Graph *reward_transform_p(std::vector<double> rewards);

        Graph dph_reward_transform(std::vector<int> rewards) {
            struct ptd_graph *res = ptd_graph_dph_reward_transform(c_graph(), &rewards[0]);

            if (res == NULL) {
                throw std::runtime_error((char *) ptd_err);
            }

            return Graph(res);
        }

        Graph *dph_reward_transform_p(std::vector<int> rewards) {
            struct ptd_graph *res = ptd_graph_dph_reward_transform(c_graph(), &rewards[0]);

            if (res == NULL) {
                throw std::runtime_error((char *) ptd_err);
            }

            return new Graph(res);
        }

        std::vector<double> normalize() {
            double *rewards = ptd_normalize_graph(this->c_graph());

            std::vector<double> res(this->c_graph()->vertices_length);

            for (size_t i = 0; i < res.size(); ++i) {
                res[i] = rewards[i];
            }

            free(rewards);

            notify_change();

            return res;
        }

        std::vector<double> dph_normalize() {
            double *rewards = ptd_dph_normalize_graph(this->c_graph());

            std::vector<double> res(this->c_graph()->vertices_length);

            for (size_t i = 0; i < res.size(); ++i) {
                res[i] = rewards[i];
            }

            free(rewards);

            notify_change();

            return res;
        }

        void notify_change() {
            ptd_notify_change(c_graph());
            ptd_probability_distribution_context_destroy(this->rf_graph->ph_context);
            this->rf_graph->ph_context = NULL;
            ptd_probability_distribution_context_destroy(this->rf_graph->ph_context_markov);
            this->rf_graph->ph_context_markov = NULL;
            ptd_dph_probability_distribution_context_destroy(this->rf_graph->dph_context);
            this->rf_graph->dph_context = NULL;
            ptd_dph_probability_distribution_context_destroy(this->rf_graph->dph_context_markov);
            this->rf_graph->dph_context_markov = NULL;
        }

        double defect() {
            return ptd_defect(c_graph());
        }

        Graph clone() {
            struct ptd_clone_res r = ptd_clone_graph(c_graph(), c_avl_tree());

            return Graph(r.graph, r.avl_tree);
        }


        Graph *clone_p() {
            struct ptd_clone_res r = ptd_clone_graph(c_graph(), c_avl_tree());

            return new Graph(r.graph, r.avl_tree);
        }

        double pdf(float time, int granularity = 0) {
            if (this->rf_graph->ph_context == NULL || this->rf_graph->granularity != granularity) {
                if (this->rf_graph->ph_context != NULL) {
                    ptd_probability_distribution_context_destroy(this->rf_graph->ph_context);
                }
                this->rf_graph->ph_context = ptd_probability_distribution_context_create(c_graph(), granularity);

                if (this->rf_graph->ph_context == NULL) {
                    throw std::runtime_error((char *) ptd_err);
                }

                _pdf.clear();
                _cdf.clear();
                _pdf.push_back(this->rf_graph->ph_context->pdf);
                _cdf.push_back(this->rf_graph->ph_context->cdf);
                this->rf_graph->granularity = granularity;
            }

            while (time >= this->rf_graph->ph_context->time) {
                ptd_probability_distribution_step(
                        this->rf_graph->ph_context
                );
                _pdf.push_back(this->rf_graph->ph_context->pdf);
                _cdf.push_back(this->rf_graph->ph_context->cdf);
            }

            return _pdf[this->rf_graph->ph_context->granularity * time];
        }

        double cdf(float time, int granularity = 0) {
            pdf(time, granularity);

            return _cdf[this->rf_graph->ph_context->granularity * time];
        }

        double dph_pmf(int jumps) {
            if (this->rf_graph->dph_context == NULL) {
                this->rf_graph->dph_context = ptd_dph_probability_distribution_context_create(c_graph());

                if (this->rf_graph->dph_context == NULL) {
                    throw std::runtime_error((char *) ptd_err);
                }
                _dph_pmf.clear();
                _dph_cdf.clear();
                _dph_pmf.push_back(this->rf_graph->dph_context->pmf);
                _dph_cdf.push_back(this->rf_graph->dph_context->cdf);
            }

            if (jumps > this->rf_graph->dph_context->jumps) {
                for (int i = this->rf_graph->dph_context->jumps; i < jumps; ++i) {
                    ptd_dph_probability_distribution_step(
                            this->rf_graph->dph_context
                    );
                    _dph_pmf.push_back(this->rf_graph->dph_context->pmf);
                    _dph_cdf.push_back(this->rf_graph->dph_context->cdf);
                }
            }

            return _dph_pmf[jumps];
        }

        double dph_cdf(int jumps) {
            dph_pmf(jumps);

            return _dph_cdf[jumps];
        }

        std::vector<double> stop_probability(double time, int granularity = 0) {
            if (this->rf_graph->ph_context_markov == NULL
                || this->rf_graph->granularity_markov != granularity
                || this->rf_graph->ph_context_markov->time -
                   ((double) 1.0) / this->rf_graph->ph_context_markov->granularity >
                   time) {
                if (this->rf_graph->ph_context_markov != NULL) {
                    ptd_probability_distribution_context_destroy(this->rf_graph->ph_context_markov);
                }
                this->rf_graph->ph_context_markov = ptd_probability_distribution_context_create(c_graph(), granularity);

                if (this->rf_graph->ph_context_markov == NULL) {
                    throw std::runtime_error((char *) ptd_err);
                }

                this->rf_graph->granularity_markov = granularity;
            }

            while (time > this->rf_graph->ph_context_markov->time) {
                ptd_probability_distribution_step(
                        this->rf_graph->ph_context_markov
                );
            }

            std::vector<double> ret(this->rf_graph->ph_context_markov->graph->vertices_length);

            for (size_t i = 0; i < this->rf_graph->ph_context_markov->graph->vertices_length; ++i) {
                ret[i] = (double) this->rf_graph->ph_context_markov->probability_at[i];
            }

            return ret;
        }

        std::vector<double> accumulated_visiting_time(double time, int granularity = 0) {
            if (this->rf_graph->ph_context_markov == NULL
                || this->rf_graph->granularity_markov != granularity
                || this->rf_graph->ph_context_markov->time -
                   ((double) 1.0) / this->rf_graph->ph_context_markov->granularity >
                   time) {
                if (this->rf_graph->ph_context_markov != NULL) {
                    ptd_probability_distribution_context_destroy(this->rf_graph->ph_context_markov);
                }
                this->rf_graph->ph_context_markov = ptd_probability_distribution_context_create(c_graph(), granularity);

                if (this->rf_graph->ph_context_markov == NULL) {
                    throw std::runtime_error((char *) ptd_err);
                }

                this->rf_graph->granularity_markov = granularity;
            }

            while (time > this->rf_graph->ph_context_markov->time) {
                ptd_probability_distribution_step(
                        this->rf_graph->ph_context_markov
                );
            }

            std::vector<double> ret(this->rf_graph->ph_context_markov->graph->vertices_length);

            for (size_t i = 0; i < this->rf_graph->ph_context_markov->graph->vertices_length; ++i) {
                ret[i] = (double) this->rf_graph->ph_context_markov->accumulated_visits[i];
                ret[i] /= this->rf_graph->ph_context_markov->granularity;
            }

            return ret;
        }

        std::vector<double> dph_stop_probability(int jumps) {
            if (this->rf_graph->dph_context_markov == NULL
                || this->rf_graph->dph_context_markov->jumps > jumps) {
                if (this->rf_graph->dph_context_markov != NULL) {
                    ptd_dph_probability_distribution_context_destroy(this->rf_graph->dph_context_markov);
                }
                this->rf_graph->dph_context_markov = ptd_dph_probability_distribution_context_create(c_graph());

                if (this->rf_graph->dph_context_markov == NULL) {
                    throw std::runtime_error((char *) ptd_err);
                }
            }

            while (jumps > this->rf_graph->dph_context_markov->jumps) {
                ptd_dph_probability_distribution_step(
                        this->rf_graph->dph_context_markov
                );
            }

            std::vector<double> ret(this->rf_graph->dph_context_markov->graph->vertices_length);

            for (size_t i = 0; i < this->rf_graph->dph_context_markov->graph->vertices_length; ++i) {
                ret[i] = (double) this->rf_graph->dph_context_markov->probability_at[i];
            }

            return ret;
        }

        std::vector<double> dph_accumulated_visits(int jumps) {
            if (this->rf_graph->dph_context_markov == NULL
                || this->rf_graph->dph_context_markov->jumps > jumps) {
                if (this->rf_graph->dph_context_markov != NULL) {
                    ptd_dph_probability_distribution_context_destroy(this->rf_graph->dph_context_markov);
                }
                this->rf_graph->dph_context_markov = ptd_dph_probability_distribution_context_create(c_graph());

                if (this->rf_graph->dph_context_markov == NULL) {
                    throw std::runtime_error((char *) ptd_err);
                }
            }

            while (jumps > this->rf_graph->dph_context_markov->jumps) {
                ptd_dph_probability_distribution_step(
                        this->rf_graph->dph_context_markov
                );
            }

            std::vector<double> ret(this->rf_graph->dph_context_markov->graph->vertices_length);

            for (size_t i = 0; i < this->rf_graph->dph_context_markov->graph->vertices_length; ++i) {
                ret[i] = (double) this->rf_graph->dph_context_markov->accumulated_visits[i];
            }

            return ret;
        }

        std::vector<double> dph_expected_visits(int jumps) {
            if (this->rf_graph->dph_context_markov == NULL
                || this->rf_graph->dph_context_markov->jumps > jumps) {
                if (this->rf_graph->dph_context_markov != NULL) {
                    ptd_dph_probability_distribution_context_destroy(this->rf_graph->dph_context_markov);
                }
                this->rf_graph->dph_context_markov = ptd_dph_probability_distribution_context_create(c_graph());

                if (this->rf_graph->dph_context_markov == NULL) {
                    throw std::runtime_error((char *) ptd_err);
                }
            }

            while (jumps > this->rf_graph->dph_context_markov->jumps) {
                ptd_dph_probability_distribution_step(
                        this->rf_graph->dph_context_markov
                );
            }

            std::vector<double> ret(this->rf_graph->dph_context_markov->graph->vertices_length);

            for (size_t i = 0; i < this->rf_graph->dph_context_markov->graph->vertices_length; ++i) {
                ret[i] = (double) this->rf_graph->dph_context_markov->accumulated_visits[i];
            }

            return ret;
        }

    public:
        Graph &operator=(const Graph &o) {
            if (this == &o) {
                return *this;
            }

            (*this).rf_graph = (struct rf_graph *) malloc(sizeof(*this->rf_graph));

            *this->rf_graph->references -= 1;

            if (*this->rf_graph->references == 0) {
                ptd_avl_tree_destroy(this->rf_graph->tree);
                ptd_graph_destroy(this->rf_graph->graph);
                free(this->rf_graph->references);
                ptd_probability_distribution_context_destroy(this->rf_graph->ph_context);
                ptd_dph_probability_distribution_context_destroy(this->rf_graph->dph_context);
                this->rf_graph->dph_context = NULL;
            }

            free(this->rf_graph);

            //this->rf_graph = o.rf_graph;
            this->rf_graph = (struct rf_graph *) malloc(sizeof(*this->rf_graph));
            this->rf_graph->references = o.rf_graph->references;
            this->rf_graph->graph = o.rf_graph->graph;
            this->rf_graph->tree = o.rf_graph->tree;
            this->rf_graph->ph_context = o.rf_graph->ph_context;
            this->rf_graph->dph_context = o.rf_graph->dph_context;
            *(this->rf_graph->references) += 1;

            return *this;
        }

        struct ptd_graph *c_graph() {
            return rf_graph->graph;
        }

        struct ptd_avl_tree *c_avl_tree() {
            return rf_graph->tree;
        }

    private:
        Graph(struct rf_graph *rf_graph) {
            this->rf_graph = rf_graph;
            rf_graph->references++;
        }

        struct rf_graph *rf_graph;

        std::vector<double> _pdf;
        std::vector<double> _cdf;
        std::vector<double> _dph_pmf;
        std::vector<double> _dph_cdf;

        friend class VertexLinkedList;

        friend class Vertex;
    };

    class Vertex {
    private:
        Vertex(Graph &graph, int *state) : graph(graph) {
            this->vertex = ptd_vertex_create_state(graph.rf_graph->graph, state);

            if (this->vertex == NULL) {
                throw std::runtime_error("Failed to create ptd_vertex\n");
            }
        }

    public:
        Vertex(Graph &graph, struct ptd_vertex *vertex) : graph(graph) {
            this->vertex = vertex;
        }

        Vertex(const Vertex &o) : graph(o.graph) {
            this->vertex = o.vertex;
        }

        ~Vertex() {
        }

        void add_edge(Vertex &to, double weight);

        void add_edge_parameterized(Vertex &to, double weight, std::vector<double> edge_state);

        std::vector<int> state();

        std::vector<Edge> edges();

        std::vector<ParameterizedEdge> parameterized_edges();

        bool operator==(const Vertex &other) const {
            return vertex == other.vertex;
        }

        Vertex &operator=(const Vertex &o) {
            vertex = o.vertex;
            graph = Graph(o.graph.rf_graph);

            return *this;
        }

        struct ptd_vertex *c_vertex() {
            return vertex;
        }

        double rate() {
            return ptd_vertex_rate(vertex);
        }

    private:
        Graph &graph;

        struct ptd_vertex *vertex;

        friend class Graph;
    };

    struct Edge {
    private:
        Edge(struct ptd_vertex *vertex, struct ptd_edge *edge, Graph &graph, double weight) : graph(graph) {
            this->_weight = weight;
            this->_edge = edge;
            this->_vertex = vertex;
        }

    private:
        Graph &graph;
        struct ptd_vertex *_vertex;
        struct ptd_edge *_edge;
        double _weight;

    public:
        Vertex to() {
            return Vertex(graph, _vertex);
        }

	Vertex *to_p() {
            return new Vertex(graph, _vertex);
        }

        double weight() {
            return _weight;
        }

        void update_weight(double weight) {
            ptd_edge_update_weight(_edge, weight);
            _weight = weight;
        }

        Edge &operator=(const Edge &o) {
            _weight = o._weight;
            _vertex = o._vertex;
            graph = o.graph;

            return *this;
        }

        friend class Vertex;

        friend class ParameterizedEdge;
    };

    struct ParameterizedEdge : private Edge {
    private:
        ParameterizedEdge(
                struct ptd_vertex *vertex,
                struct ptd_edge *edge,
                Graph &graph,
                double weight,
                double *state
        ) : Edge(vertex, edge, graph, weight) {
            this->_state = state;
        }

    private:
        double *_state;
        size_t state_size;

    public:
        Vertex to() {
            return Vertex(graph, _vertex);
        }

	Vertex *to_p() {
            return new Vertex(graph, _vertex);
        }

        double weight() {
            return _weight;
        }

        std::vector<double> edge_state(size_t state_length) {
            std::vector<double> state;

            if (_state != NULL) {
                for (size_t i = 0; i < state_length; ++i) {
                    state.push_back(_state[i]);
                }
            }

            return state;
        }

        ParameterizedEdge &operator=(const ParameterizedEdge &o) {
            _weight = o._weight;
            _vertex = o._vertex;
            graph = o.graph;
            _state = o._state;

            return *this;
        }

        friend class Vertex;
    };

    class PhaseTypeDistribution {
    private:
        PhaseTypeDistribution(Graph &graph, struct ptd_phase_type_distribution *matrix) {
            this->length = matrix->length;
            this->sub_intensity_matrix = matrix->sub_intensity_matrix;
            this->initial_probability_vector = matrix->initial_probability_vector;

            for (size_t i = 0; i < matrix->length; ++i) {
                this->vertices.push_back(Vertex(graph, matrix->vertices[i]));
            }

            this->distribution = matrix;
        }

        struct ptd_phase_type_distribution *distribution;

    public:
        ~PhaseTypeDistribution() {
            ptd_phase_type_distribution_destroy(distribution);
        }

        size_t length;
        double **sub_intensity_matrix;
        double *initial_probability_vector;
        std::vector<Vertex> vertices;

        friend class Graph;
    };

    class AnyProbabilityDistributionContext {
    public:
        int is_discrete() {
            return _discrete;
        }

        virtual void step() {
            throw std::runtime_error("Not implemented");
        }

        virtual double pmf() {
            throw std::runtime_error("Not implemented");
        }

        virtual double pdf() {
            throw std::runtime_error("Not implemented");
        }

        virtual double cdf() {
            throw std::runtime_error("Not implemented");
        }

        virtual double time() {
            throw std::runtime_error("Not implemented");
        }

        virtual int jumps() {
            throw std::runtime_error("Not implemented");
        }

        virtual std::vector<long double> stop_probability() {
            throw std::runtime_error("Not implemented");
        }

        virtual std::vector<long double> accumulated_visits() {
            throw std::runtime_error("Not implemented");
        }

        virtual std::vector<long double> accumulated_visiting_time() {
            throw std::runtime_error("Not implemented");
        }

        virtual ~AnyProbabilityDistributionContext() {};

    protected:
        int _discrete;
    };

    class ProbabilityDistributionContext : private AnyProbabilityDistributionContext {
    public:
        ProbabilityDistributionContext(Graph &graph, int granularity = 0) : graph(graph) {
            context = ptd_probability_distribution_context_create(graph.c_graph(), granularity);
            _discrete = false;
        }

        void step() {
            ptd_probability_distribution_step(context);
        }

        double pdf() {
            return context->pdf;
        }

        double cdf() {
            return context->cdf;
        }

        double time() {
            return (double) context->time;
        }

        std::vector<long double> stop_probability() {
            return std::vector<long double>(context->probability_at,
                                            context->probability_at + context->graph->vertices_length);
        }

        std::vector<long double> accumulated_visiting_time() {
            std::vector<long double> res(context->accumulated_visits,
                                         context->accumulated_visits + context->graph->vertices_length);

            for (size_t i = 0; i < context->graph->vertices_length; ++i) {
                res[i] /= (long double) context->granularity;
            }

            return res;
        }

        ~ProbabilityDistributionContext() {
            ptd_probability_distribution_context_destroy(context);
        }

    private:
        Graph &graph;
        struct ptd_probability_distribution_context *context;
    };

    class DPHProbabilityDistributionContext : private AnyProbabilityDistributionContext {
    public:
        DPHProbabilityDistributionContext(Graph &graph) : graph(graph) {
            context = ptd_dph_probability_distribution_context_create(graph.c_graph());
            _discrete = true;
        }

        void step() {
            ptd_dph_probability_distribution_step(context);
        }

        double pmf() {
            return context->pmf;
        }

        double cdf() {
            return context->cdf;
        }

        int jumps() {
            return context->jumps;
        }

        std::vector<long double> stop_probability() {
            return std::vector<long double>(context->probability_at,
                                            context->probability_at + context->graph->vertices_length);
        }


        std::vector<long double> accumulated_visits() {
            return std::vector<long double>(context->accumulated_visits,
                                            context->accumulated_visits + context->graph->vertices_length);
        }

        ~DPHProbabilityDistributionContext() {
            ptd_dph_probability_distribution_context_destroy(context);
        }

    private:
        Graph &graph;
        struct ptd_dph_probability_distribution_context *context;
    };
}

#endif //PTDALGORITHMS_PTDCPP_H
