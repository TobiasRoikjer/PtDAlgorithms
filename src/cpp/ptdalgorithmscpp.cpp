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

#include <stdexcept>
#include <sstream>
#include <stack>
#include <cerrno>
#include <cstring>
#include <vector>
#include "ptdalgorithmscpp.h"

/* While it seems very strange to have this in a C file, the R code
 * has very strange linking behavior, and we therefore sometimes include
 * the same C file...
 */

#ifndef PTDALGORITHMS_PTDCPP_C
#define PTDALGORITHMS_PTDCPP_C

static void assert_same_length(std::vector<int> state, struct ptd_graph *graph) {
    if (state.size() != graph->state_length) {
        std::stringstream message;
        message << "Vector `state` argument must have same size as graph state length. Was '";
        message << state.size() << "' expected '" << graph->state_length << "'" << std::endl;

        throw std::invalid_argument(message.str());
    }
}

static int *force_same_length(std::vector<int> state, struct ptd_graph *graph) {
    int *res = (int *) calloc(graph->state_length, sizeof(int));
  
    for (size_t i = 0; i < state.size(); ++i) {
        res[i] = state[i];
    }
  
    return res;
}

ptdalgorithms::Vertex ptdalgorithms::Graph::create_vertex(std::vector<int> state) {
    int *c = force_same_length(state, c_graph());

    Vertex res = create_vertex(c);

    return res;
}


ptdalgorithms::Vertex ptdalgorithms::Graph::create_vertex(const int *state) {
    struct ptd_vertex *c_vertex = ptd_vertex_create_state(c_graph(), (int*)state);
  
    Vertex vertex = ptdalgorithms::Vertex(*this, c_vertex);
    notify_change();

    return vertex;
}

ptdalgorithms::Vertex *ptdalgorithms::Graph::create_vertex_p(std::vector<int> state) {
    int *c = force_same_length(state, c_graph());

    Vertex *res = create_vertex_p(c);

    return res;
}


ptdalgorithms::Vertex *ptdalgorithms::Graph::create_vertex_p(const int *state) {
    struct ptd_vertex *c_vertex = ptd_vertex_create_state(c_graph(), (int*)state);

    Vertex *vertex = new ptdalgorithms::Vertex(*this, c_vertex);
    notify_change();
    
    return vertex;
}

ptdalgorithms::Vertex ptdalgorithms::Graph::find_vertex(std::vector<int> state) {
    int *c = force_same_length(state, c_graph());

    Vertex res = find_vertex(c);

    free(c);

    return res;
}


ptdalgorithms::Vertex ptdalgorithms::Graph::find_vertex(const int *state) {
    struct ptd_avl_node *node = ptd_avl_tree_find(this->rf_graph->tree, state);

    if (node == NULL) {
        throw std::runtime_error(
                "No such vertex found\n"
        );
    }

    return ptdalgorithms::Vertex(*this, (struct ptd_vertex *) node->entry);
}

ptdalgorithms::Vertex *ptdalgorithms::Graph::find_vertex_p(std::vector<int> state) {
    int *c = force_same_length(state, c_graph());

    Vertex *res = find_vertex_p(c);

    free(c);

    return res;
}

ptdalgorithms::Vertex *ptdalgorithms::Graph::find_vertex_p(const int *state) {
    struct ptd_avl_node *node = ptd_avl_tree_find(this->rf_graph->tree, state);

    if (node == NULL) {
        throw std::runtime_error(
                "No such vertex found\n"
        );
    }

    return new ptdalgorithms::Vertex(*this, (struct ptd_vertex *) node->entry);
}

bool ptdalgorithms::Graph::vertex_exists(std::vector<int> state) {
    assert_same_length(state, this->rf_graph->graph);

    return vertex_exists(&state[0]);
}

bool ptdalgorithms::Graph::vertex_exists(const int *state) {
    struct ptd_avl_node *node = ptd_avl_tree_find(this->rf_graph->tree, state);

    return (node != NULL);
}

ptdalgorithms::Vertex ptdalgorithms::Graph::find_or_create_vertex(std::vector<int> state) {
    assert_same_length(state, this->rf_graph->graph);

    return find_or_create_vertex(&state[0]);
}

ptdalgorithms::Vertex ptdalgorithms::Graph::find_or_create_vertex(const int *state) {
    notify_change();

    return ptdalgorithms::Vertex(*this, ptd_find_or_create_vertex(c_graph(), c_avl_tree(), state));
}

ptdalgorithms::Vertex *ptdalgorithms::Graph::find_or_create_vertex_p(std::vector<int> state) {
    assert_same_length(state, this->rf_graph->graph);

    return find_or_create_vertex_p(&state[0]);
}

ptdalgorithms::Vertex *ptdalgorithms::Graph::find_or_create_vertex_p(const int *state) {
    notify_change();

    return new ptdalgorithms::Vertex(*this, ptd_find_or_create_vertex(c_graph(), c_avl_tree(), state));
}

ptdalgorithms::Vertex ptdalgorithms::Graph::starting_vertex() {
    return Vertex(*this, this->rf_graph->graph->starting_vertex);
}

ptdalgorithms::Vertex *ptdalgorithms::Graph::starting_vertex_p() {
    return new Vertex(*this, this->rf_graph->graph->starting_vertex);
}

std::vector<ptdalgorithms::Vertex> ptdalgorithms::Graph::vertices() {
    std::vector<Vertex> vec;

    for (size_t i = 0; i < c_graph()->vertices_length; ++i) {
        vec.push_back(Vertex(*this, c_graph()->vertices[i]));
    }

    return vec;
}

std::vector<ptdalgorithms::Vertex *> ptdalgorithms::Graph::vertices_p() {
    std::vector<Vertex *> vec;

    for (size_t i = 0; i < c_graph()->vertices_length; ++i) {
        vec.push_back(new Vertex(*this, c_graph()->vertices[i]));
    }

    return vec;
}

ptdalgorithms::Vertex ptdalgorithms::Graph::vertex_at(size_t index) {
    return Vertex(*this, c_graph()->vertices[index]);
}

ptdalgorithms::Vertex *ptdalgorithms::Graph::vertex_at_p(size_t index) {
    return new Vertex(*this, c_graph()->vertices[index]);
}

size_t ptdalgorithms::Graph::vertices_length() {
    return c_graph()->vertices_length;
}

ptdalgorithms::PhaseTypeDistribution ptdalgorithms::Graph::phase_type_distribution() {
    struct ptd_phase_type_distribution *matrix = ptd_graph_as_phase_type_distribution(this->rf_graph->graph);

    if (matrix == NULL) {
        char msg[1024];

        snprintf(msg, 1024, "Failed to make sub-intensity matrix: %s \n", std::strerror(errno));

        throw new std::runtime_error(
                msg
        );
    }

    return PhaseTypeDistribution(*this, matrix);
}

void ptdalgorithms::Vertex::add_edge(Vertex &to, double weight) {
    if (this->vertex == to.vertex) {
        throw new std::invalid_argument(
                "The edge to add is between the same vertex\n"
        );
    }

    graph.notify_change();

    ptd_graph_add_edge(this->vertex, to.vertex, weight);
}


void ptdalgorithms::Vertex::add_edge_parameterized(Vertex &to, double weight, std::vector<double> edge_state) {
    if (this->vertex == to.vertex) {
        throw new std::invalid_argument(
                "The edge to add is between the same vertex\n"
        );
    }

    double *state = (double *) calloc(edge_state.size(), sizeof(*state));

    for (size_t i = 0; i < edge_state.size(); ++i) {
        state[i] = edge_state[i];
    }

    graph.notify_change();

    ptd_graph_add_edge_parameterized(
            this->vertex, to.vertex, weight, state
    );
}


std::vector<int> ptdalgorithms::Vertex::state() {
    return std::vector<int>(
            this->vertex->state,
            this->vertex->state + this->vertex->graph->state_length
    );
}

std::vector<ptdalgorithms::Edge> ptdalgorithms::Vertex::edges() {
    std::vector<Edge> vector;

    for (size_t i = 0; i < this->vertex->edges_length; ++i) {
        Edge edge_i(
                this->vertex->edges[i]->to,
                this->vertex->edges[i],
                graph,
                this->vertex->edges[i]->weight
        );

        vector.push_back(edge_i);
    }

    return vector;
}

std::vector<ptdalgorithms::ParameterizedEdge> ptdalgorithms::Vertex::parameterized_edges() {
    std::vector<ParameterizedEdge> vector;

    for (size_t i = 0; i < this->vertex->edges_length; ++i) {
        if (this->vertex->edges[i]->parameterized) {
            ParameterizedEdge edge_i(
                    this->vertex->edges[i]->to,
                    this->vertex->edges[i],
                    graph,
                    this->vertex->edges[i]->weight,
                    ((struct ptd_edge_parameterized *) this->vertex->edges[i])->state

            );

            vector.push_back(edge_i);
        } else {
            ParameterizedEdge edge_i(
                    this->vertex->edges[i]->to,
                    this->vertex->edges[i],
                    graph,
                    this->vertex->edges[i]->weight,
                    NULL
            );

            vector.push_back(edge_i);
        }
    }

    return vector;
}

ptdalgorithms::Graph ptdalgorithms::Graph::reward_transform(std::vector<double> rewards) {
    struct ptd_graph *res = ptd_graph_reward_transform(this->c_graph(), &rewards[0]);

    if (res == NULL) {
        throw std::runtime_error((char *) ptd_err);
    }

    return Graph(res);
}

ptdalgorithms::Graph *ptdalgorithms::Graph::reward_transform_p(std::vector<double> rewards) {
  struct ptd_graph *res = ptd_graph_reward_transform(this->c_graph(), &rewards[0]);
  
  if (res == NULL) {
    throw std::runtime_error((char *) ptd_err);
  }
  
  return new Graph(res);
}

void ptdalgorithms::Graph::update_weights_parameterized(std::vector<double> scalars) {
    ptd_graph_update_weight_parameterized(
            this->c_graph(),
            &scalars[0],
            scalars.size()
    );

    notify_change();
}

#endif