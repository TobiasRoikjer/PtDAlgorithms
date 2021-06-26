#include <stdexcept>
#include <sstream>
#include <stack>
#include <cerrno>
#include <cstring>
#include "ptdalgorithmscpp.h"

static void assert_same_length(vector<int> state, struct ptd_ph_graph *graph) {
    if (state.size() != graph->state_length) {
        stringstream message;
        message << "Vector `state` argument must have same size as graph state length. Was '";
        message << state.size() << "' expected '" << graph->state_length << "'" << std::endl;

        throw std::invalid_argument(message.str());
    }
}

ptdalgorithms::Vertex ptdalgorithms::Graph::create_vertex(vector<int> state) {
    assert_same_length(state, this->rf_graph->graph);

    int *c_state = (int *) calloc(state.size(), sizeof(*c_state));
    std::copy(state.begin(), state.end(), c_state);

    Vertex vertex = ptdalgorithms::Vertex(*this, c_state);

    if (vertex_exists(state)) {
        char buffer[1024];
        ptd_ph_vertex_to_s(vertex.vertex, buffer, 1024);
        char message[2048];
        snprintf(message, 2048, "Vertex with state '%s' already exists\n", buffer);

        throw runtime_error(
                message
        );
    }

    if (ptd_avl_tree_find_or_insert(this->rf_graph->tree, c_state, vertex.vertex) == NULL) {
        throw runtime_error(
                "Failed to insert into AVL tree\n"
        );
    }

    return vertex;
}

ptdalgorithms::Vertex *ptdalgorithms::Graph::create_vertex_p(vector<int> state) {
    assert_same_length(state, this->rf_graph->graph);

    int *c_state = (int *) calloc(state.size(), sizeof(*c_state));
    std::copy(state.begin(), state.end(), c_state);

    Vertex *vertex = new ptdalgorithms::Vertex(*this, c_state);

    if (vertex_exists(state)) {
        char buffer[1024];
        ptd_ph_vertex_to_s(vertex->vertex, buffer, 1024);
        char message[2048];
        snprintf(message, 2048, "Vertex with state '%s' already exists\n", buffer);

        throw runtime_error(
                message
        );
    }

    if (ptd_avl_tree_find_or_insert(this->rf_graph->tree, c_state, vertex->vertex) == NULL) {
        throw runtime_error(
                "Failed to insert into AVL tree\n"
        );
    }

    return vertex;
}

ptdalgorithms::Vertex ptdalgorithms::Graph::find_vertex(vector<int> state) {
    assert_same_length(state, this->rf_graph->graph);

    int *c_state = (int *) calloc(state.size(), sizeof(*c_state));
    std::copy(state.begin(), state.end(), c_state);

    struct ptd_avl_node *node = ptd_avl_tree_find(this->rf_graph->tree, c_state);

    if (node == NULL) {
        throw runtime_error(
                "No such vertex found\n"
        );
    }

    free(c_state);

    return ptdalgorithms::Vertex(*this, (struct ptd_ph_vertex*) node->entry);
}

ptdalgorithms::Vertex *ptdalgorithms::Graph::find_vertex_p(vector<int> state) {
    assert_same_length(state, this->rf_graph->graph);

    int *c_state = (int *) calloc(state.size(), sizeof(*c_state));
    std::copy(state.begin(), state.end(), c_state);

    struct ptd_avl_node *node = ptd_avl_tree_find(this->rf_graph->tree, c_state);

    if (node == NULL) {
        throw runtime_error(
                "No such vertex found\n"
        );
    }

    free(c_state);

    return new ptdalgorithms::Vertex(*this, (struct ptd_ph_vertex*) node->entry);
}

bool ptdalgorithms::Graph::vertex_exists(std::vector<int> state) {
    assert_same_length(state, this->rf_graph->graph);

    int *c_state = (int *) calloc(state.size(), sizeof(*c_state));
    std::copy(state.begin(), state.end(), c_state);

    struct ptd_avl_node *node = ptd_avl_tree_find(this->rf_graph->tree, c_state);

    free(c_state);

    return (node != NULL);
}

ptdalgorithms::Vertex ptdalgorithms::Graph::find_or_create_vertex(vector<int> state) {
    if (vertex_exists(state)) {
        return find_vertex(state);
    } else {
        return this->create_vertex(state);
    }
}

ptdalgorithms::Vertex *ptdalgorithms::Graph::find_or_create_vertex_p(vector<int> state) {
    if (vertex_exists(state)) {
        return find_vertex_p(state);
    } else {
        return this->create_vertex_p(state);
    }
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

std::vector<ptdalgorithms::Vertex*> ptdalgorithms::Graph::vertices_p() {
    std::vector<Vertex*> vec;

    for (size_t i = 0; i < c_graph()->vertices_length; ++i) {
        vec.push_back(new Vertex(*this, c_graph()->vertices[i]));
    }

    return vec;
}

ptdalgorithms::PhaseTypeDistribution ptdalgorithms::Graph::phase_type_distribution() {
    ptd_phase_type_distribution_t *matrix = ptd_graph_as_phase_type_distribution(this->rf_graph->graph);

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

    ptd_ph_graph_add_edge(this->vertex, to.vertex, weight);
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
                graph,
                this->vertex->edges[i]->weight
        );

        vector.push_back(edge_i);
    }

    return vector;
}

std::vector<double> ptdalgorithms::Graph::expected_visits() {
    double *ptr = ptd_ph_graph_expected_visits(this->c_graph());
    std::vector<double> res;
    res.assign(ptr, ptr + this->c_graph()->vertices_length);
    free(ptr);

    return res;
}

std::vector<double> ptdalgorithms::Graph::expected_waiting_time() {
    double *ptr = ptd_ph_graph_expected_waiting_time(this->c_graph());
    std::vector<double> res;
    res.assign(ptr, ptr + this->c_graph()->vertices_length);
    free(ptr);

    return res;
}

std::vector<double> ptdalgorithms::Graph::moment_rewards(std::vector<double> rewards) {
    double *ptr = ptd_ph_graph_moment_rewards(this->c_graph(), &rewards[0]);
    std::vector<double> res;
    res.assign(ptr, ptr + this->c_graph()->vertices_length);
    free(ptr);

    return res;
}

void ptdalgorithms::Graph::reward_transform(std::vector<double> rewards) {
    ptd_ph_graph_reward_transform(this->c_graph(), &rewards[0]);
}
