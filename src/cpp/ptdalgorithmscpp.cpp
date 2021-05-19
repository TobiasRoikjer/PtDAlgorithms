#include <stdexcept>
#include <sstream>
#include <stack>
#include <cerrno>
#include <cstring>
#include "../c/io.h"
#include "ptdalgorithmscpp.h"

static void assert_same_length(vector<size_t> state, struct ptd_graph *graph) {
    if (state.size() != graph->state_length) {
        stringstream message;
        message << "Vector `state` argument must have same size as graph state length. Was '";
        message << state.size() << "' expected '" << graph->state_length << "'" << std::endl;

        throw std::invalid_argument(message.str());
    }
}

ptdalgorithms::Vertex ptdalgorithms::VertexLinkedList::next(void) {
    if (this->is_first) {
        if (this->current == NULL || this->current->vertex == NULL) {
            throw std::runtime_error("No next vertex (is NULL). Did has_next return true?");
        }

        Vertex res(graph, this->current->vertex);

        this->is_first = false;

        return res;
    } else {
        if (this->current == NULL) {
            throw std::runtime_error("No next vertex (is NULL). Did has_next return true?");
        }

        if (this->current->next == NULL || this->current->next->vertex == NULL) {
            throw std::runtime_error("No next vertex (is NULL). Did has_next return true?");
        }

        this->current = this->current->next;

        Vertex res(graph, this->current->vertex);

        return res;
    }
}

ptdalgorithms::Vertex ptdalgorithms::Graph::create_vertex(vector<size_t> state) {
    assert_same_length(state, this->rf_graph->graph);

    size_t *c_state = (size_t *) calloc(state.size(), sizeof(*c_state));
    std::copy(state.begin(), state.end(), c_state);

    Vertex vertex = ptdalgorithms::Vertex(*this, c_state);

    if (ptd_avl_tree_vertex_insert(this->rf_graph->tree, c_state, vertex.vertex)) {
        throw runtime_error(
                "Failed to insert into AVL tree\n"
        );
    }

    return vertex;
}

ptdalgorithms::Vertex *ptdalgorithms::Graph::create_vertex_p(vector<size_t> state) {
    assert_same_length(state, this->rf_graph->graph);

    size_t *c_state = (size_t *) calloc(state.size(), sizeof(*c_state));
    std::copy(state.begin(), state.end(), c_state);

    Vertex *vertex = new ptdalgorithms::Vertex(*this, c_state);

    if (ptd_avl_tree_vertex_insert(this->rf_graph->tree, c_state, vertex->vertex)) {
        throw runtime_error(
                "Failed to insert into AVL tree\n"
        );
    }

    return vertex;
}

ptdalgorithms::Vertex ptdalgorithms::Graph::find_vertex(vector<size_t> state) {
    assert_same_length(state, this->rf_graph->graph);

    size_t *c_state = (size_t *) calloc(state.size(), sizeof(*c_state));
    std::copy(state.begin(), state.end(), c_state);

    struct ptd_vertex *vertex = ptd_avl_tree_vertex_find(this->rf_graph->tree, c_state);

    if (vertex == NULL) {
        throw runtime_error(
                "No such vertex found\n"
        );
    }

    free(c_state);

    return ptdalgorithms::Vertex(*this, vertex);
}

ptdalgorithms::Vertex *ptdalgorithms::Graph::find_vertex_p(vector<size_t> state) {
    assert_same_length(state, this->rf_graph->graph);

    size_t *c_state = (size_t *) calloc(state.size(), sizeof(*c_state));
    std::copy(state.begin(), state.end(), c_state);

    struct ptd_vertex *vertex = ptd_avl_tree_vertex_find(this->rf_graph->tree, c_state);

    if (vertex == NULL) {
        throw runtime_error(
                "No such vertex found\n"
        );
    }

    free(c_state);

    return new ptdalgorithms::Vertex(*this, vertex);
}

bool ptdalgorithms::Graph::vertex_exists(std::vector<size_t> state) {
    assert_same_length(state, this->rf_graph->graph);

    size_t *c_state = (size_t *) calloc(state.size(), sizeof(*c_state));
    std::copy(state.begin(), state.end(), c_state);

    struct ptd_vertex *vertex = ptd_avl_tree_vertex_find(this->rf_graph->tree, c_state);

    free(c_state);

    return (vertex != NULL);
}

ptdalgorithms::Vertex ptdalgorithms::Graph::find_or_create_vertex(vector<size_t> state) {
    if (vertex_exists(state)) {
        return find_vertex(state);
    } else {
        return this->create_vertex(state);
    }
}

ptdalgorithms::Vertex *ptdalgorithms::Graph::find_or_create_vertex_p(vector<size_t> state) {
    if (vertex_exists(state)) {
        return find_vertex_p(state);
    } else {
        return this->create_vertex_p(state);
    }
}

ptdalgorithms::Vertex ptdalgorithms::Graph::start_vertex() {
    return Vertex(*this, this->rf_graph->graph->start_vertex);
}

ptdalgorithms::Vertex *ptdalgorithms::Graph::start_vertex_p() {
    return new Vertex(*this, this->rf_graph->graph->start_vertex);
}

std::vector<ptdalgorithms::Vertex> ptdalgorithms::Graph::vertices() {
    std::vector<Vertex> vec;
    std::stack<struct avl_node *> s;

    s.push((struct avl_node *) rf_graph->tree->root);

    while (!s.empty()) {
        struct avl_node *v = s.top();
        s.pop();

        if (v == NULL) {
            continue;
        }

        vec.push_back(Vertex(*this, ((struct ptd_vertex *) v->entry)));
        s.push(v->left);
        s.push(v->right);
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

static int (*cpp_visit_func)(ptdalgorithms::Graph &graph, ptdalgorithms::Vertex &vertex);

static ptdalgorithms::Graph *cpp_graph;

static int visit_from_cpp(struct ptd_vertex *vertex) {
    ptdalgorithms::Vertex t(*cpp_graph, vertex);

    return cpp_visit_func(*cpp_graph, t);
}

void ptdalgorithms::Graph::visit_vertices(
        int (*visit_func)(ptdalgorithms::Graph &graph, ptdalgorithms::Vertex &vertex),
        bool include_start
) {
    cpp_visit_func = visit_func;
    cpp_graph = this;
    int res = ptd_visit_vertices(this->rf_graph->graph, visit_from_cpp, include_start);

    if (res != 0) {
        char msg[1024];

        snprintf(msg, 1024, "Failed to visit internal_vertices, visit function returned non-zero: %s \n", std::strerror(errno));

        throw new std::runtime_error(
                msg
        );
    }
}

void ptdalgorithms::Vertex::add_edge(Vertex &to, long double weight) {
    if (this->vertex == to.vertex) {
        throw new std::invalid_argument(
                "The edge to add is between the same vertex\n"
        );
    }

    ptd_add_edge(this->vertex, to.vertex, weight);
}

std::vector<size_t> ptdalgorithms::Vertex::state() {
    return std::vector<size_t>(
            this->vertex->state,
            this->vertex->state + this->vertex->graph->state_length
    );
}

std::vector<ptdalgorithms::Edge> ptdalgorithms::Vertex::edges() {
    std::vector<Edge> vector;

    for (size_t i = 0; i < this->vertex->edges_length; ++i) {
        Edge edge_i(
                this->vertex->edges[i].to,
                graph,
                this->vertex->edges[i].weight
        );

        vector.push_back(edge_i);
    }

    return vector;
}

/*
// [[Rcpp::export]]
SEXP start_vertex(SEXP phase_type_graph) {
    Rcpp::XPtr<PTDGraph> graph(phase_type_graph);

    return Rcpp::XPtr<PTDVertex>(
            new PTDVertex(
                    graph->graph->start_vertex
            )
    );
}

// [[Rcpp::export]]
List internal_vertices(SEXP phase_type_graph) {
    Rcpp::XPtr<PTDGraph> graph(phase_type_graph);

    ptd_label_vertices(graph->graph);
    queue<struct ptd_vertex *> q = ptd_enqueue_vertices(graph->graph);

    // Remove start vertex
    q.pop();

    List list(q.size());

    while (!q.empty()) {
        struct ptd_vertex *vertex = q.front();
        q.pop();

        list[vertex->index - 1] = Rcpp::XPtr<PTDVertex>(
                new PTDVertex(
                        vertex
                )
        );
    }

    return list;
}

Rcpp::Function *custom_visit_function;

int custom_visit(struct ptd_vertex *vertex) {
    fprintf(stderr, "foo\n");
    SEXP v = Rcpp::XPtr<PTDVertex>(
            new PTDVertex(
                    vertex
            )
    );
    fprintf(stderr, "baz\n");

    (*custom_visit_function)(
            v
    );
    fprintf(stderr, "bar\n");

    return 0;
}

// [[Rcpp::export]]
void visit_vertices(SEXP phase_type_graph, Rcpp::Function visit_function) {
    Rcpp::XPtr<PTDGraph> graph(phase_type_graph);
    fprintf(stderr, "W");
    custom_visit_function = &visit_function;

    ptd_visit_vertices(graph->graph, custom_visit);
}

*/