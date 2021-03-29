#include <stdexcept>
#include <sstream>
#include "../c/io.h"
#include "../../api/c/ptdalgorithms.h"
#include "../../api/cpp/ptdalgorithmscpp.h"

static void assert_same_length(vector<size_t> state, ptd_graph_t *graph) {
    if (state.size() != graph->state_length) {
        stringstream message;
        message << "Vector `state` argument must have same size as graph state length. Was '";
        message << state.size() << "' expected '" << graph->state_length << "'" << std::endl;

        throw std::invalid_argument(message.str());
    }
}

ptdalgorithms::Vertex ptdalgorithms::Graph::create_vertex(vector<size_t> state) {
    assert_same_length(state, this->graph);

    size_t *c_state = (size_t *) calloc(state.size(), sizeof(*c_state));
    std::copy(state.begin(), state.end(), c_state);

    Vertex vertex = ptdalgorithms::Vertex(*this, c_state);

    if (ptd_avl_tree_vertex_insert(this->tree, c_state, vertex.vertex)) {
        throw runtime_error(
                "Failed to insert into AVL tree\n"
        );
    }

    return vertex;
}

ptdalgorithms::Vertex ptdalgorithms::Graph::find_vertex(vector<size_t> state) const {
    assert_same_length(state, this->graph);

    size_t *c_state = (size_t *) calloc(state.size(), sizeof(*c_state));
    std::copy(state.begin(), state.end(), c_state);

    ptd_vertex_t *vertex = ptd_avl_tree_vertex_find(this->tree, c_state);

    if (vertex == NULL) {
        throw runtime_error(
                "No such vertex found\n"
        );
    }

    free(c_state);

    return ptdalgorithms::Vertex(*this, vertex);
}

bool ptdalgorithms::Graph::vertex_exists(std::vector<size_t> state) {
    assert_same_length(state, this->graph);

    size_t *c_state = (size_t *) calloc(state.size(), sizeof(*c_state));
    std::copy(state.begin(), state.end(), c_state);

    ptd_vertex_t *vertex = ptd_avl_tree_vertex_find(this->tree, c_state);

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

ptdalgorithms::Vertex ptdalgorithms::Graph::start_vertex() const {
    return Vertex(*this, this->graph->start_vertex);
}

std::vector<ptdalgorithms::Vertex> ptdalgorithms::Graph::vertices() const {
    std::vector<ptd_edge_t *> vec;
    std::stack<avl_vec_vertex_t *> s;

    s.push((avl_vec_vertex_t *) avl_tree->root);

    while (!s.empty()) {
        avl_vec_vertex_t *v = s.top();
        s.pop();

        if (v == NULL) {
            continue;
        }

        vec.push_back((ptd_edge_t *) v->entry);
        s.push(v->left);
        s.push(v->right);
    }


    return Vertex(*this, this->graph->start_vertex);
}

void ptdalgorithms::Vertex::add_edge(const Vertex &to, long double weight) {
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
            this->vertex->state + this->graph.graph->state_length
    );
}

std::vector<ptdalgorithms::Edge> ptdalgorithms::Vertex::edges() {
    std::vector<Edge> vector;

    for (size_t i = 0; i < this->vertex->edges_length; ++i) {
        Edge edge_i(
                Vertex(this->graph, this->vertex->edges[i].to),
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
List vertices(SEXP phase_type_graph) {
    Rcpp::XPtr<PTDGraph> graph(phase_type_graph);

    ptd_label_vertices(graph->graph);
    queue<ptd_vertex_t *> q = ptd_enqueue_vertices(graph->graph);

    // Remove start vertex
    q.pop();

    List list(q.size());

    while (!q.empty()) {
        ptd_vertex_t *vertex = q.front();
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

int custom_visit(ptd_vertex_t *vertex) {
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