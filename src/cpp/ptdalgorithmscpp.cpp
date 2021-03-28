#include "../c/io.h"
#include "../../api/c/ptdalgorithms.h"
#include "../../api/cpp/ptdalgorithmscpp.h"

void ptdalgorithms::Graph::create_vertex() {
    Vertex vertex(this);
}

/*

class PhaseTypeGraph {
public:
    PhaseTypeGraph(vertex_t *start, size_t rewards_length) {
        this->start = start;
        this->rewards_length = rewards_length;
    }

    ~PhaseTypeGraph() {
        graph_free(this->start);
    }

    vertex_t *start;
    size_t rewards_length;
};

class PhaseTypeVertex {
public:
    PhaseTypeVertex(vertex_t *vertex) {
        this->vertex = vertex;
    }

    vertex_t *vertex;
};

SEXP get_first_list_entry(SEXP e, std::string message) {
    if (Rf_isList(e) || Rf_isNewList(e)) {
        List list = Rcpp::as<List>(e);

        if (list.size() != 1) {
            char message[1024];

            snprintf(
                    message,
                    1024,
                    "Failed: When finding %s of a list-type of vertices, list can only contain 1 vertex, it contained %zu. Did you use [] as lookup instead of [[]]?",
                    message,
                    (size_t)list.size()
            );

            throw std::runtime_error(
                    message
            );
        }

        e = list[0];
    }

    return e;
}

// [[Rcpp::export]]
IntegerVector state(SEXP phase_type_vertex) {
    phase_type_vertex = get_first_list_entry(phase_type_vertex, "state");
    Rcpp::XPtr<PTDVertex> vertex(phase_type_vertex);
    size_t *state = vertex->vertex->state;
    size_t len = vertex->vertex->graph->state_length;
    IntegerVector res(len);

    for (size_t i = 0; i < len; i++) {
        res(i) = state[i];
    }

    return res;
}


// [[Rcpp::export]]
NumericVector rewards(SEXP phase_type_vertex) {
    phase_type_vertex = get_first_list_entry(phase_type_vertex, "rewards");
    Rcpp::XPtr<PTDVertex> vertex(phase_type_vertex);
    long double *rewards = vertex->vertex->rewards;
    size_t len = vertex->vertex->graph->rewards_length;
    NumericVector res(len);

    for (size_t i = 0; i < len; i++) {
        res(i) = rewards[i];
    }

    return res;
}


// [[Rcpp::export]]
List edges(SEXP phase_type_vertex) {
    phase_type_vertex = get_first_list_entry(phase_type_vertex, "edges");

    Rcpp::XPtr<PTDVertex> vertex(phase_type_vertex);
    List edges(vertex->vertex->edges_length);

    for (size_t i = 0; i < vertex->vertex->edges_length; i++) {
        SEXP child = Rcpp::XPtr<PTDVertex>(
                new PTDVertex(
                        vertex->vertex->edges[i].to
                )
        );

        double weight = vertex->vertex->edges[i].weight;

        edges[i] = List::create(Named("weight") = weight, _["child"] = child);
    }

    return edges;
}

// [[Rcpp::export]]
SEXP create_graph(size_t state_length, size_t reward_length) {
    ptd_graph_t *graph = ptd_graph_create(
            state_length,
            reward_length
    );

    if (graph == NULL) {
        throw new std::runtime_error(
                "memory allocation failed"
        );
    }

    return Rcpp::XPtr<PTDGraph>(
            new PTDGraph(
                    graph
            )
    );
}

// [[Rcpp::export]]
void add_edge(SEXP phase_type_vertex_from, SEXP phase_type_vertex_to, double weight) {
    phase_type_vertex_from = get_first_list_entry(phase_type_vertex_from, "edge from");
    phase_type_vertex_to = get_first_list_entry(phase_type_vertex_to, "edge to");

    Rcpp::XPtr<PTDVertex> from(phase_type_vertex_from);
    Rcpp::XPtr<PTDVertex> to(phase_type_vertex_to);

    if (from->vertex == to->vertex) {
        throw new std::runtime_error(
                "The edge to add is from the same vertex to the same vertex."
        );
    }

    ptd_add_edge(from->vertex, to->vertex, weight);
}

// [[Rcpp::export]]
SEXP create_vertex(SEXP phase_type_graph, IntegerVector state, NumericVector rewards) {
    Rcpp::XPtr<PTDGraph> graph(phase_type_graph);
    ptd_vertex_t *vertex = ptd_vertex_create(graph->graph);

    if (vertex == NULL) {
        throw new std::runtime_error(
                "memory allocation failed"
        );
    }

    for (size_t i = 0; i < vertex->graph->state_length; i++) {
        vertex->state[i] = state[i];
    }

    for (size_t i = 0; i < vertex->graph->rewards_length; i++) {
        vertex->rewards[i] = rewards[i];
    }

    //std::copy(state.begin(), state.end(), vertex->state);
    //std::copy(rewards.begin(), rewards.end(), vertex->rewards);

    ptd_avl_tree_vertex_insert(graph->tree, vertex->state, vertex);

    return Rcpp::XPtr<PTDVertex>(
            new PTDVertex(
                    vertex
            )
    );
}

int print_func(ptd_vertex_t *vertex) {
    return 0;
}

// [[Rcpp::export]]
SEXP find_vertex(SEXP phase_type_graph, IntegerVector state) {
    Rcpp::XPtr<PTDGraph> graph(phase_type_graph);

    if (state.size() != graph->graph->state_length) {
        throw new std::runtime_error(
                "state size must be equal to state length TODO"
        );
    }

    vec_entry_t *c_state = (vec_entry_t*) calloc(graph->graph->state_length, sizeof(*c_state));

    for (size_t i = 0; i < graph->graph->state_length; i++) {
        c_state[i] = state[i];
    }

    ptd_vertex_t *vertex = ptd_avl_tree_vertex_find(graph->tree, c_state);

    if (vertex == NULL) {
        return List::get_na();
    }

    return Rcpp::XPtr<PTDVertex>(
            new PTDVertex(
                    vertex
            )
    );
}


// [[Rcpp::export]]
SEXP find_or_create_vertex(SEXP phase_type_graph, IntegerVector state, NumericVector rewards) {
    SEXP vertex = find_vertex(phase_type_graph, state);

    if (vertex == List::get_na()) {
        return create_vertex(phase_type_graph, state, rewards);
    } else {
        return vertex;
    }
}

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