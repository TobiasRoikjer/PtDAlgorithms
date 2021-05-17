#include "../api/cpp/ptdalgorithmscpp.h"
#include <iostream>
#include <cassert>

size_t n;

int kingman_visit(ptdalgorithms::Graph &graph, ptdalgorithms::Vertex &v) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            vector<size_t> state = v.state();
            double rate = 0;

            if (i == j) {
                if (state[i] < 2) {
                    continue;
                }

                rate = state[i] * (state[i] - 1) / 2;
            } else {
                if (state[i] < 1 || state[j] < 1) {
                    continue;
                }

                rate = state[i] * state[j];
            }

            vector<size_t> child_state(state);
            child_state[i]--;
            child_state[j]--;
            child_state[i + j + 1]++;

            ptdalgorithms::Vertex to = graph.find_or_create_vertex(child_state);

            v.add_edge(to, rate);
        }
    }

    return 0;
}

long double sum;
size_t reward_index;

int compute_probability(ptdalgorithms::Graph &graph, ptdalgorithms::Vertex &v) {
    long double my_prob;

    if (v == graph.start_vertex()) {
        v.setData(malloc(sizeof(long double)));
        *((long double *) v.getData()) = 1;
    }

    my_prob = *((long double *) v.getData());

    vector<ptdalgorithms::Edge> edges = v.edges();

    for (size_t i = 0; i < edges.size(); ++i) {
        long double trans_prob = edges[i].weight() / v.rate();

        if (edges[i].to().getData() == NULL) {
            edges[i].to().setData(malloc(sizeof(long double)));
            *((long double *) edges[i].to().getData()) = 0;
        }

        *((long double *) v.edges()[i].to().getData()) += my_prob * trans_prob;
    }

    return 0;
}

int compute_expected(ptdalgorithms::Graph &graph, ptdalgorithms::Vertex &v) {
    long double my_prob = *((long double *) v.getData());
    long double my_expected;

    if (v.rate() != 0) {
        my_expected = my_prob * v.state()[reward_index] / v.rate();
    } else {
        my_expected = 0;
    }

    *((long double *) v.getData()) = my_expected;

    return 0;
}

int reduce(ptdalgorithms::Graph &graph, ptdalgorithms::Vertex &v) {
    long double data = *((long double *) v.getData());
    sum += data;

    return 0;
}

int p(ptdalgorithms::Graph &graph, ptdalgorithms::Vertex &v) {
    fprintf(stderr, "I am vertex: ");
    fprintf(stderr, "%zu %zu %zu %zu %zu %zu (data=%Lf)\n",
            v.state()[0],
            v.state()[1],
            v.state()[2],
            v.state()[3],
            v.state()[4],
            v.state()[5],
            *((long double *) v.getData())
    );

    return 0;
}

int pi(ptdalgorithms::Graph &graph, ptdalgorithms::Vertex &v) {
    char buf[1024] = {0};
    ptd_vertex_to_s(v.c_vertex(), buf, 1024);
    fprintf(stderr, "%s has index %zu\n", buf, v.c_vertex()->index);

    return 0;
}

int pc(struct ptd_vertex *v) {
    char buffer[1024];
    ptd_vertex_to_s(v, buffer, 1024);
    fprintf(stderr, "I am vertex: %s\n", buffer);

    for (size_t i = 0; i < v->edges_length; ++i) {
        struct ptd_vertex *c = v->edges[i].to;

        char buffer2[1024];
        ptd_vertex_to_s(c, buffer2, 1024);
        fprintf(stderr, "Child (%Lf): %s\n", v->edges[i].weight, buffer2);
    }

    fprintf(stderr, "\n");

    return 0;
}


double reward_by_index(struct ptd_vertex *vertex) {
    return vertex->state[reward_index];
}

double reward_by_1(struct ptd_vertex *vertex) {
    return 1;
}

bool always_true(struct ptd_vertex *v) {
    return true;
}

int main(void) {
    struct ptd_graph *two_loci = ptd_model_two_island_two_loci_recomb(
            5,
            5,
            1,
            1,
            1,
            1,
            1
    );

//    ptd_visit_vertices(two_loci, pc, true);


    ptd_strongly_connected_components_t *scc =
            ptd_find_strongly_connected_components(two_loci, always_true);

    fprintf(stderr, "We have %zu scc with %zu internal_vertices\n", scc->components_length, two_loci->vertices_length);

    fprintf(stderr, "EXP %f\n", ptd_circular_exp(two_loci, reward_by_1));
    fprintf(stderr, "We have %zu scc with %zu internal_vertices\n", scc->components_length, two_loci->vertices_length);
    return 0;

    for (size_t l = 0; l < scc->components_length; ++l) {
        ptd_strongly_connected_component_t *sc = scc->components[l];
        fprintf(stderr, "Component %zu has %zu internal_vertices\n", l, sc->internal_vertices_length);

        for (size_t i = 0; i < sc->internal_vertices_length; ++i) {
            char buffer[1024];
            ptd_vertex_to_s(sc->internal_vertices[i], buffer, 1024);
            fprintf(stderr, "Entry idx %zu: %s\n",sc->internal_vertices[i]->index,  buffer);
        }
    }

    ptd_graph_destroy(two_loci);

    return 0;


    /*vector<size_t> state1;
    state1.push_back(1);
    state1.push_back(2);
    vector<size_t> state2;
    state2.push_back(3);
    state2.push_back(4);
    ptdalgorithms::Graph graph(2);
    ptdalgorithms::Vertex vertex = graph.create_vertex(state1);
    ptdalgorithms::Vertex vertex2 = graph.find_or_create_vertex(state2);
    vertex.add_edge(vertex2, 4);
    ptdalgorithms::Graph g = graph;

    ptdalgorithms::Graph graph2(2);
    g = graph2;
    ptdalgorithms::Vertex e = vertex.edges()[0].to;
    std::cout << "Edges: " << e.state()[0] << "-" << e.state()[1] << "  " << vertex.edges()[0].weight << endl;
    ptdalgorithms::Vertex vertex3 = graph.find_or_create_vertex(state2);
    assert(vertex2 == vertex3);
    assert(vertex2 == e);

    vector<ptdalgorithms::Vertex> list = graph.internal_vertices();

    assert(list.size() == 2);
    assert(list[0].state() == state1);
    assert(list[1].state() == state2);*/

    n = 4;
    ptdalgorithms::Graph kingman(n);

    vector<size_t> initial_state;
    initial_state.push_back(n);

    for (size_t k = 1; k < n; ++k) {
        initial_state.push_back(0);
    }

    ptdalgorithms::Vertex initial = kingman.create_vertex(initial_state);
    kingman.start_vertex().add_edge(initial, 1);

    time_t seconds, seconds2, seconds3, seconds4, seconds5;

    time(&seconds);
    kingman.visit_vertices(kingman_visit);
    time(&seconds2);

    reward_index = 0;

    fprintf(stderr, "%zu Vertices: %zu\n", seconds2 - seconds, kingman.c_graph()->vertices_length);
    time(&seconds3);
    kingman.index_topological();
    kingman.visit_vertices(pi);
    kingman.visit_vertices(compute_probability, true);
    kingman.index_invert();
    kingman.visit_vertices(compute_expected, true);
    sum = 0;
    kingman.visit_vertices(reduce, true);
    fprintf(stderr, "%zu %Lf\n", seconds3 - seconds2, sum);


    time(&seconds4);
    ptdalgorithms::Graph kingman_c = ptdalgorithms::Models::kingman(4);
    time(&seconds5);
    fprintf(stderr, "C kingman %zu\n", seconds5 - seconds4);
    long double exp = 0;
    reward_index = 0;
    ptd_expected_value(&exp, kingman_c.c_graph(), reward_by_index);
    fprintf(stderr, "Kingman exp %Lf\n", exp);
    long double cov = 0;
    reward_index = 0;
    ptd_covariance(&cov, kingman_c.c_graph(), reward_by_index, reward_by_index);
    fprintf(stderr, "Kingman cov %Lf\n", cov);
    return 0;

    kingman_c.index_topological();

    ptdalgorithms::PhaseTypeDistribution sub = kingman_c.phase_type_distribution();

    for (size_t i = 0; i < sub.length; ++i) {
        fprintf(stderr, "%Lf ", sub.initial_probability_vector[i]);
    }

    fprintf(stderr, "\n\n");

    for (size_t i = 0; i < sub.length; ++i) {
        for (size_t j = 0; j < sub.length; ++j) {
            fprintf(stderr, "%Lf ", sub.sub_intensity_matrix[i][j]);
        }

        fprintf(stderr, "\n");
    }

}