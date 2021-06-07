#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <numeric>
#include <string.h>
#include "../api/cpp/ptdalgorithmscpp.h"

using namespace ptdalgorithms;

void assert(bool a) {
    if (!a) {
        exit(1);
    }
}

void test_basic_ptd_ph_graph() {
    Graph g(4);

    vector<int> state(4);
    state[0] = 0;
    Vertex v1 = g.create_vertex(state);
    Vertex v2 = g.create_vertex(state);
    Vertex v3 = g.create_vertex(state);
    Vertex v4 = g.create_vertex(state);
    Vertex v5 = g.create_vertex(state);

    assert(g.vertices().size() == 5 + 1);
    assert(g.state_length() == 4);
    assert(g.vertices()[0] == g.starting_vertex());
    assert(g.vertices()[1] == v1);
    assert(g.vertices()[2] == v2);
    assert(g.vertices()[3] == v3);
    assert(g.vertices()[4] == v4);
    assert(g.vertices()[5] == v5);
}

void test_basic_ptd_ph_graph_edges() {
    Graph g(4);

    vector<int> state(4);
    state[0] = 0;
    Vertex v1 = g.create_vertex(state);
    Vertex v2 = g.create_vertex(state);
    Vertex v3 = g.create_vertex(state);
    Vertex v4 = g.create_vertex(state);
    Vertex v5 = g.create_vertex(state);

    g.starting_vertex().add_edge(v1, 0.4);
    g.starting_vertex().add_edge(v2, 0.6);
    v3.add_edge(v2, 1.0);
    v3.add_edge(v1, 2.0);
    v3.add_edge(v4, 2.0);

    assert(g.starting_vertex().edges().size() == 2);
    assert(g.starting_vertex().edges()[0].weight() == 0.4);
    assert(g.starting_vertex().edges()[1].weight() == 0.6);
    assert(v3.edges().size() == 3);
    assert(v3.edges()[0].weight() == 1.0);
    assert(v3.edges()[0].to() == v2);
    assert(v3.edges()[1].weight() == 2.0);
    assert(v3.edges()[1].to() == v1);
    assert(v3.edges()[2].weight() == 2.0);
    assert(v3.edges()[2].to() == v4);
}

int main(int argc, char **argv) {
    test_basic_ptd_ph_graph();
    test_basic_ptd_ph_graph_edges();
}
