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

    Graph g2 = g;
    Graph g3(g);
    Vertex v10 = v3;
    Vertex v11(v3);

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


void test_expected_entry_visits() {
    struct ptd_ph_graph *graph = ptd_ph_graph_create(4);

    struct ptd_ph_vertex *S = graph->starting_vertex;
    S->state[0] = 0;
    fprintf(stderr, "S: %p\n", (void *) S);
    struct ptd_ph_vertex *T = ptd_ph_vertex_create(graph);
    T->state[0] = 1;
    fprintf(stderr, "T: %p\n", (void *) T);

    struct ptd_ph_vertex *A = ptd_ph_vertex_create(graph);
    A->state[0] = 2;
    fprintf(stderr, "A: %p\n", (void *) A);
    struct ptd_ph_vertex *B = ptd_ph_vertex_create(graph);
    B->state[0] = 3;
    fprintf(stderr, "B: %p\n", (void *) B);
    struct ptd_ph_vertex *C = ptd_ph_vertex_create(graph);
    C->state[0] = 4;
    fprintf(stderr, "C: %p\n", (void *) C);
    struct ptd_ph_vertex *D = ptd_ph_vertex_create(graph);
    D->state[0] = 5;
    fprintf(stderr, "D: %p\n", (void *) D);
    struct ptd_ph_vertex *E = ptd_ph_vertex_create(graph);
    E->state[0] = 6;
    fprintf(stderr, "E: %p\n", (void *) E);
    struct ptd_ph_vertex *F = ptd_ph_vertex_create(graph);
    F->state[0] = 7;
    fprintf(stderr, "F: %p\n", (void *) F);
    struct ptd_ph_vertex *G = ptd_ph_vertex_create(graph);
    G->state[0] = 8;
    fprintf(stderr, "G: %p\n", (void *) G);
    struct ptd_ph_vertex *H = ptd_ph_vertex_create(graph);
    H->state[0] = 9;
    fprintf(stderr, "H: %p\n", (void *) H);
    struct ptd_ph_vertex *I = ptd_ph_vertex_create(graph);
    I->state[0] = 10;
    fprintf(stderr, "I: %p\n", (void *) I);
    struct ptd_ph_vertex *J = ptd_ph_vertex_create(graph);
    J->state[0] = 11;
    fprintf(stderr, "J %p\n", (void *) J);
    struct ptd_ph_vertex *K = ptd_ph_vertex_create(graph);
    K->state[0] = 12;
    fprintf(stderr, "K %p\n", (void *) K);
    struct ptd_ph_vertex *L = ptd_ph_vertex_create(graph);
    L->state[0] = 13;
    fprintf(stderr, "L %p\n", (void *) L);


    ptd_ph_graph_add_edge(S, A, 1);

    ptd_ph_graph_add_edge(A, B, 3);
    ptd_ph_graph_add_edge(B, C, 4);
    ptd_ph_graph_add_edge(C, A, 4);
    ptd_ph_graph_add_edge(B, D, 2);
    ptd_ph_graph_add_edge(D, E, 5);
    ptd_ph_graph_add_edge(E, F, 5);
    ptd_ph_graph_add_edge(E, I, 15);
    ptd_ph_graph_add_edge(F, G, 1);
    ptd_ph_graph_add_edge(G, H, 1);
    ptd_ph_graph_add_edge(H, F, 1);
    ptd_ph_graph_add_edge(H, G, 1);

    ptd_ph_graph_add_edge(H, T, 1);
    ptd_ph_graph_add_edge(I, L, 1);
    ptd_ph_graph_add_edge(L, T, 1);

    ptd_ph_graph_add_edge(C, J, 1);

    ptd_ph_graph_add_edge(J, K, 1);
    ptd_ph_graph_add_edge(K, J, 1);

    ptd_ph_graph_add_edge(K, E, 1);

    ptd_ph_graph_add_edge(K, G, 1);

    Graph g(graph);
    vector<double> exp_visits = g.expected_visits();

    assert(abs(exp_visits[S->index] - 1) < 0.01);
    assert(abs(exp_visits[T->index] - 1) < 0.01);
    assert(abs(exp_visits[A->index] - 2.14) < 0.01);
    assert(abs(exp_visits[B->index] - 2.14) < 0.01);
    assert(abs(exp_visits[C->index] - 1.42) < 0.01);
    assert(abs(exp_visits[D->index] - 0.71) < 0.01);
    assert(abs(exp_visits[E->index] - 0.85) < 0.01);
    assert(abs(exp_visits[J->index] - 0.42) < 0.01);
    assert(abs(exp_visits[K->index] - 0.42) < 0.01);
    assert(abs(exp_visits[G->index] - 1.07) < 0.01);
    assert(abs(exp_visits[H->index] - 1.07) < 0.01);
}

int main(int argc, char **argv) {
    test_basic_ptd_ph_graph();
    test_basic_ptd_ph_graph_edges();
    test_expected_entry_visits();
}
