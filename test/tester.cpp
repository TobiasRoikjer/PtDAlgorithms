#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <numeric>
#include <string.h>
#include "../api/c/ptdalgorithms.h"
#include <math.h>


void assert(bool a) {
    if (!a) {
        exit(1);
    }
}


void test_it_can_insert_avl() {
    struct ptd_graph *graph = ptd_graph_create(4);
    assert(graph != NULL);
    struct ptd_avl_tree *tree = ptd_avl_tree_create(4);
    assert(tree != NULL);

    struct ptd_vertex *S = graph->starting_vertex;
    fprintf(stderr, "S: %p\n", (void *) S);
    struct ptd_vertex *T = ptd_vertex_create(graph);
    fprintf(stderr, "T: %p\n", (void *) T);

    struct ptd_vertex *A = ptd_vertex_create(graph);
    fprintf(stderr, "A: %p\n", (void *) A);
    struct ptd_vertex *B = ptd_vertex_create(graph);
    fprintf(stderr, "B: %p\n", (void *) B);
    struct ptd_vertex *C = ptd_vertex_create(graph);
    fprintf(stderr, "C: %p\n", (void *) C);
    struct ptd_vertex *D = ptd_vertex_create(graph);
    fprintf(stderr, "D: %p\n", (void *) D);
    struct ptd_vertex *E = ptd_vertex_create(graph);
    fprintf(stderr, "E: %p\n", (void *) E);
    struct ptd_vertex *F = ptd_vertex_create(graph);
    fprintf(stderr, "F: %p\n", (void *) F);
    struct ptd_vertex *G = ptd_vertex_create(graph);
    fprintf(stderr, "G: %p\n", (void *) G);
    struct ptd_vertex *H = ptd_vertex_create(graph);
    fprintf(stderr, "H: %p\n", (void *) H);
    struct ptd_vertex *I = ptd_vertex_create(graph);
    fprintf(stderr, "I: %p\n", (void *) I);
    struct ptd_vertex *J = ptd_vertex_create(graph);
    fprintf(stderr, "J %p\n", (void *) J);
    struct ptd_vertex *K = ptd_vertex_create(graph);
    fprintf(stderr, "K %p\n", (void *) K);
    struct ptd_vertex *L = ptd_vertex_create(graph);
    fprintf(stderr, "L %p\n", (void *) L);

    S->state[0] = 2;
    A->state[0] = 4;
    B->state[0] = 3;
    C->state[0] = 3;
    C->state[1] = 2;
    D->state[0] = 3;
    D->state[1] = 1;
    E->state[0] = 1;

    assert(ptd_avl_tree_find(tree, S->state) == NULL);

    ptd_avl_tree_find_or_insert(tree, S->state, S);
    assert(ptd_avl_tree_find(tree, S->state)->entry == S);

    assert(ptd_avl_tree_find(tree, A->state) == NULL);
    ptd_avl_tree_find_or_insert(tree, A->state, A);
    assert(ptd_avl_tree_find(tree, A->state)->entry == A);

    ptd_avl_tree_find_or_insert(tree, B->state, B);
    assert(ptd_avl_tree_find(tree, B->state)->entry == B);
    ptd_avl_tree_find_or_insert(tree, C->state, C);
    assert(ptd_avl_tree_find(tree, C->state)->entry == C);
    assert(ptd_avl_tree_find(tree, D->state) == NULL);
    ptd_avl_tree_find_or_insert(tree, D->state, D);
    assert(ptd_avl_tree_find(tree, D->state)->entry == D);
    ptd_avl_tree_find_or_insert(tree, E->state, E);
    assert(ptd_avl_tree_find(tree, E->state)->entry == E);

    ptd_avl_tree_destroy(tree);
//destroy_vertices();
    ptd_graph_destroy(graph);
}

#include <vector>
#include <cmath>
#include <time.h>

using namespace std;

void test_avl_is_balanced() {
    struct ptd_graph *graph = ptd_graph_create(4);
    struct ptd_avl_tree *tree = ptd_avl_tree_create(4);
    vector<struct ptd_vertex *> vertices;

    for (size_t i = 0; i < 1024; ++i) {
        struct ptd_vertex *vertex = ptd_vertex_create(graph);
        assert(vertex != NULL);
        vertex->state[0] = rand();
        ptd_avl_tree_find_or_insert(tree, vertex->state, vertex);
        assert(ptd_avl_tree_find(tree, vertex->state)->entry == vertex);
        vertices.push_back(vertex);
    }

    for (size_t k = 0; k < vertices.size(); ++k) {
        ptd_vertex_destroy(vertices[k]);
    }

    fprintf(stderr, "Depth: %zu\n", ptd_avl_tree_max_depth(tree->root));
    assert(ptd_avl_tree_max_depth(tree->root) < 16);
    assert(ptd_avl_tree_max_depth(tree->root) > 2);
    ptd_avl_tree_destroy(tree);
    ptd_graph_destroy(graph);
}


void test_basic_graph() {
    struct ptd_directed_graph *g = (struct ptd_directed_graph *) calloc(1, sizeof(*g));
    struct ptd_directed_vertex *v1 = (struct ptd_directed_vertex *) calloc(1, sizeof(*v1));
    struct ptd_directed_vertex *v2 = (struct ptd_directed_vertex *) calloc(1, sizeof(*v2));
    struct ptd_directed_vertex *v3 = (struct ptd_directed_vertex *) calloc(1, sizeof(*v3));
    struct ptd_directed_vertex *v4 = (struct ptd_directed_vertex *) calloc(1, sizeof(*v4));
    struct ptd_directed_vertex *v5 = (struct ptd_directed_vertex *) calloc(1, sizeof(*v5));

    ptd_directed_vertex_add(g, v1);
    ptd_directed_vertex_add(g, v2);
    ptd_directed_vertex_add(g, v3);
    ptd_directed_vertex_add(g, v4);
    ptd_directed_vertex_add(g, v5);

    assert(g->vertices_length == 5);
    assert(g->vertices[0] == v1);
    assert(g->vertices[1] == v2);
    assert(g->vertices[2] == v3);
    assert(g->vertices[3] == v4);
    assert(g->vertices[4] == v5);

    ptd_directed_graph_destroy(g);
}

void test_basic_ptd_graph() {
    struct ptd_graph *g = ptd_graph_create(4);
    struct ptd_vertex *v1 = ptd_vertex_create(g);
    struct ptd_vertex *v2 = ptd_vertex_create(g);
    struct ptd_vertex *v3 = ptd_vertex_create(g);
    struct ptd_vertex *v4 = ptd_vertex_create(g);
    struct ptd_vertex *v5 = ptd_vertex_create(g);

    assert(g->vertices_length == 5 + 1);
    assert(g->state_length == 4);
    assert(g->vertices[0] == g->starting_vertex);
    assert(g->vertices[1] == v1);
    assert(g->vertices[2] == v2);
    assert(g->vertices[3] == v3);
    assert(g->vertices[4] == v4);
    assert(g->vertices[5] == v5);

    ptd_graph_destroy(g);
}

void test_basic_ptd_graph_edges() {
    struct ptd_graph *g = ptd_graph_create(4);
    struct ptd_vertex *v1 = ptd_vertex_create(g);
    struct ptd_vertex *v2 = ptd_vertex_create(g);
    struct ptd_vertex *v3 = ptd_vertex_create(g);
    struct ptd_vertex *v4 = ptd_vertex_create(g);

    ptd_graph_add_edge(g->starting_vertex, v1, 0.4);
    ptd_graph_add_edge(g->starting_vertex, v2, 0.6);
    ptd_graph_add_edge(v3, v2, 1.0);
    ptd_graph_add_edge(v3, v1, 2.0);
    ptd_graph_add_edge(v3, v4, 2.0);

    assert(g->starting_vertex->edges_length == 2);
    assert(g->starting_vertex->edges[0]->weight == 0.4);
    assert(g->starting_vertex->edges[1]->weight == 0.6);
    assert(v3->edges_length == 3);
    assert(v3->edges[0]->weight == 1.0);
    assert(v3->edges[0]->to == v2);
    assert(v3->edges[1]->weight == 2.0);
    assert(v3->edges[1]->to == v1);
    assert(v3->edges[2]->weight == 2.0);
    assert(v3->edges[2]->to == v4);


    ptd_graph_destroy(g);
}

void test_basic_ptd_graph_scc() {
    struct ptd_graph *graph = ptd_graph_create(4);

    struct ptd_vertex *S = graph->starting_vertex;
    fprintf(stderr, "S: %p\n", (void *) S);
    struct ptd_vertex *T = ptd_vertex_create(graph);
    fprintf(stderr, "T: %p\n", (void *) T);

    struct ptd_vertex *A = ptd_vertex_create(graph);
    fprintf(stderr, "A: %p\n", (void *) A);
    struct ptd_vertex *B = ptd_vertex_create(graph);
    fprintf(stderr, "B: %p\n", (void *) B);
    struct ptd_vertex *C = ptd_vertex_create(graph);
    fprintf(stderr, "C: %p\n", (void *) C);
    struct ptd_vertex *D = ptd_vertex_create(graph);
    fprintf(stderr, "D: %p\n", (void *) D);
    struct ptd_vertex *E = ptd_vertex_create(graph);
    fprintf(stderr, "E: %p\n", (void *) E);
    struct ptd_vertex *F = ptd_vertex_create(graph);
    fprintf(stderr, "F: %p\n", (void *) F);
    struct ptd_vertex *G = ptd_vertex_create(graph);
    fprintf(stderr, "G: %p\n", (void *) G);
    struct ptd_vertex *H = ptd_vertex_create(graph);
    fprintf(stderr, "H: %p\n", (void *) H);
    struct ptd_vertex *I = ptd_vertex_create(graph);
    fprintf(stderr, "I: %p\n", (void *) I);
    struct ptd_vertex *J = ptd_vertex_create(graph);
    fprintf(stderr, "J %p\n", (void *) J);
    struct ptd_vertex *K = ptd_vertex_create(graph);
    fprintf(stderr, "K %p\n", (void *) K);
    struct ptd_vertex *L = ptd_vertex_create(graph);
    fprintf(stderr, "L %p\n", (void *) L);


    ptd_graph_add_edge(A, B, 3);
    ptd_graph_add_edge(B, C, 4);
    ptd_graph_add_edge(C, A, 4);
    ptd_graph_add_edge(B, D, 2);
    ptd_graph_add_edge(D, E, 5);
    ptd_graph_add_edge(E, F, 5);
    ptd_graph_add_edge(E, I, 15);
    ptd_graph_add_edge(F, E, 1);
    ptd_graph_add_edge(F, G, 1);
    ptd_graph_add_edge(G, H, 1);
    ptd_graph_add_edge(H, F, 1);
    ptd_graph_add_edge(H, G, 1);

    ptd_graph_add_edge(H, T, 1);
    ptd_graph_add_edge(I, L, 1);
    ptd_graph_add_edge(L, T, 1);

    ptd_graph_add_edge(C, J, 1);

    ptd_graph_add_edge(J, K, 1);
    ptd_graph_add_edge(K, J, 1);

    ptd_graph_add_edge(K, E, 1);

    assert(A->edges[0]->to == B);
    assert(graph->vertices[0] == S);
    assert(graph->vertices[1] == T);
    assert(graph->vertices[2] == A);

    struct ptd_scc_graph *scc_graph = ptd_find_strongly_connected_components(
            graph
    );

    for (size_t i = 0; i < scc_graph->vertices_length; ++i) {
        fprintf(stderr, "\nComponents %zu %p:\n", i, (void *) scc_graph->vertices[i]);

        for (size_t j = 0; j < scc_graph->vertices[i]->internal_vertices_length; ++j) {
            fprintf(stderr, "\tVertex %zu %p\n", j, (void *) (scc_graph->vertices[i]->internal_vertices[j]));

            for (size_t k = 0; k < scc_graph->vertices[i]->internal_vertices_length; ++k) {
                /*fprintf(stderr, "\t\tTo internal vertex %zu %p exp %f\n", k,
                        (void *) scc_graph->vertices[i]->internal_vertices[k],
                        scc_graph->vertices[i]->internal_expected_visits[j][k]
                );*/
            }


        }

        for (size_t j = 0; j < scc_graph->vertices[i]->edges_length; ++j) {
            fprintf(stderr, "\tEdge to scc %zu %p\n", j, (void *) (scc_graph->vertices[i]->edges[j]->to));
        }
    }

    ptd_scc_graph_destroy(scc_graph);
    assert(A->edges[0]->to == B);
    assert(graph->vertices[0] == S);
    assert(graph->vertices[1] == T);
    assert(graph->vertices[2] == A);

    ptd_graph_destroy(graph);
}

void test_acyclic_expected_visits() {
    /*struct ptd_graph *graph = ptd_graph_create(4);

    struct ptd_vertex *S = graph->starting_vertex;
    fprintf(stderr, "S: %p\n", (void *) S);
    struct ptd_vertex *T = ptd_vertex_create(graph);
    fprintf(stderr, "T: %p\n", (void *) T);

    struct ptd_vertex *A = ptd_vertex_create(graph);
    fprintf(stderr, "A: %p\n", (void *) A);
    struct ptd_vertex *B = ptd_vertex_create(graph);
    fprintf(stderr, "B: %p\n", (void *) B);
    struct ptd_vertex *C = ptd_vertex_create(graph);
    fprintf(stderr, "C: %p\n", (void *) C);
    struct ptd_vertex *D = ptd_vertex_create(graph);
    fprintf(stderr, "D: %p\n", (void *) D);
    struct ptd_vertex *E = ptd_vertex_create(graph);
    fprintf(stderr, "E: %p\n", (void *) E);
    struct ptd_vertex *F = ptd_vertex_create(graph);
    fprintf(stderr, "F: %p\n", (void *) F);
    struct ptd_vertex *G = ptd_vertex_create(graph);
    fprintf(stderr, "G: %p\n", (void *) G);
    struct ptd_vertex *H = ptd_vertex_create(graph);
    fprintf(stderr, "H: %p\n", (void *) H);
    struct ptd_vertex *I = ptd_vertex_create(graph);
    fprintf(stderr, "I: %p\n", (void *) I);
    struct ptd_vertex *J = ptd_vertex_create(graph);
    fprintf(stderr, "J %p\n", (void *) J);
    struct ptd_vertex *K = ptd_vertex_create(graph);
    fprintf(stderr, "K %p\n", (void *) K);
    struct ptd_vertex *L = ptd_vertex_create(graph);
    fprintf(stderr, "L %p\n", (void *) L);

    ptd_graph_add_edge(S, A, 0.5);
    ptd_graph_add_edge(S, C, 0.5);
    ptd_graph_add_edge(A, C, 4);
    ptd_graph_add_edge(A, B, 2);
    ptd_graph_add_edge(B, D, 10);
    ptd_graph_add_edge(D, F, 5);
    ptd_graph_add_edge(C, F, 15);
    ptd_graph_add_edge(C, G, 15);
    ptd_graph_add_edge(G, F, 15);
    ptd_graph_add_edge(F, T, 1);

    double *exp = ptd_graph_acyclic_visit_probability(graph);
    double *exp2 = ptd_graph_expected_visits(graph);

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        fprintf(stderr, "%p: %f\n", (void *) graph->vertices[i], exp[i]);
        assert(abs(exp[i] - exp2[i]) < 0.00001);
    }

    free(exp2);

    assert(exp[0] == 1.0f); // S
    assert(exp[1] <= 1.1f); // T
    assert(exp[1] >= 0.9f); // T
    assert(exp[2] <= 0.51f); // A
    assert(exp[2] >= 0.49f); // A
    assert(exp[3] <= 0.17); // B
    assert(exp[3] >= 0.15); // B
    assert(exp[4] <= 0.85f); // C
    assert(exp[4] >= 0.80f); // C
    assert(exp[5] <= 0.17f); // D
    assert(exp[5] >= 0.16f); // D
    assert(exp[7] <= 1.1f); // F
    assert(exp[7] >= 0.9f); // F

    double *rewards = (double *) calloc(graph->vertices_length, sizeof(*rewards));

    for (size_t m = 0; m < graph->vertices_length; ++m) {
        rewards[m] = (m % 2) * 2;
    }

    double *new_rewards = ptd_graph_acyclic_moment_rewards(graph, rewards);

    struct ptd_scc_graph *scc_graph = ptd_find_strongly_connected_components(graph);

    double *new_rewards2 = ptd_graph_cyclic_moment_rewards(scc_graph, rewards);

    for (size_t m = 0; m < graph->vertices_length; ++m) {
        assert(abs(new_rewards2[m] - new_rewards[m]) < 0.001);
    }

    free(new_rewards2);
    free(new_rewards);
    free(rewards);

    free(exp);

    ptd_scc_graph_destroy(scc_graph);

    ptd_graph_destroy(graph);*/
}

void test_cyclic_expected_entry_visits() {
    /*struct ptd_graph *graph = ptd_graph_create(4);

    struct ptd_vertex *S = graph->starting_vertex;
    S->state[0] = 0;
    fprintf(stderr, "S: %p\n", (void *) S);
    struct ptd_vertex *T = ptd_vertex_create(graph);
    T->state[0] = 1;
    fprintf(stderr, "T: %p\n", (void *) T);

    struct ptd_vertex *A = ptd_vertex_create(graph);
    A->state[0] = 2;
    fprintf(stderr, "A: %p\n", (void *) A);
    struct ptd_vertex *B = ptd_vertex_create(graph);
    B->state[0] = 3;
    fprintf(stderr, "B: %p\n", (void *) B);
    struct ptd_vertex *C = ptd_vertex_create(graph);
    C->state[0] = 4;
    fprintf(stderr, "C: %p\n", (void *) C);
    struct ptd_vertex *D = ptd_vertex_create(graph);
    D->state[0] = 5;
    fprintf(stderr, "D: %p\n", (void *) D);
    struct ptd_vertex *E = ptd_vertex_create(graph);
    E->state[0] = 6;
    fprintf(stderr, "E: %p\n", (void *) E);
    struct ptd_vertex *F = ptd_vertex_create(graph);
    F->state[0] = 7;
    fprintf(stderr, "F: %p\n", (void *) F);
    struct ptd_vertex *G = ptd_vertex_create(graph);
    G->state[0] = 8;
    fprintf(stderr, "G: %p\n", (void *) G);
    struct ptd_vertex *H = ptd_vertex_create(graph);
    H->state[0] = 9;
    fprintf(stderr, "H: %p\n", (void *) H);
    struct ptd_vertex *I = ptd_vertex_create(graph);
    I->state[0] = 10;
    fprintf(stderr, "I: %p\n", (void *) I);
    struct ptd_vertex *J = ptd_vertex_create(graph);
    J->state[0] = 11;
    fprintf(stderr, "J %p\n", (void *) J);
    struct ptd_vertex *K = ptd_vertex_create(graph);
    K->state[0] = 12;
    fprintf(stderr, "K %p\n", (void *) K);
    struct ptd_vertex *L = ptd_vertex_create(graph);
    L->state[0] = 13;
    fprintf(stderr, "L %p\n", (void *) L);


    ptd_graph_add_edge(S, A, 1);

    ptd_graph_add_edge(A, B, 3);
    ptd_graph_add_edge(B, C, 4);
    ptd_graph_add_edge(C, A, 4);
    ptd_graph_add_edge(B, D, 2);
    ptd_graph_add_edge(D, E, 5);
    ptd_graph_add_edge(E, F, 5);
    ptd_graph_add_edge(E, I, 15);
    ptd_graph_add_edge(F, G, 1);
    ptd_graph_add_edge(G, H, 1);
    ptd_graph_add_edge(H, F, 1);
    ptd_graph_add_edge(H, G, 1);

    ptd_graph_add_edge(H, T, 1);
    ptd_graph_add_edge(I, L, 1);
    ptd_graph_add_edge(L, T, 1);

    ptd_graph_add_edge(C, J, 1);

    ptd_graph_add_edge(J, K, 1);
    ptd_graph_add_edge(K, J, 1);

    ptd_graph_add_edge(K, E, 1);

    assert(A->edges[0]->to == B);
    assert(graph->vertices[0] == S);
    assert(graph->vertices[1] == T);
    assert(graph->vertices[2] == A);

    struct ptd_scc_graph *scc_graph = ptd_find_strongly_connected_components(
            graph
    );

    double *probs = ptd_graph_cyclic_entry_probability(scc_graph);
    double *scc_probs = ptd_scc_graph_entry_probability(scc_graph);

    for (size_t i = 0; i < scc_graph->vertices_length; ++i) {
        fprintf(stderr, "\nComponents %zu %p PROB %f:\n", i, (void *) scc_graph->vertices[i], scc_probs[i]);

        for (size_t j = 0; j < scc_graph->vertices[i]->internal_vertices_length; ++j) {
            fprintf(stderr, "\tVertex %zu %p PROB %f\n", j, (void *) (scc_graph->vertices[i]->internal_vertices[j]),
                    probs[scc_graph->vertices[i]->internal_vertices[j]->index]);

            for (size_t k = 0; k < scc_graph->vertices[i]->external_vertices_length; ++k) {
                /*fprintf(stderr, "\t\tTo external vertex %zu %p exp %f\n", k,
                        (void *) scc_graph->vertices[i]->external_vertices[k],
                        scc_graph->vertices[i]->external_expected_visits[j][k]
                );*//*
            }
        }
    }

    free(scc_probs);
    ptd_scc_graph_destroy(scc_graph);
    assert(A->edges[0]->to == B);
    assert(graph->vertices[0] == S);
    assert(graph->vertices[1] == T);
    assert(graph->vertices[2] == A);

    assert(abs(probs[A->index] - 1.0) < 0.01);
    assert(abs(probs[J->index] - 0.28) < 0.01);
    assert(abs(probs[D->index] - (1 - 0.28)) < 0.01);
    assert(abs(probs[E->index] - 1) < 0.01);
    assert(abs(probs[F->index] - 0.25) < 0.01);
    assert(abs(probs[I->index] - 0.75) < 0.01);
    assert(abs(probs[L->index] - 0.75) < 0.01);
    assert(abs(probs[G->index] - 0) < 0.01);
    assert(abs(probs[B->index] - 0) < 0.01);
    assert(abs(probs[C->index] - 0) < 0.01);
    assert(abs(probs[T->index] - 1) < 0.01);
    free(probs);

    ptd_graph_add_edge(K, G, 1);

    scc_graph = ptd_find_strongly_connected_components(
            graph
    );

    scc_probs = ptd_scc_graph_entry_probability(scc_graph);

    assert(abs(scc_probs[0] - 1.0) < 0.01);
    assert(abs(scc_probs[1] - 0.35) < 0.01);
    assert(abs(scc_probs[2] - 0.64) < 0.01);
    assert(abs(scc_probs[3] - 0.64) < 0.01);
    assert(abs(scc_probs[4] - 0.85) < 0.01);
    assert(abs(scc_probs[5] - 0.28) < 0.01);
    assert(abs(scc_probs[6] - 0.71) < 0.01);
    assert(abs(scc_probs[7] - 1) < 0.01);
    assert(abs(scc_probs[8] - 1) < 0.01);

    double *exp_visits = ptd_graph_cyclic_expected_visits(scc_graph);
    double *exp_visits2 = ptd_graph_expected_visits(graph);

    for (size_t l = 0; l < graph->vertices_length; ++l) {
        fprintf(stderr, "Vertex %p, expected visits %f\n", (void *) graph->vertices[l],
                exp_visits[graph->vertices[l]->index]);
        assert(abs(exp_visits[l] - exp_visits2[l]) < 0.00001);
    }

    free(exp_visits2);

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


    double *rewards = (double *) calloc(scc_graph->graph->vertices_length, sizeof(*rewards));

    for (size_t m = 0; m < scc_graph->graph->vertices_length; ++m) {
        rewards[m] = 1;
    }


    double *new_rewards = ptd_graph_cyclic_moment_rewards(scc_graph, rewards);

    for (size_t m = 0; m < scc_graph->graph->vertices_length; ++m) {
        fprintf(stderr, "%p : %f\n", (void *) scc_graph->graph->vertices[m],
                new_rewards[m]);
    }

    free(new_rewards);
    free(rewards);
    free(exp_visits);

    free(scc_probs);

    rewards = (double *) calloc(scc_graph->graph->vertices_length, sizeof(*rewards));

    for (size_t m = 0; m < scc_graph->graph->vertices_length; ++m) {
        rewards[m] = (m % 2) * 2;
    }

    rewards[0] = 1;
    rewards[1] = 0;
    ptd_scc_graph_destroy(scc_graph);


    ptd_graph_add_edge(J, F, 1);
    ptd_graph_add_edge(J, E, 1);

    scc_graph = ptd_find_strongly_connected_components(
            graph
    );

    double *visits = ptd_graph_cyclic_expected_visits(scc_graph);
    double e = 0;

    for (size_t m = 0; m < scc_graph->graph->vertices_length; ++m) {
        if (graph->vertices[m]->edges_length != 0 && graph->starting_vertex != graph->vertices[m]) {
            e += visits[m] * rewards[m] / ptd_vertex_rate(graph->vertices[m]);
        }
    }

    fprintf(stderr, "e: %f\n", e);

    free(visits);

    ptd_graph_reward_transform(graph, rewards);

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        ptd_vertex *vertex = graph->vertices[i];

        fprintf(stderr, "Vertex %i %p:\n", vertex->state[0], (void *) vertex);

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            fprintf(stderr, "\tChild %i %p (%f)\n", vertex->edges[j]->to->state[0], (void *) vertex->edges[j]->to,
                    vertex->edges[j]->weight);
        }
    }

    ptd_scc_graph_destroy(scc_graph);
    scc_graph = ptd_find_strongly_connected_components(
            graph
    );
    visits = ptd_graph_cyclic_expected_visits(scc_graph);
    double e2 = 0;

    for (size_t m = 0; m < scc_graph->graph->vertices_length; ++m) {
        if (graph->vertices[m]->edges_length != 0 && graph->starting_vertex != graph->vertices[m]) {
            e2 += visits[m] * 1 / ptd_vertex_rate(graph->vertices[m]);
        }
    }

    fprintf(stderr, "e2: %f\n", e2);
    free(visits);

    free(rewards);
    ptd_scc_graph_destroy(scc_graph);

    ptd_graph_destroy(graph);
    assert(abs(e - e2) < 0.0001);*/
}

void test_expected() {
    struct ptd_graph *graph = ptd_graph_create(4);

    struct ptd_vertex *S = graph->starting_vertex;
    S->state[0] = 'S';
    fprintf(stderr, "S: %p\n", (void *) S);

    struct ptd_vertex *A = ptd_vertex_create(graph);
    A->state[0] = 'A';
    fprintf(stderr, "A: %p\n", (void *) A);
    struct ptd_vertex *B = ptd_vertex_create(graph);
    B->state[0] = 'B';
    fprintf(stderr, "B: %p\n", (void *) B);
    struct ptd_vertex *C = ptd_vertex_create(graph);
    C->state[0] = 'C';
    fprintf(stderr, "C: %p\n", (void *) C);
    struct ptd_vertex *D = ptd_vertex_create(graph);
    D->state[0] = 'D';
    fprintf(stderr, "D: %p\n", (void *) D);
    struct ptd_vertex *E = ptd_vertex_create(graph);
    E->state[0] = 'E';
    fprintf(stderr, "E: %p\n", (void *) E);


    ptd_graph_add_edge(S, A, 1);

    ptd_graph_add_edge(A, B, 3);
    ptd_graph_add_edge(B, C, 4);
    ptd_graph_add_edge(C, A, 4);
    ptd_graph_add_edge(B, D, 2);
    ptd_graph_add_edge(D, E, 5);


    double *wt = ptd_expected_waiting_time(graph, NULL);


    long double sum = 0;

    srand(1234);

    for (size_t j = 0; j < 1000000; ++j) {
        sum += ptd_random_sample(graph, NULL) / 1000000.0f;
    }

    assert(fabs(sum - wt[0])    < 0.01);


    free(wt);
    ptd_graph_destroy(graph);
}


void test_reward_transform() {
    struct ptd_graph *graph = ptd_graph_create(4);

    struct ptd_vertex *S = graph->starting_vertex;
    S->state[0] = 'S';
    fprintf(stderr, "S: %p\n", (void *) S);
    /*struct ptd_vertex *T = ptd_vertex_create(graph);
    T->state[0] = 1;
    fprintf(stderr, "T: %p\n", (void *) T);*/

    struct ptd_vertex *A = ptd_vertex_create(graph);
    A->state[0] = 'A';
    fprintf(stderr, "A: %p\n", (void *) A);
    struct ptd_vertex *B = ptd_vertex_create(graph);
    B->state[0] = 'B';
    fprintf(stderr, "B: %p\n", (void *) B);
    struct ptd_vertex *C = ptd_vertex_create(graph);
    C->state[0] = 'C';
    fprintf(stderr, "C: %p\n", (void *) C);
    struct ptd_vertex *D = ptd_vertex_create(graph);
    D->state[0] = 'D';
    fprintf(stderr, "D: %p\n", (void *) D);
    struct ptd_vertex *E = ptd_vertex_create(graph);
    E->state[0] = 'E';
    fprintf(stderr, "E: %p\n", (void *) E);
    struct ptd_vertex *F = ptd_vertex_create(graph);
    F->state[0] = 'F';
    /*fprintf(stderr, "F: %p\n", (void *) F);
    struct ptd_vertex *G = ptd_vertex_create(graph);
    G->state[0] = 8;
    fprintf(stderr, "G: %p\n", (void *) G);
    struct ptd_vertex *H = ptd_vertex_create(graph);
    H->state[0] = 9;
    fprintf(stderr, "H: %p\n", (void *) H);
    struct ptd_vertex *I = ptd_vertex_create(graph);
    I->state[0] = 10;
    fprintf(stderr, "I: %p\n", (void *) I);*/
    struct ptd_vertex *J = ptd_vertex_create(graph);
    J->state[0] = 'J';
    fprintf(stderr, "J %p\n", (void *) J);
    struct ptd_vertex *K = ptd_vertex_create(graph);
    K->state[0] = 'K';
    fprintf(stderr, "K %p\n", (void *) K);
    /*struct ptd_vertex *L = ptd_vertex_create(graph);
    L->state[0] = 13;
    fprintf(stderr, "L %p\n", (void *) L);*/


    ptd_graph_add_edge(S, A, 1);

    ptd_graph_add_edge(A, B, 3);
    ptd_graph_add_edge(B, C, 4);
    ptd_graph_add_edge(C, A, 4);
    ptd_graph_add_edge(B, D, 2);
    ptd_graph_add_edge(D, E, 5);
    /*ptd_graph_add_edge(E, F, 5);
    ptd_graph_add_edge(E, I, 15);
    ptd_graph_add_edge(F, G, 1);
    ptd_graph_add_edge(G, H, 1);
    ptd_graph_add_edge(H, F, 1);
    ptd_graph_add_edge(H, G, 1);

    ptd_graph_add_edge(H, T, 1);
    ptd_graph_add_edge(I, L, 1);
    ptd_graph_add_edge(L, T, 1);*/

    ptd_graph_add_edge(C, J, 1);

    ptd_graph_add_edge(J, K, 1);
    ptd_graph_add_edge(K, J, 1);

    ptd_graph_add_edge(K, E, 1);

    assert(A->edges[0]->to == B);
    assert(graph->vertices[0] == S);

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        fprintf(stderr, "Vertex %c:\n", (char) graph->vertices[i]->state[0]);

        for (size_t j = 0; j < graph->vertices[i]->edges_length; ++j) {
            fprintf(stderr, "\tTo %c (%f):\n", (char) graph->vertices[i]->edges[j]->to->state[0],
                    graph->vertices[i]->edges[j]->weight);
        }
    }

    double *rewards = (double *) calloc(graph->vertices_length, sizeof(double));

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        if (graph->vertices[i] == C || graph->vertices[i] == J) {
            rewards[i] = 0;
        } else {
            rewards[i] =  i + 1;
        }
    }

    double *wt = ptd_expected_waiting_time(graph, rewards);

    fprintf(stderr, "wt1: %f\n", wt[0]);

    struct ptd_graph *g = graph;
    graph = ptd_graph_reward_transform(graph, rewards);
    ptd_graph_destroy(g);

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        fprintf(stderr, "Vertex %c:\n", (char) graph->vertices[i]->state[0]);

        for (size_t j = 0; j < graph->vertices[i]->edges_length; ++j) {
            fprintf(stderr, "\tTo %c (%f):\n", (char) graph->vertices[i]->edges[j]->to->state[0],
                    graph->vertices[i]->edges[j]->weight);
        }
    }


    double *wt2 = ptd_expected_waiting_time(graph, NULL);
    for (size_t i = 0; i < graph->vertices_length; ++i) {
        //wt2 += e2[i];
    }

    fprintf(stderr, "wt2: %f\n", wt2[0]);

    assert(abs(wt[0] - wt2[0]) < 0.01);


    free(wt);
    free(wt2);
    ptd_graph_destroy(graph);
    free(rewards);
}

void test_reward_parameterized() {
    struct ptd_graph *graph = ptd_graph_create(4);

    struct ptd_vertex *S = graph->starting_vertex;
    S->state[0] = 'S';
    fprintf(stderr, "S: %p\n", (void *) S);
    /*struct ptd_vertex *T = ptd_vertex_create(graph);
    T->state[0] = 1;
    fprintf(stderr, "T: %p\n", (void *) T);*/

    struct ptd_vertex *A = ptd_vertex_create(graph);
    A->state[0] = 'A';
    fprintf(stderr, "A: %p\n", (void *) A);
    struct ptd_vertex *B = ptd_vertex_create(graph);
    B->state[0] = 'B';
    fprintf(stderr, "B: %p\n", (void *) B);
    struct ptd_vertex *C = ptd_vertex_create(graph);
    C->state[0] = 'C';
    fprintf(stderr, "C: %p\n", (void *) C);
    struct ptd_vertex *D = ptd_vertex_create(graph);
    D->state[0] = 'D';
    fprintf(stderr, "D: %p\n", (void *) D);
    struct ptd_vertex *E = ptd_vertex_create(graph);
    E->state[0] = 'E';
    fprintf(stderr, "E: %p\n", (void *) E);
    struct ptd_vertex *F = ptd_vertex_create(graph);
    F->state[0] = 'F';
    /*fprintf(stderr, "F: %p\n", (void *) F);
    struct ptd_vertex *G = ptd_vertex_create(graph);
    G->state[0] = 8;
    fprintf(stderr, "G: %p\n", (void *) G);
    struct ptd_vertex *H = ptd_vertex_create(graph);
    H->state[0] = 9;
    fprintf(stderr, "H: %p\n", (void *) H);
    struct ptd_vertex *I = ptd_vertex_create(graph);
    I->state[0] = 10;
    fprintf(stderr, "I: %p\n", (void *) I);*/
    struct ptd_vertex *J = ptd_vertex_create(graph);
    J->state[0] = 'J';
    fprintf(stderr, "J %p\n", (void *) J);
    struct ptd_vertex *K = ptd_vertex_create(graph);
    K->state[0] = 'K';
    fprintf(stderr, "K %p\n", (void *) K);
    /*struct ptd_vertex *L = ptd_vertex_create(graph);
    L->state[0] = 13;
    fprintf(stderr, "L %p\n", (void *) L);*/


    ptd_graph_add_edge(S, A, 1);

    ptd_graph_add_edge(A, B, 3);
    ptd_graph_add_edge(B, C, 4);
    ptd_graph_add_edge(C, A, 4);
    ptd_graph_add_edge(B, D, 2);
    ptd_graph_add_edge(D, E, 5);
    /*ptd_graph_add_edge(E, F, 5);
    ptd_graph_add_edge(E, I, 15);
    ptd_graph_add_edge(F, G, 1);
    ptd_graph_add_edge(G, H, 1);
    ptd_graph_add_edge(H, F, 1);
    ptd_graph_add_edge(H, G, 1);

    ptd_graph_add_edge(H, T, 1);
    ptd_graph_add_edge(I, L, 1);
    ptd_graph_add_edge(L, T, 1);*/

    ptd_graph_add_edge(C, J, 1);

    ptd_graph_add_edge(J, K, 1);
    ptd_graph_add_edge(K, J, 1);

    ptd_graph_add_edge(K, E, 1);

    assert(A->edges[0]->to == B);
    assert(graph->vertices[0] == S);

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        fprintf(stderr, "Vertex %c (%i):\n", (char) graph->vertices[i]->state[0], (int) graph->vertices[i]->index);

        for (size_t j = 0; j < graph->vertices[i]->edges_length; ++j) {
            fprintf(stderr, "\tTo %c (%f):\n", (char) graph->vertices[i]->edges[j]->to->state[0],
                    graph->vertices[i]->edges[j]->weight);
        }
    }

    double *rewards = (double *) calloc(graph->vertices_length, sizeof(double));

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        if (graph->vertices[i] == C || graph->vertices[i] == J) {
            rewards[i] = 0;
        } else {
            rewards[i] = i + 1;
        }
    }

    struct ptd_desc_reward_compute *comp = ptd_graph_ex_absorbation_time_comp_graph(graph);
    struct ptd_desc_reward_compute *comp2 = ptd_graph_ex_absorbation_time_comp_graph(graph);
    struct ptd_desc_reward_compute_parameterized *compp = ptd_graph_ex_absorbation_time_comp_graph_parameterized(graph);
    struct ptd_desc_reward_compute *comp3 = ptd_graph_build_ex_absorbation_time_comp_graph_parameterized(compp);
    struct ptd_desc_reward_compute *comp4 = ptd_graph_build_ex_absorbation_time_comp_graph_parameterized(compp);


    for (size_t i = 0; i < graph->vertices_length; ++i) {
        fprintf(stderr, "Vertex %c (%i):\n", (char) graph->vertices[i]->state[0], (int) graph->vertices[i]->index);

        for (size_t j = 0; j < graph->vertices[i]->edges_length; ++j) {
            fprintf(stderr, "\tTo %c (%f):\n", (char) graph->vertices[i]->edges[j]->to->state[0],
                    graph->vertices[i]->edges[j]->weight);
        }
    }


    assert(comp->length == comp2->length);

    for (size_t i = 0; i < comp->length; ++i) {
        assert(comp->commands[i].from == comp2->commands[i].from);
        assert(comp->commands[i].to == comp2->commands[i].to);
        assert(abs(comp->commands[i].multiplier - comp2->commands[i].multiplier) < 0.001);
    }


    assert(comp->length == comp3->length);

    for (size_t i = 0; i < comp->length; ++i) {
        assert(comp->commands[i].from == comp3->commands[i].from);
        assert(comp->commands[i].to == comp3->commands[i].to);
        assert(abs(comp->commands[i].multiplier - comp3->commands[i].multiplier) < 0.001);
    }


    assert(comp->length == comp4->length);

    for (size_t i = 0; i < comp->length; ++i) {
        assert(comp->commands[i].from == comp4->commands[i].from);
        assert(comp->commands[i].to == comp4->commands[i].to);
        assert(abs(comp->commands[i].multiplier - comp4->commands[i].multiplier) < 0.001);
    }

    A->edges[0]->weight = 8;
    B->edges[0]->weight = 9;
    B->edges[1]->weight = 10;
    C->edges[0]->weight = 11;
    C->edges[1]->weight = 12;
    struct ptd_desc_reward_compute *compA = ptd_graph_ex_absorbation_time_comp_graph(graph);
    struct ptd_desc_reward_compute *compA3 = ptd_graph_build_ex_absorbation_time_comp_graph_parameterized(compp);

    assert(compA->length == compA3->length);

    for (size_t i = 0; i < compA->length; ++i) {
        assert(compA->commands[i].from == compA3->commands[i].from);
        assert(compA->commands[i].to == compA3->commands[i].to);
        assert(abs(compA->commands[i].multiplier - compA3->commands[i].multiplier) < 0.001);
    }


    ptd_graph_destroy(graph);
    free(rewards);
}

void test_reward_parameterized2() {
    struct ptd_graph *graph = ptd_graph_create(4);

    struct ptd_vertex *S = graph->starting_vertex;
    S->state[0] = 'S';
    fprintf(stderr, "S: %p\n", (void *) S);
    /*struct ptd_vertex *T = ptd_vertex_create(graph);
    T->state[0] = 1;
    fprintf(stderr, "T: %p\n", (void *) T);*/

    struct ptd_vertex *A = ptd_vertex_create(graph);
    A->state[0] = 'A';
    fprintf(stderr, "A: %p\n", (void *) A);
    struct ptd_vertex *B = ptd_vertex_create(graph);
    B->state[0] = 'B';
    fprintf(stderr, "B: %p\n", (void *) B);
    struct ptd_vertex *C = ptd_vertex_create(graph);
    C->state[0] = 'C';
    fprintf(stderr, "C: %p\n", (void *) C);
    struct ptd_vertex *D = ptd_vertex_create(graph);
    D->state[0] = 'D';
    fprintf(stderr, "D: %p\n", (void *) D);
    struct ptd_vertex *E = ptd_vertex_create(graph);
    E->state[0] = 'E';
    fprintf(stderr, "E: %p\n", (void *) E);
    struct ptd_vertex *F = ptd_vertex_create(graph);
    F->state[0] = 'F';
    /*fprintf(stderr, "F: %p\n", (void *) F);
    struct ptd_vertex *G = ptd_vertex_create(graph);
    G->state[0] = 8;
    fprintf(stderr, "G: %p\n", (void *) G);
    struct ptd_vertex *H = ptd_vertex_create(graph);
    H->state[0] = 9;
    fprintf(stderr, "H: %p\n", (void *) H);
    struct ptd_vertex *I = ptd_vertex_create(graph);
    I->state[0] = 10;
    fprintf(stderr, "I: %p\n", (void *) I);*/
    struct ptd_vertex *J = ptd_vertex_create(graph);
    J->state[0] = 'J';
    fprintf(stderr, "J %p\n", (void *) J);
    struct ptd_vertex *K = ptd_vertex_create(graph);
    K->state[0] = 'K';
    fprintf(stderr, "K %p\n", (void *) K);
    /*struct ptd_vertex *L = ptd_vertex_create(graph);
    L->state[0] = 13;
    fprintf(stderr, "L %p\n", (void *) L);*/


    ptd_graph_add_edge_parameterized(S, A, 1, NULL);

    ptd_graph_add_edge(A, B, 3);
    ptd_graph_add_edge(B, C, 4);
    ptd_graph_add_edge(C, A, 4);
    ptd_graph_add_edge(B, D, 2);
    ptd_graph_add_edge(D, E, 5);
    /*ptd_graph_add_edge(E, F, 5);
    ptd_graph_add_edge(E, I, 15);
    ptd_graph_add_edge(F, G, 1);
    ptd_graph_add_edge(G, H, 1);
    ptd_graph_add_edge(H, F, 1);
    ptd_graph_add_edge(H, G, 1);

    ptd_graph_add_edge(H, T, 1);
    ptd_graph_add_edge(I, L, 1);
    ptd_graph_add_edge(L, T, 1);*/

    ptd_graph_add_edge(C, J, 1);

    ptd_graph_add_edge(J, K, 1);
    ptd_graph_add_edge(K, J, 1);

    ptd_graph_add_edge(K, E, 1);

    assert(A->edges[0]->to == B);
    assert(graph->vertices[0] == S);

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        fprintf(stderr, "Vertex %c (%i):\n", (char) graph->vertices[i]->state[0], (int) graph->vertices[i]->index);

        for (size_t j = 0; j < graph->vertices[i]->edges_length; ++j) {
            fprintf(stderr, "\tTo %c (%f):\n", (char) graph->vertices[i]->edges[j]->to->state[0],
                    graph->vertices[i]->edges[j]->weight);
        }
    }

    free(ptd_expected_waiting_time(graph, NULL));

    ptd_graph_destroy(graph);
}


void test_is_acyclic() {
    struct ptd_graph *graph = ptd_graph_create(4);

    struct ptd_vertex *S = graph->starting_vertex;
    fprintf(stderr, "S: %p\n", (void *) S);
    struct ptd_vertex *T = ptd_vertex_create(graph);
    fprintf(stderr, "T: %p\n", (void *) T);

    struct ptd_vertex *A = ptd_vertex_create(graph);
    fprintf(stderr, "A: %p\n", (void *) A);
    struct ptd_vertex *B = ptd_vertex_create(graph);
    fprintf(stderr, "B: %p\n", (void *) B);
    struct ptd_vertex *C = ptd_vertex_create(graph);
    fprintf(stderr, "C: %p\n", (void *) C);
    struct ptd_vertex *D = ptd_vertex_create(graph);
    fprintf(stderr, "D: %p\n", (void *) D);
    struct ptd_vertex *E = ptd_vertex_create(graph);
    fprintf(stderr, "E: %p\n", (void *) E);
    struct ptd_vertex *F = ptd_vertex_create(graph);
    fprintf(stderr, "F: %p\n", (void *) F);
    struct ptd_vertex *G = ptd_vertex_create(graph);
    fprintf(stderr, "G: %p\n", (void *) G);
    struct ptd_vertex *H = ptd_vertex_create(graph);
    fprintf(stderr, "H: %p\n", (void *) H);
    struct ptd_vertex *I = ptd_vertex_create(graph);
    fprintf(stderr, "I: %p\n", (void *) I);
    struct ptd_vertex *J = ptd_vertex_create(graph);
    fprintf(stderr, "J %p\n", (void *) J);
    struct ptd_vertex *K = ptd_vertex_create(graph);
    fprintf(stderr, "K %p\n", (void *) K);
    struct ptd_vertex *L = ptd_vertex_create(graph);
    fprintf(stderr, "L %p\n", (void *) L);

    ptd_graph_add_edge(S, A, 0.5);
    ptd_graph_add_edge(S, C, 0.5);
    ptd_graph_add_edge(A, C, 4);
    ptd_graph_add_edge(A, B, 2);
    ptd_graph_add_edge(B, D, 10);
    ptd_graph_add_edge(D, F, 5);
    ptd_graph_add_edge(C, F, 15);
    ptd_graph_add_edge(C, G, 15);
    ptd_graph_add_edge(G, F, 15);
    ptd_graph_add_edge(F, T, 1);

    assert(ptd_graph_is_acyclic(graph) == true);
    ptd_graph_add_edge(A, T, 1);
    assert(ptd_graph_is_acyclic(graph) == true);
    ptd_graph_add_edge(A, G, 1);
    assert(ptd_graph_is_acyclic(graph) == true);
    ptd_graph_add_edge(F, A, 1);
    assert(ptd_graph_is_acyclic(graph) == false);

    ptd_graph_destroy(graph);
}

void test_phase_type() {
    struct ptd_graph *graph = ptd_graph_create(4);

    struct ptd_vertex *S = graph->starting_vertex;
    fprintf(stderr, "S: %p\n", (void *) S);
    struct ptd_vertex *T = ptd_vertex_create(graph);
    fprintf(stderr, "T: %p\n", (void *) T);

    struct ptd_vertex *A = ptd_vertex_create(graph);
    fprintf(stderr, "A: %p\n", (void *) A);
    struct ptd_vertex *B = ptd_vertex_create(graph);
    fprintf(stderr, "B: %p\n", (void *) B);
    struct ptd_vertex *C = ptd_vertex_create(graph);
    fprintf(stderr, "C: %p\n", (void *) C);
    struct ptd_vertex *D = ptd_vertex_create(graph);
    fprintf(stderr, "D: %p\n", (void *) D);
    struct ptd_vertex *E = ptd_vertex_create(graph);
    fprintf(stderr, "E: %p\n", (void *) E);
    struct ptd_vertex *F = ptd_vertex_create(graph);
    fprintf(stderr, "F: %p\n", (void *) F);
    struct ptd_vertex *G = ptd_vertex_create(graph);
    fprintf(stderr, "G: %p\n", (void *) G);
    struct ptd_vertex *H = ptd_vertex_create(graph);
    fprintf(stderr, "H: %p\n", (void *) H);
    struct ptd_vertex *I = ptd_vertex_create(graph);
    fprintf(stderr, "I: %p\n", (void *) I);
    struct ptd_vertex *J = ptd_vertex_create(graph);
    fprintf(stderr, "J %p\n", (void *) J);
    struct ptd_vertex *K = ptd_vertex_create(graph);
    fprintf(stderr, "K %p\n", (void *) K);
    struct ptd_vertex *L = ptd_vertex_create(graph);
    fprintf(stderr, "L %p\n", (void *) L);

    ptd_graph_add_edge(S, A, 0.5);
    ptd_graph_add_edge(A, C, 4);
    ptd_graph_add_edge(A, B, 2);
    ptd_graph_add_edge(B, D, 10);
    ptd_graph_add_edge(D, F, 5);
    ptd_graph_add_edge(C, F, 15);
    ptd_graph_add_edge(C, G, 15);
    ptd_graph_add_edge(G, F, 15);
    ptd_graph_add_edge(F, T, 1);

    struct ptd_phase_type_distribution *ph = ptd_graph_as_phase_type_distribution(graph);

    for (size_t i = 0; i < ph->length; ++i) {
        for (size_t j = 0; j < ph->length; ++j) {
            fprintf(stderr, "%f ", ph->sub_intensity_matrix[i][j]);
        }

        fprintf(stderr, "\n");
    }

    assert(abs(ph->sub_intensity_matrix[0][0] - (-6)) < 0.0001);
    assert(abs(ph->sub_intensity_matrix[0][1] - 2) < 0.0001);
    assert(abs(ph->sub_intensity_matrix[0][2] - 4) < 0.0001);
    assert(abs(ph->sub_intensity_matrix[1][0] - 0) < 0.0001);
    assert(abs(ph->sub_intensity_matrix[2][0] - 0) < 0.0001);
    assert(abs(ph->sub_intensity_matrix[3][0] - 0) < 0.0001);
    assert(abs(ph->sub_intensity_matrix[4][0] - 0) < 0.0001);
    assert(abs(ph->sub_intensity_matrix[5][0] - 0) < 0.0001);

    ptd_phase_type_distribution_destroy(ph);
    ptd_graph_destroy(graph);
}

void test_kingman() {
    size_t n = 16;
    struct ptd_graph *kingman_graph = ptd_graph_create(n);
    struct ptd_avl_tree *avl_tree = ptd_avl_tree_create(n);
    int *istate = (int *) calloc(kingman_graph->state_length, sizeof(int));
    istate[0] = (int) n;
    ptd_graph_add_edge(kingman_graph->starting_vertex, ptd_vertex_create_state(kingman_graph, istate), 1);

    for (size_t k = 1; k < kingman_graph->vertices_length; k++) {
        struct ptd_vertex *vertex = kingman_graph->vertices[k];
        int *state = (int *) calloc(kingman_graph->state_length, sizeof(int));
        memcpy(state, vertex->state, kingman_graph->state_length * sizeof(int));

        for (size_t i = 0; i < kingman_graph->state_length; ++i) {
            for (size_t j = i; j < kingman_graph->state_length; ++j) {
                double weight;

                if (i == j) {
                    if (state[i] < 2) {
                        continue;
                    }

                    weight = state[i] * (state[i] - 1) / 2;
                } else {
                    if (state[i] < 1 || state[j] < 1) {
                        continue;
                    }

                    weight = state[i] * state[j];
                }

                state[i]--;
                state[j]--;
                state[(i + j + 2) - 1]++;

                struct ptd_vertex *child;
                struct ptd_avl_node *avl_node = ptd_avl_tree_find(avl_tree, state);

                if (avl_node == NULL) {
                    int *child_state = (int *) calloc(kingman_graph->state_length, sizeof(int));

                    memcpy(child_state, state, kingman_graph->state_length * sizeof(int));

                    child = ptd_vertex_create_state(kingman_graph, child_state);

                    ptd_avl_tree_find_or_insert(avl_tree, child_state, child);
                } else {
                    child = (struct ptd_vertex *) avl_node->entry;
                }

                state[i]++;
                state[j]++;
                state[(i + j + 2) - 1]--;

                ptd_graph_add_edge(vertex, child, weight);
            }
        }

        free(state);
    }

    fprintf(stderr, "%zu\n", kingman_graph->vertices_length);

    double *e = ptd_expected_waiting_time(kingman_graph, NULL);

    fprintf(stderr, "EXP %f\n", (float) e[0]);

    free(e);

    double *rewards = (double *) calloc(kingman_graph->vertices_length, sizeof(*rewards));

    for (size_t l = 0; l < kingman_graph->vertices_length; ++l) {
        rewards[l] = kingman_graph->vertices[l]->state[6] * kingman_graph->vertices[l]->state[3];
    }

    double *w;

    w = ptd_expected_waiting_time(kingman_graph, rewards);

    double a = w[0];
    free(w);

    struct ptd_graph *g;
    g = ptd_graph_reward_transform(kingman_graph, rewards);

    w = ptd_expected_waiting_time(g, NULL);
    double b = w[0];

    free(w);

    assert(abs(b - a) < 0.0001);


    free(rewards);

    ptd_avl_tree_destroy(avl_tree);
    ptd_graph_destroy(kingman_graph);
    ptd_graph_destroy(g);
}

struct conf {
    int locus1;
    int locus2;
    int population;
};

struct conf index_to_conf(int s, int i) {
    int d = s + 1;
    int p = i / (d * d);
    int a = (i - p * d * d) / d;
    int b = (i - p * d * d) % d;

    return ((struct conf) {.locus1 = a, .locus2=b, .population=p + 1});
}

int conf_to_index(int s, int locus1, int locus2, int population) {
    int d = s + 1;
    int i = (population - 1) * d * d + locus1 * d + locus2;

    return (i);
}

void test_2p2l() {
    int s = 3;
    int p = 2;

    size_t n = (size_t) p * (s + 1) * (s + 1);
    struct ptd_graph *graph = ptd_graph_create(n);
    struct ptd_avl_tree *avl_tree = ptd_avl_tree_create(n);

    struct ptd_vertex *first_vertex = ptd_vertex_create(graph);
    first_vertex->state[conf_to_index(s, 1, 1, 1)] = s;
    ptd_graph_add_edge(graph->starting_vertex, first_vertex, 1);

    for (int index = 1; index < graph->vertices_length; index++) {
        struct ptd_vertex *vertex = graph->vertices[index];
        int count = 0;

        for (size_t i = 0; i < n; ++i) {
            count += vertex->state[i];
        }

        if (count <= 1) {
            // Only one lineage, stop
            continue;
        }

        for (int i = 0; i < n; ++i) {
            struct conf conf_i = index_to_conf(s, i);

            // coalescence
            for (int j = i; j < n; ++j) {
                struct conf conf_j = index_to_conf(s, j);

                if (conf_i.population != conf_j.population) {
                    continue;
                }

                double rate;

                if (i == j) {
                    if (vertex->state[i] < 2) {
                        continue;
                    }

                    rate = vertex->state[i] * (vertex->state[i] - 1) / 2;
                } else {
                    if (vertex->state[i] < 1 || vertex->state[j] < 1) {
                        continue;
                    }

                    rate = vertex->state[i] * vertex->state[j];
                }

                int *child_state = (int *) calloc(n, sizeof(int));
                memcpy(child_state, vertex->state, sizeof(*child_state) * n);
                child_state[i] -= 1;
                child_state[j] -= 1;

                int k = conf_to_index(s, conf_i.locus1 + conf_j.locus1, conf_i.locus2 + conf_j.locus2,
                                      conf_i.population);
                child_state[k] += 1;

                struct ptd_vertex *child_vertex;

                struct ptd_avl_node *avl_node = ptd_avl_tree_find(avl_tree, child_state);

                if (avl_node == NULL) {
                    child_vertex = ptd_vertex_create_state(graph, child_state);

                    ptd_avl_tree_find_or_insert(avl_tree, child_state, child_vertex);
                } else {
                    free(child_state);

                    child_vertex = (struct ptd_vertex *) avl_node->entry;
                }

                ptd_graph_add_edge(vertex, child_vertex, rate);
            }

            if (vertex->state[i] > 0 && conf_i.locus1 > 0 && conf_i.locus2 > 0) {
                // recombination
                double rate = 1;

                int k = conf_to_index(s, conf_i.locus1, 0, conf_i.population);
                int l = conf_to_index(s, 0, conf_i.locus2, conf_i.population);
                int *child_state = (int *) calloc(n, sizeof(int));
                memcpy(child_state, vertex->state, sizeof(*child_state) * n);
                child_state[i] -= 1;
                child_state[k] += 1;
                child_state[l] += 1;

                struct ptd_vertex *child_vertex;

                struct ptd_avl_node *avl_node = ptd_avl_tree_find(avl_tree, child_state);

                if (avl_node == NULL) {
                    child_vertex = ptd_vertex_create_state(graph, child_state);

                    ptd_avl_tree_find_or_insert(avl_tree, child_state, child_vertex);
                } else {
                    free(child_state);

                    child_vertex = (struct ptd_vertex *) avl_node->entry;
                }

                ptd_graph_add_edge(vertex, child_vertex, rate);
            }

            if (vertex->state[i] > 0) {
                // migration
                double rate = 0.1;

                int *child_state = (int *) calloc(n, sizeof(int));
                memcpy(child_state, vertex->state, sizeof(*child_state) * n);

                int m;

                if (conf_i.population == 1) {
                    m = 2;
                } else {
                    m = 1;
                }

                int k = conf_to_index(s, conf_i.locus1, conf_i.locus2, m);
                child_state[i] -= 1;
                child_state[k] += 1;

                struct ptd_vertex *child_vertex;

                struct ptd_avl_node *avl_node = ptd_avl_tree_find(avl_tree, child_state);

                if (avl_node == NULL) {
                    child_vertex = ptd_vertex_create_state(graph, child_state);

                    ptd_avl_tree_find_or_insert(avl_tree, child_state, child_vertex);
                } else {
                    free(child_state);

                    child_vertex = (struct ptd_vertex *) avl_node->entry;
                }

                ptd_graph_add_edge(vertex, child_vertex, rate);
            }
        }
    }

    fprintf(stderr, "Done creating graph. It has %zu vertices\n", graph->vertices_length);

    double *rewards = (double *) calloc(graph->vertices_length, sizeof(*rewards));

    // NULL: (A), B, (C), D, E, F, G, H, J
    for (size_t i = 0; i < graph->vertices_length; ++i) {
        rewards[i] = graph->vertices[i]->state[2];
        /*if (i == 1 || i == 3 || i == 7) {
            rewards[i] = 0;
        } else {
            rewards[i] = 1;
        }*/
    }

    /*double *e = ptd_expected_waiting_time(graph, );
    double wt = 0;

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        wt += e[i] * rewards[i]; //graph->vertices[i]->state[2];
        graph->vertices[i]->state[0] = (int) ('A' + i);
    }

    fprintf(stderr, "wt1: %f\n", wt);


    ptd_graph_reward_transform(graph, rewards);


    double wt2 = 0;
    double *e2 = ptd_graph_expected_waiting_time(graph);

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        wt2 += e2[i];
    }

    fprintf(stderr, "wt2: %f\n", wt2);

    assert(abs(wt - wt2) < 0.01);


    free(e);
    free(e2);
    free(rewards);*/

    /*double *rewards = (double*) calloc(graph->vertices_length, sizeof(*rewards));
    int relevant_index1 = conf_to_index(s, 3, 1, 1);
    int relevant_index2 = conf_to_index(s, 3, 1, 2);


    for (size_t i = 0; i < graph->vertices_length; ++i) {
        rewards[i] = graph->vertices[i]->state[relevant_index1] != 0 || graph->vertices[i]->state[relevant_index2] != 0 ? 1 : 0;
    }

    ptd_graph_reward_transform(graph, rewards);

    double *e = ptd_graph_expected_waiting_time(graph);
    double wt = 0;

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        wt += e[i];
    }

    fprintf(stderr, "Exp waiting time: %f\n", wt);

    free(e);*/

    ptd_avl_tree_destroy(avl_tree);
    ptd_graph_destroy(graph);
}

void test_2p2lTIME() {
    int s = 5;
    int p = 2;

    size_t n = (size_t) p * (s + 1) * (s + 1);
    struct ptd_graph *graph = ptd_graph_create(n);
    struct ptd_avl_tree *avl_tree = ptd_avl_tree_create(n);

    struct ptd_vertex *first_vertex = ptd_vertex_create(graph);
    first_vertex->state[conf_to_index(s, 1, 1, 1)] = s;
    ptd_graph_add_edge(graph->starting_vertex, first_vertex, 1);

    for (int index = 1; index < graph->vertices_length; index++) {
        struct ptd_vertex *vertex = graph->vertices[index];
        int count = 0;

        for (size_t i = 0; i < n; ++i) {
            count += vertex->state[i];
        }

        if (count <= 1) {
            // Only one lineage, stop
            continue;
        }

        for (int i = 0; i < n; ++i) {
            struct conf conf_i = index_to_conf(s, i);

            // coalescence
            for (int j = i; j < n; ++j) {
                struct conf conf_j = index_to_conf(s, j);

                if (conf_i.population != conf_j.population) {
                    continue;
                }

                double rate;

                if (i == j) {
                    if (vertex->state[i] < 2) {
                        continue;
                    }

                    rate = vertex->state[i] * (vertex->state[i] - 1) / 2;
                } else {
                    if (vertex->state[i] < 1 || vertex->state[j] < 1) {
                        continue;
                    }

                    rate = vertex->state[i] * vertex->state[j];
                }

                int *child_state = (int *) calloc(n, sizeof(int));
                memcpy(child_state, vertex->state, sizeof(*child_state) * n);
                child_state[i] -= 1;
                child_state[j] -= 1;

                int k = conf_to_index(s, conf_i.locus1 + conf_j.locus1, conf_i.locus2 + conf_j.locus2,
                                      conf_i.population);
                child_state[k] += 1;

                struct ptd_vertex *child_vertex;

                struct ptd_avl_node *avl_node = ptd_avl_tree_find(avl_tree, child_state);

                if (avl_node == NULL) {
                    child_vertex = ptd_vertex_create_state(graph, child_state);

                    ptd_avl_tree_find_or_insert(avl_tree, child_state, child_vertex);
                } else {
                    free(child_state);

                    child_vertex = (struct ptd_vertex *) avl_node->entry;
                }

                ptd_graph_add_edge(vertex, child_vertex, rate);
            }

            if (vertex->state[i] > 0 && conf_i.locus1 > 0 && conf_i.locus2 > 0) {
                // recombination
                double rate = 1;

                int k = conf_to_index(s, conf_i.locus1, 0, conf_i.population);
                int l = conf_to_index(s, 0, conf_i.locus2, conf_i.population);
                int *child_state = (int *) calloc(n, sizeof(int));
                memcpy(child_state, vertex->state, sizeof(*child_state) * n);
                child_state[i] -= 1;
                child_state[k] += 1;
                child_state[l] += 1;

                struct ptd_vertex *child_vertex;

                struct ptd_avl_node *avl_node = ptd_avl_tree_find(avl_tree, child_state);

                if (avl_node == NULL) {
                    child_vertex = ptd_vertex_create_state(graph, child_state);

                    ptd_avl_tree_find_or_insert(avl_tree, child_state, child_vertex);
                } else {
                    free(child_state);

                    child_vertex = (struct ptd_vertex *) avl_node->entry;
                }

                ptd_graph_add_edge(vertex, child_vertex, rate);
            }

            if (vertex->state[i] > 0) {
                // migration
                double rate = 0.1;

                int *child_state = (int *) calloc(n, sizeof(int));
                memcpy(child_state, vertex->state, sizeof(*child_state) * n);

                int m;

                if (conf_i.population == 1) {
                    m = 2;
                } else {
                    m = 1;
                }

                int k = conf_to_index(s, conf_i.locus1, conf_i.locus2, m);
                child_state[i] -= 1;
                child_state[k] += 1;

                struct ptd_vertex *child_vertex;

                struct ptd_avl_node *avl_node = ptd_avl_tree_find(avl_tree, child_state);

                if (avl_node == NULL) {
                    child_vertex = ptd_vertex_create_state(graph, child_state);

                    ptd_avl_tree_find_or_insert(avl_tree, child_state, child_vertex);
                } else {
                    free(child_state);

                    child_vertex = (struct ptd_vertex *) avl_node->entry;
                }

                ptd_graph_add_edge(vertex, child_vertex, rate);
            }
        }
    }

    fprintf(stderr, "Done creating graph. It has %zu vertices\n", graph->vertices_length);

    double *rewards = (double *) calloc(graph->vertices_length, sizeof(*rewards));

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        rewards[i] = 1;
    }

    struct ptd_desc_reward_compute *comp =
            ptd_graph_ex_absorbation_time_comp_graph(graph);

    for (size_t j = 0; j < comp->length; ++j) {
        struct ptd_reward_increase command = comp->commands[j];


        fprintf(stderr, "%c <- %f + %f * %f = ", graph->vertices[command.from]->state[0],
                rewards[command.from],
                rewards[command.to],
                command.multiplier
        );
        rewards[command.from] += rewards[command.to] * command.multiplier;
        fprintf(stderr, "%f (%c)\n", rewards[command.from], graph->vertices[command.to]->state[0]);
    }

    free(comp->commands);
    free(comp);


    struct ptd_phase_type_distribution *ph = ptd_graph_as_phase_type_distribution(graph);

    for (size_t i = 0; i < ph->length; ++i) {
        for (size_t j = 0; j < ph->length; ++j) {
            fprintf(stderr, "%.4f,", ph->sub_intensity_matrix[i][j]);
        }

        fprintf(stderr, "\n");
    }

    ptd_phase_type_distribution_destroy(ph);


    for (size_t i = 0; i < graph->vertices_length; ++i) {
        fprintf(stderr, "REWARD %zu is %f\n", i, rewards[i]);
    }

    free(rewards);

    /*double *rewards = (double*) calloc(graph->vertices_length, sizeof(*rewards));
    int relevant_index1 = conf_to_index(s, 3, 1, 1);
    int relevant_index2 = conf_to_index(s, 3, 1, 2);


    for (size_t i = 0; i < graph->vertices_length; ++i) {
        rewards[i] = graph->vertices[i]->state[relevant_index1] != 0 || graph->vertices[i]->state[relevant_index2] != 0 ? 1 : 0;
    }

    ptd_graph_reward_transform(graph, rewards);

    double *e = ptd_graph_expected_waiting_time(graph);
    double wt = 0;

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        wt += e[i];
    }

    fprintf(stderr, "Exp waiting time: %f\n", wt);

    free(e);*/

    ptd_avl_tree_destroy(avl_tree);
    ptd_graph_destroy(graph);
}

void test_2p2lTIME2() {
    int s = 3;
    int p = 2;

    size_t n = (size_t) p * (s + 1) * (s + 1);
    struct ptd_graph *graph = ptd_graph_create(n);
    struct ptd_avl_tree *avl_tree = ptd_avl_tree_create(n);

    struct ptd_vertex *first_vertex = ptd_vertex_create(graph);
    first_vertex->state[conf_to_index(s, 1, 1, 1)] = s;
    ptd_graph_add_edge(graph->starting_vertex, first_vertex, 1);

    for (int index = 1; index < graph->vertices_length; index++) {
        struct ptd_vertex *vertex = graph->vertices[index];
        int count = 0;

        for (size_t i = 0; i < n; ++i) {
            count += vertex->state[i];
        }

        if (count <= 1) {
            // Only one lineage, stop
            continue;
        }

        for (int i = 0; i < n; ++i) {
            struct conf conf_i = index_to_conf(s, i);

            // coalescence
            for (int j = i; j < n; ++j) {
                struct conf conf_j = index_to_conf(s, j);

                if (conf_i.population != conf_j.population) {
                    continue;
                }

                double rate;

                if (i == j) {
                    if (vertex->state[i] < 2) {
                        continue;
                    }

                    rate = vertex->state[i] * (vertex->state[i] - 1) / 2;
                } else {
                    if (vertex->state[i] < 1 || vertex->state[j] < 1) {
                        continue;
                    }

                    rate = vertex->state[i] * vertex->state[j];
                }

                int *child_state = (int *) calloc(n, sizeof(int));
                memcpy(child_state, vertex->state, sizeof(*child_state) * n);
                child_state[i] -= 1;
                child_state[j] -= 1;

                int k = conf_to_index(s, conf_i.locus1 + conf_j.locus1, conf_i.locus2 + conf_j.locus2,
                                      conf_i.population);
                child_state[k] += 1;

                struct ptd_vertex *child_vertex;

                struct ptd_avl_node *avl_node = ptd_avl_tree_find(avl_tree, child_state);

                if (avl_node == NULL) {
                    child_vertex = ptd_vertex_create_state(graph, child_state);

                    ptd_avl_tree_find_or_insert(avl_tree, child_state, child_vertex);
                } else {
                    free(child_state);

                    child_vertex = (struct ptd_vertex *) avl_node->entry;
                }

                ptd_graph_add_edge(vertex, child_vertex, rate);
            }

            if (vertex->state[i] > 0 && conf_i.locus1 > 0 && conf_i.locus2 > 0) {
                // recombination
                double rate = 1;

                int k = conf_to_index(s, conf_i.locus1, 0, conf_i.population);
                int l = conf_to_index(s, 0, conf_i.locus2, conf_i.population);
                int *child_state = (int *) calloc(n, sizeof(int));
                memcpy(child_state, vertex->state, sizeof(*child_state) * n);
                child_state[i] -= 1;
                child_state[k] += 1;
                child_state[l] += 1;

                struct ptd_vertex *child_vertex;

                struct ptd_avl_node *avl_node = ptd_avl_tree_find(avl_tree, child_state);

                if (avl_node == NULL) {
                    child_vertex = ptd_vertex_create_state(graph, child_state);

                    ptd_avl_tree_find_or_insert(avl_tree, child_state, child_vertex);
                } else {
                    free(child_state);

                    child_vertex = (struct ptd_vertex *) avl_node->entry;
                }

                ptd_graph_add_edge(vertex, child_vertex, rate);
            }

            if (vertex->state[i] > 0) {
                // migration
                double rate = 0.1;

                int *child_state = (int *) calloc(n, sizeof(int));
                memcpy(child_state, vertex->state, sizeof(*child_state) * n);

                int m;

                if (conf_i.population == 1) {
                    m = 2;
                } else {
                    m = 1;
                }

                int k = conf_to_index(s, conf_i.locus1, conf_i.locus2, m);
                child_state[i] -= 1;
                child_state[k] += 1;

                struct ptd_vertex *child_vertex;

                struct ptd_avl_node *avl_node = ptd_avl_tree_find(avl_tree, child_state);

                if (avl_node == NULL) {
                    child_vertex = ptd_vertex_create_state(graph, child_state);

                    ptd_avl_tree_find_or_insert(avl_tree, child_state, child_vertex);
                } else {
                    free(child_state);

                    child_vertex = (struct ptd_vertex *) avl_node->entry;
                }

                ptd_graph_add_edge(vertex, child_vertex, rate);
            }
        }
    }

    fprintf(stderr, "Done creating graph. It has %zu vertices\n", graph->vertices_length);

    graph->parameterized = true;
    ptd_expected_waiting_time(graph, NULL);

    ptd_avl_tree_destroy(avl_tree);
    ptd_graph_destroy(graph);
}

void test_2pmigTIME() {
#define popsize 10

    struct mig {
        int pop1[(popsize) * (popsize)];
        int pop2[(popsize) * (popsize)];
    };

    size_t state_length = sizeof(mig) / sizeof(int);

    struct ptd_graph *graphA = ptd_graph_create(state_length);
    struct ptd_avl_tree *avl_treeA = ptd_avl_tree_create(state_length);

    struct mig *startA = (struct mig *) calloc(1, sizeof(*startA));
    startA->pop1[1] = popsize - 1;
    startA->pop1[0 + popsize] = popsize - 1;

    ptd_find_or_create_vertex(graphA, avl_treeA, (int *) startA);
    free(startA);

    for (size_t index = 1; index < graphA->vertices_length; ++index) {
        struct ptd_vertex *vertex = graphA->vertices[index];
        struct mig *state = (struct mig *) vertex->state;

        int count = 0;

        for (size_t i = 0; i < popsize; ++i) {
            for (size_t j = 0; j < popsize; ++j) {
                count += state->pop1[i + j * popsize];
                count += state->pop2[i + j * popsize];
            }
        }

        if (count == 0 || count == 1) {
            // Only one lineage left, absorb
            continue;
        }

        assert(count <= popsize * 2);

        for (size_t i = 0; i < popsize; ++i) {
            for (size_t j = 0; j < popsize; ++j) {
                if (state->pop1[i + j * popsize] == 0) {
                    continue;
                }
                for (size_t i2 = 0; i2 < popsize; ++i2) {
                    for (size_t j2 = 0; j2 < popsize; ++j2) {
                        if (i + j * popsize < i2 + j2 * popsize) {
                            continue;
                        }
                        if (state->pop1[i2 + j2 * popsize] == 0) {
                            continue;
                        }
                        if (i == i2 && j == j2) {
                            if (state->pop1[i2 + j2 * popsize] == 1) {
                                continue;
                            }
                        }

                        double weight;
                        if (i == i2 && j == j2) {
                            weight = state->pop1[i + j * popsize] * (state->pop1[i2 + j2 * popsize] - 1) / 2;
                        } else {
                            weight = state->pop1[i + j * popsize] * state->pop1[i2 + j2 * popsize];
                        }

                        size_t new_indexi = i + i2;
                        size_t new_indexj = j + j2;

                        struct mig *child_state = (struct mig *) malloc(sizeof(*child_state));
                        memcpy(child_state, state, sizeof(struct mig));
                        child_state->pop1[i + j * popsize]--;
                        child_state->pop1[i2 + j2 * popsize]--;
                        child_state->pop1[new_indexi + new_indexj * popsize]++;
                        struct ptd_vertex *child = ptd_find_or_create_vertex(graphA, avl_treeA, (int *) child_state);
                        ptd_graph_add_edge(vertex, child, weight);
                        free(child_state);
                    }
                }
            }
        }


        if (graphA->vertices_length % 100000 == 0) {
            fprintf(stderr, "AGraph has %zu vertices\n", graphA->vertices_length);
        }
    }

    clock_t start_time = time(NULL);

    struct ptd_graph *graph = ptd_graph_create(state_length);
    struct ptd_avl_tree *avl_tree = ptd_avl_tree_create(state_length);

    struct mig *start = (struct mig *) calloc(1, sizeof(*start));
    start->pop1[1] = popsize - 1;
    start->pop2[0 + popsize] = popsize - 1;

    struct ptd_vertex *first_vertex = ptd_vertex_create_state(graph, (int *) start);
    ptd_graph_add_edge(graph->starting_vertex, first_vertex, 1);
    struct mig child_state;

    for (size_t index = 1; index < graph->vertices_length; ++index) {
        struct ptd_vertex *vertex = graph->vertices[index];
        struct mig *state = (struct mig *) vertex->state;
        double migration_rate = 0.1;

        int count = 0;

        for (size_t i = 0; i < popsize; ++i) {
            for (size_t j = 0; j < popsize; ++j) {
                count += state->pop1[i + j * popsize];
                count += state->pop2[i + j * popsize];
            }
        }

        if (count == 0 || count == 1) {
            // Only one lineage left, absorb
            continue;
        }

        assert(count <= popsize * 2);

        for (size_t i = 0; i < popsize; ++i) {
            for (size_t j = 0; j < popsize; ++j) {
                if (state->pop1[i + j * popsize] == 0) {
                    continue;
                }
                memcpy(&child_state, state, sizeof(struct mig));
                child_state.pop1[i + j * popsize]--;
                child_state.pop2[i + j * popsize]++;
                struct ptd_vertex *child = ptd_find_or_create_vertex(graph, avl_tree, (int *) &child_state);
                ptd_graph_add_edge(vertex, child, migration_rate);
            }
        }

        for (size_t i = 0; i < popsize; ++i) {
            for (size_t j = 0; j < popsize; ++j) {
                if (state->pop2[i + j * popsize] == 0) {
                    continue;
                }
                memcpy(&child_state, state, sizeof(struct mig));
                child_state.pop2[i + j * popsize]--;
                child_state.pop1[i + j * popsize]++;
                struct ptd_vertex *child = ptd_find_or_create_vertex(graph, avl_tree, (int *) &child_state);
                ptd_graph_add_edge(vertex, child, migration_rate);
            }
        }


        for (size_t i = 0; i < popsize; ++i) {
            for (size_t j = 0; j < popsize; ++j) {
                if (state->pop1[i + j * popsize] == 0) {
                    continue;
                }
                for (size_t i2 = 0; i2 < popsize; ++i2) {
                    for (size_t j2 = 0; j2 < popsize; ++j2) {
                        if (i + j * popsize < i2 + j2 * popsize) {
                            continue;
                        }

                        if (state->pop1[i2 + j2 * popsize] == 0) {
                            continue;
                        }
                        if (i == i2 && j == j2) {
                            if (state->pop1[i2 + j2 * popsize] == 1) {
                                continue;
                            }
                        }

                        double weight;
                        if (i == i2 && j == j2) {
                            weight = state->pop1[i + j * popsize] * (state->pop1[i2 + j2 * popsize] - 1) / 2;
                        } else {
                            weight = state->pop1[i + j * popsize] * state->pop1[i2 + j2 * popsize];
                        }

                        size_t new_indexi = i + i2;
                        size_t new_indexj = j + j2;

                        memcpy(&child_state, state, sizeof(struct mig));
                        child_state.pop1[i + j * popsize]--;
                        child_state.pop1[i2 + j2 * popsize]--;
                        child_state.pop1[new_indexi + new_indexj * popsize]++;
                        struct ptd_vertex *child = ptd_find_or_create_vertex(graph, avl_tree, (int *) &child_state);
                        ptd_graph_add_edge(vertex, child, weight);
                    }
                }
            }
        }

        for (size_t i = 0; i < popsize; ++i) {
            for (size_t j = 0; j < popsize; ++j) {
                if (state->pop2[i + j * popsize] == 0) {
                    continue;
                }
                for (size_t i2 = 0; i2 < popsize; ++i2) {
                    for (size_t j2 = 0; j2 < popsize; ++j2) {
                        if (i + j * popsize < i2 + j2 * popsize) {
                            continue;
                        }

                        if (state->pop2[i2 + j2 * popsize] == 0) {
                            continue;
                        }
                        if (i == i2 && j == j2) {
                            if (state->pop2[i2 + j2 * popsize] == 1) {
                                continue;
                            }
                        }

                        size_t new_indexi = i + i2;
                        size_t new_indexj = j + j2;

                        double weight;
                        if (i == i2 && j == j2) {
                            weight = state->pop2[i + j * popsize] * (state->pop2[i2 + j2 * popsize] - 1) / 2;
                        } else {
                            weight = state->pop2[i + j * popsize] * state->pop2[i2 + j2 * popsize];
                        }

                        memcpy(&child_state, state, sizeof(struct mig));
                        child_state.pop2[i + j * popsize]--;
                        child_state.pop2[i2 + j2 * popsize]--;

                        child_state.pop2[new_indexi + new_indexj * popsize]++;
                        struct ptd_vertex *child = ptd_find_or_create_vertex(graph, avl_tree, (int *) &child_state);

                        ptd_graph_add_edge(vertex, child, weight);
                    }
                }
            }
        }

        if (graph->vertices_length % 10000 == 0) {
            fprintf(stderr, "Graph has %zu vertices\n", graph->vertices_length);
        }
    }

    size_t nedges = 0;

    for (size_t k = 0; k < graph->vertices_length; ++k) {
        nedges += graph->vertices[k]->edges_length;
    }

    clock_t end_time = time(NULL);

    fprintf(stderr, "Done creating graph. It has %zu vertices and %zu edges. Took %i seconds\n", graph->vertices_length,
            nedges, (int) (end_time - start_time));

    start_time = time(NULL);

    /* "Traditional" (Kern and Hey) method */

    /* Compute the stopping time for each vertex */
    struct ptd_probability_distribution_context *context =
            ptd_probability_distribution_context_create(graph, 1000);

    assert(context != 0);

    for (double t = 0; t <= 1.5; t += 1.0 / context->granularity) {
        ptd_probability_distribution_step(context);
    }

    long double *im_stop_probability = context->probability_at;

    end_time = time(NULL);

    fprintf(stderr, "Done computing im_stop_probability (cdf %f). Took %i seconds\n", context->cdf,
            (int) (end_time - start_time));

    start_time = time(NULL);

    /* Example, take singletons left doubletons right */
    double *singleton_doubleton_uncut =
            (double *) calloc(graph->vertices_length, sizeof(*singleton_doubleton_uncut));

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        struct ptd_vertex *vertex = graph->vertices[i];
        struct mig *state = (struct mig *) vertex->state;
        singleton_doubleton_uncut[i] = state->pop1[1 + 2 * popsize] + state->pop2[1 + 2 * popsize];
    }



    double *im_singdoub_waiting_time = ptd_expected_waiting_time(graph, singleton_doubleton_uncut);
    double uncut_im_singdoub_waiting_time = im_singdoub_waiting_time[0];


    end_time = time(NULL);

    fprintf(stderr, "Done expectation. Took %i seconds\n", (int) (end_time - start_time));

    start_time = time(NULL);

    long double tail_singdoub = 0;

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        tail_singdoub += im_stop_probability[i] * im_singdoub_waiting_time[i];
    }

    /* Expectation of X-Y = E[X] - E[Y]. Total expected singleton pop1 / doubleton pop2 rewarded "time" */
    long double isolated_singdoub = uncut_im_singdoub_waiting_time - tail_singdoub;

    /* How much in ancestral? */
    double *singleton_doubleton_A =
            (double *) calloc(graphA->vertices_length, sizeof(*singleton_doubleton_A));

    long double *starting_probabilities_A = (long double *) calloc(
            graphA->vertices_length,
            sizeof(*starting_probabilities_A)
    );

    for (size_t i = 1; i < graph->vertices_length; ++i) {
        struct ptd_vertex *im_vertex = graph->vertices[i];
        struct mig *IM_state = (struct mig*) im_vertex->state;
        struct mig *A_state = (struct mig*) calloc(1, sizeof(*A_state));

        for (size_t j = 0; j < (popsize) * (popsize); ++j) {
            A_state->pop1[j] = IM_state->pop1[j] + IM_state->pop2[j];
            A_state->pop2[j] = 0;
        }

        struct ptd_avl_node *node = ptd_avl_tree_find(avl_treeA, (int*)A_state);


        struct ptd_vertex *corresponding_A_vertex = (struct ptd_vertex *) node->entry;

        starting_probabilities_A[corresponding_A_vertex->index] += im_stop_probability[i];

        free(A_state);
    }

    for (size_t i = 0; i < graphA->vertices_length; ++i) {
        if (i == 0){
            singleton_doubleton_A[i] = 0;
        } else {
            struct mig *state = (struct mig *) graphA->vertices[i]->state;

            if (starting_probabilities_A[i] != 0) {
                ptd_graph_add_edge(graphA->starting_vertex, graphA->vertices[i], (double) starting_probabilities_A[i]);
            }

            singleton_doubleton_A[i] = state->pop1[1 + 2 * popsize];
        }
    }
    double *A_singdoub_waiting_time = ptd_expected_waiting_time(graphA, singleton_doubleton_A);
    double A_singdoub = A_singdoub_waiting_time[0];


    long double total_singdoub = isolated_singdoub + A_singdoub;

    fprintf(stderr, "Total singdoub = %f - %Lf + %f = %Lf\n",
            uncut_im_singdoub_waiting_time, tail_singdoub,
            A_singdoub,
            total_singdoub);

    free(im_singdoub_waiting_time);
    free(A_singdoub_waiting_time);
    free(singleton_doubleton_A);
    free(singleton_doubleton_uncut);
    free(starting_probabilities_A);



    ptd_probability_distribution_context_destroy(context);

    ptd_avl_tree_destroy(avl_tree);
    ptd_graph_destroy(graph);
    ptd_avl_tree_destroy(avl_treeA);
    ptd_graph_destroy(graphA);


    if (false) {
        for (size_t i1 = 0; i1 < 100000; i1 += 1) {
            /*expected += context->pdf / context->granularity * i1;

            fprintf(stderr, "%f\n", expected);
            ptd_probability_distribution_step(context);
            continue;*/

            if (i1 % 100 == 0) {
                fprintf(stderr, "CDF at time %f=%f (time spent: %i seconds)\n", (float) context->time, context->cdf,
                        (int) (time(NULL) - start_time));

                long double *prob_at_samples = (long double *) calloc(popsize * 2 + 2, sizeof(*prob_at_samples));

                for (size_t index = 1; index < graph->vertices_length; ++index) {
                    struct ptd_vertex *vertex = graph->vertices[index];
                    struct mig *state = (struct mig *) vertex->state;

                    size_t count = 0;

                    for (size_t i = 0; i < popsize; ++i) {
                        for (size_t j = 0; j < popsize; ++j) {
                            count += state->pop1[i + j * popsize];
                            count += state->pop2[i + j * popsize];
                        }
                    }

                    prob_at_samples[count] += context->probability_at[vertex->index];
                }

                fprintf(stderr, "Probability to stand at samples:\n");

                for (size_t k = 0; k < popsize * 2 + 1; ++k) {
                    fprintf(stderr, "\t%i\t%.7Lf\n", (int) k, prob_at_samples[k]);
                }

                fprintf(stderr, "\n");

                free(prob_at_samples);
            }

            ptd_probability_distribution_step(context);
        }
    }
}

void test_2pTIME() {
    int s = 7;
    int p = 2;

    size_t n = (size_t) p * (s + 1) * (s + 1);
    struct ptd_graph *graph = ptd_graph_create(n);
    struct ptd_avl_tree *avl_tree = ptd_avl_tree_create(n);

    struct ptd_vertex *first_vertex = ptd_vertex_create(graph);
    first_vertex->state[conf_to_index(s, 1, 1, 1)] = s;
    ptd_graph_add_edge(graph->starting_vertex, first_vertex, 1);

    for (int index = 1; index < graph->vertices_length; index++) {
        struct ptd_vertex *vertex = graph->vertices[index];
        int count = 0;

        for (size_t i = 0; i < n; ++i) {
            count += vertex->state[i];
        }

        if (count <= 1) {
            // Only one lineage, stop
            continue;
        }

        for (int i = 0; i < n; ++i) {
            struct conf conf_i = index_to_conf(s, i);

            // coalescence
            for (int j = i; j < n; ++j) {
                struct conf conf_j = index_to_conf(s, j);

                if (conf_i.population != conf_j.population) {
                    continue;
                }

                double rate;

                if (i == j) {
                    if (vertex->state[i] < 2) {
                        continue;
                    }

                    rate = vertex->state[i] * (vertex->state[i] - 1) / 2;
                } else {
                    if (vertex->state[i] < 1 || vertex->state[j] < 1) {
                        continue;
                    }

                    rate = vertex->state[i] * vertex->state[j];
                }

                int *child_state = (int *) calloc(n, sizeof(int));
                memcpy(child_state, vertex->state, sizeof(*child_state) * n);
                child_state[i] -= 1;
                child_state[j] -= 1;

                int k = conf_to_index(s, conf_i.locus1 + conf_j.locus1, conf_i.locus2 + conf_j.locus2,
                                      conf_i.population);
                child_state[k] += 1;

                struct ptd_vertex *child_vertex;

                struct ptd_avl_node *avl_node = ptd_avl_tree_find(avl_tree, child_state);

                if (avl_node == NULL) {
                    child_vertex = ptd_vertex_create_state(graph, child_state);

                    ptd_avl_tree_find_or_insert(avl_tree, child_state, child_vertex);
                } else {
                    free(child_state);

                    child_vertex = (struct ptd_vertex *) avl_node->entry;
                }

                ptd_graph_add_edge(vertex, child_vertex, rate);
            }

            if (vertex->state[i] > 0 && conf_i.locus1 > 0 && conf_i.locus2 > 0) {
                // recombination
                double rate = 1;

                int k = conf_to_index(s, conf_i.locus1, 0, conf_i.population);
                int l = conf_to_index(s, 0, conf_i.locus2, conf_i.population);
                int *child_state = (int *) calloc(n, sizeof(int));
                memcpy(child_state, vertex->state, sizeof(*child_state) * n);
                child_state[i] -= 1;
                child_state[k] += 1;
                child_state[l] += 1;

                struct ptd_vertex *child_vertex;

                struct ptd_avl_node *avl_node = ptd_avl_tree_find(avl_tree, child_state);

                if (avl_node == NULL) {
                    child_vertex = ptd_vertex_create_state(graph, child_state);

                    ptd_avl_tree_find_or_insert(avl_tree, child_state, child_vertex);
                } else {
                    free(child_state);

                    child_vertex = (struct ptd_vertex *) avl_node->entry;
                }

                ptd_graph_add_edge(vertex, child_vertex, rate);
            }

            if (vertex->state[i] > 0) {
                // migration
                double rate = 0.1;

                int *child_state = (int *) calloc(n, sizeof(int));
                memcpy(child_state, vertex->state, sizeof(*child_state) * n);

                int m;

                if (conf_i.population == 1) {
                    m = 2;
                } else {
                    m = 1;
                }

                int k = conf_to_index(s, conf_i.locus1, conf_i.locus2, m);
                child_state[i] -= 1;
                child_state[k] += 1;

                struct ptd_vertex *child_vertex;

                struct ptd_avl_node *avl_node = ptd_avl_tree_find(avl_tree, child_state);

                if (avl_node == NULL) {
                    child_vertex = ptd_vertex_create_state(graph, child_state);

                    ptd_avl_tree_find_or_insert(avl_tree, child_state, child_vertex);
                } else {
                    free(child_state);

                    child_vertex = (struct ptd_vertex *) avl_node->entry;
                }

                ptd_graph_add_edge(vertex, child_vertex, rate);
            }
        }
    }

    fprintf(stderr, "Done creating graph. It has %zu vertices\n", graph->vertices_length);

    double *r = (double *) calloc(graph->vertices_length, sizeof(*r));


    for (size_t i1 = 0; i1 < graph->vertices_length; ++i1) {
        r[i1] = 1000.0;
    }

    graph = ptd_graph_reward_transform(graph, r);

    struct ptd_dph_probability_distribution_context *context =
            ptd_dph_probability_distribution_context_create(graph);


    for (size_t i1 = 0; i1 < 10000; ++i1) {
        fprintf(stderr, "CDF %i=%f\n", (int) i1, context->cdf);

        ptd_dph_probability_distribution_step(context);
    }

    ptd_avl_tree_destroy(avl_tree);
    ptd_graph_destroy(graph);
}

void test_desc_exp() {
    struct ptd_graph *graph = ptd_graph_create(4);

    struct ptd_vertex *S = graph->starting_vertex;
    fprintf(stderr, "S: %p\n", (void *) S);
    struct ptd_vertex *T = ptd_vertex_create(graph);
    fprintf(stderr, "T: %p\n", (void *) T);

    struct ptd_vertex *C = ptd_vertex_create(graph);
    fprintf(stderr, "C: %p\n", (void *) C);
    struct ptd_vertex *B = ptd_vertex_create(graph);
    fprintf(stderr, "B: %p\n", (void *) B);
    struct ptd_vertex *D = ptd_vertex_create(graph);
    fprintf(stderr, "D: %p\n", (void *) D);
    struct ptd_vertex *A = ptd_vertex_create(graph);
    fprintf(stderr, "A: %p\n", (void *) A);
    struct ptd_vertex *F = ptd_vertex_create(graph);
    fprintf(stderr, "F: %p\n", (void *) F);


    S->state[0] = 'S';
    A->state[0] = 'A';
    B->state[0] = 'B';
    C->state[0] = 'C';
    D->state[0] = 'D';
    F->state[0] = 'F';
    T->state[0] = 'T';

    ptd_graph_add_edge(F, B, 2.3);
    ptd_graph_add_edge(F, D, 0.7);
    ptd_graph_add_edge(S, A, 1);
    ptd_graph_add_edge(A, B, 0.2);
    ptd_graph_add_edge(A, F, 0.2);
    ptd_graph_add_edge(A, D, 1.6);
    ptd_graph_add_edge(B, C, 5);
    ptd_graph_add_edge(C, A, 0.5);
    ptd_graph_add_edge(C, D, 1.5);
    ptd_graph_add_edge(D, T, 1);

    double *rewards = (double *) calloc(graph->vertices_length, sizeof(*rewards));

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        rewards[i] = 1;
    }

    struct ptd_desc_reward_compute *comp =
            ptd_graph_ex_absorbation_time_comp_graph(graph);

    for (size_t j = 0; j < comp->length; ++j) {
        struct ptd_reward_increase command = comp->commands[j];


        fprintf(stderr, "%c <- %f + %f * %f = ", graph->vertices[command.from]->state[0],
                rewards[command.from],
                rewards[command.to],
                command.multiplier
        );
        rewards[command.from] += rewards[command.to] * command.multiplier;
        fprintf(stderr, "%f (%c)\n", rewards[command.from], graph->vertices[command.to]->state[0]);
    }

    free(comp->commands);
    free(comp);


    for (size_t i = 0; i < graph->vertices_length; ++i) {
        fprintf(stderr, "REWARD %c is %f\n", graph->vertices[i]->state[0], rewards[i]);
    }

    assert(abs(1.68 - rewards[0]) < 0.1);
    assert(abs(0 - rewards[1]) < 0.1);
    assert(abs(1.67 - rewards[2]) < 0.1);
    assert(abs(1.87 - rewards[3]) < 0.1);
    assert(abs(1.00 - rewards[4]) < 0.1);
    assert(abs(1.68 - rewards[5]) < 0.1);
    assert(abs(2.00 - rewards[6]) < 0.1);

    double *rewardsinput = (double *) calloc(graph->vertices_length, sizeof(*rewardsinput));

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        rewardsinput[i] = 1;
    }

    double *n = ptd_expected_waiting_time(graph, rewardsinput);

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        assert(abs(rewards[i] - n[i]) < 0.01);
    }

    free(n);
    free(rewardsinput);
    free(rewards);
    ptd_graph_destroy(graph);
}

void test_desc_exp2() {
    struct ptd_graph *graph = ptd_graph_create(2);

    struct ptd_vertex *S = graph->starting_vertex;
    fprintf(stderr, "S: %p\n", (void *) S);
    struct ptd_vertex *T = ptd_vertex_create(graph);
    fprintf(stderr, "T: %p\n", (void *) T);

    struct ptd_vertex *C = ptd_vertex_create(graph);
    fprintf(stderr, "C: %p\n", (void *) C);
    struct ptd_vertex *B = ptd_vertex_create(graph);
    fprintf(stderr, "B: %p\n", (void *) B);
    struct ptd_vertex *D = ptd_vertex_create(graph);
    fprintf(stderr, "D: %p\n", (void *) D);
    struct ptd_vertex *A = ptd_vertex_create(graph);
    fprintf(stderr, "A: %p\n", (void *) A);
    struct ptd_vertex *F = ptd_vertex_create(graph);
    fprintf(stderr, "F: %p\n", (void *) F);


    S->state[0] = 'S';
    A->state[0] = 'A';
    B->state[0] = 'B';
    C->state[0] = 'C';
    D->state[0] = 'D';
    F->state[0] = 'F';
    T->state[0] = 'T';

    ptd_graph_add_edge(S, A, 1);
    ptd_graph_add_edge(A, B, 1);
    ptd_graph_add_edge(A, T, 2);
    ptd_graph_add_edge(B, C, 1);
    ptd_graph_add_edge(B, D, 2);
    ptd_graph_add_edge(B, A, 1);
    ptd_graph_add_edge(B, F, 4);
    ptd_graph_add_edge(C, B, 1);
    ptd_graph_add_edge(C, T, 4);
    ptd_graph_add_edge(D, F, 1);
    ptd_graph_add_edge(D, T, 2);
    ptd_graph_add_edge(F, D, 1);
    ptd_graph_add_edge(F, T, 4);

    double *rewards = (double *) calloc(graph->vertices_length, sizeof(*rewards));

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        rewards[i] = 1;
    }

    struct ptd_desc_reward_compute *comp =
            ptd_graph_ex_absorbation_time_comp_graph(graph);

    for (size_t j = 0; j < comp->length; ++j) {
        struct ptd_reward_increase command = comp->commands[j];


        fprintf(stderr, "%c <- %f + %f * %f = ", graph->vertices[command.from]->state[0],
                rewards[command.from],
                rewards[command.to],
                command.multiplier
        );
        rewards[command.from] += rewards[command.to] * command.multiplier;
        fprintf(stderr, "%f (%c)\n", rewards[command.from], graph->vertices[command.to]->state[0]);
    }

    free(comp->commands);
    free(comp);


    for (size_t i = 0; i < graph->vertices_length; ++i) {
        fprintf(stderr, "REWARD %c is %f\n", graph->vertices[i]->state[0], rewards[i]);
    }

    assert(abs(1.68 - rewards[0]) < 0.1);
    assert(abs(0 - rewards[1]) < 0.1);
    assert(abs(1.67 - rewards[2]) < 0.1);
    assert(abs(1.87 - rewards[3]) < 0.1);
    assert(abs(1.00 - rewards[4]) < 0.1);
    assert(abs(1.68 - rewards[5]) < 0.1);
    assert(abs(2.00 - rewards[6]) < 0.1);

    double *rewardsinput = (double *) calloc(graph->vertices_length, sizeof(*rewardsinput));

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        rewardsinput[i] = 1;
    }

    double *n = ptd_expected_waiting_time(graph, rewardsinput);

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        assert(abs(rewards[i] - n[i]) < 0.01);
    }

    free(n);
    free(rewardsinput);
    free(rewards);
    ptd_graph_destroy(graph);
}

/*
struct ptd_vertex *find_or_create_child(struct ptd_graph *graph, struct ptd_avl_tree *avl_tree, int *child_state) {
    struct ptd_vertex *child;
    struct ptd_avl_node *avl_node = ptd_avl_tree_find(avl_tree, child_state);

    if (avl_node == NULL) {
        child = ptd_vertex_create_state(graph, child_state);

        ptd_avl_tree_find_or_insert(avl_tree, child_state, child);
    } else {
        child = (struct ptd_vertex *) avl_node->entry;
    }

    return child;
}


void test_rabbit() {
    size_t state_size = 2;
    size_t starting_rabbits = 100;
    struct ptd_graph *graph = ptd_graph_create(state_size);
    struct ptd_avl_tree *avl_tree = ptd_avl_tree_create(state_size);
    int *initial_state = (int *) calloc(graph->state_length, sizeof(int));
    initial_state[0] = (int) starting_rabbits;
    ptd_graph_add_edge(
            graph->starting_vertex,
            ptd_vertex_create_state(graph, initial_state),
            1
    );

    for (size_t k = 1; k < graph->vertices_length; k++) {
        struct ptd_vertex *vertex = graph->vertices[k];
        int *state = (int *) calloc(graph->state_length, sizeof(int));
        memcpy(state, vertex->state, graph->state_length * sizeof(int));

        struct ptd_vertex *child;
        struct ptd_avl_node *avl_node = ptd_avl_tree_find(avl_tree, vertex->state);

        if (avl_node == NULL) {
            ptd_avl_tree_find_or_insert()
            int *child_state = (int *) calloc(graph->state_length, sizeof(int));
            memcpy(child_state, state, graph->state_length * sizeof(int));
            child_state[]

            child = ptd_vertex_create_state(graph, child_state);

            ptd_avl_tree_find_or_insert(avl_tree, child_state, child);
        } else {
            child = (struct ptd_vertex *) avl_node->entry;
        }

        state[i]++;
        state[j]++;
        state[(i + j + 2) - 1]--;

        ptd_graph_add_edge(vertex, child, weight);

        free(state);
    }

    ptd_avl_tree_destroy(avl_tree);
    ptd_graph_destroy(graph);
}*/

void test_random_sample() {
    struct ptd_graph *graph = ptd_graph_create(4);

    struct ptd_vertex *S = graph->starting_vertex;
    S->state[0] = 'S';
    struct ptd_vertex *T = ptd_vertex_create(graph);
    T->state[0] = 1;
    struct ptd_vertex *A = ptd_vertex_create(graph);
    A->state[0] = 'A';
    struct ptd_vertex *B = ptd_vertex_create(graph);
    B->state[0] = 'B';
    struct ptd_vertex *C = ptd_vertex_create(graph);
    C->state[0] = 'C';
    struct ptd_vertex *D = ptd_vertex_create(graph);
    D->state[0] = 'D';
    struct ptd_vertex *E = ptd_vertex_create(graph);
    E->state[0] = 'E';
    struct ptd_vertex *F = ptd_vertex_create(graph);
    F->state[0] = 'F';
    struct ptd_vertex *J = ptd_vertex_create(graph);
    J->state[0] = 'J';
    struct ptd_vertex *K = ptd_vertex_create(graph);
    K->state[0] = 'K';

    ptd_graph_add_edge(S, A, 1);

    ptd_graph_add_edge(A, B, 3);
    ptd_graph_add_edge(B, C, 4);
    ptd_graph_add_edge(C, A, 4);
    ptd_graph_add_edge(B, D, 2);
    ptd_graph_add_edge(D, E, 5);
    ptd_graph_add_edge(E, F, 5);
    ptd_graph_add_edge(F, E, 3);
    ptd_graph_add_edge(E, J, 3);

    ptd_graph_add_edge(C, J, 1);

    ptd_graph_add_edge(J, K, 1);
    ptd_graph_add_edge(K, J, 1);

    ptd_graph_add_edge(K, E, 1);
    ptd_graph_add_edge(K, T, 1);

    double *rewards = (double *) calloc(graph->vertices_length, sizeof(double));
    double *rewardsM = (double *) calloc(graph->vertices_length * 3, sizeof(double));

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        if (graph->vertices[i] == C || graph->vertices[i] == J) {
            rewards[i] = 0;
            rewardsM[i * 3] = 3;
            rewardsM[i * 3 + 1] = 0;
            rewardsM[i * 3 + 2] = 5;
        } else {
            rewards[i] = i + 1;
            rewardsM[i * 3] = 3;
            rewardsM[i * 3 + 1] = i + 1;
            rewardsM[i * 3 + 2] = i;
        }
    }

    double *wt = ptd_expected_waiting_time(graph, NULL);
    double e = wt[0];

    free(wt);

    long double sum = 0;

    srand(1234);

    for (size_t j = 0; j < 100000; ++j) {
        sum += ptd_random_sample(graph, NULL) / 100000.0f;
    }

    assert((sum - e) < 0.01);

    wt = ptd_expected_waiting_time(graph, rewards);
    e = wt[0];

    free(wt);

    sum = 0;

    srand(1234);

    for (size_t j = 0; j < 10000; ++j) {
        sum += ptd_random_sample(graph, rewards);
    }

    assert((sum / (10000.0) - e) < 0.01);

    sum = 0;

    srand(1234);

    for (size_t j = 0; j < 10000; ++j) {
        long double *res = ptd_mph_random_sample(graph, rewardsM, 3);
        sum += res[1];
        free(res);
    }

    assert((sum / (10000.0) - e) < 0.01);

    {
        struct ptd_graph *graphD = ptd_graph_create(4);

        struct ptd_vertex *SD = graphD->starting_vertex;
        S->state[0] = 'S';
        struct ptd_vertex *TD = ptd_vertex_create(graphD);
        T->state[0] = 1;
        struct ptd_vertex *AD = ptd_vertex_create(graphD);
        A->state[0] = 'A';
        struct ptd_vertex *BD = ptd_vertex_create(graphD);
        B->state[0] = 'B';

        ptd_graph_add_edge(SD, AD, 1);

        ptd_graph_add_edge(AD, BD, 0.8);
        ptd_graph_add_edge(BD, TD, 0.5);

        long double re1 = 0, re1v = 0;

        srand(1234);

        for (size_t j = 0; j < 10000; ++j) {
            long double r = ptd_dph_random_sample(graphD, NULL);
            assert(!isnan(r));
            re1 += r;
            re1v += r * r;
        }

        double *eres = ptd_expected_waiting_time(graphD, NULL);
        double *eres2 = ptd_expected_waiting_time(graphD, eres);

        double eexp = eres[0];
        double evar = eres2[0];

        free(eres);
        free(eres2);


        double *new_rewards = ptd_dph_normalize_graph(graphD);
        long double re2 = 0, re2v = 0;

        srand(1234);

        for (size_t j = 0; j < 10000; ++j) {
            long double r = ptd_dph_random_sample(graphD, new_rewards);
            re2 += r;
            re2v += r * r;
        }

        assert((re1 / 10000.0 - re2 / 10000.0) < 0.01);
        assert((re1 / 10000.0 - eexp) < 0.01);
        assert((re1v / 10000.0 - re2v / 10000.0) < 0.01);
        assert((re1v / 10000.0 - (evar * 2 - eexp)) < 0.01);

        free(new_rewards);
        ptd_graph_destroy(graphD);
    }

    {

        double *drw = (double *) calloc(4, sizeof(double));
        drw[0] = 0;
        drw[1] = 2;
        drw[2] = 3;
        drw[3] = 0;
        struct ptd_graph *graphD = ptd_graph_create(4);

        struct ptd_vertex *SD = graphD->starting_vertex;
        S->state[0] = 'S';
        struct ptd_vertex *BD = ptd_vertex_create(graphD);
        B->state[0] = 'B';
        struct ptd_vertex *AD = ptd_vertex_create(graphD);
        A->state[0] = 'A';
        struct ptd_vertex *TD = ptd_vertex_create(graphD);
        T->state[0] = 1;

        ptd_graph_add_edge(SD, AD, 1);

        ptd_graph_add_edge(AD, BD, 0.8);
        ptd_graph_add_edge(BD, TD, 0.5);

        long double re1 = 0, re1v = 0;

        srand(1234);

        for (size_t j = 0; j < 10000; ++j) {
            long double r = ptd_dph_random_sample(graphD, drw);
            assert(!isnan(r));
            re1 += r;
            re1v += r * r;
        }

        double *eres = ptd_expected_waiting_time(graphD, drw);
        double eexp = eres[0];
        eres[0] *= drw[0];
        eres[1] *= drw[1];
        eres[2] *= drw[2];
        eres[3] *= drw[3];

        double *eres2 = ptd_expected_waiting_time(graphD, eres);

        double evar = eres2[0];

        free(eres);
        free(eres2);


        double *new_rewards = ptd_dph_normalize_graph(graphD);
        new_rewards[0] *= drw[0];
        new_rewards[1] *= drw[1];
        new_rewards[2] *= drw[2];
        new_rewards[3] *= drw[3];
        long double re2 = 0, re2v = 0;

        srand(1234);

        for (size_t j = 0; j < 10000; ++j) {
            long double r = ptd_dph_random_sample(graphD, new_rewards);
            re2 += r;
            re2v += r * r;
        }

        assert((re1 / 10000.0 - re2 / 10000.0) < 0.01);
        assert((re1 / 10000.0 - eexp) < 0.01);
        assert((re1v / 10000.0 - re2v / 10000.0) < 0.01);
        assert((re1v / 10000.0 - (evar * 2 - eexp)) < 0.01);

        free(new_rewards);
        free(drw);
        ptd_graph_destroy(graphD);
    }

    ptd_graph_destroy(graph);

    free(rewardsM);
    free(rewards);
}

void test_pmf() {
    struct ptd_graph *graphD = ptd_graph_create(4);

    struct ptd_vertex *SD = graphD->starting_vertex;
    struct ptd_vertex *TD = ptd_vertex_create(graphD);
    struct ptd_vertex *TD2 = ptd_vertex_create(graphD);
    struct ptd_vertex *BD = ptd_vertex_create(graphD);
    struct ptd_vertex *AD = ptd_vertex_create(graphD);

    ptd_graph_add_edge(SD, AD, 0.75);
    ptd_graph_add_edge(SD, TD2, 0.25);

    ptd_graph_add_edge(AD, BD, 0.8);
    ptd_graph_add_edge(BD, TD, 0.5);
    ptd_graph_add_edge(BD, TD2, 0.2);
    struct ptd_dph_probability_distribution_context *context =
            ptd_dph_probability_distribution_context_create(graphD);

    assert(fabs(context->pmf - 0.25) <= 0.0001);
    assert(fabs(context->cdf - 0.25) <= 0.0001);
    assert(context->jumps == 0);

    ptd_dph_probability_distribution_step(context);

    assert(context->jumps == 1);
    assert(fabs(context->pmf - 0) <= 0.0001);
    assert(fabs(context->cdf - 0.25) <= 0.0001);

    ptd_dph_probability_distribution_step(context);

    assert(context->jumps == 2);
    assert(fabs(context->pmf - 0.75 * 0.8 * 0.7) <= 0.0001);
    assert(fabs(context->cdf - 0.25 - 0.75 * 0.8 * 0.7) <= 0.0001);

    ptd_dph_probability_distribution_step(context);

    assert(context->jumps == 3);

    double hits = 0;
    srand(1234);

    for (size_t i = 0; i < 100000; ++i) {
        if (ptd_dph_random_sample(graphD, NULL) == 3) {
            hits += 1.0 / 100000.0;
        }
    }

    assert(fabs(context->pmf - hits) <= 0.001);
    assert(fabs(context->cdf - 0.25 - 0.75 * 0.8 * 0.7 - hits) <= 0.001);

    ptd_dph_probability_distribution_step(context);

    assert(context->jumps == 4);

    double hits2 = 0;
    srand(1234);

    for (size_t i = 0; i < 100000; ++i) {
        if (ptd_dph_random_sample(graphD, NULL) == 4) {
            hits2 += 1.0 / 100000.0;
        }
    }

    assert(fabs(context->pmf - hits2) <= 0.001);
    assert(fabs(context->cdf - 0.25 - 0.75 * 0.8 * 0.7 - hits - hits2) <= 0.001);

    ptd_dph_probability_distribution_context_destroy(context);

    ptd_graph_destroy(graphD);
}

void test_rabbit() {
    size_t state_size = 2;
    int starting_rabbits = 2;

    struct ptd_graph *graph = ptd_graph_create(state_size);
    struct ptd_avl_tree *avl_tree = ptd_avl_tree_create(state_size);
    int *initial_state = (int*)calloc(graph->state_length, sizeof(*initial_state));
    int *child_state = (int*)calloc(graph->state_length, sizeof(*initial_state));
    initial_state[0] = starting_rabbits;
    ptd_graph_add_edge(
            graph->starting_vertex,
            ptd_vertex_create_state(graph, initial_state),
            1
    );

    for (size_t k = 1; k < graph->vertices_length; k++) {
        struct ptd_vertex *vertex = graph->vertices[k];
        int *state = vertex->state;

        if (state[0] > 0) {
            // Rabbit jump left to right
            memcpy(child_state, vertex->state, graph->state_length * sizeof(int));
            child_state[0] -= 1;
            child_state[1] += 1;
            ptd_graph_add_edge(
                    vertex,
                    ptd_find_or_create_vertex(graph, avl_tree, child_state),
                    1
            );

            memcpy(child_state, vertex->state, graph->state_length * sizeof(int));
            child_state[0] = 0;
            ptd_graph_add_edge(
                    vertex,
                    ptd_find_or_create_vertex(graph, avl_tree, child_state),
                    2
            );
        }

        if (state[1] > 0) {
            // Rabbit jump right to left
            memcpy(child_state, vertex->state, graph->state_length * sizeof(int));
            child_state[1] -= 1;
            child_state[0] += 1;
            ptd_graph_add_edge(
                    vertex,
                    ptd_find_or_create_vertex(graph, avl_tree, child_state),
                    1
            );

            memcpy(child_state, vertex->state, graph->state_length * sizeof(int));
            child_state[1] = 0;
            ptd_graph_add_edge(
                    vertex,
                    ptd_find_or_create_vertex(graph, avl_tree, child_state),
                    4
            );
        }
    }

    double *rw = (double*)calloc(graph->vertices_length, sizeof(*rw));

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        rw[i] = graph->vertices[i]->state[1];
    }

    struct ptd_graph *g = ptd_graph_reward_transform(graph, rw);

    free(child_state);
    assert(abs(ptd_defect(g) - 0.666) < 0.01);
}

int main(int argc, char **argv) {
    test_rabbit();
    // test_desc_exp();
    // test_desc_exp2();
    //test_2pmigTIME();
    //return 0;
    test_expected();
    test_reward_transform();
    //test_reward_parameterized();
    //test_reward_parameterized2();
    test_kingman();
    //return 0;
    test_pmf();
    test_random_sample();
    test_basic_graph();
    test_basic_ptd_graph();
    test_basic_ptd_graph_edges();
    test_basic_ptd_graph_scc();
    test_acyclic_expected_visits();
    test_cyclic_expected_entry_visits();
    test_is_acyclic();
    test_rabbit();
    return 0;
    //test_2p2lTIME();
    //test_2p2lTIME2();
    // test_phase_type();
    //test_2p2l();


    return 0;
}
