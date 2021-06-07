#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <numeric>
#include <string.h>
#include "../api/c/ptdalgorithms.h"


void assert(bool a) {
    if (!a) {
        exit(1);
    }
}


void test_it_can_insert_avl() {
    struct ptd_ph_graph *graph = ptd_ph_graph_create(4);
    assert(graph != NULL);
    ptd_avl_tree_t *tree = ptd_avl_tree_create(4);
    assert(tree != NULL);

    struct ptd_ph_vertex *S = graph->starting_vertex;
    fprintf(stderr, "S: %p\n", (void *) S);
    struct ptd_ph_vertex *T = ptd_ph_vertex_create(graph);
    fprintf(stderr, "T: %p\n", (void *) T);

    struct ptd_ph_vertex *A = ptd_ph_vertex_create(graph);
    fprintf(stderr, "A: %p\n", (void *) A);
    struct ptd_ph_vertex *B = ptd_ph_vertex_create(graph);
    fprintf(stderr, "B: %p\n", (void *) B);
    struct ptd_ph_vertex *C = ptd_ph_vertex_create(graph);
    fprintf(stderr, "C: %p\n", (void *) C);
    struct ptd_ph_vertex *D = ptd_ph_vertex_create(graph);
    fprintf(stderr, "D: %p\n", (void *) D);
    struct ptd_ph_vertex *E = ptd_ph_vertex_create(graph);
    fprintf(stderr, "E: %p\n", (void *) E);
    struct ptd_ph_vertex *F = ptd_ph_vertex_create(graph);
    fprintf(stderr, "F: %p\n", (void *) F);
    struct ptd_ph_vertex *G = ptd_ph_vertex_create(graph);
    fprintf(stderr, "G: %p\n", (void *) G);
    struct ptd_ph_vertex *H = ptd_ph_vertex_create(graph);
    fprintf(stderr, "H: %p\n", (void *) H);
    struct ptd_ph_vertex *I = ptd_ph_vertex_create(graph);
    fprintf(stderr, "I: %p\n", (void *) I);
    struct ptd_ph_vertex *J = ptd_ph_vertex_create(graph);
    fprintf(stderr, "J %p\n", (void *) J);
    struct ptd_ph_vertex *K = ptd_ph_vertex_create(graph);
    fprintf(stderr, "K %p\n", (void *) K);
    struct ptd_ph_vertex *L = ptd_ph_vertex_create(graph);
    fprintf(stderr, "L %p\n", (void *) L);

    S->state[0] = 2;
    A->state[0] = 4;
    B->state[0] = 3;
    C->state[0] = 3;
    C->state[1] = 2;
    D->state[0] = 3;
    D->state[1] = 1;
    E->state[0] = 1;

    assert(ptd_avl_tree_vertex_find(tree, S->state) == NULL);

    ptd_avl_tree_vertex_insert(tree, S->state, S);
    assert(ptd_avl_tree_vertex_find(tree, S->state) == S);

    assert(ptd_avl_tree_vertex_find(tree, A->state) == NULL);
    ptd_avl_tree_vertex_insert(tree, A->state, A);
    assert(ptd_avl_tree_vertex_find(tree, A->state) == A);

    ptd_avl_tree_vertex_insert(tree, B->state, B);
    assert(ptd_avl_tree_vertex_find(tree, B->state) == B);
    ptd_avl_tree_vertex_insert(tree, C->state, C);
    assert(ptd_avl_tree_vertex_find(tree, C->state) == C);
    assert(ptd_avl_tree_vertex_find(tree, D->state) == NULL);
    ptd_avl_tree_vertex_insert(tree, D->state, D);
    assert(ptd_avl_tree_vertex_find(tree, D->state) == D);
    ptd_avl_tree_vertex_insert(tree, E->state, E);
    assert(ptd_avl_tree_vertex_find(tree, E->state) == E);

    ptd_avl_tree_vertex_destroy(tree);
//destroy_vertices();
    ptd_ph_graph_destroy(graph);
}


void test_avl_is_balanced() {
    struct ptd_ph_graph *graph = ptd_ph_graph_create(4);
    ptd_avl_tree_t *tree = ptd_avl_tree_create(4);
    vector<struct ptd_ph_vertex *> vertices;

    for (size_t i = 0; i < 1024; ++i) {
        struct ptd_ph_vertex *vertex = ptd_ph_vertex_create(graph);
        assert(vertex != NULL);
        vertex->state[0] = rand();
        ptd_avl_tree_vertex_insert(tree, vertex->state, vertex);
        assert(ptd_avl_tree_vertex_find(tree, vertex->state) == vertex);
        vertices.push_back(vertex);
    }

    for (size_t k = 0; k < vertices.size(); ++k) {
        ptd_ph_vertex_destroy(vertices[k]);
    }

    fprintf(stderr, "Depth: %zu\n", ptd_avl_tree_max_depth(tree->root));
    assert(ptd_avl_tree_max_depth(tree->root) < 16);
    assert(ptd_avl_tree_max_depth(tree->root) > 2);
    ptd_avl_tree_vertex_destroy(tree);
    ptd_ph_graph_destroy(graph);
}

void test_avl_is_balanced_and_updated() {
    struct ptd_ph_graph *graph = ptd_ph_graph_create(4);
    ptd_avl_tree_t *tree = ptd_avl_tree_create(4);

    vector<struct ptd_ph_vertex *> vertices;

    for (size_t i = 0; i < 100; ++i) {
        for (size_t j = 0; j < 100; ++j) {
            struct ptd_ph_vertex *vertex = ptd_ph_vertex_create(graph);
            assert(vertex != NULL);
            vertex->state[0] = (int)j;
            assert(ptd_avl_tree_edge_insert_or_increment(tree, vertex->state, vertex, 1) == 0);
            if (i == 0) {
                assert(ptd_avl_tree_edge_find(tree, vertex->state)->to == vertex);
            } else {
                assert(ptd_avl_tree_edge_find(tree, vertex->state)->to != vertex);
            }

            assert(ptd_avl_tree_edge_find(tree, vertex->state)->weight == i + 1);
            vertices.push_back(vertex);
        }
    }

    for (size_t k = 0; k < vertices.size(); ++k) {
        ptd_ph_vertex_destroy(vertices[k]);
    }

    fprintf(stderr, "Depth: %zu\n", ptd_avl_tree_max_depth(tree->root));
    assert(ptd_avl_tree_max_depth(tree->root) < 16);
    assert(ptd_avl_tree_max_depth(tree->root) > 2);
    ptd_avl_tree_edge_destroy(tree);
    ptd_ph_graph_destroy(graph);
}


void test_avl_is_balanced_and_removed() {
    struct ptd_ph_graph *graph = ptd_ph_graph_create(4);
    ptd_avl_tree_t *tree = ptd_avl_tree_create(4);

    vector<struct ptd_ph_vertex *> vertices;


    for (size_t i = 0; i < 100; ++i) {
        struct ptd_ph_vertex *vertex = ptd_ph_vertex_create(graph);
        assert(vertex != NULL);
        vertex->state[0] = (int) i;
        assert(ptd_avl_tree_edge_insert_or_increment(tree, vertex->state, vertex, 1) == 0);
        assert(ptd_avl_tree_edge_find(tree, vertex->state)->to == vertex);
        assert(ptd_avl_tree_edge_find(tree, vertex->state)->weight == 1);
        ptd_avl_tree_edge_remove(tree, vertex->state);
        vertices.push_back(vertex);
    }

    struct ptd_ph_vertex *vertex[100];

    for (size_t i = 0; i < 100; ++i) {
        vertex[i] = ptd_ph_vertex_create(graph);
        assert(vertex[i] != NULL);
        vertex[i]->state[0] = (size_t) i;
        assert(ptd_avl_tree_edge_insert_or_increment(tree, vertex[i]->state, vertex[i], 1) == 0);
        assert(ptd_avl_tree_edge_find(tree, vertex[i]->state)->to == vertex[i]);
        assert(ptd_avl_tree_edge_find(tree, vertex[i]->state)->weight == 1);
        vertices.push_back(vertex[i]);
    }

    for (size_t i = 0; i < 100; ++i) {
        ptd_avl_tree_edge_remove(tree, vertex[i]->state);
        for (size_t j = i + 1; j < 100; ++j) {
            // TODO: Remove does not work!
            assert(ptd_avl_tree_edge_find(tree, vertex[j]->state) != NULL);
            assert(ptd_avl_tree_edge_find(tree, vertex[j]->state)->to == vertex[j]);
            assert(ptd_avl_tree_edge_find(tree, vertex[j]->state)->weight == 1);
        }
    }

    for (size_t k = 0; k < vertices.size(); ++k) {
        ptd_ph_vertex_destroy(vertices[k]);
    }

    fprintf(stderr, "Depth: %zu\n", ptd_avl_tree_max_depth(tree->root));
    assert(ptd_avl_tree_max_depth(tree->root) < 16);
    assert(ptd_avl_tree_max_depth(tree->root) > 2);
    ptd_avl_tree_edge_destroy(tree);
    ptd_ph_graph_destroy(graph);
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

void test_basic_ptd_ph_graph() {
    struct ptd_ph_graph *g = ptd_ph_graph_create(4);
    struct ptd_ph_vertex *v1 = ptd_ph_vertex_create(g);
    struct ptd_ph_vertex *v2 = ptd_ph_vertex_create(g);
    struct ptd_ph_vertex *v3 = ptd_ph_vertex_create(g);
    struct ptd_ph_vertex *v4 = ptd_ph_vertex_create(g);
    struct ptd_ph_vertex *v5 = ptd_ph_vertex_create(g);

    assert(g->vertices_length == 5 + 1);
    assert(g->state_length == 4);
    assert(g->vertices[0] == g->starting_vertex);
    assert(g->vertices[1] == v1);
    assert(g->vertices[2] == v2);
    assert(g->vertices[3] == v3);
    assert(g->vertices[4] == v4);
    assert(g->vertices[5] == v5);

    ptd_ph_graph_destroy(g);
}

void test_basic_ptd_ph_graph_edges() {
    struct ptd_ph_graph *g = ptd_ph_graph_create(4);
    struct ptd_ph_vertex *v1 = ptd_ph_vertex_create(g);
    struct ptd_ph_vertex *v2 = ptd_ph_vertex_create(g);
    struct ptd_ph_vertex *v3 = ptd_ph_vertex_create(g);
    struct ptd_ph_vertex *v4 = ptd_ph_vertex_create(g);

    ptd_ph_graph_add_edge(g->starting_vertex, v1, 0.4);
    ptd_ph_graph_add_edge(g->starting_vertex, v2, 0.6);
    ptd_ph_graph_add_edge(v3, v2, 1.0);
    ptd_ph_graph_add_edge(v3, v1, 2.0);
    ptd_ph_graph_add_edge(v3, v4, 2.0);

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


    ptd_ph_graph_destroy(g);
}

void test_basic_ptd_ph_graph_scc() {
    struct ptd_ph_graph *graph = ptd_ph_graph_create(4);

    struct ptd_ph_vertex *S = graph->starting_vertex;
    fprintf(stderr, "S: %p\n", (void *) S);
    struct ptd_ph_vertex *T = ptd_ph_vertex_create(graph);
    fprintf(stderr, "T: %p\n", (void *) T);

    struct ptd_ph_vertex *A = ptd_ph_vertex_create(graph);
    fprintf(stderr, "A: %p\n", (void *) A);
    struct ptd_ph_vertex *B = ptd_ph_vertex_create(graph);
    fprintf(stderr, "B: %p\n", (void *) B);
    struct ptd_ph_vertex *C = ptd_ph_vertex_create(graph);
    fprintf(stderr, "C: %p\n", (void *) C);
    struct ptd_ph_vertex *D = ptd_ph_vertex_create(graph);
    fprintf(stderr, "D: %p\n", (void *) D);
    struct ptd_ph_vertex *E = ptd_ph_vertex_create(graph);
    fprintf(stderr, "E: %p\n", (void *) E);
    struct ptd_ph_vertex *F = ptd_ph_vertex_create(graph);
    fprintf(stderr, "F: %p\n", (void *) F);
    struct ptd_ph_vertex *G = ptd_ph_vertex_create(graph);
    fprintf(stderr, "G: %p\n", (void *) G);
    struct ptd_ph_vertex *H = ptd_ph_vertex_create(graph);
    fprintf(stderr, "H: %p\n", (void *) H);
    struct ptd_ph_vertex *I = ptd_ph_vertex_create(graph);
    fprintf(stderr, "I: %p\n", (void *) I);
    struct ptd_ph_vertex *J = ptd_ph_vertex_create(graph);
    fprintf(stderr, "J %p\n", (void *) J);
    struct ptd_ph_vertex *K = ptd_ph_vertex_create(graph);
    fprintf(stderr, "K %p\n", (void *) K);
    struct ptd_ph_vertex *L = ptd_ph_vertex_create(graph);
    fprintf(stderr, "L %p\n", (void *) L);


    ptd_ph_graph_add_edge(A, B, 3);
    ptd_ph_graph_add_edge(B, C, 4);
    ptd_ph_graph_add_edge(C, A, 4);
    ptd_ph_graph_add_edge(B, D, 2);
    ptd_ph_graph_add_edge(D, E, 5);
    ptd_ph_graph_add_edge(E, F, 5);
    ptd_ph_graph_add_edge(E, I, 15);
    ptd_ph_graph_add_edge(F, E, 1);
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

    assert(A->edges[0]->to == B);
    assert(graph->vertices[0] == S);
    assert(graph->vertices[1] == T);
    assert(graph->vertices[2] == A);

    struct ptd_ph_scc_graph *scc_graph = ptd_ph_find_strongly_connected_components(
            graph
    );

    for (size_t i = 0; i < scc_graph->vertices_length; ++i) {
        fprintf(stderr, "\nComponents %zu %p:\n", i, (void *) scc_graph->vertices[i]);

        for (size_t j = 0; j < scc_graph->vertices[i]->internal_vertices_length; ++j) {
            fprintf(stderr, "\tVertex %zu %p\n", j, (void *) (scc_graph->vertices[i]->internal_vertices[j]));

            for (size_t k = 0; k < scc_graph->vertices[i]->internal_vertices_length; ++k) {
                fprintf(stderr, "\t\tTo internal vertex %zu %p exp %f\n", k,
                        (void *) scc_graph->vertices[i]->internal_vertices[k],
                        scc_graph->vertices[i]->internal_expected_visits[j][k]
                );
            }


            for (size_t k = 0; k < scc_graph->vertices[i]->external_vertices_length; ++k) {
                fprintf(stderr, "\t\tTo external vertex %zu %p exp %f\n", k,
                        (void *) scc_graph->vertices[i]->external_vertices[k],
                        scc_graph->vertices[i]->external_expected_visits[j][k]
                );
            }
        }

        for (size_t j = 0; j < scc_graph->vertices[i]->edges_length; ++j) {
            fprintf(stderr, "\tEdge to scc %zu %p\n", j, (void *) (scc_graph->vertices[i]->edges[j]->to));
        }

        for (size_t j = 0; j < scc_graph->vertices[i]->external_vertices_length; ++j) {
            fprintf(stderr, "\tEdge to vertex %zu %p\n", j, (void *) (scc_graph->vertices[i]->external_vertices[j]));
        }
    }

    ptd_ph_scc_graph_destroy(scc_graph);
    assert(A->edges[0]->to == B);
    assert(graph->vertices[0] == S);
    assert(graph->vertices[1] == T);
    assert(graph->vertices[2] == A);

    ptd_ph_graph_destroy(graph);
}

void test_acyclic_expected_visits() {
    struct ptd_ph_graph *graph = ptd_ph_graph_create(4);

    struct ptd_ph_vertex *S = graph->starting_vertex;
    fprintf(stderr, "S: %p\n", (void *) S);
    struct ptd_ph_vertex *T = ptd_ph_vertex_create(graph);
    fprintf(stderr, "T: %p\n", (void *) T);

    struct ptd_ph_vertex *A = ptd_ph_vertex_create(graph);
    fprintf(stderr, "A: %p\n", (void *) A);
    struct ptd_ph_vertex *B = ptd_ph_vertex_create(graph);
    fprintf(stderr, "B: %p\n", (void *) B);
    struct ptd_ph_vertex *C = ptd_ph_vertex_create(graph);
    fprintf(stderr, "C: %p\n", (void *) C);
    struct ptd_ph_vertex *D = ptd_ph_vertex_create(graph);
    fprintf(stderr, "D: %p\n", (void *) D);
    struct ptd_ph_vertex *E = ptd_ph_vertex_create(graph);
    fprintf(stderr, "E: %p\n", (void *) E);
    struct ptd_ph_vertex *F = ptd_ph_vertex_create(graph);
    fprintf(stderr, "F: %p\n", (void *) F);
    struct ptd_ph_vertex *G = ptd_ph_vertex_create(graph);
    fprintf(stderr, "G: %p\n", (void *) G);
    struct ptd_ph_vertex *H = ptd_ph_vertex_create(graph);
    fprintf(stderr, "H: %p\n", (void *) H);
    struct ptd_ph_vertex *I = ptd_ph_vertex_create(graph);
    fprintf(stderr, "I: %p\n", (void *) I);
    struct ptd_ph_vertex *J = ptd_ph_vertex_create(graph);
    fprintf(stderr, "J %p\n", (void *) J);
    struct ptd_ph_vertex *K = ptd_ph_vertex_create(graph);
    fprintf(stderr, "K %p\n", (void *) K);
    struct ptd_ph_vertex *L = ptd_ph_vertex_create(graph);
    fprintf(stderr, "L %p\n", (void *) L);

    ptd_ph_graph_add_edge(S, A, 0.5);
    ptd_ph_graph_add_edge(S, C, 0.5);
    ptd_ph_graph_add_edge(A, C, 4);
    ptd_ph_graph_add_edge(A, B, 2);
    ptd_ph_graph_add_edge(B, D, 10);
    ptd_ph_graph_add_edge(D, F, 5);
    ptd_ph_graph_add_edge(C, F, 15);
    ptd_ph_graph_add_edge(C, G, 15);
    ptd_ph_graph_add_edge(G, F, 15);
    ptd_ph_graph_add_edge(F, T, 1);

    double *exp = ptd_ph_graph_acyclic_visit_probability(graph);

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        fprintf(stderr, "%p: %f\n", (void *) graph->vertices[i], exp[i]);
    }

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

    double *new_rewards = ptd_ph_graph_acyclic_moment_rewards(graph, rewards);

    struct ptd_ph_scc_graph *scc_graph = ptd_ph_find_strongly_connected_components(graph);

    double *new_rewards2 = ptd_ph_graph_cyclic_moment_rewards(scc_graph, rewards);

    for (size_t m = 0; m < graph->vertices_length; ++m) {
        assert(abs(new_rewards2[m] - new_rewards[m]) < 0.001);
    }

    free(new_rewards2);
    free(new_rewards);
    free(rewards);

    free(exp);

    ptd_ph_scc_graph_destroy(scc_graph);

    ptd_ph_graph_destroy(graph);
}

void test_cyclic_expected_entry_visits() {
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

    assert(A->edges[0]->to == B);
    assert(graph->vertices[0] == S);
    assert(graph->vertices[1] == T);
    assert(graph->vertices[2] == A);

    struct ptd_ph_scc_graph *scc_graph = ptd_ph_find_strongly_connected_components(
            graph
    );

    double *probs = ptd_ph_graph_cyclic_entry_probability(scc_graph);
    double *scc_probs = ptd_ph_scc_graph_entry_probability(scc_graph);

    for (size_t i = 0; i < scc_graph->vertices_length; ++i) {
        fprintf(stderr, "\nComponents %zu %p PROB %f:\n", i, (void *) scc_graph->vertices[i], scc_probs[i]);

        for (size_t j = 0; j < scc_graph->vertices[i]->internal_vertices_length; ++j) {
            fprintf(stderr, "\tVertex %zu %p PROB %f\n", j, (void *) (scc_graph->vertices[i]->internal_vertices[j]),
                    probs[scc_graph->vertices[i]->internal_vertices[j]->index]);

            for (size_t k = 0; k < scc_graph->vertices[i]->external_vertices_length; ++k) {
                fprintf(stderr, "\t\tTo external vertex %zu %p exp %f\n", k,
                        (void *) scc_graph->vertices[i]->external_vertices[k],
                        scc_graph->vertices[i]->external_expected_visits[j][k]
                );
            }
        }
    }

    free(scc_probs);
    ptd_ph_scc_graph_destroy(scc_graph);
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

    ptd_ph_graph_add_edge(K, G, 1);

    scc_graph = ptd_ph_find_strongly_connected_components(
            graph
    );

    scc_probs = ptd_ph_scc_graph_entry_probability(scc_graph);

    assert(abs(scc_probs[0] - 1.0) < 0.01);
    assert(abs(scc_probs[1] - 0.35) < 0.01);
    assert(abs(scc_probs[2] - 0.64) < 0.01);
    assert(abs(scc_probs[3] - 0.64) < 0.01);
    assert(abs(scc_probs[4] - 0.85) < 0.01);
    assert(abs(scc_probs[5] - 0.28) < 0.01);
    assert(abs(scc_probs[6] - 0.71) < 0.01);
    assert(abs(scc_probs[7] - 1) < 0.01);
    assert(abs(scc_probs[8] - 1) < 0.01);

    double *exp_visits = ptd_ph_graph_cyclic_expected_visits(scc_graph);

    for (size_t l = 0; l < graph->vertices_length; ++l) {
        fprintf(stderr, "Vertex %p, expected visits %f\n", (void *) graph->vertices[l],
                exp_visits[graph->vertices[l]->index]);
    }

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


    double *new_rewards = ptd_ph_graph_cyclic_moment_rewards(scc_graph, rewards);

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
    ptd_ph_scc_graph_destroy(scc_graph);


    ptd_ph_graph_add_edge(J, F, 1);
    ptd_ph_graph_add_edge(J, E, 1);

    scc_graph = ptd_ph_find_strongly_connected_components(
            graph
    );

    double *visits = ptd_ph_graph_cyclic_expected_visits(scc_graph);
    double e = 0;

    for (size_t m = 0; m < scc_graph->graph->vertices_length; ++m) {
        if (graph->vertices[m]->edges_length != 0 && graph->starting_vertex != graph->vertices[m]) {
            e += visits[m] * rewards[m] / ptd_ph_vertex_rate(graph->vertices[m]);
        }
    }

    fprintf(stderr, "e: %f\n", e);

    free(visits);

    ptd_ph_graph_reward_transform(graph, rewards);

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        ptd_ph_vertex *vertex = graph->vertices[i];

        fprintf(stderr, "Vertex %i %p:\n", vertex->state[0], (void *) vertex);

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            fprintf(stderr, "\tChild %i %p (%f)\n", vertex->edges[j]->to->state[0],  (void *) vertex->edges[j]->to, vertex->edges[j]->weight);
        }
    }

    ptd_ph_scc_graph_destroy(scc_graph);
    scc_graph = ptd_ph_find_strongly_connected_components(
            graph
    );
    visits = ptd_ph_graph_cyclic_expected_visits(scc_graph);
    double e2 = 0;

    for (size_t m = 0; m < scc_graph->graph->vertices_length; ++m) {
        if (graph->vertices[m]->edges_length != 0 && graph->starting_vertex != graph->vertices[m]) {
            e2 += visits[m] * 1 / ptd_ph_vertex_rate(graph->vertices[m]);
        }
    }

    fprintf(stderr, "e2: %f\n", e2);
    free(visits);

    free(rewards);
    ptd_ph_scc_graph_destroy(scc_graph);

    ptd_ph_graph_destroy(graph);
    assert(abs(e-e2) < 0.0001);
}

int main(int argc, char **argv) {
    test_basic_graph();
    test_basic_ptd_ph_graph();
    test_basic_ptd_ph_graph_edges();
    test_basic_ptd_ph_graph_scc();
    test_acyclic_expected_visits();
    test_cyclic_expected_entry_visits();
}
