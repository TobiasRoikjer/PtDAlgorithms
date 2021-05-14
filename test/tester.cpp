#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <numeric>
#include <string.h>
#include "../api/c/ptdalgorithms.h"
#include "../src/c/io.h"

#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>

static gsl_matrix *
matrix_invert(gsl_matrix *matrix, size_t size) {
    gsl_matrix *inverse = gsl_matrix_alloc(size, size);

    int sign;
    gsl_permutation *p = gsl_permutation_alloc(size);

    gsl_linalg_LU_decomp(matrix, p, &sign);

    gsl_linalg_LU_invert(matrix, p, inverse);

    gsl_permutation_free(p);

    return inverse;
}

void assert(bool a) {
    if (!a) {
        exit(1);
    }
}

ptd_vertex_t *S, *A, *B, *C, *D, *E, *F, *G, *H, *I, *J, *K, *L, *T;

void create_vertices(ptd_graph_t *graph) {
    S = graph->start_vertex;
    T = ptd_vertex_create(graph);

    A = ptd_vertex_create(graph);
    B = ptd_vertex_create(graph);
    C = ptd_vertex_create(graph);
    D = ptd_vertex_create(graph);
    E = ptd_vertex_create(graph);
    F = ptd_vertex_create(graph);
    G = ptd_vertex_create(graph);
    H = ptd_vertex_create(graph);
    I = ptd_vertex_create(graph);
    J = ptd_vertex_create(graph);
    K = ptd_vertex_create(graph);
    L = ptd_vertex_create(graph);

    A->state[0] = 1;
    B->state[0] = 2;
    C->state[0] = 3;
    D->state[0] = 4;
    E->state[0] = 5;
    F->state[0] = 6;
    G->state[0] = 7;
    H->state[0] = 8;
    I->state[0] = 9;
    J->state[0] = 10;
    K->state[0] = 11;
    L->state[0] = 12;
}

void destroy_vertices() {
    ptd_vertex_destroy(A);
    ptd_vertex_destroy(B);
    ptd_vertex_destroy(C);
    ptd_vertex_destroy(D);
    ptd_vertex_destroy(E);
    ptd_vertex_destroy(F);
    ptd_vertex_destroy(G);
    ptd_vertex_destroy(H);
    ptd_vertex_destroy(I);
    ptd_vertex_destroy(L);

    A = NULL;
    B = NULL;
    C = NULL;
    D = NULL;
    E = NULL;
    F = NULL;
    G = NULL;
    H = NULL;
    I = NULL;
    L = NULL;
}

char *vertex_name(ptd_vertex_t *vertex) {
    if (vertex == S) {
        return "S";
    }

    if (vertex == T) {
        return "T";
    }

    if (vertex == A) {
        return "A";
    }

    if (vertex == B) {
        return "B";
    }

    if (vertex == C) {
        return "C";
    }

    if (vertex == D) {
        return "D";
    }

    if (vertex == E) {
        return "E";
    }

    if (vertex == F) {
        return "F";
    }

    if (vertex == G) {
        return "G";
    }

    if (vertex == H) {
        return "H";
    }

    if (vertex == I) {
        return "I";
    }

    if (vertex == K) {
        return "K";
    }

    if (vertex == J) {
        return "J";
    }

    if (vertex == L) {
        return "L";
    }

    return "Unknown";
}

size_t vertex_len;

void print_vertex_state(ptd_vertex_t *vertex) {
    fprintf(stderr, "(");

    for (size_t i = 0; i < vertex_len; ++i) {
        fprintf(stderr, "%zu, ", vertex->state[i]);
    }

    fprintf(stderr, ")");
}

int print_func(ptd_vertex_t *vertex) {
    fprintf(stderr, "Vertex %s ", vertex_name(vertex));

    print_vertex_state(vertex);

    fprintf(stderr, ")\n");

    for (size_t i = 0; i < vertex->edges_length; ++i) {
        fprintf(stderr, "\t %s %Lf  (", vertex_name(vertex->edges[i].to), vertex->edges[i].weight);
        print_vertex_state(vertex->edges[i].to);
        fprintf(stderr, ")\n");
    }

    fprintf(stderr, "\n");

    return 0;
}

void draw_graph(ptd_graph_t *graph) {
    vertex_len = graph->state_length;
    ptd_visit_vertices(graph, print_func, true);
}

bool keep_all(ptd_vertex_t *vertex) {
    return true;
}

double reward_identity(ptd_vertex_t *vertex) {
    return 1;
}

void test_can_find_scc() {
    ptd_graph_t *graph = ptd_graph_create(4);

    create_vertices(graph);

    ptd_add_edge(S, H, 1);

    ptd_add_edge(A, B, 1);
    ptd_add_edge(B, C, 1);
    ptd_add_edge(C, A, 1);

    ptd_add_edge(D, B, 1);
    ptd_add_edge(D, C, 1);
    ptd_add_edge(D, E, 1);
    ptd_add_edge(E, D, 1);
    ptd_add_edge(E, F, 1);

    ptd_add_edge(F, C, 1);
    ptd_add_edge(F, G, 1);
    ptd_add_edge(G, F, 1);

    ptd_add_edge(H, G, 1);
    ptd_add_edge(H, E, 1);

    draw_graph(graph);

    ptd_strongly_connected_components_t *sccs = ptd_find_strongly_connected_components(graph, keep_all);

    for (size_t i = 0; i < sccs->components_length; ++i) {
        fprintf(stderr, "\nComponents %zu:\n", i);

        for (size_t j = 0; j < sccs->components[i]->vertices_length; ++j) {
            fprintf(stderr, "Vertex %zu %s\n", j, vertex_name((sccs->components[i]->vertices[j])));
        }
    }

    ptd_strongly_connected_components_destroy(sccs);
    ptd_graph_destroy(graph);
}

bool remove_D_and_H_and_I(ptd_vertex_t *vertex) {
    return !(vertex == D || vertex == H || vertex == I);

}

void test_can_find_scc2() {
    ptd_graph_t *graph = ptd_graph_create(4);

    create_vertices(graph);

    ptd_add_edge(S, A, 1);

    ptd_add_edge(A, B, 1);
    ptd_add_edge(B, C, 1);
    ptd_add_edge(C, A, 1);
    ptd_add_edge(B, D, 1);
    ptd_add_edge(D, B, 1);
    ptd_add_edge(D, E, 1);
    ptd_add_edge(E, F, 1);
    ptd_add_edge(E, I, 1);
    ptd_add_edge(F, E, 1);
    ptd_add_edge(F, G, 1);
    ptd_add_edge(G, D, 1);
    ptd_add_edge(G, H, 1);
    ptd_add_edge(H, F, 1);

    ptd_add_edge(H, T, 1);
    ptd_add_edge(I, T, 1);

    draw_graph(graph);

    ptd_strongly_connected_components_t *sccs = ptd_find_strongly_connected_components(
            graph, remove_D_and_H_and_I
    );

    for (size_t i = 0; i < sccs->components_length; ++i) {
        fprintf(stderr, "\nComponents %zu:\n", i);

        for (size_t j = 0; j < sccs->components[i]->vertices_length; ++j) {
            fprintf(stderr, "Vertex %zu %s\n", j, vertex_name((sccs->components[i]->vertices[j])));
        }
    }

    ptd_strongly_connected_components_destroy(sccs);
    ptd_graph_destroy(graph);
}

double reward_one(ptd_vertex_t *vertex) {
    return 1;
}

int set_data_as_int(ptd_vertex_t *vertex) {
    vertex->data = malloc(sizeof(double));
    *((double *) vertex->data) = 0;

    return 0;
}

void test_can_find_scc2_graph() {
    ptd_graph_t *graph = ptd_graph_create(4);

    create_vertices(graph);

    ptd_add_edge(S, A, 1);

    ptd_add_edge(A, B, 3);
    ptd_add_edge(B, C, 4);
    ptd_add_edge(C, A, 4);
    ptd_add_edge(B, D, 2);
    ptd_add_edge(D, E, 5);
    ptd_add_edge(E, F, 5);
    ptd_add_edge(E, I, 15);
    ptd_add_edge(F, E, 1);
    ptd_add_edge(F, G, 1);
    ptd_add_edge(G, H, 1);
    ptd_add_edge(H, F, 1);
    ptd_add_edge(H, G, 1);

    ptd_add_edge(H, T, 1);
    ptd_add_edge(I, L, 1);
    ptd_add_edge(L, T, 1);

    ptd_add_edge(C, J, 1);

    ptd_add_edge(J, K, 1);
    ptd_add_edge(K, J, 1);

    ptd_add_edge(K, E, 1);



    long double **mat;
    ptd_vertex_t **vs;
    size_t length;

    ptd_phase_type_distribution_t *ptd = ptd_graph_as_phase_type_distribution(graph);
    mat = ptd->sub_intensity_matrix;
    vs = ptd->vertices;
    length = ptd->length;



    ptd_desc_multipliers_t *multipliers = ptd_cyclic_descendant_multipliers(graph);

    ptd_cyclic_desc(graph, multipliers,  reward_identity);

    draw_graph(graph);
    for (size_t j = 0; j < length; ++j) {
        fprintf(stderr, "%s ", vertex_name(vs[j]));
    }

    fprintf(stderr, "\n");
    fprintf(stderr, "MATRIX:\n");
    gsl_matrix *full = gsl_matrix_alloc(length, length);

    for (size_t k = 0; k < length; ++k) {
        for (size_t j = 0; j < length; ++j) {
            fprintf(stderr, "%.2Lf ", mat[k][j]);
            gsl_matrix_set(full, k, j, (double) mat[k][j]);
        }

        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "INVERTED:\n");
    gsl_matrix *inv = matrix_invert(full, length);

    for (size_t k = 0; k < length; ++k) {
        for (size_t j = 0; j < length; ++j) {
            fprintf(stderr, "%.2f ", -gsl_matrix_get(inv, k, j));
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");



    for (size_t k = 0; k < length; ++k) {

        double desc = 0;

        for (size_t j = 0; j < length; ++j) {
            desc += gsl_matrix_get(inv, k, j) / vs[j]->rate;
        }
        fprintf(stderr, "P: Vertex %zu (%s) has %f\n", vs[k]->index, vertex_name(vs[k]), desc);

    }

    ptd_graph_destroy(graph);
}

void test_can_find_scc2_rand_graph() {
    ptd_graph_t *graph = ptd_graph_create(4);
    srand(1234);

    ptd_vertex_t **vertices = (ptd_vertex_t**) calloc(50, sizeof(*vertices));
    double *fulld = (double*) calloc(50, sizeof(*fulld));

    ptd_vertex_t *abs = ptd_vertex_create(graph);

    for (size_t i = 0; i < 50; ++i) {
        vertices[i] = ptd_vertex_create(graph);
        vertices[i]->state[0] = i + 1;
        ptd_add_edge(vertices[i], abs, 2);
    }

    ptd_add_edge(graph->start_vertex, vertices[0], 1);

    for (size_t i = 0; i < 100; ++i) {
        int x = rand() % 50;
        int y = rand() % 50;

        if (x == y) {
            continue;
        }

        assert(ptd_add_edge(vertices[x], vertices[y], rand() % 2) == 0);
    }

    long double **mat;
    ptd_vertex_t **vs;
    size_t length;

    ptd_phase_type_distribution_t *ptd = ptd_graph_as_phase_type_distribution(graph);
    mat = ptd->sub_intensity_matrix;
    vs = ptd->vertices;
    length = ptd->length;

    ptd_desc_multipliers_t *multipliers = ptd_cyclic_descendant_multipliers(graph);

    double *algod = ptd_cyclic_desc(graph, multipliers,  reward_identity);

    gsl_matrix *full = gsl_matrix_alloc(length, length);

    ptd_strongly_connected_components_t *sccs = ptd_find_strongly_connected_components(graph, keep_all);

    fprintf(stderr, "Sccs: %zu\n", sccs->components_length);

    ptd_strongly_connected_components_destroy(sccs);

    for (size_t k = 0; k < length; ++k) {
        for (size_t j = 0; j < length; ++j) {
            gsl_matrix_set(full, k, j, (double) mat[k][j]);
        }
    }

    gsl_matrix *inv = matrix_invert(full, length);

    for (size_t k = 0; k < length; ++k) {
        double desc = 0;

        for (size_t j = 0; j < length; ++j) {
            desc += -gsl_matrix_get(inv, k, j) / vs[j]->rate;
        }

        fulld[vs[k]->index] = desc;
    }

    queue<ptd_vertex_t *> q = ptd_enqueue_vertices(graph);

    while (!q.empty()) {
        ptd_vertex_t *v = q.front();
        q.pop();

        fprintf(stdout, "Vertex %zu: %f vs %f\n", v->index, fulld[v->index], algod[v->index]);

        assert(fabs(fulld[v->index] - algod[v->index]) < 0.01);
    }

    ptd_phase_type_distribution_destroy(ptd);
    gsl_matrix_free(inv);
    gsl_matrix_free(full);
    free(algod);
    free(fulld);

    for (size_t i = 0; i < 50; ++i) {
        ptd_vertex_destroy(vertices[i]);
    }

    ptd_vertex_destroy(abs);

    free(vertices);

    ptd_graph_destroy(graph);
}

void test_it_can_insert_avl() {
    ptd_graph_t *graph = ptd_graph_create(4);
    assert(graph != NULL);
    ptd_avl_tree_t *tree = ptd_avl_tree_create(4);
    assert(tree != NULL);
    create_vertices(graph);

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
    ptd_graph_destroy(graph);
}


void test_avl_is_balanced() {
    ptd_graph_t *graph = ptd_graph_create(4);
    ptd_avl_tree_t *tree = ptd_avl_tree_create(4);
    vector<ptd_vertex_t *> vertices;

    for (size_t i = 0; i < 1024; ++i) {
        ptd_vertex_t *vertex = ptd_vertex_create(graph);
        assert(vertex != NULL);
        vertex->state[0] = (size_t) rand();
        ptd_avl_tree_vertex_insert(tree, vertex->state, vertex);
        assert(ptd_avl_tree_vertex_find(tree, vertex->state) == vertex);
        vertices.push_back(vertex);
    }

    for (size_t k = 0; k < vertices.size(); ++k) {
        ptd_vertex_destroy(vertices[k]);
    }

    fprintf(stderr, "Depth: %zu\n", ptd_avl_tree_max_depth(tree->root));
    assert(ptd_avl_tree_max_depth(tree->root) < 16);
    assert(ptd_avl_tree_max_depth(tree->root) > 2);
    ptd_avl_tree_vertex_destroy(tree);
    ptd_graph_destroy(graph);
}

void test_avl_is_balanced_and_updated() {
    ptd_graph_t *graph = ptd_graph_create(4);
    ptd_avl_tree_t *tree = ptd_avl_tree_create(4);

    vector<ptd_vertex_t *> vertices;

    for (size_t i = 0; i < 100; ++i) {
        for (size_t j = 0; j < 100; ++j) {
            ptd_vertex_t *vertex = ptd_vertex_create(graph);
            assert(vertex != NULL);
            vertex->state[0] = (size_t) j;
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
        ptd_vertex_destroy(vertices[k]);
    }

    fprintf(stderr, "Depth: %zu\n", ptd_avl_tree_max_depth(tree->root));
    assert(ptd_avl_tree_max_depth(tree->root) < 16);
    assert(ptd_avl_tree_max_depth(tree->root) > 2);
    ptd_avl_tree_edge_destroy(tree);
    ptd_graph_destroy(graph);
}


void test_avl_is_balanced_and_removed() {
    ptd_graph_t *graph = ptd_graph_create(4);
    ptd_avl_tree_t *tree = ptd_avl_tree_create(4);

    vector<ptd_vertex_t *> vertices;


    for (size_t i = 0; i < 100; ++i) {
        ptd_vertex_t *vertex = ptd_vertex_create(graph);
        assert(vertex != NULL);
        vertex->state[0] = (size_t) i;
        assert(ptd_avl_tree_edge_insert_or_increment(tree, vertex->state, vertex, 1) == 0);
        assert(ptd_avl_tree_edge_find(tree, vertex->state)->to == vertex);
        assert(ptd_avl_tree_edge_find(tree, vertex->state)->weight == 1);
        ptd_avl_tree_edge_remove(tree, vertex->state);
        vertices.push_back(vertex);
    }

    ptd_vertex_t *vertex[100];

    for (size_t i = 0; i < 100; ++i) {
        vertex[i] = ptd_vertex_create(graph);
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
        ptd_vertex_destroy(vertices[k]);
    }

    fprintf(stderr, "Depth: %zu\n", ptd_avl_tree_max_depth(tree->root));
    assert(ptd_avl_tree_max_depth(tree->root) < 16);
    assert(ptd_avl_tree_max_depth(tree->root) > 2);
    ptd_avl_tree_edge_destroy(tree);
    ptd_graph_destroy(graph);
}

static ptd_graph_t *kingman_graph;
static ptd_avl_tree_t *avl_tree;

int make_kingman(ptd_vertex_t *vertex) {
    vec_entry_t *state = vertex->state;

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

            vec_entry_t *child_state = (vec_entry_t *) calloc(kingman_graph->state_length, sizeof(vec_entry_t));

            if (child_state == NULL) {
                return -1;
            }

            memcpy(child_state, state, kingman_graph->state_length * sizeof(vec_entry_t));
            child_state[i]--;
            child_state[j]--;
            child_state[(i + j + 2) - 1]++;

            ptd_vertex_t *child = ptd_avl_tree_vertex_find(avl_tree, child_state);

            if (child == NULL) {
                child = ptd_vertex_create_state(kingman_graph, child_state);

                if (ptd_avl_tree_vertex_insert(avl_tree, child_state, child)) {
                    return -1;
                }
            } else {
                free(child_state);
            }

            ptd_add_edge(vertex, child, weight);
        }
    }

    return 0;
}

void test_it_can_construct_kingman() {
    avl_tree = ptd_avl_tree_create(5);
    kingman_graph = ptd_graph_create(5);
    ptd_vertex_t *initial = ptd_vertex_create(kingman_graph);
    initial->state[0] = 5;
    ptd_add_edge(kingman_graph->start_vertex, initial, 1);
    ptd_visit_vertices(kingman_graph, make_kingman, false);
    draw_graph(kingman_graph);
    ptd_avl_tree_vertex_destroy(avl_tree);
    ptd_graph_destroy(kingman_graph);
}

double reward_identity(const ptd_vertex_t *vertex) {
    return 1;
}

void test_it_can_reward_transform() {
    ptd_graph_t *graph = ptd_graph_create(4);

    create_vertices(graph);

    ptd_add_edge(S, A, 1);

    ptd_add_edge(A, B, 1);
    ptd_add_edge(B, C, 2);
    ptd_add_edge(C, A, 3);
    ptd_add_edge(B, D, 4);
    ptd_add_edge(D, B, 5);
    ptd_add_edge(D, T, 6);
    ptd_add_edge(D, E, 7);
    ptd_add_edge(E, A, 8.2);

    draw_graph(graph);

    ptd_reward_transform(graph, reward_identity);
    draw_graph(graph);

    ptd_graph_destroy(graph);
}

int main(int argc, char **argv) {
//    test_can_find_scc();
//    test_can_find_scc2();
    //test_can_find_scc2_graph();
    test_can_find_scc2_rand_graph();

//    test_it_can_insert_avl();
//    test_avl_is_balanced();
    //   test_avl_is_balanced_and_updated();
    // test_avl_is_balanced_and_removed();
    //  test_it_can_construct_kingman();
    //test_it_can_reward_transform();
}

vertex_t *s, *a, *b, *c, *d, *e, *f, *t;

double reward_by_index(vertex_t *vertex) {
    if (vertex == s) {
        return 1;
    }

    if (vertex == a) {
        return 0;
    }

    if (vertex == b) {
        return 1;
    }

    if (vertex == c) {
        return 1;
    }

    if (vertex == d) {
        return 4;
    }

    if (vertex == e) {
        return 5;
    }

    if (vertex == f) {
        return 1;
    }

    if (vertex == t) {
        return 1;
    }

    exit(1);
}


int main2(int argc, char **argv) {
    vec_entry_t *states = (vec_entry_t *) calloc(10, sizeof(vec_entry_t));
    vec_entry_t *statea = (vec_entry_t *) calloc(10, sizeof(vec_entry_t));
    vec_entry_t *stateb = (vec_entry_t *) calloc(10, sizeof(vec_entry_t));
    vec_entry_t *statec = (vec_entry_t *) calloc(10, sizeof(vec_entry_t));
    vec_entry_t *stated = (vec_entry_t *) calloc(10, sizeof(vec_entry_t));
    vec_entry_t *statee = (vec_entry_t *) calloc(10, sizeof(vec_entry_t));
    vec_entry_t *statef = (vec_entry_t *) calloc(10, sizeof(vec_entry_t));
    vec_entry_t *statet = (vec_entry_t *) calloc(10, sizeof(vec_entry_t));

    *states = 1;
    *statea = 2;
    *stateb = 3;
    *statec = 4;
    *stated = 5;
    *statee = 6;
    *statef = 7;
    *statet = 8;

    s = vertex_init(states, vector<double>(), 0);
    a = vertex_init(statea, vector<double>(), 0);
    b = vertex_init(stateb, vector<double>(), 0);
    c = vertex_init(statec, vector<double>(), 0);
    d = vertex_init(stated, vector<double>(), 0);
    e = vertex_init(statee, vector<double>(), 0);
    f = vertex_init(statef, vector<double>(), 0);
    t = vertex_init(statet, vector<double>(), 0);

    vertex_add_edge(s, a, 0.3);
    vertex_add_edge(s, b, 0.7);
    vertex_add_edge(a, b, 3);
    vertex_add_edge(b, a, 1);
    vertex_add_edge(b, c, 2);
    vertex_add_edge(b, d, 5);
    vertex_add_edge(c, b, 6);
    vertex_add_edge(c, t, 6);
    vertex_add_edge(d, c, 7);

    fprintf(stderr, "Vertex a %p\n", (void *) a);
    fprintf(stderr, "Vertex b %p\n", (void *) b);
    fprintf(stderr, "Vertex c %p\n", (void *) c);
    fprintf(stderr, "Vertex d %p\n", (void *) d);
    fprintf(stderr, "Vertex s %p\n", (void *) s);
    fprintf(stderr, "Vertex t %p\n", (void *) t);

    struct graph_info graph_info = get_graph_info(s);
    fprintf(stderr, "Vertices %zu, edges %zu\n", graph_info.vertices, graph_info.edges);

    reward_transform(s, reward_by_index);

    vertex_t **vertices;
    double **mat;
    size_t size;

    graph_as_mat(&mat, &size, &vertices, s);

    for (size_t i = 2; i < size; ++i) {
        fprintf(stdout, "%f ", mat[1][i]);
    }

    fprintf(stdout, "\n\n");

    for (size_t i = 2; i < size; ++i) {
        for (size_t j = 2; j < size; ++j) {
            fprintf(stdout, "%f ", mat[i][j]);
        }

        fprintf(stdout, "\n");
    }

    for (size_t i = 0; i < size; ++i) {
        free(mat[i]);
    }

    free(mat);
    free(vertices);

    graph_free(s);

    return 0;
}

/*
int main(int argv, char **argc) {
    vertex_t *a = vertex_init();
    vertex_t *b = vertex_init();
    vertex_t *c = vertex_init();
    vertex_t *d = vertex_init();
    vertex_t *e = vertex_init();
    vertex_t *f = vertex_init();

    vertex_t *me = vertex_init();
    ll_insert(me, b, 1);
    ll_insert(me, c, 1);
    ll_insert(me, a, 1);
    ll_insert(me, f, 1);
    ll_insert(me, e, 1);
    ll_insert(me, d, 1);

    assert(me->edges->vertex == a);
    assert(me->edges->next->vertex == b);
    assert(me->edges->next->next->vertex == c);
    assert(me->edges->next->next->next->vertex == d);
    assert(me->edges->next->next->next->next->vertex == e);
    assert(me->edges->next->next->next->next->next->vertex == f);

    me = vertex_init();
    ll_insert_p(b, me);
    ll_insert_p(c, me);
    ll_insert_p(a, me);
    ll_insert_p(f, me);
    ll_insert_p(e, me);
    ll_insert_p(d, me);

    assert(me->parents->vertex == a);
    assert(me->parents->next->vertex == b);
    assert(me->parents->next->next->vertex == c);
    assert(me->parents->next->next->next->vertex == d);
    assert(me->parents->next->next->next->next->vertex == e);
    assert(me->parents->next->next->next->next->next->vertex == f);

    me = vertex_init();
    ll_insert_p_existing(b, me);
    ll_insert_p_existing(b, me);
    ll_insert_p_existing(c, me);
    ll_insert_p_existing(c, me);
    ll_insert_p_existing(a, me);
    ll_insert_p_existing(a, me);
    ll_insert_p_existing(f, me);
    ll_insert_p_existing(f, me);
    ll_insert_p_existing(e, me);
    ll_insert_p_existing(e, me);
    ll_insert_p_existing(d, me);
    ll_insert_p_existing(d, me);
    ll_insert_p_existing(b, me);
    ll_insert_p_existing(b, me);
    ll_insert_p_existing(c, me);
    ll_insert_p_existing(c, me);
    ll_insert_p_existing(a, me);
    ll_insert_p_existing(a, me);
    ll_insert_p_existing(f, me);
    ll_insert_p_existing(f, me);
    ll_insert_p_existing(e, me);
    ll_insert_p_existing(e, me);
    ll_insert_p_existing(d, me);
    ll_insert_p_existing(d, me);

    assert(me->parents->vertex == a);
    assert(me->parents->next->vertex == b);
    assert(me->parents->next->next->vertex == c);
    assert(me->parents->next->next->next->vertex == d);
    assert(me->parents->next->next->next->next->vertex == e);
    assert(me->parents->next->next->next->next->next->vertex == f);

    a = vertex_init();
    b = vertex_init();
    c = vertex_init();
    d = vertex_init();
    e = vertex_init();
    f = vertex_init();
    me = vertex_init();
    vertex_t *other = vertex_init();
    ll_insert(me, c, 1);
    ll_insert(me, d, 1);

    ll_insert(other, a, 1);
    ll_insert(other, d, 1);
    ll_insert(other, f, 1);
    ll_insert(other, b, 1);
    ll_insert(other, e, 1);

    ll_insert_or_inc_list(me, other->edges);

    assert(me->edges->vertex == a);
    assert(me->edges->next->vertex == b);
    assert(me->edges->next->next->vertex == c);
    assert(me->edges->next->next->next->vertex == d);
    assert(me->edges->next->next->next->next->vertex == e);
    assert(me->edges->next->next->next->next->next->vertex == f);

    assert(me->edges->weight == 1);
    assert(me->edges->next->weight == 1);
    assert(me->edges->next->next->weight == 1);
    assert(me->edges->next->next->next->weight == 2);
    assert(me->edges->next->next->next->next->weight == 1);
    assert(me->edges->next->next->next->next->next->weight == 1);

    ll_insert_or_inc_list(me, other->edges);

    assert(me->edges->vertex == a);
    assert(me->edges->next->vertex == b);
    assert(me->edges->next->next->vertex == c);
    assert(me->edges->next->next->next->vertex == d);
    assert(me->edges->next->next->next->next->vertex == e);
    assert(me->edges->next->next->next->next->next->vertex == f);

    assert(me->edges->weight == 2);
    assert(me->edges->next->weight == 2);
    assert(me->edges->next->next->weight == 1);
    assert(me->edges->next->next->next->weight == 3);
    assert(me->edges->next->next->next->next->weight == 2);
    assert(me->edges->next->next->next->next->next->weight == 2);

    assert(a->parents->vertex == me);
    assert(a->parents->next->vertex == other);

    assert(b->parents->vertex == me);
    assert(b->parents->next->vertex == other);

    assert(c->parents->vertex == me);

    assert(d->parents->vertex == me);
    assert(d->parents->next->vertex == other);

    assert(e->parents->vertex == me);
    assert(e->parents->next->vertex == other);

    assert(f->parents->vertex == me);
    assert(f->parents->next->vertex == other);

    return 0;
}*/

vector<pair<double, vector<size_t> > > visit_function(vector<size_t> state) {
    vector<pair<double, vector<size_t> > > children;

    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = i; j < 4; ++j) {
            double rate;

            if (i == j) {
                if (state[i] <= 1) {
                    continue;
                }

                rate = state[i] * (state[i] - 1) / 2;
            } else {
                if (state[i] == 0 || state[j] == 0) {
                    continue;
                }

                rate = state[i] * state[j];
            }

            const size_t inc_pos = (i + j + 2) - 1;

            vector<size_t> child_state(state);
            child_state[i]--;
            child_state[j]--;
            child_state[inc_pos]++;
            children.push_back(pair<double, vector<size_t> >(rate, child_state));
        }
    }

    return children;
}

vector<pair<double, vector<size_t> > > initial_states(void) {
    vector<pair<double, vector<size_t> > > children;
    vector<size_t> child_state;
    child_state.push_back(4);
    child_state.push_back(0);
    child_state.push_back(0);
    child_state.push_back(0);
    children.push_back(pair<double, vector<size_t> >(1.0f, child_state));

    return children;
}

vector<double> rewards(vector<size_t> input) {
    vector<double> rew(input.begin(), input.end());

    return rew;
}

size_t n1, n2, nrow, ncol, t1, t2, matrix_size, state_size;
double m;

void mig_set_state_properties(size_t num_pop1, size_t num_pop2,
                              size_t types_pop1, size_t types_pop2, double mig) {
    n1 = num_pop1;
    n2 = num_pop2;
    nrow = types_pop1 + 1;
    ncol = types_pop2 + 1;
    matrix_size = nrow * ncol;
    state_size = matrix_size * 2;
    t1 = types_pop1;
    t2 = types_pop2;
    m = mig;
}

size_t mig_state_index(size_t p1, size_t p2, size_t population) {
    if (p1 > t1) {
        p1 = t1;
    }

    if (p2 > t2) {
        p2 = t2;
    }

    return nrow * (p2) + p1 + (population == 1 ? 0 : matrix_size);
}

vector<pair<double, vector<size_t> > > mig_visit_function(vector<size_t> state) {
    vector<pair<double, vector<size_t> > > children;

    if (accumulate(state.begin(), state.end(), 0) == 1) {
        return children;
    }

    for (size_t p = 1; p <= 2; p++) {
        for (size_t i1 = 0; i1 <= t1; i1++) {
            for (size_t j1 = 0; j1 <= t2; j1++) {
                size_t idx1 = mig_state_index(i1, j1, p);

                {
                    size_t idx2 = mig_state_index(i1, j1, 2 - p + 1);
                    // Migration from the islands
                    if (state[idx1] >= 1) {
                        double rate = m * state[idx1];

                        vector<size_t> child_state(state);

                        child_state[idx1]--;
                        child_state[idx2]++;

                        children.push_back(
                                pair<double, vector<size_t> >(rate, child_state)
                        );
                    }
                }

                if (state[idx1] != 0) {
                    // Coalescence
                    for (size_t i2 = 0; i2 <= t1; i2++) {
                        for (size_t j2 = 0; j2 <= t2; j2++) {
                            double rate;
                            size_t idx2 = mig_state_index(i2, j2, p);

                            if (i1 == i2 && j1 == j2) {
                                if (state[idx1] <= 1) {
                                    continue;
                                }

                                rate = (double) (state[idx1] * (state[idx1] - 1)) / 2;
                            } else {
                                if (state[idx1] == 0 || state[idx2] == 0) {
                                    continue;
                                }

                                rate = (double) state[idx1] * state[idx2];
                            }

                            vector<size_t> child_state(state);
                            size_t idx3 = mig_state_index(i1 + i2, j1 + j2, p);

                            child_state[idx1]--;
                            child_state[idx2]--;
                            child_state[idx3]++;

                            children.push_back(
                                    pair<double, vector<size_t> >(rate, child_state)
                            );
                        }
                    }
                }
            }
        }
    }

    return (children);
}

vector<pair<double, vector<size_t> > > mig_initial_states(void) {
    vector<pair<double, vector<size_t> > > children;
    vector<size_t> state = vector<size_t>(state_size, 0);

    state[mig_state_index(1, 0, 1)] = n1;
    state[mig_state_index(0, 1, 2)] = n2;

    children.push_back(pair<double, vector<size_t> >(1.0f, state));

    return children;
}

vector<double> mig_rewards(vector<size_t> input) {
    vector<double> rew(input.begin(), input.end());

    return rew;
}

void set_bool_for_parents(vertex_t *vertex) {
    if (vertex->boolean) {
        return;
    }

    vertex->boolean = true;

    llp_t *parent = vertex->parents->next;

    while (parent != NULL) {
        set_bool_for_parents(parent->parent);

        parent = parent->next;
    }
}

int main15(int argc, char **argv) {
    vertex_t *g;
    mph_cov_exp_all(g, 4);

    queue<vertex_t *> qu = enqueue_vertices(g);
    vector<vertex_t *> vec;

    while (!qu.empty()) {
        vec.push_back(qu.front());
        qu.pop();
    }


    mig_set_state_properties(5, 5, 5, 1, 0.1);

    vertex_t *graph = generate_state_space(
            state_size,
            mig_visit_function,
            mig_initial_states,
            mig_rewards
    );

    queue<vertex_t *> queue = enqueue_vertices(graph);

    while (!queue.empty()) {
        vertex_t *vertex = queue.front();
        queue.pop();

        vertex->boolean = false;
    }

    queue = enqueue_vertices(graph);

    while (!queue.empty()) {
        vertex_t *vertex = queue.front();
        queue.pop();

        if (vertex->state != NULL &&
            (vertex->state[mig_state_index(5, 5, 1)] != 0)) {
            set_bool_for_parents(vertex);
        }
    }

    size_t nt = 0, nf = 0;

    queue = enqueue_vertices(graph);

    while (!queue.empty()) {
        vertex_t *vertex = queue.front();
        queue.pop();

        if (vertex->boolean) {
            nt++;
        } else {
            nf++;
        }
    }

    fprintf(stderr, "nt %zu, nf %zu\n", nt, nf);

    struct graph_info info = get_graph_info(graph);
    fprintf(stderr, "Vertices %zu edges %zu\n", info.vertices, info.edges);

    return 0;
    vertex_t **vertices;
    double **mat;
    size_t size;

    graph_as_mat(&mat, &size, &vertices, graph);

    for (size_t i = 2; i < size; ++i) {
        fprintf(stdout, "%f ", mat[1][i]);
    }

    fprintf(stdout, "\n\n");

    for (size_t i = 2; i < size; ++i) {
        for (size_t j = 2; j < size; ++j) {
            fprintf(stdout, "%f ", mat[i][j]);
        }

        fprintf(stdout, "\n");
    }

    for (size_t i = 0; i < size; ++i) {
        free(mat[i]);
    }

    free(mat);
    free(vertices);

    graph_free(graph);

    return 0;
}
