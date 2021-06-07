#ifndef PTDALGORITHMS_PTD_H
#define PTDALGORITHMS_PTD_H

#include <vector>
#include <utility>
#include <queue>
#include <stdlib.h>
#include <stdio.h>
#include <stdexcept>

#define DIE_ERROR(error_code, error, ...) do {     \
char error_formatted[1024];                        \
char error_formatted_line[1024];                   \
                                                   \
snprintf(error_formatted,                          \
         sizeof(error_formatted),                  \
         error, ##__VA_ARGS__);                    \
snprintf(error_formatted_line,                     \
         sizeof(error_formatted_line),             \
         "%s @ %s (%d)", error_formatted,          \
         __FILE__, __LINE__);                      \
                                                   \
fprintf(stderr, "%s\n", error_formatted_line);     \
/*exit(error_code);*/throw std::runtime_error(error_formatted_line); \
} while(0)

#define DEBUG_PRINT(message, ...) do {             \
char formatted[2048];                              \
                                                   \
snprintf(formatted,                                \
         sizeof(formatted),                        \
         message, ##__VA_ARGS__);                  \
                                                   \
fprintf(stderr, "%s", formatted);                  \
} while(0)

using namespace std;

struct avl_node {
    struct avl_node *left;
    struct avl_node *right;
    struct avl_node *parent;
    signed short balance;
    char *key;
    void *entry;
};

struct ptd_directed_graph;
struct ptd_directed_edge;
struct ptd_directed_vertex;

struct ptd_directed_graph {
    size_t vertices_length;
    struct ptd_directed_vertex **vertices;
    struct ptd_directed_vertex *starting_vertex;
};

struct ptd_directed_edge {
    struct ptd_directed_vertex *to;
};

struct ptd_directed_vertex {
    size_t edges_length;
    struct ptd_directed_edge **edges;
    struct ptd_directed_graph *graph;
    size_t index;
};

int ptd_directed_graph_add_edge(struct ptd_directed_vertex *vertex, struct ptd_directed_edge *edge);

void ptd_directed_graph_destroy(struct ptd_directed_graph *graph);

int ptd_directed_vertex_add(struct ptd_directed_graph *graph, struct ptd_directed_vertex *vertex);

void ptd_directed_vertex_destroy(struct ptd_directed_vertex *vertex);

struct ptd_ph_graph;
struct ptd_ph_edge;
struct ptd_ph_vertex;

struct ptd_ph_graph {
    size_t vertices_length;
    struct ptd_ph_vertex **vertices;
    struct ptd_ph_vertex *starting_vertex;
    size_t state_length;
};

struct ptd_ph_edge {
    struct ptd_ph_vertex *to;
    double weight;
};

struct ptd_ph_vertex {
    size_t edges_length;
    struct ptd_ph_edge **edges;
    struct ptd_ph_graph *graph;
    size_t index;
    int *state;
};

struct ptd_ph_graph *ptd_ph_graph_create(size_t state_length);

void ptd_ph_graph_destroy(struct ptd_ph_graph *graph);

struct ptd_ph_vertex *ptd_ph_vertex_create(struct ptd_ph_graph *graph);

struct ptd_ph_vertex *ptd_ph_vertex_create_state(
        struct ptd_ph_graph *graph,
        int *state
);

double ptd_ph_vertex_rate(struct ptd_ph_vertex *vertex);

void ptd_ph_vertex_destroy(struct ptd_ph_vertex *vertex);

struct ptd_ph_edge *ptd_ph_graph_add_edge(
        struct ptd_ph_vertex *from,
        struct ptd_ph_vertex *to,
        double weight
);

struct ptd_ph_vertex **ptd_ph_graph_topological_sort(struct ptd_ph_graph *graph);

double *ptd_ph_graph_acyclic_visit_probability(struct ptd_ph_graph *graph);

double *ptd_ph_graph_acyclic_moment_rewards(struct ptd_ph_graph *graph, double *rewards);

struct ptd_ph_scc_graph;
struct ptd_ph_scc_edge;
struct ptd_ph_scc_vertex;

struct ptd_ph_scc_graph {
    size_t vertices_length;
    struct ptd_ph_scc_vertex **vertices;
    struct ptd_ph_scc_vertex *starting_vertex;
    struct ptd_ph_graph *graph;
};

struct ptd_ph_scc_edge {
    struct ptd_ph_scc_vertex *to;
};

struct ptd_ph_scc_vertex {
    size_t edges_length;
    struct ptd_ph_scc_edge **edges;
    struct ptd_ph_scc_graph *graph;
    size_t index;
    size_t internal_vertices_length;
    struct ptd_ph_vertex **internal_vertices;
    size_t external_vertices_length;
    struct ptd_ph_vertex **external_vertices;
    double **internal_expected_visits;
    double **external_expected_visits;
};

struct ptd_ph_scc_graph *ptd_ph_find_strongly_connected_components(struct ptd_ph_graph *graph);

void ptd_ph_scc_graph_destroy(struct ptd_ph_scc_graph *scc_graph);

double *ptd_ph_graph_cyclic_entry_probability(struct ptd_ph_scc_graph *scc_graph);

double *ptd_ph_scc_graph_entry_probability(struct ptd_ph_scc_graph *scc_graph);

double *ptd_ph_graph_cyclic_expected_visits(struct ptd_ph_scc_graph *scc_graph);

double *ptd_ph_graph_cyclic_moment_rewards(
        struct ptd_ph_scc_graph *scc_graph, double *rewards
);

int ptd_ph_graph_reward_transform(struct ptd_ph_graph *graph, double *rewards);

typedef struct ptd_avl_tree {
    void *root;
    size_t vec_length;
} ptd_avl_tree_t;

ptd_avl_tree_t *ptd_avl_tree_create(size_t vec_length);

void ptd_avl_tree_vertex_destroy(ptd_avl_tree_t *avl_tree);

void ptd_avl_tree_vertex_destroy_free(ptd_avl_tree_t *avl_tree);

void ptd_avl_tree_edge_destroy(ptd_avl_tree_t *avl_tree);

int ptd_avl_tree_vertex_insert(ptd_avl_tree_t *avl_tree, const int *key, const struct ptd_ph_vertex *vertex);

struct ptd_ph_vertex *ptd_avl_tree_vertex_find(const ptd_avl_tree_t *avl_tree, const int *key);

int ptd_avl_tree_edge_insert_or_increment(ptd_avl_tree_t *avl_tree, const int *key, struct ptd_ph_vertex *vertex,
                                          double weight);

int ptd_avl_tree_edge_remove(ptd_avl_tree_t *avl_tree, const int *key);

struct ptd_ph_edge *ptd_avl_tree_edge_find(const ptd_avl_tree_t *avl_tree, const int *key);

//TODO: remove?
size_t ptd_avl_tree_max_depth(void *avl_vec_vertex);

typedef struct {
    size_t length;
    double *initial_probability_vector;
    double **sub_intensity_matrix;
    struct ptd_ph_vertex **vertices;
    size_t memory_allocated;
} ptd_phase_type_distribution_t;

ptd_phase_type_distribution_t *ptd_graph_as_phase_type_distribution(struct ptd_ph_graph *graph);

void ptd_phase_type_distribution_destroy(ptd_phase_type_distribution_t *ptd);

int ptd_ph_vertex_to_s(struct ptd_ph_vertex *vertex, char *buffer, size_t buffer_length);

/*
 * Models
 */

struct ptd_ph_graph *ptd_model_kingman(size_t n);

struct ptd_ph_graph *ptd_model_two_island_two_loci_recomb(
        size_t n1,
        size_t n2,
        size_t effective_population_size1,
        size_t effective_population_size2,
        double migration_rate1,
        double migration_rate2,
        double recombination_rate
);

#endif //PTDALGORITHMS_PTD_H
