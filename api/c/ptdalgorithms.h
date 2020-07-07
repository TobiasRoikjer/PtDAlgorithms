/*
 * MIT License
 *
 * Copyright (c) 2021 Tobias Røikjer
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef PTDALGORITHMS_PTD_H
#define PTDALGORITHMS_PTD_H

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ptd_avl_node {
    struct ptd_avl_node *left;
    struct ptd_avl_node *right;
    struct ptd_avl_node *parent;
    signed short balance;
    int *key;
    void *entry;
};

struct ptd_avl_tree {
    struct ptd_avl_node *root;
    size_t key_length;
};

struct ptd_graph;
struct ptd_edge;
struct ptd_vertex;

struct ptd_scc_graph;
struct ptd_scc_edge;
struct ptd_scc_vertex;

extern volatile char ptd_err[4096];

#ifndef PTD_DEBUG_1_INDEX
#define PTD_DEBUG_1_INDEX 0
#endif

struct ptd_avl_tree *ptd_avl_tree_create(size_t key_length);

void ptd_avl_tree_destroy(struct ptd_avl_tree *avl_tree);

struct ptd_avl_node *ptd_avl_tree_find_or_insert(struct ptd_avl_tree *avl_tree, const int *key, const void *entry);

struct ptd_avl_node *ptd_avl_tree_find(const struct ptd_avl_tree *avl_tree, const int *key);

struct ptd_vertex *ptd_avl_tree_find_vertex(const struct ptd_avl_tree *avl_tree, const int *key);

size_t ptd_avl_tree_max_depth(void *avl_vec_vertex);

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

struct ptd_graph {
    size_t vertices_length;
    struct ptd_vertex **vertices;
    struct ptd_vertex *starting_vertex;
    size_t state_length;
    bool parameterized;
    struct ptd_desc_reward_compute *reward_compute_graph;
    struct ptd_desc_reward_compute_parameterized *parameterized_reward_compute_graph;
    bool was_dph;
};

struct ptd_edge {
    struct ptd_vertex *to;
    double weight;
    bool parameterized;
};

struct ptd_edge_parameterized {
    struct ptd_vertex *to;
    double weight;
    bool parameterized;
    double *state;
    bool should_free_state;
};


struct ptd_vertex {
    size_t edges_length;
    struct ptd_edge **edges;
    struct ptd_graph *graph;
    size_t index;
    int *state;
};

struct ptd_graph *ptd_graph_create(size_t state_length);

void ptd_graph_destroy(struct ptd_graph *graph);

struct ptd_vertex *ptd_vertex_create(struct ptd_graph *graph);

struct ptd_vertex *ptd_vertex_create_state(
        struct ptd_graph *graph,
        int *state
);

struct ptd_vertex *ptd_find_or_create_vertex(
        struct ptd_graph *graph,
        struct ptd_avl_tree *avl_tree,
        const int *child_state
);

double ptd_vertex_rate(struct ptd_vertex *vertex);

void ptd_vertex_destroy(struct ptd_vertex *vertex);

struct ptd_edge *ptd_graph_add_edge(
        struct ptd_vertex *from,
        struct ptd_vertex *to,
        double weight
);

struct ptd_edge_parameterized *ptd_graph_add_edge_parameterized(
        struct ptd_vertex *from,
        struct ptd_vertex *to,
        double weight,
        double *edge_state
);

void ptd_edge_update_weight(
        struct ptd_edge *edge,
        double weight
);

void ptd_notify_change(
        struct ptd_graph *graph
);

void ptd_edge_update_weight_parameterized(
        struct ptd_edge *edge,
        double *scalars,
        size_t scalars_length
);

void ptd_graph_update_weight_parameterized(
        struct ptd_graph *graph,
        double *scalars,
        size_t scalars_length
);

double *ptd_normalize_graph(struct ptd_graph *graph);

double *ptd_dph_normalize_graph(struct ptd_graph *graph);

double *ptd_expected_waiting_time(struct ptd_graph *graph, double *rewards);

double *ptd_dph_expected_waiting_time(struct ptd_graph *graph, double *rewards);

bool ptd_graph_is_acyclic(struct ptd_graph *graph);

struct ptd_vertex **ptd_graph_topological_sort(struct ptd_graph *graph);

struct ptd_graph *ptd_graph_reward_transform(struct ptd_graph *graph, double *rewards);

struct ptd_graph *ptd_graph_dph_reward_transform(struct ptd_graph *graph, int *rewards);

long double ptd_random_sample(struct ptd_graph *graph, double *rewards);

long double *ptd_mph_random_sample(struct ptd_graph *graph, double *rewards, size_t vertex_rewards_length);

long double ptd_dph_random_sample(struct ptd_graph *graph, double *rewards);

long double *ptd_mdph_random_sample(struct ptd_graph *graph, double *rewards, size_t vertex_rewards_length);

struct ptd_vertex *ptd_random_sample_stop_vertex(struct ptd_graph *graph, double time);

struct ptd_vertex *ptd_dph_random_sample_stop_vertex(struct ptd_graph *graph, int jumps);

double ptd_defect(struct ptd_graph *graph);

int ptd_validate_graph(const struct ptd_graph *graph);

struct ptd_clone_res {
    struct ptd_avl_tree *avl_tree;
    struct ptd_graph *graph;
};

struct ptd_clone_res ptd_clone_graph(struct ptd_graph *graph, struct ptd_avl_tree *avl_tree);

struct ptd_desc_reward_compute {
    size_t length;
    struct ptd_reward_increase *commands;
};


struct ptd_reward_increase {
    size_t from;
    size_t to;
    double multiplier;
};

struct ptd_comp_graph_parameterized {
    size_t from;
    size_t to;
    double multiplier;
    double *multiplierptr;
    int type;
    double *fromT;
    double *toT;
};

struct ptd_desc_reward_compute_parameterized {
    size_t length;
    struct ptd_comp_graph_parameterized *commands;
    void *mem;
    void *memr;
};

struct ptd_desc_reward_compute *ptd_graph_ex_absorbation_time_comp_graph(struct ptd_graph *graph);

struct ptd_desc_reward_compute_parameterized *
ptd_graph_ex_absorbation_time_comp_graph_parameterized(struct ptd_graph *graph);

struct ptd_desc_reward_compute *
ptd_graph_build_ex_absorbation_time_comp_graph_parameterized(struct ptd_desc_reward_compute_parameterized *compute);

void ptd_parameterized_reward_compute_graph_destroy(
        struct ptd_desc_reward_compute_parameterized *compute_graph
);


struct ptd_scc_graph {
    size_t vertices_length;
    struct ptd_scc_vertex **vertices;
    struct ptd_scc_vertex *starting_vertex;
    struct ptd_graph *graph;
};

struct ptd_scc_edge {
    struct ptd_scc_vertex *to;
};

struct ptd_scc_vertex {
    size_t edges_length;
    struct ptd_scc_edge **edges;
    struct ptd_scc_graph *graph;
    size_t index;
    size_t internal_vertices_length;
    struct ptd_vertex **internal_vertices;
};

int ptd_precompute_reward_compute_graph(struct ptd_graph *graph);

struct ptd_scc_graph *ptd_find_strongly_connected_components(struct ptd_graph *graph);

struct ptd_scc_vertex **ptd_scc_graph_topological_sort(struct ptd_scc_graph *graph);

void ptd_scc_graph_destroy(struct ptd_scc_graph *scc_graph);

struct ptd_phase_type_distribution {
    size_t length;
    double *initial_probability_vector;
    double **sub_intensity_matrix;
    struct ptd_vertex **vertices;
    size_t memory_allocated;
};

struct ptd_phase_type_distribution *ptd_graph_as_phase_type_distribution(struct ptd_graph *graph);

void ptd_phase_type_distribution_destroy(struct ptd_phase_type_distribution *ptd);

int ptd_vertex_to_s(struct ptd_vertex *vertex, char *buffer, size_t buffer_length);

struct ptd_probability_distribution_context {
    double pdf;
    double cdf;
    long double *probability_at;
    long double *accumulated_visits;
    struct ptd_graph *graph;
    void *priv;
    long double time;
    int granularity;
};

struct ptd_probability_distribution_context *ptd_probability_distribution_context_create(
        struct ptd_graph *graph,
        int granularity
);

void ptd_probability_distribution_context_destroy(
        struct ptd_probability_distribution_context *context
);

int ptd_probability_distribution_step(
        struct ptd_probability_distribution_context *context
);

struct ptd_dph_probability_distribution_context {
    double pmf;
    double cdf;
    long double *probability_at;
    long double *accumulated_visits;
    struct ptd_graph *graph;
    void *priv;
    size_t priv2;
    double priv3;
    int jumps;
};

struct ptd_dph_probability_distribution_context *ptd_dph_probability_distribution_context_create(
        struct ptd_graph *graph
);

void ptd_dph_probability_distribution_context_destroy(
        struct ptd_dph_probability_distribution_context *context
);

int ptd_dph_probability_distribution_step(
        struct ptd_dph_probability_distribution_context *context
);

#ifndef PTD_INTEGRATE_EXCEPTIONS
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
exit(error_code);   \
} while(0)
#else
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
throw std::runtime_error(error_formatted_line); \
} while(0)
#endif

#define DEBUG_PRINT(message, ...) do {             \
char formatted[2048];                              \
                                                   \
snprintf(formatted,                                \
         sizeof(formatted),                        \
         message, ##__VA_ARGS__);                  \
                                                   \
fprintf(stderr, "%s", formatted);                  \
} while(0)

#ifdef __cplusplus
}
#endif

#endif //PTDALGORITHMS_PTD_H
