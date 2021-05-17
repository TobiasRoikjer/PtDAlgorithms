#ifndef PTDALGORITHMS_PTD_H
#define PTDALGORITHMS_PTD_H

#include <vector>
#include <utility>
#include <queue>
#include <stdlib.h>
#include <stdio.h>
#include <stdexcept>

typedef struct {
    double ***cov;
    double **exp;
    double **defect;
} cov_exp_return;


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

typedef size_t vec_entry_t;

struct llc;
struct llp;

static size_t id = 0;


struct avl_node {
    struct avl_node *left;
    struct avl_node *right;
    struct avl_node *parent;
    signed short balance;
    char *key;
    void *entry;
};

typedef struct vertex {
    vertex(vec_entry_t *state, vector<double> rewards, size_t state_length) {
        this->vertex_index = id;
        id++;
        this->state = state;
        this->rewards = rewards;
        this->edges = NULL;
        this->parents = NULL;
        this->visited = false;
        this->rate = 0;
        this->prob = 0;
        this->reset = true;
    }

    ~vertex() {
        free(this->state);
    }

    vec_entry_t *state;
    struct llc *edges;
    size_t nedges;
    struct llp *parents;
    size_t nparents;

    bool visited;
    double rate;
    double prob;
    vector<double> rewards;

    vector<double> exp;
    vector<double> desc;
    vector<double> no_defect;

    size_t vertex_index;
    int integer;
    bool boolean;
    bool reset;
} vertex_t;

typedef struct llc {
    struct vertex *child;
    struct llp *llp;
    double weight;
} llc_t;

typedef struct llp {
    struct vertex *parent;
    struct llc *llc;
    struct llp *prev;
    struct llp *next;
} llp_t;

vertex_t *vertex_init(vec_entry_t *state, vector<double> rewards, size_t state_length);

void vertex_destroy(vertex_t *vertex);

void vertex_add_edge(vertex_t *from, vertex_t *to, double weight);

void graph_free(vertex_t *graph);

vertex_t *get_abs_vertex(vertex_t *graph);

queue<vertex_t *> topological_queue(vertex_t *graph);

queue<vertex_t *> enqueue_vertices(vertex_t *graph);

int label_vertex_index(size_t *largest_index, vertex_t *graph);

struct graph_info {
    size_t vertices;
    size_t edges;
};

struct pdf_values {
    long double lambda;
    long double k;
    size_t n;
    int c;
};

struct vertex_pdf {
    long double defect_prob;
    vector<struct pdf_values> *parts;
};

struct graph_info get_graph_info(vertex_t *graph);

cov_exp_return mph_cov_exp_all(vertex_t *graph, size_t m);

int reward_transform(vertex_t *graph, double (*reward_func)(vertex_t *));

double fold(vertex_t *graph, double (*vertex_func)(vertex_t *));

void set_graph_rewards(vertex_t *graph, vector<double> (*set_rewards_func)(vector<double>));

void reduce_graph(vertex_t *graph);

void calculate_prob(vertex_t *graph, size_t *size, double **probs);

void calculate_var(vertex_t *graph, size_t *size, double **vars);

void pdf(vertex_t *graph, double (*reward_func)(vertex_t *), struct vertex_pdf *out_vertex_pdfs);

vertex_t *generate_state_space(
        size_t state_length,
        vector<pair<double, vector<size_t> > >(*visit_function)(vector<size_t>),
        vector<pair<double, vector<size_t> > >(*initial_states)(void),
        vector<double>(*rewards)(vector<size_t>)
);

int gen_kingman_graph(vertex_t **graph, size_t n, size_t m);

int gen_kingman_graph_rw(vertex_t **graph, size_t n, size_t rw);

/*
// TODO: Not actual
typedef struct expanding_array {
    vector<char> vector;
} expanding_array_t;

expanding_array_t *expanding_array_create(size_t entry_size);

expanding_array_t *expanding_array_create_from_array(size_t entry_size);

void expanding_array_destroy(expanding_array_t *expanding_array);

size_t expanding_array_length(const expanding_array_t *expanding_array);

size_t expanding_array_entry_size(const expanding_array_t *expanding_array);

void *expanding_array_entries(const expanding_array_t *expanding_array);
 */

struct ptd_vertex_linked_list_item {
    struct ptd_vertex *vertex;
    struct ptd_vertex_linked_list_item *previous;
    struct ptd_vertex_linked_list_item *next;
};

struct ptd_vertex_linked_list {
    size_t vertices_count;
    struct ptd_vertex_linked_list_item *first;
    struct ptd_vertex_linked_list_item *tail;
};

struct ptd_vertex {
    struct ptd_graph *graph;
    size_t edges_length;
    size_t edges_limit;
    struct ptd_edge *edges;
    long double rate;

    size_t index;
    bool visited;
    bool reset;

    struct ptd_vertex_linked_list_item *item_in_graph_list;
    vec_entry_t *state;
    void *data;
};

struct ptd_edge {
    struct ptd_vertex *to;
    long double weight;
};

struct ptd_graph {
    size_t state_length;
    struct ptd_vertex *start_vertex;
    struct ptd_vertex_linked_list *vertices_list;

    bool is_indexed;
};

struct ptd_vertex *ptd_vertex_create_state(struct ptd_graph *graph, vec_entry_t *state);

struct ptd_graph *ptd_graph_create(size_t state_length);

void ptd_graph_vertices_destroy(struct ptd_graph *graph);

void ptd_graph_destroy(struct ptd_graph *graph);

struct ptd_vertex *ptd_vertex_create(struct ptd_graph *graph);

void ptd_vertex_destroy(struct ptd_vertex *vertex);

int ptd_visit_vertices(struct ptd_graph *graph, int (*visit_func)(struct ptd_vertex *), bool include_start);

int ptd_add_edge(struct ptd_vertex *from, struct ptd_vertex *to, long double weight);

int ptd_remove_edge(struct ptd_vertex *from, struct ptd_vertex *to);

int ptd_reward_transform(struct ptd_graph *graph, double (*reward_func)(const struct ptd_vertex *));

queue<struct ptd_vertex *> ptd_enqueue_vertices(struct ptd_graph *graph);

// TODO: This does not belong here
int ptd_label_vertices(struct ptd_graph *graph);

// TODO: static
typedef struct strongly_connected_component ptd_strongly_connected_component_t;

typedef struct strongly_connected_component {
    size_t internal_vertices_length;
    struct ptd_vertex **internal_vertices;
    size_t external_sccs_length;
    ptd_strongly_connected_component_t **external_sccs;
    size_t external_vertices_length;
    struct ptd_vertex **external_vertices;

    bool visited;
    size_t index;
} ptd_strongly_connected_component_t;

typedef struct strongly_connected_components {
    size_t components_length;
    ptd_strongly_connected_component_t **components;
} ptd_strongly_connected_components_t;

/*typedef struct scc_vertex {
    ptd_strongly_connected_component_t *scc;
    size_t scc_edges_length;
    scc_vertex **scc_edges;
    size_t local_edges_length;
    struct ptd_vertex **local_edges;
    size_t topo;
    bool visited;
} ptd_ordered_scc_t;

typedef struct {
    size_t ordered_components_length;
    ptd_ordered_scc_t **ordered_components;
} ptd_ordered_sccs_t;*/

ptd_strongly_connected_components_t *
ptd_find_strongly_connected_components(struct ptd_graph *graph);

void ptd_strongly_connected_components_destroy(ptd_strongly_connected_components_t *sccs);

int
ptd_order_strongly_connected_components(
        ptd_strongly_connected_components_t *sccs
);

typedef struct ptd_vertex_group ptd_vertex_group_t;

struct ptd_vertex_group {
    struct ptd_graph *graph;
    size_t edges_length;
    size_t edges_limit;
    struct ptd_edge *edges;
    long double rate;

    size_t index;
    bool visited;
    bool reset;

    vec_entry_t *state;
};

struct ptd_graph *
ptd_convert_strongly_connected_components_to_group(struct ptd_graph *graph, ptd_strongly_connected_components_t *sccs);

typedef struct ptd_avl_tree {
    void *root;
    size_t vec_length;
} ptd_avl_tree_t;

ptd_avl_tree_t *ptd_avl_tree_create(size_t vec_length);

void ptd_avl_tree_vertex_destroy(ptd_avl_tree_t *avl_tree);

void ptd_avl_tree_vertex_destroy_free(ptd_avl_tree_t *avl_tree);

void ptd_avl_tree_edge_destroy(ptd_avl_tree_t *avl_tree);

int ptd_avl_tree_vertex_insert(ptd_avl_tree_t *avl_tree, const vec_entry_t *key, const struct ptd_vertex *vertex);

struct ptd_vertex *ptd_avl_tree_vertex_find(const ptd_avl_tree_t *avl_tree, const vec_entry_t *key);

int ptd_avl_tree_edge_insert_or_increment(ptd_avl_tree_t *avl_tree, const vec_entry_t *key, struct ptd_vertex *vertex,
                                          long double weight);

int ptd_avl_tree_edge_remove(ptd_avl_tree_t *avl_tree, const vec_entry_t *key);

struct ptd_edge *ptd_avl_tree_edge_find(const ptd_avl_tree_t *avl_tree, const vec_entry_t *key);

//TODO: remove?
size_t ptd_avl_tree_max_depth(void *avl_vec_vertex);

typedef struct {
    size_t length;
    long double *initial_probability_vector;
    long double **sub_intensity_matrix;
    struct ptd_vertex **vertices;
    size_t memory_allocated;
} ptd_phase_type_distribution_t;

ptd_phase_type_distribution_t *ptd_graph_as_phase_type_distribution(struct ptd_graph *graph);

ptd_phase_type_distribution_t *
ptd_find_local_matrix(ptd_strongly_connected_components_t *in);

void ptd_phase_type_distribution_destroy(ptd_phase_type_distribution_t *ptd);

int ptd_index_topological(struct ptd_graph *graph);

int ptd_index_invert(struct ptd_graph *graph);

int ptd_vertex_to_s(struct ptd_vertex *vertex, char *buffer, size_t buffer_length);

/*
 * Models
 */

struct ptd_graph *ptd_model_kingman(size_t n);

struct ptd_graph *ptd_model_two_island_two_loci_recomb(
        size_t n1,
        size_t n2,
        size_t effective_population_size1,
        size_t effective_population_size2,
        double migration_rate1,
        double migration_rate2,
        double recombination_rate
);

/*
 * Visit algorithms
 */

int ptd_visit_alloc(struct ptd_graph *graph, size_t size, char *value);

int ptd_visit_assign_probability(struct ptd_graph *graph);

int ptd_visit_assign_expectation(struct ptd_graph *graph, double (*reward)(struct ptd_vertex *));

long double ptd_visit_reduce_sum_long_double(struct ptd_graph *graph);


/*
 * Properties
 */

int ptd_expected_value(long double *expected, struct ptd_graph *graph, double (*reward)(struct ptd_vertex *));

int ptd_covariance(
        long double *covariance, struct ptd_graph *graph,
        double (*reward_1)(struct ptd_vertex *), double (*reward_2)(struct ptd_vertex *)
);


typedef struct {
    struct ptd_vertex *vertex;
    double multiplier;
    bool external;
} ptd_desc_multiplier_t;

typedef struct {
    size_t *desc_length;
    ptd_desc_multiplier_t **desc_multipliers;
    ptd_strongly_connected_components_t *ordered;
} ptd_desc_multipliers_t;

double ptd_circular_exp(struct ptd_graph *graph, double (*reward)(struct ptd_vertex *));

ptd_desc_multipliers_t *ptd_cyclic_descendant_multipliers(struct ptd_graph *graph);
double *ptd_cyclic_desc(struct ptd_graph *graph, ptd_desc_multipliers_t *multipliers, double (*reward)(struct ptd_vertex *));

#endif //PTDALGORITHMS_PTD_H
