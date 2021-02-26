#ifndef PTDALGORITHMS_PHASE_H
#define PTDALGORITHMS_PHASE_H

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
/*char formatted[2048];                              \
                                                   \
snprintf(formatted,                                \
         sizeof(formatted),                        \
         message, ##__VA_ARGS__);                  \
                                                   \
fprintf(stderr, "%s", formatted); */                 \
} while(0)

using namespace std;

typedef size_t vec_entry_t;

struct llc;
struct llp;

static size_t id = 0;

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
    size_t integer;
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

#endif //PTDALGORITHMS_PHASE_H
