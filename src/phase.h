#ifndef PTDALGORITHMS_PHASE_H
#define PTDALGORITHMS_PHASE_H

#include <vector>
#include <utility>
#include <queue>

typedef struct {
    double ***cov;
    double **exp;
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
exit(error_code);                                  \
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
        this->rewards = std::move(rewards);
        this->edges = NULL;
        this->parents = NULL;
    }

    ~vertex() {
        free(this->state);
    }

    vec_entry_t *state;
    struct llc *edges;
    size_t nedges;
    struct llp *parents;
    size_t nparents;

    bool visited = false;
    double rate = 0;
    double prob = 0;
    vector<double> rewards;

    vector<double> exp;
    vector<double> desc;

    size_t vertex_index;
    size_t reset_int = 0;
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

queue<vertex_t *> enqueue_vertices(vertex_t *graph);
int label_vertex_index(size_t *largest_index, vertex_t *graph);

cov_exp_return mph_cov_exp_all(vertex_t *graph, size_t m);

void graph_free(vertex_t *graph);

int gen_kingman_graph(vertex_t **graph, size_t n, size_t m);

int reward_transform(vertex_t *graph, double (*reward_func)(vertex_t *));

#endif //PTDALGORITHMS_PHASE_H
