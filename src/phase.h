#ifndef PTDALGORITHMS_PHASE_H
#define PTDALGORITHMS_PHASE_H


#include <vector>

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

static size_t id = 0;

typedef struct vertex {
    vertex(vec_entry_t *state, size_t state_length) {
        this->vertex_index = id;
        id++;
        this->state = state;
        this->rewards = vector<double>(state, state + state_length);
    }

    ~vertex() {
        free(this->state);
    }

    vec_entry_t *state;
    vector<struct vertex*> children = vector<struct vertex*>();
    vector<double> weights = vector<double>();

    vector<struct vertex*> parents = vector<struct vertex*>();
    vector<double> weights_parent = vector<double>();

    bool visited = false;
    double rate = 0;
    double prob = 0;
    vector<double> rewards;

    vector<double> exp;
    vector<double> desc;

    size_t vertex_index;
    size_t reset_int = 0;
} vertex_t;


cov_exp_return mph_cov_exp_all(vertex_t *graph, size_t m);
void graph_free(vertex_t *graph);
int gen_kingman_graph(vertex_t **graph, size_t n, size_t m);

#endif //PTDALGORITHMS_PHASE_H
