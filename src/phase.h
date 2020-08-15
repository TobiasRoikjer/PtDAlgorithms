#include <utility>

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
    vertex(vec_entry_t *state, vector<double> rewards, size_t state_length) {
        this->vertex_index = id;
        id++;
        this->state = state;
        this->rewards = std::move(rewards);
    }

    ~vertex() {
        free(this->state);
    }

    typedef struct edge {
        struct vertex* vertex;
        double weight;

        bool operator==(const struct vertex* b) {
            return (vertex == b);
        }

        bool operator<(const struct vertex* b) {
            return (vertex->vertex_index < b->vertex_index);
        }

        bool operator<(const struct edge& b) const {
            return (vertex->vertex_index < b.vertex->vertex_index);
        }
    } edge_t;

    bool operator<(const struct edge& b) const {
        return (this->vertex_index < b.vertex->vertex_index);
    }

    vec_entry_t *state;
    vector<edge_t> children;
    vector<edge_t> parents;

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
int graph_as_mat(double ***weights, vertex_t ***vertices, size_t *out_size, vertex_t *graph);
int reward_transform(vertex_t *graph, double (*reward_func)(vertex_t*));
void print_graph_list(FILE *stream, vertex_t *graph, bool indexed, size_t vec_length, size_t vec_spacing);

#endif //PTDALGORITHMS_PHASE_H
