#include <iostream>
#include <vector>
#include <cstring>
#include <stack>
#include <queue>
#include "phase.h"

using namespace std;


void print_graph_list(FILE *stream, vertex_t *graph,
                      bool indexed,
                      size_t vec_length, size_t vec_spacing);

vertex_t *vertex_init(vec_entry_t *state, vector<double> rewards, size_t state_length) {
    vertex_t *vertex = new vertex_t(state, rewards, state_length);
    llp_t *llp_dummy = (llp_t *) calloc(1, sizeof(llp_t));

    vertex->edges = NULL;
    vertex->parents = llp_dummy;
    vertex->nedges = 0;
    vertex->nparents = 0;

    return vertex;
}

void ll_destroy(llp_t *llp) {
    if (llp->next != NULL) {
        ll_destroy(llp->next);
    }

    fprintf(stderr, "ll_destroy freeing %p\n", (void*)llp);
    free(llp);
}

void ll_remove(llp_t *llp) {
    llp_t *prev = llp->prev;
    llp_t *next = llp->next;

    if (prev != NULL) {
        prev->next = next;
    }

    if (next != NULL) {
        next->prev = prev;
    }


    fprintf(stderr, "ll_remove prev, next is %p, %p\n",
            (void*)llp->prev,(void*)llp->next);

    fprintf(stderr, "ll_remove freeing %p\n", (void*)llp);

    free(llp);
}

/*
 * Assumes that my children will have removed me from their
 * parents list independently of this call
 */
void vertex_destroy(vertex_t *vertex) {
    ll_destroy(vertex->parents);
    fprintf(stderr, "Free alloc %p\n", vertex->edges);
    free(vertex->edges);
    delete (vertex);
}

/*
 * Destroys all parent links from my children
 */
void vertex_destroy_parents(vertex_t *vertex) {
    fprintf(stderr, "I am vertex %p, I will destroy the links from my children to me. I assume that I am no-ones child\n", (void*)vertex);

    for (size_t i = 0; i < vertex->nedges; ++i) {
        if (vertex->edges[i].child == NULL) {
            DIE_ERROR(1, "NULL child\n");
        }

        if (vertex->edges[i].weight <= 0.0001) {
            DIE_ERROR(1, "zero weight\n");
        }

        if (vertex->edges[i].llp == NULL) {
            DIE_ERROR(1, "NULL parent edge\n");
        }


        fprintf(stderr, "My child %zu is %p, the llp is %p. I will remove that as I am no longer its parent\n",
                i,(void*)(vertex->edges[i].child), (void*)vertex->edges[i].llp);

        ll_remove(vertex->edges[i].llp);
    }

    if (vertex->parents == NULL) {
        DIE_ERROR(1, "No parent link\n");
    }
}

inline void add_edge(vertex_t *from, vertex_t *to, size_t index, double weight) {
    llc_t *entry = &(from->edges[index]);

    entry->weight = weight;
    entry->child = to;
    llp_t *llp = (llp_t *) malloc(sizeof(llp_t));
    llp->llc = entry;
    llp->parent = from;
    entry->llp = llp;

    llp_t *next = to->parents->next;
    llp->next = next;
    llp->prev = to->parents;

    if (next != NULL) {
        next->prev = llp;
    }

    to->parents->next = llp;

    from->nedges++;
    to->nparents++;
    from->rate += weight;
}

void vertex_add_edge(vertex_t *from, vertex_t *to, double weight) {
    from->edges = (llc_t *) reallocarray(from->edges, from->nedges + 1, sizeof(llc_t));
    fprintf(stderr, "New alloc %p\n", from->edges);

    for (size_t i = 0; i < from->nedges; ++i) {
        from->edges[i].llp->llc = &(from->edges[i]);
    }

    add_edge(from, to, from->nedges, weight);
}

static void print_vector_spacing(FILE *stream, vec_entry_t *v, size_t nmemb, size_t spacing) {
    if (v == NULL) {
        fprintf(stream, "(NULL)");
        return;
    }

    fprintf(stream, "(");

    for (size_t i = 0; i < nmemb; i++) {
        if (i % spacing == 0) {
            fprintf(stream, "-");
        }

        fprintf(stream, "%zu", v[i]);
    }
    fprintf(stream, ")");
}

void _reset_graph_visited(vertex_t *vertex, size_t reset_int) {
    if (vertex->reset_int == reset_int) {
        return;
    }

    vertex->reset_int = reset_int;

    for (size_t i = 0; i < vertex->nedges; ++i) {
        vertex_t *child = vertex->edges[i].child;
        _reset_graph_visited(child, reset_int);
    }

    vertex->visited = false;
}

static void reset_graph_visited(vertex_t *vertex) {
    _reset_graph_visited(vertex, vertex->reset_int + 1);
}

typedef struct avl_vec_vertex {
    struct avl_vec_vertex *left;
    struct avl_vec_vertex *right;
    struct avl_vec_vertex *parent;
    signed short balance;
    vec_entry_t *key;
    vertex_t *entry;
} avl_vec_vertex_t;

static inline int radix_cmp(const vec_entry_t *a, const vec_entry_t *b,
                            const size_t vec_length) {
    return (memcmp(a, b, sizeof(vec_entry_t) * vec_length));
}

/* Example:
*     A            A            A
*   B   (left)    B  (right)   D
* C       ->    D      ->    C   B
*   D         C
* In this case:
*  C: child
*  B: parent
*  D: child_right
*/
avl_vec_vertex_t *rotate_left_right(avl_vec_vertex_t *parent, avl_vec_vertex_t *child) {
    avl_vec_vertex_t *child_right_left, *child_right_right;
    avl_vec_vertex_t *child_right = child->right;
    child_right_left = child_right->left;
    child->right = child_right_left;

    if (child_right_left != nullptr) {
        child_right_left->parent = child;
    }

    child_right->left = child;
    child->parent = child_right;
    child_right_right = child_right->right;
    parent->left = child_right_right;

    if (child_right_right != nullptr) {
        child_right_right->parent = parent;
    }

    child_right->right = parent;
    parent->parent = child_right;

    if (child_right->balance > 0) {
        parent->balance = -1;
        child->balance = 0;
    } else if (child_right->balance == 0) {
        parent->balance = 0;
        child->balance = 0;
    } else {
        parent->balance = 0;
        child->balance = +1;
    }

    child_right->balance = 0;

    return child_right;
}

/* Example:
*  A          A            A
*   B  (right)  B   (left)   D
*     C   ->      D    -> B    C
*   D               C
* In this case:
*  C: child
*  B: parent
*  D: child_left
*/
avl_vec_vertex_t *rotate_right_left(avl_vec_vertex_t *parent, avl_vec_vertex_t *child) {
    avl_vec_vertex_t *child_left_right, *child_left_left;
    avl_vec_vertex_t *child_left = child->left;

    child_left_right = child_left->right;

    child->left = child_left_right;

    if (child_left_right != nullptr) {
        child_left_right->parent = child;
    }

    child_left->right = child;

    child->parent = child_left;
    child_left_left = child_left->left;
    parent->right = child_left_left;

    if (child_left_left != nullptr) {
        child_left_left->parent = parent;
    }

    child_left->left = parent;
    parent->parent = child_left;

    if (child_left->balance > 0) {
        parent->balance = -1;
        child->balance = 0;
    } else if (child_left->balance == 0) {
        parent->balance = 0;
        child->balance = 0;
    } else {
        parent->balance = 0;
        child->balance = 1;
    }

    child_left->balance = 0;

    return child_left;
}

/*
* Example:
*  A              B
*    B   (left) A   C
*      C   ->
*/
avl_vec_vertex_t *rotate_left(avl_vec_vertex_t *parent, avl_vec_vertex_t *child) {
    avl_vec_vertex_t *child_left;

    child_left = child->left;
    parent->right = child_left;

    if (child_left != nullptr) {
        child_left->parent = parent;
    }

    child->left = parent;
    parent->parent = child;

    if (child->balance == 0) {
        parent->balance = 1;
        child->balance = -1;
    } else {
        parent->balance = 0;
        child->balance = 0;
    }

    return child;
}

/*
* Example:
*      A            B
*    B    (right) C   A
*  C        ->
*/
avl_vec_vertex_t *rotate_right(avl_vec_vertex_t *parent, avl_vec_vertex_t *child) {
    avl_vec_vertex_t *child_right;

    child_right = child->right;
    parent->left = child_right;

    if (child_right != nullptr) {
        child_right->parent = parent;
    }

    child->right = parent;
    parent->parent = child;

    if (child->balance == 0) {
        parent->balance = +1;
        child->balance = -1;
    } else {
        parent->balance = 0;
        child->balance = 0;
    }

    return child;
}

int avl_vec_vertex_create(avl_vec_vertex_t **vertex, vec_entry_t *key, vertex_t *entry, avl_vec_vertex_t *parent) {
    if ((*vertex = (avl_vec_vertex_t *) malloc(sizeof(avl_vec_vertex_t))) == nullptr) {
        return 1;
    }

    (*vertex)->key = key;
    (*vertex)->entry = entry;
    (*vertex)->left = nullptr;
    (*vertex)->right = nullptr;
    (*vertex)->parent = parent;
    (*vertex)->balance = 0;

    return 0;
}

void avl_vec_vertex_destroy(avl_vec_vertex_t *vertex) {
    if (vertex == nullptr) {
        return;
    }

    avl_vec_vertex_destroy(vertex->left);
    avl_vec_vertex_destroy(vertex->right);

    free(vertex);
}

static void avl_free(avl_vec_vertex_t *vertex) {
    if (vertex == nullptr) {
        return;
    }

    avl_free(vertex->left);
    avl_free(vertex->right);
    free(vertex);
}

const avl_vec_vertex_t *avl_vec_find(const avl_vec_vertex_t *rootptr, const vec_entry_t *key, const size_t vec_length) {
    if (rootptr == nullptr) {
        return nullptr;
    }

    const avl_vec_vertex_t *vertex = rootptr;

    while (true) {
        int res = radix_cmp(key, vertex->key, vec_length);
        if (res < 0) {
            if (vertex->left == nullptr) {
                return nullptr;
            } else {
                vertex = vertex->left;
            }
        } else if (res > 0) {
            if (vertex->right == nullptr) {
                return nullptr;
            } else {
                vertex = vertex->right;
            }
        } else {
            return vertex;
        }
    }
}

int find_or_insert_vec(avl_vec_vertex_t **out, avl_vec_vertex_t *rootptr, vec_entry_t *key, vertex_t *entry,
                       const size_t vec_length) {
    if (avl_vec_vertex_create(out, key, entry, nullptr)) {
        return -1;
    }

    if (rootptr == nullptr) {
        return 1;
    }

    avl_vec_vertex_t *vertex = rootptr;

    while (true) {
        int res = radix_cmp(key, vertex->key, vec_length);
        if (res < 0) {
            if (vertex->left == nullptr) {
                vertex->left = *out;
                break;
            } else {
                vertex = vertex->left;
            }
        } else if (res > 0) {
            if (vertex->right == nullptr) {
                vertex->right = *out;
                break;
            } else {
                vertex = vertex->right;
            }
        } else {
            *out = vertex;
        }
    }

    (*out)->parent = vertex;

    return 0;
}

int avl_vec_insert(avl_vec_vertex_t **root, vec_entry_t *key, vertex_t *entry, const size_t vec_length) {
    avl_vec_vertex_t *child;

    if (*root == nullptr) {
        if (avl_vec_vertex_create(root, key, entry, nullptr)) {
            return 1;
        }

        return 0;
    }

    int res = find_or_insert_vec(&child, *root, key, entry, vec_length);

    if (res == -1) {
        return 1;
    }

    if (res == 0) {
        return 0;
    }

    avl_vec_vertex_t *pivot, *rotated_parent;

    for (avl_vec_vertex_t *parent = child->parent; parent != nullptr; parent = child->parent) {
        if (child == parent->right) {
            if (parent->balance > 0) {
                pivot = parent->parent;

                if (child->balance < 0) {
                    rotated_parent = rotate_right_left(parent, child);
                } else {
                    rotated_parent = rotate_left(parent, child);
                }
            } else {
                if (parent->balance < 0) {
                    parent->balance = 0;

                    return 0;
                }

                parent->balance = 1;
                child = parent;

                continue;
            }
        } else {
            if (parent->balance < 0) {
                pivot = parent->parent;

                if (child->balance > 0) {
                    rotated_parent = rotate_left_right(parent, child);
                } else {
                    rotated_parent = rotate_right(parent, child);
                }
            } else {
                if (parent->balance > 0) {
                    parent->balance = 0;

                    return 0;
                }

                parent->balance = -1;
                child = parent;
                continue;
            }
        }

        rotated_parent->parent = pivot;

        if (pivot != nullptr) {
            if (parent == pivot->left) {
                pivot->left = rotated_parent;
            } else {
                pivot->right = rotated_parent;
            }

            return 0;
        } else {
            *root = rotated_parent;
        }
    }

    return 0;
}

static size_t avl_vec_get_size(avl_vec_vertex_t *vertex) {
    if (vertex == nullptr) {
        return 0;
    }

    return 1 + avl_vec_get_size(vertex->left) + avl_vec_get_size(vertex->right);
}

static queue<vertex_t *> enqueue_vertices(vertex_t *graph) {
    queue<vertex_t *> ret;
    queue<vertex_t *> queue;
    reset_graph_visited(graph);

    queue.push(graph);

    while (!queue.empty()) {
        vertex_t *vertex = queue.front();
        queue.pop();

        if (vertex->visited) {
            continue;
        }

        vertex->visited = true;
        ret.push(vertex);

        for (size_t i = 0; i < vertex->nedges; ++i) {
            vertex_t *child = vertex->edges[i].child;
            queue.push(child);
        }
    }

    return ret;
}

static vertex_t *get_abs_vertex(vertex_t *graph) {
    queue<vertex_t *> queue = enqueue_vertices(graph);
    vertex_t *abs_vertex = nullptr;

    while (!queue.empty()) {
        vertex_t *vertex = queue.front();
        queue.pop();

        if (vertex->nedges == 0) {
            if (abs_vertex == nullptr) {
                abs_vertex = vertex;
            } else {
                DIE_ERROR(1, "Found multiple absorbing vertices\n");
            }
        }
    }

    if (abs_vertex == nullptr) {
        DIE_ERROR(1, "No absorbing vertex found!\n");
    }

    return abs_vertex;
}

static int kingman_visit_vertex(vertex_t **out_initial_vertex,
                                vec_entry_t *initial_state,
                                vertex_t *abs_vertex,
                                const size_t m) {
    avl_vec_vertex_t *bst = nullptr;

    queue<vertex_t *> vertices_to_visit;
    vertex_t *initial_vertex = vertex_init(initial_state,
                                           vector<double>(initial_state, initial_state + m),
                                           m);
    vertices_to_visit.push(initial_vertex);
    vec_entry_t *v = (vec_entry_t *) malloc(sizeof(vec_entry_t) * m);
    avl_vec_vertex_t *bst_vertex;

    size_t end = 0;

    while (!vertices_to_visit.empty()) {
        vertex_t *vertex = vertices_to_visit.front();
        vertices_to_visit.pop();
        memcpy(v, vertex->state, sizeof(vec_entry_t) * m);
        size_t n_remaining = 0;
        size_t start = -1;

        for (vec_entry_t i = 0; i < m; i++) {
            n_remaining += v[i];

            if (v[i] > 0) {
                end = i;

                if (start == (size_t) -1) {
                    start = i;
                }
            }
        }

        for (vec_entry_t i = start; i <= end; i++) {
            if (v[i] == 0) {
                continue;
            }

            for (vec_entry_t j = i; j <= end; j++) {
                if (((i == j && v[i] >= 2) || (i != j && v[i] > 0 && v[j] > 0))) {
                    vertex->nedges++;
                }
            }

            vertex->edges = (llc_t *) calloc(vertex->nedges, sizeof(llc_t));

            size_t pos = 0;

            for (vec_entry_t j = i; j <= end; j++) {
                if (((i == j && v[i] >= 2) || (i != j && v[i] > 0 && v[j] > 0))) {
                    double t = i == j ? v[i] * (v[i] - 1) / 2 : v[i] * v[j];

                    const size_t inc_pos = min((i + j + 2) - 1, m - 1);
                    v[i]--;
                    v[j]--;
                    v[inc_pos]++;

                    vertex_t *new_vertex;
                    const bool only_tail = (v[m - 1] == n_remaining - 1);

                    if (only_tail) {
                        new_vertex = abs_vertex;
                    } else {
                        bst_vertex = (avl_vec_vertex_t *) avl_vec_find(bst, v, m);

                        if (bst_vertex == nullptr) {
                            vec_entry_t *new_state = (vec_entry_t *) malloc(sizeof(vec_entry_t) * m);
                            memcpy(new_state, v, sizeof(vec_entry_t) * m);
                            vertex_t *to = vertex_init(new_state,
                                                       vector<double>(new_state, new_state + m),
                                                       m);

                            avl_vec_insert(&bst, new_state, to, m);
                            vertices_to_visit.push(to);
                            new_vertex = to;
                        } else {
                            new_vertex = bst_vertex->entry;
                        }
                    }

                    v[i]++;
                    v[j]++;
                    v[inc_pos]--;

                    add_edge(vertex, new_vertex, pos, t);
                    pos++;
                }
            }
        }
    }

    free(v);
    *out_initial_vertex = initial_vertex;
    avl_free(bst);
    return 0;
}

int gen_kingman_graph(vertex_t **graph, size_t n, size_t m) {
    if (m < n) {
        m += 1;
    }

    vec_entry_t *initial = (vec_entry_t *) calloc(m, sizeof(vec_entry_t));
    initial[0] = n;

    vec_entry_t *mrca = (vec_entry_t *) calloc(m, sizeof(vec_entry_t));
    mrca[m - 1] = 1;

    vertex_t *absorbing_vertex = vertex_init(mrca,
                                             vector<double>(mrca, mrca + m),
                                             m);

    vertex_t *state_graph;

    kingman_visit_vertex(&state_graph,
                         initial, absorbing_vertex,
                         m);

    vec_entry_t *start_state = (vec_entry_t *) calloc(m, sizeof(vec_entry_t));

    vertex_t *start = vertex_init(start_state, vector<double>(start_state, start_state + m), m);

    start->nedges = 1;
    start->edges = (llc_t *) calloc(1, sizeof(llc_t));
    add_edge(start, state_graph, 0, 1);

    *graph = start;

    return 0;
}

static void _print_graph_list(FILE *stream, vertex_t *vertex,
                              bool indexed,
                              size_t vec_length, size_t vec_spacing) {
    if (vertex->visited) {
        return;
    }

    vertex->visited = true;

    fprintf(stream, "vertex: ");
    print_vector_spacing(stream, vertex->state,
                         vec_length, vec_spacing);
    if (indexed) {
        fprintf(stream, " (%zu)", vertex->vertex_index);
    }
    fprintf(stream, " %p", (void*)vertex);
    fprintf(stream, ":\n");

    for (size_t i = 0; i < vertex->nedges; ++i) {
        llc_t child_edge = vertex->edges[i];
        vertex_t *child = child_edge.child;
        fprintf(stream, "\t");

        fprintf(stream, " %p ", (void*)child);
        fprintf(stream, "(%p) ", (void*)&child_edge);
        fprintf(stream, " (w%f) ", child_edge.weight);
        //print_vector_spacing(stream, child->state,
        //                     vec_length, vec_spacing);

        if (indexed) {
            fprintf(stream, " (%zu)", vertex->vertex_index);
        }

        //fprintf(stream, "\n");
    }


    fprintf(stream, "\n");

    llp_t *parent = vertex->parents->next;

    while (parent != NULL) {
        fprintf(stream, "\t");

        fprintf(stream, "P %p ", (void*)parent->parent);
        fprintf(stream, "(%p) ", (void*)&parent);
        fprintf(stream, " (llc%p) ", (void*)parent->llc);
        //print_vector_spacing(stream, parent->parent->state,
         //                    vec_length, vec_spacing);

        if (indexed) {
            fprintf(stream, " (%zu)", vertex->vertex_index);
        }

        //fprintf(stream, "\n");
        parent = parent->next;
    }

    fprintf(stream, "\n");

    for (size_t i = 0; i < vertex->nedges; ++i) {
        vertex_t *child = vertex->edges[i].child;
        _print_graph_list(stream, child,
                          indexed,
                          vec_length, vec_spacing);
    }
}

void print_graph_list(FILE *stream, vertex_t *graph,
                      bool indexed,
                      size_t vec_length, size_t vec_spacing) {
    reset_graph_visited(graph);
    fprintf(stream, "-- Graph list --\n");
    _print_graph_list(stream, graph, indexed, vec_length, vec_spacing);
    fflush(stream);
}

void mph_cov_assign_vertex_all(vertex_t *vertex, size_t m) {
    if (vertex->visited) {
        return;
    }

    vertex->visited = true;

    if (vertex->parents->next == NULL) {
        // Starting vertex
        vertex->prob = 1.0f;
    } else {
        vertex->prob = 0.0f;
    }

    llp_t *parent = vertex->parents->next;

    while (parent != NULL) {
        mph_cov_assign_vertex_all(parent->parent, m);

        vertex->prob += parent->llc->weight / parent->parent->rate * parent->parent->prob;

        parent = parent->next;
    }

    vertex->exp = vector<double>();

    for (size_t j = 0; j < m; ++j) {
        if (vertex->rate != 0) {
            vertex->exp.push_back(vertex->prob * vertex->rewards[j] / vertex->rate);
        } else {
            vertex->exp.push_back(0);
        }
    }
}

void mph_cov_assign_desc_all(vertex_t *vertex, size_t m) {
    if (vertex->visited) {
        return;
    }

    vertex->visited = true;

    vertex->desc = vector<double>(m);

    for (size_t i = 0; i < vertex->nedges; ++i) {
        vertex_t *child = vertex->edges[i].child;
        mph_cov_assign_desc_all(child, m);
    }

    for (size_t i = 0; i < vertex->nedges; ++i) {
        llc_t child_edge = vertex->edges[i];

        for (size_t j = 0; j < m; ++j) {
            vertex->desc[j] += child_edge.weight / vertex->rate * child_edge.child->desc[j];
        }
    }

    vector<double> exp(m);

    for (size_t j = 0; j < m; ++j) {
        if (vertex->rate != 0) {
            exp[j] = vertex->rewards[j] / vertex->rate;
        } else {
            exp[j] = 0;
        }
    }

    for (size_t j = 0; j < m; ++j) {
        vertex->desc[j] += exp[j];
    }
}

double **cov;
double *expectation;

void _mph_cov_all(vertex_t *vertex, size_t m) {
    if (vertex->visited) {
        return;
    }

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            cov[i][j] += vertex->desc[i] * vertex->exp[j];
            cov[i][j] += vertex->desc[j] * vertex->exp[i];
        }
    }

    vertex->visited = true;
    for (size_t i = 0; i < vertex->nedges; ++i) {
        vertex_t *child = vertex->edges[i].child;
        _mph_cov_all(child, m);
    }
}

cov_exp_return mph_cov_exp_all(vertex_t *graph, size_t m) {
    vertex_t *abs = get_abs_vertex(graph);
    reset_graph_visited(graph);
    mph_cov_assign_vertex_all(abs, m);
    reset_graph_visited(graph);
    mph_cov_assign_desc_all(graph, m);

    cov = (double **) calloc(m, sizeof(double *));

    for (size_t i = 0; i < m; ++i) {
        cov[i] = (double *) calloc(m, sizeof(double));
    }

    reset_graph_visited(graph);
    _mph_cov_all(graph, m);

    expectation = (double *) calloc(m, sizeof(double));

    for (size_t i = 0; i < m; ++i) {
        expectation[i] = graph->desc[i];

        for (size_t j = 0; j <= i; ++j) {
            cov[i][j] -= graph->desc[i] *
                         graph->desc[j];
        }
    }

    return (cov_exp_return) {.cov = &cov, .exp = &expectation};
}

void graph_free(vertex_t *graph) {
    queue<vertex_t *> to_free = enqueue_vertices(graph);

    while (!to_free.empty()) {
        vertex_t *vertex = to_free.front();
        to_free.pop();
        vertex_destroy(vertex);
    }
}

double calculate_rate(vertex_t *vertex) {
    double rate = 0;

    for (size_t i = 0; i < vertex->nedges; ++i) {
        vertex_t *child = vertex->edges[i].child;
        llc_t child_edge = vertex->edges[i];
        rate += child_edge.weight;
    }

    return rate;
}

size_t times = 0;
void ensure_valid_graph(vertex_t *graph) {
    return;
    queue<vertex_t *> queue = enqueue_vertices(graph);

    while (!queue.empty()) {
        vertex_t *vertex = queue.front();
        queue.pop();

        for (size_t i = 0; i < vertex->nedges; ++i) {
            if (vertex->edges[i].llp->llc != &(vertex->edges[i])) {
                DIE_ERROR(1, "Invalid link\n");
            }
        }

        llp_t *parent = vertex->parents->next;

        while(parent != NULL) {
            if (times > 6) {
                2+2;
                //print_graph_list(stderr, graph, false, 1, 1);
                //exit(1);
            }
            times++;
            if (parent->llc->llp != parent) {
                DIE_ERROR(1, "Invalid link\n");
            }

            parent = parent->next;
        }
    }
}

int reward_transform(vertex_t *graph, double (*reward_func)(vertex_t *)) {
    queue<vertex_t *> queue = enqueue_vertices(graph);

    while (!queue.empty()) {
        vertex_t *vertex = queue.front();
        queue.pop();

        if (vertex->nedges == 0) {
            // Absorbing vertex
            continue;
        }

        if (vertex->nparents == 0) {
            // Starting vertex
            continue;
        }
        size_t times = 0;

        double reward = reward_func(vertex);

        if (reward == 0) {
            fprintf(stderr, "==Vertex=== %p\n", (void*)vertex);
            vertex_destroy_parents(vertex);

            llc_t *children = vertex->edges;
            size_t nchildren = vertex->nedges;
            vertex->edges = NULL;
            vertex->nedges = 0;

            // Take all my edges and add to my parent instead.
            llp_t *parent_edge = vertex->parents->next;

            while (parent_edge != NULL) {
                fprintf(stderr, "Pedged parent: %p\n", (void*)parent_edge->parent);
                vertex_t *parent = parent_edge->parent;
                llc_t *parent_old_children = parent->edges;
                size_t parent_nchildren = parent->nedges;
                ensure_valid_graph(parent);

                llc_t *new_parent_children = (llc_t *) calloc(
                        parent_nchildren + nchildren,
                        sizeof(llc_t)
                );

                parent->edges = new_parent_children;

                double parent_weight = parent_edge->llc->weight;

                size_t i = 0, j = 0, k;

                for (k = 0; i < parent_nchildren || j < nchildren;) {
                    if (times > 4) {
                        2+2;
                        //exit(0);
                    }
                    times++;


                    if (j < nchildren &&
                        children[j].child == parent) {
                        double prob = children[j].weight / vertex->rate;
                        parent->rate -= parent_weight * prob;
                        j++;
                        continue;
                    }

                    if (i < parent_nchildren &&
                        parent_old_children[i].child == vertex) {
                        i++;
                        continue;
                    }

                    if (i >= parent_nchildren) {
                        double prob = children[j].weight / vertex->rate;
                        add_edge(parent, children[j].child, k, prob * parent_weight);
                        parent->rate -= prob * parent_weight;
                        j++;
                    } else if (j >= nchildren ||
                               parent_old_children[i].child < children[j].child) {
                        new_parent_children[k] = parent_old_children[i];
                        new_parent_children[k].llp->llc = &(new_parent_children[k]);
                        i++;
                    } else if (parent_old_children[i].child > children[j].child) {
                        // TODO we must remove the old parents
                        //      ll_remove(children[j].llp);
                        double prob = children[j].weight / vertex->rate;
                        add_edge(parent, children[j].child, k, prob * parent_weight);
                        parent->rate -= prob * parent_weight;
                        j++;
                    } else {
                        // ==
                        double prob = children[j].weight / vertex->rate;
                        new_parent_children[k] = parent_old_children[i];
                        new_parent_children[k].weight += prob * parent_weight;
                        new_parent_children[k].llp->llc = &(new_parent_children[k]);
                        i++;
                        j++;
                    }

                    for (size_t l = 0; l < k; ++l) {
                        if (new_parent_children[l].child == NULL) {
                            DIE_ERROR(1, "NULL child\n");
                        }

                        if (new_parent_children[l].weight <= 0.001) {
                            DIE_ERROR(1, "Zero weight\n");
                        }
                    }

                    // THE RATE IS WRONG!!!

                    k++;
                }

                parent->nedges = k;

                free(parent_old_children);
                parent_edge = parent_edge->next;
                ensure_valid_graph(graph);
            }

            free(children);
            vertex_destroy(vertex);
            ensure_valid_graph(graph);
        } else {
            for (size_t i = 0; i < vertex->nedges; ++i) {
                vertex->edges[i].weight /= reward;
            }

            vertex->rate /= reward;
        }
    }

    return 0;
}

/*
 * Also ensures that the absorbing vertex has index 0
 */
int label_vertex_index(size_t *largest_index, vertex_t *graph) {
    vertex_t *abs_vertex;
    abs_vertex = get_abs_vertex(graph);
    reset_graph_visited(graph);
    queue<vertex_t *> queue;
    size_t index = 0;

    // The absorbing vertex should have index 0
    queue.push(abs_vertex);
    queue.push(graph);

    while (!queue.empty()) {
        vertex_t *vertex = queue.front();
        queue.pop();

        if (vertex->visited) {
            continue;
        }

        vertex->visited = true;

        vertex->vertex_index = index++;

        for (size_t i = 0; i < vertex->nedges; ++i) {
            llc_t child = vertex->edges[i];
            queue.push(child.child);
        }
    }

    if (largest_index != nullptr) {
        *largest_index = index - 1;
    }

    return 0;
}


/*
 * Assumes visited state reset and indexed vertices.
 */
void insert_into_weight_mat(double **weights, vertex_t **vertices, vertex_t *vertex) {
    if (vertex->visited) {
        return;
    }

    vertex->visited = true;


    for (size_t i = 0; i < vertex->nedges; ++i) {
        llc_t child = vertex->edges[i];

        weights[vertex->vertex_index][child.child->vertex_index] = child.weight;
        weights[vertex->vertex_index][vertex->vertex_index] -= child.weight;
    }

    vertices[vertex->vertex_index] = vertex;


    for (size_t i = 0; i < vertex->nedges; ++i) {
        llc_t child = vertex->edges[i];
        insert_into_weight_mat(weights, vertices, child.child);
    }
}

int graph_as_mat(double ***weights, vertex_t ***vertices, size_t *out_size, vertex_t *graph) {
    size_t largest_index;
    label_vertex_index(&largest_index, graph);
    reset_graph_visited(graph);
    size_t size = largest_index + 1;
    *out_size = size;
    *weights = (double **) malloc(sizeof(double *) * size);

    for (size_t i = 0; i < size; ++i) {
        (*weights)[i] = (double *) calloc(size, sizeof(double));
    }

    *vertices = (vertex_t **) calloc(size, sizeof(vertex_t *));

    insert_into_weight_mat(*weights, *vertices, graph);
    return 0;
}

