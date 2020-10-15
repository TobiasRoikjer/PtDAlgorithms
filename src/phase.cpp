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

static void verify_next(ll_t *ll) {
    return;
    ll_t *s = ll;
    vertex_t *seen[50];
    size_t seen_count = 0;

    while (s != NULL) {
        if (s->prev != NULL) {
            if (s->prev->next != s) {
                DIE_ERROR(1, "No <->\n");
            }
        }


        if (s->next != NULL) {
            if (s->next->prev != s) {
                DIE_ERROR(1, "No <->\n");
            }
        }


        for (size_t i = 0; i < seen_count; ++i) {
            if (seen[i] == s->vertex) {
                DIE_ERROR(1, "ALREADY SEEN \n");
            }
        }

        seen[seen_count] = s->vertex;
        seen_count++;
        s = s->next;
    }
}

static void ll_destroy(ll_t *ll) {
    if (ll != NULL) {
        ll_destroy(ll->next);
    }

    free(ll);
}

static void ll_remove(llc_t *ll) {
    ll->llp->parent->rate -= ll->weight;
    ll->prev->next = ll->next;

    if (ll->next != NULL) {
        ll->next->prev = ll->prev;
    }

    ll->llp->prev->next = ll->llp->next;

    if (ll->llp->next != NULL) {
        ll->llp->next->prev = ll->llp->prev;
    }

    verify_next((ll_t *) ll->prev);
    verify_next((ll_t *) ll->next);

    verify_next((ll_t *) ll->llp->prev);
    verify_next((ll_t *) ll->llp->next);

    free(ll->llp);
    free(ll);
}

struct ll_init ll_init(vertex_t *from, vertex_t *to, double weight) {
    llc_t *llc = (llc_t *) malloc(sizeof(llc_t));

    llc->child = to;
    llc->weight = weight;
    llc->next = NULL;
    llc->prev = NULL;

    llp_t *llp = (llp_t *) malloc(sizeof(llp_t));

    llp->parent = from;
    llp->next = NULL;
    llp->prev = NULL;

    llc->llp = llp;
    llp->llc = llc;

    return {.llc = llc, .llp = llp};
}

void ll_insert(vertex_t *from, vertex_t *to, double weight) {
    verify_next((ll_t*)from->edges);
    verify_next((ll_t*)to->parents);
    struct ll_init ll = ll_init(from, to, weight);
    llc_t *current = from->edges;

    while (true) {
        if (current->next == NULL
            || current->next->child > to) {
            llc_t *next = current->next;
            current->next = ll.llc;
            (ll.llc)->prev = current;
            (ll.llc)->next = next;

            if (next != NULL) {
                next->prev = ll.llc;
            }

            break;
        }

        current = current->next;
    }

    from->rate += weight;

    llp_t *next = to->parents->next;
    to->parents->next = ll.llp;
    (ll.llp)->next = next;
    (ll.llp)->prev = to->parents;

    if (next != NULL) {
        next->prev = ll.llp;
    }

    verify_next((ll_t*)from->edges);
    verify_next((ll_t*)to->parents);
}

vertex_t *vertex_init(vec_entry_t *state, vector<double> rewards, size_t state_length) {
    vertex_t *vertex = new vertex_t(state, rewards, state_length);
    llc_t *llc_dummy = (llc_t *) calloc(1, sizeof(llc_t));
    llp_t *llp_dummy = (llp_t *) calloc(1, sizeof(llp_t));

    vertex->edges = llc_dummy;
    vertex->parents = llp_dummy;

    return vertex;
}

void vertex_destroy(vertex_t *vertex) {
    // TODO: Delete edges?
    delete (vertex);
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

    llc_t *child = vertex->edges->next;

    while (child != NULL) {
        _reset_graph_visited(child->child, reset_int);
        child = child->next;
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
        llc_t *child = vertex->edges->next;

        while (child != NULL) {
            queue.push(child->child);
            child = child->next;
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

        if (vertex->edges->next == NULL) {
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

                    ll_insert(vertex, new_vertex, t);
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

    ll_insert(start, state_graph, 1);

    *graph = start;
    print_graph_list(stderr, start, true, n, n);

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
    fprintf(stream, ":\n");


    llc_t *child = vertex->edges->next;

    while (child != NULL) {
        fprintf(stream, "\t");
        fprintf(stream, "(%f) ", child->weight);
        print_vector_spacing(stream, child->child->state,
                             vec_length, vec_spacing);

        if (indexed) {
            fprintf(stream, " (%zu)", vertex->vertex_index);
        }

        //fprintf(stream, "\n");
        child = child->next;
    }


    fprintf(stream, "\n");

    llp_t *parent = vertex->parents->next;

    while (parent != NULL) {
        fprintf(stream, "\t");
        fprintf(stream, "P (%f) ", parent->llc->weight);
        print_vector_spacing(stream, parent->parent->state,
                             vec_length, vec_spacing);

        if (indexed) {
            fprintf(stream, " (%zu)", vertex->vertex_index);
        }

        //fprintf(stream, "\n");
        parent = parent->next;
    }

    fprintf(stream, "\n");

    child = vertex->edges->next;

    while (child != NULL) {
        _print_graph_list(stream, child->child,
                          indexed,
                          vec_length, vec_spacing);
        child = child->next;
    }
}

void print_graph_list(FILE *stream, vertex_t *graph,
                      bool indexed,
                      size_t vec_length, size_t vec_spacing) {
    return;
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

    llc_t *child = vertex->edges->next;

    while (child != NULL) {
        mph_cov_assign_desc_all(child->child, m);
        child = child->next;
    }

    child = vertex->edges->next;

    while (child != NULL) {
        for (size_t j = 0; j < m; ++j) {
            vertex->desc[j] += child->weight / vertex->rate * child->child->desc[j];
        }

        child = child->next;
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
    llc_t *child = vertex->edges->next;

    while (child != NULL) {
        _mph_cov_all(child->child, m);
        child = child->next;
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

    llc_t *child = vertex->edges->next;

    while (child != NULL) {
        rate += child->weight;
        child = child->next;
    }

    return rate;
}

// TODO: use parents->next and children->next everywhere

int reward_transform(vertex_t *graph, double (*reward_func)(vertex_t *)) {
    queue<vertex_t *> queue = enqueue_vertices(graph);
    print_graph_list(stderr, graph, false, 1, 1);
    while (!queue.empty()) {
        vertex_t *vertex = queue.front();
        queue.pop();

        //fprintf(stderr, "============ VERTEX  %zu  =============\n",
        //        *vertex->state);

        if (vertex->edges->next == NULL) {
            // Absorbing vertex
            continue;
        }

        if (vertex->parents->next == NULL) {
            // Starting vertex
            continue;
        }

        double reward = reward_func(vertex);


        if (reward == 0) {

            verify_next((ll_t *) vertex->edges);
            verify_next((ll_t *) vertex->parents);
            // Take all my edges and add to my parent instead.
            llp_t *parent = vertex->parents->next;

            while (parent != NULL) {


                verify_next((ll_t *) vertex->edges);
                verify_next((ll_t *) vertex->parents);
                verify_next((ll_t *) parent->parent->edges);
                verify_next((ll_t *) parent->parent->parents);
                llc_t *parent_child = parent->parent->edges;
                llc_t *child = vertex->edges->next;

                while (child != NULL) {
                    //fprintf(stderr, "Visiting parent %zu and child %zu\n",
                    //        *parent->parent->state,
                    //        *child->child->state);


                    verify_next((ll_t *) vertex->edges);
                    verify_next((ll_t *) vertex->parents);
                    verify_next((ll_t *) child->child->edges);
                    verify_next((ll_t *) child->child->parents);
                    verify_next((ll_t *) parent->parent->edges);
                    verify_next((ll_t *) parent->parent->parents);

                    double prob = child->weight / vertex->rate;

                    if (vertex->rate == 0) {
                        prob = 0;
                    }

                    if (child->child == parent->parent) {
                        continue;
                    }

                    if (parent_child->next == NULL ||
                        child->child < parent_child->next->child) {
                        //fprintf(stderr, "A\n");

                        verify_next((ll_t *) child->child->edges);
                        verify_next((ll_t *) child->child->parents);
                        // We must insert new
                        struct ll_init ll = ll_init(parent->parent, child->child, prob * parent->llc->weight);

                        llc_t *next_child = parent_child->next;
                        parent_child->next = ll.llc;
                        (ll.llc)->next = next_child;
                        (ll.llc)->prev = parent_child;

                        if (next_child != NULL) {
                            next_child->prev = ll.llc;
                        }


                        verify_next((ll_t *) child->child->edges);
                        verify_next((ll_t *) child->child->parents);

                        llp_t *next_parent = child->child->parents->next;
                        child->child->parents->next = ll.llp;
                        (ll.llp)->next = next_parent;
                        (ll.llp)->prev = child->child->parents;

                        if (next_parent != NULL) {
                            next_parent->prev = ll.llp;
                        }

                        verify_next((ll_t *) vertex->edges);
                        verify_next((ll_t *) vertex->parents);
                        verify_next((ll_t *) child->child->edges);
                        verify_next((ll_t *) child->child->parents);
                        verify_next((ll_t *) parent->parent->edges);
                        verify_next((ll_t *) parent->parent->parents);

                        child = child->next;
                    } else if (child->child > parent_child->next->child) {
                        verify_next((ll_t *) child->child->edges);
                        verify_next((ll_t *) child->child->parents);
                        //fprintf(stderr, "B\n");
                        parent_child = parent_child->next;
                    } else {
                        verify_next((ll_t *) child->child->edges);
                        verify_next((ll_t *) child->child->parents);
                        //fprintf(stderr, "C prob %f weight %f\n", prob, parent->llc->weight);
                        // ==
                        parent_child->next->weight += prob * parent->llc->weight;

                        verify_next((ll_t *) child->child->edges);
                        verify_next((ll_t *) child->child->parents);
                        child = child->next;
                        parent_child = parent_child->next;
                        if (child != NULL) {
                            verify_next((ll_t *) child->child->edges);
                            verify_next((ll_t *) child->child->parents);
                        }
                    }

                    verify_next((ll_t *) vertex->edges);
                    verify_next((ll_t *) vertex->parents);
                    if (child != NULL) {
                        verify_next((ll_t *) child->child->edges);
                        verify_next((ll_t *) child->child->parents);
                    }
                    verify_next((ll_t *) parent->parent->edges);
                    verify_next((ll_t *) parent->parent->parents);
                    print_graph_list(stderr, graph, false, 1, 1);
                }

                //fprintf(stderr, "Pre remove parents\n");
                print_graph_list(stderr, graph, false, 1, 1);
                ll_remove(parent->llc);
                parent = vertex->parents->next;

                //fprintf(stderr, "POST remove parents\n");
                print_graph_list(stderr, graph, false, 1, 1);
            }

            llc_t *child = vertex->edges->next;

            verify_next((ll_t *) vertex->edges);
            verify_next((ll_t *) vertex->parents);

            while (child != NULL) {
                if (*vertex->state == 3) {
                    //fprintf(stderr, "SPECIAL 3 PRE CHILD RM\n");
                    print_graph_list(stderr, graph, false, 1, 1);
                }
                ll_remove(child);
                child = vertex->edges->next;
            }

            //fprintf(stderr, "POST remove children\n");
            print_graph_list(stderr, graph, false, 1, 1);
            parent = vertex->parents->next;

            while (parent != NULL) {
                // ll_remove(parent->llc);

                parent = vertex->parents->next;
            }


            //fprintf(stderr, "POST remove parents2\n");
            print_graph_list(stderr, graph, false, 1, 1);

            vertex_destroy(vertex);
        } else {
            llc_t *child = vertex->edges->next;

            while (child != NULL) {
                child->weight /= reward;
                child = child->next;
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

        llc_t *child = vertex->edges->next;

        while (child != NULL) {
            queue.push(child->child);
            child = child->next;
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

    llc_t *child = vertex->edges->next;

    while (child != NULL) {
        vertex_t *child_vertex = child->child;

        weights[vertex->vertex_index][child_vertex->vertex_index] = child->weight;
        weights[vertex->vertex_index][vertex->vertex_index] -= child->weight;

        child = child->next;
    }

    vertices[vertex->vertex_index] = vertex;

    child = vertex->edges->next;

    while (child != NULL) {
        insert_into_weight_mat(weights, vertices, child->child);
        child = child->next;
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

void vertex_add_edge(vertex_t *from, vertex_t *to, double weight) {
    ll_insert(from, to, weight);
}

