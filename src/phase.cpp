#include <iostream>
#include <vector>
#include <cstring>
#include <stack>
#include <queue>
#include <stdint.h>
#include <math.h>
#include "phase.h"

static int reward_transform_vertex(vertex_t *vertex, double (*reward_func)(vertex_t *));

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

    free(llp);
}

/*
 * Assumes that my children will have removed me from their
 * parents list independently of this call
 */
void vertex_destroy(vertex_t *vertex) {
    ll_destroy(vertex->parents);
    free(vertex->edges);
    delete (vertex);
}

/*
 * Destroys all parent links from my children
 */
void vertex_destroy_parents(vertex_t *vertex) {
    for (size_t i = 0; i < vertex->nedges; ++i) {
        ll_remove(vertex->edges[i].llp);
    }
}


inline void ensure_valid_children(vertex_t *from) {
    return;
    if (from->nedges == 0) {
        return;
    }

    for (size_t i = 0; i < from->nedges - 1; ++i) {
        if (from->edges[i].child >= from->edges[i + 1].child) {
            DIE_ERROR(1, "Child %p  at index %zu is larger than (or eq) to next child %p at index %zu\n",
                      (void *) from->edges[i].child,
                      i,
                      (void *) from->edges[i + 1].child,
                      i + 1);
        }
    }
}

inline void add_edge_no_rate(vertex_t *from, vertex_t *to, size_t index, double weight) {
    if (to == NULL) {
        DIE_ERROR(1, "to is NULL\n");
    }

    llc_t *entry = &(from->edges[index]);
    DEBUG_PRINT("Adding edge from %p to %p at index %zu with weight %f\n",
                (void *) from, (void *) to, index, weight);
    /*for (size_t i = 0; i < index; ++i) {
        if (from->edges[i].child >= to) {
            DIE_ERROR(1, "Child %p to be inserted at index %zu is smaller than (or eq) to prev child %p at index %zu\n",
                      (void *) to, index, (void *) from->edges[i].child, i);
        }
    }*/

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
}

inline void add_edge(vertex_t *from, vertex_t *to, size_t index, double weight) {
    add_edge_no_rate(from, to, index, weight);
    from->rate += weight;
}

void vertex_add_edge(vertex_t *from, vertex_t *to, double weight) {
    //TODO super slow
    llc_t *edges = from->edges;
    from->edges = (llc_t *) calloc(from->nedges + 1, sizeof(llc_t));

    size_t inc = 0;
    size_t nedges = from->nedges + 1;

    for (size_t i = 0; i < nedges; ++i) {
        if (inc == 0 &&
            (i == nedges - 1 || edges[i].child > to)) {
            add_edge(from, to, i, weight);
            inc = 1;
        } else {
            from->edges[i] = edges[i - inc];
        }
    }

    for (size_t i = 0; i < from->nedges; ++i) {
        from->edges[i].llp->llc = &(from->edges[i]);
    }

    free(edges);
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

static void reset_graph_visited(vertex_t *graph) {
    queue<vertex_t *> to_reset;
    queue<vertex_t *> to_visit;

    to_visit.push(graph);

    while (!to_visit.empty()) {
        vertex_t *vertex = to_visit.front();
        to_visit.pop();

        if (!vertex->reset) {
            continue;
        }

        vertex->reset = false;
        to_reset.push(vertex);

        vertex->visited = false;

        for (size_t i = 0; i < vertex->nedges; ++i) {
            vertex_t *child = vertex->edges[i].child;
            to_visit.push(child);
        }
    }


    while (!to_reset.empty()) {
        vertex_t *vertex = to_reset.front();
        to_reset.pop();

        vertex->reset = true;
    }
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

    if (child_right_left != NULL) {
        child_right_left->parent = child;
    }

    child_right->left = child;
    child->parent = child_right;
    child_right_right = child_right->right;
    parent->left = child_right_right;

    if (child_right_right != NULL) {
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

    if (child_left_right != NULL) {
        child_left_right->parent = child;
    }

    child_left->right = child;

    child->parent = child_left;
    child_left_left = child_left->left;
    parent->right = child_left_left;

    if (child_left_left != NULL) {
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

    if (child_left != NULL) {
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

    if (child_right != NULL) {
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
    if ((*vertex = (avl_vec_vertex_t *) malloc(sizeof(avl_vec_vertex_t))) == NULL) {
        return 1;
    }

    (*vertex)->key = key;
    (*vertex)->entry = entry;
    (*vertex)->left = NULL;
    (*vertex)->right = NULL;
    (*vertex)->parent = parent;
    (*vertex)->balance = 0;

    return 0;
}

void avl_vec_vertex_destroy(avl_vec_vertex_t *vertex) {
    if (vertex == NULL) {
        return;
    }

    avl_vec_vertex_destroy(vertex->left);
    avl_vec_vertex_destroy(vertex->right);

    free(vertex);
}

static void avl_free(avl_vec_vertex_t *vertex) {
    if (vertex == NULL) {
        return;
    }

    avl_free(vertex->left);
    avl_free(vertex->right);
    free(vertex);
}

const avl_vec_vertex_t *
avl_vec_find(const avl_vec_vertex_t *rootptr, const vec_entry_t *key, const size_t vec_length) {
    if (rootptr == NULL) {
        return NULL;
    }

    const avl_vec_vertex_t *vertex = rootptr;

    while (true) {
        int res = radix_cmp(key, vertex->key, vec_length);
        if (res < 0) {
            if (vertex->left == NULL) {
                return NULL;
            } else {
                vertex = vertex->left;
            }
        } else if (res > 0) {
            if (vertex->right == NULL) {
                return NULL;
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
    if (avl_vec_vertex_create(out, key, entry, NULL)) {
        return -1;
    }

    if (rootptr == NULL) {
        return 1;
    }

    avl_vec_vertex_t *vertex = rootptr;

    while (true) {
        int res = radix_cmp(key, vertex->key, vec_length);
        if (res < 0) {
            if (vertex->left == NULL) {
                vertex->left = *out;
                break;
            } else {
                vertex = vertex->left;
            }
        } else if (res > 0) {
            if (vertex->right == NULL) {
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

    if (*root == NULL) {
        if (avl_vec_vertex_create(root, key, entry, NULL)) {
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

    for (avl_vec_vertex_t *parent = child->parent; parent != NULL; parent = child->parent) {
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

        if (pivot != NULL) {
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
    if (vertex == NULL) {
        return 0;
    }

    return 1 + avl_vec_get_size(vertex->left) + avl_vec_get_size(vertex->right);
}

queue<vertex_t *> enqueue_vertices(vertex_t *graph) {
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
        DEBUG_PRINT("Adding %p to queue\n", (void *) vertex);
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
    vertex_t *abs_vertex = NULL;

    while (!queue.empty()) {
        vertex_t *vertex = queue.front();
        queue.pop();

        if (vertex->nedges == 0) {
            if (abs_vertex == NULL) {
                abs_vertex = vertex;
            } else {
                DIE_ERROR(1, "Found multiple absorbing vertices (both %zu and %zu)\n",
                          abs_vertex->vertex_index, vertex->vertex_index);
            }
        }
    }

    if (abs_vertex == NULL) {
        DIE_ERROR(1, "No absorbing vertex found!\n");
    }

    return abs_vertex;
}

struct edge {
    vertex_t *vertex;
    double weight;
};

int edgecmp(const void *a, const void *b) {
    uintptr_t vertex_position_a = (uintptr_t) ((struct edge *) a)->vertex;
    uintptr_t vertex_position_b = (uintptr_t) ((struct edge *) b)->vertex;
    ptrdiff_t diff = vertex_position_a - vertex_position_b;

    if (diff < 0) {
        return -1;
    } else if (diff > 0) {
        return 1;
    } else {
        return 0;
    }
}

vertex_t *generate_state_space(
        size_t state_length,
        vector<pair<double, vector<size_t> > >(*visit_function)(vector<size_t>),
        vector<pair<double, vector<size_t> > >(*initial_states)(void),
        vector<double>(*rewards)(vector<size_t>)
) {
    DEBUG_PRINT("Start generate state space\n");

    avl_vec_vertex_t *bst = NULL;
    queue<pair<vertex_t *, vector<pair<double, vector<size_t> > > > > vertices_to_visit;
    queue<vertex_t *> q;
    vertex_t *start_vertex = vertex_init(
            NULL,
            rewards(vector<size_t>(state_length, 0)),
            0
    );

    start_vertex->vertex_index = 1;

    vertex_t *absorbing_vertex = vertex_init(
            NULL,
            rewards(vector<size_t>(state_length, 0)),
            0
    );

    absorbing_vertex->vertex_index = 0;

    size_t index = 2;

    vector<edge> initial_edges;
    vector<pair<double, vector<size_t> > > initial =
            initial_states();

    for (vector<pair<double, vector<size_t> > >::iterator it =
            initial.begin(); it != initial.end(); ++it) {
        double weight = it->first;
        vector<size_t> state(it->second);

        if (state.size() != state_length) {
            fprintf(stderr, "States must have given length '%zu', output state had length '%zu'",
                    state_length, state.size());
            goto error;
        }

        vertex_t *child;
        avl_vec_vertex_t *bst_entry = (avl_vec_vertex_t *) avl_vec_find(bst, &state[0], state_length);

        if (bst_entry != NULL) {
            fprintf(stderr, "Two edges with the same state have already been added to the initial state\n");
            goto error;
        }

        vector<pair<double, vector<size_t> > > child_children = visit_function(state);

        if (child_children.empty()) {
            child = absorbing_vertex;
        } else {
            size_t *vec_entry = (size_t *) calloc(
                    state.size(),
                    sizeof(size_t)
            );

            memcpy(vec_entry, &state[0], state.size() * sizeof(size_t));

            child = vertex_init(
                    vec_entry,
                    rewards(state),
                    state_length
            );

            child->vertex_index = index++;

            avl_vec_insert(&bst, vec_entry, child, state_length);

            pair<vertex_t *, vector<pair<double, vector<size_t> > > > pair(child, child_children);
            vertices_to_visit.push(pair);
        }

        struct edge edge;
        edge.weight = weight;
        edge.vertex = child;

        initial_edges.push_back(edge);
    }

    if (initial_edges.size() != 0) {
        start_vertex->edges = (llc_t *) calloc(
                initial_edges.size(), sizeof(llc_t)
        );

        struct edge *e = &initial_edges[0];
        qsort(e, initial_edges.size(), sizeof(struct edge), edgecmp);

        for (size_t i = 0; i < initial_edges.size(); ++i) {
            add_edge(start_vertex, e[i].vertex, i, e[i].weight);
        }

        for (size_t i = 0; i < initial_edges.size(); ++i) {
            start_vertex->edges[i].llp->llc = &(start_vertex->edges[i]);
        }
    } else {
        fprintf(stderr, "Start vertex must have one edge");
        goto error;
    }

    while (!vertices_to_visit.empty()) {
        pair<vertex_t *, vector<pair<double, vector<size_t> > > > p = vertices_to_visit.front();
        vertices_to_visit.pop();

        vector<edge> edges;
        vertex_t *visiting_vertex = p.first;
        vector<pair<double, vector<size_t> > > visiting_children = p.second;

        for (vector<pair<double, vector<size_t> > >::iterator it =
                visiting_children.begin(); it != visiting_children.end(); ++it) {
            double child_weight = it->first;
            vector<size_t> child_state = it->second;

            if (child_state.size() != state_length) {
                fprintf(stderr, "States must have given length '%zu', output state had length '%zu'",
                        state_length, child_state.size());
                goto error;
            }

            vertex_t *child;
            avl_vec_vertex_t *bst_entry = (avl_vec_vertex_t *) avl_vec_find(bst, &child_state[0], state_length);

            if (bst_entry == NULL) {
                vector<pair<double, vector<size_t> > > children = visit_function(child_state);

                if (children.empty()) {
                    child = absorbing_vertex;
                } else {
                    size_t *vec_entry = (size_t *) calloc(
                            child_state.size(),
                            sizeof(size_t)
                    );

                    memcpy(vec_entry, &child_state[0], child_state.size() * sizeof(size_t));

                    child = vertex_init(
                            vec_entry,
                            rewards(child_state),
                            state_length
                    );

                    avl_vec_insert(&bst, vec_entry, child, state_length);

                    pair<vertex_t *, vector<pair<double, vector<size_t> > > > pair(child, children);
                    vertices_to_visit.push(pair);
                }
            } else {
                child = bst_entry->entry;
            }

            struct edge edge;
            edge.weight = child_weight;
            edge.vertex = child;

            edges.push_back(edge);
        }

        if (edges.size() != 0) {
            // The list of edges may contain the same vertex twice
            size_t n_edges = 0;

            struct edge *e = &edges[0];
            qsort(e, edges.size(), sizeof(struct edge), edgecmp);
            vertex_t *prev_vertex = NULL;

            for (size_t i = 0; i < edges.size(); ++i) {
                vertex_t *child = edges[i].vertex;

                if (prev_vertex != child) {
                    n_edges++;
                }

                prev_vertex = child;
            }

            visiting_vertex->edges = (llc_t *) calloc(
                    n_edges, sizeof(llc_t)
            );

            prev_vertex = NULL;
            size_t k = 0;

            for (size_t i = 0; i < edges.size(); ++i) {
                vertex_t *child = e[i].vertex;

                if (prev_vertex != child) {
                    add_edge(visiting_vertex, child, k, e[i].weight);
                    k++;
                }

                prev_vertex = child;
            }

            for (size_t i = 0; i < n_edges; ++i) {
                // TODO: needed?
                visiting_vertex->edges[i].llp->llc = &(visiting_vertex->edges[i]);
            }
        }
    }

    DEBUG_PRINT("Freeing BST\n");

    avl_free(bst);

    q = enqueue_vertices(start_vertex);
    DEBUG_PRINT("Queue size: %zu\n", q.size());

    while (!q.empty()) {
        DEBUG_PRINT("Checking children of %p\n",
                    (void *) q.front());
        ensure_valid_children(q.front());
        q.pop();
    }

    return start_vertex;

    error:
    DEBUG_PRINT("FAILURE\n");
    avl_free(bst);
    return NULL;
}

static int kingman_visit_vertex(vertex_t **out_initial_vertex,
                                vec_entry_t *initial_state,
                                vertex_t *abs_vertex,
                                const size_t m) {
    avl_vec_vertex_t *bst = NULL;

    queue<vertex_t *> vertices_to_visit;
    vertex_t *initial_vertex = vertex_init(initial_state,
                                           vector<double>(initial_state, initial_state + m),
                                           m);
    vertices_to_visit.push(initial_vertex);
    vec_entry_t *v = (vec_entry_t *) malloc(sizeof(vec_entry_t) * m);
    avl_vec_vertex_t *bst_vertex;

    size_t end = 0;

    initial_vertex->vertex_index = 2;
    size_t index = 3;

    while (!vertices_to_visit.empty()) {
        vertex_t *vertex = vertices_to_visit.front();
        vertices_to_visit.pop();
        memcpy(v, vertex->state, sizeof(vec_entry_t) * m);
        size_t n_remaining = 0;
        size_t start = -1;

        vector<edge> edges;

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

                        if (bst_vertex == NULL) {
                            vec_entry_t *new_state = (vec_entry_t *) malloc(sizeof(vec_entry_t) * m);
                            memcpy(new_state, v, sizeof(vec_entry_t) * m);
                            vertex_t *to = vertex_init(
                                    new_state,
                                    vector<double>(new_state, new_state + m),
                                    m
                            );

                            to->vertex_index = index++;

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

                    struct edge edge;
                    edge.weight = t;
                    edge.vertex = new_vertex;

                    edges.push_back(edge);
                }
            }
        }

        if (edges.size() != 0) {
            vertex->edges = (llc_t *) calloc(edges.size(), sizeof(llc_t));
            struct edge *e = &edges[0];
            qsort(e, edges.size(), sizeof(struct edge), edgecmp);

            for (size_t i = 0; i < edges.size(); ++i) {
                add_edge(vertex, e[i].vertex, i, e[i].weight);
            }

            for (size_t i = 0; i < edges.size(); ++i) {
                vertex->edges[i].llp->llc = &(vertex->edges[i]);
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

    vertex_t *absorbing_vertex = vertex_init(
            mrca,
            vector<double>(mrca, mrca + m),
            m
    );

    absorbing_vertex->vertex_index = 0;

    vertex_t *state_graph;

    kingman_visit_vertex(
            &state_graph,
            initial, absorbing_vertex,
            m
    );

    vec_entry_t *start_state = (vec_entry_t *) calloc(m, sizeof(vec_entry_t));

    vertex_t *start = vertex_init(start_state, vector<double>(start_state, start_state + m), m);

    start->edges = (llc_t *) calloc(1, sizeof(llc_t));
    add_edge(start, state_graph, 0, 1);

    start->vertex_index = 1;

    // TODO: remove this. Replace entire function with generic one
    get_abs_vertex(start)->vertex_index = 0;

    *graph = start;

    return 0;
}

size_t reward_state = 0;

static double reward_by_state(vertex_t *vertex) {
    return vertex->rewards[reward_state];
}

static int kingman_visit_vertex_rw(vertex_t **out_initial_vertex,
                                   vec_entry_t *initial_state,
                                   vertex_t *abs_vertex,
                                   const size_t rw) {
    reward_state = rw;
    size_t m = rw + 2;
    avl_vec_vertex_t *bst = NULL;

    queue<vertex_t *> vertices_to_visit;

    vertex_t *initial_vertex = vertex_init(initial_state,
                                           vector<double>(initial_state, initial_state + m),
                                           m);
    vertices_to_visit.push(initial_vertex);
    vec_entry_t *v = (vec_entry_t *) malloc(sizeof(vec_entry_t) * m);
    avl_vec_vertex_t *bst_vertex;

    size_t end = 0;
    size_t previous_layer = 0;
    vector<vertex_t *> unrewarded_vertices;
    vector<vertex_t *> current_layer_vertices;

    while (!vertices_to_visit.empty()) {
        vertex_t *vertex = vertices_to_visit.front();
        vertices_to_visit.pop();


        memcpy(v, vertex->state, sizeof(vec_entry_t) * m);
        size_t n_remaining = 0;
        size_t start = -1;

        vector<edge> edges;

        for (vec_entry_t i = 0; i < m; i++) {
            n_remaining += v[i];

            if (v[i] > 0) {
                end = i;

                if (start == (size_t) -1) {
                    start = i;
                }
            }
        }

        if (n_remaining != previous_layer) {
            previous_layer = n_remaining;

            for (vector<vertex_t *>::iterator it =
                    unrewarded_vertices.begin(); it != unrewarded_vertices.end(); ++it) {
                reward_transform_vertex(*it, &reward_by_state);
            }

            unrewarded_vertices = current_layer_vertices;
            current_layer_vertices.clear();
        }

        current_layer_vertices.push_back(vertex);

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

                        if (bst_vertex == NULL) {
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

                    struct edge edge;
                    edge.weight = t;
                    edge.vertex = new_vertex;

                    edges.push_back(edge);
                }
            }
        }

        if (edges.size() != 0) {
            vertex->edges = (llc_t *) calloc(edges.size(), sizeof(llc_t));
            struct edge *e = &edges[0];
            qsort(e, edges.size(), sizeof(struct edge), edgecmp);

            for (size_t i = 0; i < edges.size(); ++i) {
                add_edge(vertex, e[i].vertex, i, e[i].weight);
            }

            for (size_t i = 0; i < edges.size(); ++i) {
                vertex->edges[i].llp->llc = &(vertex->edges[i]);
            }
        }
    }

    for (vector<vertex_t *>::iterator it =
            unrewarded_vertices.begin(); it != unrewarded_vertices.end(); ++it) {
        reward_transform_vertex(*it, &reward_by_state);
    }

    for (vector<vertex_t *>::iterator it =
            current_layer_vertices.begin(); it != current_layer_vertices.end(); ++it) {
        reward_transform_vertex(*it, &reward_by_state);
    }

    unrewarded_vertices.clear();
    current_layer_vertices.clear();

    free(v);
    *out_initial_vertex = initial_vertex;
    avl_free(bst);
    return 0;
}

int gen_kingman_graph_rw(vertex_t **graph, size_t n, size_t rw) {
    size_t m = rw + 2;

    vec_entry_t *initial = (vec_entry_t *) calloc(m, sizeof(vec_entry_t));
    initial[0] = n;

    vec_entry_t *mrca = (vec_entry_t *) calloc(m, sizeof(vec_entry_t));
    mrca[m - 1] = 1;

    vertex_t *absorbing_vertex = vertex_init(mrca,
                                             vector<double>(mrca, mrca + m),
                                             m);

    vertex_t *state_graph;

    kingman_visit_vertex_rw(&state_graph,
                            initial, absorbing_vertex,
                            rw);

    vec_entry_t *start_state = (vec_entry_t *) calloc(m, sizeof(vec_entry_t));

    vertex_t *start = vertex_init(start_state, vector<double>(start_state, start_state + m), m);

    start->edges = (llc_t *) calloc(1, sizeof(llc_t));
    add_edge(start, state_graph, 0, 1);

    reward_transform_vertex(state_graph, &reward_by_state);

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
    fprintf(stream, " %p", (void *) vertex);
    fprintf(stream, ":\n");

    for (size_t i = 0; i < vertex->nedges; ++i) {
        llc_t child_edge = vertex->edges[i];
        vertex_t *child = child_edge.child;
        fprintf(stream, "\t");

        fprintf(stream, " %p ", (void *) child);


        print_vector_spacing(stream, child->state,
                             vec_length, vec_spacing);

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

        fprintf(stream, "P %p ", (void *) parent->parent);

        print_vector_spacing(stream, parent->parent->state,
                             vec_length, vec_spacing);
        //fprintf(stream, "(%p) ", (void*)&parent);
        //fprintf(stream, " (llc%p) ", (void*)parent->llc);
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
    //reset_graph_visited(graph);
    fprintf(stream, "-- Graph list --\n");
    _print_graph_list(stream, graph, indexed, vec_length, vec_spacing);
    fflush(stream);
}

void assign_probability(vertex_t *graph) {
    queue<vertex_t *> q = enqueue_vertices(graph);

    while (!q.empty()) {
        vertex_t *vertex = q.front();
        q.pop();

        vertex->prob = 0;
    }

    // Start has probability 1
    graph->prob = 1;

    q = enqueue_vertices(graph);

    while (!q.empty()) {
        vertex_t *vertex = q.front();
        q.pop();

        for (size_t i = 0; i < vertex->nedges; ++i) {
            llc_t edge = vertex->edges[i];
            vertex_t *child = edge.child;

            child->prob += vertex->prob * edge.weight / vertex->rate;
        }
    }
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

    vertex->no_defect = vector<double>(m, 0);

    for (size_t i = 0; i < vertex->nedges; ++i) {
        llc_t child_edge = vertex->edges[i];

        for (size_t j = 0; j < m; ++j) {
            double probability = child_edge.weight / vertex->rate;
            vertex->desc[j] += probability * child_edge.child->desc[j];
            vertex->no_defect[j] += probability * child_edge.child->no_defect[j];
        }
    }

    for (size_t j = 0; j < m; ++j) {
        if (vertex->rewards[j] != 0) {
            vertex->no_defect[j] = 1;
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
double *defect;

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
    defect = (double *) calloc(m, sizeof(double));

    for (size_t i = 0; i < m; ++i) {
        expectation[i] = graph->desc[i];
        defect[i] = 1 - graph->no_defect[i];

        for (size_t j = 0; j <= i; ++j) {
            cov[i][j] -= graph->desc[i] *
                         graph->desc[j];
        }
    }

    return (cov_exp_return) {.cov = &cov, .exp = &expectation, .defect = &defect};
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
        llc_t child_edge = vertex->edges[i];
        rate += child_edge.weight;
    }

    return rate;
}

int reward_transform(vertex_t *graph, double (*reward_func)(vertex_t *)) {
    DEBUG_PRINT("Reward transforming\n");
    queue<vertex_t *> queue = enqueue_vertices(graph);

    if (graph->nparents != 0) {
        DIE_ERROR(1, "Expected start vertex to have 0 children\n");
    }

    while (!queue.empty()) {
        vertex_t *vertex = queue.front();
        queue.pop();

        reward_transform_vertex(vertex, reward_func);
    }

    label_vertex_index(NULL, graph);

    return 0;
}

int reward_transform_vertex(vertex_t *vertex, double (*reward_func)(vertex_t *)) {
    DEBUG_PRINT("Visiting vertex %p\n",
                (void *) vertex);

    ensure_valid_children(vertex);

    if (vertex->nedges == 0) {
        // Absorbing vertex
        DEBUG_PRINT("Vertex %p is absorbing\n",
                    (void *) vertex);
        return 0;
    }

    if (vertex->nparents == 0) {
        // Starting vertex
        DEBUG_PRINT("Vertex %p is start\n",
                    (void *) vertex);
        return 0;
    }

    double reward = reward_func(vertex);

    if (reward == 0) {
        DEBUG_PRINT("\t zero reward\n");
        ensure_valid_children(vertex);
        llc_t *children = vertex->edges;
        size_t nchildren = vertex->nedges;


        // Take all my edges and add to my parent instead.
        llp_t *parent_edge = vertex->parents->next;

        while (parent_edge != NULL) {
            vertex_t *parent = parent_edge->parent;
            llc_t *parent_old_children = parent->edges;
            size_t parent_nchildren = parent->nedges;
            ensure_valid_children(parent);

            llc_t *new_parent_children = (llc_t *) calloc(
                    parent_nchildren + nchildren,
                    sizeof(llc_t)
            );

            DEBUG_PRINT("Me (rw0) Vertex has %zu children\n",
                        nchildren);

            for (size_t n = 0; n < nchildren; ++n) {
                DEBUG_PRINT("\tChild %p with edge %.*f\n",
                            (void *) children[n].child, 3,
                            children[n].weight);
            }


            DEBUG_PRINT("Parent Vertex has %zu children\n",
                        parent_nchildren);

            for (size_t n = 0; n < parent_nchildren; ++n) {
                DEBUG_PRINT("\tChild %p with edge %.*f\n",
                            (void *) parent_old_children[n].child, 3,
                            parent_old_children[n].weight);
            }

            parent->edges = new_parent_children;

            double parent_weight = parent_edge->llc->weight;

            size_t i = 0, j = 0, k;

            llc_t child_i, child_j;

            for (k = 0; i < parent_nchildren || j < nchildren;) {
                //TODO add null after all children in array
                // loop over that instead
                if (i >= parent_nchildren) {
                    while (j < nchildren) {
                        // We have more children, but our parent does not have
                        // any more.
                        child_j = children[j];
                        double prob = children[j].weight / vertex->rate;

                        if (child_j.child == parent) {
                            // A 'self-loop' is removed
                            parent->rate -= parent_weight * prob;
                            j++;
                            continue;
                        }

                        DEBUG_PRINT("CASE i>=: (i %zu, parent_nchildren %zu, j %zu, nchildren %zu)\n",
                                    i, parent_nchildren, j, nchildren);
                        add_edge_no_rate(parent, child_j.child, k, prob * parent_weight);
                        j++;
                        k++;
                    }

                    break;
                }

                if (j >= nchildren) {
                    while (i < parent_nchildren) {
                        child_i = parent_old_children[i];

                        if (child_i.child == vertex) {
                            i++;
                            continue;
                        }

                        new_parent_children[k] = child_i;
                        new_parent_children[k].llp->llc = &(new_parent_children[k]);
                        i++;
                        k++;
                    }

                    break;
                }

                child_i = parent_old_children[i];
                child_j = children[j];

                if (child_j.child == parent) {
                    double prob = children[j].weight / vertex->rate;
                    parent->rate -= parent_weight * prob;
                    j++;
                    continue;
                }

                if (child_i.child == vertex) {
                    i++;
                    continue;
                }

                if (j >= nchildren ||
                    child_i.child < child_j.child) {
                    DEBUG_PRINT("Case A (i %zu, j %zu, k %zu, nchildren %zu, child_i %p, child_j %p)\n",
                                i, j, k, nchildren, (void *) child_i.child, (void *) child_j.child);
                    new_parent_children[k] = child_i;
                    new_parent_children[k].llp->llc = &(new_parent_children[k]);
                    i++;
                } else if (child_i.child > child_j.child) {
                    DEBUG_PRINT("Case B (i %zu, j %zu, k %zu, child_i %p, child_j %p)\n",
                                i, j, k, (void *) child_i.child, (void *) child_j.child);
                    double prob = child_j.weight / vertex->rate;
                    add_edge(parent, child_j.child, k, prob * parent_weight);
                    parent->rate -= prob * parent_weight;
                    j++;
                } else {
                    // ==
                    DEBUG_PRINT("Case C (i %zu, j %zu, k %zu, child_i %p, child_j %p)\n",
                                i, j, k, (void *) child_i.child, (void *) child_j.child);
                    double prob = child_j.weight / vertex->rate;
                    new_parent_children[k] = child_i;
                    new_parent_children[k].weight += prob * parent_weight;
                    new_parent_children[k].llp->llc = &(new_parent_children[k]);
                    i++;
                    j++;
                }

                k++;
            }


            parent->nedges = k;
            free(parent_old_children);
            parent_edge = parent_edge->next;
        }

        vertex_destroy_parents(vertex);
        vertex_destroy(vertex);

    } else {
        DEBUG_PRINT("\t non-zero reward\n");

        for (size_t i = 0; i < vertex->nedges; ++i) {
            vertex->edges[i].weight /= reward;
        }

        vertex->rate /= reward;
    }

    return 0;
}

/*
 * Also ensures that the absorbing vertex has index 0
 */
int label_vertex_index(size_t *largest_index, vertex_t *graph) {
    vertex_t *abs_vertex = get_abs_vertex(graph);
    queue<vertex_t *> q = enqueue_vertices(graph);
    size_t index = 1;

    while (!q.empty()) {
        vertex_t *vertex = q.front();
        q.pop();

        if (vertex != abs_vertex) {
            vertex->vertex_index = index++;
        } else {
            DEBUG_PRINT("%p is absorbing vertex\n", (void *) vertex);

            vertex->vertex_index = 0;
        }

        DEBUG_PRINT("Labelled %p as %zu\n", (void *) vertex, vertex->vertex_index);
    }

    if (largest_index != NULL) {
        *largest_index = index - 1;
    }

    return 0;
}

struct graph_info get_graph_info(vertex_t *graph) {
    size_t vertices, edges;

    queue<vertex_t *> queue = enqueue_vertices(graph);
    vertices = (size_t) queue.size();
    edges = 0;

    while (!queue.empty()) {
        vertex_t *vertex = queue.front();
        queue.pop();
        edges += vertex->nedges;
    }

    struct graph_info graph_info;

    graph_info.vertices = vertices;
    graph_info.edges = edges;

    return graph_info;
}

void set_graph_rewards(vertex_t *graph, vector<double> (*set_rewards_func)(vector<double>)) {
    queue<vertex_t *> queue = enqueue_vertices(graph);

    while (!queue.empty()) {
        vertex_t *vertex = queue.front();
        queue.pop();
        vertex->rewards = set_rewards_func(vertex->rewards);
    }
}


void assign_vertex_dist(vertex_t *vertex) {
    if (vertex->visited) {
        return;
    }

    vertex->visited = true;

    for (size_t i = 0; i < vertex->nedges; ++i) {
        vertex_t *child = vertex->edges[i].child;
        assign_vertex_dist(child);
    }

    size_t max = 0;

    for (size_t i = 0; i < vertex->nedges; ++i) {
        vertex_t *child = vertex->edges[i].child;

        if (child->integer > max) {
            max = child->integer;
        }
    }

    vertex->integer = max + 1;
}

int same_vertex_cmp(const void *a, const void *b) {
    vertex_t *va = *((vertex_t **) (a));
    vertex_t *vb = *((vertex_t **) (b));

    if (va->integer != vb->integer) {
        return (int) va->integer - (int) vb->integer;
    }

    if (va->nedges != vb->nedges) {
        return (int) va->nedges - (int) vb->nedges;
    }

    if (va->rate != vb->rate) {
        return (int) va->rate - (int) vb->rate;
    }

    return 0;
}


void reduce_graph(vertex_t *graph) {
    reset_graph_visited(graph);
    assign_vertex_dist(graph);

    struct graph_info info = get_graph_info(graph);
    size_t graph_size = info.vertices;

    queue<vertex_t *> queue2;
    queue<vertex_t *> queue = enqueue_vertices(graph);
    vertex_t **vertices = (vertex_t **) calloc(graph_size, sizeof(vertex_t *));

    size_t k = 0;

    while (!queue.empty()) {
        vertex_t *vertex = queue.front();
        queue.pop();

        vertices[k] = vertex;
        k++;
    }

    qsort(vertices, graph_size, sizeof(vertex_t *), same_vertex_cmp);
    bool *removed = (bool *) calloc(graph_size, sizeof(bool));

    for (size_t i = 0; i < graph_size; ++i) {
        for (size_t j = i + 1; j < graph_size; ++j) {
            DEBUG_PRINT("comparing %i and %i (%i %i)\n",
                        removed[i] ? -1 : (int) (vertices[i]->vertex_index),
                        removed[j] ? -1 : (int) (vertices[j]->vertex_index),
                        removed[i], removed[j]);
            if (removed[i] || removed[j]) {
                continue;
            }

            //TODO: Do not loop over all twice, replace with linear loop
            vertex_t *vertex_i = vertices[i];
            vertex_t *vertex_j = vertices[j];

            if (vertex_i->integer != vertex_j->integer ||
                vertex_i->nedges != vertex_j->nedges) {
                continue;
            }

            DEBUG_PRINT("Vertices %zu and %zu are the same group with depth %zu and nedges %zu\n",
                        vertex_i->vertex_index, vertex_j->vertex_index,
                        vertex_i->integer, vertex_i->nedges);

            bool equal = true;

            for (size_t l = 0; l < vertex_i->nedges; ++l) {
                if (fabs(vertex_i->edges[l].weight - vertex_j->edges[l].weight) > 0.001 ||
                    vertex_i->edges[l].child != vertex_j->edges[l].child) {
                    equal = false;
                    break;
                }
            }

            if (!equal) {
                DEBUG_PRINT("But NOT equal:\n");
            }
            DEBUG_PRINT("Vertex %zu has %zu children\n", vertex_i->vertex_index,
                        vertex_i->nedges);

            for (size_t n = 0; n < vertex_i->nedges; ++n) {
                DEBUG_PRINT("\tChild %zu with edge %.*f\n",
                            vertex_i->edges[n].child->vertex_index,
                            3, vertex_i->edges[n].weight);
            }
            DEBUG_PRINT("Vertex %zu has %zu children\n", vertex_j->vertex_index,
                        vertex_j->nedges);

            for (size_t n = 0; n < vertex_j->nedges; ++n) {
                DEBUG_PRINT("\tChild %zu with edge %.*f\n",
                            vertex_j->edges[n].child->vertex_index,
                            3, vertex_j->edges[n].weight);
            }

            if (equal) {
                DEBUG_PRINT("And they are equal\n");
                vertex_t **parents = (vertex_t **) calloc(
                        vertex_j->nparents,
                        sizeof(vertex_t *)
                );

                double *weights = (double *) calloc(
                        vertex_j->nparents,
                        sizeof(double)
                );

                llp_t *parent_edge = vertex_j->parents->next;
                size_t p = 0;

                while (parent_edge != NULL) {
                    vertex_t *parent = parent_edge->parent;
                    parents[p] = parent;
                    weights[p] = parent_edge->llc->weight;
                    parent_edge = parent_edge->next;
                    p++;
                }

                for (size_t l = 0; l < p; ++l) {
                    vertex_t *parent = parents[l];
                    double weight_j = weights[l];
                    llc_t *parent_children = parent->edges;
                    size_t parent_nchildren = parent->nedges;

                    DEBUG_PRINT("Redirecting the parent %zu with an edge to vertex j (%zu) to vertex i (%zu)\n",
                                parent->vertex_index,
                                vertex_j->vertex_index, vertex_i->vertex_index);

                    DEBUG_PRINT("Parent had old %zu children:\n", parent_nchildren);

                    for (size_t n = 0; n < parent_nchildren; ++n) {
                        DEBUG_PRINT("\tChild %zu with edge %.*f\n",
                                    parent_children[n].child->vertex_index,
                                    3, parent_children[n].weight);
                    }

                    // Remove all children's parent link
                    for (size_t m = 0; m < parent_nchildren; ++m) {
                        llc_t llc = parent_children[m];
                        llc.child->nparents--;

                        ll_remove(llc.llp);
                    }

                    llc_t *new_parent_children = (llc_t *) calloc(
                            parent_nchildren + 1,
                            sizeof(llc_t)
                    );

                    parent->edges = new_parent_children;
                    bool has_vertex_i = false;
                    size_t index = 0;
                    parent->rate = 0;
                    parent->nedges = 0;

                    for (size_t m = 0; m < parent_nchildren; ++m) {
                        llc_t edge = parent_children[m];

                        if (edge.child == vertex_j) {
                            // If vertex_j was the last child, we must add vertex_i
                            // Unless it has been seen/added already
                            if (m + 1 == parent_nchildren && !has_vertex_i) {
                                add_edge(parent, vertex_i, index, weight_j);
                                index++;
                                has_vertex_i = true;
                            }

                            continue;
                        }

                        if (!has_vertex_i && parent_children[m].child > vertex_i) {
                            add_edge(parent, vertex_i, index, weight_j);
                            index++;
                            has_vertex_i = true;
                            m--;
                            continue;
                        }

                        if (edge.child == vertex_i) {
                            add_edge(parent, edge.child, index, edge.weight + weight_j);
                            index++;
                            has_vertex_i = true;
                            continue;
                        }

                        add_edge(parent, edge.child, index, edge.weight);
                        index++;
                    }

                    DEBUG_PRINT("Parent NEW %zu children:\n",
                                parent->nedges);

                    for (size_t n = 0; n < parent->nedges; ++n) {
                        DEBUG_PRINT("\tChild %zu with edge %.*f\n",
                                    parent->edges[n].child->vertex_index,
                                    3, parent->edges[n].weight);
                    }

                    free(parent_children);
                }

                free(weights);
                free(parents);

                // Remove all my parents link
                for (size_t m = 0; m < vertex_j->nedges; ++m) {
                    llc_t llc = vertex_j->edges[m];
                    llc.child->nparents--;

                    ll_remove(llc.llp);
                }

                vertex_destroy(vertex_j);
                removed[j] = true;

                /*queue2 = enqueue_vertices(graph);

                while (!queue2.empty()) {
                    vertex_t *vertex = queue2.front();
                    queue2.pop();

                    if (vertex == vertex_j) {
                        //TODO: Remove this check
                        DIE_ERROR(1, "Should have removed vertex, it still exists\n");
                    }
                }*/
            }
        }
    }

    free(removed);
    free(vertices);
}

double fold(vertex_t *graph, double (*vertex_func)(vertex_t *)) {
    assign_probability(graph);

    double fold = 0;

    queue<vertex_t *> q = enqueue_vertices(graph);

    // Remove start vertex
    q.pop();

    while (!q.empty()) {
        vertex_t *vertex = q.front();
        q.pop();

        // Do not invoke absorbing vertex
        if (vertex->nedges == 0) {
            continue;
        }

        fprintf(stderr, "Vertex prob %f function res %f\n",
                vertex->prob, vertex_func(vertex));

        fold += vertex->prob * vertex_func(vertex);
    }

    return fold;
}

void _calculate_prob(vertex_t *vertex, double *probs) {
    if (vertex->visited) {
        return;
    }

    vertex->visited = true;

    if (vertex->nparents == 0) {
        // Start vertex
        probs[vertex->vertex_index] = 1;

        return;
    }

    llp_t *parent = vertex->parents->next;

    while (parent != NULL) {
        _calculate_prob(parent->parent, probs);

        probs[vertex->vertex_index] +=
                parent->llc->weight / parent->parent->rate *
                probs[parent->parent->vertex_index];

        parent = parent->next;
    }
}

void calculate_prob(vertex_t *graph, size_t *size, double **probs) {
    vertex_t *abs = get_abs_vertex(graph);
    graph_info info = get_graph_info(graph);

    (*probs) = (double *) calloc(info.vertices, sizeof(double));
    *size = info.vertices;

    reset_graph_visited(graph);
    _calculate_prob(abs, *probs);
}

void _calculate_var(vertex_t *vertex, double *vars) {
    if (vertex->visited) {
        return;
    }

    vertex->visited = true;

    if (vertex->nedges == 0) {
        // abs vertex
        vars[vertex->vertex_index] = 0;

        return;
    }

    for (size_t i = 0; i < vertex->nedges; ++i) {
        llc_t edge = vertex->edges[i];

        _calculate_var(edge.child, vars);
        double prob = edge.weight / vertex->rate;

        vars[vertex->vertex_index] += prob * prob * (vars[edge.child->vertex_index]);
    }

    if (vertex->nparents != 0 && vertex->nedges != 0) {
        vars[vertex->vertex_index] += 2 / (vertex->rate * vertex->rate);
    }
}

void calculate_var(vertex_t *graph, size_t *size, double **vars) {
    graph_info info = get_graph_info(graph);

    (*vars) = (double *) calloc(info.vertices, sizeof(double));
    *size = info.vertices;

    reset_graph_visited(graph);
    _calculate_var(graph, *vars);
}

struct vertex_pdf *vertex_pdfs;

size_t fac(size_t n) {
    if (n == 0) {
        return 1;
    }

    return n * fac(n - 1);
}

int sign(long double n) {
    if (n < 0) {
        return -1;
    } else if (n > 0) {
        return 1;
    } else {
        // NOTE:
        return 1;
    }
}

#define EPSILON 0.0001

int cmp_pdf_part(const void *a, const void *b) {
    struct pdf_values *aa = (struct pdf_values *) a;
    struct pdf_values *bb = (struct pdf_values *) b;

    if (aa->n != bb->n) {
        return (int) (aa->n - bb->n);
    }

    if (abs(aa->lambda - bb->lambda) < EPSILON) {
        return 0;
    } else if (aa->lambda > bb->lambda) {
        return 1;
    } else {
        return -1;
    }
}

// TODO: Use fractions not decimal...
void _pdf(vertex_t *vertex, double (*reward_func)(vertex_t *)) {
    if (vertex->visited) {
        return;
    }

    vertex->visited = true;

    for (size_t f = 0; f < vertex->nedges; ++f) {
        llc_t edge = vertex->edges[f];

        _pdf(edge.child, reward_func);
    }

    if (vertex->nedges == 0) {
        vertex_pdfs[vertex->vertex_index].defect_prob = 1;
        vertex_pdfs[vertex->vertex_index].parts = new vector<struct pdf_values>();

        DEBUG_PRINT("I am vertex %zu I am the absorbing\n", vertex->vertex_index);
        return;
    }

    long double reward;

    if (vertex->nparents == 0) {
        reward = 0;
    } else {
        reward = reward_func(vertex);
    }

    DEBUG_PRINT("I am vertex %zu I have rewards %f %f %f %f\n", vertex->vertex_index,
                vertex->rewards[0], vertex->rewards[1], vertex->rewards[2], vertex->rewards[3]);


    if (reward == 0) {
        DEBUG_PRINT("My reward is zero\n");
        vertex_pdfs[vertex->vertex_index].parts =
                new vector<struct pdf_values>();

        vector<struct pdf_values> *parts =
                vertex_pdfs[vertex->vertex_index].parts;

        long double defect_prob = 0;

        for (size_t f = 0; f < vertex->nedges; ++f) {
            llc_t edge = vertex->edges[f];

            vector<struct pdf_values> *partsz =
                    vertex_pdfs[edge.child->vertex_index].parts;

            long double prob = edge.weight / vertex->rate;

            for (size_t p = 0; p < partsz->size(); ++p) {

                DEBUG_PRINT("PUSHING RC\n");
                (*parts).push_back((struct pdf_values) {
                        .lambda = (*partsz)[p].lambda,
                        .k = log(prob) + (*partsz)[p].k,
                        .n = (*partsz)[p].n,
                        .c = (*partsz)[p].c
                });
            }

            defect_prob += prob * vertex_pdfs[edge.child->vertex_index].defect_prob;
        }


        vertex_pdfs[vertex->vertex_index].defect_prob = defect_prob;
    } else {
        DEBUG_PRINT("My reward is not zero\n");
        vertex_pdfs[vertex->vertex_index].parts =
                new vector<struct pdf_values>();

        vector<struct pdf_values> *parts =
                vertex_pdfs[vertex->vertex_index].parts;

        for (size_t f = 0; f < vertex->nedges; ++f) {
            llc_t edge = vertex->edges[f];
            DEBUG_PRINT("\n====\n");

            long double prob = edge.weight / vertex->rate;
            long double mu = vertex->rate / reward;

            vector<struct pdf_values> *partsz =
                    vertex_pdfs[edge.child->vertex_index].parts;

            DEBUG_PRINT("I have a child %zu with %zu parts and defect of %Lf\n", edge.child->vertex_index,
                        (*partsz).size(), vertex_pdfs[edge.child->vertex_index].defect_prob);


            long double defect_probz = vertex_pdfs[edge.child->vertex_index].defect_prob;

            /*for (size_t i = 0; i < (*partsz).size(); ++i) {
                DEBUG_PRINT("\n== part %zu =\n", i);
                int czi = (*partsz)[i].c;
                long double kzi = (*partsz)[i].k;
                size_t nzi = (*partsz)[i].n;
                long double lambdazi = (*partsz)[i].lambda;

                DEBUG_PRINT("Child part %zu has lambdazi %Lf kzi %Lf (%Lf) nzi %zu\n",
                            i, lambdazi, kzi, expl(kzi), nzi);

                if (abs(lambdazi - (-mu)) < EPSILON) {
                    DEBUG_PRINT("The rates are the same (mu=%Lf)\n", mu);
                    long double newk = log(mu) + kzi - log(nzi);
                    long double k = log(prob) + newk;

                    if (abs(prob * newk) > EPSILON) {
                        DEBUG_PRINT("PUSHING F c %i lambda  %Lf k %Lf (exp %Lf) n %zu\n", czi, -mu, k, expl(k),
                                    nzi + 1);

                        parts->push_back((struct pdf_values) {
                                .lambda = -mu,
                                .k = k,
                                .n = nzi + 1,
                                .c = czi
                        });
                    }
                } else {
                    DEBUG_PRINT("The rates are NOT the same (%Lf != %Lf)\n", lambdazi, -mu);
                    long double a;
                    int c;

                    DEBUG_PRINT("powl(-lambdazi - mu, nzi) > 0, powl(%Lf - %Lf, %zu) = %Lf\n",
                                -lambdazi, mu, nzi, powl(-lambdazi - mu, nzi));

                    if (powl(-lambdazi - mu, nzi) > 0) {
                        a = log(mu) + kzi + log(fac(nzi - 1)) - log(powl(-lambdazi - mu, nzi));
                        c = czi;
                    } else {
                        a = log(mu) + kzi + log(fac(nzi - 1)) - log(powl(lambdazi + mu, nzi));
                        c = -czi;
                    }

                    //DEBUG_PRINT("My a is constructed from   %Lf * (%Lf * %zu)/ (pow(-%Lf - %Lf, %zu)) = %Lf\n",
                    //           mu, kzi, fac(nzi - 1), lambdazi, mu, nzi, a);
                    // Add the first part
                    DEBUG_PRINT("abs(log(prob) + a) = abs(log(%Lf) + %Lf)=%Lf", prob, a, abs(log(prob) + a));
                    if (abs(log(prob) + a) > EPSILON) {
                        //DEBUG_PRINT("The prob is %Lf giving a*prob=%Lf\n", prob, prob * a);
                        DEBUG_PRINT("PUSHING A c %i lambda  %Lf k %Lf (exp %Lf) n %zu\n",
                                    c,
                                    -mu, log(prob) + a, exp(log(prob) + a), (size_t) 1);

                        parts->push_back((struct pdf_values) {
                                .lambda = -mu,
                                .k = log(prob) + a,
                                .n = 1,
                                .c = c
                        });
                    }

                    for (int j = 0; j <= nzi - 1; ++j) {
                        int snzi = (int) nzi;
                        long double b;
                        int c2;

                        DEBUG_PRINT(
                                "powl(-lambdazi - mu, j - nzi) > 0, powl(%Lf - %Lf, %i - %i) = powl(%Lf, %i)= %Lf\n",
                                -lambdazi, mu, j, snzi, -lambdazi - mu, j - snzi, powl(-lambdazi - mu, j - snzi));

                        if (powl(-lambdazi - mu, j - snzi) > 0) {
                            b = log(mu) + kzi + log(fac(nzi - 1)) + log(powl(-lambdazi - mu, j - snzi)) -
                                log(fac((size_t) j));
                            c2 = -czi;
                            DEBUG_PRINT(
                                    "COMPUTING B by log(mu) + kzi + log(fac(nzi - 1)) + log(-lambdazi - mu) * (j - snzi) - log(fac((size_t)j)) = log(%Lf) + %Lf + log(fac(%zu - 1)) - log(-%Lf - %Lf) * (%i - %i) - log(fac((size_t)%i)) = %Lf (%Lf)\n",
                                    mu, kzi, nzi, lambdazi, mu, j, snzi, j, b, expl(b)
                            );
                        } else {
                            b = log(mu) + kzi + log(fac(nzi - 1)) + log(powl(lambdazi + mu, (j - snzi))) -
                                log(fac((size_t) j));
                            c2 = czi;
                            DEBUG_PRINT(
                                    "COMPUTING B by log(mu) + kzi + log(fac(nzi - 1)) + log(lambdazi + mu) * (j - snzi) - log(fac((size_t)j)) = log(%Lf) + %Lf + log(fac(%zu - 1)) - log(%Lf + %Lf) * (%i - %i) - log(fac((size_t)%i)) = %Lf (%Lf)\n",
                                    mu, kzi, nzi, lambdazi, mu, j, snzi, j, b, expl(b)
                            );
                        }

                        long double newk = log(prob) + b;
                        long double newlambda = (lambdazi);
                        size_t newn = (size_t) j + 1;

                        DEBUG_PRINT("log(prob) + b = log(%Lf) + %Lf = %Lf\n",
                                    prob, b, log(prob) + b);


                        if (abs(newk) > EPSILON) {
                            DEBUG_PRINT("PUSHING B c %i lambda  %Lf k %Lf (exp %Lf) n %zu\n", c2, newlambda, newk,
                                        exp(newk), newn);
                            parts->push_back((struct pdf_values) {
                                    .lambda = newlambda,
                                    .k = newk,
                                    .n = newn,
                                    .c=c2
                            });
                        }
                    }
                }
            }*/

            for (size_t i = 0; i < (*partsz).size(); ++i) {
                long double kzi2 = ((*partsz)[i].k);
                long double kzi = (*partsz)[i].c * expl((*partsz)[i].k);
                size_t nzi = (*partsz)[i].n;
                long double lambdazi = (*partsz)[i].lambda;
                int czi = (*partsz)[i].c;

                DEBUG_PRINT("Child part %zu has lambdazi %Lf kzi %Lf nzi %zu\n", i, lambdazi, kzi, nzi);

                if (abs(lambdazi - (-mu)) < EPSILON) {
                    DEBUG_PRINT("The rates are the same (mu=%Lf)\n", mu);
                    /* OLD{
                        long double newk = mu * kzi * 1 / (nzi);
                        DEBUG_PRINT("long double newk = mu * kzi * 1 / (nzi)=%Lf * %Lf * 1 / (%zu)=%Lf\n",
                                    mu, kzi, nzi, newk);
                        DEBUG_PRINT("abs(prob * newk)=abs(%Lf * %Lf)= %Lf\n",
                                    prob, newk, abs(prob * newk));
                        DEBUG_PRINT("OKAY SO IF WE TAKE THE EXP:%i %Lf\n",
                                    -sign(mu * kzi * 1 / (nzi)), expl(abs(log(mu) + kzi2 - log(nzi))));

                        if (abs(prob * newk) > EPSILON) {
                            DEBUG_PRINT("PUSHING F lambda OLD  %Lf k %Lf n %zu\n", -mu, prob * newk, nzi + 1);

                            parts->push_back((struct pdf_values) {
                                    .lambda = -mu, .k = log(prob * abs(newk)), .n = nzi + 1, .c=sign(newk)
                            });
                        }
                    }*/
                    if (prob > EPSILON) {
                        long double newk = log(mu) + kzi2 - log((long double) nzi);
                        long double k = log(prob) + newk;
                        DEBUG_PRINT("long double newk = log(mu) + kzi2 - log(nzi)="
                                    "log(%Lf) + %Lf - log(%zu)=(%Lf) + %Lf - %Lf=%Lf\n",
                                    mu, kzi2, nzi, log(mu), kzi2, log((long double) nzi), newk);

                        DEBUG_PRINT("log(prob) + newk=log(%Lf)+ %Lf= %Lf\n",
                                    prob, newk, log(prob) + newk);

                        DEBUG_PRINT("PUSHING F c %i lambda  %Lf k %Lf (exp %Lf) n %zu\n",
                                    czi, -mu, k, expl(k),
                                    nzi + 1);

                        parts->push_back((struct pdf_values) {
                                .lambda = -mu,
                                .k = k,
                                .n = nzi + 1,
                                .c = czi
                        });
                    }
                } else {
                    DEBUG_PRINT("The rates are NOT the same\n");
                    long double NEWlambda, NEWk;
                    size_t NEWn;
                    long double aa;
                    int cc;
                    {
                        DEBUG_PRINT("The rates are NOT the same (%Lf != %Lf)\n", lambdazi, -mu);
                        long double a;
                        int c;

                        DEBUG_PRINT("powl(-lambdazi - mu, nzi) > 0, powl(%Lf - %Lf, %zu) = %Lf\n",
                                    -lambdazi, mu, nzi, powl(-lambdazi - mu, nzi));

                        if (powl(-lambdazi - mu, nzi) > 0) {
                            DEBUG_PRINT("PATH 1\n");
                            a = log(mu) + kzi2 + log(fac(nzi - 1)) - log(powl(-lambdazi - mu, nzi));
                            c = czi;
                        } else {
                            DEBUG_PRINT("PATH 2\n");
                            a = log(mu) + kzi2 + log(fac(nzi - 1)) - log(powl(lambdazi + mu, nzi));
                            c = -czi;
                        }

                        //DEBUG_PRINT("My a is constructed from   %Lf * (%Lf * %zu)/ (pow(-%Lf - %Lf, %zu)) = %Lf\n",
                        //           mu, kzi, fac(nzi - 1), lambdazi, mu, nzi, a);
                        // Add the first part
                        DEBUG_PRINT("abs(log(prob) + a) = abs(log(%Lf) + %Lf)=%Lf\n", prob, a, abs(log(prob) + a));
                        {
                            //DEBUG_PRINT("The prob is %Lf giving a*prob=%Lf\n", prob, prob * a);
                            DEBUG_PRINT("PUSHING A c %i lambda  %Lf k %Lf (exp %Lf) n %zu\n",
                                        c,
                                        -mu, log(prob) + a, exp(log(prob) + a), (size_t) 1);

                            NEWlambda = -mu;
                            NEWk = c*expl(log(prob) + a);
                            NEWn = 1;
                            parts->push_back((struct pdf_values) {
                                    .lambda = -mu,
                                    .k = log(prob) + a,
                                    .n = 1,
                                    .c = c
                            });
                            aa = a;
                            cc = c;
                        }
                    }

                    long double a = mu * (kzi * fac(nzi - 1)) / pow(-lambdazi - mu, nzi);

                    DEBUG_PRINT("My a is constructed from   %Lf * (%Lf * %zu)/ (pow(-%Lf - %Lf, %zu)) = %Lf\n",
                                mu, kzi, fac(nzi - 1), lambdazi, mu, nzi, a);
                    // Add the first part
                    if (abs(prob * a) > EPSILON) {
                        DEBUG_PRINT("The prob is %Lf giving a*prob=%Lf\n", prob, prob * a);
                        DEBUG_PRINT("PUSHING A OLD lambda  %Lf k %Lf n %zu\n", -mu, prob * a, (size_t) 1);

                        /*parts->push_back((struct pdf_values) {
                                .lambda = -mu, .k = log(abs(prob * a)), .n = 1, .c=sign(prob * a)
                        });*/

                        if (abs(-mu - NEWlambda) > EPSILON) {
                            DIE_ERROR(1, "??");
                        }

                        if (prob * a - NEWk > EPSILON) {
                            DIE_ERROR(1, "??");
                        }
                    }




                    for (int j = 0; j <= nzi - 1; ++j) {
                        {
                            int snzi = (int) nzi;
                            long double b;
                            int c2;

                            DEBUG_PRINT(
                                    "powl(-lambdazi - mu, j - nzi) > 0, powl(%Lf - %Lf, %i - %i) = powl(%Lf, %i)= %Lf\n",
                                    -lambdazi, mu, j, snzi, -lambdazi - mu, j - snzi, powl(-lambdazi - mu, j - snzi));

                            if (powl(-lambdazi - mu, j - snzi) > 0) {
                                b = log(mu) + kzi2 + log(fac(nzi - 1)) + log(powl(-lambdazi - mu, j - snzi)) -
                                    log(fac((size_t) j));
                                c2 = -czi;
                                DEBUG_PRINT(
                                        "COMPUTING B by log(mu) + kzi + log(fac(nzi - 1)) + log(-lambdazi - mu) * (j - snzi) - log(fac((size_t)j)) = log(%Lf) + %Lf + log(fac(%zu - 1)) - log(-%Lf - %Lf) * (%i - %i) - log(fac((size_t)%i)) = %Lf (%Lf)\n",
                                        mu, kzi2, nzi, lambdazi, mu, j, snzi, j, b, expl(b)
                                );
                            } else {
                                b = log(mu) + kzi2 + log(fac(nzi - 1)) + log(powl(lambdazi + mu, (j - snzi))) -
                                    log(fac((size_t) j));
                                c2 = czi;
                                DEBUG_PRINT(
                                        "COMPUTING B by log(mu) + kzi + log(fac(nzi - 1)) + log(lambdazi + mu) * (j - snzi) - log(fac((size_t)j)) = log(%Lf) + %Lf + log(fac(%zu - 1)) - log(%Lf + %Lf) * (%i - %i) - log(fac((size_t)%i)) = %Lf (%Lf)\n",
                                        mu, kzi2, nzi, lambdazi, mu, j, snzi, j, b, expl(b)
                                );
                            }

                            long double newk = log(prob) + b;
                            long double newlambda = (lambdazi);
                            size_t newn = (size_t) j + 1;

                            DEBUG_PRINT("log(prob) + b = log(%Lf) + %Lf = %Lf\n",
                                        prob, b, log(prob) + b);

                            NEWlambda = 9999;
                            NEWk = 999;

                            //if (abs(newk) > EPSILON) {
                            {
                                NEWlambda = newlambda;
                                NEWk = c2*expl(newk);

                                DEBUG_PRINT("PUSHING B c %i lambda  %Lf k %Lf (exp %Lf) n %zu\n", c2, newlambda, newk,
                                            exp(newk), newn);
                                parts->push_back((struct pdf_values) {
                                        .lambda = newlambda,
                                        .k = newk,
                                        .n = newn,
                                        .c=c2
                                });
                            }
                        }
                        // OLD
                        long double newk = prob * -a * (1 / fac((size_t)j)) * pow(-lambdazi - mu, j);
                        long double newlambda = (lambdazi);
                        size_t newn = (size_t)j + 1;

                        if (abs(newk) > EPSILON) {
                            DEBUG_PRINT("PUSHING B OLD lambda  %Lf k %Lf n %zu\n", newlambda, newk, newn);
                            /*parts->push_back(
                                    (struct pdf_values) {
                                        .lambda = newlambda, .k = log(abs(newk)), .n = newn, .c=sign(
                                            newk)});*/



                            if (abs(newlambda - NEWlambda) > EPSILON) {
                                DIE_ERROR(1, "??");
                            }

                            if (newk - NEWk > EPSILON) {
                                DIE_ERROR(1, "??");
                            }
                        }

                    }
                }

            }
            /* OLD
            // Defect addition
            if (abs(prob * defect_probz * mu) > EPSILON) {
                long double k = log(prob) + log(defect_probz) + log(mu);
                DEBUG_PRINT("PUSHING E c %i lambda  %Lf k %Lf (exp %Lf) n %zu\n", 1, -mu, k, expl(k), (size_t) 1);
                parts->push_back((struct pdf_values) {
                        .lambda = -mu,
                        .k = k,
                        .n = 1,
                        .c = 1
                });
            }*/
            // Defect addition
            if (abs(prob * defect_probz) > EPSILON) {
                DEBUG_PRINT("PUSHING E lambda  %Lf k %Lf n %zu\n", -mu, prob * defect_probz * mu, (size_t) 1);
                parts->push_back((struct pdf_values) {
                        .lambda = -mu, .k = log(prob) + log(defect_probz) + log(mu), .n = 1, .c =1
                });
            }
        }

        vertex_pdfs[vertex->vertex_index].defect_prob = 0;
    }


    // Combine the same lambda/n
    // TODO: Do this already before...
/*
    struct pdf_values *values = (struct pdf_values *) calloc(
            vertex_pdfs[vertex->vertex_index].parts->size(),
            sizeof(struct pdf_values)
    );

    for (size_t i = 0; i < vertex_pdfs[vertex->vertex_index].parts->size(); ++i) {
        values[i] = (*vertex_pdfs[vertex->vertex_index].parts)[i];
    }

    qsort(values, vertex_pdfs[vertex->vertex_index].parts->size(), sizeof(struct pdf_values), cmp_pdf_part);

    DEBUG_PRINT("ORDER:\n");
    for (size_t i = 0; i < vertex_pdfs[vertex->vertex_index].parts->size(); ++i) {
        DEBUG_PRINT("%Lf %Lf %zu\n", values[i].lambda, values[i].k, values[i].n);
    }

    size_t len = vertex_pdfs[vertex->vertex_index].parts->size();

    delete (vertex_pdfs[vertex->vertex_index].parts);

    vertex_pdfs[vertex->vertex_index].parts = new vector<struct pdf_values>();
    long double k = 0;
    size_t prev_n = 0;
    long double prev_lambda = 0;

    for (size_t i = 0; i < len; ++i) {
        DEBUG_PRINT("(%zu) %Lf %Lf %zu\n", i, values[i].lambda, values[i].k, values[i].n);
        if (prev_n == 0 || (prev_n == values[i].n && abs(prev_lambda - values[i].lambda) < EPSILON)) {
            DEBUG_PRINT("Increasing k by %Lf to %Lf\n", values[i].k, values[i].k + k);
            k += values[i].k;
        } else {
            DEBUG_PRINT("ADDING?  %Lf %Lf %zu\n", prev_lambda, k, prev_n);

            if (abs(k) > EPSILON) {
                vertex_pdfs[vertex->vertex_index].parts->push_back(
                        (struct pdf_values) {.lambda = prev_lambda, .k = k, .n = prev_n});
            }

            k = 0;
            k += values[i].k;
        }

        prev_lambda = values[i].lambda;
        prev_n = values[i].n;
    }

    DEBUG_PRINT("ADDING?2  %Lf %Lf %zu\n", prev_lambda, k, prev_n);

    if (prev_n != 0) {
        if (abs(k) > EPSILON) {
            vertex_pdfs[vertex->vertex_index].parts->push_back(
                    (struct pdf_values) {.lambda = prev_lambda, .k = k, .n = prev_n});
        }
    }

    DEBUG_PRINT("NEW ORDER:\n");
    for (size_t i = 0; i < vertex_pdfs[vertex->vertex_index].parts->size(); ++i) {
        DEBUG_PRINT("%Lf %Lf %zu\n",
                    (*vertex_pdfs[vertex->vertex_index].parts)[i].lambda,
                    (*vertex_pdfs[vertex->vertex_index].parts)[i].k,
                    (*vertex_pdfs[vertex->vertex_index].parts)[i].n);
    }*/


    long double tt = 0.5;
    long double prob = 0;

    for (size_t i = 0; i < vertex_pdfs[vertex->vertex_index].parts->size(); ++i) {
        DEBUG_PRINT("MY PARTS c %i k %Lf (%Lf) lambda %Lf n %zu\n",
                    (*vertex_pdfs[vertex->vertex_index].parts)[i].c,
                    (*vertex_pdfs[vertex->vertex_index].parts)[i].k,
                    expl((*vertex_pdfs[vertex->vertex_index].parts)[i].k),
                    ((*vertex_pdfs[vertex->vertex_index].parts)[i].lambda),
                    (*vertex_pdfs[vertex->vertex_index].parts)[i].n);
        prob += (*vertex_pdfs[vertex->vertex_index].parts)[i].c *
                expl((*vertex_pdfs[vertex->vertex_index].parts)[i].k) *
                pow(tt, (*vertex_pdfs[vertex->vertex_index].parts)[i].n - 1) *
                expl((*vertex_pdfs[vertex->vertex_index].parts)[i].lambda * tt);
    }

    DEBUG_PRINT("RESULT %Lf\t %Lf\n", tt, prob);
    fprintf(stderr, "I am vertex %zu I have rewards %f %f %f %f %f\n", vertex->vertex_index,
            vertex->rewards[0], vertex->rewards[1], vertex->rewards[2], vertex->rewards[3], vertex->rewards[4]);

    fprintf(stderr, "My prob %Lf\n", prob);
    DEBUG_PRINT("\n============================\n\n");
}

void pdf(vertex_t *graph, double (*reward_func)(vertex_t *), struct vertex_pdf *out_vertex_pdfs) {
    graph_info info = get_graph_info(graph);

    // TODO: Should always be indexed
    label_vertex_index(NULL, graph);
    reset_graph_visited(graph);

    vertex_pdfs = (struct vertex_pdf *) calloc(info.vertices, sizeof(struct vertex_pdf));

    _pdf(graph, reward_func);

    // TODO: Clone
    out_vertex_pdfs->defect_prob = vertex_pdfs[1].defect_prob;
    out_vertex_pdfs->parts = vertex_pdfs[1].parts;
}