#include <iostream>
#include <vector>
#include <cstring>
#include <stack>
#include <queue>
#include <stdint.h>
#include <stddef.h>
#include <math.h>
#include <set>

extern void *create_matrix(long double **mat, size_t length);

extern void *matrix_invert(void *matrix, size_t size);

extern double matrix_get(void *matrix, size_t i, size_t j);

#ifdef PTD_USE_GSL

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

void *create_matrix(long double **mat, size_t length) {
    gsl_matrix *full = gsl_matrix_alloc(length, length);

    for (size_t k = 0; k < length; ++k) {
        for (size_t j = 0; j < length; ++j) {
            gsl_matrix_set(full, k, j, (double) mat[k][j]);
        }
    }

    return full;
}

void *matrix_invert(void *matrix, size_t size) {
    gsl_matrix *gsl_mat = (gsl_matrix *) matrix;
    gsl_matrix *inverse = gsl_matrix_alloc(size, size);

    int sign;
    gsl_permutation *p = gsl_permutation_alloc(size);

    gsl_linalg_LU_decomp(gsl_mat, p, &sign);

    gsl_linalg_LU_invert(gsl_mat, p, inverse);

    gsl_permutation_free(p);

    return inverse;
}

double matrix_get(void *matrix, size_t i, size_t j) {
    gsl_matrix *gsl_mat = (gsl_matrix *) matrix;

    return gsl_matrix_get(gsl_mat, i, j);
}

#else // PTD_USE_GSL

#endif // PTD_USE_GSL

#include "ptdalgorithms.h"


extern char *vertex_name(struct ptd_vertex *vertex);

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
    for (size_t i = 0; i < from->nedges; ++i) {
        if (from->edges[i].child == to) {
            from->edges[i].weight += weight;
            from->rate += weight;
            return;
        }
    }

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
struct avl_node *rotate_left_right(struct avl_node *parent, struct avl_node *child) {
    struct avl_node *child_right_left, *child_right_right;
    struct avl_node *child_right = child->right;
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
struct avl_node *rotate_right_left(struct avl_node *parent, struct avl_node *child) {
    struct avl_node *child_left_right, *child_left_left;
    struct avl_node *child_left = child->left;

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
struct avl_node *rotate_left(struct avl_node *parent, struct avl_node *child) {
    struct avl_node *child_left;

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
struct avl_node *rotate_right(struct avl_node *parent, struct avl_node *child) {
    struct avl_node *child_right;

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

struct avl_node *avl_vec_vertex_create(char *key, void *entry, struct avl_node *parent) {
    struct avl_node *vertex;

    if ((vertex = (struct avl_node *) malloc(sizeof(*vertex))) == NULL) {
        return NULL;
    }

    vertex->key = key;
    vertex->entry = entry;
    vertex->left = NULL;
    vertex->right = NULL;
    vertex->parent = parent;
    vertex->balance = 0;

    return vertex;
}

void avl_vec_vertex_destroy(struct avl_node *vertex) {
    if (vertex == NULL) {
        return;
    }

    avl_vec_vertex_destroy(vertex->left);
    avl_vec_vertex_destroy(vertex->right);

    free(vertex);
}

static void avl_free(struct avl_node *vertex) {
    if (vertex == NULL) {
        return;
    }

    avl_free(vertex->left);
    avl_free(vertex->right);
    free(vertex);
}

const struct avl_node *
avl_vec_find(const struct avl_node *rootptr, const char *key, const size_t vec_length) {
    if (rootptr == NULL) {
        return NULL;
    }

    const struct avl_node *vertex = rootptr;

    while (true) {
        int res = memcmp(key, vertex->key, vec_length);

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

int find_or_insert_vec(struct avl_node **out, struct avl_node *rootptr, char *key, void *entry,
                       const size_t vec_length) {
    if ((*out = avl_vec_vertex_create(key, entry, NULL)) == NULL) {
        return -1;
    }

    if (rootptr == NULL) {
        return 1;
    }

    struct avl_node *vertex = rootptr;

    while (true) {
        int res = memcmp(key, vertex->key, vec_length);

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
            return 0;
        }
    }

    (*out)->parent = vertex;

    return 0;
}

int avl_rebalance_tree(struct avl_node **root, struct avl_node *child) {
    struct avl_node *pivot, *rotated_parent;

    for (struct avl_node *parent = child->parent; parent != NULL; parent = child->parent) {
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

int avl_vec_insert(struct avl_node **root, char *key, void *entry, const size_t vec_length) {
    struct avl_node *child;

    if (*root == NULL) {
        if ((*root = avl_vec_vertex_create(key, entry, NULL)) == NULL) {
            return -1;
        }

        return 0;
    }

    int res = find_or_insert_vec(&child, *root, key, entry, vec_length);

    if (res == -1) {
        return -1;
    }

    if (res == 0) {
        return 0;
    }

    avl_rebalance_tree(root, child);

    return 0;
}


static size_t avl_vec_get_size(struct avl_node *vertex) {
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

vertex_t *get_abs_vertex(vertex_t *graph) {
    queue<vertex_t *> queue = enqueue_vertices(graph);
    vertex_t *abs_vertex = NULL;

    while (!queue.empty()) {
        vertex_t *vertex = queue.front();
        queue.pop();

        if (vertex->nedges == 0) {
            if (abs_vertex == NULL) {
                abs_vertex = vertex;
            } else {
                DIE_ERROR(1, "Found multiple absorbing internal_vertices (both %zu and %zu)\n",
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

    struct avl_node *bst = NULL;
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
        struct avl_node *bst_entry = (struct avl_node *) avl_vec_find(bst, (char *) &state[0],
                                                                      state_length * sizeof(vec_entry_t));

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

            avl_vec_insert(&bst, (char *) vec_entry, child, state_length * sizeof(vec_entry_t));

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
            struct avl_node *bst_entry = (struct avl_node *) avl_vec_find(bst, (char *) &child_state[0],
                                                                          state_length * sizeof(vec_entry_t));

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

                    avl_vec_insert(&bst, (char *) vec_entry, child, state_length * sizeof(vec_entry_t));

                    pair<vertex_t *, vector<pair<double, vector<size_t> > > > pair(child, children);
                    vertices_to_visit.push(pair);
                }
            } else {
                child = (vertex_t *) bst_entry->entry;
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


size_t reward_state = 0;

static double reward_by_state(vertex_t *vertex) {
    return vertex->rewards[reward_state];
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

            DEBUG_PRINT("Vertices %zu and %zu are the same group with depth %i and nedges %zu\n",
                        vertex_i->vertex_index, vertex_j->vertex_index,
                        vertex_i->integer, vertex_i->nedges);

            bool equal = true;

            for (size_t l = 0; l < vertex_i->nedges; ++l) {
                if (fabsl(vertex_i->edges[l].weight - vertex_j->edges[l].weight) > 0.001 ||
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

    if (fabsl(aa->lambda - bb->lambda) < EPSILON) {
        return 0;
    } else if (aa->lambda > bb->lambda) {
        return 1;
    } else {
        return -1;
    }
}

vector<struct pdf_values> *combine_densities(vector<struct pdf_values> *parts) {
    struct pdf_values *new_values = (struct pdf_values *) calloc(
            parts->size(),
            sizeof(struct pdf_values)
    );

    for (size_t i = 0; i < parts->size(); ++i) {
        new_values[i] = (*parts)[i];
    }

    qsort(new_values, parts->size(), sizeof(struct pdf_values), cmp_pdf_part);
    size_t len = parts->size();

    vector<struct pdf_values> *new_parts = new vector<struct pdf_values>();

    if (len > 0) {
        size_t prev_n = new_values[0].n;
        long double prev_lambda = new_values[0].lambda;
        size_t g = 0;
        vector<vector<size_t> > groups;
        groups.push_back(vector<size_t>());

        for (size_t i = 0; i < len; ++i) {
            if (prev_n == new_values[i].n && fabsl(prev_lambda - new_values[i].lambda) < EPSILON) {
                groups[g].push_back(i);
            } else {
                groups.push_back(vector<size_t>());
                g++;
                groups[g].push_back(i);
                prev_n = new_values[i].n;
                prev_lambda = new_values[i].lambda;
            }
        }

        for (size_t g = 0; g < groups.size(); ++g) {
            long double kp = 0;
            long double kn = 0;

            for (size_t j = 0; j < groups[g].size(); ++j) {
                size_t i = groups[g][j];

                if (new_values[i].c == -1) {
                    kn += exp2l(new_values[i].k);
                } else {
                    kp += exp2l(new_values[i].k);
                }
            }

            long double k = log2l(fabsl(kp - kn));

            new_parts->push_back(
                    (struct pdf_values) {
                            .lambda = new_values[groups[g][0]].lambda,
                            .k = k,
                            .n = new_values[groups[g][0]].n,
                            .c = sign(kp - kn)
                    });
        }
    }

    return new_parts;
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

    vector<struct pdf_values> allchildparts;
    vector<struct pdf_values> *allchildparts2;

    for (size_t f = 0; f < vertex->nedges; ++f) {
        llc_t edge = vertex->edges[f];

        vector<struct pdf_values> *partsz =
                vertex_pdfs[edge.child->vertex_index].parts;

        long double prob = edge.weight / vertex->rate;

        for (size_t p = 0; p < partsz->size(); ++p) {
            allchildparts.push_back((struct pdf_values) {
                    .lambda = (*partsz)[p].lambda,
                    .k = log2l(prob) + (*partsz)[p].k,
                    .n = (*partsz)[p].n,
                    .c = (*partsz)[p].c
            });
        }
    }

    allchildparts2 = combine_densities(&allchildparts);
    allchildparts = *allchildparts2;

    vertex_pdfs[vertex->vertex_index].parts =
            new vector<struct pdf_values>();
    vector<struct pdf_values> *vertex_parts =
            vertex_pdfs[vertex->vertex_index].parts;

    if (reward == 0) {
        DEBUG_PRINT("My reward is zero\n");

        for (size_t p = 0; p < allchildparts.size(); ++p) {
            vertex_parts->push_back((struct pdf_values) {
                    .lambda = allchildparts[p].lambda,
                    .k = allchildparts[p].k,
                    .n = allchildparts[p].n,
                    .c = allchildparts[p].c
            });
        }

        long double defect_prob = 0;

        for (size_t f = 0; f < vertex->nedges; ++f) {
            llc_t edge = vertex->edges[f];
            long double prob = edge.weight / vertex->rate;
            defect_prob += prob * vertex_pdfs[edge.child->vertex_index].defect_prob;
        }

        vertex_pdfs[vertex->vertex_index].defect_prob = defect_prob;
    } else {
        DEBUG_PRINT("My reward is not zero\n");

        for (size_t f = 0; f < vertex->nedges; ++f) {
            llc_t edge = vertex->edges[f];
            long double prob = edge.weight / vertex->rate;
            long double mu = vertex->rate / reward;
            long double defect_probz = vertex_pdfs[edge.child->vertex_index].defect_prob;

            if (fabsl(prob * defect_probz) > EPSILON) {
                vertex_parts->push_back((struct pdf_values) {
                        .lambda = -mu,
                        .k = log2l(prob) + log2l(defect_probz) + log2l(mu),
                        .n = 1,
                        .c =1
                });
            }
        }

        vertex_pdfs[vertex->vertex_index].defect_prob = 0;

        for (size_t i = 0; i < allchildparts.size(); ++i) {
            long double mu = vertex->rate / reward;
            long double kzi = allchildparts[i].k;
            size_t nzi = allchildparts[i].n;
            long double lambdazi = allchildparts[i].lambda;
            int czi = allchildparts[i].c;

            if (fabsl(lambdazi - (-mu)) < EPSILON) {
                DEBUG_PRINT("The rates are the same (mu=%Lf)\n", mu);
                long double newk = log2l(mu) + kzi - log2l((long double) nzi);

                vertex_parts->push_back((struct pdf_values) {
                        .lambda = -mu,
                        .k = newk,
                        .n = nzi + 1,
                        .c = czi
                });
            } else {
                DEBUG_PRINT("The rates are NOT the same (%Lf != %Lf)\n", lambdazi, -mu);
                long double a;
                int c;

                if (powl(-lambdazi - mu, nzi) > 0) {
                    int s = sign(powl(-lambdazi - mu, nzi));
                    a = log2l(mu) + kzi + log2l(fac(nzi - 1)) - log2l(fabsl(-lambdazi - mu)) * nzi * s;
                    c = czi;
                } else {
                    int s = sign(powl(lambdazi + mu, nzi));
                    a = log2l(mu) + kzi + log2l(fac(nzi - 1)) - log2l(fabsl(lambdazi + mu)) * nzi * s;
                    c = -czi;
                }
                vertex_parts->push_back((struct pdf_values) {
                        .lambda = -mu,
                        .k = a,
                        .n = 1,
                        .c = c
                });

                for (int j = 0; j < nzi; ++j) {
                    int snzi = (int) nzi;
                    long double b;
                    int c2;
                    if (powl(-lambdazi - mu, j - snzi) > 0) {
                        int s = sign(powl(-lambdazi - mu, j - snzi));
                        b = log2l(mu) + kzi + log2l(fac(nzi - 1)) +
                            s * log2l(fabsl(-lambdazi - mu)) * (j - snzi) -
                            log2l(fac((size_t) j));
                        c2 = -czi;
                    } else {
                        int s = sign(powl(lambdazi + mu, j - snzi));
                        b = log2l(mu) + kzi + log2l(fac(nzi - 1)) + s * log2l(lambdazi + mu) * (j - snzi) -
                            log2l(fac((size_t) j));
                        c2 = czi;
                    }

                    vertex_parts->push_back((struct pdf_values) {
                            .lambda = lambdazi,
                            .k = b,
                            .n = (size_t) j + 1,
                            .c = c2
                    });
                }
            }
        }
    }

    vector<struct pdf_values> *new_parts = combine_densities(vertex_pdfs[vertex->vertex_index].parts);
    delete (vertex_pdfs[vertex->vertex_index].parts);
    vertex_pdfs[vertex->vertex_index].parts = new_parts;
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

queue<vertex_t *> topological_queue(vertex_t *graph) {
    queue<vertex_t *> topological_queue;
    queue<vertex_t *> q;

    q = enqueue_vertices(graph);

    while (!q.empty()) {
        vertex_t *vertex = q.front();
        q.pop();

        vertex->integer = (int) vertex->nparents;
    }

    reset_graph_visited(graph);
    q.push(graph);

    while (!q.empty()) {
        vertex_t *vertex = q.front();
        q.pop();

        topological_queue.push(vertex);

        for (size_t i = 0; i < vertex->nedges; ++i) {
            vertex_t *child = vertex->edges[i].child;

            child->integer -= 1;

            if (child->integer == 0 && !child->visited) {
                child->visited = true;
                q.push(child);
            }
        }
    }

    return topological_queue;
}

static void ptd_reset_graph_visited(struct ptd_graph *graph) {
    queue<struct ptd_vertex *> to_reset;
    queue<struct ptd_vertex *> to_visit;

    to_visit.push(graph->start_vertex);

    while (!to_visit.empty()) {
        struct ptd_vertex *vertex = to_visit.front();
        to_visit.pop();

        if (!vertex->reset) {
            continue;
        }

        vertex->reset = false;
        to_reset.push(vertex);

        vertex->visited = false;

        for (size_t i = 0; i < vertex->edges_length; ++i) {
            struct ptd_vertex *child = vertex->edges[i].to;
            to_visit.push(child);
        }
    }

    while (!to_reset.empty()) {
        struct ptd_vertex *vertex = to_reset.front();
        to_reset.pop();

        vertex->reset = true;
    }
}

queue<struct ptd_vertex *> ptd_enqueue_vertices(struct ptd_graph *graph) {
    ptd_reset_graph_visited(graph);
    // TODO: add failures to these
    queue<struct ptd_vertex *> ret;
    queue<struct ptd_vertex *> queue;

    queue.push(graph->start_vertex);

    while (!queue.empty()) {
        struct ptd_vertex *vertex = queue.front();
        queue.pop();

        if (vertex->visited) {
            continue;
        }

        vertex->visited = true;
        ret.push(vertex);

        for (size_t i = 0; i < vertex->edges_length; ++i) {
            struct ptd_edge edge = vertex->edges[i];
            queue.push(edge.to);
        }
    }

    return ret;
}

int cmp_vertex_index(const void *a, const void *b) {
    if ((*((const struct ptd_vertex **) a))->index > (*((const struct ptd_vertex **) b))->index) {
        return 1;
    } else if ((*((const struct ptd_vertex **) a))->index < (*((const struct ptd_vertex **) b))->index) {
        return -1;
    } else {
        return 0;
    }
}

queue<struct ptd_vertex *> ptd_enqueue_vertices_sorted(struct ptd_graph *graph) {
    if (!graph->is_indexed) {
        ptd_label_vertices(graph);
    }

    queue<struct ptd_vertex *> unsorted = ptd_enqueue_vertices(graph);
    queue<struct ptd_vertex *> sorted_queue;
    size_t len = unsorted.size();
    struct ptd_vertex **sorted = (struct ptd_vertex **) calloc(len, sizeof(*sorted));
    size_t index = 0;

    while (!unsorted.empty()) {
        struct ptd_vertex *vertex = unsorted.front();
        unsorted.pop();

        sorted[index] = vertex;
        index++;
    }

    qsort(sorted, index, sizeof(*sorted), cmp_vertex_index);

    for (size_t i = 0; i < len; ++i) {
        sorted_queue.push(sorted[i]);
    }

    free(sorted);

    return sorted_queue;
}

static struct ptd_vertex_linked_list_item *add_to_ll(
        struct ptd_graph *graph,
        struct ptd_vertex *vertex) {
    struct ptd_vertex_linked_list_item *item = (struct ptd_vertex_linked_list_item *) malloc(
            sizeof(*item)
    );

    if (item == NULL) {
        return NULL;
    }

    item->next = NULL;
    item->previous = graph->vertices_list->tail;
    item->vertex = vertex;

    if (graph->vertices_list->first == NULL) {
        graph->vertices_list->first = item;
        graph->vertices_list->tail = item;
        graph->vertices_list->vertices_count = 1;
    } else {
        graph->vertices_list->tail->next = item;
        graph->vertices_list->tail = item;
        graph->vertices_list->vertices_count++;
    }

    return item;
}

static void remove_from_ll(
        struct ptd_graph *graph,
        struct ptd_vertex *vertex
) {
    if (vertex->item_in_graph_list->previous != NULL) {
        vertex->item_in_graph_list->previous->next = vertex->item_in_graph_list->next;
    }

    if (vertex->item_in_graph_list->next != NULL) {
        vertex->item_in_graph_list->next->previous = vertex->item_in_graph_list->previous;
    }

    if (graph->vertices_list->first == vertex->item_in_graph_list) {
        graph->vertices_list->first = vertex->item_in_graph_list->next;
    }

    if (graph->vertices_list->tail == vertex->item_in_graph_list) {
        graph->vertices_list->tail = vertex->item_in_graph_list->previous;
    }

    free(vertex->item_in_graph_list);
    graph->vertices_list->vertices_count--;
}

struct ptd_graph *ptd_graph_create(size_t state_length) {
    struct ptd_graph *graph;
    struct ptd_vertex *start;

    if ((graph = (struct ptd_graph *) malloc(sizeof(struct ptd_graph))) == NULL) {
        return NULL;
    }

    graph->state_length = state_length;

    if ((graph->vertices_list = (struct ptd_vertex_linked_list *) malloc(
            sizeof(*graph->vertices_list)
    )) == NULL) {
        // TODO: proper cleanup, e.g. call graph free allow for nulls
        free(graph);
        return NULL;
    }

    graph->vertices_list->first = NULL;
    graph->vertices_list->tail = NULL;

    if ((start = ptd_vertex_create(graph)) == NULL) {
        free(graph);
        return NULL;
    }

    graph->start_vertex = start;

    graph->is_indexed = true;

    return graph;
}

void ptd_graph_vertices_destroy(struct ptd_graph *graph) {
    for (struct ptd_vertex_linked_list_item *item = graph->vertices_list->first;
         item->next != NULL;) {
        ptd_vertex_destroy(item->next->vertex);
    }
}

void ptd_graph_destroy(struct ptd_graph *graph) {
    ptd_vertex_destroy(graph->start_vertex);

    struct ptd_vertex_linked_list_item *item = graph->vertices_list->first;
    struct ptd_vertex_linked_list_item *next;

    while (item != NULL) {
        next = item->next;
        free(item);
        item = next;
    }

    free(graph->vertices_list);

    memset(graph, 0, sizeof(*graph));
    free(graph);
}

queue<struct ptd_vertex *> *vertices_to_visit;

struct ptd_vertex *ptd_vertex_create_state(struct ptd_graph *graph, vec_entry_t *state) {
    struct ptd_vertex *vertex;

    vertex = (struct ptd_vertex *) malloc(sizeof(*vertex));

    if (vertex == NULL) {
        return NULL;
    }

    vertex->edges_limit = 8;
    vertex->edges = (struct ptd_edge *) calloc(8, sizeof(*vertex->edges));

    if (vertex->edges == NULL) {
        return NULL;
    }

    vertex->edges_length = 0;

    vertex->state = state;

    vertex->rate = 0;
    vertex->visited = false;
    vertex->reset = true;
    vertex->graph = graph;
    vertex->data = NULL;
    vertex->index = 0;

    vertex->item_in_graph_list = add_to_ll(graph, vertex);

    if (vertices_to_visit != NULL) {
        vertices_to_visit->push(vertex);
    }

    return vertex;
}

struct ptd_vertex *ptd_vertex_create(struct ptd_graph *graph) {
    vec_entry_t *state = (vec_entry_t *) calloc(graph->state_length, sizeof(*state));

    if (state == NULL) {
        return NULL;
    }

    return ptd_vertex_create_state(graph, state);
}

void ptd_vertex_destroy(struct ptd_vertex *vertex) {
    free(vertex->edges);
    vertex->edges = NULL;
    free(vertex->state);
    vertex->state = NULL;
    vertex->edges_length = 0;
    vertex->edges_limit = 0;
    free(vertex->data);
    vertex->data = NULL;
    remove_from_ll(vertex->graph, vertex);

    free(vertex);
}

ptd_avl_tree_t *ptd_avl_tree_create(size_t vec_length) {
    ptd_avl_tree_t *avl_tree = (ptd_avl_tree_t *) malloc(sizeof(ptd_avl_tree_t));

    if (avl_tree == NULL) {
        return NULL;
    }

    avl_tree->root = NULL;
    avl_tree->vec_length = vec_length;

    return avl_tree;
}

void _ptd_avl_tree_vertex_destroy(struct avl_node *avl_vertex) {
    if (avl_vertex == NULL) {
        return;
    }

    _ptd_avl_tree_vertex_destroy(avl_vertex->left);
    _ptd_avl_tree_vertex_destroy(avl_vertex->right);

    avl_vertex->left = NULL;
    avl_vertex->right = NULL;
    avl_vertex->entry = NULL;
    free(avl_vertex);
}

void ptd_avl_tree_vertex_destroy(ptd_avl_tree_t *avl_tree) {
    _ptd_avl_tree_vertex_destroy((struct avl_node *) avl_tree->root);
    avl_tree->root = NULL;
    free(avl_tree);
}

void _ptd_avl_tree_vertex_destroy_free(struct avl_node *avl_vertex) {
    if (avl_vertex == NULL) {
        return;
    }

    _ptd_avl_tree_vertex_destroy_free(avl_vertex->left);
    _ptd_avl_tree_vertex_destroy_free(avl_vertex->right);

    ptd_vertex_destroy((struct ptd_vertex *) avl_vertex->entry);
    avl_vertex->left = NULL;
    avl_vertex->right = NULL;
    avl_vertex->entry = NULL;
    free(avl_vertex);
}

void ptd_avl_tree_vertex_destroy_free(ptd_avl_tree_t *avl_tree) {
    _ptd_avl_tree_vertex_destroy_free((struct avl_node *) avl_tree->root);
    avl_tree->root = NULL;
    free(avl_tree);
}


void _ptd_avl_tree_edge_destroy(struct avl_node *avl_vertex) {
    if (avl_vertex == NULL) {
        return;
    }

    _ptd_avl_tree_edge_destroy(avl_vertex->left);
    _ptd_avl_tree_edge_destroy(avl_vertex->right);

    avl_vertex->left = NULL;
    avl_vertex->right = NULL;
    free(avl_vertex->entry);
    avl_vertex->entry = NULL;
    free(avl_vertex);
}

void ptd_avl_tree_edge_destroy(ptd_avl_tree_t *avl_tree) {
    _ptd_avl_tree_edge_destroy((struct avl_node *) avl_tree->root);
    avl_tree->root = NULL;
    free(avl_tree);
}

size_t ptd_avl_tree_max_depth(void *avl_vec_vertex) {
    if ((struct avl_node *) avl_vec_vertex == NULL) {
        return 0;
    }

    return max(
            ptd_avl_tree_max_depth((void *) ((struct avl_node *) avl_vec_vertex)->left) + 1,
            ptd_avl_tree_max_depth((void *) ((struct avl_node *) avl_vec_vertex)->left) + 1
    );
}


int ptd_avl_tree_vertex_insert(ptd_avl_tree_t *avl_tree, const vec_entry_t *key, const struct ptd_vertex *vertex) {
    int res;

    struct avl_node *root = (struct avl_node *) avl_tree->root;
    res = avl_vec_insert(&root, (char *) key, (void *) vertex, avl_tree->vec_length * sizeof(vec_entry_t));

    if (res != 0) {
        return res;
    }

    avl_tree->root = root;

    return 0;
}

struct ptd_vertex *ptd_avl_tree_vertex_find(const ptd_avl_tree_t *avl_tree, const vec_entry_t *key) {
    const struct avl_node *avl_vertex = avl_vec_find(
            (struct avl_node *) avl_tree->root,
            (char *) key,
            avl_tree->vec_length * sizeof(vec_entry_t)
    );

    /*fprintf(stderr, "finding vertex ");
    fprintf(stderr, "%zu %zu %zu %zu",
            key[0],
            key[1],
            key[2],
            key[3]
    );*/


    if (avl_vertex == NULL) {
        //fprintf(stderr, " NOT EXISTS\n");
        return NULL;
    }

    //fprintf(stderr, " EXISTS\n");

    return (struct ptd_vertex *) avl_vertex->entry;
}

int ptd_avl_tree_edge_insert_or_increment(ptd_avl_tree_t *avl_tree, const vec_entry_t *key, struct ptd_vertex *vertex,
                                          long double weight) {
    struct avl_node *root = (struct avl_node *) avl_tree->root;
    struct avl_node *child;

    if (root == NULL) {
        struct ptd_edge *edge = (struct ptd_edge *) malloc(sizeof(*edge));
        edge->weight = weight;
        edge->to = vertex;

        if ((root = avl_vec_vertex_create((char *) key, (void *) edge, NULL)) == NULL) {
            return -1;
        }

        avl_tree->root = root;

        return 0;
    }

    struct avl_node *parent = root;

    while (true) {
        int res = memcmp(parent->key, key, avl_tree->vec_length * sizeof(vec_entry_t));

        if (res < 0) {
            if (parent->left == NULL) {
                struct ptd_edge *edge = (struct ptd_edge *) malloc(sizeof(*edge));
                edge->weight = weight;
                edge->to = vertex;

                child = avl_vec_vertex_create((char *) key, edge, parent);

                if (child == NULL) {
                    return -1;
                }

                parent->left = child;

                break;
            } else {
                parent = parent->left;
            }
        } else if (res > 0) {
            if (parent->right == NULL) {
                struct ptd_edge *edge = (struct ptd_edge *) malloc(sizeof(*edge));
                edge->weight = weight;
                edge->to = vertex;

                child = avl_vec_vertex_create((char *) key, edge, parent);

                if (child == NULL) {
                    return -1;
                }

                parent->right = child;

                break;
            } else {
                parent = parent->right;
            }
        } else {
            ((struct ptd_edge *) parent->entry)->weight += weight;

            return 0;
        }
    }

    avl_rebalance_tree(&root, child);
    avl_tree->root = root;

    return 0;


    return 0;
}

/*struct avl_node * _ptd_avl_tree_edge_remove(struct avl_node *root, char *key, size_t length) {
    if (root == NULL) {
        return root;
    }

    int res = memcmp(root->key, key, length);

    if (res < 0) {
        root->left = _ptd_avl_tree_edge_remove(root->left, key, length);
    } else if (res > 0) {
        root->right = _ptd_avl_tree_edge_remove(root->right, key, length);
    } else {
        // node with only one child or no child
        if (root->left == NULL || root->right == NULL){
            struct avl_node *temp = NULL;

            if (temp == root->left)
                temp = root->right;
            else
                temp = root->left;

            // No child case
            if (temp == NULL) {
                temp = root;
                root = NULL;
            }
            else // One child case
                root = temp; // Copy the contents of
            // the non-empty child
        }
        else
        {
            Node temp = minValueNode(root.right);

            // Copy the inorder successor's data to this node
            root.key = temp.key;

            // Delete the inorder successor
            root.right = deleteNode(root.right, temp.key);
        }
    }

    // If the tree had only one node then return
    if (root == null)
        return root;

    // STEP 2: UPDATE HEIGHT OF THE CURRENT NODE
    root.height = max(height(root.left),
                      height(root.right)) + 1;

    // STEP 3: GET THE BALANCE FACTOR
    // OF THIS NODE (to check whether
    // this node became unbalanced)
    int balance = getBalance(root);

    // If this node becomes unbalanced,
    // then there are 4 cases
    // Left Left Case
    if (balance > 1 && getBalance(root.left) >= 0)
        return rightRotate(root);

    // Left Right Case
    if (balance > 1 && getBalance(root.left) < 0)
    {
        root.left = leftRotate(root.left);
        return rightRotate(root);
    }

    // Right Right Case
    if (balance < -1 && getBalance(root.right) <= 0)
        return leftRotate(root);

    // Right Left Case
    if (balance < -1 && getBalance(root.right) > 0)
    {
        root.right = rightRotate(root.right);
        return leftRotate(root);
    }

    return root;
}

int ptd_avl_tree_edge_remove(ptd_avl_tree_t *avl_tree, const vec_entry_t *key) {

}
*/
int ptd_avl_tree_edge_remove(ptd_avl_tree_t *avl_tree, const vec_entry_t *key) {
    struct avl_node *root = (struct avl_node *) avl_tree->root;

    if (root == NULL) {
        return 0;
    }

    struct avl_node *parent = root;

    enum DIR {
        NONE, LEFT, RIGHT
    };
    DIR dir = NONE;

    while (true) {
        int res = memcmp(parent->key, key, avl_tree->vec_length * sizeof(vec_entry_t));

        if (res < 0) {
            if (parent->left == NULL) {
                return 0;
            } else {
                parent = parent->left;
                dir = LEFT;
            }
        } else if (res > 0) {
            if (parent->right == NULL) {
                return 0;
            } else {
                parent = parent->right;
                dir = RIGHT;
            }
        } else {
            break;
        }
    }

    if (parent->parent != NULL) {
        if (dir == LEFT) {
            parent->parent->left = NULL;
        }

        avl_rebalance_tree((struct avl_node **) &avl_tree->root, parent->parent);
    } else {
        if (parent->left != NULL) {
            avl_tree->root = parent->left;
        } else if (parent->right != NULL) {
            avl_tree->root = parent->right;
        } else {
            avl_tree->root = NULL;
            return 0;
        }

        avl_rebalance_tree((struct avl_node **) &avl_tree->root, (struct avl_node *) avl_tree->root);
    }

    return 0;
}

struct ptd_edge *ptd_avl_tree_edge_find(const ptd_avl_tree_t *avl_tree, const vec_entry_t *key) {
    struct avl_node *parent = (struct avl_node *) avl_tree->root;

    while (true) {
        if (parent == NULL) {
            return NULL;
        }

        int res = memcmp(parent->key, (char *) key, avl_tree->vec_length * sizeof(vec_entry_t));

        if (res < 0) {
            parent = parent->left;
        } else if (res > 0) {
            parent = parent->right;
        } else {
            return (struct ptd_edge *) parent->entry;
        }
    }
}

int ptd_visit_vertices(struct ptd_graph *graph, int (*visit_func)(struct ptd_vertex *), bool include_start) {
    vertices_to_visit = new queue<struct ptd_vertex *>;
    *vertices_to_visit = ptd_enqueue_vertices_sorted(graph);

    while (!vertices_to_visit->empty()) {
        struct ptd_vertex *vertex = vertices_to_visit->front();
        vertices_to_visit->pop();

        if (!include_start && vertex == graph->start_vertex) {
            continue;
        }

        char buf[1024] = {0};
        ptd_vertex_to_s(vertex, buf, 1024);
        //fprintf(stderr, "Now visiting %s with index %zu\n", buf, vertex->index);
        int res = visit_func(vertex);

        if (res != 0) {
            return res;
        }
    }

    delete vertices_to_visit;
    vertices_to_visit = NULL;

    return 0;
}

int ptd_add_edge(struct ptd_vertex *from, struct ptd_vertex *to, long double weight) {
    /*fprintf(stderr, "Adding edge between ");
    fprintf(stderr, "%zu %zu %zu %zu ",
            from->state[0],
            from->state[1],
            from->state[2],
            from->state[3]
    );

    fprintf(stderr, "AND %zu %zu %zu %zu\n",
            to->state[0],
            to->state[1],
            to->state[2],
            to->state[3]
    );*/

    // TODO: improve the speed of this...
    for (size_t i = 0; i < from->edges_length; ++i) {
        if (from->edges[i].to == to) {
            from->edges[i].weight += weight;
            from->rate += weight;
            return 0;
        }
    }


    if (from->edges_length + 1 >= from->edges_limit) {
        from->edges_limit *= 2;
        from->edges = (struct ptd_edge *) realloc(
                from->edges,
                from->edges_limit * sizeof(struct ptd_edge)
        );

        if (from->edges == NULL) {
            return -1;
        }
    }

    from->edges[from->edges_length].to = to;
    from->edges[from->edges_length].weight = weight;
    from->edges_length++;
    from->rate += weight;
    from->graph->is_indexed = false;

    return 0;
}


stack<struct ptd_ph_vertex *> *scc_stack2;
vector<struct ptd_ph_scc_vertex *> *scc_components2;
size_t scc_index2;
size_t *scc_indices2;
size_t *low_links2;
bool *scc_on_stack2;
static bool *visited;
ptd_strongly_connected_components_t *sccs2;

int strongconnect2(struct ptd_ph_vertex *vertex) {
    scc_indices2[vertex->index] = scc_index2;
    low_links2[vertex->index] = scc_index2;
    visited[vertex->index] = true;
    scc_index2++;
    scc_stack2->push(vertex);
    scc_on_stack2[vertex->index] = true;

    for (size_t i = 0; i < vertex->edges_length; ++i) {
        struct ptd_ph_edge *edge = vertex->edges[i];

        if (!visited[edge->to->index]) {
            int res = strongconnect2(edge->to);

            if (res != 0) {
                return res;
            }

            low_links2[vertex->index] = min(
                    low_links2[vertex->index],
                    low_links2[edge->to->index]
            );
        } else if (scc_on_stack2[edge->to->index]) {
            low_links2[vertex->index] = min(
                    low_links2[vertex->index],
                    scc_indices2[edge->to->index]
            );
        }
    }

    if (low_links2[vertex->index] == scc_indices2[vertex->index]) {
        struct ptd_ph_vertex *w;
        vector<struct ptd_ph_vertex *> list;

        do {
            w = scc_stack2->top();
            scc_stack2->pop();
            scc_on_stack2[w->index] = false;

            list.push_back(w);
        } while (w != vertex);

        ptd_ph_scc_vertex *scc = (ptd_ph_scc_vertex *) malloc(sizeof(*scc));

        if (scc == NULL) {
            return -1;
        }

        scc->internal_edges_length = list.size();
        scc->internal_edges = (struct ptd_ph_edge **) calloc(
                scc->internal_edges_length,
                sizeof(*(scc->internal_edges))
        );

        for (size_t i = 0; i < scc->internal_edges_length; ++i) {
            scc->internal_edges[i] = (struct ptd_ph_edge *) malloc(sizeof(*(scc->internal_edges[i])));
            scc->internal_edges[i]->to = list.at(i);
        }

        scc_components2->push_back(scc);
    }

    return 0;
}

struct ptd_ph_scc_graph *ptd_ph_find_strongly_connected_components(struct ptd_ph_graph *graph) {
    struct ptd_ph_scc_graph *scc_graph = (struct ptd_ph_scc_graph *) malloc(
            sizeof(*scc_graph)
    );

    scc_stack2 = new stack<struct ptd_ph_vertex *>;

    scc_index2 = 0;
    scc_indices2 = (size_t *) calloc(graph->vertices_length, sizeof(size_t));
    low_links2 = (size_t *) calloc(graph->vertices_length, sizeof(size_t));
    scc_on_stack2 = (bool *) calloc(graph->vertices_length, sizeof(bool));
    visited = (bool *) calloc(graph->vertices_length, sizeof(bool));
    scc_components2 = new vector<ptd_ph_scc_vertex *>();

    size_t vertex_order;

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        struct ptd_ph_vertex *vertex = graph->vertices[i];

        if (!visited[i]) {
            if (strongconnect2(vertex) != 0) {
                return NULL;
            }
        }
    }

    size_t non_empty_components = 0;

    for (size_t i = 0; i < scc_components2->size(); ++i) {
        struct ptd_ph_scc_vertex *c = scc_components2->at(i);

        if (c->internal_edges_length != 0) {
            non_empty_components++;
        }
    }

    scc_graph->vertices_length = non_empty_components;
    scc_graph->vertices = (struct ptd_ph_scc_vertex **) calloc(
            scc_graph->vertices_length,
            sizeof(*(scc_graph->vertices))
    );

    size_t index = 0;

    for (size_t i = 0; i < scc_graph->vertices_length; ++i) {
        struct ptd_ph_scc_vertex *scc = scc_components2->at(i);

        if (scc->internal_edges_length != 0) {
            scc_graph->vertices[index] = scc_components2->at(i);
            index++;
        } else {
            free(scc->internal_edges);
            free(scc);
        }
    }

    struct ptd_ph_scc_vertex **sccs_for_vertices = (struct ptd_ph_scc_vertex **) calloc(
            graph->vertices_length,
            sizeof(*sccs_for_vertices)
    );

    for (size_t i = 0; i < scc_graph->vertices_length; ++i) {
        struct ptd_ph_scc_vertex *scc = scc_graph->vertices[i];
        scc->index = i;
        set<struct ptd_ph_scc_vertex *> external_sccs;
        set<struct ptd_ph_vertex *> external_vertices;

        for (size_t j = 0; j < scc->internal_edges_length; ++j) {
            struct ptd_ph_edge *edge = scc->internal_edges[j];
            struct ptd_ph_vertex *vertex = edge->to;

            sccs_for_vertices[vertex->index] = scc;
        }
    }

    for (size_t i = 0; i < scc_graph->vertices_length; ++i) {
        struct ptd_ph_scc_vertex *scc = scc_graph->vertices[i];
        set<struct ptd_ph_scc_vertex *> external_sccs;
        set<struct ptd_ph_vertex *> external_vertices;

        for (size_t j = 0; j < scc->internal_edges_length; ++j) {
            struct ptd_ph_edge *edge = scc->internal_edges[j];
            struct ptd_ph_vertex *vertex = edge->to;

            for (size_t k = 0; k < vertex->edges_length; ++k) {
                struct ptd_ph_vertex *child = vertex->edges[k]->to;
                struct ptd_ph_scc_vertex *child_scc = sccs_for_vertices[child->index];

                if (child_scc != scc) {
                    external_sccs.insert(child_scc);
                    external_vertices.insert(child);
                }
            }
        }

        scc->edges_length = external_sccs.size();
        scc->edges = (struct ptd_ph_scc_edge **) calloc(
                scc->edges_length,
                sizeof(*(scc->edges))
        );

        scc->external_edges_length = external_vertices.size();
        scc->external_edges = (struct ptd_ph_edge **) calloc(
                scc->external_edges_length,
                sizeof(*scc->external_edges)
        );

        size_t set_index;

        set_index = 0;

        for (set<struct ptd_ph_scc_vertex *>::iterator itr = external_sccs.begin();
             itr != external_sccs.end();
             itr++) {
            scc->edges[set_index] = (struct ptd_ph_scc_edge*) malloc(sizeof(*(scc->edges[set_index])));
            scc->edges[set_index]->to = *itr;
            set_index++;
        }

        set_index = 0;

        for (set<struct ptd_ph_vertex *>::iterator itr = external_vertices.begin();
             itr != external_vertices.end();
             itr++) {
            scc->external_edges[set_index] = (struct ptd_ph_edge*) malloc(sizeof(*(scc->external_edges[set_index])));
            scc->external_edges[set_index]->to = *itr;
            set_index++;
        }
    }

    free(scc_indices2);
    free(low_links2);
    free(scc_on_stack2);
    free(visited);
    delete scc_components2;

    free(sccs_for_vertices);

    return scc_graph;
}

void ptd_ph_scc_graph_destroy(struct ptd_ph_scc_graph *scc_graph) {
    for (size_t i = 0; i < scc_graph->vertices_length; ++i) {
        struct ptd_ph_scc_vertex *scc = scc_graph->vertices[i];

        for (size_t j = 0; j < scc->edges_length; ++j) {
            free(scc->edges[j]);
        }

        free(scc->edges);

        for (size_t j = 0; j < scc->internal_edges_length; ++j) {
            free(scc->internal_edges[j]);
        }

        free(scc->internal_edges);

        for (size_t j = 0; j < scc->external_edges_length; ++j) {
            free(scc->external_edges[j]);
        }

        free(scc->external_edges);

        free(scc);
    }

    free(scc_graph->vertices);
    free(scc_graph);
}

stack<struct ptd_vertex *> *scc_stack;
vector<ptd_strongly_connected_component_t *> *scc_components;
size_t scc_index;
size_t *scc_indices;
size_t *low_links;
bool *scc_on_stack;
ptd_strongly_connected_components_t *sccs;

int strongconnect(struct ptd_vertex *vertex) {
    scc_indices[vertex->index] = scc_index;
    low_links[vertex->index] = scc_index;
    vertex->visited = true;
    scc_index++;
    scc_stack->push(vertex);
    scc_on_stack[vertex->index] = true;

    for (size_t i = 0; i < vertex->edges_length; ++i) {
        struct ptd_edge edge = vertex->edges[i];

        if (!edge.to->visited) {
            int res = strongconnect(edge.to);

            if (res != 0) {
                return res;
            }

            low_links[vertex->index] = min(
                    low_links[vertex->index],
                    low_links[edge.to->index]
            );
        } else if (scc_on_stack[edge.to->index]) {
            low_links[vertex->index] = min(
                    low_links[vertex->index],
                    scc_indices[edge.to->index]
            );
        }
    }

    if (low_links[vertex->index] == scc_indices[vertex->index]) {
        struct ptd_vertex *w;
        vector<struct ptd_vertex *> list;

        do {
            w = scc_stack->top();
            scc_stack->pop();
            scc_on_stack[w->index] = false;

            list.push_back(w);
        } while (w != vertex);

        ptd_strongly_connected_component_t *scc = (ptd_strongly_connected_component_t *) malloc(sizeof(*scc));

        if (scc == NULL) {
            return -1;
        }

        scc->internal_vertices_length = list.size();
        scc->internal_vertices = (struct ptd_vertex **) calloc(scc->internal_vertices_length,
                                                               sizeof(struct ptd_vertex *));

        for (size_t i = 0; i < scc->internal_vertices_length; ++i) {
            scc->internal_vertices[i] = list.at(i);
        }

        scc_components->push_back(scc);
    }

    return 0;
}

int ptd_label_vertices(struct ptd_graph *graph) {
    ptd_reset_graph_visited(graph);

    queue<struct ptd_vertex *> q = ptd_enqueue_vertices(graph);
    size_t index = 1;

    while (!q.empty()) {
        struct ptd_vertex *vertex = q.front();
        q.pop();

        if (vertex->edges_length == 0) {
            vertex->index = 0;
        } else {
            vertex->index = index;
            index++;
        }
    }

    graph->is_indexed = true;

    return 0;
}

ptd_strongly_connected_components_t *
ptd_find_strongly_connected_components(struct ptd_graph *graph) {
    ptd_strongly_connected_components_t *components = (ptd_strongly_connected_components_t *) malloc(
            sizeof(ptd_strongly_connected_components_t)
    );

    if (components == NULL) {
        return NULL;
    }

    void **old_data = (void **) calloc(graph->vertices_list->vertices_count, sizeof(void *));

    ptd_label_vertices(graph);
    scc_stack = new stack<struct ptd_vertex *>;
    queue<struct ptd_vertex *> q = ptd_enqueue_vertices(graph);
    queue<struct ptd_vertex *> q_reset;
    ptd_reset_graph_visited(graph);
    scc_index = 0;
    scc_indices = (size_t *) calloc(q.size() + 2, sizeof(size_t));
    low_links = (size_t *) calloc(q.size() + 2, sizeof(size_t));
    scc_on_stack = (bool *) calloc(q.size() + 2, sizeof(bool));
    scc_components = new vector<ptd_strongly_connected_component_t *>();

    size_t vertex_order;

    vertex_order = 0;

    while (!q.empty()) {
        struct ptd_vertex *vertex = q.front();
        q.pop();
        q_reset.push(vertex);
        old_data[vertex_order] = vertex->data;
        vertex_order++;

        if (!vertex->visited) {
            if (strongconnect(vertex) != 0) {
                return NULL;
            }
        }
    }

    size_t non_empty_components = 0;

    for (size_t i = 0; i < scc_components->size(); ++i) {
        ptd_strongly_connected_component_t *c = scc_components->at(i);

        if (c->internal_vertices_length != 0) {
            non_empty_components++;
        }
    }

    components->components_length = non_empty_components;
    components->components = (ptd_strongly_connected_component_t **) calloc(
            components->components_length,
            sizeof(ptd_strongly_connected_component_t *)
    );

    size_t index = 0;

    for (size_t i = 0; i < scc_components->size(); ++i) {
        ptd_strongly_connected_component_t *scc = scc_components->at(i);

        if (scc->internal_vertices_length != 0) {
            ptd_strongly_connected_component_t *scc = scc_components->at(i);
            components->components[index] = scc;
            index++;
        } else {
            free(scc->internal_vertices);
            free(scc);
        }
    }

    for (size_t i = 0; i < components->components_length; ++i) {
        ptd_strongly_connected_component_t *scc = components->components[i];

        for (size_t j = 0; j < scc->internal_vertices_length; ++j) {
            struct ptd_vertex *vertex = scc->internal_vertices[j];
            vertex->data = (ptd_strongly_connected_component_t *) scc;
        }
    }

    for (size_t i = 0; i < components->components_length; ++i) {
        ptd_strongly_connected_component_t *scc = components->components[i];
        set<ptd_strongly_connected_component_t *> external_sccs;
        set<struct ptd_vertex *> external_vertices;

        for (size_t j = 0; j < scc->internal_vertices_length; ++j) {
            struct ptd_vertex *vertex = scc->internal_vertices[j];

            for (size_t k = 0; k < vertex->edges_length; ++k) {
                struct ptd_vertex *child = vertex->edges[k].to;

                if ((ptd_strongly_connected_component_t *) child->data != scc) {
                    external_sccs.insert((ptd_strongly_connected_component_t *) child->data);
                    external_vertices.insert(child);
                }
            }
        }

        scc->external_sccs_length = external_sccs.size();
        scc->external_sccs = (ptd_strongly_connected_component_t **) calloc(scc->external_sccs_length,
                                                                            sizeof(*scc->external_sccs));

        scc->external_vertices_length = external_vertices.size();
        scc->external_vertices = (struct ptd_vertex **) calloc(scc->external_vertices_length,
                                                               sizeof(*scc->external_vertices));

        size_t set_index;

        set_index = 0;

        for (set<ptd_strongly_connected_component_t *>::iterator itr = external_sccs.begin();
             itr != external_sccs.end();
             itr++) {
            scc->external_sccs[set_index] = *itr;
            set_index++;
        }

        set_index = 0;

        for (set<struct ptd_vertex *>::iterator itr = external_vertices.begin();
             itr != external_vertices.end();
             itr++) {
            scc->external_vertices[set_index] = *itr;
            set_index++;
        }
    }

    vertex_order = 0;

    while (!q_reset.empty()) {
        struct ptd_vertex *vertex = q_reset.front();
        q_reset.pop();

        vertex->data = old_data[vertex_order];
        vertex_order++;
    }

    free(old_data);
    delete scc_components;
    free(scc_on_stack);
    scc_on_stack = NULL;
    delete scc_stack;
    free(low_links);
    low_links = NULL;
    free(scc_indices);
    scc_indices = NULL;

    return components;
}

void ptd_strongly_connected_components_destroy(ptd_strongly_connected_components_t *sccs) {
    for (size_t i = 0; i < sccs->components_length; ++i) {
        free(sccs->components[i]->internal_vertices);
        free(sccs->components[i]->external_sccs);
        free(sccs->components[i]->external_vertices);
        memset(sccs->components[i], 0, sizeof(*sccs->components[i]));

        free(sccs->components[i]);
        sccs->components[i] = NULL;
    }

    free(sccs->components);
    sccs->components = NULL;

    free(sccs);
}

ptd_strongly_connected_components_t *ptd_scc_index_topological(
        ptd_strongly_connected_component_t **in, size_t length
) {
    ptd_strongly_connected_component_t **res = (ptd_strongly_connected_component_t **) calloc(
            length, sizeof(*res)
    );

    queue<ptd_strongly_connected_component_t *> q;

    for (size_t i = 0; i < length; ++i) {
        ptd_strongly_connected_component_t *vertex = in[i];
        vertex->index = 0;
        vertex->visited = false;
    }

    for (size_t i = 0; i < length; ++i) {
        ptd_strongly_connected_component_t *vertex = in[i];

        for (size_t j = 0; j < vertex->external_sccs_length; ++j) {
            ptd_strongly_connected_component_t *child = vertex->external_sccs[j];

            child->index++;
        }
    }

    q.push(in[length - 1]);
    size_t idx = 0;

    while (!q.empty()) {
        ptd_strongly_connected_component_t *vertex = q.front();
        q.pop();

        res[idx] = vertex;
        idx++;

        for (size_t i = 0; i < vertex->external_sccs_length; ++i) {
            ptd_strongly_connected_component_t *child = vertex->external_sccs[i];

            child->index--;

            if (child->index == 0 && !child->visited) {
                child->visited = true;
                q.push(child);
            }
        }
    }

    for (size_t i = 0; i < length; ++i) {
        ptd_strongly_connected_component_t *vertex = in[i];
        vertex->index = 0;
        vertex->visited = false;
    }

    ptd_strongly_connected_components_t *r = (ptd_strongly_connected_components_t *) malloc(sizeof(*r));
    r->components = res;
    r->components_length = length;

    return r;
}

int
ptd_order_strongly_connected_components(
        ptd_strongly_connected_components_t *sccs
) {
    struct ptd_graph *graph = sccs->components[0]->internal_vertices[0]->graph;

    if (!graph->is_indexed) {
        ptd_label_vertices(graph);
    }

    ptd_strongly_connected_components_t *sorted = ptd_scc_index_topological(sccs->components, sccs->components_length);
    sccs->components = sorted->components;

    return 0;
}

ptd_phase_type_distribution_t *ptd_find_local_matrix(ptd_strongly_connected_component_t *in) {
    ptd_phase_type_distribution_t *res = (ptd_phase_type_distribution_t *) malloc(sizeof(*res));
    size_t full_length = in->internal_vertices_length + in->external_vertices_length;

    long double **mat = (long double **) calloc(full_length, sizeof(long double *));
    struct ptd_vertex **vertices = (struct ptd_vertex **) calloc(full_length, sizeof(struct ptd_vertex *));
    void **old_data = (void **) calloc(full_length, sizeof(*old_data));
    size_t *old_indices = (size_t *) calloc(full_length, sizeof(*old_indices));

    for (size_t i = 0; i < full_length; ++i) {
        mat[i] = (long double *) calloc(full_length, sizeof(long double));
    }

    for (size_t i = 0; i < in->internal_vertices_length; ++i) {
        old_indices[i] = in->internal_vertices[i]->index;
        vertices[i] = in->internal_vertices[i];
        old_data[i] = in->internal_vertices[i]->data;
        vertices[i]->data = malloc(sizeof(size_t));
        *((size_t *) (vertices)[i]->data) = i;
    }

    for (size_t i = 0; i < in->external_vertices_length; ++i) {
        size_t index = in->internal_vertices_length + i;
        (vertices)[index] = in->external_vertices[i];
        old_indices[index] = (vertices)[index]->index;
        (vertices)[index]->index = index;
        old_data[index] = (vertices)[index]->data;
        (vertices)[index]->data = malloc(sizeof(size_t));
        *((size_t *) (vertices)[index]->data) = index;
    }

    for (size_t i = 0; i < in->internal_vertices_length; ++i) {
        struct ptd_vertex *vertex = (vertices)[i];

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            size_t child_index = *((size_t *) vertex->edges[j].to->data);
            (mat)[i][child_index] += vertex->edges[j].weight / vertex->rate;
            (mat)[i][i] -= vertex->edges[j].weight / vertex->rate;
        }
    }

    for (size_t i = 0; i < in->external_vertices_length; ++i) {
        size_t index = in->internal_vertices_length + i;
        (mat)[index][index] = -1;
    }

    res->length = full_length;
    res->vertices = vertices;
    res->sub_intensity_matrix = mat;
    res->initial_probability_vector = NULL;
    res->memory_allocated = full_length;

    for (size_t i = 0; i < full_length; ++i) {
        free((vertices)[i]->data);
        (vertices)[i]->data = old_data[i];
        (vertices)[i]->index = old_indices[i];
    }

    free(old_indices);
    free(old_data);

    return res;
}

ptd_vertex_group_t *ptd_convert_strongly_connected_component_to_group(ptd_strongly_connected_component_t *scc) {
    DIE_ERROR(1, "Not implemented\n");
}

int ptd_remove_edge(struct ptd_vertex *from, struct ptd_vertex *to) {
    DIE_ERROR(1, "Not implemented\n");
}

struct ptd_graph *
ptd_convert_strongly_connected_components_to_group(struct ptd_graph *graph, ptd_strongly_connected_components_t *sccs) {
    ptd_label_vertices(graph);
    struct ptd_graph *ret = ptd_graph_create(0);

    if (graph == NULL) {
        return NULL;
    }

    struct ptd_vertex **vertices_scc;
    vertices_scc = (struct ptd_vertex **) calloc(graph->vertices_list->vertices_count, sizeof(*vertices_scc));

    if (vertices_scc == NULL) {
        return NULL;
    }

    for (size_t i = 0; i < sccs->components_length; ++i) {
        ptd_strongly_connected_component_t *scc = sccs->components[i];

        // Make the scc into a vertex
        struct ptd_vertex *group = ptd_vertex_create(ret);

        if (group == NULL) {
            return NULL;
        }

        for (size_t j = 0; j < scc->internal_vertices_length; ++j) {
            struct ptd_vertex *vertex = scc->internal_vertices[j];

            vertices_scc[vertex->index] = group;
        }
    }

    set<struct ptd_vertex *> edges;
/*
    for (size_t j = 0; j < scc->internal_vertices_length; ++j) {
        struct ptd_vertex *vertex = scc->internal_vertices[j];

        for (size_t k = 0; k < vertex->edges_length; ++k) {
            edges.insert(vertex->edges[k].to);
        }
    }

    for (size_t k = 0; k < edges.size(); ++k) {
        ptd_add_edge()
    }

    for (size_t i = 0; i < sccs->components_length; ++i) {
        set<vertex_t *> edges_set;

        edges.push_back(set);
    }
*/
    return NULL;
}

double (*reward_function)(struct ptd_vertex *);

static bool keep_zero_rewarded(struct ptd_vertex *vertex) {
    return (reward_function(vertex) == 0);
}

static ptd_avl_tree_t **avl_edges;
static ptd_avl_tree_t **parents;
static struct ptd_vertex **vertices;
static vec_entry_t **states;
static double *rewards;

static inline vector<struct ptd_edge *> avl_tree_as_list(ptd_avl_tree_t *avl_tree) {
    vector<struct ptd_edge *> vec;
    stack<struct avl_node *> s;

    s.push((struct avl_node *) avl_tree->root);

    while (!s.empty()) {
        struct avl_node *v = s.top();
        s.pop();

        if (v == NULL) {
            continue;
        }

        vec.push_back((struct ptd_edge *) v->entry);
        s.push(v->left);
        s.push(v->right);
    }

    return vec;
}

//extern int print_func(struct ptd_vertex *vertex);

// TODO: We should always reset before doing anything else everywhere

int ptd_reward_transform(struct ptd_graph *graph, double (*reward_func)(const struct ptd_vertex *)) {
    ptd_reset_graph_visited(graph);
    queue<struct ptd_vertex *> q = ptd_enqueue_vertices(graph);

    ptd_reset_graph_visited(graph);
    ptd_label_vertices(graph);
    size_t n = q.size();

    avl_edges = (ptd_avl_tree_t **) calloc(n, sizeof(*avl_edges));

    if (avl_edges == NULL) {
        return 1;
    }

    parents = (ptd_avl_tree_t **) calloc(n, sizeof(*parents));
    vertices = (struct ptd_vertex **) calloc(n, sizeof(*vertices));
    states = (vec_entry_t **) calloc(n, sizeof(*states));
    rewards = (double *) calloc(n, sizeof(*rewards));

    while (!q.empty()) {
        struct ptd_vertex *vertex = q.front();
        q.pop();

        vertices[vertex->index] = vertex;
        avl_edges[vertex->index] = ptd_avl_tree_create(1);

        if (avl_edges[vertex->index] == NULL) {
            return -1;
        }

        parents[vertex->index] = ptd_avl_tree_create(1);

        if (parents[vertex->index] == NULL) {
            // TODO: Free rest
            return -1;
        }

        states[vertex->index] = (vec_entry_t *) calloc(1, sizeof(*(states[vertex->index])));
        states[vertex->index][0] = vertex->index;

        double reward;

        if (vertex->index == 0 || vertex->edges_length == 0) {
            reward = 1;
        } else {
            reward = reward_func(vertex);
        }

        rewards[vertex->index] = reward;
    }

    for (size_t i = 0; i < n; ++i) {
        struct ptd_vertex *vertex = vertices[i];

        // We make the weight be the probability of transitioning.
        // Store the new reward for later
        rewards[i] /= vertex->rate;

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            struct ptd_vertex *child_vertex = vertex->edges[j].to;
            size_t child_index = child_vertex->index;

            ptd_avl_tree_edge_insert_or_increment(
                    avl_edges[i],
                    states[child_index],
                    child_vertex,
                    vertex->edges[j].weight / vertex->rate
            );

            ptd_avl_tree_edge_insert_or_increment(
                    parents[child_index],
                    states[i],
                    vertex,
                    vertex->edges[j].weight / vertex->rate
            );
        }
    }

    for (size_t i = 0; i < n; ++i) {
        double reward = rewards[i];

        if (reward == 0) {
            vector<struct ptd_edge *> outgoing_edges = avl_tree_as_list(avl_edges[i]);
            vector<struct ptd_edge *> ingoing_edges = avl_tree_as_list(parents[i]);

            for (size_t p = 0; p < ingoing_edges.size(); ++p) {
                struct ptd_edge *ingoing_edge = ingoing_edges[p];
                size_t parent_index = ingoing_edge->to->index;

                for (size_t c = 0; c < outgoing_edges.size(); ++c) {
                    struct ptd_edge *outgoing_edge = outgoing_edges[c];
                    size_t child_index = outgoing_edge->to->index;

                    if (child_index != parent_index) {
                        // Add new child to parent
                        // Weight has been set to probability earlier
                        ptd_avl_tree_edge_insert_or_increment(
                                avl_edges[parent_index],
                                states[child_index],
                                outgoing_edge->to,
                                ingoing_edge->weight * outgoing_edge->weight
                        );

                        ptd_avl_tree_edge_insert_or_increment(
                                parents[child_index],
                                states[parent_index],
                                ingoing_edge->to,
                                ingoing_edge->weight * outgoing_edge->weight
                        );
                    }
                }

                ptd_avl_tree_edge_remove(avl_edges[parent_index], states[i]);
            }

            for (size_t c = 0; c < outgoing_edges.size(); ++c) {
                struct ptd_edge *outgoing_edge = outgoing_edges[c];
                size_t child_index = outgoing_edge->to->index;

                ptd_avl_tree_edge_remove(parents[child_index], states[i]);
            }
        }
    }

    for (size_t i = 0; i < n; ++i) {
        struct ptd_vertex *vertex = vertices[i];
        vertex->rate = 0;
        vertex->edges_length = 0;
        double reward = rewards[i];

        if (reward == 0) {
            continue;
        }

        vector<struct ptd_edge *> outgoing_edges = avl_tree_as_list(avl_edges[i]);


        fprintf(stderr, "\n=============\n\nGoing to vertex with %zu edges: ", outgoing_edges.size());
        //print_func(vertex);

        for (size_t j = 0; j < outgoing_edges.size(); ++j) {
            fprintf(stderr, "Adding edge between: ");
            //print_func(vertex);
            fprintf(stderr, "and: ");
            //print_func(outgoing_edges[j]->to);
            ptd_add_edge(
                    vertex,
                    outgoing_edges[j]->to,
                    outgoing_edges[j]->weight / reward
            );
        }

        fprintf(stderr, "Resulting in: ");
        //print_func(vertex);
    }

    return 0;
}

void ptd_phase_type_distribution_destroy(ptd_phase_type_distribution_t *ptd) {
    for (size_t i = 0; i < ptd->memory_allocated; ++i) {
        free(ptd->sub_intensity_matrix[i]);

        ptd->sub_intensity_matrix[i] = NULL;
    }

    free(ptd->vertices);
    free(ptd->sub_intensity_matrix);
    free(ptd->initial_probability_vector);

    ptd->vertices = NULL;
    ptd->sub_intensity_matrix = NULL;
    ptd->initial_probability_vector = NULL;

    ptd->memory_allocated = 0;
    ptd->length = 0;

    free(ptd);
}

ptd_phase_type_distribution_t *ptd_graph_as_phase_type_distribution(struct ptd_graph *graph) {
    if (!graph->is_indexed) {
        ptd_label_vertices(graph);
    }

    queue<struct ptd_vertex *> q = ptd_enqueue_vertices(graph);

    ptd_phase_type_distribution_t *res = (ptd_phase_type_distribution_t *) malloc(sizeof(*res));

    if (res == NULL) {
        return NULL;
    }

    res->length = 0;

    size_t size = (size_t) q.size();

    res->memory_allocated = size;
    res->vertices = (struct ptd_vertex **) calloc(size, sizeof(struct ptd_vertex *));

    if (res->vertices == NULL) {
        free(res);
        return NULL;
    }

    res->initial_probability_vector = (long double *) calloc(size, sizeof(long double));

    if (res->initial_probability_vector == NULL) {
        free(res->vertices);
        free(res);
        return NULL;
    }

    res->sub_intensity_matrix = (long double **) calloc(size, sizeof(long double *));

    if (res->sub_intensity_matrix == NULL) {
        free(res->initial_probability_vector);
        free(res->vertices);
        free(res);
        return NULL;
    }

    for (size_t i = 0; i < size; ++i) {
        res->sub_intensity_matrix[i] = (long double *) calloc(size, sizeof(long double));

        if ((res->sub_intensity_matrix)[i] == NULL) {
            for (size_t j = 0; j < i; ++j) {
                free(res->sub_intensity_matrix[j]);
            }

            free(res->sub_intensity_matrix);
            free(res->initial_probability_vector);
            free(res->vertices);
            free(res);
            return NULL;
        }

    }

    while (!q.empty()) {
        struct ptd_vertex *vertex = q.front();
        q.pop();

        if (vertex->index == 0) {
            continue;
        }

        if (vertex == graph->start_vertex) {
            for (size_t i = 0; i < vertex->edges_length; ++i) {
                struct ptd_edge edge = vertex->edges[i];

                if (edge.to->index != 0) {
                    res->initial_probability_vector[edge.to->index - 1 - 1] = edge.weight;
                }
            }

            continue;
        }

        res->vertices[vertex->index - 2] = vertex;
        res->length++;

        for (size_t i = 0; i < vertex->edges_length; ++i) {
            struct ptd_edge edge = vertex->edges[i];

            if (edge.to->index != 0) {
                res->sub_intensity_matrix[vertex->index - 2][edge.to->index - 2] += edge.weight;
            }

            res->sub_intensity_matrix[vertex->index - 2][vertex->index - 2] -= edge.weight;
        }
    }

    return res;
}

static int *topo_values;

int ptd_index_topological(struct ptd_graph *graph) {
    queue<struct ptd_vertex *> topological_queue;
    queue<struct ptd_vertex *> q;

    q = ptd_enqueue_vertices(graph);
    topo_values = (int *) calloc(q.size(), sizeof(*topo_values));

    if (topo_values == NULL) {
        return -1;
    }

    size_t index = 0;

    while (!q.empty()) {
        struct ptd_vertex *vertex = q.front();
        q.pop();

        vertex->index = index;

        index++;
    }

    q = ptd_enqueue_vertices(graph);

    while (!q.empty()) {
        struct ptd_vertex *vertex = q.front();
        q.pop();

        for (size_t i = 0; i < vertex->edges_length; ++i) {
            topo_values[vertex->edges[i].to->index]++;
        }
    }

    ptd_reset_graph_visited(graph);

    q.push(graph->start_vertex);

    while (!q.empty()) {
        struct ptd_vertex *vertex = q.front();
        q.pop();

        topological_queue.push(vertex);

        for (size_t i = 0; i < vertex->edges_length; ++i) {
            struct ptd_vertex *child = vertex->edges[i].to;

            topo_values[child->index] -= 1;

            if (topo_values[child->index] == 0 && !child->visited) {
                child->visited = true;
                q.push(child);
            }
        }
    }

    index = 1;

    while (!topological_queue.empty()) {
        struct ptd_vertex *vertex = topological_queue.front();
        topological_queue.pop();

        if (vertex->edges_length == 0) {
            vertex->index = 0;
            continue;
        }

        vertex->index = index;
        index++;
    }

    free(topo_values);
    graph->is_indexed = true;

    return 0;
}

int ptd_index_invert(struct ptd_graph *graph) {
    queue<struct ptd_vertex *> q = ptd_enqueue_vertices(graph);
    size_t largest_index = 0;

    while (!q.empty()) {
        struct ptd_vertex *vertex = q.front();
        q.pop();

        if (vertex->index > largest_index) {
            largest_index = vertex->index;
        }
    }

    q = ptd_enqueue_vertices(graph);

    while (!q.empty()) {
        struct ptd_vertex *vertex = q.front();
        q.pop();

        if (vertex->edges_length == 0) {
            vertex->index = 0;
            continue;
        }

        if (vertex == graph->start_vertex) {
            vertex->index = 1;
            continue;
        }

        vertex->index = largest_index - vertex->index + 2;
    }

    return 0;
}

int ptd_vertex_to_s(struct ptd_vertex *vertex, char *buffer, size_t buffer_length) {
    memset(buffer, '\0', buffer_length);

    char *build = (char *) calloc(buffer_length, sizeof(char));

    for (size_t i = 0; i < vertex->graph->state_length; ++i) {
        if (i == 0) {
            snprintf(build, buffer_length, "%s%zu", buffer, vertex->state[i]);
        } else {
            snprintf(build, buffer_length, "%s %zu", buffer, vertex->state[i]);
        }

        strncpy(buffer, build, buffer_length);
    }

    free(build);

    return 0;
}

/*
 * Models
 */


static ptd_avl_tree_t *avl_tree = NULL;
static struct ptd_graph *kingman_graph = NULL;

static int make_kingman(struct ptd_vertex *vertex) {
    vec_entry_t *state = (vec_entry_t *) calloc(kingman_graph->state_length, sizeof(vec_entry_t));
    memcpy(state, vertex->state, kingman_graph->state_length * sizeof(vec_entry_t));

    for (size_t i = 0; i < kingman_graph->state_length; ++i) {
        for (size_t j = i; j < kingman_graph->state_length; ++j) {
            double weight;

            if (i == j) {
                if (state[i] < 2) {
                    continue;
                }

                weight = state[i] * (state[i] - 1) / 2;
            } else {
                if (state[i] < 1 || state[j] < 1) {
                    continue;
                }

                weight = state[i] * state[j];
            }

            state[i]--;
            state[j]--;
            state[(i + j + 2) - 1]++;

            struct ptd_vertex *child = ptd_avl_tree_vertex_find(avl_tree, state);

            if (child == NULL) {
                vec_entry_t *child_state = (vec_entry_t *) calloc(kingman_graph->state_length, sizeof(vec_entry_t));

                if (child_state == NULL) {
                    return -1;
                }

                memcpy(child_state, state, kingman_graph->state_length * sizeof(vec_entry_t));

                child = ptd_vertex_create_state(kingman_graph, child_state);

                if (ptd_avl_tree_vertex_insert(avl_tree, child_state, child)) {
                    return -1;
                }
            }

            state[i]++;
            state[j]++;
            state[(i + j + 2) - 1]--;

            ptd_add_edge(vertex, child, weight);
        }
    }

    free(state);

    return 0;
}

struct ptd_graph *ptd_model_kingman(size_t n) {
    avl_tree = ptd_avl_tree_create(n);

    if (avl_tree == NULL) {
        return NULL;
    }

    kingman_graph = ptd_graph_create(n);

    if (kingman_graph == NULL) {
        ptd_avl_tree_vertex_destroy(avl_tree);
        return NULL;
    }

    struct ptd_vertex *initial = ptd_vertex_create(kingman_graph);

    if (initial == NULL) {
        ptd_graph_destroy(kingman_graph);
        ptd_avl_tree_vertex_destroy(avl_tree);
    }

    initial->state[0] = n;

    if (ptd_add_edge(kingman_graph->start_vertex, initial, 1) != 0) {
        ptd_graph_destroy(kingman_graph);
        ptd_avl_tree_vertex_destroy(avl_tree);

        return NULL;
    }

    if (ptd_visit_vertices(kingman_graph, make_kingman, false) != 0) {
        ptd_graph_destroy(kingman_graph);
        ptd_avl_tree_vertex_destroy(avl_tree);

        return NULL;
    }

    ptd_avl_tree_vertex_destroy(avl_tree);

    return kingman_graph;
}

/*
 * State has two islands and two loci. For two loci we can have the states
 * A-A, A-, -A, for the state of recombination. Since we have two islands,
 * this gives us 6 states.
 *
 * State 0,1,2: island 1
 * State 3,4,5: island 2
 */

static struct ptd_graph *two_island_two_loci_recomb_graph;

static double recomb_rate;
static size_t N[2];
static double mig[2];

#define IDX_COMB 0
#define IDX_L 1
#define IDX_R 2

static struct ptd_vertex *add_child(vec_entry_t *state) {
    struct ptd_vertex *child = ptd_avl_tree_vertex_find(avl_tree, state);

    if (child == NULL) {
        vec_entry_t *child_state = (vec_entry_t *) calloc(
                two_island_two_loci_recomb_graph->state_length,
                sizeof(vec_entry_t)
        );

        if (child_state == NULL) {
            return NULL;
        }

        memcpy(
                child_state, state,
                two_island_two_loci_recomb_graph->state_length * sizeof(vec_entry_t)
        );

        child = ptd_vertex_create_state(
                two_island_two_loci_recomb_graph, child_state
        );

        if (ptd_avl_tree_vertex_insert(avl_tree, child_state, child)) {
            return NULL;
        }
    }

    return child;
}

static int visit_two_island_two_loci_recomb(struct ptd_vertex *vertex) {
    // If there is only one left
    if (vertex->state[0] + vertex->state[1] +
        vertex->state[2] + vertex->state[3] +
        vertex->state[4] + vertex->state[5] <= 1) {
        return 0;
    }

    vec_entry_t *state = (vec_entry_t *) calloc(
            two_island_two_loci_recomb_graph->state_length,
            sizeof(vec_entry_t)
    );

    memcpy(
            state, vertex->state,
            two_island_two_loci_recomb_graph->state_length * sizeof(vec_entry_t)
    );

    // Migration
    {
        for (size_t island = 0; island < 2; ++island) {
            for (size_t i = 0; i < 3; ++i) {
                size_t pos = island * 3 + i;
                size_t other_pos = (1 - island) * 3 + i;

                if (state[pos] > 0) {
                    state[pos]--;
                    state[other_pos]++;

                    struct ptd_vertex *child = add_child(state);

                    state[pos]++;
                    state[other_pos]--;

                    ptd_add_edge(vertex, child, mig[island] * N[island]);
                }
            }
        }
    }

    // Recombination
    {
        for (size_t island = 0; island < 2; ++island) {
            double rate = recomb_rate * N[island];

            // From combined to split
            if (state[island * 3 + IDX_COMB] > 0) {
                state[island * 3 + IDX_COMB]--;
                state[island * 3 + IDX_L]++;
                state[island * 3 + IDX_R]++;

                struct ptd_vertex *child = add_child(state);

                state[island * 3 + IDX_COMB]++;
                state[island * 3 + IDX_L]--;
                state[island * 3 + IDX_R]--;

                ptd_add_edge(vertex, child, state[island * 3 + IDX_COMB] * rate);

            }

            // From split to combined
            if (state[island * 3 + IDX_L] > 0 && state[island * 3 + IDX_R] > 0) {
                state[island * 3 + IDX_COMB]++;
                state[island * 3 + IDX_L]--;
                state[island * 3 + IDX_R]--;

                struct ptd_vertex *child = add_child(state);

                state[island * 3 + IDX_COMB]--;
                state[island * 3 + IDX_L]++;
                state[island * 3 + IDX_R]++;

                ptd_add_edge(
                        vertex, child,
                        state[island * 3 + IDX_R] * state[island * 3 + IDX_L] * rate
                );
            }
        }
    }

    // Coalescence
    {
        for (size_t island = 0; island < 2; ++island) {
            // Combined
            if (state[island * 3 + IDX_COMB] > 1) {
                size_t combinations = state[island * 3 + IDX_COMB] * (state[island * 3 + IDX_COMB] - 1) / 2;

                state[island * 3 + IDX_COMB]--;
                struct ptd_vertex *child = add_child(state);
                state[island * 3 + IDX_COMB]++;

                ptd_add_edge(vertex, child, combinations * N[island]);
            }

            // Left and combined
            if (state[island * 3 + IDX_COMB] > 0 && state[island * 3 + IDX_L] > 0) {
                size_t combinations = state[island * 3 + IDX_COMB] * state[island * 3 + IDX_L];

                state[island * 3 + IDX_L]--;

                struct ptd_vertex *child = add_child(state);

                state[island * 3 + IDX_L]++;

                ptd_add_edge(vertex, child, combinations * N[island]);
            }

            // Right and combined
            if (state[island * 3 + IDX_COMB] > 0 && state[island * 3 + IDX_R] > 0) {
                size_t combinations = state[island * 3 + IDX_COMB] * state[island * 3 + IDX_R];

                state[island * 3 + IDX_R]--;
                struct ptd_vertex *child = add_child(state);
                state[island * 3 + IDX_R]++;

                ptd_add_edge(vertex, child, combinations * N[island]);
            }

            // Left and left
            if (state[island * 3 + IDX_L] > 1) {
                size_t combinations = state[island * 3 + IDX_L] * (state[island * 3 + IDX_L] - 1) / 2;

                state[island * 3 + IDX_L]--;
                struct ptd_vertex *child = add_child(state);
                state[island * 3 + IDX_L]++;

                ptd_add_edge(vertex, child, combinations * N[island]);
            }

            // Right and right
            if (state[island * 3 + IDX_R] > 1) {
                size_t combinations = state[island * 3 + IDX_R] * (state[island * 3 + IDX_R] - 1) / 2;

                state[island * 3 + IDX_R]--;
                struct ptd_vertex *child = add_child(state);
                state[island * 3 + IDX_R]++;

                ptd_add_edge(vertex, child, combinations * N[island]);
            }
        }
    }

    free(state);

    return 0;
}


struct ptd_graph *ptd_model_two_island_two_loci_recomb(
        size_t n1,
        size_t n2,
        size_t effective_population_size1,
        size_t effective_population_size2,
        double migration_rate1,
        double migration_rate2,
        double recombination_rate
) {
    N[0] = effective_population_size1;
    N[1] = effective_population_size2;
    mig[0] = migration_rate1;
    mig[1] = migration_rate2;
    recomb_rate = recombination_rate;

    two_island_two_loci_recomb_graph = ptd_graph_create(6);
    avl_tree = ptd_avl_tree_create(6);

    struct ptd_vertex *initial = ptd_vertex_create(two_island_two_loci_recomb_graph);
    initial->state[0] = n1;
    initial->state[3] = n2;

    ptd_add_edge(
            two_island_two_loci_recomb_graph->start_vertex,
            initial,
            1
    );

    ptd_visit_vertices(
            two_island_two_loci_recomb_graph,
            visit_two_island_two_loci_recomb,
            false
    );

    ptd_avl_tree_vertex_destroy(avl_tree);

    return two_island_two_loci_recomb_graph;
}

/*
 * Visit algorithms
 */

size_t alloc_size;
char *alloc_value;

static int visit_alloc(struct ptd_vertex *vertex) {
    vertex->data = malloc(alloc_size);
    memcpy((char *) vertex->data, alloc_value, alloc_size);

    return (vertex->data == NULL);
}

int ptd_visit_alloc(struct ptd_graph *graph, size_t size, char *value) {
    alloc_size = size;
    alloc_value = value;

    return ptd_visit_vertices(graph, visit_alloc, true);
}

static int visit_probability(struct ptd_vertex *vertex) {
    long double my_prob;

    if (vertex == vertex->graph->start_vertex) {
        *((long double *) vertex->data) = 1;
    }

    my_prob = *((long double *) vertex->data);

    struct ptd_edge *edges = vertex->edges;

    for (size_t i = 0; i < vertex->edges_length; ++i) {
        struct ptd_edge edge = edges[i];

        long double trans_prob = edge.weight / vertex->rate;

        *((long double *) edge.to->data) += my_prob * trans_prob;
    }

    return 0;
}

int ptd_visit_assign_probability(struct ptd_graph *graph) {
    // TODO: do only if graph is acyclic
    ptd_index_topological(graph);

    return ptd_visit_vertices(graph, visit_probability, true);
}

static int compute_expected(struct ptd_vertex *vertex) {
    long double my_prob = *((long double *) vertex->data);

    long double my_expected;

    if (vertex->rate != 0) {
        my_expected = my_prob * reward_function(vertex) / vertex->rate;
    } else {
        my_expected = 0;
    }

    *((long double *) vertex->data) = my_expected;

    return 0;
}

int ptd_visit_assign_expectation(struct ptd_graph *graph, double (*reward)(struct ptd_vertex *)) {
    reward_function = reward;
    return ptd_visit_vertices(graph, compute_expected, true);
}

long double reduce_ld_sum;

static int reduce_ld(struct ptd_vertex *vertex) {
    long double data = *((long double *) vertex->data);
    reduce_ld_sum += data;

    return 0;
}

long double ptd_visit_reduce_sum_long_double(struct ptd_graph *graph) {
    reduce_ld_sum = 0;
    ptd_visit_vertices(graph, reduce_ld, true);

    return reduce_ld_sum;
}

/*
 * Properties
 */
int ptd_expected_value(long double *expected, struct ptd_graph *graph, double (*reward)(struct ptd_vertex *)) {
    long double *start_value = (long double *) malloc(sizeof(long double));
    *start_value = 0;
    ptd_visit_alloc(graph, sizeof(long double), (char *) start_value);
    ptd_visit_assign_probability(graph);
    ptd_visit_assign_expectation(graph, reward);
    *expected = ptd_visit_reduce_sum_long_double(graph);

    return 0;
}

struct cov {
    long double prob;
    long double desc;
};

static int visit_desc(struct ptd_vertex *vertex) {
    long double desc = 0;

    struct ptd_edge *edges = vertex->edges;

    for (size_t i = 0; i < vertex->edges_length; ++i) {
        struct ptd_edge edge = edges[i];

        long double trans_prob = edge.weight / vertex->rate;

        desc += trans_prob * ((struct cov *) edge.to->data)->desc;
    }


    char buf[1024] = {0};
    ptd_vertex_to_s(vertex, buf, 1024);
    fprintf(stderr, "I am %s, I have B %Lf\n", buf, desc);

    if (vertex->rate != 0) {
        desc += reward_function(vertex) / vertex->rate;
    }

    memset(buf, 0, sizeof(buf) * sizeof(char));

    ptd_vertex_to_s(vertex, buf, 1024);
    fprintf(stderr, "I am %s, I have A %Lf\n", buf, desc);

    ((struct cov *) vertex->data)->desc = desc;

    return 0;
}

long double cov_sum;

static int visit_cov(struct ptd_vertex *vertex) {
    struct cov data = *((struct cov *) vertex->data);
    cov_sum += data.prob * reward_function(vertex) * data.desc;

    return 0;
}

int ptd_covariance(
        long double *covariance, struct ptd_graph *graph,
        double (*reward_1)(struct ptd_vertex *), double (*reward_2)(struct ptd_vertex *)
) {
    struct cov *start_value = (struct cov *) malloc(sizeof(struct cov));
    start_value->prob = 0;
    start_value->desc = 0;
    cov_sum = 0;
    ptd_visit_alloc(graph, sizeof(struct cov), (char *) start_value);
    ptd_visit_assign_probability(graph);
    ptd_index_topological(graph);
    ptd_index_invert(graph);
    reward_function = reward_2;
    ptd_visit_vertices(graph, visit_desc, true);
    return 0;
    reward_function = reward_1;
    ptd_visit_vertices(graph, visit_cov, true);
    reward_function = reward_1;
    ptd_visit_vertices(graph, visit_desc, true);
    reward_function = reward_2;
    ptd_visit_vertices(graph, visit_cov, true);

    long double exp1;
    long double exp2;
    ptd_expected_value(&exp1, graph, reward_1);
    ptd_expected_value(&exp2, graph, reward_2);

    cov_sum -= exp1 * exp2;

    *covariance = cov_sum;

    return 0;
}

static int set_data_as_int(struct ptd_vertex *vertex) {
    vertex->data = malloc(sizeof(long double));
    *((long double *) vertex->data) = 0;

    return 0;
}

static bool keep_all(struct ptd_vertex *vertex) {
    return true;
}

double ptd_circular_exp(struct ptd_graph *graph, double (*reward)(struct ptd_vertex *)) {
    ptd_visit_vertices(graph, set_data_as_int, true);

    *((long double *) graph->start_vertex->data) = 1;

    ptd_strongly_connected_components_t *sccs = ptd_find_strongly_connected_components(
            graph
    );

    ptd_order_strongly_connected_components(sccs);
    ptd_strongly_connected_component_t **vertices = sccs->components;

    for (size_t i = 0; i < sccs->components_length; ++i) {
        DEBUG_PRINT("At scc %zu size %zu\n", i, sccs->components[i]->internal_vertices_length);
        ptd_strongly_connected_component_t *v = vertices[i];

        long double **mat;
        struct ptd_vertex **vs;
        size_t length;

        DEBUG_PRINT("As ptd...\n");
        ptd_phase_type_distribution_t *ptd = ptd_find_local_matrix(v);
        DEBUG_PRINT("As ptd\n");
        mat = ptd->sub_intensity_matrix;
        vs = ptd->vertices;
        length = ptd->length;

        fprintf(stderr, "A");
        long double *ipv = (long double *) calloc(length, sizeof(*ipv));

        for (size_t k = 0; k < length; ++k) {
            ipv[k] = *((long double *) vs[k]->data);
            *((long double *) vs[k]->data) = 0;
        }

        fprintf(stderr, "B\n");

        ptd->initial_probability_vector = ipv;

        if (length == 1) {

            fprintf(stderr, "F1\n");
            ptd_phase_type_distribution_destroy(ptd);
            fprintf(stderr, "F2\n");
            continue;
        }

        DEBUG_PRINT("Making mat...\n");
        void *full = create_matrix(mat, length);

        DEBUG_PRINT("Making mat %zu\n", length);
        DEBUG_PRINT("Making inv...\n");
        void *inv = matrix_invert(full, length);
        DEBUG_PRINT("inv\n");

        for (size_t k = 0; k < length; ++k) {
            for (size_t j = 0; j < length; ++j) {
                *((long double *) vs[j]->data) += ipv[k] * -matrix_get(inv, k, j);
            }
        }
        //gsl_matrix_free(full);
        //gsl_matrix_free(inv);

        ptd_phase_type_distribution_destroy(ptd);
    }

    *((long double *) graph->start_vertex->data) = 0;
    ptd_visit_assign_expectation(graph, reward);

    /*{
        ptd_label_vertices(graph);

        ptd_phase_type_distribution_t *ptd = ptd_graph_as_phase_type_distribution(graph);
        gsl_matrix *full = gsl_matrix_alloc(ptd->length, ptd->length);

        for (size_t k = 0; k < ptd->length; ++k) {
            for (size_t j = 0; j < ptd->length; ++j) {
                gsl_matrix_set(full, k, j, (double) ptd->sub_intensity_matrix[k][j]);
            }
        }

        gsl_matrix *inv = matrix_invert(full, ptd->length);
        gsl_matrix_free(full);

        fprintf(stderr, "TRUE RES: ");
        double e = 0;

        for (size_t k = 0; k < ptd->length; ++k) {
            e += -matrix_get(inv, 0, k);
        }
        fprintf(stderr, "%f", e);
        fprintf(stderr, "\n");
    }*/

    return ptd_visit_reduce_sum_long_double(graph);


    ptd_strongly_connected_components_destroy(sccs);
}

ptd_desc_multipliers_t *ptd_cyclic_descendant_multipliers(struct ptd_graph *graph) {
    if (!graph->is_indexed) {
        ptd_label_vertices(graph);
    }

    ptd_desc_multiplier_t **desc_multipliers = (ptd_desc_multiplier_t **) calloc(
            graph->vertices_list->vertices_count,
            sizeof(*desc_multipliers)
    );

    size_t *desc_length = (size_t *) calloc(graph->vertices_list->vertices_count, sizeof(*desc_length));

    ptd_visit_vertices(graph, set_data_as_int, true);

    *((double *) graph->start_vertex->data) = 1;

    ptd_strongly_connected_components_t *sccs = ptd_find_strongly_connected_components(
            graph
    );

    ptd_order_strongly_connected_components(sccs);
    ptd_strongly_connected_component_t **vertices = sccs->components;

    for (size_t i = 0; i < sccs->components_length; ++i) {
        ptd_strongly_connected_component_t *v = vertices[i];

        long double **mat;
        struct ptd_vertex **vs;
        size_t length;

        ptd_phase_type_distribution_t *ptd = ptd_find_local_matrix(v);
        mat = ptd->sub_intensity_matrix;
        vs = ptd->vertices;
        length = ptd->length;

        if (length == 1) {
            ptd_phase_type_distribution_destroy(ptd);
            continue;
        }

        void *full = create_matrix(mat, length);

        void *inv = matrix_invert(full, length);

        for (size_t k = 0; k < v->internal_vertices_length; ++k) {
            if (vs[k]->index != 0) {
                desc_length[vs[k]->index] = length;
                desc_multipliers[vs[k]->index] = (ptd_desc_multiplier_t *) calloc(
                        length, sizeof(*(desc_multipliers[vs[k]->index]))
                );

                for (size_t j = 0; j < length; ++j) {
                    desc_multipliers[vs[k]->index][j].external = j >= v->internal_vertices_length;
                    desc_multipliers[vs[k]->index][j].vertex = vs[j];
                    if (!desc_multipliers[vs[k]->index][j].external) {
                        desc_multipliers[vs[k]->index][j].multiplier =
                                -matrix_get(inv, k, j) / (double) desc_multipliers[vs[k]->index][j].vertex->rate;
                    } else {
                        desc_multipliers[vs[k]->index][j].multiplier =
                                -matrix_get(inv, k, j);
                    }
                }
            } else {
                desc_length[vs[k]->index] = 0;
            }
        }

        //gsl_matrix_free(inv);
        //gsl_matrix_free(full);
        ptd_phase_type_distribution_destroy(ptd);
    }

    ptd_desc_multipliers_t *res = (ptd_desc_multipliers_t *) malloc(sizeof(*res));
    res->desc_length = desc_length;
    res->desc_multipliers = desc_multipliers;
    res->ordered = sccs;

    return res;
}

double *
ptd_cyclic_desc(struct ptd_graph *graph, ptd_desc_multipliers_t *multipliers, double (*reward)(struct ptd_vertex *)) {
    ptd_desc_multiplier_t **desc_multipliers = multipliers->desc_multipliers;
    size_t *desc_length = multipliers->desc_length;
    ptd_strongly_connected_components_t *ordered = multipliers->ordered;

    double *desc_value = (double *) calloc(graph->vertices_list->vertices_count, sizeof(*desc_value));

    for (size_t ii = 0; ii < ordered->components_length; ++ii) {
        size_t i = ordered->components_length - 1 - ii;
        size_t length = ordered->components[i]->internal_vertices_length;
        struct ptd_vertex **vertices = ordered->components[i]->internal_vertices;

        for (size_t k = 0; k < length; ++k) {
            if (vertices[k]->index != 0) {
                for (size_t j = 0; j < desc_length[vertices[k]->index]; ++j) {

                    if (desc_multipliers[vertices[k]->index][j].vertex->index == 0) {
                        continue;
                    }

                    if (!desc_multipliers[vertices[k]->index][j].external) {
                        double multiplier = desc_multipliers[vertices[k]->index][j].multiplier;
                        double reward_value = reward(desc_multipliers[vertices[k]->index][j].vertex);
                        double rate = (double) desc_multipliers[vertices[k]->index][j].vertex->rate;
                        desc_value[vertices[k]->index] += multiplier * reward_value / rate;
                    } else {
                        double multiplier = desc_multipliers[vertices[k]->index][j].multiplier;
                        desc_value[vertices[k]->index] +=
                                multiplier * desc_value[desc_multipliers[vertices[k]->index][j].vertex->index];
                    }
                }
            }
        }
    }

    desc_value[1] = 0;
/*
    for (size_t i = 0; i < graph->internal_vertices_length; ++i) {
        free(desc_multipliers[i]);
    }

    free(desc_multipliers);
    ptd_ordered_sccs_destroy(ordered);
    ptd_strongly_connected_components_destroy(sccs);
    free(desc_length);*/

    return desc_value;
}


/*** NEW ***/

void ptd_directed_graph_destroy(struct ptd_directed_graph *graph) {
    for (size_t i = 0; i < graph->vertices_length; ++i) {
        ptd_directed_vertex_destroy(graph->vertices[i]);
    }

    free(graph->vertices);
    graph->vertices = NULL;
    free(graph);
}

int ptd_directed_vertex_add(struct ptd_directed_graph *graph, struct ptd_directed_vertex *vertex) {
    bool is_power_of_2 = (graph->vertices_length & (graph->vertices_length - 1)) == 0;

    if (is_power_of_2) {
        size_t new_length = graph->vertices_length == 0 ? 1 : graph->vertices_length * 2;

        if ((graph->vertices = (struct ptd_directed_vertex **) realloc(
                graph->vertices, new_length *
                                 sizeof(struct ptd_directed_vertex *))
            ) == NULL) {
            return -1;
        }
    }

    vertex->graph = graph;

    graph->vertices[graph->vertices_length] = vertex;
    vertex->index = graph->vertices_length;
    graph->vertices_length++;

    return 0;
}

int ptd_directed_graph_add_edge(struct ptd_directed_vertex *vertex, struct ptd_directed_edge *edge) {
    /*bool is_power_of_2 = (vertex->edges_length & (vertex->edges_length - 1)) == 0;

    if (is_power_of_2) {
        size_t new_length = vertex->edges_length == 0 ? 1 : vertex->edges_length * 2;

        if ((vertex->edges = realloc(vertex->edges, new_length * edge_size)) == NULL) {
            return -1;
        }
    }

    size_t ptr_location = vertex->edges_length * edge_size;
    ((char*)(vertex->edges))[ptr_location] = edge;*/
    bool is_power_of_2 = (vertex->edges_length & (vertex->edges_length - 1)) == 0;

    if (is_power_of_2) {
        size_t new_length = vertex->edges_length == 0 ? 1 : vertex->edges_length * 2;

        if ((vertex->edges = (struct ptd_directed_edge **) realloc(
                vertex->edges,
                new_length * sizeof(struct ptd_directed_edge *))
            ) == NULL) {
            return -1;
        }
    }

    vertex->edges[vertex->edges_length] = edge;
    vertex->edges_length++;

    return 0;
}

void ptd_directed_vertex_destroy(struct ptd_directed_vertex *vertex) {
    for (size_t i = 0; i < vertex->edges_length; ++i) {
        free(vertex->edges[i]);
    }

    free(vertex->edges);
    vertex->edges = NULL;
    free(vertex);
}

struct ptd_ph_graph *ptd_ph_graph_create(size_t state_length) {
    struct ptd_ph_graph *graph = (struct ptd_ph_graph *) malloc(sizeof(*graph));
    graph->vertices_length = 0;
    graph->state_length = state_length;
    graph->vertices = NULL;
    graph->starting_vertex = ptd_ph_vertex_create(graph);

    return graph;
}

void ptd_ph_graph_destroy(struct ptd_ph_graph *graph) {
    for (size_t i = 0; i < graph->vertices_length; ++i) {
        ptd_ph_vertex_destroy(graph->vertices[i]);
    }

    free(graph->vertices);
    memset(graph, 0, sizeof(*graph));
    free(graph);
}

struct ptd_ph_vertex *ptd_ph_vertex_create(struct ptd_ph_graph *graph) {
    int *state = (int *) calloc(graph->state_length, sizeof(*state));

    return ptd_ph_vertex_create_state(graph, state);
}

struct ptd_ph_vertex *ptd_ph_vertex_create_state(struct ptd_ph_graph *graph, int *state) {
    struct ptd_ph_vertex *vertex = (struct ptd_ph_vertex *) malloc(sizeof(*vertex));
    vertex->graph = graph;
    vertex->edges_length = 0;
    vertex->state = state;
    vertex->edges = NULL;
    ptd_directed_vertex_add(
            (struct ptd_directed_graph *) graph,
            (struct ptd_directed_vertex *) vertex
    );

    return vertex;
}

void ptd_ph_vertex_destroy(struct ptd_ph_vertex *vertex) {
    for (size_t i = 0; i < vertex->edges_length; ++i) {
        free(vertex->edges[i]);
    }

    free(vertex->edges);
    free(vertex->state);
    memset(vertex, 0, sizeof(*vertex));
    free(vertex);
}

struct ptd_ph_edge *ptd_ph_graph_add_edge(
        struct ptd_ph_vertex *from,
        struct ptd_ph_vertex *to,
        double weight
) {
    struct ptd_ph_edge *edge = (struct ptd_ph_edge *) malloc(sizeof(*edge));

    edge->to = to;
    edge->weight = weight;

    ptd_directed_graph_add_edge(
            (struct ptd_directed_vertex *) from,
            (struct ptd_directed_edge *) edge
    );

    return edge;
}
