#include <iostream>
#include <vector>
#include <cstring>
#include <stack>
#include <queue>
#include <stdint.h>
#include <stddef.h>
#include <math.h>
#include <set>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include "ptdalgorithms.h"

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
avl_node_t *rotate_left_right(avl_node_t *parent, avl_node_t *child) {
    avl_node_t *child_right_left, *child_right_right;
    avl_node_t *child_right = child->right;
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
avl_node_t *rotate_right_left(avl_node_t *parent, avl_node_t *child) {
    avl_node_t *child_left_right, *child_left_left;
    avl_node_t *child_left = child->left;

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
avl_node_t *rotate_left(avl_node_t *parent, avl_node_t *child) {
    avl_node_t *child_left;

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
avl_node_t *rotate_right(avl_node_t *parent, avl_node_t *child) {
    avl_node_t *child_right;

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

avl_node_t *avl_vec_vertex_create(char *key, void *entry, avl_node_t *parent) {
    avl_node_t *vertex;

    if ((vertex = (avl_node_t *) malloc(sizeof(*vertex))) == NULL) {
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

void avl_vec_vertex_destroy(avl_node_t *vertex) {
    if (vertex == NULL) {
        return;
    }

    avl_vec_vertex_destroy(vertex->left);
    avl_vec_vertex_destroy(vertex->right);

    free(vertex);
}

static void avl_free(avl_node_t *vertex) {
    if (vertex == NULL) {
        return;
    }

    avl_free(vertex->left);
    avl_free(vertex->right);
    free(vertex);
}

const avl_node_t *
avl_vec_find(const avl_node_t *rootptr, const char *key, const size_t vec_length) {
    if (rootptr == NULL) {
        return NULL;
    }

    const avl_node_t *vertex = rootptr;

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

int find_or_insert_vec(avl_node_t **out, avl_node_t *rootptr, char *key, void *entry,
                       const size_t vec_length) {
    if ((*out = avl_vec_vertex_create(key, entry, NULL)) == NULL) {
        return -1;
    }

    if (rootptr == NULL) {
        return 1;
    }

    avl_node_t *vertex = rootptr;

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

int avl_rebalance_tree(avl_node_t **root, avl_node_t *child) {
    avl_node_t *pivot, *rotated_parent;

    for (avl_node_t *parent = child->parent; parent != NULL; parent = child->parent) {
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

int avl_vec_insert(avl_node_t **root, char *key, void *entry, const size_t vec_length) {
    avl_node_t *child;

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


static size_t avl_vec_get_size(avl_node_t *vertex) {
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

    avl_node_t *bst = NULL;
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
        avl_node_t *bst_entry = (avl_node_t *) avl_vec_find(bst, (char *) &state[0],
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
            avl_node_t *bst_entry = (avl_node_t *) avl_vec_find(bst, (char *) &child_state[0],
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

            DEBUG_PRINT("Vertices %zu and %zu are the same group with depth %zu and nedges %zu\n",
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

static void ptd_reset_graph_visited(ptd_graph_t *graph) {
    queue<ptd_vertex_t *> to_reset;
    queue<ptd_vertex_t *> to_visit;

    to_visit.push(graph->start_vertex);

    while (!to_visit.empty()) {
        ptd_vertex_t *vertex = to_visit.front();
        to_visit.pop();

        if (!vertex->reset) {
            continue;
        }

        vertex->reset = false;
        to_reset.push(vertex);

        vertex->visited = false;

        for (size_t i = 0; i < vertex->edges_length; ++i) {
            ptd_vertex_t *child = vertex->edges[i].to;
            to_visit.push(child);
        }
    }

    while (!to_reset.empty()) {
        ptd_vertex_t *vertex = to_reset.front();
        to_reset.pop();

        vertex->reset = true;
    }
}

queue<ptd_vertex_t *> ptd_enqueue_vertices(ptd_graph_t *graph) {
    ptd_reset_graph_visited(graph);
    // TODO: add failures to these
    queue<ptd_vertex_t *> ret;
    queue<ptd_vertex_t *> queue;

    queue.push(graph->start_vertex);

    while (!queue.empty()) {
        ptd_vertex_t *vertex = queue.front();
        queue.pop();

        if (vertex->visited) {
            continue;
        }

        vertex->visited = true;
        ret.push(vertex);

        for (size_t i = 0; i < vertex->edges_length; ++i) {
            ptd_edge_t edge = vertex->edges[i];
            queue.push(edge.to);
        }
    }

    return ret;
}

int cmp_vertex_index(const void *a, const void *b) {
    if ((*((const ptd_vertex_t **) a))->index > (*((const ptd_vertex_t **) b))->index) {
        return 1;
    } else if ((*((const ptd_vertex_t **) a))->index < (*((const ptd_vertex_t **) b))->index) {
        return -1;
    } else {
        return 0;
    }
}

queue<ptd_vertex_t *> ptd_enqueue_vertices_sorted(ptd_graph_t *graph) {
    if (!graph->is_indexed) {
        ptd_label_vertices(graph);
    }

    queue<ptd_vertex_t *> unsorted = ptd_enqueue_vertices(graph);
    queue<ptd_vertex_t *> sorted_queue;
    size_t len = unsorted.size();
    ptd_vertex_t **sorted = (ptd_vertex_t **) calloc(len, sizeof(*sorted));
    size_t index = 0;

    while (!unsorted.empty()) {
        ptd_vertex_t *vertex = unsorted.front();
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

ptd_graph_t *ptd_graph_create(size_t state_length) {
    ptd_graph_t *graph;
    ptd_vertex_t *start;

    if ((graph = (ptd_graph_t *) malloc(sizeof(ptd_graph_t))) == NULL) {
        return NULL;
    }

    graph->state_length = state_length;

    if ((start = ptd_vertex_create(graph)) == NULL) {
        free(graph);
        return NULL;
    }

    graph->start_vertex = start;
    graph->vertices_length = 2;
    graph->is_indexed = true;

    return graph;
}

void ptd_graph_destroy(ptd_graph_t *graph) {
    queue<ptd_vertex_t *> q = ptd_enqueue_vertices(graph);

    // TODO: Add failures

    //if (q == NULL) {
    // We have failed to allocate or push to a queue for all the vertices.
    // This means that this freeing function cannot succeed. For now, we
    // will just ignore this error.
    //}

    while (!q.empty()) {
        ptd_vertex_t *vertex = q.front();
        q.pop();

        ptd_vertex_destroy(vertex);
    }

    memset(graph, 0, sizeof(*graph));
    free(graph);
}

queue<ptd_vertex_t *> *vertices_to_visit;

ptd_vertex_t *ptd_vertex_create_state(ptd_graph_t *graph, vec_entry_t *state) {
    ptd_vertex_t *vertex;

    vertex = (ptd_vertex_t *) malloc(sizeof(*vertex));

    if (vertex == NULL) {
        return NULL;
    }

    vertex->edges_limit = 8;
    vertex->edges = (ptd_edge_t *) calloc(8, sizeof(*vertex->edges));

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

    graph->vertices_length++;

    if (vertices_to_visit != NULL) {
        vertices_to_visit->push(vertex);
    }

    return vertex;
}

ptd_vertex_t *ptd_vertex_create(ptd_graph_t *graph) {
    vec_entry_t *state = (vec_entry_t *) calloc(graph->state_length, sizeof(*state));

    if (state == NULL) {
        return NULL;
    }

    return ptd_vertex_create_state(graph, state);
}

void ptd_vertex_destroy(ptd_vertex_t *vertex) {
    free(vertex->edges);
    vertex->edges = NULL;
    free(vertex->state);
    vertex->state = NULL;
    vertex->edges_length = 0;
    vertex->edges_limit = 0;
    free(vertex->data);
    vertex->data = NULL;

    if (vertex->graph != NULL) {
        vertex->graph->vertices_length--;
    }

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

void _ptd_avl_tree_vertex_destroy(avl_node_t *avl_vertex) {
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
    _ptd_avl_tree_vertex_destroy((avl_node_t *) avl_tree->root);
    avl_tree->root = NULL;
    free(avl_tree);
}

void _ptd_avl_tree_vertex_destroy_free(avl_node_t *avl_vertex) {
    if (avl_vertex == NULL) {
        return;
    }

    _ptd_avl_tree_vertex_destroy_free(avl_vertex->left);
    _ptd_avl_tree_vertex_destroy_free(avl_vertex->right);

    ptd_vertex_destroy((ptd_vertex_t *) avl_vertex->entry);
    avl_vertex->left = NULL;
    avl_vertex->right = NULL;
    avl_vertex->entry = NULL;
    free(avl_vertex);
}

void ptd_avl_tree_vertex_destroy_free(ptd_avl_tree_t *avl_tree) {
    _ptd_avl_tree_vertex_destroy_free((avl_node_t *) avl_tree->root);
    avl_tree->root = NULL;
    free(avl_tree);
}


void _ptd_avl_tree_edge_destroy(avl_node_t *avl_vertex) {
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
    _ptd_avl_tree_edge_destroy((avl_node_t *) avl_tree->root);
    avl_tree->root = NULL;
    free(avl_tree);
}

size_t ptd_avl_tree_max_depth(void *avl_vec_vertex) {
    if ((avl_node_t *) avl_vec_vertex == NULL) {
        return 0;
    }

    return max(
            ptd_avl_tree_max_depth((void *) ((avl_node_t *) avl_vec_vertex)->left) + 1,
            ptd_avl_tree_max_depth((void *) ((avl_node_t *) avl_vec_vertex)->left) + 1
    );
}


int ptd_avl_tree_vertex_insert(ptd_avl_tree_t *avl_tree, const vec_entry_t *key, const ptd_vertex_t *vertex) {
    int res;

    avl_node_t *root = (avl_node_t *) avl_tree->root;
    res = avl_vec_insert(&root, (char *) key, (void *) vertex, avl_tree->vec_length * sizeof(vec_entry_t));

    if (res != 0) {
        return res;
    }

    avl_tree->root = root;

    return 0;
}

ptd_vertex_t *ptd_avl_tree_vertex_find(const ptd_avl_tree_t *avl_tree, const vec_entry_t *key) {
    const avl_node_t *avl_vertex = avl_vec_find(
            (avl_node_t *) avl_tree->root,
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

    return (ptd_vertex_t *) avl_vertex->entry;
}

int ptd_avl_tree_edge_insert_or_increment(ptd_avl_tree_t *avl_tree, const vec_entry_t *key, ptd_vertex_t *vertex,
                                          long double weight) {
    avl_node_t *root = (avl_node_t *) avl_tree->root;
    avl_node_t *child;

    if (root == NULL) {
        ptd_edge_t *edge = (ptd_edge_t *) malloc(sizeof(*edge));
        edge->weight = weight;
        edge->to = vertex;

        if ((root = avl_vec_vertex_create((char *) key, (void *) edge, NULL)) == NULL) {
            return -1;
        }

        avl_tree->root = root;

        return 0;
    }

    avl_node_t *parent = root;

    while (true) {
        int res = memcmp(parent->key, key, avl_tree->vec_length * sizeof(vec_entry_t));

        if (res < 0) {
            if (parent->left == NULL) {
                ptd_edge_t *edge = (ptd_edge_t *) malloc(sizeof(*edge));
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
                ptd_edge_t *edge = (ptd_edge_t *) malloc(sizeof(*edge));
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
            ((ptd_edge_t *) parent->entry)->weight += weight;

            return 0;
        }
    }

    avl_rebalance_tree(&root, child);
    avl_tree->root = root;

    return 0;


    return 0;
}

/*avl_node_t * _ptd_avl_tree_edge_remove(avl_node_t *root, char *key, size_t length) {
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
            avl_node_t *temp = NULL;

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
    avl_node_t *root = (avl_node_t *) avl_tree->root;

    if (root == NULL) {
        return 0;
    }

    avl_node_t *parent = root;

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

        avl_rebalance_tree((avl_node_t **) &avl_tree->root, parent->parent);
    } else {
        if (parent->left != NULL) {
            avl_tree->root = parent->left;
        } else if (parent->right != NULL) {
            avl_tree->root = parent->right;
        } else {
            avl_tree->root = NULL;
            return 0;
        }

        avl_rebalance_tree((avl_node_t **) &avl_tree->root, (avl_node_t *) avl_tree->root);
    }

    return 0;
}

ptd_edge_t *ptd_avl_tree_edge_find(const ptd_avl_tree_t *avl_tree, const vec_entry_t *key) {
    avl_node_t *parent = (avl_node_t *) avl_tree->root;

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
            return (ptd_edge_t *) parent->entry;
        }
    }
}

int ptd_visit_vertices(ptd_graph_t *graph, int (*visit_func)(ptd_vertex_t *), bool include_start) {
    vertices_to_visit = new queue<ptd_vertex_t *>;
    *vertices_to_visit = ptd_enqueue_vertices_sorted(graph);

    while (!vertices_to_visit->empty()) {
        ptd_vertex_t *vertex = vertices_to_visit->front();
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

int ptd_add_edge(ptd_vertex_t *from, ptd_vertex_t *to, long double weight) {
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
        from->edges = (ptd_edge_t *) realloc(
                from->edges,
                from->edges_limit * sizeof(ptd_edge_t)
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

stack<ptd_vertex_t *> *scc_stack;
vector<ptd_strongly_connected_component_t *> *scc_components;
size_t scc_index;
size_t *scc_indices;
size_t *low_links;
bool *scc_on_stack;
ptd_strongly_connected_components_t *sccs;

int strongconnect(ptd_vertex_t *vertex, bool (*is_included_func)(ptd_vertex_t *)) {
    scc_indices[vertex->index] = scc_index;
    low_links[vertex->index] = scc_index;
    vertex->visited = true;
    scc_index++;
    scc_stack->push(vertex);
    scc_on_stack[vertex->index] = true;

    if (is_included_func(vertex)) {
        for (size_t i = 0; i < vertex->edges_length; ++i) {
            ptd_edge_t edge = vertex->edges[i];

            if (!edge.to->visited) {
                int res = strongconnect(edge.to, is_included_func);

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
    }

    if (low_links[vertex->index] == scc_indices[vertex->index]) {
        ptd_vertex_t *w;
        vector<ptd_vertex_t *> list;

        do {
            w = scc_stack->top();
            scc_stack->pop();
            scc_on_stack[w->index] = false;

            if (is_included_func(w)) {
                list.push_back(w);
            }
        } while (w != vertex);

        ptd_strongly_connected_component_t *scc = (ptd_strongly_connected_component_t *) malloc(sizeof(*scc));

        if (scc == NULL) {
            return -1;
        }

        scc->vertices_length = list.size();
        scc->vertices = (ptd_vertex_t **) calloc(scc->vertices_length, sizeof(ptd_vertex_t *));

        for (size_t i = 0; i < scc->vertices_length; ++i) {
            scc->vertices[i] = list.at(i);
        }

        scc_components->push_back(scc);
    }

    return 0;
}

int ptd_label_vertices(ptd_graph_t *graph) {
    ptd_reset_graph_visited(graph);

    queue<ptd_vertex_t *> q = ptd_enqueue_vertices(graph);
    size_t index = 1;

    while (!q.empty()) {
        ptd_vertex_t *vertex = q.front();
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
ptd_find_strongly_connected_components(ptd_graph_t *graph, bool (*is_included_func)(ptd_vertex_t *)) {
    ptd_strongly_connected_components_t *components =
            (ptd_strongly_connected_components_t *) malloc(sizeof(ptd_strongly_connected_components_t));

    if (components == NULL) {
        return NULL;
    }

    ptd_label_vertices(graph);
    scc_stack = new stack<ptd_vertex_t *>;
    queue<ptd_vertex_t *> q = ptd_enqueue_vertices(graph);
    ptd_reset_graph_visited(graph);
    scc_index = 0;
    scc_indices = (size_t *) calloc(q.size() + 2, sizeof(size_t));
    low_links = (size_t *) calloc(q.size() + 2, sizeof(size_t));
    scc_on_stack = (bool *) calloc(q.size() + 2, sizeof(bool));
    scc_components = new vector<ptd_strongly_connected_component_t *>();

    fprintf(stderr, "Variable sizes %zu\n", q.size() + 2);

    while (!q.empty()) {
        ptd_vertex_t *vertex = q.front();
        q.pop();

        if (!vertex->visited) {
            if (strongconnect(vertex, is_included_func) != 0) {
                return NULL;
            }
        }
    }

    size_t non_empty_components = 0;

    for (size_t i = 0; i < scc_components->size(); ++i) {
        ptd_strongly_connected_component_t *c = scc_components->at(i);

        if (c->vertices_length != 0) {
            non_empty_components++;
        }
    }

    components->components_length = non_empty_components;
    components->components = (ptd_strongly_connected_component_t **)
            calloc(
                    components->components_length,
                    sizeof(ptd_strongly_connected_component_t *)
            );

    size_t index = 0;

    for (size_t i = 0; i < scc_components->size(); ++i) {
        ptd_strongly_connected_component_t *c = scc_components->at(i);

        if (c->vertices_length != 0) {
            components->components[index] = (scc_components->at(i));
            index++;
        } else {
            free(c->vertices);
            free(c);
        }
    }

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
        free(sccs->components[i]->vertices);
        sccs->components[i]->vertices = NULL;

        free(sccs->components[i]);
        sccs->components[i] = NULL;
    }

    free(sccs->components);
    sccs->components = NULL;

    free(sccs);
}

ptd_scc_vertex_t **ptd_scc_index_topological(ptd_scc_vertex_t **in, size_t length) {
    ptd_scc_vertex_t **res = (ptd_scc_vertex_t **) calloc(length, sizeof(*res));
    queue<ptd_scc_vertex_t *> q;


    for (size_t i = 0; i < length; ++i) {
        ptd_scc_vertex_t *vertex = in[i];
        vertex->topo = 0;
        vertex->visited = false;
    }

    for (size_t i = 0; i < length; ++i) {
        ptd_scc_vertex_t *vertex = in[i];

        for (size_t j = 0; j < vertex->scc_edges_length; ++j) {
            ptd_scc_vertex_t *child = vertex->scc_edges[j];

            child->topo++;
        }
    }

    q.push(in[length - 1]);
    size_t idx = 0;

    while (!q.empty()) {
        ptd_scc_vertex_t *vertex = q.front();
        q.pop();

        res[idx] = vertex;
        idx++;

        for (size_t i = 0; i < vertex->scc_edges_length; ++i) {
            ptd_scc_vertex_t *child = vertex->scc_edges[i];

            child->topo -= 1;

            if (child->topo == 0 && !child->visited) {
                child->visited = true;
                q.push(child);
            }
        }
    }

    return res;
}

ptd_scc_vertex_t **
ptd_order_strongly_connected_components(
        ptd_strongly_connected_components_t *sccs
) {
    ptd_graph_t *graph = sccs->components[0]->vertices[0]->graph;

    if (!graph->is_indexed) {
        ptd_label_vertices(graph);
    }

    void **old_data = (void **) calloc(graph->vertices_length, sizeof(void *));
    ptd_scc_vertex_t **res = (ptd_scc_vertex_t **) calloc(sccs->components_length, sizeof(*res));

    for (size_t i = 0; i < sccs->components_length; ++i) {
        ptd_strongly_connected_component_t *scc = sccs->components[i];
        res[i] = (ptd_scc_vertex_t *) malloc(sizeof(*(res[i])));
        res[i]->scc = scc;
    }

    size_t p = 0;

    for (size_t i = 0; i < sccs->components_length; ++i) {
        ptd_strongly_connected_component_t *scc = sccs->components[i];

        for (size_t j = 0; j < scc->vertices_length; ++j) {
            old_data[p] = scc->vertices[j]->data;
            scc->vertices[j]->data = res[i];
            p++;
        }
    }

    for (size_t i = 0; i < sccs->components_length; ++i) {
        set<ptd_scc_vertex_t *> scc_edges;
        set<ptd_vertex_t *> local_edges;
        ptd_strongly_connected_component_t *scc = sccs->components[i];

        for (size_t j = 0; j < scc->vertices_length; ++j) {
            ptd_vertex_t *vertex = scc->vertices[j];

            for (size_t k = 0; k < vertex->edges_length; ++k) {
                ptd_vertex_t *child = vertex->edges[k].to;

                if (((ptd_scc_vertex_t *) (child->data))->scc == scc) {
                    continue;
                }

                scc_edges.insert((ptd_scc_vertex_t *) child->data);
                local_edges.insert(child);
            }
        }

        res[i]->scc_edges_length = scc_edges.size();
        res[i]->scc_edges = (ptd_scc_vertex_t **) calloc(sizeof(*(res[i]->scc_edges)), res[i]->scc_edges_length);

        size_t l = 0;
        for (
                set<ptd_scc_vertex_t *>::iterator itr = scc_edges.begin();
                itr != scc_edges.end();
                itr++
                ) {
            res[i]->scc_edges[l] = *itr;
            l++;
        }

        res[i]->local_edges_length = local_edges.size();
        res[i]->local_edges = (ptd_vertex_t **) calloc(sizeof(*(res[i]->local_edges)), res[i]->local_edges_length);

        l = 0;
        for (
                set<ptd_vertex_t *>::iterator itr = local_edges.begin();
                itr != local_edges.end();
                itr++
                ) {
            res[i]->local_edges[l] = *itr;
            l++;
        }
    }

    for (size_t i = 0; i < sccs->components_length; ++i) {
        ptd_strongly_connected_component_t *scc = sccs->components[i];

        for (size_t j = 0; j < scc->vertices_length; ++j) {
            scc->vertices[j]->data = NULL;
        }
    }

    ptd_scc_vertex_t **sorted = ptd_scc_index_topological(res, sccs->components_length);
    free(res);

    p = 0;

    for (size_t i = 0; i < sccs->components_length; ++i) {
        ptd_strongly_connected_component_t *scc = sccs->components[i];

        for (size_t j = 0; j < scc->vertices_length; ++j) {
            scc->vertices[j]->data = old_data[p];
            p++;
        }
    }

    return sorted;
}

ptd_phase_type_distribution_t *ptd_find_local_matrix(ptd_scc_vertex_t *in) {
    ptd_phase_type_distribution_t *res = (ptd_phase_type_distribution_t *) malloc(sizeof(*res));
    size_t full_length = in->local_edges_length + in->scc->vertices_length;

    long double **mat = (long double **) calloc(full_length, sizeof(long double *));
    ptd_vertex_t **vertices = (ptd_vertex_t **) calloc(full_length, sizeof(ptd_vertex_t *));
    void **old_data = (void **) calloc(full_length, sizeof(void *));

    for (size_t i = 0; i < full_length; ++i) {
        mat[i] = (long double *) calloc(full_length, sizeof(long double));
    }

    for (size_t i = 0; i < in->scc->vertices_length; ++i) {
        vertices[i] = in->scc->vertices[i];
        old_data[i] = in->scc->vertices[i]->data;
        vertices[i]->data = malloc(sizeof(size_t));
        *((size_t *) (vertices)[i]->data) = i;
    }

    for (size_t i = 0; i < in->local_edges_length; ++i) {
        size_t index = in->scc->vertices_length + i;
        (vertices)[index] = in->local_edges[i];
        (vertices)[index]->index = index;
        old_data[index] = (vertices)[index]->data;
        (vertices)[index]->data = malloc(sizeof(size_t));
        *((size_t *) (vertices)[index]->data) = index;
    }

    for (size_t i = 0; i < in->scc->vertices_length; ++i) {
        ptd_vertex_t *vertex = (vertices)[i];

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            size_t child_index = *((size_t *) vertex->edges[j].to->data);
            (mat)[i][child_index] += vertex->edges[j].weight / vertex->rate;
            (mat)[i][i] -= vertex->edges[j].weight / vertex->rate;
        }
    }

    for (size_t i = 0; i < in->local_edges_length; ++i) {
        size_t index = in->scc->vertices_length + i;
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
    }

    free(old_data);

    return res;
}

ptd_vertex_group_t *ptd_convert_strongly_connected_component_to_group(ptd_strongly_connected_component_t *scc) {
    DIE_ERROR(1, "Not implemented\n");
}

int ptd_remove_edge(ptd_vertex_t *from, ptd_vertex_t *to) {
    DIE_ERROR(1, "Not implemented\n");
}

ptd_graph_t *
ptd_convert_strongly_connected_components_to_group(ptd_graph_t *graph, ptd_strongly_connected_components_t *sccs) {
    ptd_label_vertices(graph);
    ptd_graph_t *ret = ptd_graph_create(0);

    if (graph == NULL) {
        return NULL;
    }

    ptd_vertex_t **vertices_scc;
    vertices_scc = (ptd_vertex_t **) calloc(graph->vertices_length, sizeof(*vertices_scc));

    if (vertices_scc == NULL) {
        return NULL;
    }

    for (size_t i = 0; i < sccs->components_length; ++i) {
        ptd_strongly_connected_component_t *scc = sccs->components[i];

        // Make the scc into a vertex
        ptd_vertex_t *group = ptd_vertex_create(ret);

        if (group == NULL) {
            return NULL;
        }

        for (size_t j = 0; j < scc->vertices_length; ++j) {
            ptd_vertex_t *vertex = scc->vertices[j];

            vertices_scc[vertex->index] = group;
        }
    }

    set<ptd_vertex_t *> edges;
/*
    for (size_t j = 0; j < scc->vertices_length; ++j) {
        ptd_vertex_t *vertex = scc->vertices[j];

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

double (*reward_function)(ptd_vertex_t *);

static bool keep_zero_rewarded(ptd_vertex_t *vertex) {
    return (reward_function(vertex) == 0);
}

static ptd_avl_tree_t **edges;
static ptd_avl_tree_t **parents;
static ptd_vertex_t **vertices;
static vec_entry_t **states;
static double *rewards;

static inline vector<ptd_edge_t *> avl_tree_as_list(ptd_avl_tree_t *avl_tree) {
    vector<ptd_edge_t *> vec;
    stack<avl_node_t *> s;

    s.push((avl_node_t *) avl_tree->root);

    while (!s.empty()) {
        avl_node_t *v = s.top();
        s.pop();

        if (v == NULL) {
            continue;
        }

        vec.push_back((ptd_edge_t *) v->entry);
        s.push(v->left);
        s.push(v->right);
    }

    return vec;
}

//extern int print_func(ptd_vertex_t *vertex);

// TODO: We should always reset before doing anything else everywhere

int ptd_reward_transform(ptd_graph_t *graph, double (*reward_func)(const ptd_vertex_t *)) {
    ptd_reset_graph_visited(graph);
    queue<ptd_vertex_t *> q = ptd_enqueue_vertices(graph);

    ptd_reset_graph_visited(graph);
    ptd_label_vertices(graph);
    size_t n = q.size();

    edges = (ptd_avl_tree_t **) calloc(n, sizeof(*edges));

    if (edges == NULL) {
        return 1;
    }

    parents = (ptd_avl_tree_t **) calloc(n, sizeof(*parents));
    vertices = (ptd_vertex_t **) calloc(n, sizeof(*vertices));
    states = (vec_entry_t **) calloc(n, sizeof(*states));
    rewards = (double *) calloc(n, sizeof(*rewards));

    while (!q.empty()) {
        ptd_vertex_t *vertex = q.front();
        q.pop();

        vertices[vertex->index] = vertex;
        edges[vertex->index] = ptd_avl_tree_create(1);

        if (edges[vertex->index] == NULL) {
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
        ptd_vertex_t *vertex = vertices[i];

        // We make the weight be the probability of transitioning.
        // Store the new reward for later
        rewards[i] /= vertex->rate;

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            ptd_vertex_t *child_vertex = vertex->edges[j].to;
            size_t child_index = child_vertex->index;

            ptd_avl_tree_edge_insert_or_increment(
                    edges[i],
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
            vector<ptd_edge_t *> outgoing_edges = avl_tree_as_list(edges[i]);
            vector<ptd_edge_t *> ingoing_edges = avl_tree_as_list(parents[i]);

            for (size_t p = 0; p < ingoing_edges.size(); ++p) {
                ptd_edge_t *ingoing_edge = ingoing_edges[p];
                size_t parent_index = ingoing_edge->to->index;

                for (size_t c = 0; c < outgoing_edges.size(); ++c) {
                    ptd_edge_t *outgoing_edge = outgoing_edges[c];
                    size_t child_index = outgoing_edge->to->index;

                    if (child_index != parent_index) {
                        // Add new child to parent
                        // Weight has been set to probability earlier
                        ptd_avl_tree_edge_insert_or_increment(
                                edges[parent_index],
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

                ptd_avl_tree_edge_remove(edges[parent_index], states[i]);
            }

            for (size_t c = 0; c < outgoing_edges.size(); ++c) {
                ptd_edge_t *outgoing_edge = outgoing_edges[c];
                size_t child_index = outgoing_edge->to->index;

                ptd_avl_tree_edge_remove(parents[child_index], states[i]);
            }
        }
    }

    for (size_t i = 0; i < n; ++i) {
        ptd_vertex_t *vertex = vertices[i];
        vertex->rate = 0;
        vertex->edges_length = 0;
        double reward = rewards[i];

        if (reward == 0) {
            continue;
        }

        vector<ptd_edge_t *> outgoing_edges = avl_tree_as_list(edges[i]);


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

void ptd_phase_type_distribution_free(ptd_phase_type_distribution_t *ptd) {
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

ptd_phase_type_distribution_t *ptd_graph_as_phase_type_distribution(ptd_graph_t *graph) {
    if (!graph->is_indexed) {
        ptd_label_vertices(graph);
    }

    queue<ptd_vertex_t *> q = ptd_enqueue_vertices(graph);

    ptd_phase_type_distribution_t *res = (ptd_phase_type_distribution_t *) malloc(sizeof(*res));

    if (res == NULL) {
        return NULL;
    }

    res->length = 0;

    size_t size = (size_t) q.size();

    res->memory_allocated = size;
    res->vertices = (ptd_vertex_t **) calloc(size, sizeof(ptd_vertex_t *));

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
        ptd_vertex_t *vertex = q.front();
        q.pop();

        if (vertex->index == 0) {
            continue;
        }

        if (vertex == graph->start_vertex) {
            for (size_t i = 0; i < vertex->edges_length; ++i) {
                ptd_edge_t edge = vertex->edges[i];

                if (edge.to->index != 0) {
                    res->initial_probability_vector[edge.to->index - 1 - 1] = edge.weight;
                }
            }

            continue;
        }

        res->vertices[vertex->index - 2] = vertex;
        res->length++;

        for (size_t i = 0; i < vertex->edges_length; ++i) {
            ptd_edge_t edge = vertex->edges[i];

            if (edge.to->index != 0) {
                res->sub_intensity_matrix[vertex->index - 2][edge.to->index - 2] += edge.weight;
            }

            res->sub_intensity_matrix[vertex->index - 2][vertex->index - 2] -= edge.weight;
        }
    }

    return res;
}

static int *topo_values;

int ptd_index_topological(ptd_graph_t *graph) {
    queue<ptd_vertex_t *> topological_queue;
    queue<ptd_vertex_t *> q;

    q = ptd_enqueue_vertices(graph);
    topo_values = (int *) calloc(q.size(), sizeof(*topo_values));

    if (topo_values == NULL) {
        return -1;
    }

    size_t index = 0;

    while (!q.empty()) {
        ptd_vertex_t *vertex = q.front();
        q.pop();

        vertex->index = index;

        index++;
    }

    q = ptd_enqueue_vertices(graph);

    while (!q.empty()) {
        ptd_vertex_t *vertex = q.front();
        q.pop();

        for (size_t i = 0; i < vertex->edges_length; ++i) {
            topo_values[vertex->edges[i].to->index]++;
        }
    }

    ptd_reset_graph_visited(graph);

    q.push(graph->start_vertex);

    while (!q.empty()) {
        ptd_vertex_t *vertex = q.front();
        q.pop();

        topological_queue.push(vertex);

        for (size_t i = 0; i < vertex->edges_length; ++i) {
            ptd_vertex_t *child = vertex->edges[i].to;

            topo_values[child->index] -= 1;

            if (topo_values[child->index] == 0 && !child->visited) {
                child->visited = true;
                q.push(child);
            }
        }
    }

    index = 1;

    while (!topological_queue.empty()) {
        ptd_vertex_t *vertex = topological_queue.front();
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

int ptd_index_invert(ptd_graph_t *graph) {
    queue<ptd_vertex_t *> q = ptd_enqueue_vertices(graph);
    size_t largest_index = 0;

    while (!q.empty()) {
        ptd_vertex_t *vertex = q.front();
        q.pop();

        if (vertex->index > largest_index) {
            largest_index = vertex->index;
        }
    }

    q = ptd_enqueue_vertices(graph);

    while (!q.empty()) {
        ptd_vertex_t *vertex = q.front();
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

int ptd_vertex_to_s(ptd_vertex_t *vertex, char *buffer, size_t buffer_length) {
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
static ptd_graph_t *kingman_graph = NULL;

static int make_kingman(ptd_vertex_t *vertex) {
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

            ptd_vertex_t *child = ptd_avl_tree_vertex_find(avl_tree, state);

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

ptd_graph_t *ptd_model_kingman(size_t n) {
    avl_tree = ptd_avl_tree_create(n);

    if (avl_tree == NULL) {
        return NULL;
    }

    kingman_graph = ptd_graph_create(n);

    if (kingman_graph == NULL) {
        ptd_avl_tree_vertex_destroy(avl_tree);
        return NULL;
    }

    ptd_vertex_t *initial = ptd_vertex_create(kingman_graph);

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

static ptd_graph_t *two_island_two_loci_recomb_graph;

static double recomb_rate;
static size_t N[2];
static double mig[2];

#define IDX_COMB 0
#define IDX_L 1
#define IDX_R 2

static ptd_vertex_t *add_child(vec_entry_t *state) {
    ptd_vertex_t *child = ptd_avl_tree_vertex_find(avl_tree, state);

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

static int visit_two_island_two_loci_recomb(ptd_vertex_t *vertex) {
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

                    ptd_vertex_t *child = add_child(state);

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

                ptd_vertex_t *child = add_child(state);

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

                ptd_vertex_t *child = add_child(state);

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
                ptd_vertex_t *child = add_child(state);
                state[island * 3 + IDX_COMB]++;

                ptd_add_edge(vertex, child, combinations * N[island]);
            }

            // Left and combined
            if (state[island * 3 + IDX_COMB] > 0 && state[island * 3 + IDX_L] > 0) {
                size_t combinations = state[island * 3 + IDX_COMB] * state[island * 3 + IDX_L];

                state[island * 3 + IDX_L]--;

                ptd_vertex_t *child = add_child(state);

                state[island * 3 + IDX_L]++;

                ptd_add_edge(vertex, child, combinations * N[island]);
            }

            // Right and combined
            if (state[island * 3 + IDX_COMB] > 0 && state[island * 3 + IDX_R] > 0) {
                size_t combinations = state[island * 3 + IDX_COMB] * state[island * 3 + IDX_R];

                state[island * 3 + IDX_R]--;
                ptd_vertex_t *child = add_child(state);
                state[island * 3 + IDX_R]++;

                ptd_add_edge(vertex, child, combinations * N[island]);
            }

            // Left and left
            if (state[island * 3 + IDX_L] > 1) {
                size_t combinations = state[island * 3 + IDX_L] * (state[island * 3 + IDX_L] - 1) / 2;

                state[island * 3 + IDX_L]--;
                ptd_vertex_t *child = add_child(state);
                state[island * 3 + IDX_L]++;

                ptd_add_edge(vertex, child, combinations * N[island]);
            }

            // Right and right
            if (state[island * 3 + IDX_R] > 1) {
                size_t combinations = state[island * 3 + IDX_R] * (state[island * 3 + IDX_R] - 1) / 2;

                state[island * 3 + IDX_R]--;
                ptd_vertex_t *child = add_child(state);
                state[island * 3 + IDX_R]++;

                ptd_add_edge(vertex, child, combinations * N[island]);
            }
        }
    }

    free(state);

    return 0;
}


ptd_graph_t *ptd_model_two_island_two_loci_recomb(
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

    ptd_vertex_t *initial = ptd_vertex_create(two_island_two_loci_recomb_graph);
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

static int visit_alloc(ptd_vertex_t *vertex) {
    vertex->data = malloc(alloc_size);
    memcpy((char *) vertex->data, alloc_value, alloc_size);

    return (vertex->data == NULL);
}

int ptd_visit_alloc(ptd_graph_t *graph, size_t size, char *value) {
    alloc_size = size;
    alloc_value = value;

    return ptd_visit_vertices(graph, visit_alloc, true);
}

static int visit_probability(ptd_vertex_t *vertex) {
    long double my_prob;

    if (vertex == vertex->graph->start_vertex) {
        *((long double *) vertex->data) = 1;
    }

    my_prob = *((long double *) vertex->data);

    ptd_edge_t *edges = vertex->edges;

    for (size_t i = 0; i < vertex->edges_length; ++i) {
        ptd_edge_t edge = edges[i];

        long double trans_prob = edge.weight / vertex->rate;

        *((long double *) edge.to->data) += my_prob * trans_prob;
    }

    return 0;
}

int ptd_visit_assign_probability(ptd_graph_t *graph) {
    // TODO: do only if graph is acyclic
    ptd_index_topological(graph);

    return ptd_visit_vertices(graph, visit_probability, true);
}

static int compute_expected(ptd_vertex_t *vertex) {
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

int ptd_visit_assign_expectation(ptd_graph_t *graph, double (*reward)(ptd_vertex_t *)) {
    reward_function = reward;
    return ptd_visit_vertices(graph, compute_expected, true);
}

long double reduce_ld_sum;

static int reduce_ld(ptd_vertex_t *vertex) {
    long double data = *((long double *) vertex->data);
    reduce_ld_sum += data;

    return 0;
}

long double ptd_visit_reduce_sum_long_double(ptd_graph_t *graph) {
    reduce_ld_sum = 0;
    ptd_visit_vertices(graph, reduce_ld, true);

    return reduce_ld_sum;
}

/*
 * Properties
 */
int ptd_expected_value(long double *expected, ptd_graph_t *graph, double (*reward)(ptd_vertex_t *)) {
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

static int visit_desc(ptd_vertex_t *vertex) {
    long double desc = 0;

    ptd_edge_t *edges = vertex->edges;

    for (size_t i = 0; i < vertex->edges_length; ++i) {
        ptd_edge_t edge = edges[i];

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

static int visit_cov(ptd_vertex_t *vertex) {
    struct cov data = *((struct cov *) vertex->data);
    cov_sum += data.prob * reward_function(vertex) * data.desc;

    return 0;
}

int ptd_covariance(
        long double *covariance, ptd_graph_t *graph,
        double (*reward_1)(ptd_vertex_t *), double (*reward_2)(ptd_vertex_t *)
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

static int set_data_as_int(ptd_vertex_t *vertex) {
    vertex->data = malloc(sizeof(long double));
    *((long double *) vertex->data) = 0;

    return 0;
}

static bool keep_all(ptd_vertex_t *vertex) {
    return true;
}

static gsl_matrix *
matrix_invert(gsl_matrix *matrix, size_t size) {
    gsl_matrix *inverse = gsl_matrix_alloc(size, size);

    int sign;
    gsl_permutation *p = gsl_permutation_alloc(size);

    gsl_linalg_LU_decomp(matrix, p, &sign);

    gsl_linalg_LU_invert(matrix, p, inverse);

    gsl_permutation_free(p);

    return inverse;
}

double ptd_circular_exp(ptd_graph_t *graph, double (*reward)(ptd_vertex_t *)) {
    ptd_visit_vertices(graph, set_data_as_int, true);

    *((long double *) graph->start_vertex->data) = 1;

    ptd_strongly_connected_components_t *sccs = ptd_find_strongly_connected_components(
            graph, keep_all
    );

    ptd_scc_vertex_t **vertices = ptd_order_strongly_connected_components(sccs);

    for (size_t i = 0; i < sccs->components_length; ++i) {
        ptd_scc_vertex_t *v = vertices[i];

        long double **mat;
        ptd_vertex_t **vs;
        size_t length;

        ptd_phase_type_distribution_t *ptd = ptd_find_local_matrix(v);
        mat = ptd->sub_intensity_matrix;
        vs = ptd->vertices;
        length = ptd->length;

        long double *ipv = (long double *) calloc(length, sizeof(*ipv));

        for (size_t k = 0; k < length; ++k) {
            ipv[k] = *((long double *) vs[k]->data);
            *((long double *) vs[k]->data) = 0;
        }

        ptd->initial_probability_vector = ipv;

        if (length == 1) {
            ptd_phase_type_distribution_free(ptd);
            continue;
        }
        gsl_matrix *full = gsl_matrix_alloc(length, length);

        for (size_t k = 0; k < length; ++k) {
            for (size_t j = 0; j < length; ++j) {
                gsl_matrix_set(full, k, j, (double) mat[k][j]);
            }
        }

        gsl_matrix *inv = matrix_invert(full, length);

        for (size_t k = 0; k < length; ++k) {
            for (size_t j = 0; j < length; ++j) {
                *((long double *) vs[j]->data) += ipv[k] * -gsl_matrix_get(inv, k, j);
            }
        }
        gsl_matrix_free(full);
        gsl_matrix_free(inv);

        ptd_phase_type_distribution_free(ptd);
    }

    *((long double*)graph->start_vertex->data) = 0;
    ptd_visit_assign_expectation(graph, reward);

    return ptd_visit_reduce_sum_long_double(graph);


    ptd_strongly_connected_components_destroy(sccs);
}