/*
 * MIT License
 *
 * Copyright (c) 2021 Tobias Røikjer
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <stdint.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include "ptdalgorithms.h"

volatile char ptd_err[4096] = {'\0'};

/*
 * Utility data structures
 */

struct ptd_ll {
    void *value;
    struct ptd_ll *next;
};

struct ptd_vector {
    size_t entries;
    void **arr;
};

static struct ptd_vector *vector_create();

static int vector_add(struct ptd_vector *vector, void *entry);

static void *vector_get(struct ptd_vector *vector, size_t index);

static size_t vector_length(struct ptd_vector *vector);

static void vector_destroy(struct ptd_vector *vector);

struct ptd_queue {
    struct ptd_ll *ll;
    struct ptd_ll *tail;
};

static struct ptd_queue *queue_create();

static void queue_destroy(struct ptd_queue *queue);

static int queue_enqueue(struct ptd_queue *queue, void *entry);

static void *queue_dequeue(struct ptd_queue *queue);

static int queue_empty(struct ptd_queue *queue);

struct ptd_stack {
    struct ptd_ll *ll;
};

static struct ptd_stack *stack_create();

static void stack_destroy(struct ptd_stack *stack);

static int stack_push(struct ptd_stack *stack, void *entry);

static void *stack_pop(struct ptd_stack *stack);

static int stack_empty(struct ptd_stack *stack);

struct ll_of_a {
    struct ll_of_a *next;
    double *mem;
    size_t current_mem_index;
    double *current_mem_position;
};

/*
 * AVL tree
 */

static void _ptd_avl_tree_destroy(struct ptd_avl_node *avl_vertex);

struct ptd_avl_tree *ptd_avl_tree_create(size_t key_length) {
    struct ptd_avl_tree *avl_tree = (struct ptd_avl_tree *) malloc(sizeof(struct ptd_avl_tree));

    if (avl_tree == NULL) {
        return NULL;
    }

    avl_tree->root = NULL;
    avl_tree->key_length = key_length;

    return avl_tree;
}

void ptd_avl_tree_destroy(struct ptd_avl_tree *avl_tree) {
    _ptd_avl_tree_destroy((struct ptd_avl_node *) avl_tree->root);
    avl_tree->root = NULL;
    free(avl_tree);
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
struct ptd_avl_node *rotate_left_right(struct ptd_avl_node *parent, struct ptd_avl_node *child) {
    struct ptd_avl_node *child_right_left, *child_right_right;
    struct ptd_avl_node *child_right = child->right;
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
struct ptd_avl_node *rotate_right_left(struct ptd_avl_node *parent, struct ptd_avl_node *child) {
    struct ptd_avl_node *child_left_right, *child_left_left;
    struct ptd_avl_node *child_left = child->left;

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
struct ptd_avl_node *rotate_left(struct ptd_avl_node *parent, struct ptd_avl_node *child) {
    struct ptd_avl_node *child_left;

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
struct ptd_avl_node *rotate_right(struct ptd_avl_node *parent, struct ptd_avl_node *child) {
    struct ptd_avl_node *child_right;

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

struct ptd_avl_node *ptd_avl_node_create(const int *key, const void *entry, struct ptd_avl_node *parent) {
    struct ptd_avl_node *vertex;

    if ((vertex = (struct ptd_avl_node *) malloc(sizeof(*vertex))) == NULL) {
        return NULL;
    }

    vertex->key = (int *) key;
    vertex->entry = (void *) entry;
    vertex->left = NULL;
    vertex->right = NULL;
    vertex->parent = parent;
    vertex->balance = 0;

    return vertex;
}

static void ptd_avl_node_destroy(struct ptd_avl_node *vertex) {
    if (vertex == NULL) {
        return;
    }

    ptd_avl_node_destroy(vertex->left);
    ptd_avl_node_destroy(vertex->right);

    free(vertex);
}

static void avl_free(struct ptd_avl_node *vertex) {
    if (vertex == NULL) {
        return;
    }

    avl_free(vertex->left);
    avl_free(vertex->right);
    free(vertex);
}

const struct ptd_avl_node *
avl_vec_find(const struct ptd_avl_node *rootptr, const char *key, const size_t vec_length) {
    if (rootptr == NULL) {
        return NULL;
    }

    const struct ptd_avl_node *vertex = rootptr;

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

int find_or_insert_vec(struct ptd_avl_node **out, struct ptd_avl_node *rootptr, int *key, void *entry,
                       const size_t vec_length) {
    if ((*out = ptd_avl_node_create(key, entry, NULL)) == NULL) {
        return -1;
    }

    if (rootptr == NULL) {
        return 1;
    }

    struct ptd_avl_node *vertex = rootptr;

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
            free(*out);
            *out = vertex;
            return 0;
        }
    }

    (*out)->parent = vertex;

    return 0;
}

int avl_rebalance_tree(struct ptd_avl_node **root, struct ptd_avl_node *child) {
    struct ptd_avl_node *pivot, *rotated_parent;

    for (struct ptd_avl_node *parent = child->parent; parent != NULL; parent = child->parent) {
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


static size_t avl_vec_get_size(struct ptd_avl_node *vertex) {
    if (vertex == NULL) {
        return 0;
    }

    return 1 + avl_vec_get_size(vertex->left) + avl_vec_get_size(vertex->right);
}


static void _ptd_avl_tree_destroy(struct ptd_avl_node *avl_vertex) {
    if (avl_vertex == NULL) {
        return;
    }

    _ptd_avl_tree_destroy(avl_vertex->left);
    _ptd_avl_tree_destroy(avl_vertex->right);

    avl_vertex->left = NULL;
    avl_vertex->right = NULL;
    avl_vertex->entry = NULL;
    free(avl_vertex);
}

#define _ptd_max(a, b) a >= b ? a : b
#define _ptd_min(a, b) a <= b ? a : b

size_t ptd_avl_tree_max_depth(void *avl_vec_vertex) {
    if ((struct ptd_avl_node *) avl_vec_vertex == NULL) {
        return 0;
    }

    return _ptd_max(
                   ptd_avl_tree_max_depth((void *) ((struct ptd_avl_node *) avl_vec_vertex)->left) + 1,
                   ptd_avl_tree_max_depth((void *) ((struct ptd_avl_node *) avl_vec_vertex)->left) + 1
           );
}


struct ptd_avl_node *ptd_avl_tree_find_or_insert(struct ptd_avl_tree *avl_tree, const int *key, const void *entry) {
    struct ptd_avl_node *new_node = ptd_avl_node_create(key, entry, NULL);

    if (new_node == NULL) {
        return NULL;
    }

    if (avl_tree->root == NULL) {
        avl_tree->root = new_node;

        return new_node;
    }

    struct ptd_avl_node *vertex = avl_tree->root;

    while (true) {
        int res = memcmp(key, vertex->key, sizeof(int) * avl_tree->key_length);

        if (res < 0) {
            if (vertex->left == NULL) {
                vertex->left = new_node;
                break;
            } else {
                vertex = vertex->left;
            }
        } else if (res > 0) {
            if (vertex->right == NULL) {
                vertex->right = new_node;
                break;
            } else {
                vertex = vertex->right;
            }
        } else {
            free(new_node);
            return vertex;
        }
    }

    new_node->parent = vertex;

    avl_rebalance_tree(&avl_tree->root, new_node);

    return new_node;
}

struct ptd_avl_node *ptd_avl_tree_find(const struct ptd_avl_tree *avl_tree, const int *key) {
    struct ptd_avl_node *vertex = avl_tree->root;

    while (true) {
        if (vertex == NULL) {
            return NULL;
        }

        int res = memcmp(key, vertex->key, sizeof(int) * avl_tree->key_length);

        if (res < 0) {
            vertex = vertex->left;
        } else if (res > 0) {
            vertex = vertex->right;
        } else {
            return vertex;
        }
    }
}

int ptd_precompute_reward_compute_graph(struct ptd_graph *graph) {
    if (graph->was_dph) {
        graph->was_dph = false;

        if (graph->reward_compute_graph != NULL) {
            free(graph->reward_compute_graph->commands);
            free(graph->reward_compute_graph);
        }

        if (graph->parameterized_reward_compute_graph != NULL) {
            ptd_parameterized_reward_compute_graph_destroy(
                    graph->parameterized_reward_compute_graph
            );
        }

        graph->reward_compute_graph = NULL;
        graph->parameterized_reward_compute_graph = NULL;
    }

    if (graph->reward_compute_graph == NULL) {
        if (graph->parameterized) {
            if (graph->parameterized_reward_compute_graph == NULL) {
                DEBUG_PRINT("INFO: building parameterized compute graph...\n");
                graph->parameterized_reward_compute_graph =
                        ptd_graph_ex_absorbation_time_comp_graph_parameterized(graph);
            }

            if (graph->reward_compute_graph != NULL) {
                free(graph->reward_compute_graph->commands);
                free(graph->reward_compute_graph);
            }

            DEBUG_PRINT("INFO: building reward compute graph from parameterized compute graph...\n");
            graph->reward_compute_graph =
                    ptd_graph_build_ex_absorbation_time_comp_graph_parameterized(
                            graph->parameterized_reward_compute_graph
                    );
        } else {
            DEBUG_PRINT("INFO: building reward compute graph...\n");
            graph->reward_compute_graph = ptd_graph_ex_absorbation_time_comp_graph(graph);

            if (graph->reward_compute_graph == NULL) {
                return -1;
            }
        }
    }

    return 0;
}


static struct ptd_stack *scc_stack2 = NULL;
static struct ptd_vector *scc_components2 = NULL;
static size_t scc_index2 = 0;
static size_t *scc_indices2 = NULL;
static size_t *low_links2 = NULL;
static bool *scc_on_stack2 = NULL;
static bool *visited = NULL;

static int strongconnect2(struct ptd_vertex *vertex) {
    scc_indices2[vertex->index] = scc_index2;
    low_links2[vertex->index] = scc_index2;
    visited[vertex->index] = true;
    scc_index2++;
    stack_push(scc_stack2, vertex);
    scc_on_stack2[vertex->index] = true;

    for (size_t i = 0; i < vertex->edges_length; ++i) {
        struct ptd_edge *edge = vertex->edges[i];

        if (!visited[edge->to->index]) {
            int res = strongconnect2(edge->to);

            if (res != 0) {
                return res;
            }

            low_links2[vertex->index] = _ptd_min(
                                                low_links2[vertex->index],
                                                low_links2[edge->to->index]
                                        );
        } else if (scc_on_stack2[edge->to->index]) {
            low_links2[vertex->index] = _ptd_min(
                                                low_links2[vertex->index],
                                                scc_indices2[edge->to->index]
                                        );
        }
    }

    if (low_links2[vertex->index] == scc_indices2[vertex->index]) {
        struct ptd_vertex *w;
        struct ptd_vector *list = vector_create();

        do {
            if (stack_empty(scc_stack2)) {
                DIE_ERROR(1, "Stack is empty.\n");
            }
            w = (struct ptd_vertex *) stack_pop(scc_stack2);
            scc_on_stack2[w->index] = false;

            vector_add(list, w);
        } while (w != vertex);

        struct ptd_scc_vertex *scc = (struct ptd_scc_vertex *) malloc(sizeof(*scc));

        if (scc == NULL) {
            return -1;
        }

        scc->internal_vertices_length = vector_length(list);
        scc->internal_vertices = (struct ptd_vertex **) calloc(
                scc->internal_vertices_length,
                sizeof(*(scc->internal_vertices))
        );

        for (size_t i = 0; i < scc->internal_vertices_length; ++i) {
            scc->internal_vertices[i] = (struct ptd_vertex *) vector_get(list, i);
        }

        vector_add(scc_components2, scc);
        vector_destroy(list);
    }

    return 0;
}

static int scc_edge_cmp(const void *a, const void *b) {
    struct ptd_scc_edge *ea = *((struct ptd_scc_edge **) a);
    struct ptd_scc_edge *eb = *((struct ptd_scc_edge **) b);

    if (ea->to->index < eb->to->index) {
        return -1;
    } else if (ea->to->index > eb->to->index) {
        return 1;
    } else {
        return 0;
    }
}

struct ptd_scc_graph *ptd_find_strongly_connected_components(struct ptd_graph *graph) {
    struct ptd_scc_graph *scc_graph = (struct ptd_scc_graph *) malloc(
            sizeof(*scc_graph)
    );

    scc_graph->graph = graph;

    scc_stack2 = NULL;
    scc_components2 = NULL;
    scc_index2 = 0;
    scc_indices2 = NULL;
    low_links2 = NULL;
    scc_on_stack2 = NULL;
    visited = NULL;

    scc_stack2 = stack_create();

    scc_index2 = 0;
    scc_indices2 = (size_t *) calloc(graph->vertices_length * 10, sizeof(size_t));
    low_links2 = (size_t *) calloc(graph->vertices_length * 10, sizeof(size_t));
    scc_on_stack2 = (bool *) calloc(graph->vertices_length * 10, sizeof(bool));
    visited = (bool *) calloc(graph->vertices_length * 10, sizeof(bool));
    scc_components2 = vector_create();

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        struct ptd_vertex *vertex = graph->vertices[i];

        if (!visited[i]) {
            if (strongconnect2(vertex) != 0) {
                return NULL;
            }
        }
    }

    size_t non_empty_components = 0;

    for (size_t i = 0; i < vector_length(scc_components2); ++i) {
        struct ptd_scc_vertex *c =
                (struct ptd_scc_vertex *) vector_get(scc_components2, i);

        if (c->internal_vertices_length != 0) {
            non_empty_components++;
        }
    }

    scc_graph->vertices_length = non_empty_components;
    scc_graph->vertices = (struct ptd_scc_vertex **) calloc(
            scc_graph->vertices_length,
            sizeof(*(scc_graph->vertices))
    );

    size_t index = 0;

    for (size_t i = 0; i < scc_graph->vertices_length; ++i) {
        struct ptd_scc_vertex *scc =
                (struct ptd_scc_vertex *) vector_get(scc_components2, i);

        if (scc->internal_vertices_length != 0) {
            scc_graph->vertices[index] =
                    (struct ptd_scc_vertex *) vector_get(scc_components2, i);
            scc_graph->vertices[index]->index = index;
            index++;
        } else {
            free(scc->internal_vertices);
            free(scc);
        }
    }

    struct ptd_scc_vertex **sccs_for_vertices = (struct ptd_scc_vertex **) calloc(
            graph->vertices_length,
            sizeof(*sccs_for_vertices)
    );

    for (size_t i = 0; i < scc_graph->vertices_length; ++i) {
        struct ptd_scc_vertex *scc = scc_graph->vertices[i];
        scc->index = i;

        for (size_t j = 0; j < scc->internal_vertices_length; ++j) {
            struct ptd_vertex *vertex = scc->internal_vertices[j];

            sccs_for_vertices[vertex->index] = scc;
        }
    }

    scc_graph->starting_vertex = sccs_for_vertices[graph->starting_vertex->index];

    for (size_t i = 0; i < scc_graph->vertices_length; ++i) {
        struct ptd_scc_vertex *scc = scc_graph->vertices[i];
        struct ptd_avl_tree *external_sccs = ptd_avl_tree_create(1);

        for (size_t j = 0; j < scc->internal_vertices_length; ++j) {
            struct ptd_vertex *vertex = scc->internal_vertices[j];

            for (size_t k = 0; k < vertex->edges_length; ++k) {
                struct ptd_vertex *child = vertex->edges[k]->to;
                struct ptd_scc_vertex *child_scc = sccs_for_vertices[child->index];

                if (child_scc != scc) {
                    ptd_avl_tree_find_or_insert(external_sccs, (int *) &(child_scc->index), child_scc);
                }
            }
        }

        struct ptd_vector *external_sccs_vector = vector_create();
        struct ptd_stack *tree_stack;
        tree_stack = stack_create();

        if (external_sccs->root != NULL) {
            stack_push(tree_stack, external_sccs->root);
        }

        while (!stack_empty(tree_stack)) {
            struct ptd_avl_node *node = (struct ptd_avl_node *) stack_pop(tree_stack);
            vector_add(external_sccs_vector, node->entry);

            if (node->left != NULL) {
                stack_push(tree_stack, node->left);
            }

            if (node->right != NULL) {
                stack_push(tree_stack, node->right);
            }
        }

        scc->edges_length = vector_length(external_sccs_vector);
        scc->edges = (struct ptd_scc_edge **) calloc(
                scc->edges_length,
                sizeof(*(scc->edges))
        );

        size_t set_index;

        set_index = 0;

        for (size_t l = 0; l < vector_length(external_sccs_vector); ++l) {
            scc->edges[set_index] = (struct ptd_scc_edge *) malloc(sizeof(*(scc->edges[set_index])));
            scc->edges[set_index]->to = (struct ptd_scc_vertex *) vector_get(external_sccs_vector, l);
            set_index++;
        }

        qsort(scc->edges, scc->edges_length, sizeof(*(scc->edges)), scc_edge_cmp);

        vector_destroy(external_sccs_vector);
        stack_destroy(tree_stack);
        ptd_avl_tree_destroy(external_sccs);
    }

    free(scc_indices2);
    free(low_links2);
    free(scc_on_stack2);
    free(visited);
    vector_destroy(scc_components2);
    stack_destroy(scc_stack2);

    free(sccs_for_vertices);


    scc_stack2 = NULL;
    scc_components2 = NULL;
    scc_index2 = 0;
    scc_indices2 = NULL;
    low_links2 = NULL;
    scc_on_stack2 = NULL;
    visited = NULL;

    return scc_graph;
}

void ptd_scc_graph_destroy(struct ptd_scc_graph *scc_graph) {
    if (scc_graph == NULL) {
        return;
    }

    for (size_t i = 0; i < scc_graph->vertices_length; ++i) {
        struct ptd_scc_vertex *scc = scc_graph->vertices[i];

        for (size_t j = 0; j < scc->edges_length; ++j) {
            free(scc->edges[j]);
        }


        free(scc->edges);
        free(scc->internal_vertices);
        free(scc);
    }

    free(scc_graph->vertices);
    free(scc_graph);
}

double *ptd_normalize_graph(struct ptd_graph *graph) {
    double *res = (double *) calloc(graph->vertices_length, sizeof(*res));

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        struct ptd_vertex *vertex = graph->vertices[i];
        double rate = 0;

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            rate += vertex->edges[j]->weight;
        }

        if (rate == 0) {
            res[i] = 1.0;
        } else {
            res[i] = 1.0 / rate;
        }

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            vertex->edges[j]->weight /= rate;
        }
    }

    return res;
}

double *ptd_dph_normalize_graph(struct ptd_graph *graph) {
    size_t old_length = graph->vertices_length;
    double *res = (double *) calloc(old_length * 2, sizeof(*res));

    for (size_t i = 0; i < old_length; ++i) {
        res[i] = 1;

        struct ptd_vertex *vertex = graph->vertices[i];
        double rate = 0;

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            rate += vertex->edges[j]->weight;
        }

        if (rate == 0 || graph->starting_vertex == vertex) {
            continue;
        }

        if (1 - rate > 0.0000001) {
            struct ptd_vertex *auxiliary_vertex = ptd_vertex_create(graph);
            ptd_graph_add_edge(vertex, auxiliary_vertex, 1 - rate);
            ptd_graph_add_edge(auxiliary_vertex, vertex, 1);
            res[auxiliary_vertex->index] = 0;
        }
    }

    return res;
}

struct ptd_phase_type_distribution *ptd_graph_as_phase_type_distribution(struct ptd_graph *graph) {
    struct ptd_phase_type_distribution *res = (struct ptd_phase_type_distribution *) malloc(sizeof(*res));

    if (res == NULL) {
        return NULL;
    }

    res->length = 0;

    size_t size = graph->vertices_length;

    res->memory_allocated = size;
    res->vertices = (struct ptd_vertex **) calloc(size, sizeof(struct ptd_vertex *));

    if (res->vertices == NULL) {
        free(res);
        return NULL;
    }

    res->initial_probability_vector = (double *) calloc(size, sizeof(double));

    if (res->initial_probability_vector == NULL) {
        free(res->vertices);
        free(res);
        return NULL;
    }

    res->sub_intensity_matrix = (double **) calloc(size, sizeof(double *));

    if (res->sub_intensity_matrix == NULL) {
        free(res->initial_probability_vector);
        free(res->vertices);
        free(res);
        return NULL;
    }

    for (size_t i = 0; i < size; ++i) {
        res->sub_intensity_matrix[i] = (double *) calloc(size, sizeof(double));

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

    struct ptd_scc_graph *scc = ptd_find_strongly_connected_components(graph);
    struct ptd_vertex **vertices =
            (struct ptd_vertex **) calloc(graph->vertices_length, sizeof(*vertices));
    struct ptd_scc_vertex **v = ptd_scc_graph_topological_sort(scc);
    size_t idx = 0;

    for (size_t i = 0; i < scc->vertices_length; ++i) {
        for (size_t j = 0; j < v[i]->internal_vertices_length; ++j) {
            vertices[idx] = v[i]->internal_vertices[j];
            vertices[idx]->index = idx;
            idx++;
        }
    }

    size_t *indices = (size_t *) calloc(size, sizeof(*indices));
    size_t index = 0;

    for (size_t k = 0; k < graph->vertices_length; ++k) {
        struct ptd_vertex *vertex = vertices[k];

        if (graph->starting_vertex != vertex && vertex->edges_length != 0) {
            indices[vertex->index] = index;
            res->vertices[index] = vertex;
            index++;
        }
    }

    res->length = index;

    for (size_t k = 0; k < graph->vertices_length; ++k) {
        struct ptd_vertex *vertex = vertices[k];

        if (vertex->edges_length == 0) {
            continue;
        }

        if (vertex == graph->starting_vertex) {
            double rate = 0;

            for (size_t i = 0; i < vertex->edges_length; ++i) {
                struct ptd_edge *edge = vertex->edges[i];

                rate += edge->weight;
            }

            for (size_t i = 0; i < vertex->edges_length; ++i) {
                struct ptd_edge *edge = vertex->edges[i];

                if (edge->to->edges_length != 0) {
                    res->initial_probability_vector[indices[edge->to->index]] = edge->weight / rate;
                }
            }

            continue;
        }

        for (size_t i = 0; i < vertex->edges_length; ++i) {
            struct ptd_edge *edge = vertex->edges[i];

            if (edge->to->edges_length != 0) {
                res->sub_intensity_matrix[indices[vertex->index]][indices[edge->to->index]] += edge->weight;
            }

            res->sub_intensity_matrix[indices[vertex->index]][indices[vertex->index]] -= edge->weight;
        }
    }

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        graph->vertices[i]->index = i;
    }

    free(v);
    ptd_scc_graph_destroy(scc);
    free(indices);
    free(vertices);

    return res;
}

void ptd_phase_type_distribution_destroy(struct ptd_phase_type_distribution *ptd) {
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

int ptd_vertex_to_s(struct ptd_vertex *vertex, char *buffer, size_t buffer_length) {
    memset(buffer, '\0', buffer_length);

    char *build = (char *) calloc(buffer_length, sizeof(char));

    for (size_t i = 0; i < vertex->graph->state_length; ++i) {
        if (i == 0) {
            snprintf(build, buffer_length, "%s%i", buffer, vertex->state[i]);
        } else {
            snprintf(build, buffer_length, "%s %i", buffer, vertex->state[i]);
        }

        strncpy(buffer, build, buffer_length);
    }

    free(build);

    return 0;
}

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

struct ptd_graph *ptd_graph_create(size_t state_length) {
    struct ptd_graph *graph = (struct ptd_graph *) malloc(sizeof(*graph));
    graph->vertices_length = 0;
    graph->state_length = state_length;
    graph->vertices = NULL;
    graph->reward_compute_graph = NULL;
    graph->parameterized_reward_compute_graph = NULL;
    graph->starting_vertex = ptd_vertex_create(graph);
    graph->parameterized = false;
    graph->was_dph = false;

    return graph;
}

void ptd_parameterized_reward_compute_graph_destroy(
        struct ptd_desc_reward_compute_parameterized *compute_graph
) {
    struct ll_of_a *mem = (struct ll_of_a *) compute_graph->mem;

    while (mem != NULL) {
        struct ll_of_a *memp = mem;
        mem = mem->next;
        free(memp->mem);
        free(memp);
    }

    free(compute_graph->memr);
    free(compute_graph->commands);
    free(compute_graph);
}

void ptd_graph_destroy(struct ptd_graph *graph) {
    for (size_t i = 0; i < graph->vertices_length; ++i) {
        ptd_vertex_destroy(graph->vertices[i]);
    }

    free(graph->vertices);

    if (graph->reward_compute_graph != NULL) {
        free(graph->reward_compute_graph->commands);
        free(graph->reward_compute_graph);
    }

    if (graph->parameterized_reward_compute_graph != NULL) {
        ptd_parameterized_reward_compute_graph_destroy(
                graph->parameterized_reward_compute_graph
        );
    }

    graph->reward_compute_graph = NULL;
    graph->parameterized_reward_compute_graph = NULL;
    memset(graph, 0, sizeof(*graph));
    free(graph);
}

struct ptd_vertex *ptd_vertex_create(struct ptd_graph *graph) {
    int *state = (int *) calloc(graph->state_length, sizeof(*state));

    return ptd_vertex_create_state(graph, state);
}

struct ptd_vertex *ptd_vertex_create_state(struct ptd_graph *graph, int *state) {
    struct ptd_vertex *vertex = (struct ptd_vertex *) malloc(sizeof(*vertex));
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

double ptd_vertex_rate(struct ptd_vertex *vertex) {
    double rate = 0;

    for (size_t i = 0; i < vertex->edges_length; ++i) {
        rate += vertex->edges[i]->weight;
    }

    return rate;
}

void ptd_vertex_destroy(struct ptd_vertex *vertex) {
    for (size_t i = 0; i < vertex->edges_length; ++i) {
        if (vertex->edges[i]->parameterized) {
            if (((struct ptd_edge_parameterized *) vertex->edges[i])->should_free_state) {
                free(((struct ptd_edge_parameterized *) vertex->edges[i])->state);
            }
        }

        free(vertex->edges[i]);
    }

    free(vertex->edges);
    free(vertex->state);
    memset(vertex, 0, sizeof(*vertex));
    free(vertex);
}


static inline int edge_cmp(const void *a, const void *b) {
    if ((*((struct ptd_edge **) a))->to < (*((struct ptd_edge **) b))->to) {
        return -1;
    } else {
        return 1;
    }
}

int ptd_validate_graph(const struct ptd_graph *graph) {
    struct ptd_edge **edges_buffer = (struct ptd_edge **) calloc(graph->vertices_length, sizeof(*edges_buffer));

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        struct ptd_vertex *vertex = graph->vertices[i];

        if (vertex->edges_length >= graph->vertices_length) {
            // Definitely have a problem...
            edges_buffer = (struct ptd_edge **) realloc(edges_buffer, vertex->edges_length * sizeof(*edges_buffer));
        }

        memcpy(edges_buffer, vertex->edges, vertex->edges_length * sizeof(*edges_buffer));
        qsort(edges_buffer, vertex->edges_length, sizeof(*edges_buffer), edge_cmp);

        for (size_t j = 1; j < vertex->edges_length; ++j) {
            if (vertex->edges[j]->to == vertex->edges[j - 1]->to) {
                struct ptd_vertex *from = vertex;
                struct ptd_vertex *to = vertex->edges[j]->to;
                size_t debug_index_from = from->index;
                size_t debug_index_to = to->index;

                if (PTD_DEBUG_1_INDEX) {
                    debug_index_from++;
                    debug_index_to++;
                }

                char state[1024] = {'\0'};
                char state_to[1024] = {'\0'};
                char starting_vertex[] = " (starting vertex)";

                if (from != from->graph->starting_vertex) {
                    starting_vertex[0] = '\0';
                }

                ptd_vertex_to_s(from, state, 1023);
                ptd_vertex_to_s(to, state_to, 1023);

                snprintf(
                        (char *) ptd_err,
                        sizeof(ptd_err),
                        "Multiple edges to the same vertex!. From vertex with index %i%s (state %s)."
                        " To vertex with index %i (state %s)\n",
                        (int) debug_index_from, starting_vertex, state,
                        (int) debug_index_to, state_to
                );

                free(edges_buffer);
                return 1;
            }
        }
    }

    free(edges_buffer);

    return 0;
}

struct ptd_edge *ptd_graph_add_edge(
        struct ptd_vertex *from,
        struct ptd_vertex *to,
        double weight
) {
    if (weight < 0) {
        size_t debug_index = from->index;

        if (PTD_DEBUG_1_INDEX) {
            debug_index++;
        }

        char state[1024] = {'\0'};
        char starting_vertex[] = " (starting vertex)";

        if (from != from->graph->starting_vertex) {
            starting_vertex[0] = '\0';
        }

        ptd_vertex_to_s(from, state, 1023);

        snprintf(
                (char *) ptd_err,
                sizeof(ptd_err),
                "Tried to add edge with non-positive weight '%f'. Vertex index %i%s (state %s). Weight must be strictly larger than 0.\n",
                weight, (int) debug_index, starting_vertex, state
        );

        return NULL;
    }

    if (from == to) {
        size_t debug_index = from->index;

        if (PTD_DEBUG_1_INDEX) {
            debug_index++;
        }

        char state[1024] = {'\0'};
        char starting_vertex[] = " (starting vertex)";

        if (from != from->graph->starting_vertex) {
            starting_vertex[0] = '\0';
        }

        ptd_vertex_to_s(from, state, 1023);

        snprintf(
                (char *) ptd_err,
                sizeof(ptd_err),
                "Tried to add edge to itself. Vertex index %i%s (state %s). Self-loops are not allowed, discrete self-loops are set as the missing out-going weight.\n",
                (int) debug_index, starting_vertex, state
        );

        return NULL;
    }

    /*for (size_t i = 0; i < from->edges_length; ++i) {
        if (from->edges[i]->to == to) {
            size_t debug_index = from->index;

            if (PTD_DEBUG_1_INDEX) {
                debug_index++;
            }
            size_t debug_index_to = to->index;

            if (PTD_DEBUG_1_INDEX) {
                debug_index_to++;
            }

            char state[1024] = {'\0'};
            char starting_vertex[] = " (starting vertex)";

            if (from != from->graph->starting_vertex) {
                starting_vertex[0] = '\0';
            }

            ptd_vertex_to_s(from, state, 1023);

            char state_to[1024] = {'\0'};

            ptd_vertex_to_s(to, state_to, 1023);

            snprintf(
                    (char *) ptd_err,
                    sizeof(ptd_err),
                    "Tried to add to a vertex with an already existing edge. Vertex index %i%s (state %s), to %i (state %s).\n",
                    (int) debug_index, starting_vertex, state, (int) debug_index_to, state_to
            );

            return NULL;
        }
    }*/

    struct ptd_edge *edge = (struct ptd_edge *) malloc(sizeof(*edge));

    edge->to = to;
    edge->weight = weight;
    edge->parameterized = false;

    ptd_directed_graph_add_edge(
            (struct ptd_directed_vertex *) from,
            (struct ptd_directed_edge *) edge
    );

    if (from->graph->reward_compute_graph != NULL) {
        free(from->graph->reward_compute_graph->commands);
        free(from->graph->reward_compute_graph);
    }

    if (from->graph->parameterized_reward_compute_graph != NULL) {
        ptd_parameterized_reward_compute_graph_destroy(
                from->graph->parameterized_reward_compute_graph
        );
    }

    from->graph->reward_compute_graph = NULL;
    from->graph->parameterized_reward_compute_graph = NULL;

    return edge;
}

struct ptd_edge_parameterized *ptd_graph_add_edge_parameterized(
        struct ptd_vertex *from,
        struct ptd_vertex *to,
        double weight,
        double *edge_state
) {
    //from->graph->parameterized = true;

    struct ptd_edge_parameterized *edge = (struct ptd_edge_parameterized *) malloc(sizeof(*edge));

    edge->to = to;
    edge->weight = weight;
    edge->parameterized = true;
    edge->state = edge_state;
    edge->should_free_state = true;

    ptd_directed_graph_add_edge(
            (struct ptd_directed_vertex *) from,
            (struct ptd_directed_edge *) edge
    );

    if (from->graph->reward_compute_graph != NULL) {
        free(from->graph->reward_compute_graph->commands);
        free(from->graph->reward_compute_graph);
    }

    if (from->graph->parameterized_reward_compute_graph != NULL) {
        ptd_parameterized_reward_compute_graph_destroy(
                from->graph->parameterized_reward_compute_graph
        );
    }

    from->graph->reward_compute_graph = NULL;
    from->graph->parameterized_reward_compute_graph = NULL;

    return edge;
}

void ptd_notify_change(
        struct ptd_graph *graph
) {
    if (graph->reward_compute_graph != NULL) {
        free(graph->reward_compute_graph->commands);
        free(graph->reward_compute_graph);
        graph->reward_compute_graph = NULL;
    }
}

void ptd_edge_update_weight(
        struct ptd_edge *edge,
        double weight
) {
    edge->weight = weight;

    if (edge->to->graph->reward_compute_graph != NULL) {
        free(edge->to->graph->reward_compute_graph->commands);
        edge->to->graph->reward_compute_graph = NULL;
    }
}

void ptd_edge_update_weight_parameterized(
        struct ptd_edge *edge,
        double *scalars,
        size_t scalars_length
) {
    double weight = 0;

    for (size_t i = 0; i < scalars_length; ++i) {
        weight += scalars[i] * ((struct ptd_edge_parameterized *) edge)->state[i];
    }

    edge->weight = weight;

    if (edge->to->graph->reward_compute_graph != NULL) {
        free(edge->to->graph->reward_compute_graph->commands);
        edge->to->graph->reward_compute_graph = NULL;
    }
}

void ptd_graph_update_weight_parameterized(
        struct ptd_graph *graph,
        double *scalars,
        size_t scalars_length
) {
    for (size_t i = 0; i < graph->vertices_length; ++i) {
        for (size_t j = 0; j < graph->vertices[i]->edges_length; ++j) {
            if (graph->vertices[i]->edges[j]->parameterized) {
                ptd_edge_update_weight_parameterized(
                        graph->vertices[i]->edges[j], scalars, scalars_length
                );
            }
        }
    }
}


struct ptd_directed_vertex **ptd_directed_graph_topological_sort(struct ptd_directed_graph *graph) {
    struct ptd_directed_vertex **res = (struct ptd_directed_vertex **) calloc(
            graph->vertices_length, sizeof(*res)
    );

    bool *visited = (bool *) calloc(graph->vertices_length, sizeof(*visited));
    size_t *nparents = (size_t *) calloc(graph->vertices_length, sizeof(*nparents));

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        struct ptd_directed_vertex *vertex = graph->vertices[i];

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            struct ptd_directed_vertex *child = vertex->edges[j]->to;

            nparents[child->index]++;
        }
    }

    bool has_pushed_all_others = false;
    struct ptd_queue *q = queue_create();
    queue_enqueue(q, graph->vertices[0]);
    size_t topo_index = 0;

    while (!queue_empty(q)) {
        struct ptd_directed_vertex *vertex = (struct ptd_directed_vertex *) queue_dequeue(q);

        res[topo_index] = vertex;
        visited[vertex->index] = true;
        topo_index++;

        for (size_t i = 0; i < vertex->edges_length; ++i) {
            struct ptd_directed_vertex *child = vertex->edges[i]->to;

            nparents[child->index]--;

            if (nparents[child->index] == 0 && !visited[child->index]) {
                visited[child->index] = true;
                queue_enqueue(q, child);
            }
        }

        if (queue_empty(q) && !has_pushed_all_others) {
            for (size_t i = 0; i < graph->vertices_length; ++i) {
                struct ptd_directed_vertex *independent_vertex = graph->vertices[i];

                if (nparents[independent_vertex->index] == 0 && !visited[independent_vertex->index]) {
                    queue_enqueue(q, independent_vertex);
                }
            }

            has_pushed_all_others = true;
        }
    }

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        if (nparents[i] != 0) {
            free(nparents);
            free(visited);
            free(res);
            queue_destroy(q);

            return NULL;
        }
    }

    free(nparents);
    free(visited);
    queue_destroy(q);

    return res;
}

bool ptd_graph_is_acyclic(struct ptd_graph *graph) {
    struct ptd_vertex **sorted = ptd_graph_topological_sort(graph);

    bool is_acyclic = (sorted != NULL);

    free(sorted);

    return is_acyclic;
}

struct ptd_vertex **ptd_graph_topological_sort(struct ptd_graph *graph) {
    return (struct ptd_vertex **) ptd_directed_graph_topological_sort((struct ptd_directed_graph *) graph);
}

struct ptd_scc_vertex **ptd_scc_graph_topological_sort(struct ptd_scc_graph *graph) {
    return (struct ptd_scc_vertex **) ptd_directed_graph_topological_sort((struct ptd_directed_graph *) graph);
}

struct ll_c {
    struct ll_c *next;
    struct ll_c *prev;
    double weight;
    struct ptd_vertex *c;
    struct ll_p *ll_p;
};

struct ll_p {
    struct ll_p *next;
    struct ll_p *prev;
    struct ptd_vertex *p;
    struct ll_c *ll_c;
};


struct ll_c2 {
    struct ll_c2 *next;
    struct ll_c2 *prev;
    double *weight;
    struct ptd_vertex *c;
    struct ll_p2 *ll_p;
};

struct ll_p2 {
    struct ll_p2 *next;
    struct ll_p2 *prev;
    struct ptd_vertex *p;
    struct ll_c2 *ll_c;
};

#define REWARD_EPSILON 0.000001

struct arr_c {
    double prob;
    struct ptd_vertex *to;
    size_t arr_p_index;
};

struct arr_p {
    struct ptd_vertex *p;
    size_t arr_c_index;
};

static inline int arr_c_cmp(const void *a, const void *b) {
    if ((*((struct arr_c *) a)).to < (*((struct arr_c *) b)).to) {
        return -1;
    } else {
        return 1;
    }
}

struct ptd_graph *_ptd_graph_reward_transform(struct ptd_graph *graph, double *__rewards, size_t **new_indices_r) {
    double *rewards = (double *) calloc(graph->vertices_length, sizeof(*rewards));

    struct ptd_vertex *dummy__ptd_min = (struct ptd_vertex *) 1, *dummy__ptd_max = 0;

    struct ptd_vertex **vertices = (struct ptd_vertex **) calloc(graph->vertices_length, sizeof(*vertices));
    size_t *original_indices = (size_t *) calloc(graph->vertices_length, sizeof(*original_indices));

    size_t vertices_length = graph->vertices_length;

    struct ptd_scc_graph *scc = ptd_find_strongly_connected_components(graph);
    struct ptd_scc_vertex **v = ptd_scc_graph_topological_sort(scc);

    size_t idx = 0;

    for (size_t sii = 0; sii < scc->vertices_length; ++sii) {
        for (size_t j = 0; j < v[sii]->internal_vertices_length; ++j) {
            if (v[sii]->internal_vertices[j]->edges_length == 0) {
                continue;
            }

            original_indices[idx] = v[sii]->internal_vertices[j]->index;
            v[sii]->internal_vertices[j]->index = idx;
            vertices[idx] = v[sii]->internal_vertices[j];
            idx++;
        }
    }

    for (size_t sii = 0; sii < scc->vertices_length; ++sii) {
        for (size_t j = 0; j < v[sii]->internal_vertices_length; ++j) {
            if (v[sii]->internal_vertices[j]->edges_length != 0) {
                continue;
            }

            original_indices[idx] = v[sii]->internal_vertices[j]->index;
            v[sii]->internal_vertices[j]->index = idx;
            vertices[idx] = v[sii]->internal_vertices[j];
            idx++;
        }
    }

    for (size_t i = 0; i < vertices_length; ++i) {
        if (__rewards[original_indices[i]] <= REWARD_EPSILON) {
            rewards[i] = 0;
        } else {
            rewards[i] = __rewards[original_indices[i]];
        }

        if (graph->starting_vertex == vertices[i] || vertices[i]->edges_length == 0) {
            rewards[i] = 1;
        }
    }

    struct arr_p **vertex_parents;
    size_t *vertex_parents_length;
    struct arr_c **vertex_edges;
    size_t *vertex_edges_length;
    double *old_rates = (double *) calloc(vertices_length, sizeof(*old_rates));

    for (size_t i = 0; i < vertices_length; ++i) {
        struct ptd_vertex *vertex = vertices[i];

        if (vertex >= dummy__ptd_max) {
            dummy__ptd_max = vertex + 1;
        }

        if (vertex <= dummy__ptd_min) {
            dummy__ptd_min = vertex - 1;
        }

        double rate = 0;

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            rate += vertex->edges[j]->weight;
        }

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            vertex->edges[j]->weight /= rate;
        }

        if (rewards[i] != 0) {
            rewards[i] /= rate;
        }

        old_rates[i] = rate;
    }

    vertex_parents = (struct arr_p **) calloc(vertices_length, sizeof(*vertex_parents));
    vertex_parents_length = (size_t *) calloc(vertices_length, sizeof(*vertex_parents_length));
    size_t *vertex_parents_alloc_length = (size_t *) calloc(vertices_length, sizeof(*vertex_parents_alloc_length));
    vertex_edges = (struct arr_c **) calloc(vertices_length, sizeof(*vertex_edges));
    vertex_edges_length = (size_t *) calloc(vertices_length, sizeof(*vertex_edges_length));
    size_t *vertex_edges_alloc_length = (size_t *) calloc(vertices_length, sizeof(*vertex_edges_alloc_length));

    for (size_t i = 0; i < vertices_length; ++i) {
        vertex_edges_alloc_length[i] = 64;
        struct ptd_vertex *vertex = vertices[i];

        while (vertex->edges_length + 2 >= vertex_edges_alloc_length[i]) {
            vertex_edges_alloc_length[i] *= 2;
        }

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            vertex_parents_length[vertex->edges[j]->to->index]++;
        }

        vertex_edges[i] = (struct arr_c *) calloc(vertex_edges_alloc_length[i], sizeof(*(vertex_edges[i])));
        vertex_edges_length[i] = vertex->edges_length + 2;
    }

    for (size_t i = 0; i < vertices_length; ++i) {
        vertex_parents_alloc_length[i] = 64;

        while (vertex_parents_length[i] >= vertex_parents_alloc_length[i]) {
            vertex_parents_alloc_length[i] *= 2;
        }

        vertex_parents[i] = (struct arr_p *) calloc(vertex_parents_alloc_length[i], sizeof(*(vertex_parents[i])));
        vertex_parents_length[i] = 0;
    }

    for (size_t i = 0; i < vertices_length; ++i) {
        struct ptd_vertex *vertex = vertices[i];

        vertex_edges[i][0].to = dummy__ptd_min;
        vertex_edges[i][0].prob = 0;
        vertex_edges[i][0].arr_p_index = (unsigned int) ((int) -1);

        double rate = 0;

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            rate += vertex->edges[j]->weight;
        }

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            vertex_edges[i][j + 1].to = vertex->edges[j]->to;
            vertex_edges[i][j + 1].prob = vertex->edges[j]->weight / rate;
        }

        vertex_edges[i][vertex->edges_length + 1].prob = 0;
        vertex_edges[i][vertex->edges_length + 1].to = dummy__ptd_max;
        vertex_edges[i][vertex->edges_length + 1].arr_p_index = (unsigned int) ((int) -1);

        qsort(vertex_edges[i], vertex_edges_length[i], sizeof(*(vertex_edges[i])), arr_c_cmp);
    }


    for (size_t i = 0; i < vertices_length; ++i) {
        struct ptd_vertex *vertex = vertices[i];

        for (size_t j = 1; j < vertex_edges_length[i] - 1; ++j) {
            struct arr_c *child = &(vertex_edges[i][j]);
            size_t k = child->to->index;
            child->arr_p_index = vertex_parents_length[k];
            vertex_parents[k][vertex_parents_length[k]].p = vertex;
            vertex_parents[k][vertex_parents_length[k]].arr_c_index = j;
            vertex_parents_length[k]++;
        }
    }

    struct arr_c *old_edges_buffer =
            (struct arr_c *) calloc(vertices_length + 2, sizeof(*old_edges_buffer));

    for (size_t i = 0; i < vertices_length; ++i) {
        if (rewards[i] != 0) {
            continue;
        }

        struct ptd_vertex *me = vertices[i];
        struct arr_c *my_children = vertex_edges[i];
        size_t my_parents_length = vertex_parents_length[i];
        size_t my_edges_length = vertex_edges_length[i];


        for (size_t p = 0; p < my_parents_length; ++p) {
            struct arr_p me_to_parent = vertex_parents[i][p];
            struct ptd_vertex *parent_vertex = me_to_parent.p;

            size_t parent_vertex_index = parent_vertex->index;
            struct arr_c parent_to_me = vertex_edges[parent_vertex_index][me_to_parent.arr_c_index];

            size_t parent_edges_length = vertex_edges_length[parent_vertex_index];

            bool should_resize = false;
            size_t new_parent_edges_alloc_length = my_edges_length + parent_edges_length;

            while (new_parent_edges_alloc_length >= vertex_edges_alloc_length[parent_vertex_index]) {
                vertex_edges_alloc_length[parent_vertex_index] *= 2;
                should_resize = true;
            }

            if (should_resize) {
                vertex_edges[parent_vertex_index] = (struct arr_c *) realloc(
                        vertex_edges[parent_vertex_index],
                        vertex_edges_alloc_length[parent_vertex_index] * sizeof(*(vertex_edges[parent_vertex_index]))
                );
            }

            vertex_edges_length[parent_vertex_index] = 0;

            double parent_weight_to_me = parent_to_me.prob;
            double new_parent_total_prob = 0;

            memcpy(
                    old_edges_buffer, vertex_edges[parent_vertex_index],
                    sizeof(struct arr_c) * parent_edges_length
            );

            struct arr_c *new_parent_children = vertex_edges[parent_vertex_index];

            size_t child_index = 0;
            size_t parent_child_index = 0;

            while (child_index < my_edges_length || parent_child_index < parent_edges_length) {
                struct arr_c me_to_child = my_children[child_index];
                struct ptd_vertex *me_to_child_v = me_to_child.to;
                struct arr_c parent_to_child = old_edges_buffer[parent_child_index];
                struct ptd_vertex *parent_to_child_v = parent_to_child.to;
                double me_to_child_p = me_to_child.prob;

                if (me_to_child_v == parent_vertex) {
                    double prob = parent_weight_to_me * me_to_child_p;
                    rewards[parent_vertex_index] *= 1 / (1 - prob);

                    child_index++;
                    continue;
                }

                if (parent_to_child_v == me) {
                    parent_child_index++;
                    continue;
                }

                if (me_to_child_v == parent_to_child_v) {
                    new_parent_children[vertex_edges_length[parent_vertex_index]].to = parent_to_child_v;
                    new_parent_children[vertex_edges_length[parent_vertex_index]].prob =
                            parent_to_child.prob + me_to_child_p * parent_weight_to_me;

                    new_parent_children[vertex_edges_length[parent_vertex_index]].arr_p_index = parent_to_child.arr_p_index;

                    if (parent_to_child_v != dummy__ptd_min && parent_to_child_v != dummy__ptd_max) {
                        size_t current_parent_index = parent_to_child.arr_p_index;
                        vertex_parents[parent_to_child_v->index][current_parent_index].arr_c_index = vertex_edges_length[parent_vertex_index];

                    }

                    new_parent_total_prob += new_parent_children[vertex_edges_length[parent_vertex_index]].prob;
                    vertex_edges_length[parent_vertex_index]++;

                    child_index++;
                    parent_child_index++;
                } else if (me_to_child_v < parent_to_child_v) {
                    size_t child_parents_length = vertex_parents_length[me_to_child_v->index];

                    if (child_parents_length >= vertex_parents_alloc_length[me_to_child_v->index]) {
                        vertex_parents_alloc_length[me_to_child_v->index] *= 2;
                        vertex_parents[me_to_child_v->index] = (struct arr_p *) realloc(
                                vertex_parents[me_to_child_v->index],
                                vertex_parents_alloc_length[me_to_child_v->index] *
                                sizeof(*(vertex_parents[me_to_child_v->index]))
                        );
                    }

                    vertex_parents[me_to_child_v->index][child_parents_length].arr_c_index = vertex_edges_length[parent_vertex_index];
                    vertex_parents[me_to_child_v->index][child_parents_length].p = parent_vertex;

                    new_parent_children[vertex_edges_length[parent_vertex_index]].to = me_to_child_v;
                    new_parent_children[vertex_edges_length[parent_vertex_index]].prob =
                            me_to_child_p * parent_weight_to_me;
                    new_parent_children[vertex_edges_length[parent_vertex_index]].arr_p_index = child_parents_length;
                    new_parent_total_prob += me_to_child_p * parent_weight_to_me;

                    vertex_edges_length[parent_vertex_index]++;
                    vertex_parents_length[me_to_child_v->index]++;

                    child_index++;
                } else {
                    new_parent_children[vertex_edges_length[parent_vertex_index]] = parent_to_child;
                    vertex_parents[parent_to_child_v->index][parent_to_child.arr_p_index].arr_c_index = vertex_edges_length[parent_vertex_index];
                    new_parent_total_prob += parent_to_child.prob;
                    vertex_edges_length[parent_vertex_index]++;

                    parent_child_index++;
                }
            }


            // Make sure parent has rate of 1
            for (size_t j = 0; j < vertex_edges_length[parent_vertex_index]; ++j) {
                new_parent_children[j].prob /= new_parent_total_prob;
            }

            vertex_edges_length[parent_vertex_index] = vertex_edges_length[parent_vertex_index];
        }

        for (size_t j = 1; j < my_edges_length - 1; ++j) {
            struct arr_c me_to_child = my_children[j];
            struct ptd_vertex *me_to_child_v = me_to_child.to;
            size_t index_to_remove = me_to_child.arr_p_index;
            size_t index_to_move = vertex_parents_length[me_to_child_v->index] - 1;
            vertex_parents[me_to_child_v->index][index_to_remove] =
                    vertex_parents[me_to_child_v->index][index_to_move];
            vertex_parents_length[me_to_child_v->index]--;
            struct arr_p child_to_move_parent = vertex_parents[me_to_child_v->index][index_to_remove];
            vertex_edges[child_to_move_parent.p->index][child_to_move_parent.arr_c_index].arr_p_index = index_to_remove;
        }
    }

    struct ptd_graph *new_graph = ptd_graph_create(graph->state_length);
    size_t *new_indicesGtoN = (size_t *) calloc(vertices_length, sizeof(*new_indicesGtoN));
    size_t *new_indicesNtoG = (size_t *) calloc(vertices_length, sizeof(*new_indicesNtoG));
    size_t *new_indicesNtoO = (size_t *) calloc(vertices_length, sizeof(*new_indicesNtoO));
    new_indicesGtoN[graph->starting_vertex->index] = 0;
    new_indicesNtoG[0] = graph->starting_vertex->index;
    new_indicesNtoO[0] = 0;
    size_t new_idx = 1;
    memcpy(graph->starting_vertex->state, new_graph->starting_vertex->state, graph->state_length * sizeof(int));

    for (size_t i = 0; i < vertices_length; ++i) {
        if (vertices[i] == graph->starting_vertex) {
            continue;
        }

        if (rewards[i] == 0) {
            continue;
        }

        struct ptd_vertex *vertex = ptd_vertex_create(new_graph);
        memcpy(vertex->state, vertices[i]->state, graph->state_length * sizeof(int));
        new_indicesGtoN[i] = new_idx;
        new_indicesNtoG[new_idx] = i;
        new_indicesNtoO[new_idx] = original_indices[i];
        new_idx++;
    }

    for (size_t i = 0; i < vertices_length; ++i) {
        if (rewards[i] == 0) {
            continue;
        }

        for (size_t j = 1; j < vertex_edges_length[i] - 1; ++j) {
            ptd_graph_add_edge(
                    new_graph->vertices[new_indicesGtoN[i]],
                    new_graph->vertices[new_indicesGtoN[vertex_edges[i][j].to->index]],
                    vertex_edges[i][j].prob / rewards[i]
            );
        }
    }

    *(new_indices_r) = new_indicesNtoO;

    free(new_indicesGtoN);
    free(new_indicesNtoG);

    for (size_t i = 0; i < vertices_length; ++i) {
        struct ptd_vertex *vertex = vertices[i];

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            vertex->edges[j]->weight *= old_rates[i];
        }
    }

    for (size_t i = 0; i < vertices_length; ++i) {
        graph->vertices[i]->index = i;
    }

    for (size_t i = 0; i < vertices_length; ++i) {
        free(vertex_edges[i]);
        free(vertex_parents[i]);
    }

    free(old_rates);
    free(vertex_parents_length);
    free(vertex_parents_alloc_length);
    free(vertex_parents);
    free(vertex_edges);
    free(vertex_edges_length);
    free(vertex_edges_alloc_length);
    free(original_indices);
    free(vertices);
    free(old_edges_buffer);
    free(v);
    ptd_scc_graph_destroy(scc);
    free(rewards);


    return new_graph;
}

struct ptd_graph *ptd_graph_reward_transform(struct ptd_graph *graph, double *rewards) {
    if (ptd_validate_graph(graph)) {
        return NULL;
    }

    size_t *new_indices;
    struct ptd_graph *res = _ptd_graph_reward_transform(graph, rewards, &new_indices);

    free(new_indices);

    return res;
}

struct ptd_graph *ptd_graph_dph_reward_transform(struct ptd_graph *_graph, int *rewards) {
    if (ptd_validate_graph(_graph)) {
        return NULL;
    }

    for (size_t i = 0; i < _graph->vertices_length; ++i) {
        if (rewards[i] <= REWARD_EPSILON) {
            continue;
        }

        struct ptd_vertex *vertex = _graph->vertices[i];

        if (vertex->edges_length == 0) {
            continue;
        }
        if (vertex == _graph->starting_vertex) {
            continue;
        }

        double rate = 0;

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            rate += vertex->edges[j]->weight;
        }

        if (rate > 1.0001) {
            size_t debug_index = vertex->index;

            if (PTD_DEBUG_1_INDEX) {
                debug_index++;
            }

            char state[1024] = {'\0'};
            char starting_vertex[] = " (starting vertex)";

            if (vertex != _graph->starting_vertex) {
                starting_vertex[0] = '\0';
            }

            ptd_vertex_to_s(vertex, state, 1023);

            snprintf(
                    (char *) ptd_err,
                    sizeof(ptd_err),
                    "Expected vertex with index %i%s (state %s) to have outgoing rate <= 1. Is '%f'. Are you sure this is a discrete phase-type distribution?\n",
                    (int) debug_index, starting_vertex, state, (float) rate
            );

            return NULL;
        }
    }

    double *zero_rewards = (double *) calloc(_graph->vertices_length, sizeof(*zero_rewards));

    for (size_t i = 0; i < _graph->vertices_length; ++i) {
        if (rewards[i] == 0) {
            zero_rewards[i] = 0;
        } else {
            zero_rewards[i] = 1;
        }
    }

    zero_rewards[0] = 1;

    size_t *new_graph_indices;
    struct ptd_graph *graph = _ptd_graph_reward_transform(_graph, zero_rewards, &new_graph_indices);

    struct ptd_vertex **vertices = (struct ptd_vertex **) calloc(
            graph->vertices_length, sizeof(*vertices)
    );

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        vertices[i] = graph->vertices[i];
    }

    free(zero_rewards);

    int *non_zero_rewards = (int *) calloc(
            graph->vertices_length, sizeof(*non_zero_rewards)
    );

    for (size_t i = 1; i < graph->vertices_length; ++i) {
        size_t old_index = new_graph_indices[i];

        non_zero_rewards[i] = rewards[old_index];
    }

    free(vertices);

    size_t old_length = graph->vertices_length;

    for (size_t i = 0; i < old_length; ++i) {
        struct ptd_vertex *vertex = graph->vertices[i];

        if (vertex->edges_length == 0) {
            continue;
        }

        if (non_zero_rewards[i] == 1) {
            continue;
        }

        if (vertex == graph->starting_vertex) {
            continue;
        }

        double rate = 0;

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            rate += vertex->edges[j]->weight;
        }

        if (rate > 1.0001) {
            size_t debug_index = vertex->index;

            if (PTD_DEBUG_1_INDEX) {
                debug_index++;
            }

            char state[1024] = {'\0'};
            char starting_vertex[] = " (starting vertex)";

            if (vertex != graph->starting_vertex) {
                starting_vertex[0] = '\0';
            }

            ptd_vertex_to_s(vertex, state, 1023);

            snprintf(
                    (char *) ptd_err,
                    sizeof(ptd_err),
                    "Expected vertex with index %i%s (state %s) to have outgoing rate <= 1. Is '%f'. Are you sure this is a discrete phase-type distribution?\n",
                    (int) debug_index, starting_vertex, state, (float) rate
            );

            free(non_zero_rewards);

            return NULL;
        }

        struct ptd_vertex **auxiliary_vertices = (struct ptd_vertex **) calloc(
                (size_t) non_zero_rewards[i],
                sizeof(*auxiliary_vertices)
        );

        auxiliary_vertices[0] = vertex;

        for (int k = 1; k < non_zero_rewards[i]; ++k) {
            auxiliary_vertices[k] = ptd_vertex_create(graph);
        }

        size_t edges_length = vertex->edges_length;

        for (size_t j = 0; j < edges_length; ++j) {
            ptd_graph_add_edge(
                    auxiliary_vertices[non_zero_rewards[i] - 1],
                    vertex->edges[j]->to,
                    vertex->edges[j]->weight
            );

            free(vertex->edges[j]);
        }

        vertex->edges_length = 0;

        for (int k = 0; k < non_zero_rewards[i] - 1; ++k) {
            ptd_graph_add_edge(
                    auxiliary_vertices[k],
                    auxiliary_vertices[k + 1],
                    1
            );
        }

        if (1 - rate > REWARD_EPSILON) {
            ptd_graph_add_edge(
                    auxiliary_vertices[non_zero_rewards[i] - 1],
                    vertex,
                    1 - rate
            );
        }

        free(auxiliary_vertices);
    }

    free(new_graph_indices);
    free(non_zero_rewards);

    return graph;
}


static struct ptd_reward_increase *add_command(
        struct ptd_reward_increase *cmd,
        size_t from,
        size_t to,
        double weight,
        size_t index
) {
    bool is_power_of_2 = (index & (index - 1)) == 0;

    if (is_power_of_2) {
        size_t new_length = index == 0 ? 1 : index * 2;

        cmd = (struct ptd_reward_increase *) realloc(
                cmd, new_length *
                     sizeof(*cmd)
        );
    }

    if (from != to) {
//        fprintf(stderr, "ADD COMMAND %zu += %zu * %f\n", from, to, weight);
        cmd[index].from = from;
        cmd[index].to = to;
        cmd[index].multiplier = weight;
    } else {
        //      fprintf(stderr, "ADD COMMAND %zu *= %f\n", from, weight);
        cmd[index].from = from;
        cmd[index].to = to;
        cmd[index].multiplier = weight - 1;
    }

    return cmd;
}

enum command_types {
    PP = 3,
    P = 1,
    INV = 2,
    ZERO = 6,
    DIVIDE = 5,
    ONE_MINUS = 4,
    NEW_ADD = 0
};

static struct ptd_comp_graph_parameterized *add_command_param_pp(
        struct ptd_comp_graph_parameterized *cmd,
        double *from,
        double *to,
        double *weight,
        size_t index
) {
    bool is_power_of_2 = (index & (index - 1)) == 0;

    if (is_power_of_2) {
        size_t new_length = index == 0 ? 1 : index * 2;

        cmd = (struct ptd_comp_graph_parameterized *) realloc(
                cmd, new_length *
                     sizeof(*cmd)
        );
    }

    cmd[index].type = PP;

    if (from != to) {
        cmd[index].fromT = from;
        cmd[index].toT = to;
        cmd[index].multiplierptr = weight;
    } else {
        cmd[index].fromT = from;
        cmd[index].toT = to;
        cmd[index].multiplierptr = weight - 1;
    }

    return cmd;
}


static struct ptd_comp_graph_parameterized *add_command_param_p(
        struct ptd_comp_graph_parameterized *cmd,
        double *from,
        double *to,
        double weight,
        size_t index
) {
    bool is_power_of_2 = (index & (index - 1)) == 0;

    if (is_power_of_2) {
        size_t new_length = index == 0 ? 1 : index * 2;

        cmd = (struct ptd_comp_graph_parameterized *) realloc(
                cmd, new_length *
                     sizeof(*cmd)
        );
    }

    cmd[index].type = P;

    if (from != to) {
        cmd[index].fromT = from;
        cmd[index].toT = to;
        cmd[index].multiplier = weight;
    } else {
        cmd[index].fromT = from;
        cmd[index].toT = to;
        cmd[index].multiplier = weight - 1;
    }

    return cmd;
}

static struct ptd_comp_graph_parameterized *add_command_param_inverse(
        struct ptd_comp_graph_parameterized *cmd,
        double *from,
        size_t index
) {
    bool is_power_of_2 = (index & (index - 1)) == 0;

    if (is_power_of_2) {
        size_t new_length = index == 0 ? 1 : index * 2;

        cmd = (struct ptd_comp_graph_parameterized *) realloc(
                cmd, new_length *
                     sizeof(*cmd)
        );
    }

    cmd[index].type = INV;
    cmd[index].fromT = from;

    return cmd;
}

static struct ptd_comp_graph_parameterized *add_command_param_zero(
        struct ptd_comp_graph_parameterized *cmd,
        double *from,
        size_t index
) {
    bool is_power_of_2 = (index & (index - 1)) == 0;

    if (is_power_of_2) {
        size_t new_length = index == 0 ? 1 : index * 2;

        cmd = (struct ptd_comp_graph_parameterized *) realloc(
                cmd, new_length *
                     sizeof(*cmd)
        );
    }

    cmd[index].type = ZERO;
    cmd[index].fromT = from;

    return cmd;
}


static struct ptd_comp_graph_parameterized *add_command_param_p_divide(
        struct ptd_comp_graph_parameterized *cmd,
        double *from,
        double *to,
        size_t index
) {
    bool is_power_of_2 = (index & (index - 1)) == 0;

    if (is_power_of_2) {
        size_t new_length = index == 0 ? 1 : index * 2;

        cmd = (struct ptd_comp_graph_parameterized *) realloc(
                cmd, new_length *
                     sizeof(*cmd)
        );
    }

    cmd[index].type = DIVIDE;
    cmd[index].fromT = from;
    cmd[index].toT = to;

    return cmd;
}


static struct ptd_comp_graph_parameterized *add_command_param_one__ptd_minus(
        struct ptd_comp_graph_parameterized *cmd,
        double *from,
        size_t index
) {
    bool is_power_of_2 = (index & (index - 1)) == 0;

    if (is_power_of_2) {
        size_t new_length = index == 0 ? 1 : index * 2;

        cmd = (struct ptd_comp_graph_parameterized *) realloc(
                cmd, new_length *
                     sizeof(*cmd)
        );
    }

    cmd[index].type = ONE_MINUS;
    cmd[index].fromT = from;

    return cmd;
}


static struct ptd_comp_graph_parameterized *add_command_param(
        struct ptd_comp_graph_parameterized *cmd,
        size_t from,
        size_t to,
        double *weight,
        size_t index
) {
    bool is_power_of_2 = (index & (index - 1)) == 0;

    if (is_power_of_2) {
        size_t new_length = index == 0 ? 1 : index * 2;

        cmd = (struct ptd_comp_graph_parameterized *) realloc(
                cmd, new_length *
                     sizeof(*cmd)
        );
    }

    cmd[index].type = NEW_ADD;

    cmd[index].from = from;
    cmd[index].to = to;
    cmd[index].multiplierptr = weight;

    return cmd;
}

struct ptd_desc_reward_compute *ptd_graph_ex_absorbation_time_comp_graph(struct ptd_graph *graph) {
    if (ptd_validate_graph(graph)) {
        return NULL;
    }

    struct ptd_vertex *dummy__ptd_min = (struct ptd_vertex *) 1, *dummy__ptd_max = 0;

    struct ptd_vertex **vertices = (struct ptd_vertex **) calloc(graph->vertices_length, sizeof(*vertices));
    size_t *original_indices = (size_t *) calloc(graph->vertices_length, sizeof(*original_indices));

    struct ptd_reward_increase *commands = NULL;
    size_t command_index = 0;
    size_t vertices_length = graph->vertices_length;

    struct ptd_scc_graph *scc = ptd_find_strongly_connected_components(graph);
    struct ptd_scc_vertex **v = ptd_scc_graph_topological_sort(scc);

    size_t idx = 0;

    for (size_t sii = 0; sii < scc->vertices_length; ++sii) {
        for (size_t j = 0; j < v[sii]->internal_vertices_length; ++j) {
            if (v[sii]->internal_vertices[j]->edges_length == 0) {
                continue;
            }

            original_indices[idx] = v[sii]->internal_vertices[j]->index;
            v[sii]->internal_vertices[j]->index = idx;
            vertices[idx] = v[sii]->internal_vertices[j];
            idx++;
        }
    }

    for (size_t sii = 0; sii < scc->vertices_length; ++sii) {
        for (size_t j = 0; j < v[sii]->internal_vertices_length; ++j) {
            if (v[sii]->internal_vertices[j]->edges_length != 0) {
                continue;
            }

            original_indices[idx] = v[sii]->internal_vertices[j]->index;
            v[sii]->internal_vertices[j]->index = idx;
            vertices[idx] = v[sii]->internal_vertices[j];
            idx++;
        }
    }

    struct arr_p **vertex_parents;
    size_t *vertex_parents_length;
    struct arr_c **vertex_edges;
    size_t *vertex_edges_length;

    for (size_t i = 0; i < vertices_length; ++i) {
        struct ptd_vertex *vertex = vertices[i];

        if (vertex >= dummy__ptd_max) {
            dummy__ptd_max = vertex + 1;
        }

        if (vertex <= dummy__ptd_min) {
            dummy__ptd_min = vertex - 1;
        }

        double rate = 0;

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            rate += vertex->edges[j]->weight;
        }

        // Add the "real" rate as our first reward

        if (graph->starting_vertex == vertex || vertex->edges_length == 0) {
            commands = add_command(
                    commands,
                    original_indices[i],
                    original_indices[i],
                    0,
                    command_index++
            );
        } else {
            commands = add_command(
                    commands,
                    original_indices[i],
                    original_indices[i],
                    1 / rate,
                    command_index++
            );
        }
    }

    vertex_parents = (struct arr_p **) calloc(vertices_length, sizeof(*vertex_parents));
    vertex_parents_length = (size_t *) calloc(vertices_length, sizeof(*vertex_parents_length));
    size_t *vertex_parents_alloc_length = (size_t *) calloc(vertices_length, sizeof(*vertex_parents_alloc_length));
    vertex_edges = (struct arr_c **) calloc(vertices_length, sizeof(*vertex_edges));
    vertex_edges_length = (size_t *) calloc(vertices_length, sizeof(*vertex_edges_length));
    size_t *vertex_edges_alloc_length = (size_t *) calloc(vertices_length, sizeof(*vertex_edges_alloc_length));

    for (size_t i = 0; i < vertices_length; ++i) {
        vertex_edges_alloc_length[i] = 64;
        struct ptd_vertex *vertex = vertices[i];

        while (vertex->edges_length + 2 >= vertex_edges_alloc_length[i]) {
            vertex_edges_alloc_length[i] *= 2;
        }

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            vertex_parents_length[vertex->edges[j]->to->index]++;
        }

        vertex_edges[i] = (struct arr_c *) calloc(vertex_edges_alloc_length[i], sizeof(*(vertex_edges[i])));
        vertex_edges_length[i] = vertex->edges_length + 2;
    }

    for (size_t i = 0; i < vertices_length; ++i) {
        vertex_parents_alloc_length[i] = 64;

        while (vertex_parents_length[i] >= vertex_parents_alloc_length[i]) {
            vertex_parents_alloc_length[i] *= 2;
        }

        vertex_parents[i] = (struct arr_p *) calloc(vertex_parents_alloc_length[i], sizeof(*(vertex_parents[i])));
        vertex_parents_length[i] = 0;
    }

    for (size_t i = 0; i < vertices_length; ++i) {
        struct ptd_vertex *vertex = vertices[i];

        vertex_edges[i][0].to = dummy__ptd_min;
        vertex_edges[i][0].prob = 0;
        vertex_edges[i][0].arr_p_index = (unsigned int) ((int) -1);

        double rate = 0;

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            rate += vertex->edges[j]->weight;
        }

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            vertex_edges[i][j + 1].to = vertex->edges[j]->to;
            vertex_edges[i][j + 1].prob = vertex->edges[j]->weight / rate;
        }

        vertex_edges[i][vertex->edges_length + 1].prob = 0;
        vertex_edges[i][vertex->edges_length + 1].to = dummy__ptd_max;
        vertex_edges[i][vertex->edges_length + 1].arr_p_index = (unsigned int) ((int) -1);

        qsort(vertex_edges[i], vertex_edges_length[i], sizeof(*(vertex_edges[i])), arr_c_cmp);
    }


    for (size_t i = 0; i < vertices_length; ++i) {
        struct ptd_vertex *vertex = vertices[i];

        for (size_t j = 1; j < vertex_edges_length[i] - 1; ++j) {
            struct arr_c *child = &(vertex_edges[i][j]);
            size_t k = child->to->index;
            child->arr_p_index = vertex_parents_length[k];
            vertex_parents[k][vertex_parents_length[k]].p = vertex;
            vertex_parents[k][vertex_parents_length[k]].arr_c_index = j;
            vertex_parents_length[k]++;
        }
    }

    struct arr_c *old_edges_buffer =
            (struct arr_c *) calloc(vertices_length + 2, sizeof(*old_edges_buffer));

    for (size_t i = 0; i < vertices_length; ++i) {
        struct ptd_vertex *me = vertices[i];
        struct arr_c *my_children = vertex_edges[i];
        size_t my_parents_length = vertex_parents_length[i];
        size_t my_edges_length = vertex_edges_length[i];


        for (size_t p = 0; p < my_parents_length; ++p) {
            struct arr_p me_to_parent = vertex_parents[i][p];
            struct ptd_vertex *parent_vertex = me_to_parent.p;

            size_t parent_vertex_index = parent_vertex->index;
            struct arr_c parent_to_me = vertex_edges[parent_vertex_index][me_to_parent.arr_c_index];

            size_t parent_edges_length = vertex_edges_length[parent_vertex_index];

            if (parent_vertex_index < i) {
                continue;
            }

            bool should_resize = false;
            size_t new_parent_edges_alloc_length = my_edges_length + parent_edges_length;

            while (new_parent_edges_alloc_length >= vertex_edges_alloc_length[parent_vertex_index]) {
                vertex_edges_alloc_length[parent_vertex_index] *= 2;
                should_resize = true;
            }

            if (should_resize) {
                vertex_edges[parent_vertex_index] = (struct arr_c *) realloc(
                        vertex_edges[parent_vertex_index],
                        vertex_edges_alloc_length[parent_vertex_index] * sizeof(*(vertex_edges[parent_vertex_index]))
                );
            }

            vertex_edges_length[parent_vertex_index] = 0;

            double parent_weight_to_me = parent_to_me.prob;
            double new_parent_total_prob = 0;

            if (memcpy(
                    old_edges_buffer, vertex_edges[parent_vertex_index],
                    sizeof(struct arr_c) * parent_edges_length
            ) != old_edges_buffer) {
                return NULL;
            }

            struct arr_c *new_parent_children = vertex_edges[parent_vertex_index];

            commands = add_command(
                    commands,
                    original_indices[parent_vertex_index],
                    original_indices[i],
                    parent_weight_to_me,
                    command_index++
            );

            size_t child_index = 0;
            size_t parent_child_index = 0;

            while (child_index < my_edges_length || parent_child_index < parent_edges_length) {
                struct arr_c me_to_child = my_children[child_index];
                struct ptd_vertex *me_to_child_v = me_to_child.to;
                struct arr_c parent_to_child = old_edges_buffer[parent_child_index];
                struct ptd_vertex *parent_to_child_v = parent_to_child.to;
                double me_to_child_p = me_to_child.prob;

                if (me_to_child_v == parent_vertex) {
                    double prob = parent_weight_to_me * me_to_child_p;
                    commands = add_command(
                            commands,
                            original_indices[parent_vertex->index],
                            original_indices[parent_vertex->index],
                            1 / (1 - prob),
                            command_index++
                    );

                    child_index++;
                    continue;
                }

                if (parent_to_child_v == me) {
                    parent_child_index++;
                    continue;
                }

                if (me_to_child_v == parent_to_child_v) {
                    new_parent_children[vertex_edges_length[parent_vertex_index]].to = parent_to_child_v;
                    new_parent_children[vertex_edges_length[parent_vertex_index]].prob =
                            parent_to_child.prob + me_to_child_p * parent_weight_to_me;
                    new_parent_children[vertex_edges_length[parent_vertex_index]].arr_p_index = parent_to_child.arr_p_index;
                    if (parent_to_child_v != dummy__ptd_min && parent_to_child_v != dummy__ptd_max) {
                        size_t current_parent_index = parent_to_child.arr_p_index;
                        vertex_parents[parent_to_child_v->index][current_parent_index].arr_c_index = vertex_edges_length[parent_vertex_index];

                    }
                    new_parent_total_prob += new_parent_children[vertex_edges_length[parent_vertex_index]].prob;
                    vertex_edges_length[parent_vertex_index]++;

                    child_index++;
                    parent_child_index++;
                } else if (me_to_child_v < parent_to_child_v) {
                    size_t child_parents_length = vertex_parents_length[me_to_child_v->index];

                    if (child_parents_length >= vertex_parents_alloc_length[me_to_child_v->index]) {
                        vertex_parents_alloc_length[me_to_child_v->index] *= 2;
                        vertex_parents[me_to_child_v->index] = (struct arr_p *) realloc(
                                vertex_parents[me_to_child_v->index],
                                vertex_parents_alloc_length[me_to_child_v->index] *
                                sizeof(*(vertex_parents[me_to_child_v->index]))
                        );
                    }

                    vertex_parents[me_to_child_v->index][child_parents_length].arr_c_index = vertex_edges_length[parent_vertex_index];
                    vertex_parents[me_to_child_v->index][child_parents_length].p = parent_vertex;

                    new_parent_children[vertex_edges_length[parent_vertex_index]].to = me_to_child_v;
                    new_parent_children[vertex_edges_length[parent_vertex_index]].prob =
                            me_to_child_p * parent_weight_to_me;
                    new_parent_children[vertex_edges_length[parent_vertex_index]].arr_p_index = child_parents_length;
                    new_parent_total_prob += me_to_child_p * parent_weight_to_me;

                    vertex_edges_length[parent_vertex_index]++;
                    vertex_parents_length[me_to_child_v->index]++;

                    child_index++;
                } else {
                    new_parent_children[vertex_edges_length[parent_vertex_index]] = parent_to_child;
                    vertex_parents[parent_to_child_v->index][parent_to_child.arr_p_index].arr_c_index = vertex_edges_length[parent_vertex_index];
                    new_parent_total_prob += parent_to_child.prob;
                    vertex_edges_length[parent_vertex_index]++;

                    parent_child_index++;
                }
            }


            // Make sure parent has rate of 1
            for (size_t j = 0; j < vertex_edges_length[parent_vertex_index]; ++j) {
                new_parent_children[j].prob /= new_parent_total_prob;
            }

            //free(vertex_edges[parent->p->index]);
            //vertex_edges[parent->p->index] = new_parent_children;
            vertex_edges_length[parent_vertex_index] = vertex_edges_length[parent_vertex_index];
        }
    }

    for (size_t ii = 0; ii < vertices_length; ++ii) {
        size_t i = vertices_length - ii - 1;
        struct ptd_vertex *vertex = vertices[i];


        for (size_t j = 1; j < vertex_edges_length[i] - 1; ++j) {
            struct arr_c child = vertex_edges[i][j];
            commands = add_command(
                    commands,
                    original_indices[vertex->index],
                    original_indices[child.to->index],
                    child.prob,
                    command_index++
            );
        }
    }

    for (size_t i = 0; i < vertices_length; ++i) {
        graph->vertices[i]->index = i;
    }

    for (size_t i = 0; i < vertices_length; ++i) {
        free(vertex_edges[i]);
        free(vertex_parents[i]);
    }

    free(vertex_parents_length);
    free(vertex_parents_alloc_length);
    free(vertex_parents);
    free(vertex_edges);
    free(vertex_edges_length);
    free(vertex_edges_alloc_length);
    free(original_indices);
    free(vertices);
    free(old_edges_buffer);
    free(v);
    ptd_scc_graph_destroy(scc);

    commands = add_command(
            commands,
            0,
            0,
            NAN,
            command_index
    );

    struct ptd_desc_reward_compute *res = (struct ptd_desc_reward_compute *) malloc(sizeof(*res));
    res->length = command_index;
    res->commands = commands;

    return res;
}


struct ll_c2_a {
    struct ll_c2_a *next;
    struct ll_c2 *mem;
};

static struct ll_c2_a **ll_c2_alloced;
static size_t ll_c2_alloced__ptd_max = 1024;
static size_t *ll_c2_alloced_index;

static void ll_c2_alloc_init(size_t length) {
    ll_c2_alloced_index = (size_t *) calloc(length, sizeof(*ll_c2_alloced_index));
    ll_c2_alloced = (struct ll_c2_a **) calloc(length, sizeof(*ll_c2_alloced));

    for (size_t i = 0; i < length; ++i) {
        ll_c2_alloced[i] = (struct ll_c2_a *) malloc(sizeof(*(ll_c2_alloced[i])));
        ll_c2_alloced[i]->next = NULL;
        ll_c2_alloced[i]->mem = (struct ll_c2 *) calloc(ll_c2_alloced__ptd_max, sizeof(struct ll_c2));
        ll_c2_alloced_index[i] = 0;
    }
}

static void ll_c2_alloc_init_free(size_t length) {
    free(ll_c2_alloced_index);
    free(ll_c2_alloced);
}

static struct ll_c2 *ll_c2_alloc(size_t index) {
    if (ll_c2_alloced_index[index] >= ll_c2_alloced__ptd_max) {
        struct ll_c2_a *old = ll_c2_alloced[index];
        ll_c2_alloced[index] = (struct ll_c2_a *) malloc(sizeof(*(ll_c2_alloced[index])));
        ll_c2_alloced[index]->next = old;
        ll_c2_alloced[index]->mem = (struct ll_c2 *) calloc(ll_c2_alloced__ptd_max, sizeof(struct ll_c2));
        ll_c2_alloced_index[index] = 0;
    }

    return &(ll_c2_alloced[index]->mem[ll_c2_alloced_index[index]++]);
}

static void ll_c2_free(size_t index) {
    struct ll_c2_a *old = ll_c2_alloced[index];

    while (old != NULL) {
        free(old->mem);
        struct ll_c2_a *next = old->next;
        free(old);
        old = next;
    }
}

struct ll_p2_a {
    struct ll_p2_a *next;
    struct ll_p2 *mem;
};

static struct ll_p2_a **ll_p2_alloced;
static size_t ll_p2_alloced__ptd_max = 1024;
static size_t *ll_p2_alloced_index;

static void ll_p2_alloc_init(size_t length) {
    ll_p2_alloced_index = (size_t *) calloc(length, sizeof(*ll_p2_alloced_index));
    ll_p2_alloced = (struct ll_p2_a **) calloc(length, sizeof(*ll_p2_alloced));

    for (size_t i = 0; i < length; ++i) {
        ll_p2_alloced[i] = (struct ll_p2_a *) malloc(sizeof(*(ll_p2_alloced[i])));
        ll_p2_alloced[i]->next = NULL;
        ll_p2_alloced[i]->mem = (struct ll_p2 *) calloc(ll_p2_alloced__ptd_max, sizeof(struct ll_p2));
        ll_p2_alloced_index[i] = 0;
    }
}

static void ll_p2_alloc_init_free(size_t length) {
    free(ll_p2_alloced);
    free(ll_p2_alloced_index);
}

static struct ll_p2 *ll_p2_alloc(size_t index) {
    if (ll_p2_alloced_index[index] >= ll_p2_alloced__ptd_max) {
        struct ll_p2_a *old = ll_p2_alloced[index];
        ll_p2_alloced[index] = (struct ll_p2_a *) malloc(sizeof(*(ll_p2_alloced[index])));
        ll_p2_alloced[index]->next = old;
        ll_p2_alloced[index]->mem = (struct ll_p2 *) calloc(ll_p2_alloced__ptd_max, sizeof(struct ll_p2));
        ll_p2_alloced_index[index] = 0;
    }

    return &(ll_p2_alloced[index]->mem[ll_p2_alloced_index[index]++]);
}

static void ll_p2_free(size_t index) {
    struct ll_p2_a *old = ll_p2_alloced[index];

    while (old != NULL) {
        free(old->mem);
        struct ll_p2_a *next = old->next;
        free(old);
        old = next;
    }
}

static int t = 0;

static struct ll_of_a *add_mem(struct ll_of_a *current_mem_ll, double what) {
    struct ll_of_a *n;

    if (current_mem_ll == NULL || current_mem_ll->current_mem_index >= 32768) {
        n = (struct ll_of_a *) malloc(sizeof(*n));
        n->next = current_mem_ll;
        n->mem = (double *) calloc(32768, sizeof(double));
        n->current_mem_index = 0;
        n->current_mem_position = n->mem;
        t++;
    } else {
        n = current_mem_ll;
    }

    n->mem[n->current_mem_index] = what;
    n->current_mem_position = &(n->mem[n->current_mem_index]);
    n->current_mem_index++;

    return n;
}

struct ptd_desc_reward_compute_parameterized *ptd_graph_ex_absorbation_time_comp_graph_parameterized(
        struct ptd_graph *graph
) {
    struct ptd_vertex *dummy__ptd_min = 0, *dummy__ptd_max = 0;

    struct ll_of_a *current_mem_ll = NULL;
    current_mem_ll = add_mem(current_mem_ll, 0);
    double *SIMPLE_ZERO = current_mem_ll->current_mem_position;

    struct ll_c2 **edges;
    struct ll_p2 **parents;

    struct ptd_vertex **vertices = (struct ptd_vertex **) calloc(graph->vertices_length, sizeof(*vertices));
    size_t *original_indices = (size_t *) calloc(graph->vertices_length, sizeof(*original_indices));
    edges = (struct ll_c2 **) calloc(graph->vertices_length, sizeof(*edges));
    parents = (struct ll_p2 **) calloc(graph->vertices_length, sizeof(*parents));
    ll_c2_alloc_init(1);
    ll_p2_alloc_init(1);
    struct ptd_comp_graph_parameterized *commands = NULL;
    size_t command_index = 0;
    size_t vertices_length = graph->vertices_length;


    struct ptd_scc_graph *scc = ptd_find_strongly_connected_components(graph);
    struct ptd_scc_vertex **v = ptd_scc_graph_topological_sort(scc);
    size_t idx = 0;

    for (size_t sii = 0; sii < scc->vertices_length; ++sii) {
        for (size_t j = 0; j < v[sii]->internal_vertices_length; ++j) {
            if (v[sii]->internal_vertices[j]->edges_length == 0) {
                continue;
            }

            original_indices[idx] = v[sii]->internal_vertices[j]->index;
            v[sii]->internal_vertices[j]->index = idx;
            vertices[idx] = v[sii]->internal_vertices[j];
            idx++;
        }
    }

    for (size_t sii = 0; sii < scc->vertices_length; ++sii) {
        for (size_t j = 0; j < v[sii]->internal_vertices_length; ++j) {
            if (v[sii]->internal_vertices[j]->edges_length != 0) {
                continue;
            }

            original_indices[idx] = v[sii]->internal_vertices[j]->index;
            v[sii]->internal_vertices[j]->index = idx;
            vertices[idx] = v[sii]->internal_vertices[j];
            idx++;
        }
    }

    double **rates = (double **) calloc(graph->vertices_length, sizeof(*rates));

    for (size_t i = 0; i < vertices_length; ++i) {
        struct ptd_vertex *vertex = vertices[i];

        if (vertex >= dummy__ptd_max) {
            dummy__ptd_max = vertex + 1;
        }

        if (vertex <= dummy__ptd_min) {
            dummy__ptd_min = vertex - 1;
        }

        current_mem_ll = add_mem(current_mem_ll, 0);
        rates[i] = current_mem_ll->current_mem_position;
        commands = add_command_param_zero(
                commands,
                rates[i],
                command_index++
        );

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            commands = add_command_param_p(
                    commands,
                    rates[i],
                    &(vertex->edges[j]->weight),
                    1,
                    command_index++
            );
        }

        commands = add_command_param_inverse(
                commands,
                rates[i],
                command_index++
        );

        // Add the "real" rate as our first reward

        if (graph->starting_vertex == vertex || vertex->edges_length == 0) {
            commands = add_command_param(
                    commands,
                    original_indices[i],
                    original_indices[i],
                    SIMPLE_ZERO,
                    command_index++
            );
        } else {
            commands = add_command_param(
                    commands,
                    original_indices[i],
                    original_indices[i],
                    rates[i],
                    command_index++
            );
        }
    }

    for (size_t i = 0; i < vertices_length; ++i) {
        struct ptd_vertex *vertex = vertices[i];

        struct ll_c2 *dummy_first = ll_c2_alloc(0);
        dummy_first->next = NULL;
        dummy_first->prev = NULL;
        dummy_first->weight = 0;
        dummy_first->c = dummy__ptd_min;
        dummy_first->ll_p = NULL;
        edges[i] = dummy_first;

        struct ll_c2 *last = dummy_first;

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            struct ll_p2 *n = ll_p2_alloc(0);

            n->next = parents[vertex->edges[j]->to->index];
            n->p = vertex;
            n->prev = NULL;

            if (parents[vertex->edges[j]->to->index] != NULL) {
                parents[vertex->edges[j]->to->index]->prev = n;
            }

            parents[vertex->edges[j]->to->index] = n;

            struct ll_c2 *nc = ll_c2_alloc(0);
            nc->next = NULL;

            nc->prev = last;
            last->next = nc;

            current_mem_ll = add_mem(current_mem_ll, 0);

            commands = add_command_param_zero(
                    commands,
                    current_mem_ll->current_mem_position,
                    command_index++
            );

            commands = add_command_param_pp(
                    commands,
                    current_mem_ll->current_mem_position,
                    &(vertex->edges[j]->weight),
                    rates[i],
                    command_index++
            );

            nc->weight = current_mem_ll->current_mem_position;

            nc->c = vertex->edges[j]->to;
            nc->ll_p = n;
            n->ll_c = nc;
            last = nc;
        }

        struct ll_c2 *dummy_last = ll_c2_alloc(0);
        dummy_last->next = NULL;
        dummy_last->prev = last;
        dummy_last->weight = 0;
        dummy_last->c = dummy__ptd_max;
        dummy_last->ll_p = NULL;
        last->next = dummy_last;
    }

    int ri = 0;

    for (size_t i = 0; i < vertices_length; ++i) {
        struct ptd_vertex *vertex = vertices[i];

        ri++;

        struct ll_p2 *parent = parents[i];

        struct ll_c2 *c = edges[i];
        size_t n_edges = 0;

        while (c != NULL) {
            n_edges += 1;
            c = c->next;
        }

        struct ll_c2 *children_arr = (struct ll_c2 *) calloc(n_edges, sizeof(*children_arr));
        c = edges[i];
        size_t l = 0;

        while (c != NULL) {
            children_arr[l] = *c;
            l++;
            c = c->next;
        }

        while (parent != NULL) {
            if (parent->p->index < i) {
                parent = parent->next;
                continue;
            }

            l = 0;
            struct ll_c2 *parent_child = edges[parent->p->index];
            double *parent_weight_to_me = parent->ll_c->weight;

            commands = add_command_param(
                    commands,
                    original_indices[parent->p->index],
                    original_indices[i],
                    parent_weight_to_me,
                    command_index++
            );

            while (children_arr[l].c != dummy__ptd_max) {
                double *prob = children_arr[l].weight;
                struct ptd_vertex *child_vertex = children_arr[l].c;
                struct ptd_vertex *parent_vertex = parent->p;
                struct ptd_vertex *parent_child_vertex = parent_child->c;

                if (child_vertex == parent_vertex) {
                    current_mem_ll = add_mem(current_mem_ll, 0);
                    double *p = current_mem_ll->current_mem_position;

                    commands = add_command_param_zero(
                            commands,
                            p,
                            command_index++
                    );

                    commands = add_command_param_pp(
                            commands,
                            p,
                            parent_weight_to_me,
                            prob,
                            command_index++
                    );

                    commands = add_command_param_one__ptd_minus(
                            commands,
                            p,
                            command_index++
                    );

                    commands = add_command_param_inverse(
                            commands,
                            p,
                            command_index++
                    );

                    commands = add_command_param(
                            commands,
                            original_indices[parent_vertex->index],
                            original_indices[parent_vertex->index],
                            p,
                            command_index++
                    );

                    l++;
                    continue;
                }

                if (parent_child_vertex == vertex) {
                    parent_child = parent_child->next;
                    continue;
                }

                if (child_vertex == parent_child_vertex) {
                    if (child_vertex != dummy__ptd_min) {
                        current_mem_ll = add_mem(current_mem_ll, 0);
                        double *p = current_mem_ll->current_mem_position;

                        commands = add_command_param_zero(
                                commands,
                                p,
                                command_index++
                        );

                        commands = add_command_param_pp(
                                commands,
                                p,
                                parent_weight_to_me,
                                prob,
                                command_index++
                        );

                        commands = add_command_param_p(
                                commands,
                                parent_child->weight,
                                p,
                                1,
                                command_index++
                        );
                    }

                    l++;
                    parent_child = parent_child->next;
                } else if (child_vertex < parent_child_vertex) {
                    current_mem_ll = add_mem(current_mem_ll, 0);
                    double *p = current_mem_ll->current_mem_position;
                    commands = add_command_param_zero(
                            commands,
                            p,
                            command_index++
                    );

                    commands = add_command_param_pp(
                            commands,
                            p,
                            parent_weight_to_me,
                            prob,
                            command_index++
                    );

                    struct ll_c2 *to = ll_c2_alloc(0);
                    to->c = child_vertex;

                    current_mem_ll = add_mem(current_mem_ll, 0);
                    commands = add_command_param_zero(
                            commands,
                            current_mem_ll->current_mem_position,
                            command_index++
                    );
                    to->weight = current_mem_ll->current_mem_position;

                    commands = add_command_param_p(
                            commands,
                            to->weight,
                            p,
                            1,
                            command_index++
                    );
                    to->next = parent_child;
                    to->prev = parent_child->prev;


                    struct ll_p2 *ll_p = ll_p2_alloc(0);
                    ll_p->next = parents[child_vertex->index];
                    parents[child_vertex->index]->prev = ll_p;
                    parents[child_vertex->index] = ll_p;
                    ll_p->prev = NULL;
                    ll_p->p = parent_vertex;

                    ll_p->ll_c = to;
                    to->ll_p = ll_p;

                    to->next = parent_child;
                    to->prev = parent_child->prev;
                    parent_child->prev->next = to;
                    parent_child->prev = to;

                    l++;
                } else {
                    parent_child = parent_child->next;
                }
            }

            struct ll_c2 *edge_to_me = parent->ll_c;
            edge_to_me->prev->next = edge_to_me->next;
            edge_to_me->next->prev = edge_to_me->prev;

            // Make sure parent has rate of 1
            current_mem_ll = add_mem(current_mem_ll, 0);
            double *rate = current_mem_ll->current_mem_position;
            commands = add_command_param_zero(
                    commands,
                    rate,
                    command_index++
            );

            parent_child = edges[parent->p->index]->next;

            while (parent_child->c != dummy__ptd_max) {
                commands = add_command_param_p(
                        commands,
                        rate,
                        parent_child->weight,
                        1,
                        command_index++
                );

                parent_child = parent_child->next;
            }

            parent_child = edges[parent->p->index]->next;

            while (parent_child->c != dummy__ptd_max) {
                commands = add_command_param_p_divide(
                        commands,
                        parent_child->weight,
                        rate,
                        command_index++
                );

                parent_child = parent_child->next;
            }

            parent_child = edges[parent->p->index]->next;

            while (parent_child->c != dummy__ptd_max) {
                parent_child = parent_child->next;
            }

            parent = parent->next;
        }

        struct ll_c2 *child = edges[i]->next;

        while (child->c != dummy__ptd_max) {
            if (child->ll_p->prev != NULL) {
                if (child->ll_p->next != NULL) {
                    child->ll_p->next->prev = child->ll_p->prev;
                    child->ll_p->prev->next = child->ll_p->next;
                } else {
                    child->ll_p->prev->next = NULL;
                }
            } else {
                if (child->ll_p->next != NULL) {
                    child->ll_p->next->prev = NULL;
                }

                parents[child->c->index] = child->ll_p->next;
            }

            child = child->next;
        }

        free(children_arr);
    }

    for (size_t ii = 0; ii < vertices_length; ++ii) {
        size_t i = vertices_length - ii - 1;
        struct ptd_vertex *vertex = vertices[i];

        struct ll_c2 *child = edges[vertex->index]->next;

        while (child->c != dummy__ptd_max) {
            commands = add_command_param(
                    commands,
                    original_indices[vertex->index],
                    original_indices[child->c->index],
                    child->weight,
                    command_index++
            );
            child = child->next;
        }
    }

    for (size_t i = 0; i < vertices_length; ++i) {
        graph->vertices[i]->index = i;
    }


    free(original_indices);
    free(vertices);
    free(edges);
    free(parents);
    free(v);
    ptd_scc_graph_destroy(scc);
    ll_c2_free(0);
    ll_p2_free(0);
    ll_c2_alloc_init_free(1);
    ll_p2_alloc_init_free(1);

    commands = add_command_param(
            commands,
            0,
            0,
            NULL,
            command_index
    );

    struct ptd_desc_reward_compute_parameterized *res = (struct ptd_desc_reward_compute_parameterized *) malloc(
            sizeof(*res)
    );
    res->length = command_index;
    res->commands = commands;
    res->mem = current_mem_ll;
    res->memr = rates;

    return res;
}

struct ptd_desc_reward_compute *ptd_graph_build_ex_absorbation_time_comp_graph_parameterized(
        struct ptd_desc_reward_compute_parameterized *compute
) {
    struct ptd_reward_increase *commands = NULL;
    size_t command_index = 0;
    enum command_types {
        PP = 3,
        P = 1,
        INV = 2,
        ZERO = 6,
        DIVIDE = 5,
        ONE_MINUS = 4,
        NEW_ADD = 0
    };
    for (size_t i = 0; i < compute->length; ++i) {
        struct ptd_comp_graph_parameterized command = compute->commands[i];

        switch (command.type) {
            case NEW_ADD:
                commands = add_command(
                        commands,
                        command.from,
                        command.to,
                        *command.multiplierptr,
                        command_index++
                );
                break;
            case P:
                *(command.fromT) = *(command.fromT) + *command.toT * command.multiplier;
                break;
            case PP:
                *(command.fromT) = *(command.fromT) + *command.toT * *command.multiplierptr;
                break;
            case INV:
                *(command.fromT) = 1 / *(command.fromT);
                break;
            case ONE_MINUS:
                *(command.fromT) = 1 - *command.fromT;
                break;
            case DIVIDE:
                *(command.fromT) /= *command.toT;
                break;
            case ZERO:
                *command.fromT = 0;
                break;
            default:
                DIE_ERROR(1, "Unknown command\n");
        }
    }

    struct ptd_desc_reward_compute *res = (struct ptd_desc_reward_compute *) malloc(sizeof(*res));
    res->length = command_index;
    res->commands = commands;

    return res;
}


double *ptd_expected_waiting_time(struct ptd_graph *graph, double *rewards) {
    if (ptd_precompute_reward_compute_graph(graph)) {
        return NULL;
    }

    double *result = (double *) calloc(graph->vertices_length, sizeof(*result));

    if (rewards != NULL) {
        // TODO: fix this if reward is nan...
        memcpy(result, rewards, sizeof(*result) * graph->vertices_length);
    } else {
        for (size_t j = 0; j < graph->vertices_length; ++j) {
            result[j] = 1;
        }
    }

    for (size_t j = 0; j < graph->reward_compute_graph->length; ++j) {
        struct ptd_reward_increase command = graph->reward_compute_graph->commands[j];

        result[command.from] += result[command.to] * command.multiplier;

        //TODO: if inf, give error stating that there is an infinite loop
    }

    return result;
}

long double ptd_random_sample(struct ptd_graph *graph, double *rewards) {
    long double outcome = 0;

    struct ptd_vertex *vertex = graph->starting_vertex;

    while (vertex->edges_length != 0) {
        long double draw_wait = (long double) rand() / (long double) RAND_MAX;

        double rate = 0;

        for (size_t i = 0; i < vertex->edges_length; ++i) {
            long double edge_weight = vertex->edges[i]->weight;
            rate += edge_weight;
        }

        long double waiting_time = -logl(draw_wait + 0.0000001) / rate;

        if (rewards != NULL) {
            waiting_time *= rewards[vertex->index];
        }

        if (vertex == graph->starting_vertex) {
            waiting_time = 0;
        }

        outcome += waiting_time;

        long double draw_direction = (long double) rand() / (long double) RAND_MAX;
        long double weight_sum = 0;
        size_t edge_index = 0;

        for (size_t i = 0; i < vertex->edges_length; ++i) {
            long double edge_weight = vertex->edges[i]->weight;
            weight_sum += edge_weight;

            if (weight_sum / rate >= draw_direction) {
                edge_index = i;
                break;
            }
        }

        vertex = vertex->edges[edge_index]->to;
    }

    return outcome;
}

long double *ptd_mph_random_sample(struct ptd_graph *graph, double *rewards, size_t vertex_rewards_length) {
    long double *outcome = (long double *) calloc(vertex_rewards_length, sizeof(*outcome));

    for (size_t j = 0; j < vertex_rewards_length; ++j) {
        outcome[j] = 0;
    }

    struct ptd_vertex *vertex = graph->starting_vertex;

    while (vertex->edges_length != 0) {
        long double draw_wait = (long double) rand() / (long double) RAND_MAX;

        double rate = 0;

        for (size_t i = 0; i < vertex->edges_length; ++i) {
            long double edge_weight = vertex->edges[i]->weight;
            rate += edge_weight;
        }

        long double waiting_time = -logl(draw_wait + 0.0000001) / rate;

        if (vertex != graph->starting_vertex) {
            for (size_t i = 0; i < vertex_rewards_length; ++i) {
                outcome[i] += waiting_time * rewards[vertex->index * vertex_rewards_length + i];
            }
        }

        long double draw_direction = (long double) rand() / (long double) RAND_MAX;
        long double weight_sum = 0;
        size_t edge_index = 0;

        for (size_t i = 0; i < vertex->edges_length; ++i) {
            long double edge_weight = vertex->edges[i]->weight;
            weight_sum += edge_weight;

            if (weight_sum / rate >= draw_direction) {
                edge_index = i;
                break;
            }
        }

        vertex = vertex->edges[edge_index]->to;
    }

    return outcome;
}


long double ptd_dph_random_sample(struct ptd_graph *graph, double *rewards) {
    long double jumps = 0;

    struct ptd_vertex *vertex = graph->starting_vertex;

    while (vertex->edges_length != 0) {
        long double draw_direction = (long double) rand() / (long double) RAND_MAX;
        long double weight_sum = 0;
        int edge_index = -1;

        double rate = 0;

        for (size_t i = 0; i < vertex->edges_length; ++i) {
            long double edge_weight = vertex->edges[i]->weight;
            rate += edge_weight;
        }

        if (rate > 1.0001) {
            size_t debug_index = vertex->index;

            if (PTD_DEBUG_1_INDEX) {
                debug_index++;
            }

            char state[1024] = {'\0'};
            char starting_vertex[] = " (starting vertex)";

            if (vertex != graph->starting_vertex) {
                starting_vertex[0] = '\0';
            }

            ptd_vertex_to_s(vertex, state, 1023);

            snprintf(
                    (char *) ptd_err,
                    sizeof(ptd_err),
                    "Expected vertex with index %i%s (state %s) to have outgoing rate <= 1. Is '%f'. Are you sure this is a discrete phase-type distribution?\n",
                    (int) debug_index, starting_vertex, state, (float) rate
            );

            return NAN;
        }

        for (int i = 0; i < (int) vertex->edges_length; ++i) {
            long double edge_weight = vertex->edges[i]->weight;
            weight_sum += edge_weight;

            if (weight_sum >= draw_direction) {
                edge_index = i;
                break;
            }
        }

        if (vertex != graph->starting_vertex) {
            if (rewards == NULL) {
                jumps += 1;
            } else {
                jumps += rewards[vertex->index];
            }
        }

        if (edge_index != -1) {
            vertex = vertex->edges[edge_index]->to;
        }
    }

    return jumps;
}

long double *ptd_mdph_random_sample(struct ptd_graph *graph, double *rewards, size_t vertex_rewards_length) {
    long double *jumps = (long double *) calloc(vertex_rewards_length, sizeof(*jumps));

    for (size_t j = 0; j < vertex_rewards_length; ++j) {
        jumps[j] = 0;
    }

    struct ptd_vertex *vertex = graph->starting_vertex;

    while (vertex->edges_length != 0) {
        long double draw_direction = (long double) rand() / (long double) RAND_MAX;
        long double weight_sum = 0;
        int edge_index = -1;
        double rate = 0;

        for (size_t i = 0; i < vertex->edges_length; ++i) {
            long double edge_weight = vertex->edges[i]->weight;
            rate += edge_weight;
        }

        if (rate > 1.0001) {
            size_t debug_index = vertex->index;

            if (PTD_DEBUG_1_INDEX) {
                debug_index++;
            }

            char state[1024] = {'\0'};
            char starting_vertex[] = " (starting vertex)";

            if (vertex != graph->starting_vertex) {
                starting_vertex[0] = '\0';
            }

            ptd_vertex_to_s(vertex, state, 1023);

            snprintf(
                    (char *) ptd_err,
                    sizeof(ptd_err),
                    "Expected vertex with index %i%s (state %s) to have outgoing rate <= 1. Is '%f'. Are you sure this is a discrete phase-type distribution?\n",
                    (int) debug_index, starting_vertex, state, (float) rate
            );

            free(jumps);

            return NULL;
        }

        for (int i = 0; i < (int) vertex->edges_length; ++i) {
            long double edge_weight = vertex->edges[i]->weight;
            weight_sum += edge_weight;

            if (weight_sum >= draw_direction) {
                edge_index = i;
                break;
            }
        }


        if (vertex != graph->starting_vertex) {
            for (size_t i = 0; i < vertex_rewards_length; ++i) {
                jumps[i] += rewards[vertex->index * vertex_rewards_length + i];
            }
        }


        if (edge_index != -1) {
            vertex = vertex->edges[edge_index]->to;
        }
    }

    return jumps;
}

struct ptd_vertex *ptd_random_sample_stop_vertex(struct ptd_graph *graph, double time) {
    double time_spent = 0;

    struct ptd_vertex *vertex = graph->starting_vertex;

    while (vertex->edges_length != 0) {
        long double draw_wait = (long double) rand() / (long double) RAND_MAX;

        double rate = 0;

        for (size_t i = 0; i < vertex->edges_length; ++i) {
            long double edge_weight = vertex->edges[i]->weight;
            rate += edge_weight;
        }

        long double waiting_time = -logl(draw_wait + 0.0000001) / rate;

        if (vertex == graph->starting_vertex) {
            waiting_time = 0;
        }

        time_spent += waiting_time;

        if (time_spent >= time && vertex != graph->starting_vertex) {
            return vertex;
        }

        long double draw_direction = (long double) rand() / (long double) RAND_MAX;
        long double weight_sum = 0;
        size_t edge_index = 0;

        for (size_t i = 0; i < vertex->edges_length; ++i) {
            long double edge_weight = vertex->edges[i]->weight;
            weight_sum += edge_weight;

            if (weight_sum / rate >= draw_direction) {
                edge_index = i;
                break;
            }
        }

        vertex = vertex->edges[edge_index]->to;
    }

    return vertex;
}

struct ptd_vertex *ptd_dph_random_sample_stop_vertex(struct ptd_graph *graph, int jumps) {
    int jumps_taken = -1;

    struct ptd_vertex *vertex = graph->starting_vertex;

    while (vertex->edges_length != 0 && jumps < jumps_taken) {
        long double draw_direction = (long double) rand() / (long double) RAND_MAX;
        long double weight_sum = 0;
        int edge_index = -1;

        double rate = 0;

        for (size_t i = 0; i < vertex->edges_length; ++i) {
            long double edge_weight = vertex->edges[i]->weight;
            rate += edge_weight;
        }

        if (rate > 1.0001) {
            size_t debug_index = vertex->index;

            if (PTD_DEBUG_1_INDEX) {
                debug_index++;
            }

            char state[1024] = {'\0'};
            char starting_vertex[] = " (starting vertex)";

            if (vertex != graph->starting_vertex) {
                starting_vertex[0] = '\0';
            }

            ptd_vertex_to_s(vertex, state, 1023);

            snprintf(
                    (char *) ptd_err,
                    sizeof(ptd_err),
                    "Expected vertex with index %i%s (state %s) to have outgoing rate <= 1. Is '%f'. Are you sure this is a discrete phase-type distribution?\n",
                    (int) debug_index, starting_vertex, state, (float) rate
            );

            return NULL;
        }

        for (int i = 0; i < (int) vertex->edges_length; ++i) {
            long double edge_weight = vertex->edges[i]->weight;
            weight_sum += edge_weight;

            if (weight_sum >= draw_direction) {
                edge_index = i;
                break;
            }
        }

        if (vertex != graph->starting_vertex) {
            jumps += 1;
        }

        if (edge_index != -1) {
            vertex = vertex->edges[edge_index]->to;
        }
    }

    return vertex;
}


struct ptd_vertex *
ptd_find_or_create_vertex(struct ptd_graph *graph, struct ptd_avl_tree *avl_tree, const int *child_state) {
    struct ptd_vertex *child;
    struct ptd_avl_node *avl_node = ptd_avl_tree_find(avl_tree, child_state);

    if (avl_node == NULL) {
        child = ptd_vertex_create(graph);
        memcpy(child->state, child_state, graph->state_length * sizeof(int));

        ptd_avl_tree_find_or_insert(avl_tree, child->state, child);
    } else {
        child = (struct ptd_vertex *) avl_node->entry;
    }

    return child;
}

struct dph_prob_increment {
    size_t from;
    size_t to;
    double *weight;
};

struct ptd_dph_probability_distribution_context *_ptd_dph_probability_distribution_context_create(
        struct ptd_graph *graph,
        bool dont_worry
) {
    if (!dont_worry) {
        for (size_t i = 0; i < graph->vertices_length; ++i) {
            double rate = 0;

            struct ptd_vertex *vertex = graph->vertices[i];

            for (size_t j = 0; j < vertex->edges_length; ++j) {
                rate += vertex->edges[j]->weight;
            }

            if (rate > 1.0001) {
                size_t debug_index = vertex->index;

                if (PTD_DEBUG_1_INDEX) {
                    debug_index++;
                }

                char state[1024] = {'\0'};
                char starting_vertex[] = " (starting vertex)";

                if (vertex != graph->starting_vertex) {
                    starting_vertex[0] = '\0';
                }

                ptd_vertex_to_s(vertex, state, 1023);

                snprintf(
                        (char *) ptd_err,
                        sizeof(ptd_err),
                        "Expected vertex with index %i%s (state %s) to have outgoing rate <= 1. Is '%f'. Are you sure this is a discrete phase-type distribution?\n",
                        (int) debug_index, starting_vertex, state, (float) rate
                );

                return NULL;
            }
        }
    }

    struct ptd_dph_probability_distribution_context *res =
            (struct ptd_dph_probability_distribution_context *) malloc(sizeof(*res));

    res->graph = graph;
    res->probability_at = (long double *) calloc(
            graph->vertices_length,
            sizeof(*(res->probability_at))
    );
    res->accumulated_visits = (long double *) calloc(
            graph->vertices_length,
            sizeof(*(res->accumulated_visits))
    );

    size_t number_of_edges = 0;

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        struct ptd_vertex *vertex = graph->vertices[i];
        res->accumulated_visits[i] = 0;
        number_of_edges += vertex->edges_length;
    }

    res->priv2 = number_of_edges;
    res->priv3 = 1;

    res->priv = calloc(number_of_edges, sizeof(struct dph_prob_increment));
    struct dph_prob_increment *inc_list = (struct dph_prob_increment *) res->priv;
    size_t inc_index = 0;

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        struct ptd_vertex *vertex = graph->vertices[i];

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            inc_list[inc_index].from = i;
            inc_list[inc_index].to = vertex->edges[j]->to->index;
            inc_list[inc_index].weight = &(vertex->edges[j]->weight);

            inc_index++;
        }
    }

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        res->probability_at[i] = 0;
    }

    res->probability_at[0] = 1;

    res->cdf = 0;
    res->pmf = 0;
    res->jumps = 0;

    ptd_dph_probability_distribution_step(res);

    res->jumps = 0;

    return res;
}

struct ptd_dph_probability_distribution_context *ptd_dph_probability_distribution_context_create(
        struct ptd_graph *graph
) {
    return _ptd_dph_probability_distribution_context_create(graph, false);
}

void ptd_dph_probability_distribution_context_destroy(struct ptd_dph_probability_distribution_context *context) {
    if (context == NULL) {
        return;
    }

    free(context->accumulated_visits);
    free(context->probability_at);
    free(context->priv);
    free(context);
}

int ptd_dph_probability_distribution_step(
        struct ptd_dph_probability_distribution_context *context
) {
    context->jumps++;
    context->pmf = 0;

    long double *old_probability_at = (long double *) calloc(
            context->graph->vertices_length, sizeof(*old_probability_at)
    );

    memcpy(
            old_probability_at,
            context->probability_at,
            sizeof(*old_probability_at) * context->graph->vertices_length
    );

    for (size_t i = 0; i < context->graph->vertices_length; ++i) {
        old_probability_at[i] = context->probability_at[i];
        struct ptd_vertex *vertex = context->graph->vertices[i];

        if (vertex->edges_length == 0) {
            context->probability_at[i] = 0;
        }
    }

    for (size_t i = 0; i < context->priv2; ++i) {
        struct dph_prob_increment inc = ((struct dph_prob_increment *) (context->priv))[i];
        long double add = old_probability_at[inc.from] * (*inc.weight) * context->priv3;
        context->probability_at[inc.to] += add;
        context->probability_at[inc.from] -= add;
    }

    for (size_t i = 0; i < context->graph->vertices_length; ++i) {
        struct ptd_vertex *vertex = context->graph->vertices[i];

        if (vertex->edges_length == 0) {
            context->pmf += context->probability_at[i];
            context->probability_at[i] = 0;
        } else {
            context->accumulated_visits[i] += context->probability_at[i];
        }
    }

    context->accumulated_visits[0] = 0;
    context->probability_at[0] = 0;

    context->cdf += context->pmf;

    free(old_probability_at);

    return 0;
}

struct ptd_probability_distribution_context *ptd_probability_distribution_context_create(
        struct ptd_graph *graph,
        int granularity
) {
    double max_rate = 1024;

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        double rate = 0;

        struct ptd_vertex *vertex = graph->vertices[i];

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            rate += vertex->edges[j]->weight;
        }

        if (rate > max_rate) {
            max_rate = rate;
        }
    }

    if (granularity == 0) {
        granularity = max_rate;
    }

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        double rate = 0;

        struct ptd_vertex *vertex = graph->vertices[i];

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            rate += vertex->edges[j]->weight;
        }

        if (rate / granularity > 1.0001) {
            size_t debug_index = vertex->index;

            if (PTD_DEBUG_1_INDEX) {
                debug_index++;
            }

            char state[1024] = {'\0'};
            char starting_vertex[] = " (starting vertex)";

            if (vertex != graph->starting_vertex) {
                starting_vertex[0] = '\0';
            }

            ptd_vertex_to_s(vertex, state, 1023);

            snprintf(
                    (char *) ptd_err,
                    sizeof(ptd_err),
                    "Expected vertex with index %i%s (state %s) to have outgoing rate divided by granularity <= 1. Rate is '%f' ('%f'). Increase the granularity\n",
                    (int) debug_index, starting_vertex, state, (float) rate, (float) (rate / granularity)
            );

            return NULL;
        }
    }

    struct ptd_probability_distribution_context *res = (struct ptd_probability_distribution_context *)
            malloc(sizeof(*res));

    struct ptd_dph_probability_distribution_context *dph_res =
            _ptd_dph_probability_distribution_context_create(graph, true);
    dph_res->priv3 = (double) 1.0 / granularity;

    long double cdf1 = dph_res->cdf * granularity;

    ptd_dph_probability_distribution_step(dph_res);

    long double cdf2 = dph_res->cdf * granularity;

    ptd_dph_probability_distribution_context_destroy(dph_res);

    dph_res = _ptd_dph_probability_distribution_context_create(graph, true);
    dph_res->priv3 = (double) 1.0 / granularity;

    res->cdf = dph_res->cdf;
    res->pdf = (double) ((cdf2 - cdf1));
    res->graph = dph_res->graph;
    res->probability_at = dph_res->probability_at;
    res->accumulated_visits = dph_res->accumulated_visits;

    res->time = 0;
    res->priv = (void *) dph_res;
    res->granularity = granularity;

    return res;
}

void ptd_probability_distribution_context_destroy(struct ptd_probability_distribution_context *context) {
    if (context == NULL) {
        return;
    }

    ptd_dph_probability_distribution_context_destroy(
            (struct ptd_dph_probability_distribution_context *) context->priv
    );

    free(context);
}

int ptd_probability_distribution_step(
        struct ptd_probability_distribution_context *context
) {
    struct ptd_dph_probability_distribution_context *dph_context =
            (struct ptd_dph_probability_distribution_context *) context->priv;

    ptd_dph_probability_distribution_step(dph_context);

    context->time = ((long double) dph_context->jumps) / context->granularity;
    context->cdf = dph_context->cdf;
    context->pdf = dph_context->pmf * context->granularity;

    return 0;
}

double ptd_defect(struct ptd_graph *graph) {
    double rate = 0;

    for (size_t i = 0; i < graph->starting_vertex->edges_length; ++i) {
        struct ptd_edge *edge = graph->starting_vertex->edges[i];

        rate += edge->weight;
    }

    double defect = 0;

    for (size_t i = 0; i < graph->starting_vertex->edges_length; ++i) {
        struct ptd_edge *edge = graph->starting_vertex->edges[i];

        if (edge->to->edges_length == 0) {
            defect += edge->weight / rate;
        }
    }

    return defect;
}

struct ptd_clone_res ptd_clone_graph(struct ptd_graph *graph, struct ptd_avl_tree *avl_tree) {
    struct ptd_graph *res = ptd_graph_create(graph->state_length);

    for (size_t i = 1; i < graph->vertices_length; ++i) {
        ptd_vertex_create(res);
    }

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        struct ptd_vertex *v = graph->vertices[i];
        struct ptd_vertex *v2 = res->vertices[i];

        if (v->state != NULL) {
            memcpy(v2->state, v->state, sizeof(int) * res->state_length);
        }
    }

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        struct ptd_vertex *v = graph->vertices[i];
        struct ptd_vertex *v2 = res->vertices[i];

        for (size_t j = 0; j < v->edges_length; ++j) {
            struct ptd_edge *e = v->edges[j];

            if (e->parameterized) {
                ptd_graph_add_edge_parameterized(
                        v2,
                        res->vertices[e->to->index],
                        e->weight,
                        ((struct ptd_edge_parameterized *) e)->state
                )->should_free_state = false;
            } else {
                ptd_graph_add_edge(v2, res->vertices[e->to->index], e->weight);
            }
        }
    }

    struct ptd_avl_tree *new_tree = ptd_avl_tree_create(avl_tree->key_length);

    struct ptd_stack *stack = stack_create();

    stack_push(stack, avl_tree->root);

    while (!stack_empty(stack)) {
        struct ptd_avl_node *v = (struct ptd_avl_node *) stack_pop(stack);

        if (v == NULL) {
            continue;
        }

        ptd_avl_tree_find_or_insert(
                new_tree,
                v->key,
                res->vertices[((struct ptd_vertex *) v->entry)->index]
        );

        stack_push(stack, v->left);
        stack_push(stack, v->right);
    }

    stack_destroy(stack);

    struct ptd_clone_res ret;
    ret.graph = res;
    ret.avl_tree = new_tree;

    return ret;
}

/*
 * Utilities
 */


static struct ptd_vector *vector_create() {
    struct ptd_vector *vector = (struct ptd_vector *) malloc(sizeof(*vector));

    vector->entries = 0;
    vector->arr = NULL;

    return vector;
}

static int vector_add(struct ptd_vector *vector, void *entry) {
    bool is_power_of_2 = (vector->entries & (vector->entries - 1)) == 0;

    if (is_power_of_2) {
        size_t new_length = vector->entries == 0 ? 1 : vector->entries * 2;

        if ((vector->arr = (void **) realloc(
                vector->arr,
                new_length * sizeof(void *))
            ) == NULL) {
            return -1;
        }
    }

    vector->arr[vector->entries] = entry;
    vector->entries++;

    return 0;
}

static void *vector_get(struct ptd_vector *vector, size_t index) {
    return vector->arr[index];
}

static size_t vector_length(struct ptd_vector *vector) {
    return vector->entries;
}

static void vector_destroy(struct ptd_vector *vector) {
    free(vector->arr);
    free(vector);
}

static struct ptd_queue *queue_create() {
    struct ptd_queue *queue = (struct ptd_queue *) malloc(sizeof(struct ptd_queue));

    queue->ll = NULL;
    queue->tail = NULL;

    return queue;
}

static void queue_destroy(struct ptd_queue *queue) {
    free(queue->ll);
    free(queue);
}

static int queue_enqueue(struct ptd_queue *queue, void *entry) {
    struct ptd_ll *new_ll = (struct ptd_ll *) malloc(sizeof(*new_ll));
    new_ll->next = NULL;
    new_ll->value = entry;

    if (queue->tail != NULL) {
        queue->tail->next = new_ll;
    } else {
        queue->ll = new_ll;
    }

    queue->tail = new_ll;

    return 0;
}

static void *queue_dequeue(struct ptd_queue *queue) {
    void *result = queue->ll->value;
    struct ptd_ll *value = queue->ll;
    queue->ll = queue->ll->next;

    if (queue->tail == value) {
        queue->tail = NULL;
    }

    free(value);

    return result;
}

static int queue_empty(struct ptd_queue *queue) {
    return (queue->tail == NULL);
}

static struct ptd_stack *stack_create() {
    struct ptd_stack *stack = (struct ptd_stack *) malloc(sizeof(struct ptd_stack));
    stack->ll = NULL;

    return stack;
}

static void stack_destroy(struct ptd_stack *stack) {
    free(stack->ll);
    free(stack);
}

static int stack_push(struct ptd_stack *stack, void *entry) {
    struct ptd_ll *new_ll = (struct ptd_ll *) malloc(sizeof(*new_ll));
    new_ll->next = stack->ll;
    new_ll->value = entry;

    stack->ll = new_ll;

    return 0;
}

static void *stack_pop(struct ptd_stack *stack) {
    void *result = stack->ll->value;
    struct ptd_ll *ll = stack->ll;

    stack->ll = stack->ll->next;
    free(ll);

    return result;
}

static int stack_empty(struct ptd_stack *stack) {
    return (stack->ll == NULL);
}
