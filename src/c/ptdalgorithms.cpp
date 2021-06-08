#include <iostream>
#include <vector>
#include <cstring>
#include <stack>
#include <queue>
#include <stdint.h>
#include <stddef.h>
#include <math.h>
#include <set>
#include "ptdalgorithms.h"

// TODO: name these ptd_...

extern void *create_matrix(double **mat, size_t length);

extern void *matrix_invert(void *matrix, size_t size);

extern double matrix_get(void *matrix, size_t i, size_t j);

extern double matrix_get(void *matrix, size_t i, size_t j);

extern void *matrix_init(size_t size);

extern void matrix_destroy(void *matrix, size_t size);

#ifdef PTD_USE_GSL

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <ptd.h>

void *create_matrix(double **mat, size_t length) {
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

void matrix_set(void *matrix, size_t i, size_t j, double x) {
    gsl_matrix *gsl_mat = (gsl_matrix *) matrix;

    gsl_matrix_set(gsl_mat, i, j, x);
}

void *matrix_init(size_t size) {
    return gsl_matrix_calloc(size, size);
}

void matrix_destroy(void *matrix, size_t size) {
    gsl_matrix_free((gsl_matrix *) matrix);
}

#endif // PTD_USE_GSL

#include "ptdalgorithms.h"

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
            free(*out);
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

    ptd_ph_vertex_destroy((struct ptd_ph_vertex *) avl_vertex->entry);
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


int ptd_avl_tree_vertex_insert(ptd_avl_tree_t *avl_tree, const int *key, const struct ptd_ph_vertex *vertex) {
    int res;

    struct avl_node *root = (struct avl_node *) avl_tree->root;
    res = avl_vec_insert(&root, (char *) key, (void *) vertex, avl_tree->vec_length * sizeof(int));

    if (res != 0) {
        return res;
    }

    avl_tree->root = root;

    return 0;
}

struct ptd_ph_vertex *ptd_avl_tree_vertex_find(const ptd_avl_tree_t *avl_tree, const int *key) {
    const struct avl_node *avl_vertex = avl_vec_find(
            (struct avl_node *) avl_tree->root,
            (char *) key,
            avl_tree->vec_length * sizeof(int)
    );


    if (avl_vertex == NULL) {
        return NULL;
    }

    return (struct ptd_ph_vertex *) avl_vertex->entry;
}

int ptd_avl_tree_edge_insert_or_increment(ptd_avl_tree_t *avl_tree, const int *key, struct ptd_ph_vertex *vertex,
                                          double weight) {
    struct avl_node *root = (struct avl_node *) avl_tree->root;
    struct avl_node *child;

    if (root == NULL) {
        struct ptd_ph_edge *edge = (struct ptd_ph_edge *) malloc(sizeof(*edge));
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
        int res = memcmp(parent->key, key, avl_tree->vec_length * sizeof(int));

        if (res < 0) {
            if (parent->left == NULL) {
                struct ptd_ph_edge *edge = (struct ptd_ph_edge *) malloc(sizeof(*edge));
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
                struct ptd_ph_edge *edge = (struct ptd_ph_edge *) malloc(sizeof(*edge));
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
            ((struct ptd_ph_edge *) parent->entry)->weight += weight;

            return 0;
        }
    }

    avl_rebalance_tree(&root, child);
    avl_tree->root = root;

    return 0;


    return 0;
}

int ptd_avl_tree_edge_remove(ptd_avl_tree_t *avl_tree, const int *key) {
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
        int res = memcmp(parent->key, key, avl_tree->vec_length * sizeof(int));

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

struct ptd_ph_edge *ptd_avl_tree_edge_find(const ptd_avl_tree_t *avl_tree, const int *key) {
    struct avl_node *parent = (struct avl_node *) avl_tree->root;

    while (true) {
        if (parent == NULL) {
            return NULL;
        }

        int res = memcmp(parent->key, (char *) key, avl_tree->vec_length * sizeof(int));

        if (res < 0) {
            parent = parent->left;
        } else if (res > 0) {
            parent = parent->right;
        } else {
            return (struct ptd_ph_edge *) parent->entry;
        }
    }
}

int ptd_ph_precompute_strongly_connected_components(struct ptd_ph_graph *graph) {
    if (graph->scc_graph == NULL) {
        graph->scc_graph = ptd_ph_find_strongly_connected_components(graph);

        if (graph->scc_graph == NULL) {
            return -1;
        }
    }

    return 0;
}

stack<struct ptd_ph_vertex *> *scc_stack2;
vector<struct ptd_ph_scc_vertex *> *scc_components2;
size_t scc_index2;
size_t *scc_indices2;
size_t *low_links2;
bool *scc_on_stack2;
static bool *visited;

static int strongconnect2(struct ptd_ph_vertex *vertex) {
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

        scc->internal_vertices_length = list.size();
        scc->internal_vertices = (struct ptd_ph_vertex **) calloc(
                scc->internal_vertices_length,
                sizeof(*(scc->internal_vertices))
        );

        for (size_t i = 0; i < scc->internal_vertices_length; ++i) {
            scc->internal_vertices[i] = list.at(i);
        }

        scc_components2->push_back(scc);
    }

    return 0;
}

struct ptd_ph_scc_graph *ptd_ph_find_strongly_connected_components(struct ptd_ph_graph *graph) {
    struct ptd_ph_scc_graph *scc_graph = (struct ptd_ph_scc_graph *) malloc(
            sizeof(*scc_graph)
    );

    scc_graph->graph = graph;

    scc_stack2 = new stack<struct ptd_ph_vertex *>;

    scc_index2 = 0;
    scc_indices2 = (size_t *) calloc(graph->vertices_length, sizeof(size_t));
    low_links2 = (size_t *) calloc(graph->vertices_length, sizeof(size_t));
    scc_on_stack2 = (bool *) calloc(graph->vertices_length, sizeof(bool));
    visited = (bool *) calloc(graph->vertices_length, sizeof(bool));
    scc_components2 = new vector<ptd_ph_scc_vertex *>();

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

        if (c->internal_vertices_length != 0) {
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

        if (scc->internal_vertices_length != 0) {
            scc_graph->vertices[index] = scc_components2->at(i);
            index++;
        } else {
            free(scc->internal_vertices);
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

        for (size_t j = 0; j < scc->internal_vertices_length; ++j) {
            struct ptd_ph_vertex *vertex = scc->internal_vertices[j];

            sccs_for_vertices[vertex->index] = scc;
        }
    }

    scc_graph->starting_vertex = sccs_for_vertices[graph->starting_vertex->index];

    for (size_t i = 0; i < scc_graph->vertices_length; ++i) {
        struct ptd_ph_scc_vertex *scc = scc_graph->vertices[i];
        set<struct ptd_ph_scc_vertex *> external_sccs;
        set<struct ptd_ph_vertex *> external_vertices;

        for (size_t j = 0; j < scc->internal_vertices_length; ++j) {
            struct ptd_ph_vertex *vertex = scc->internal_vertices[j];

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

        scc->external_vertices_length = external_vertices.size();
        scc->external_vertices = (struct ptd_ph_vertex **) calloc(
                scc->external_vertices_length,
                sizeof(*scc->external_vertices)
        );

        size_t set_index;

        set_index = 0;

        for (set<struct ptd_ph_scc_vertex *>::iterator itr = external_sccs.begin();
             itr != external_sccs.end();
             itr++) {
            scc->edges[set_index] = (struct ptd_ph_scc_edge *) malloc(sizeof(*(scc->edges[set_index])));
            scc->edges[set_index]->to = *itr;
            set_index++;
        }

        set_index = 0;

        for (set<struct ptd_ph_vertex *>::iterator itr = external_vertices.begin();
             itr != external_vertices.end();
             itr++) {
            scc->external_vertices[set_index] = *itr;
            set_index++;
        }
    }

    for (size_t i = 0; i < scc_graph->vertices_length; ++i) {
        struct ptd_ph_scc_vertex *scc = scc_graph->vertices[i];
        scc->index = i;
        size_t length = scc->internal_vertices_length + scc->external_vertices_length;
        size_t *old_indices = (size_t *) calloc(length, sizeof(*old_indices));

        for (size_t j = 0; j < scc->internal_vertices_length; ++j) {
            old_indices[j] = scc->internal_vertices[j]->index;
        }

        for (size_t j = 0; j < scc->external_vertices_length; ++j) {
            old_indices[j + scc->internal_vertices_length] = scc->external_vertices[j]->index;
        }

        for (size_t j = 0; j < scc->internal_vertices_length; ++j) {
            scc->internal_vertices[j]->index = j;
        }

        for (size_t j = 0; j < scc->external_vertices_length; ++j) {
            scc->external_vertices[j]->index = j + scc->internal_vertices_length;
        }

        void *matrix = matrix_init(length);
        // Top-left is internal->internal
        // Top-right is internal->external
        // Bottom-right is a diagonal of -1
        // Bottom-left is 0

        for (size_t j = 0; j < scc->internal_vertices_length; ++j) {
            struct ptd_ph_vertex *from = scc->internal_vertices[j];
            double rate = 0;

            for (size_t k = 0; k < from->edges_length; ++k) {
                rate += from->edges[k]->weight;
            }

            for (size_t k = 0; k < from->edges_length; ++k) {
                struct ptd_ph_edge *edge = from->edges[k];

                matrix_set(matrix, j, edge->to->index, edge->weight / rate);
            }

            matrix_set(matrix, j, j, -1);
        }

        for (size_t j = scc->internal_vertices_length; j < length; ++j) {
            matrix_set(matrix, j, j, -1);
        }

        void *inverted_matrix = matrix_invert(matrix, length);

        scc->internal_expected_visits = (double **) calloc(
                scc->internal_vertices_length,
                sizeof(*scc->internal_expected_visits)
        );

        for (size_t j = 0; j < scc->internal_vertices_length; ++j) {
            scc->internal_expected_visits[j] = (double *) calloc(
                    scc->internal_vertices_length,
                    sizeof(*scc->internal_expected_visits)
            );

            for (size_t k = 0; k < scc->internal_vertices_length; ++k) {
                scc->internal_expected_visits[j][k] = -matrix_get(inverted_matrix, j, k);
                if (scc->internal_expected_visits[j][k] < 0) {
                    fprintf(stderr, "ERROR! <0 (a)\n");


                }
            }
        }

        scc->external_expected_visits = (double **) calloc(
                scc->internal_vertices_length,
                sizeof(*scc->external_expected_visits)
        );

        for (size_t j = 0; j < scc->internal_vertices_length; ++j) {
            scc->external_expected_visits[j] = (double *) calloc(
                    scc->external_vertices_length,
                    sizeof(*scc->external_expected_visits)
            );

            for (size_t k = 0; k < scc->external_vertices_length; ++k) {
                scc->external_expected_visits[j][k] = -matrix_get(
                        inverted_matrix, j, k + scc->internal_vertices_length
                );

                if (scc->external_expected_visits[j][k] < 0) {
                    fprintf(stderr, "ERROR! <0 (b)\n");
                }
            }
        }

        for (size_t j = 0; j < scc->internal_vertices_length; ++j) {
            scc->internal_vertices[j]->index = old_indices[j];
        }

        for (size_t j = 0; j < scc->external_vertices_length; ++j) {
            scc->external_vertices[j]->index = old_indices[j + scc->internal_vertices_length];
        }

        matrix_destroy(inverted_matrix, length);
        matrix_destroy(matrix, length);
        free(old_indices);
    }

    delete scc_stack2;
    free(scc_indices2);
    free(low_links2);
    free(scc_on_stack2);
    free(visited);
    delete scc_components2;

    free(sccs_for_vertices);

    return scc_graph;
}

void ptd_ph_scc_graph_destroy(struct ptd_ph_scc_graph *scc_graph) {
    if (scc_graph == NULL) {
        return;
    }

    for (size_t i = 0; i < scc_graph->vertices_length; ++i) {
        struct ptd_ph_scc_vertex *scc = scc_graph->vertices[i];

        for (size_t j = 0; j < scc->edges_length; ++j) {
            free(scc->edges[j]);
        }

        for (size_t j = 0; j < scc->internal_vertices_length; ++j) {
            free(scc->internal_expected_visits[j]);
        }

        free(scc->internal_expected_visits);

        for (size_t j = 0; j < scc->internal_vertices_length; ++j) {
            free(scc->external_expected_visits[j]);
        }

        free(scc->external_expected_visits);

        free(scc->edges);
        free(scc->internal_vertices);
        free(scc->external_vertices);
        free(scc);
    }

    free(scc_graph->vertices);
    free(scc_graph);
}

ptd_phase_type_distribution_t *ptd_graph_as_phase_type_distribution(struct ptd_ph_graph *graph) {
    ptd_phase_type_distribution_t *res = (ptd_phase_type_distribution_t *) malloc(sizeof(*res));

    if (res == NULL) {
        return NULL;
    }

    res->length = 0;

    size_t size = graph->vertices_length;

    res->memory_allocated = size;
    res->vertices = (struct ptd_ph_vertex **) calloc(size, sizeof(struct ptd_ph_vertex *));

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

    size_t *indices = (size_t *) calloc(size, sizeof(*indices));
    size_t index = 0;

    for (size_t k = 0; k < graph->vertices_length; ++k) {
        struct ptd_ph_vertex *vertex = graph->vertices[k];

        if (graph->starting_vertex != vertex && vertex->edges_length != 0) {
            indices[vertex->index] = index;
            res->vertices[index] = vertex;
            index++;
        }
    }

    res->length = index;

    for (size_t k = 0; k < graph->vertices_length; ++k) {
        struct ptd_ph_vertex *vertex = graph->vertices[k];

        if (vertex->edges_length == 0) {
            continue;
        }

        if (vertex == graph->starting_vertex) {
            for (size_t i = 0; i < vertex->edges_length; ++i) {
                struct ptd_ph_edge *edge = vertex->edges[i];

                if (edge->to->edges_length != 0) {
                    res->initial_probability_vector[indices[edge->to->index]] = edge->weight;
                }
            }

            continue;
        }

        for (size_t i = 0; i < vertex->edges_length; ++i) {
            struct ptd_ph_edge *edge = vertex->edges[i];

            if (edge->to->edges_length != 0) {
                res->sub_intensity_matrix[indices[vertex->index]][indices[edge->to->index]] += edge->weight;
            }

            res->sub_intensity_matrix[indices[vertex->index]][indices[vertex->index]] -= edge->weight;
        }
    }

    free(indices);

    return res;
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

int ptd_ph_vertex_to_s(struct ptd_ph_vertex *vertex, char *buffer, size_t buffer_length) {
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

/*
 * Models
 */


/*static ptd_avl_tree_t *avl_tree = NULL;
static struct ptd_ph_graph *kingman_graph = NULL;

static int make_kingman(struct ptd_ph_vertex *vertex) {
    int *state = (int *) calloc(kingman_graph->state_length, sizeof(int));
    memcpy(state, vertex->state, kingman_graph->state_length * sizeof(int));

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

            struct ptd_ph_vertex *child = ptd_avl_tree_vertex_find(avl_tree, state);

            if (child == NULL) {
                int *child_state = (int *) calloc(kingman_graph->state_length, sizeof(int));

                if (child_state == NULL) {
                    return -1;
                }

                memcpy(child_state, state, kingman_graph->state_length * sizeof(int));

                child = ptd_ph_vertex_create_state(kingman_graph, child_state);

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
}*/

/*struct ptd_ph_graph *ptd_model_kingman(size_t n) {
    avl_tree = ptd_avl_tree_create(n);

    if (avl_tree == NULL) {
        return NULL;
    }

    kingman_graph = ptd_ph_graph_create(n);

    if (kingman_graph == NULL) {
        ptd_avl_tree_vertex_destroy(avl_tree);
        return NULL;
    }

    struct ptd_ph_vertex *initial = ptd_ph_vertex_create(kingman_graph);

    if (initial == NULL) {
        ptd_ph_graph_destroy(kingman_graph);
        ptd_avl_tree_vertex_destroy(avl_tree);
    }

    initial->state[0] = (int)n;

    if (ptd_ph_graph_add_edge(kingman_graph->starting_vertex, initial, 1) != 0) {
        ptd_ph_graph_destroy(kingman_graph);
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
}*/

/*
 * State has two islands and two loci. For two loci we can have the states
 * A-A, A-, -A, for the state of recombination. Since we have two islands,
 * this gives us 6 states.
 *
 * State 0,1,2: island 1
 * State 3,4,5: island 2
 */
/*
static struct ptd_graph *two_island_two_loci_recomb_graph;

static double recomb_rate;
static size_t N[2];
static double mig[2];

#define IDX_COMB 0
#define IDX_L 1
#define IDX_R 2

static struct ptd_ph_vertex *add_child(int *state) {
    struct ptd_ph_vertex *child = ptd_avl_tree_vertex_find(avl_tree, state);

    if (child == NULL) {
        int *child_state = (int *) calloc(
                two_island_two_loci_recomb_graph->state_length,
                sizeof(int)
        );

        if (child_state == NULL) {
            return NULL;
        }

        memcpy(
                child_state, state,
                two_island_two_loci_recomb_graph->state_length * sizeof(int)
        );

        child = ptd_ph_vertex_create_state(
                two_island_two_loci_recomb_graph, child_state
        );

        if (ptd_avl_tree_vertex_insert(avl_tree, child_state, child)) {
            return NULL;
        }
    }

    return child;
}

static int visit_two_island_two_loci_recomb(struct ptd_ph_vertex *vertex) {
    // If there is only one left
    if (vertex->state[0] + vertex->state[1] +
        vertex->state[2] + vertex->state[3] +
        vertex->state[4] + vertex->state[5] <= 1) {
        return 0;
    }

    int *state = (int *) calloc(
            two_island_two_loci_recomb_graph->state_length,
            sizeof(int)
    );

    memcpy(
            state, vertex->state,
            two_island_two_loci_recomb_graph->state_length * sizeof(int)
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

                    struct ptd_ph_vertex *child = add_child(state);

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

                struct ptd_ph_vertex *child = add_child(state);

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

                struct ptd_ph_vertex *child = add_child(state);

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
                struct ptd_ph_vertex *child = add_child(state);
                state[island * 3 + IDX_COMB]++;

                ptd_add_edge(vertex, child, combinations * N[island]);
            }

            // Left and combined
            if (state[island * 3 + IDX_COMB] > 0 && state[island * 3 + IDX_L] > 0) {
                size_t combinations = state[island * 3 + IDX_COMB] * state[island * 3 + IDX_L];

                state[island * 3 + IDX_L]--;

                struct ptd_ph_vertex *child = add_child(state);

                state[island * 3 + IDX_L]++;

                ptd_add_edge(vertex, child, combinations * N[island]);
            }

            // Right and combined
            if (state[island * 3 + IDX_COMB] > 0 && state[island * 3 + IDX_R] > 0) {
                size_t combinations = state[island * 3 + IDX_COMB] * state[island * 3 + IDX_R];

                state[island * 3 + IDX_R]--;
                struct ptd_ph_vertex *child = add_child(state);
                state[island * 3 + IDX_R]++;

                ptd_add_edge(vertex, child, combinations * N[island]);
            }

            // Left and left
            if (state[island * 3 + IDX_L] > 1) {
                size_t combinations = state[island * 3 + IDX_L] * (state[island * 3 + IDX_L] - 1) / 2;

                state[island * 3 + IDX_L]--;
                struct ptd_ph_vertex *child = add_child(state);
                state[island * 3 + IDX_L]++;

                ptd_add_edge(vertex, child, combinations * N[island]);
            }

            // Right and right
            if (state[island * 3 + IDX_R] > 1) {
                size_t combinations = state[island * 3 + IDX_R] * (state[island * 3 + IDX_R] - 1) / 2;

                state[island * 3 + IDX_R]--;
                struct ptd_ph_vertex *child = add_child(state);
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

    struct ptd_ph_vertex *initial = ptd_ph_vertex_create(two_island_two_loci_recomb_graph);
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
 */

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
    graph->scc_graph = NULL;
    graph->starting_vertex = ptd_ph_vertex_create(graph);

    return graph;
}

void ptd_ph_graph_destroy(struct ptd_ph_graph *graph) {
    for (size_t i = 0; i < graph->vertices_length; ++i) {
        ptd_ph_vertex_destroy(graph->vertices[i]);
    }

    free(graph->vertices);
    ptd_ph_scc_graph_destroy(graph->scc_graph);
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


    ptd_ph_scc_graph_destroy(graph->scc_graph);
    graph->scc_graph = NULL;

    return vertex;
}

double ptd_ph_vertex_rate(struct ptd_ph_vertex *vertex) {
    double rate = 0;

    for (size_t i = 0; i < vertex->edges_length; ++i) {
        rate += vertex->edges[i]->weight;
    }

    return rate;
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

    ptd_ph_scc_graph_destroy(from->graph->scc_graph);
    from->graph->scc_graph = NULL;

    return edge;
}

static struct ptd_directed_vertex **ptd_graph_topological_sort(struct ptd_directed_graph *graph) {
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
    queue<struct ptd_directed_vertex *> q;
    q.push(graph->vertices[0]);
    size_t topo_index = 0;

    while (!q.empty()) {
        struct ptd_directed_vertex *vertex = q.front();
        q.pop();

        res[topo_index] = vertex;
        visited[vertex->index] = true;
        topo_index++;

        for (size_t i = 0; i < vertex->edges_length; ++i) {
            struct ptd_directed_vertex *child = vertex->edges[i]->to;

            nparents[child->index]--;

            if (nparents[child->index] == 0 && !visited[child->index]) {
                visited[child->index] = true;
                q.push(child);
            }
        }

        if (q.empty() && !has_pushed_all_others) {
            for (size_t i = 0; i < graph->vertices_length; ++i) {
                struct ptd_directed_vertex *independent_vertex = graph->vertices[i];

                if (nparents[independent_vertex->index] == 0 && !visited[independent_vertex->index]) {
                    q.push(independent_vertex);
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

            return NULL;
        }
    }

    free(nparents);
    free(visited);

    return res;
}

bool ptd_ph_graph_is_acyclic(struct ptd_ph_graph *graph) {
    struct ptd_ph_vertex **sorted = ptd_ph_graph_topological_sort(graph);

    bool is_acyclic = (sorted != NULL);

    free(sorted);

    return is_acyclic;
}

struct ptd_ph_vertex **ptd_ph_graph_topological_sort(struct ptd_ph_graph *graph) {
    return (struct ptd_ph_vertex **) ptd_graph_topological_sort((struct ptd_directed_graph *) graph);
}

struct ptd_ph_scc_vertex **ptd_ph_scc_graph_topological_sort(struct ptd_ph_scc_graph *graph) {
    return (struct ptd_ph_scc_vertex **) ptd_graph_topological_sort((struct ptd_directed_graph *) graph);
}

double *ptd_ph_graph_acyclic_visit_probability(struct ptd_ph_graph *graph) {
    struct ptd_ph_vertex **topo = ptd_ph_graph_topological_sort(graph);

    if (topo == NULL) {
        return NULL;
    }

    double *res = (double *) calloc(graph->vertices_length, sizeof(*res));
    res[graph->starting_vertex->index] = 1.0f;

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        struct ptd_ph_vertex *vertex = topo[i];
        double prob = res[vertex->index];

        double rate = ptd_ph_vertex_rate(vertex);

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            struct ptd_ph_edge *edge = vertex->edges[j];
            res[edge->to->index] += edge->weight / rate * prob;
        }
    }

    free(topo);

    return res;
}

double *ptd_ph_graph_cyclic_entry_probability(struct ptd_ph_scc_graph *scc_graph) {
    struct ptd_ph_scc_vertex **topo = ptd_ph_scc_graph_topological_sort(scc_graph);

    if (topo == NULL) {
        return NULL;
    }

    double *res = (double *) calloc(scc_graph->graph->vertices_length, sizeof(*res));
    res[scc_graph->graph->starting_vertex->index] = 1.0f;

    for (size_t i = 0; i < scc_graph->vertices_length; ++i) {
        struct ptd_ph_scc_vertex *scc_vertex = topo[i];

        for (size_t k = 0; k < scc_vertex->internal_vertices_length; ++k) {
            struct ptd_ph_vertex *vertex = scc_vertex->internal_vertices[k];
            double prob = res[vertex->index];

            for (size_t j = 0; j < scc_vertex->external_vertices_length; ++j) {
                struct ptd_ph_vertex *to = scc_vertex->external_vertices[j];
                res[to->index] += scc_vertex->external_expected_visits[k][j] * prob;
            }
        }
    }

    free(topo);

    return res;
}

double *ptd_ph_scc_graph_entry_probability(struct ptd_ph_scc_graph *scc_graph) {
    double *res = (double *) calloc(scc_graph->vertices_length, sizeof(*res));
    double *vertex_probs = ptd_ph_graph_cyclic_entry_probability(scc_graph);

    for (size_t i = 0; i < scc_graph->vertices_length; ++i) {
        struct ptd_ph_scc_vertex *scc_vertex = scc_graph->vertices[i];

        for (size_t j = 0; j < scc_vertex->internal_vertices_length; ++j) {
            struct ptd_ph_vertex *vertex = scc_vertex->internal_vertices[j];

            res[i] += vertex_probs[vertex->index];
        }
    }

    free(vertex_probs);

    return res;
}

double *ptd_ph_graph_cyclic_expected_visits(struct ptd_ph_scc_graph *scc_graph) {
    double *res = (double *) calloc(scc_graph->graph->vertices_length, sizeof(*res));
    double *vertex_probs = ptd_ph_graph_cyclic_entry_probability(scc_graph);

    for (size_t i = 0; i < scc_graph->vertices_length; ++i) {
        struct ptd_ph_scc_vertex *scc_vertex = scc_graph->vertices[i];

        for (size_t j = 0; j < scc_vertex->internal_vertices_length; ++j) {
            struct ptd_ph_vertex *vertex = scc_vertex->internal_vertices[j];

            if (vertex_probs[vertex->index] == 0) {
                continue;
            }

            for (size_t k = 0; k < scc_vertex->internal_vertices_length; ++k) {
                struct ptd_ph_vertex *to = scc_vertex->internal_vertices[k];
                double visits = scc_vertex->internal_expected_visits[j][k];

                res[to->index] += visits * vertex_probs[vertex->index];
            }
        }
    }

    free(vertex_probs);

    return res;
}

double *ptd_ph_graph_acyclic_moment_rewards(struct ptd_ph_graph *graph, double *rewards) {
    struct ptd_ph_vertex **topo = ptd_ph_graph_topological_sort(graph);

    if (topo == NULL) {
        return NULL;
    }

    double *res = (double *) calloc(graph->vertices_length, sizeof(*res));
    double *rates = (double *) calloc(graph->vertices_length, sizeof(*rates));

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        rates[i] = ptd_ph_vertex_rate(graph->vertices[i]);
    }

    for (size_t ii = 0; ii < graph->vertices_length; ++ii) {
        size_t i = graph->vertices_length - ii - 1;

        struct ptd_ph_vertex *vertex = topo[i];

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            struct ptd_ph_vertex *to = vertex->edges[j]->to;
            double weight = vertex->edges[j]->weight;

            res[vertex->index] += weight / rates[vertex->index] * res[to->index];
        }

        if (rates[vertex->index] != 0) {
            res[vertex->index] += rewards[vertex->index] / rates[vertex->index];
        }
    }

    free(rates);
    free(topo);

    return res;
}

double *ptd_ph_graph_cyclic_moment_rewards(struct ptd_ph_scc_graph *scc_graph, double *rewards) {
    struct ptd_ph_scc_vertex **topo = ptd_ph_scc_graph_topological_sort(scc_graph);

    if (topo == NULL) {
        return NULL;
    }

    double *res = (double *) calloc(scc_graph->graph->vertices_length, sizeof(*res));
    double *rates = (double *) calloc(scc_graph->graph->vertices_length, sizeof(*rates));

    for (size_t i = 0; i < scc_graph->graph->vertices_length; ++i) {
        rates[i] = ptd_ph_vertex_rate(scc_graph->graph->vertices[i]);
    }

    for (size_t ii = 0; ii < scc_graph->vertices_length; ++ii) {
        size_t i = scc_graph->vertices_length - ii - 1;

        struct ptd_ph_scc_vertex *scc_vertex = topo[i];

        for (size_t j = 0; j < scc_vertex->internal_vertices_length; ++j) {
            struct ptd_ph_vertex *vertex = scc_vertex->internal_vertices[j];

            for (size_t k = 0; k < scc_vertex->internal_vertices_length; ++k) {
                struct ptd_ph_vertex *to = scc_vertex->internal_vertices[k];

                double visits = scc_vertex->internal_expected_visits[j][k];

                if (rates[to->index] != 0) {
                    res[vertex->index] += visits * rewards[to->index] / rates[to->index];
                }
            }

            for (size_t k = 0; k < scc_vertex->external_vertices_length; ++k) {
                struct ptd_ph_vertex *to = scc_vertex->external_vertices[k];
                double visits = scc_vertex->external_expected_visits[j][k];

                res[vertex->index] += visits * res[to->index];
            }
        }
    }

    free(rates);
    free(topo);

    return res;
}

struct edge_cmp {
    bool operator()(ptd_ph_edge *a, ptd_ph_edge *b) {
        return a->to < b->to;
    }
};

#define REWARD_EPSILON 0.000001

int ptd_ph_graph_reward_transform(struct ptd_ph_graph *graph, double *rewards) {
    vector<set<struct ptd_ph_vertex *> > parents(graph->vertices_length);
    double *rates = (double *) calloc(graph->vertices_length, sizeof(*rates));

    vector<set<struct ptd_ph_edge *, struct edge_cmp> > edges(graph->vertices_length);

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        struct ptd_ph_vertex *vertex = graph->vertices[i];

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            parents[vertex->edges[j]->to->index].insert(vertex);

            if (rewards[i] > REWARD_EPSILON && vertex != graph->starting_vertex && vertex->edges_length != 0) {
                vertex->edges[j]->weight /= rewards[i];
            }

            edges[i].insert(vertex->edges[j]);
        }

        rates[i] = ptd_ph_vertex_rate(vertex);
    }

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        struct ptd_ph_vertex *vertex = graph->vertices[i];

        if (graph->starting_vertex == vertex || rates[vertex->index] == 0) {
            continue;
        }

        if (rewards[i] <= REWARD_EPSILON) {
            vector<struct ptd_ph_vertex *> ps(parents[i].begin(), parents[i].end());

            for (size_t k = 0; k < ps.size(); ++k) {
                struct ptd_ph_vertex *parent = ps[k];

                for (set<struct ptd_ph_edge *, struct edge_cmp>::iterator it =
                        edges[vertex->index].begin();
                     it != edges[vertex->index].end();
                     ++it) {
                    struct ptd_ph_vertex *child = (*it)->to;
                    double weight = (*it)->weight;

                    if (child == parent) {
                        continue;
                    }

                    struct ptd_ph_edge fake_edge;
                    fake_edge.to = vertex;
                    struct ptd_ph_edge *edge_p_to_me = *(edges[parent->index].find(&fake_edge));

                    double prob = weight / rates[vertex->index];

                    set<struct ptd_ph_edge *, struct edge_cmp>::iterator it_edge;
                    it_edge = edges[parent->index].find(*it);

                    if (it_edge != edges[parent->index].end()) {
                        (*it_edge)->weight += prob * edge_p_to_me->weight;
                    } else {
                        struct ptd_ph_edge *new_edge = (struct ptd_ph_edge *) malloc(sizeof(*new_edge));
                        new_edge->weight = prob * edge_p_to_me->weight;
                        new_edge->to = child;
                        edges[parent->index].insert(new_edge);
                        parents[child->index].insert(parent);
                    }
                }
            }

            for (size_t k = 0; k < ps.size(); ++k) {
                struct ptd_ph_vertex *parent = ps[k];

                struct ptd_ph_edge fake_edge;
                fake_edge.to = vertex;

                set<struct ptd_ph_edge *, struct edge_cmp>::iterator it_edge;
                it_edge = edges[parent->index].find(&fake_edge);
                struct ptd_ph_edge *edge = *it_edge;
                edges[parent->index].erase(it_edge);
                free(edge);
            }

            for (set<struct ptd_ph_edge *, struct edge_cmp>::iterator it =
                    edges[vertex->index].begin();
                 it != edges[vertex->index].end();
                 ++it) {
                struct ptd_ph_vertex *child = (*it)->to;

                parents[child->index].erase(vertex);
                free(*it);
            }

            vertex->edges_length = 0;
            ptd_ph_vertex_destroy(vertex);
        }
    }

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        if (rates[i] != 0 &&
            graph->starting_vertex->index != i &&
            rewards[i] <= REWARD_EPSILON) {
            continue;
        }

        struct ptd_ph_vertex *vertex = graph->vertices[i];

        vertex->edges_length = 0;

        for (set<struct ptd_ph_edge *, struct edge_cmp>::iterator it =
                edges[vertex->index].begin();
             it != edges[vertex->index].end();
             ++it) {
            ptd_ph_graph_add_edge(vertex, (*it)->to, (*it)->weight);
        }

        for (set<struct ptd_ph_edge *, struct edge_cmp>::iterator it =
                edges[vertex->index].begin();
             it != edges[vertex->index].end();
             ++it) {
            free(*it);
        }
    }

    size_t old_length = graph->vertices_length;
    struct ptd_ph_vertex **old_vertices = graph->vertices;

    graph->vertices_length = 1;
    graph->vertices[0] = graph->starting_vertex;
    graph->starting_vertex->index = 0;

    for (size_t i = 1; i < old_length; ++i) {
        if (rewards[i] <= REWARD_EPSILON && rates[i] != 0) {
            continue;
        }

        struct ptd_ph_vertex *vertex = old_vertices[i];

        graph->vertices[graph->vertices_length] = vertex;
        vertex->index = graph->vertices_length;
        graph->vertices_length++;
    }

    free(rates);

    ptd_ph_scc_graph_destroy(graph->scc_graph);
    graph->scc_graph = NULL;

    return 0;
}

double *ptd_ph_graph_expected_visits(struct ptd_ph_graph *graph) {
    if (ptd_ph_graph_is_acyclic(graph)) {
        return ptd_ph_graph_acyclic_visit_probability(graph);
    } else {
        ptd_ph_precompute_strongly_connected_components(graph);

        return ptd_ph_graph_cyclic_expected_visits(graph->scc_graph);
    }
}

double *ptd_ph_graph_expected_waiting_time(struct ptd_ph_graph *graph) {
    double *res = ptd_ph_graph_expected_visits(graph);

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        struct ptd_ph_vertex *vertex = graph->vertices[i];
        double rate = ptd_ph_vertex_rate(vertex);

        if (rate != 0 && vertex != graph->starting_vertex) {
            res[i] /= ptd_ph_vertex_rate(vertex);
        } else {
            res[i] = 0;
        }
    }

    return res;
}

double *ptd_ph_graph_moment_rewards(struct ptd_ph_graph *graph, double *rewards) {
    if (ptd_ph_graph_is_acyclic(graph)) {
        return ptd_ph_graph_acyclic_moment_rewards(graph, rewards);
    } else {
        ptd_ph_precompute_strongly_connected_components(graph);

        return ptd_ph_graph_cyclic_moment_rewards(graph->scc_graph, rewards);
    }
}