#include <iostream>
#include <vector>
#include <cstring>
#include <stack>
#include "phase.h"

using namespace std;

void _reset_graph_visited(vertex_t *node, size_t reset_int) {
    if (node->reset_int == reset_int) {
        return;
    }

    node->reset_int = reset_int;


    for (vertex_t *child : node->children) {
        _reset_graph_visited(child, reset_int);
    }

    node->visited = false;
}

static void reset_graph_visited(vertex_t *node) {
    _reset_graph_visited(node, node->reset_int + 1);
}

typedef struct avl_vec_node {
    struct avl_vec_node *left;
    struct avl_vec_node *right;
    struct avl_vec_node *parent;
    signed short balance;
    vec_entry_t *key;
    vertex_t *entry;
} avl_vec_node_t;

static inline int radix_cmp(const vec_entry_t* a, const vec_entry_t* b,
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
avl_vec_node_t *rotate_left_right(avl_vec_node_t *parent, avl_vec_node_t *child) {
    avl_vec_node_t *child_right_left, *child_right_right;
    avl_vec_node_t *child_right = child->right;
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
avl_vec_node_t *rotate_right_left(avl_vec_node_t *parent, avl_vec_node_t *child) {
    avl_vec_node_t *child_left_right, *child_left_left;
    avl_vec_node_t *child_left = child->left;

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
avl_vec_node_t *rotate_left(avl_vec_node_t *parent, avl_vec_node_t *child) {
    avl_vec_node_t *child_left;

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
avl_vec_node_t *rotate_right(avl_vec_node_t *parent, avl_vec_node_t *child) {
    avl_vec_node_t *child_right;

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

int avl_vec_node_create(avl_vec_node_t **node, vec_entry_t *key, vertex_t *entry, avl_vec_node_t *parent) {
    if ((*node = (avl_vec_node_t*) malloc(sizeof(avl_vec_node_t))) == nullptr) {
        return 1;
    }

    (*node)->key = key;
    (*node)->entry = entry;
    (*node)->left = nullptr;
    (*node)->right = nullptr;
    (*node)->parent = parent;
    (*node)->balance = 0;

    return 0;
}

void avl_vec_node_destroy(avl_vec_node_t *node) {
    if (node == nullptr) {
        return;
    }

    avl_vec_node_destroy(node->left);
    avl_vec_node_destroy(node->right);

    free(node);
}

const avl_vec_node_t * avl_vec_find(const avl_vec_node_t *rootptr, const vec_entry_t *key, const size_t vec_length) {
    if (rootptr == nullptr) {
        return nullptr;
    }

    const avl_vec_node_t *node = rootptr;

    while (true) {
        int res = radix_cmp(key, node->key, vec_length);
        if (res < 0) {
            if (node->left == nullptr) {
                return nullptr;
            } else {
                node = node->left;
            }
        } else if (res > 0) {
            if (node->right == nullptr) {
                return nullptr;
            } else {
                node = node->right;
            }
        } else {
            return node;
        }
    }
}

int find_or_insert_vec(avl_vec_node_t **out, avl_vec_node_t *rootptr, vec_entry_t *key, vertex_t *entry, const size_t vec_length) {
    if (avl_vec_node_create(out, key, entry, nullptr)) {
        return -1;
    }

    if (rootptr == nullptr) {
        return 1;
    }

    avl_vec_node_t *node = rootptr;

    while(true) {
        int res = radix_cmp(key, node->key, vec_length);
        if (res < 0) {
            if (node->left == nullptr) {
                node->left = *out;
                break;
            } else {
                node = node->left;
            }
        } else if (res > 0) {
            if (node->right == nullptr) {
                node->right = *out;
                break;
            } else {
                node = node->right;
            }
        } else {
            *out = node;
        }
    }

    (*out)->parent = node;

    return 0;
}

int avl_vec_insert(avl_vec_node_t **root, vec_entry_t *key, vertex_t *entry, const size_t vec_length) {
    avl_vec_node_t *child;

    if (*root == nullptr) {
        if (avl_vec_node_create(root, key, entry, nullptr)) {
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

    avl_vec_node_t *pivot, *rotated_parent;

    for (avl_vec_node_t *parent = child->parent; parent != nullptr; parent = child->parent) {
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

static size_t avl_vec_get_size(avl_vec_node_t *node) {
    if (node == nullptr) {
        return 0;
    }

    return 1 + avl_vec_get_size(node->left) + avl_vec_get_size(node->right);
}

static void add_edge(vertex_t *from, vertex_t *to, double weight) {
    from->children.push_back(to);
    from->weights.push_back(weight);
    from->rate += weight;
    to->parents.push_back(from);
    to->weights_parent.push_back(weight);
}

static void avl_free(avl_vec_node_t *node) {
    if (node->left != nullptr) {
        avl_free(node->left);
    }

    if (node->right != nullptr) {
        avl_free(node->right);
    }

    free(node);
}

static void _get_abs_vertex(vertex_t **abs_vertex, vertex_t *graph) {
    if (graph->visited) {
        return;
    }

    graph->visited = true;

    if (graph->children.empty()) {
        if (*abs_vertex == nullptr || *abs_vertex == graph) {
            *abs_vertex = graph;
        } else {
            DIE_ERROR(1, "Found multiple absorbing vertices\n");
        }
    }

    for (vertex_t *child : graph->children) {
        _get_abs_vertex(abs_vertex, child);
    }
}

static vertex_t *get_abs_vertex(vertex_t *graph) {
    reset_graph_visited(graph);
    vertex_t *abs_vertex = nullptr;
    _get_abs_vertex(&abs_vertex, graph);

    return abs_vertex;
}

static int kingman_visit_vertex(vertex_t **out,
                                vec_entry_t *state,
                                avl_vec_node_t *bst,
                                size_t state_size,
                                size_t vector_length,
                                const size_t vec_nmemb,
                                size_t m) {
    avl_vec_node_t *bst_node = (avl_vec_node_t*)avl_vec_find(bst, state, vector_length);

    if (bst_node != nullptr) {
        *out = bst_node->entry;
        return 0;
    } else {
        vec_entry_t *vertex_state = (vec_entry_t*) malloc(sizeof(vec_entry_t) * vector_length);
        memcpy(vertex_state, state, sizeof(vec_entry_t) * vector_length);

        *out = new vertex_t(vertex_state, state_size);

        avl_vec_insert(&bst, vertex_state, *out, vector_length);

        vec_entry_t *v = state;

        for (vec_entry_t i = 0; i < state_size; i++) {
            for (vec_entry_t j = i; j < state_size; j++) {
                if (((i == j && v[i] >= 2) || (i != j && v[i] > 0 && v[j] > 0))) {
                    double t = i == j ? v[i] * (v[i] - 1) / 2 : v[i] * v[j];

                    v[i]--;
                    v[j]--;

                    const size_t inc_pos = min((i + j + 2) - 1, m+1);

                    v[inc_pos]++;

                    vertex_t *new_vertex;
                    kingman_visit_vertex(&new_vertex,
                                         v, bst,
                                         state_size,
                                         vector_length,
                                         vec_nmemb,
                                         m);

                    v[i]++;
                    v[j]++;
                    v[inc_pos]--;

                    add_edge(*out, new_vertex, t);
                }
            }
        }

        return 0;
    }
}

int gen_kingman_graph(vertex_t **graph, size_t n, size_t m) {
    size_t state_size = n;

    vec_entry_t *initial = (vec_entry_t*)calloc(state_size, sizeof(vec_entry_t));
    initial[0] = state_size;

    vec_entry_t *mrca = (vec_entry_t*)calloc(state_size, sizeof(vec_entry_t));
    mrca[state_size-1] = 1;

    vertex_t *absorbing_vertex = new vertex_t(mrca, state_size);

    avl_vec_node_t *BST;
    avl_vec_node_create(&BST, mrca, absorbing_vertex, nullptr);

    vertex_t *state_graph;

    kingman_visit_vertex(&state_graph, initial, BST, state_size, state_size,
                         state_size, m);

    vec_entry_t *start_state = (vec_entry_t*)calloc(state_size, sizeof(vec_entry_t));

    vertex_t *start = new vertex_t(start_state, state_size);

    add_edge(start, state_graph, 1);
    avl_free(BST);
    free(initial);
    *graph = start;

    return 0;
}


static void print_vector_spacing(FILE *stream,vec_entry_t *v, size_t nmemb, size_t spacing) {
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


static void _print_graph_list(FILE *stream, vertex_t *node,
                              bool indexed,
                              size_t vec_length, size_t vec_spacing) {
    if (node->visited) {
        return;
    }

    node->visited = true;

    fprintf(stream, "Node: ");
    print_vector_spacing(stream, node->state,
                         vec_length, vec_spacing);
    if (indexed) {
        fprintf(stream, " (%zu)", node->vertex_index);
    }
    fprintf(stream, ":\n");


    for (size_t i = 0; i < node->children.size(); ++i) {
        fprintf(stream, "\t");
        fprintf(stream, "(%f) ", node->weights[i]);
        print_vector_spacing(stream, node->children[i]->state,
                             vec_length, vec_spacing);

        fprintf(stream, "e %f d %f ", node->exp[0], node->desc[0]);

        if (indexed) {
            fprintf(stream, " (%zu)", node->vertex_index);
        }

        fprintf(stream, "\n");
    }

    fprintf(stream, "\n");
    for (size_t i = 0; i < node->children.size(); i++) {
        _print_graph_list(stream, node->children[i],
                          indexed,
                          vec_length, vec_spacing);
    }
}

void print_graph_list(FILE *stream, vertex_t *graph,
                      bool indexed,
                      size_t vec_length, size_t vec_spacing) {
    reset_graph_visited(graph);
    _print_graph_list(stream, graph, indexed, vec_length, vec_spacing);
    fflush(stream);
}

void mph_cov_assign_vertex_all(vertex_t *node, size_t m) {
    if (node->visited) {
        return;
    }

    node->visited = true;

    if (node->parents.empty()) {
        // Starting vertex
        node->prob = 1.0f;
    } else {
        node->prob = 0.0f;
    }

    for (size_t i = 0; i < node->parents.size(); i++) {
        vertex_t *parent = node->parents[i];

        mph_cov_assign_vertex_all(parent, m);

        node->prob += node->weights_parent[i]/parent->rate * parent->prob;
    }


    node->exp = vector<double>();

    for (size_t j = 0; j < m; ++j) {
        if (node->rate != 0) {
            node->exp.push_back(node->prob * node->state[j] / node->rate);
        } else {
            node->exp.push_back(0);
        }
    }
}

void mph_cov_assign_desc_all(vertex_t *node, size_t m) {
    if (node->visited) {
        return;
    }

    node->visited = true;

    node->desc = vector<double>(m);

    for (vertex_t *child : node->children) {
        mph_cov_assign_desc_all(child, m);
    }

    for (size_t i = 0; i < node->children.size(); i++) {
        for (size_t j = 0; j < m; ++j) {
            node->desc[j] += node->weights[i] / node->rate * node->children[i]->desc[j];
        }
    }

    vector<double> exp(m);

    for (size_t j = 0; j < m; ++j) {
        if (node->rate != 0) {
            exp[j] = node->state[j] / node->rate;
        } else {
            exp[j] = 0;
        }
    }

    for (size_t j = 0; j < m; ++j) {
        node->desc[j] += exp[j];
    }
}

double **cov;
double *expectation;

void _mph_cov_all(vertex_t *node, size_t m) {
    if (node->visited) {
        return;
    }

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            cov[i][j] += node->desc[i] * node->exp[j];
            cov[i][j] += node->desc[j] * node->exp[i];
        }
    }

    node->visited = true;

    for (vertex_t *child : node->children) {
        _mph_cov_all(child, m);
    }
}

cov_exp_return mph_cov_exp_all(vertex_t *graph, size_t m) {
    vertex_t *abs = get_abs_vertex(graph);
    reset_graph_visited(graph);
    mph_cov_assign_vertex_all(abs, m);
    reset_graph_visited(graph);
    mph_cov_assign_desc_all(graph, m);

    cov = (double**)calloc(m, sizeof(double*));

    for (size_t i = 0; i < m; ++i) {
        cov[i] = (double*)calloc(m, sizeof(double));
    }

    reset_graph_visited(graph);
    _mph_cov_all(graph, m);

    expectation = (double*) calloc(m, sizeof(double));

    for (size_t i = 0; i < m; ++i) {
        expectation[i] = graph->desc[i];

        for (size_t j = 0; j <= i; ++j) {
            cov[i][j] -= graph->desc[i] *
                         graph->desc[j];
        }
    }

    return (cov_exp_return) {.cov = &cov, .exp = &expectation};
}

vector<vertex_t *> *vertex_to_free;

static void mark_for_deletion(vertex_t *vertex) {
    if (vertex->visited) {
        return;
    }

    vertex->visited = true;

    vertex_to_free->push_back(vertex);

    for (auto child : vertex->children) {
        mark_for_deletion(child);
    }
}

void graph_free(vertex_t *graph) {
    vertex_to_free = new vector<vertex_t *>();
    reset_graph_visited(graph);
    mark_for_deletion(graph);

    for (auto vertex : *vertex_to_free) {
        delete vertex;
    }

    delete vertex_to_free;
}

