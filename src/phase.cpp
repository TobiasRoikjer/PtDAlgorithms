#include <iostream>
#include <vector>
#include <cstring>
#include <stack>
#include <queue>
#include <algorithm>
#include "phase.h"

using namespace std;

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

void _reset_graph_visited(vertex_t *vertex, size_t reset_int) {
    if (vertex->reset_int == reset_int) {
        return;
    }

    vertex->reset_int = reset_int;

    for (auto child : vertex->children) {
        _reset_graph_visited(child.vertex, reset_int);
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
    if ((*vertex = (avl_vec_vertex_t*) malloc(sizeof(avl_vec_vertex_t))) == nullptr) {
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

const avl_vec_vertex_t * avl_vec_find(const avl_vec_vertex_t *rootptr, const vec_entry_t *key, const size_t vec_length) {
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

int find_or_insert_vec(avl_vec_vertex_t **out, avl_vec_vertex_t *rootptr, vec_entry_t *key, vertex_t *entry, const size_t vec_length) {
    if (avl_vec_vertex_create(out, key, entry, nullptr)) {
        return -1;
    }

    if (rootptr == nullptr) {
        return 1;
    }

    avl_vec_vertex_t *vertex = rootptr;

    while(true) {
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

static void remove_edge(vertex_t *from, vertex_t *to){
    if (from == to) {
        return;
    }

    auto elem = lower_bound(from->children.begin(),
                     from->children.end(), to);

    auto elem_parent = lower_bound(to->parents.begin(),
                     to->parents.end(), from);

    from->rate -= elem->weight;
    from->children.erase(elem);
    to->parents.erase(elem_parent);
}

static void add_edge_unsorted(vertex_t *from, vertex_t *to, double weight) {
    if (from == to) {
        return;
    }

    from->children.push_back({.vertex = to, .weight = weight});
    from->rate += weight;
    to->parents.push_back({.vertex = from, .weight = weight});
}

static void add_edge(vertex_t *from, vertex_t *to, double weight) {
    if (from == to) {
        return;
    }

    from->children.insert(lower_bound(
            from->children.begin(), from->children.end(), to),
                    {.vertex = to, .weight = weight});
    from->rate += weight;
    to->parents.insert(lower_bound(
            to->parents.begin(), to->parents.end(), from),
                    {.vertex = from, .weight = weight});
}

static queue<vertex_t *> enqueue_vertices(vertex_t *graph) {
    queue<vertex_t *> ret;
    queue<vertex_t *> queue;
    reset_graph_visited(graph);

    queue.push(graph);

    while(!queue.empty()) {
        vertex_t *vertex = queue.front();
        queue.pop();

        if (vertex->visited) {
            continue;
        }

        vertex->visited = true;
        ret.push(vertex);

        for (auto child : vertex->children) {
            queue.push(child.vertex);
        }
    }

    return ret;
}

static void increase_weight(vertex_t *from, vertex_t *to, double inc_weight) {
    if (from == to) {
        return;
    }


    auto elem_child = lower_bound(from->children.begin(),
                                   from->children.end(), to);

    auto elem_parent = lower_bound(to->parents.begin(),
                            to->parents.end(), from);

    (*elem_child).weight += inc_weight;
    from->rate += inc_weight;
    (*elem_parent).weight += inc_weight;
}

static void avl_free(avl_vec_vertex_t *vertex) {
    if (vertex == nullptr) {
        return;
    }

    avl_free(vertex->left);
    avl_free(vertex->right);
    free(vertex);
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

    for (auto child : graph->children) {
        _get_abs_vertex(abs_vertex, child.vertex);
    }
}

static vertex_t *get_abs_vertex(vertex_t *graph) {
    reset_graph_visited(graph);
    vertex_t *abs_vertex = nullptr;
    _get_abs_vertex(&abs_vertex, graph);

    return abs_vertex;
}

static int kingman_visit_vertex(vertex_t **out_initial_vertex,
                                vec_entry_t *initial_state,
                                vertex_t *abs_vertex,
                                const size_t m) {
    avl_vec_vertex_t *bst = nullptr;

    queue<vertex_t *> vertices_to_visit;
    vertex_t *initial_vertex = new vertex_t(initial_state, m);
    vertices_to_visit.push(initial_vertex);
    vec_entry_t *v = (vec_entry_t *) malloc(sizeof(vec_entry_t) * m);
    avl_vec_vertex_t *bst_vertex;

    queue<vertex_t*> sorting_queue;

    while (!vertices_to_visit.empty()) {
        vertex_t *vertex = vertices_to_visit.front();
        sorting_queue.push(vertex);
        vertices_to_visit.pop();
        memcpy(v, vertex->state, sizeof(vec_entry_t) * m);
        size_t n_remaining = 0;

        for (vec_entry_t i = 0; i < m; i++) {
            n_remaining += v[i];
        }

        for (vec_entry_t i = 0; i < m; i++) {
            if (v[i] == 0) {
                continue;
            }

            for (vec_entry_t j = i; j < m; j++) {
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
                            vertex_t *to = new vertex_t(new_state, m);

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

                    add_edge_unsorted(vertex, new_vertex, t);
                }
            }
        }
    }

    while (!sorting_queue.empty()) {
        auto vertex = sorting_queue.front();
        sorting_queue.pop();

        sort(vertex->children.begin(), vertex->children.end());
        sort(vertex->parents.begin(), vertex->parents.end());
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

    vec_entry_t *initial = (vec_entry_t*)calloc(m, sizeof(vec_entry_t));
    initial[0] = n;

    vec_entry_t *mrca = (vec_entry_t*)calloc(m, sizeof(vec_entry_t));
    mrca[m-1] = 1;

    vertex_t *absorbing_vertex = new vertex_t(mrca, m);

    vertex_t *state_graph;

    kingman_visit_vertex(&state_graph,
            initial, absorbing_vertex,
            m);

    vec_entry_t *start_state = (vec_entry_t*)calloc(m, sizeof(vec_entry_t));

    vertex_t *start = new vertex_t(start_state, m);

    add_edge(start, state_graph, 1);

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
    fprintf(stream, ":\n");


    for (size_t i = 0; i < vertex->children.size(); ++i) {
        fprintf(stream, "\t");
        fprintf(stream, "(%f) ", vertex->children[i].weight);
        print_vector_spacing(stream, vertex->children[i].vertex->state,
                             vec_length, vec_spacing);

        if (indexed) {
            fprintf(stream, " (%zu)", vertex->vertex_index);
        }

        fprintf(stream, "\n");
    }

    fprintf(stream, "\n");
    for (size_t i = 0; i < vertex->children.size(); i++) {
        _print_graph_list(stream, vertex->children[i].vertex,
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

void mph_cov_assign_vertex_all(vertex_t *vertex, size_t m) {
    if (vertex->visited) {
        return;
    }

    vertex->visited = true;

    if (vertex->parents.empty()) {
        // Starting vertex
        vertex->prob = 1.0f;
    } else {
        vertex->prob = 0.0f;
    }

    for (size_t i = 0; i < vertex->parents.size(); i++) {
        vertex_t *parent = vertex->parents[i].vertex;

        mph_cov_assign_vertex_all(parent, m);

        vertex->prob += vertex->parents[i].weight/parent->rate * parent->prob;
    }


    vertex->exp = vector<double>();

    for (size_t j = 0; j < m; ++j) {
        if (vertex->rate != 0) {
            vertex->exp.push_back(vertex->prob * vertex->state[j] / vertex->rate);
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

    for (auto child : vertex->children) {
        mph_cov_assign_desc_all(child.vertex, m);
    }

    for (size_t i = 0; i < vertex->children.size(); i++) {
        for (size_t j = 0; j < m; ++j) {
            vertex->desc[j] += vertex->children[i].weight / vertex->rate * vertex->children[i].vertex->desc[j];
        }
    }

    vector<double> exp(m);

    for (size_t j = 0; j < m; ++j) {
        if (vertex->rate != 0) {
            exp[j] = vertex->state[j] / vertex->rate;
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

    for (auto child : vertex->children) {
        _mph_cov_all(child.vertex, m);
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

void graph_free(vertex_t *graph) {
    queue<vertex_t*> to_free = enqueue_vertices(graph);

    while (!to_free.empty()) {
        vertex_t *vertex = to_free.front();
        to_free.pop();
        delete vertex;
    }
}

double calculate_rate(vertex_t *vertex) {
    double rate = 0;

    for (size_t i = 0; i < vertex->children.size(); ++i) {
        rate += vertex->children[i].weight;
    }

    return rate;
}

int reward_transform(vertex_t *graph, double (*reward_func)(vertex_t*)) {
    queue<vertex_t *> queue = enqueue_vertices(graph);

    while (!queue.empty()) {
        vertex_t *vertex = queue.front();
        queue.pop();

        if (vertex->children.empty()) {
            // Absorbing vertex
            continue;
        }

        if (vertex->parents.empty()) {
            continue;
        }

        double reward = reward_func(vertex);

        if (reward == 0) {
            // Take all my edges and add to my parent instead.
            auto parents = vertex->parents;
            auto children = vertex->children;

            for (auto child : children) {
                double prob = child.weight / vertex->rate;

                for (auto parent : parents) {
                    double weight = prob * parent.weight;

                    bool found = binary_search(parent.vertex->children.begin(),
                            parent.vertex->children.end(), child);

                    if (!found) {
                        add_edge(parent.vertex, child.vertex, weight);
                    } else {
                        increase_weight(parent.vertex, child.vertex, weight);
                    }
                }
            }

            for (auto parent : parents) {
                remove_edge(parent.vertex, vertex);
            }

            for (auto child : children) {
                remove_edge(vertex, child.vertex);
            }

            delete vertex;
        } else {
            for (size_t i = 0; i < vertex->children.size(); ++i) {
                vertex->children[i].weight /= reward;
                auto elem = lower_bound(vertex->children[i].vertex->parents.begin(), vertex->children[i].vertex->parents.end(),
                                 vertex);
                (*elem).weight /= reward;
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
    queue<vertex_t*> queue;
    size_t index = 0;

    // The absorbing vertex should have index 0
    queue.push(abs_vertex);
    queue.push(graph);

    while(!queue.empty()) {
        vertex_t *vertex = queue.front();
        queue.pop();

        if (vertex->visited) {
            continue;
        }

        vertex->visited = true;

        vertex->vertex_index = index++;

        for (auto child : vertex->children) {
            queue.push(child.vertex);
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

    for (size_t i = 0; i < vertex->children.size(); ++i) {
        vertex_t *child = vertex->children[i].vertex;

        weights[vertex->vertex_index][child->vertex_index] = vertex->children[i].weight;
        weights[vertex->vertex_index][vertex->vertex_index] -= vertex->children[i].weight;
    }

    vertices[vertex->vertex_index] = vertex;

    for (auto child : vertex->children) {
        insert_into_weight_mat(weights, vertices, child.vertex);
    }
}

int graph_as_mat(double ***weights, vertex_t ***vertices, size_t *out_size, vertex_t *graph) {
    size_t largest_index;
    label_vertex_index(&largest_index, graph);
    reset_graph_visited(graph);
    size_t size = largest_index+1;
    *out_size = size;
    *weights = (double**)malloc(sizeof(double*)*size);

    for (size_t i = 0; i < size; ++i) {
        (*weights)[i] = (double*)calloc(size, sizeof(double));
    }

    *vertices = (vertex_t **)calloc(size, sizeof(vertex_t*));

    insert_into_weight_mat(*weights, *vertices, graph);
    return 0;
}

