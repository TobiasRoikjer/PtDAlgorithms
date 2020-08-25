#include <cstdlib>
#include <algorithm>
#include <cstdio>
#include "phase.h"

inline bool has_child(vertex_t *vertex, vertex_t *child) {
    for (auto c : vertex->children) {
        if (c.vertex == child) {
            return true;
        }
    }

    return false;
}

vertex_t *generate_graph(unsigned int seed,
        size_t n_states, size_t n_edges,
        size_t n_zero_rewards) {
    srand(seed);
    vertex_t *ipv = new vertex_t(nullptr, {1.0f}, 0);
    vertex_t *abs = new vertex_t(nullptr, {0.0f}, 0);
    vertex_t **vertices = (vertex_t**)calloc(n_states, sizeof(vertex_t*));

    ipv->vertex_index = 1;
    abs->vertex_index = 0;

    for (size_t i = 0; i < n_states; ++i) {
        vertices[i] = new vertex_t(nullptr, {2.0f}, 0);
        vertices[i]->vertex_index = i+2;
        add_edge(ipv, vertices[i], 1.0f);
        add_edge(vertices[i], abs, 1.0f);
    }

    vector<pair<size_t, size_t>> combinations;

    for (size_t i = 0; i < n_states; ++i) {
        for (size_t j = 0; j < n_states; ++j) {
            if (i == j) {
                continue;
            }

            combinations.push_back(pair<size_t,size_t>(i,j));
        }
    }

    random_shuffle(combinations.begin(), combinations.end());

    for (size_t i = 0; i < n_edges; ++i) {
        add_edge_unsorted(vertices[combinations[i].first],
                vertices[combinations[i].second], 1.0f);
    }

    for (size_t i = 0; i < n_states; ++i) {
        sort(vertices[i]->children.begin(), vertices[i]->children.end());
        sort(vertices[i]->parents.begin(), vertices[i]->parents.end());
    }

    for (size_t i = 0; i < n_zero_rewards; ++i) {
        vertices[i]->rewards[0] = 0.0f;
    }

    return ipv;
}
