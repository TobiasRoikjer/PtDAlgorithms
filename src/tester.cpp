#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include "phase.h"

void assert(bool a) {
    if (!a) {
        exit(1);
    }
}

vertex_t *s,*a,*b,*c,*d,*e,*f,*t;

double reward_by_index(vertex_t *vertex) {
    if (vertex == s) {
        return 1;
    }

    if (vertex == a) {
        return 0;
    }

    if (vertex == b) {
        return 0;
    }

    if (vertex == c) {
        return 0;
    }

    if (vertex == d) {
        return 0;
    }

    if (vertex == e) {
        return 1;
    }

    if (vertex == f) {
        return 1;
    }

    if (vertex == t) {
        return 1;
    }

    exit(1);
}


int main(int argv, char **argc) {
    vec_entry_t *states =(vec_entry_t*) calloc(10, sizeof(vec_entry_t));
    vec_entry_t *statea =(vec_entry_t*) calloc(10, sizeof(vec_entry_t));
    vec_entry_t *stateb =(vec_entry_t*) calloc(10, sizeof(vec_entry_t));
    vec_entry_t *statec =(vec_entry_t*) calloc(10, sizeof(vec_entry_t));
    vec_entry_t *stated =(vec_entry_t*) calloc(10, sizeof(vec_entry_t));
    vec_entry_t *statee =(vec_entry_t*) calloc(10, sizeof(vec_entry_t));
    vec_entry_t *statef =(vec_entry_t*) calloc(10, sizeof(vec_entry_t));
    vec_entry_t *statet =(vec_entry_t*) calloc(10, sizeof(vec_entry_t));

    *states = 1;
    *statea = 2;
    *stateb = 3;
    *statec = 4;
    *stated = 5;
    *statee = 6;
    *statef = 7;
    *statet = 8;

    s = vertex_init(states, vector<double>(), 0);
    a = vertex_init(statea, vector<double>(), 0);
    b = vertex_init(stateb, vector<double>(), 0);
    c = vertex_init(statec, vector<double>(), 0);
    d = vertex_init(stated, vector<double>(), 0);
    e = vertex_init(statee, vector<double>(), 0);
    f = vertex_init(statef, vector<double>(), 0);
    t = vertex_init(statet, vector<double>(), 0);

    vertex_add_edge(s, a, 0.3);
    vertex_add_edge(s, b, 0.7);
    vertex_add_edge(a, b, 3);
    vertex_add_edge(b, c, 4);
    vertex_add_edge(b, d, 5);
    vertex_add_edge(c, t, 6);
    vertex_add_edge(d, c, 7);
    reward_transform(s, reward_by_index);

    double **mat;
    vertex_t **vertices;
    size_t size;

    graph_as_mat(&mat, &vertices, &size, s);

    for (size_t i = 2; i < size; ++i) {
        fprintf(stdout, "%f ", mat[1][i]);
    }

    fprintf(stdout, "\n\n");

    for (size_t i = 2; i < size; ++i) {
        for (size_t j = 2; j < size; ++j) {
            fprintf(stdout, "%f ", mat[i][j]);
        }

        fprintf(stdout, "\n");
    }

    return 0;
}

/*
int main(int argv, char **argc) {
    vertex_t *a = vertex_init();
    vertex_t *b = vertex_init();
    vertex_t *c = vertex_init();
    vertex_t *d = vertex_init();
    vertex_t *e = vertex_init();
    vertex_t *f = vertex_init();

    vertex_t *me = vertex_init();
    ll_insert(me, b, 1);
    ll_insert(me, c, 1);
    ll_insert(me, a, 1);
    ll_insert(me, f, 1);
    ll_insert(me, e, 1);
    ll_insert(me, d, 1);

    assert(me->edges->vertex == a);
    assert(me->edges->next->vertex == b);
    assert(me->edges->next->next->vertex == c);
    assert(me->edges->next->next->next->vertex == d);
    assert(me->edges->next->next->next->next->vertex == e);
    assert(me->edges->next->next->next->next->next->vertex == f);

    me = vertex_init();
    ll_insert_p(b, me);
    ll_insert_p(c, me);
    ll_insert_p(a, me);
    ll_insert_p(f, me);
    ll_insert_p(e, me);
    ll_insert_p(d, me);

    assert(me->parents->vertex == a);
    assert(me->parents->next->vertex == b);
    assert(me->parents->next->next->vertex == c);
    assert(me->parents->next->next->next->vertex == d);
    assert(me->parents->next->next->next->next->vertex == e);
    assert(me->parents->next->next->next->next->next->vertex == f);

    me = vertex_init();
    ll_insert_p_existing(b, me);
    ll_insert_p_existing(b, me);
    ll_insert_p_existing(c, me);
    ll_insert_p_existing(c, me);
    ll_insert_p_existing(a, me);
    ll_insert_p_existing(a, me);
    ll_insert_p_existing(f, me);
    ll_insert_p_existing(f, me);
    ll_insert_p_existing(e, me);
    ll_insert_p_existing(e, me);
    ll_insert_p_existing(d, me);
    ll_insert_p_existing(d, me);
    ll_insert_p_existing(b, me);
    ll_insert_p_existing(b, me);
    ll_insert_p_existing(c, me);
    ll_insert_p_existing(c, me);
    ll_insert_p_existing(a, me);
    ll_insert_p_existing(a, me);
    ll_insert_p_existing(f, me);
    ll_insert_p_existing(f, me);
    ll_insert_p_existing(e, me);
    ll_insert_p_existing(e, me);
    ll_insert_p_existing(d, me);
    ll_insert_p_existing(d, me);

    assert(me->parents->vertex == a);
    assert(me->parents->next->vertex == b);
    assert(me->parents->next->next->vertex == c);
    assert(me->parents->next->next->next->vertex == d);
    assert(me->parents->next->next->next->next->vertex == e);
    assert(me->parents->next->next->next->next->next->vertex == f);

    a = vertex_init();
    b = vertex_init();
    c = vertex_init();
    d = vertex_init();
    e = vertex_init();
    f = vertex_init();
    me = vertex_init();
    vertex_t *other = vertex_init();
    ll_insert(me, c, 1);
    ll_insert(me, d, 1);

    ll_insert(other, a, 1);
    ll_insert(other, d, 1);
    ll_insert(other, f, 1);
    ll_insert(other, b, 1);
    ll_insert(other, e, 1);

    ll_insert_or_inc_list(me, other->edges);

    assert(me->edges->vertex == a);
    assert(me->edges->next->vertex == b);
    assert(me->edges->next->next->vertex == c);
    assert(me->edges->next->next->next->vertex == d);
    assert(me->edges->next->next->next->next->vertex == e);
    assert(me->edges->next->next->next->next->next->vertex == f);

    assert(me->edges->weight == 1);
    assert(me->edges->next->weight == 1);
    assert(me->edges->next->next->weight == 1);
    assert(me->edges->next->next->next->weight == 2);
    assert(me->edges->next->next->next->next->weight == 1);
    assert(me->edges->next->next->next->next->next->weight == 1);

    ll_insert_or_inc_list(me, other->edges);

    assert(me->edges->vertex == a);
    assert(me->edges->next->vertex == b);
    assert(me->edges->next->next->vertex == c);
    assert(me->edges->next->next->next->vertex == d);
    assert(me->edges->next->next->next->next->vertex == e);
    assert(me->edges->next->next->next->next->next->vertex == f);

    assert(me->edges->weight == 2);
    assert(me->edges->next->weight == 2);
    assert(me->edges->next->next->weight == 1);
    assert(me->edges->next->next->next->weight == 3);
    assert(me->edges->next->next->next->next->weight == 2);
    assert(me->edges->next->next->next->next->next->weight == 2);

    assert(a->parents->vertex == me);
    assert(a->parents->next->vertex == other);

    assert(b->parents->vertex == me);
    assert(b->parents->next->vertex == other);

    assert(c->parents->vertex == me);

    assert(d->parents->vertex == me);
    assert(d->parents->next->vertex == other);

    assert(e->parents->vertex == me);
    assert(e->parents->next->vertex == other);

    assert(f->parents->vertex == me);
    assert(f->parents->next->vertex == other);

    return 0;
}*/