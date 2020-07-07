#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <numeric>
#include <string.h>
#include "../api/cpp/ptdalgorithmscpp.h"

using namespace std;
using namespace ptdalgorithms;

#define assert(a) {\
    do {\
    if (!(a)) {\
        fprintf(stderr, "ASSERTION FAILED %s %d\n", __FILE__, __LINE__);\
        exit(1);\
    }\
    } while (0);\
}

void test_basic_ptd_graph() {
    Graph g(4);

    vector<int> state(4);
    state[0] = 0;
    Vertex v1 = g.create_vertex(state);
    state[0] = 1;
    Vertex v2 = g.create_vertex(state);
    state[0] = 2;
    Vertex v3 = g.create_vertex(state);
    state[0] = 3;
    Vertex v4 = g.create_vertex(state);
    state[0] = 4;
    Vertex v5 = g.create_vertex(state);

    assert(g.vertices().size() == 5 + 1);
    assert(g.state_length() == 4);
    assert(g.vertices()[0] == g.starting_vertex());
    assert(g.vertices()[1] == v1);
    assert(g.vertices()[2] == v2);
    assert(g.vertices()[3] == v3);
    assert(g.vertices()[4] == v4);
    assert(g.vertices()[5] == v5);
}

void test_basic_ptd_graph_edges() {
    Graph g(4);

    vector<int> state(4);
    state[0] = 0;
    Vertex v1 = g.create_vertex(state);
    state[0] = 1;
    Vertex v2 = g.create_vertex(state);
    state[0] = 2;
    Vertex v3 = g.create_vertex(state);
    state[0] = 3;
    Vertex v4 = g.create_vertex(state);
    state[0] = 4;
    Vertex v5 = g.create_vertex(state);

    g.starting_vertex().add_edge(v1, 0.4);
    g.starting_vertex().add_edge(v2, 0.6);
    v3.add_edge(v2, 1.0);
    v3.add_edge(v1, 2.0);
    vector<double> edge_state;
    edge_state.push_back(1);
    edge_state.push_back(2);
    v3.add_edge_parameterized(v4, 2.0, edge_state);

    Graph g2 = g;
    Graph g3(g);
    Vertex v10 = v3;
    Vertex v11(v3);

    assert(g.starting_vertex().edges().size() == 2);
    assert(g.starting_vertex().edges()[0].weight() == 0.4);
    assert(g.starting_vertex().edges()[1].weight() == 0.6);
    assert(v3.edges().size() == 3);
    assert(v3.edges()[0].weight() == 1.0);
    assert(v3.edges()[0].to() == v2);
    assert(v3.edges()[1].weight() == 2.0);
    assert(v3.edges()[1].to() == v1);
    assert(v3.edges()[2].weight() == 2.0);
    assert(v3.edges()[2].to() == v4);
    ParameterizedEdge e = v3.parameterized_edges()[2];
    assert(e.edge_state(2)[0] == 1);
    assert(e.edge_state(2)[1] == 2);

    vector<double> scalars;
    scalars.push_back(9);
    scalars.push_back(5);
    g.update_weights_parameterized(scalars);
    assert(v3.edges()[2].weight() == 19);
    Graph gg = g.clone();
}


void test_expected_entry_visits() {
    struct ptd_graph *graph = ptd_graph_create(4);

    struct ptd_vertex *S = graph->starting_vertex;
    S->state[0] = 0;
    fprintf(stderr, "S: %p\n", (void *) S);
    struct ptd_vertex *T = ptd_vertex_create(graph);
    T->state[0] = 1;
    fprintf(stderr, "T: %p\n", (void *) T);

    struct ptd_vertex *A = ptd_vertex_create(graph);
    A->state[0] = 2;
    fprintf(stderr, "A: %p\n", (void *) A);
    struct ptd_vertex *B = ptd_vertex_create(graph);
    B->state[0] = 3;
    fprintf(stderr, "B: %p\n", (void *) B);
    struct ptd_vertex *C = ptd_vertex_create(graph);
    C->state[0] = 4;
    fprintf(stderr, "C: %p\n", (void *) C);
    struct ptd_vertex *D = ptd_vertex_create(graph);
    D->state[0] = 5;
    fprintf(stderr, "D: %p\n", (void *) D);
    struct ptd_vertex *E = ptd_vertex_create(graph);
    E->state[0] = 6;
    fprintf(stderr, "E: %p\n", (void *) E);
    struct ptd_vertex *F = ptd_vertex_create(graph);
    F->state[0] = 7;
    fprintf(stderr, "F: %p\n", (void *) F);
    struct ptd_vertex *G = ptd_vertex_create(graph);
    G->state[0] = 8;
    fprintf(stderr, "G: %p\n", (void *) G);
    struct ptd_vertex *H = ptd_vertex_create(graph);
    H->state[0] = 9;
    fprintf(stderr, "H: %p\n", (void *) H);
    struct ptd_vertex *I = ptd_vertex_create(graph);
    I->state[0] = 10;
    fprintf(stderr, "I: %p\n", (void *) I);
    struct ptd_vertex *J = ptd_vertex_create(graph);
    J->state[0] = 11;
    fprintf(stderr, "J %p\n", (void *) J);
    struct ptd_vertex *K = ptd_vertex_create(graph);
    K->state[0] = 12;
    fprintf(stderr, "K %p\n", (void *) K);
    struct ptd_vertex *L = ptd_vertex_create(graph);
    L->state[0] = 13;
    fprintf(stderr, "L %p\n", (void *) L);


    ptd_graph_add_edge(S, A, 1);

    ptd_graph_add_edge(A, B, 3);
    ptd_graph_add_edge(B, C, 4);
    ptd_graph_add_edge(C, A, 4);
    ptd_graph_add_edge(B, D, 2);
    ptd_graph_add_edge(D, E, 5);
    ptd_graph_add_edge(E, F, 5);
    ptd_graph_add_edge(E, I, 15);
    ptd_graph_add_edge(F, G, 1);
    ptd_graph_add_edge(G, H, 1);
    ptd_graph_add_edge(H, F, 1);
    ptd_graph_add_edge(H, G, 1);

    ptd_graph_add_edge(H, T, 1);
    ptd_graph_add_edge(I, L, 1);
    ptd_graph_add_edge(L, T, 1);

    ptd_graph_add_edge(C, J, 1);

    ptd_graph_add_edge(J, K, 1);
    ptd_graph_add_edge(K, J, 1);

    ptd_graph_add_edge(K, E, 1);

    ptd_graph_add_edge(K, G, 1);

    Graph g(graph);
    vector<double> exp_waiting_time = g.expected_waiting_time();
    /*
    vector<double> exp_visits = g.expected_visits();

    assert(abs(exp_visits[S->index] - 1) < 0.01);
    assert(abs(exp_visits[T->index] - 1) < 0.01);
    assert(abs(exp_visits[A->index] - 2.14) < 0.01);
    assert(abs(exp_visits[B->index] - 2.14) < 0.01);
    assert(abs(exp_visits[C->index] - 1.42) < 0.01);
    assert(abs(exp_visits[D->index] - 0.71) < 0.01);
    assert(abs(exp_visits[E->index] - 0.85) < 0.01);
    assert(abs(exp_visits[J->index] - 0.42) < 0.01);
    assert(abs(exp_visits[K->index] - 0.42) < 0.01);
    assert(abs(exp_visits[G->index] - 1.07) < 0.01);
    assert(abs(exp_visits[H->index] - 1.07) < 0.01);*/

    g.normalize();

    std::vector<int> rewardsi(g.vertices().size());
    std::vector<double> rewardsd(g.vertices().size());

    for (size_t j = 0; j < g.vertices().size(); ++j) {
        rewardsi[j] = (int) (j) % 3;
        rewardsd[j] = j % 3;
    }

    long double sum = 0;

    for (size_t i = 0; i < 400000; ++i) {
        long double r = g.dph_random_sample(rewardsd) / 400000.0;
        sum += r * r;
    }


    long double sum2 = 0;

    fprintf(stderr, "\n=================================\n");

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        struct ptd_vertex *vertex = graph->vertices[i];

        fprintf(stderr, "Vertex %i %p:\n", vertex->state[0], (void *) vertex);

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            fprintf(stderr, "\tChild %i %p (%f)\n", vertex->edges[j]->to->state[0], (void *) vertex->edges[j]->to,
                    vertex->edges[j]->weight);
        }
    }

    g.dph_reward_transform(rewardsi);

    fprintf(stderr, "\n=================================\n");

    for (size_t i = 0; i < graph->vertices_length; ++i) {
        struct ptd_vertex *vertex = graph->vertices[i];

        fprintf(stderr, "Vertex %i %p:\n", vertex->state[0], (void *) vertex);

        for (size_t j = 0; j < vertex->edges_length; ++j) {
            fprintf(stderr, "\tChild %i %p (%f)\n", vertex->edges[j]->to->state[0], (void *) vertex->edges[j]->to,
                    vertex->edges[j]->weight);
        }
    }

    srand(1234);

    for (size_t i = 0; i < 400000; ++i) {
        long double r = g.dph_random_sample() / 400000.0;
        sum2 += r * r;
    }

    assert(fabs(sum2 - sum) < 0.0001);
}

void test_rabbit() {
    size_t state_vector_length = 2;
    Graph graph(state_vector_length);
    Vertex starting_vertex = graph.vertices()[0];
    vector<int> initial_state;
    initial_state.push_back(15);
    initial_state.push_back(0);
    Vertex initial_vertex = graph.create_vertex(initial_state);
    starting_vertex.add_edge(initial_vertex, 1);
    int index = 1;

    while (index < graph.vertices().size()) {
        Vertex vertex = graph.vertices()[index];
        vector<int> state = vertex.state();

        if (state[0] > 0) {
            // Rabbit jump left to right
            {
                vector<int> child_state;
                child_state.push_back(state[0] - 1);
                child_state.push_back(state[1] + 1);
                Vertex child = graph.find_or_create_vertex(child_state);
                vertex.add_edge(
                        child,
                        1
                );
            }

            // Island flooding
            {
                vector<int> child_state;
                child_state.push_back(0);
                child_state.push_back(state[1]);
                Vertex child = graph.find_or_create_vertex(child_state);
                vector<double> s;
                s.push_back(1);
                vertex.add_edge_parameterized(
                        child,
                        0,
                        s
                );
            }
        }

        if (state[1] > 0) {
            // Rabbit jump right to left
            {
                vector<int> child_state;
                child_state.push_back(state[0] + 1);
                child_state.push_back(state[1] - 1);
                Vertex child = graph.find_or_create_vertex(child_state);
                vertex.add_edge(
                        child,
                        1
                );
            }

            // Island flooding
            {
                vector<int> child_state;
                child_state.push_back(state[0]);
                child_state.push_back(0);
                Vertex child = graph.find_or_create_vertex(child_state);
                vector<double> s;
                s.push_back(2);
                vertex.add_edge_parameterized(
                        child,
                        0,
                        s
                );
            }
        }


        index++;
    }

    fprintf(stderr, "Done building\n");

    fprintf(stderr, "%f\n", graph.expected_waiting_time()[0]);
    vector<double> scalars;
    scalars.push_back(10);
    graph.update_weights_parameterized(scalars);
    fprintf(stderr, "%f\n", graph.expected_waiting_time()[0]);
    graph.random_sample();
    vector<double> rewards(graph.vertices().size() * 3);
    graph.mph_random_sample(rewards, 3);

    {
        Graph graphD(4);

        Vertex SD = graphD.starting_vertex();
        Vertex TD = graphD.create_vertex();
        Vertex AD = graphD.create_vertex();
        Vertex BD = graphD.create_vertex();

        SD.add_edge(AD, 1);
        AD.add_edge(BD, 0.8);
        BD.add_edge(TD, 0.5);

        long double re1 = 0, re1v = 0;

        srand(1234);

        for (size_t j = 0; j < 10000; ++j) {
            long double r = graphD.dph_random_sample();
            assert(!isnan(r));
            re1 += r;
            re1v += r * r;
        }

        vector<double> eres = graphD.expected_waiting_time();
        vector<double> eres2 = graphD.expected_waiting_time(eres);

        double eexp = eres[0];
        double evar = eres2[0];


        vector<double> new_rewards = graphD.dph_normalize();
        long double re2 = 0, re2v = 0;

        srand(1234);

        for (size_t j = 0; j < 10000; ++j) {
            long double r = graphD.dph_random_sample(new_rewards);
            re2 += r;
            re2v += r * r;
        }

        assert((re1 / 10000.0 - re2 / 10000.0) < 0.01);
        assert((re1 / 10000.0 - eexp) < 0.01);
        assert((re1v / 10000.0 - re2v / 10000.0) < 0.01);
        assert((re1v / 10000.0 - (evar * 2 - eexp)) < 0.01);
    }

    {
        Graph graphD(4);

        Vertex SD = graphD.starting_vertex();
        Vertex TD2 = graphD.create_vertex();
        Vertex BD = graphD.create_vertex();
        Vertex AD = graphD.create_vertex();
        Vertex TD = graphD.create_vertex();

        SD.add_edge(AD, 0.75);
        SD.add_edge(TD2, 0.25);
        AD.add_edge(BD, 0.8);
        BD.add_edge(TD, 0.5);
        BD.add_edge(TD2, 0.2);

        assert(fabs(graphD.dph_pmf(0) - 0.25) <= 0.0001);
        assert(fabs(graphD.dph_cdf(0) - 0.25) <= 0.0001);


        assert(fabs(graphD.dph_pmf(1)) <= 0.0001);
        assert(fabs(graphD.dph_cdf(1) - 0.25) <= 0.0001);

        assert(fabs(graphD.dph_pmf(2) - 0.75 * 0.8 * 0.7) <= 0.0001);
        assert(fabs(graphD.dph_cdf(2) - 0.25 - 0.75 * 0.8 * 0.7) <= 0.0001);

    }

    {
        Graph graphD(4);

        Vertex SD = graphD.starting_vertex();
        Vertex TD2 = graphD.create_vertex();
        Vertex BD = graphD.create_vertex();
        Vertex AD = graphD.create_vertex();
        Vertex TD = graphD.create_vertex();

        SD.add_edge(AD, 0.75);
        SD.add_edge(TD2, 0.25);
        AD.add_edge(BD, 0.8);
        BD.add_edge(TD, 0.5);
        BD.add_edge(TD2, 0.2);

        assert(fabs(graphD.dph_pmf(2) - 0.75 * 0.8 * 0.7) <= 0.0001);
        assert(fabs(graphD.dph_cdf(2) - 0.25 - 0.75 * 0.8 * 0.7) <= 0.0001);

    }

    {
        Graph graphD(4);

        Vertex SD = graphD.starting_vertex();
        Vertex TD2 = graphD.create_vertex();
        Vertex BD = graphD.create_vertex();
        Vertex AD = graphD.create_vertex();
        Vertex TD = graphD.create_vertex();

        SD.add_edge(AD, 0.75);
        SD.add_edge(TD2, 0.25);
        AD.add_edge(BD, 0.8);
        BD.add_edge(TD, 0.5);
        BD.add_edge(TD2, 0.2);

        graphD.dph_pmf(100);

        assert(fabs(graphD.dph_pmf(2) - 0.75 * 0.8 * 0.7) <= 0.0001);
        assert(fabs(graphD.dph_cdf(2) - 0.25 - 0.75 * 0.8 * 0.7) <= 0.0001);
        assert(fabs(graphD.dph_pmf(0) - 0.25) <= 0.0001);
        assert(fabs(graphD.dph_cdf(0) - 0.25) <= 0.0001);

    }
}

void test_2pmigTIME() {
#define popsize 10

    struct mig {
        int pop1[(popsize) * (popsize)];
        int pop2[(popsize) * (popsize)];
    };

    clock_t start_time = time(NULL);

    size_t state_length = sizeof(mig) / sizeof(int);

    Graph graphA(state_length);
    ptd_graph *g = graphA.c_graph();
    ptd_avl_tree *a = graphA.c_avl_tree();

    mig startA;
    memset(&startA, 0, sizeof(mig));
    startA.pop1[1] = popsize - 1;
    startA.pop1[0 + popsize] = popsize - 1;
    std::vector<int> in((int*)startA.pop1, (int*)startA.pop1 + state_length);
    graphA.find_or_create_vertex(in);
    mig child_state;

    for (size_t index = 1; index < g->vertices_length; ++index) {
        ptd_vertex *vertex = g->vertices[index];
        struct mig state = *((mig*)vertex->state);

        int count = 0;

        for (size_t i = 0; i < popsize; ++i) {
            for (size_t j = 0; j < popsize; ++j) {
                count += state.pop1[i + j * popsize];
                count += state.pop2[i + j * popsize];
            }
        }

        if (count == 0 || count == 1) {
            // Only one lineage left, absorb
            continue;
        }

        assert(count <= popsize * 2);

        for (size_t i = 0; i < popsize; ++i) {
            for (size_t j = 0; j < popsize; ++j) {
                if (state.pop1[i + j * popsize] == 0) {
                    continue;
                }
                for (size_t i2 = 0; i2 < popsize; ++i2) {
                    for (size_t j2 = 0; j2 < popsize; ++j2) {
                        if (i + j * popsize < i2 + j2 * popsize) {
                            continue;
                        }
                        if (state.pop1[i2 + j2 * popsize] == 0) {
                            continue;
                        }
                        if (i == i2 && j == j2) {
                            if (state.pop1[i2 + j2 * popsize] == 1) {
                                continue;
                            }
                        }

                        double weight;
                        if (i == i2 && j == j2) {
                            weight = state.pop1[i + j * popsize] * (state.pop1[i2 + j2 * popsize] - 1) / 2;
                        } else {
                            weight = state.pop1[i + j * popsize] * state.pop1[i2 + j2 * popsize];
                        }

                        size_t new_indexi = i + i2;
                        size_t new_indexj = j + j2;
                        memcpy(&child_state, &state, sizeof(mig));
                        child_state.pop1[i + j * popsize]--;
                        child_state.pop1[i2 + j2 * popsize]--;
                        child_state.pop1[new_indexi + new_indexj * popsize]++;
                        ptd_vertex *to = ptd_find_or_create_vertex(g, a, (int*)&child_state);
                        ptd_graph_add_edge(vertex, to, weight);
                    }
                }
            }
        }
    }

    size_t nedges = 0;

    clock_t end_time = time(NULL);


    for (size_t k = 0; k < graphA.vertices_length(); ++k) {
        nedges += g->vertices[k]->edges_length;
    }
    clock_t taken = end_time - start_time;
    fprintf(stderr, "Done creating graph. It has %zu vertices and %zu edges. Took %i seconds\n", graphA.vertices().size(),
            nedges, (int) (taken));
}

void test_pmf() {
    {
        Graph graph(2);
        Vertex a = graph.create_vertex();
        Vertex b = graph.create_vertex();

        graph.starting_vertex().add_edge(a, 1);
        a.add_edge(b, 1);

        assert(abs(graph.dph_pmf(0) - 0) < 0.0001);
        assert(abs(graph.dph_pmf(1) - 1) < 0.0001);
        assert(abs(graph.dph_pmf(2) - 0) < 0.0001);


        for (float time = 0; time < 1; time += 0.1) {
            assert(abs(graph.cdf(time, 1000) - (1 - exp(-time))) < 0.01);
            fprintf(stderr, "%f ", graph.cdf(time, 1000));
        }
        fprintf(stderr, "\n");
        for (float time = 0; time < 1; time += 0.1) {
            fprintf(stderr, "%f ", graph.pdf(time, 1000));
            assert(abs(graph.pdf(time, 1000) - (exp(-time))) < 0.01);
        }
        fprintf(stderr, "\n");
    }

    {
        Graph graph(2);
        Vertex a = graph.create_vertex();
        Vertex b = graph.create_vertex();

        graph.starting_vertex().add_edge(a, 1);
        a.add_edge(b, 3);

        for (float time = 0; time < 1; time += 0.1) {
            fprintf(stderr, "%f ", graph.cdf(time, 1000));
            assert(abs(graph.cdf(time, 1000) - (1 - exp(-3 * time))) < 0.01);
        }
        fprintf(stderr, "\n");
        for (float time = 0; time < 1; time += 0.1) {
            fprintf(stderr, "%f ", graph.pdf(time, 1000));
            assert(abs(graph.pdf(time, 1000) - (3*exp(-3 *time))) < 0.01);
        }
        fprintf(stderr, "\n");
    }
    {
        Graph graph(2);
        Vertex a = graph.create_vertex();
        Vertex b = graph.create_vertex();

        graph.starting_vertex().add_edge(a, 0.5);
        graph.starting_vertex().add_edge(b, 0.5);
        a.add_edge(b, 3);

        for (float time = 0; time < 1; time += 0.1) {
            fprintf(stderr, "%f ", graph.cdf(time, 1000));
            assert(abs(graph.cdf(time, 1000) - (0.5 * (1 - exp(-3 * time)) + 0.5)) < 0.01);
        }
        fprintf(stderr, "\n");
        assert(abs(graph.defect() - 0.5) < 0.01);
        for (float time = 0; time < 1; time += 0.1) {
            fprintf(stderr, "%f ", graph.pdf(time, 1000));
            assert(abs(graph.pdf(time, 1000) - (0.5*3*exp(-3 *time))) < 0.01);
        }
        fprintf(stderr, "\n");



        vector<double> aa = graph.stop_probability(1);
        vector<double> bb = graph.stop_probability(1);
    }
}

int main(int argc, char **argv) {
    //test_basic_ptd_graph();
    //test_basic_ptd_graph_edges();
    //test_2pmigTIME();
    //test_expected_entry_visits();
    //test_rabbit();
    test_pmf();
}
