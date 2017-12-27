#include "shuffler.h"

namespace shuffler {
    Shuffler::Shuffler(unsigned n, unsigned hard_stop) :graph(n), n(n), 
                                                 hard_stop(hard_stop) {
        random_device rd;
        mersenne = mt19937(rd());
        edges.resize(n, vector<bool>(n, false));
        rng_for_ends = uniform_int_distribution<unsigned>(0, 3);
    }

    void Shuffler::add_edge(unsigned v, unsigned u) {
        unsigned edge_num = tokens.size();
        permutation.push_back(edge_num);
        edges[v][u] = true;
        edges[u][v] = true;
        ends.emplace_back(v, u);
        tokens.push_back(std::move(graph.add(v, u)));
    }

    vector<unsigned> Shuffler::perm() {
        return permutation;
    }

    vector<pair<unsigned, unsigned>> Shuffler::get_edges(){
        return ends;
    }

    bool Shuffler::is_connected() {
        return graph.is_connected();
    }

    bool Shuffler::shake_it(){
        while (true) {
            --hard_stop;
            if (hard_stop == 0) {
                return false;
            }

            unsigned e1 = rng_edge(mersenne);
            unsigned e2 = rng_edge(mersenne);
            unsigned ends_direct = rng_for_ends(mersenne);
            auto edge = ends[e1];
            auto other = ends[e2];

            unsigned v = edge.first;
            unsigned u = edge.second;
            unsigned w = other.first;
            unsigned z = other.second;

            if (ends_direct / 2 == 0) {
                std::swap(v, u);
            }

            if (ends_direct % 2 == 0) {
                std::swap(w, z);
            }

            if (e1 == e2)
                continue;
            if (v == z || w == u || v == w || u == z) {
                continue;
            }
            if (edges[v][z] || edges[w][u]) {
                continue;
            }
            if (do_flip(e1, e2, v, u, w, z)) {
                return true;
            }
        }
    }

    bool Shuffler::do_flip(unsigned e1, unsigned e2, unsigned v, 
                           unsigned u, unsigned w, unsigned z) {
        graph.remove(std::move(tokens[e1]));
        graph.remove(std::move(tokens[e2]));
        tokens[e1] = std::move(graph.add(v, z));
        tokens[e2] = std::move(graph.add(w, u));
        if (graph.is_connected()) {
            edges[v][u] = false;
            edges[u][v] = false;
            edges[w][z] = false;
            edges[z][w] = false;
            edges[v][z] = true;
            edges[z][v] = true;
            edges[w][u] = true;
            edges[u][w] = true;
            ends[e1] = make_pair(v, z);
            ends[e2] = make_pair(w, u);
            swap(permutation[e1], permutation[e2]);
            return true;
        }
        graph.remove(std::move(tokens[e1]));
        graph.remove(std::move(tokens[e2]));
        tokens[e1] = std::move(graph.add(v, u));
        tokens[e2] = std::move(graph.add(w, z));
        return false;
    }
}
