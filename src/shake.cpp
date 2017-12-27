#include <Rcpp.h>
#include "shuffler.h"

using namespace Rcpp;
using shuffler::Shuffler;

// [[Rcpp::export]]
Rcpp::List shake_internal(Rcpp::List& graph,
                          Rcpp::IntegerVector size,
                          Rcpp::IntegerVector permutations,
                          Rcpp::IntegerVector hard_stop) {
    auto from = as<IntegerVector>(graph["from"]);
    auto to = as<IntegerVector>(graph["to"]);

    Shuffler shuffler(size[0], hard_stop[0]);
    for (int i = 0; i < from.size(); i++) {
        shuffler.add_edge(from(i) - 1, to(i) - 1);
    }

    if (!shuffler.is_connected()) {
        stop("Given network is disconnected");
    }

    for (unsigned i = 0; i < permutations[0]; i++) {
        if (!shuffler.shake_it()) {
            warning("Reached limit of permutation attempts(may be this graph is unique?)");
            break;
        }
    }
    auto permutation = shuffler.perm();

    auto ends = shuffler.get_edges();
    vector<unsigned> f;
    vector<unsigned> t;
    for (int i = 0; i < permutation.size(); i++){
        f.push_back(ends[permutation[i]].first + 1);
        t.push_back(ends[permutation[i]].second + 1);
    }

    IntegerVector from_new(f.begin(), f.end());
    IntegerVector to_new(t.begin(), t.end());
    Rcpp::DataFrame df = DataFrame::create(Named("from") = from_new,
                                           Named("to") = to_new);
    return df;
}
