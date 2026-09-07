#include <bits/stdc++.h>
#include "rseas_metric_k45.hpp"
using namespace std;

using Vec = rseas_metric::Vec;

namespace {

constexpr double EPS = 1e-15;

bool read_vec(istream& input, Vec& values) {
    for (double& value : values) {
        if (!(input >> value)) return false;
    }
    return true;
}

bool no_extra_values(istream& input) {
    input >> ws;
    return input.eof();
}

bool parse_k(const char* text, int& k) {
    try {
        size_t used = 0;
        const string input(text);
        k = stoi(input, &used);
        return used == input.size() && (k == 4 || k == 5);
    } catch (const exception&) {
        return false;
    }
}

double find_best_multiplier(const Vec& current,
                            const Vec& target,
                            const Vec& displacement) {
    double numerator = 0.0;
    double denominator = 0.0;

    // Minimize sum_i (current_i + u*displacement_i - target_i)^2.
    for (size_t i = 0; i < current.size(); ++i) {
        numerator += displacement[i] * (target[i] - current[i]);
        denominator += displacement[i] * displacement[i];
    }

    double u = 0.0;
    if (abs(denominator) > EPS) {
        u = numerator / denominator;
    }
    return clamp(u, 0.0, 2.0);
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <k:4|5> < scale_input\n";
        return 1;
    }

    int k = 0;
    if (!parse_k(argv[1], k)) {
        cerr << "Error: k must be 4 or 5.\n";
        return 1;
    }

    const size_t dimension = rseas_metric::dimension_for_k(k);
    Vec current(dimension, 0.0);
    Vec target(dimension, 0.0);
    Vec displacement(dimension, 0.0);

    if (!read_vec(cin, current) ||
        !read_vec(cin, target) ||
        !read_vec(cin, displacement) ||
        !no_extra_values(cin)) {
        cerr << "Error: bxk4one RMSE expected exactly "
             << (3 * dimension) << " numeric input values for k="
             << k << ".\n";
        return 1;
    }

    for (size_t i = 0; i < dimension; ++i) {
        current[i] = rseas_metric::safe_log(current[i]);
        target[i] = rseas_metric::safe_log(target[i]);

        // rpll writes log(current)-log(stage1). Negating it gives the observed
        // movement from the current graph toward the stage-one graph.
        displacement[i] = -displacement[i];
    }

    double u = find_best_multiplier(current, target, displacement);
    if (abs(u) <= EPS) u = 0.0;

    cout << fixed << setprecision(12) << u << '\n';
    return 0;
}

