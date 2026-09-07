#include <bits/stdc++.h>
#include "rseas_metric_k45.hpp"
using namespace std;

using Vec = rseas_metric::Vec;

namespace {

constexpr double EPS = 1e-12;
constexpr double U_MIN = 0.0;
constexpr double U_MAX = 2.0;

struct Result {
    double u = 0.0;
    double score = numeric_limits<double>::infinity();
};

void add_candidate(vector<double>& candidates, double u) {
    if (!isfinite(u) || u < U_MIN - EPS || u > U_MAX + EPS) {
        return;
    }
    candidates.push_back(clamp(u, U_MIN, U_MAX));
}

double score_at(double u,
                const Vec& current,
                const Vec& target,
                const Vec& displacement,
                const rseas_metric::Config& config) {
    Vec signed_errors(current.size(), 0.0);
    for (size_t i = 0; i < current.size(); ++i) {
        const double predicted = current[i] + u * displacement[i];
        signed_errors[i] = predicted - target[i];
    }
    return rseas_metric::score(signed_errors, config);
}

Result find_best_u(const Vec& current,
                   const Vec& target,
                   const Vec& displacement,
                   const rseas_metric::Config& config) {
    const size_t dimension = current.size();
    Vec intercept(dimension, 0.0);
    Vec slope(dimension, 0.0);
    for (size_t i = 0; i < dimension; ++i) {
        intercept[i] = current[i] - target[i];
        slope[i] = displacement[i];
    }

    vector<double> region_breaks{U_MIN, U_MAX};
    vector<double> candidates{U_MIN, U_MAX};

    // The active set of squared-hinge terms changes only where a residual is
    // +e or -e. Between consecutive points the objective is one quadratic.
    for (size_t i = 0; i < dimension; ++i) {
        if (abs(slope[i]) <= EPS) continue;

        const double plus_boundary =
            (config.error_range - intercept[i]) / slope[i];
        const double minus_boundary =
            (-config.error_range - intercept[i]) / slope[i];
        add_candidate(region_breaks, plus_boundary);
        add_candidate(region_breaks, minus_boundary);
    }

    sort(region_breaks.begin(), region_breaks.end());
    region_breaks.erase(unique(region_breaks.begin(), region_breaks.end(),
                               [](double a, double b) {
                                   return abs(a - b) <= EPS;
                               }),
                        region_breaks.end());
    candidates.insert(candidates.end(), region_breaks.begin(),
                      region_breaks.end());

    // Inside the error range, RSEAS is max_i |residual_i| - e. This is the
    // upper envelope of the 2*D signed residual lines.
    struct Line {
        double intercept;
        double slope;
    };
    vector<Line> lines;
    lines.reserve(2 * dimension);
    for (size_t i = 0; i < dimension; ++i) {
        lines.push_back({intercept[i], slope[i]});
        lines.push_back({-intercept[i], -slope[i]});
    }

    for (size_t i = 0; i < lines.size(); ++i) {
        for (size_t j = i + 1; j < lines.size(); ++j) {
            const double denominator = lines[i].slope - lines[j].slope;
            if (abs(denominator) <= EPS) continue;
            const double u =
                (lines[j].intercept - lines[i].intercept) / denominator;
            add_candidate(candidates, u);
        }
    }

    // In every outside interval, find the exact stationary point of the active
    // weighted quadratic. Midpoints also cover constant pieces.
    for (size_t interval = 0; interval + 1 < region_breaks.size(); ++interval) {
        const double left = region_breaks[interval];
        const double right = region_breaks[interval + 1];
        if (right - left <= EPS) continue;

        const double midpoint = 0.5 * (left + right);
        add_candidate(candidates, midpoint);

        double numerator = 0.0;
        double denominator = 0.0;
        bool any_active = false;

        for (size_t i = 0; i < dimension; ++i) {
            const double residual_at_midpoint =
                intercept[i] + slope[i] * midpoint;
            if (abs(residual_at_midpoint) <= config.error_range) {
                continue;
            }

            any_active = true;
            const double sign = residual_at_midpoint > 0.0 ? 1.0 : -1.0;
            const double q_intercept =
                sign * intercept[i] - config.error_range;
            const double q_slope = sign * slope[i];

            numerator += config.weights[i] * q_slope * q_intercept;
            denominator += config.weights[i] * q_slope * q_slope;
        }

        if (any_active && denominator > EPS) {
            const double stationary = -numerator / denominator;
            if (stationary >= left - EPS && stationary <= right + EPS) {
                add_candidate(candidates, stationary);
            }
        }
    }

    sort(candidates.begin(), candidates.end());
    candidates.erase(unique(candidates.begin(), candidates.end(),
                            [](double a, double b) {
                                return abs(a - b) <= EPS;
                            }),
                     candidates.end());

    Result best;
    for (double u : candidates) {
        const double candidate_score =
            score_at(u, current, target, displacement, config);
        if (candidate_score < best.score - EPS ||
            (abs(candidate_score - best.score) <= EPS && u < best.u)) {
            best = {u, candidate_score};
        }
    }
    return best;
}

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

bool parse_int(const char* text, int& value) {
    try {
        size_t used = 0;
        const string input(text);
        value = stoi(input, &used);
        return used == input.size();
    } catch (const exception&) {
        return false;
    }
}

bool parse_double(const char* text, double& value) {
    try {
        size_t used = 0;
        const string input(text);
        value = stod(input, &used);
        return used == input.size() && isfinite(value);
    } catch (const exception&) {
        return false;
    }
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0]
             << " <k:4|5> <error_range> <w1> ... <wD>\n";
        return 1;
    }

    int k = 0;
    if (!parse_int(argv[1], k) || (k != 4 && k != 5)) {
        cerr << "Error: k must be 4 or 5.\n";
        return 1;
    }
    const size_t dimension = rseas_metric::dimension_for_k(k);
    const int expected_argc = static_cast<int>(3 + dimension);
    if (argc != expected_argc) {
        cerr << "Usage: " << argv[0]
             << " <k:4|5> <error_range> <w1> ... <wD>\n"
             << "For k=" << k << ", exactly " << dimension
             << " weights are required.\n";
        return 1;
    }

    rseas_metric::Config config;
    config.weights.assign(dimension, 1.0);
    if (!parse_double(argv[2], config.error_range) ||
        config.error_range < 0.0) {
        cerr << "Error: error_range must be a finite nonnegative number.\n";
        return 1;
    }
    for (size_t i = 0; i < dimension; ++i) {
        if (!parse_double(argv[3 + i], config.weights[i]) ||
            config.weights[i] <= 0.0) {
            cerr << "Error: all " << dimension
                 << " weights must be positive finite numbers.\n";
            return 1;
        }
    }

    Vec current(dimension, 0.0);
    Vec target(dimension, 0.0);
    Vec displacement(dimension, 0.0);

    if (!read_vec(cin, current) ||
        !read_vec(cin, target) ||
        !read_vec(cin, displacement) ||
        !no_extra_values(cin)) {
        cerr << "Error: bxk4one RSEAS expected exactly "
             << (3 * dimension) << " numeric input values for k="
             << k << ".\n";
        return 1;
    }

    for (size_t i = 0; i < dimension; ++i) {
        current[i] = rseas_metric::safe_log(current[i]);
        target[i] = rseas_metric::safe_log(target[i]);

        // rpll writes log(current)-log(stage1). Negating it gives the observed
        // movement from current toward stage1.
        displacement[i] = -displacement[i];
    }

    Result best = find_best_u(current, target, displacement, config);
    if (abs(best.u) <= EPS) best.u = 0.0;

    cout << fixed << setprecision(12) << best.u << '\n';
    return 0;
}
