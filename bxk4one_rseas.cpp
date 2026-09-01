#include <bits/stdc++.h>
#include "rseas_metric.hpp"
using namespace std;

using Vec6 = rseas_metric::Vec6;

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
                const Vec6& current,
                const Vec6& target,
                const Vec6& displacement,
                const rseas_metric::Config& config) {
    Vec6 signed_errors{};
    for (int i = 0; i < 6; ++i) {
        const double predicted = current[i] + u * displacement[i];
        signed_errors[i] = predicted - target[i];
    }
    return rseas_metric::score(signed_errors, config);
}

Result find_best_u(const Vec6& current,
                   const Vec6& target,
                   const Vec6& displacement,
                   const rseas_metric::Config& config) {
    Vec6 intercept{};
    Vec6 slope{};
    for (int i = 0; i < 6; ++i) {
        intercept[i] = current[i] - target[i];
        slope[i] = displacement[i];
    }

    vector<double> region_breaks{U_MIN, U_MAX};
    vector<double> candidates{U_MIN, U_MAX};

    // The active set of squared-hinge terms changes only where a residual is
    // +e or -e. Between consecutive points the objective is one quadratic.
    for (int i = 0; i < 6; ++i) {
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

    // Inside the error range, RSEAS is max_i |residual_i| - e. That is the
    // upper envelope of the twelve signed residual lines, whose minimum occurs
    // at an endpoint or an intersection of two such lines.
    struct Line {
        double intercept;
        double slope;
    };
    vector<Line> lines;
    lines.reserve(12);
    for (int i = 0; i < 6; ++i) {
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

    // In each outside interval, find the exact stationary point of the active
    // weighted quadratic. Midpoints are also candidates for constant pieces.
    for (size_t interval = 0; interval + 1 < region_breaks.size(); ++interval) {
        const double left = region_breaks[interval];
        const double right = region_breaks[interval + 1];
        if (right - left <= EPS) continue;

        const double midpoint = 0.5 * (left + right);
        add_candidate(candidates, midpoint);

        double numerator = 0.0;
        double denominator = 0.0;
        bool any_active = false;

        for (int i = 0; i < 6; ++i) {
            const double residual_at_midpoint =
                intercept[i] + slope[i] * midpoint;
            if (abs(residual_at_midpoint) <= config.error_range) {
                continue;
            }

            any_active = true;
            const double sign = residual_at_midpoint > 0.0 ? 1.0 : -1.0;
            const double q_intercept = sign * intercept[i] - config.error_range;
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

bool read_vec6(istream& input, Vec6& values) {
    for (double& value : values) {
        if (!(input >> value)) return false;
    }
    return true;
}

bool parse_config(char* argv[], rseas_metric::Config& config) {
    try {
        config.error_range = stod(argv[1]);
        for (int i = 0; i < 6; ++i) {
            config.weights[i] = stod(argv[2 + i]);
        }
    } catch (const exception&) {
        return false;
    }

    if (!isfinite(config.error_range) || config.error_range < 0.0) {
        return false;
    }
    for (double weight : config.weights) {
        if (!isfinite(weight) || weight <= 0.0) {
            return false;
        }
    }
    return true;
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc != 8) {
        cerr << "Usage: " << argv[0]
             << " <error_range> <w1> <w2> <w3> <w4> <w5> <w6>\n";
        return 1;
    }

    rseas_metric::Config config;
    if (!parse_config(argv, config)) {
        cerr << "Error: error_range must be nonnegative and all six weights "
                "must be positive finite numbers.\n";
        return 1;
    }

    Vec6 current{};
    Vec6 target{};
    Vec6 displacement{};

    if (!read_vec6(cin, current) ||
        !read_vec6(cin, target) ||
        !read_vec6(cin, displacement)) {
        cerr << "Error: incomplete bxk4one input.\n";
        return 1;
    }

    for (int i = 0; i < 6; ++i) {
        current[i] = rseas_metric::safe_log(current[i]);
        target[i] = rseas_metric::safe_log(target[i]);

        // rpll writes log(current)-log(stage1). Negating gives the observed
        // movement from current toward stage1.
        displacement[i] = -displacement[i];
    }

    Result best = find_best_u(current, target, displacement, config);
    if (abs(best.u) <= EPS) best.u = 0.0;

    cout << fixed << setprecision(12) << best.u << '\n';
    return 0;
}
