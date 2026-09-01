#include <bits/stdc++.h>
using namespace std;

using Vec6 = array<double, 6>;

namespace {

constexpr double EPS = 1e-12;
constexpr double LOG_FLOOR = -15.0;

struct Curve {
    double ka = 0.0;
    double kb = 0.0;
    Vec6 A{};
    Vec6 B{};
};

Vec6 lerp(const Vec6& p, const Vec6& q, double u) {
    Vec6 result{};
    for (int i = 0; i < 6; ++i) {
        result[i] = (1.0 - u) * p[i] + u * q[i];
    }
    return result;
}

Vec6 average(const Vec6& a, const Vec6& b) {
    Vec6 result{};
    for (int i = 0; i < 6; ++i) {
        result[i] = 0.5 * (a[i] + b[i]);
    }
    return result;
}

Curve normalize_curve(Curve curve) {
    curve.ka = max(0.0, curve.ka);
    curve.kb = max(0.0, curve.kb);

    // Keep each measured displacement attached to the scale that produced it.
    if (curve.kb < curve.ka) {
        swap(curve.ka, curve.kb);
        swap(curve.A, curve.B);
    }

    // Integer rounding can map both probes to the same effective scale.
    if (curve.ka > 0.0 && abs(curve.kb - curve.ka) <= EPS) {
        curve.A = average(curve.A, curve.B);
        curve.B = curve.A;
    }

    return curve;
}

Vec6 displacement(double x, const Curve& curve) {
    const Vec6 zero{};
    if (x <= 0.0 || curve.kb <= 0.0) {
        return zero;
    }

    if (curve.ka <= 0.0) {
        return lerp(zero, curve.B, x / curve.kb);
    }

    if (x <= curve.ka) {
        return lerp(zero, curve.A, x / curve.ka);
    }

    if (curve.kb - curve.ka <= EPS) {
        return lerp(zero, curve.A, x / curve.ka);
    }

    // u > 1 deliberately extrapolates beyond the larger tested scale.
    const double u = (x - curve.ka) / (curve.kb - curve.ka);
    return lerp(curve.A, curve.B, u);
}

void add_range(vector<double>& grid, double first, double last, double step) {
    if (step <= 0.0 || last < first) return;
    for (double x = first; x <= last + EPS; x += step) {
        grid.push_back(x);
    }
}

vector<double> build_grid(const Curve& first, const Curve& second) {
    vector<double> grid;

    const double upper =
        2.0 * min(max(first.ka, first.kb),
                  max(second.ka, second.kb));

    add_range(grid, 0.0, 600.0, 24.0);
    add_range(grid, 600.0, 6000.0, 240.0);
    add_range(grid, 6000.0, 60000.0, 2400.0);

    if (upper > 60000.0) {
        add_range(grid, 60000.0, upper, 12000.0);
    }

    grid.push_back(first.ka);
    grid.push_back(first.kb);
    grid.push_back(second.ka);
    grid.push_back(second.kb);

    grid.erase(remove_if(grid.begin(), grid.end(),
                         [&](double x) {
                             return x < -EPS || x > upper + EPS;
                         }),
               grid.end());

    sort(grid.begin(), grid.end());
    grid.erase(unique(grid.begin(), grid.end(), [](double a, double b) {
        return abs(a - b) <= EPS;
    }), grid.end());

    return grid;
}

// This deliberately matches the original rpll score convention for six unit
// weights: sqrt(sum_i error_i^2).  Dividing by six would only multiply every
// candidate score by the same constant and would not change the minimizer.
double rmse_at(const Vec6& target,
               const Vec6& initial,
               double x,
               double y,
               const Curve& first,
               const Curve& second) {
    const Vec6 v1 = displacement(x, first);
    const Vec6 v2 = displacement(y, second);

    double sum = 0.0;
    for (int i = 0; i < 6; ++i) {
        const double predicted = initial[i] - v1[i] - v2[i];
        const double error = predicted - target[i];
        sum += error * error;
    }
    return sqrt(sum);
}

struct SearchResult {
    double x = 0.0;
    double y = 0.0;
    double score = numeric_limits<double>::infinity();
};

SearchResult minimize_piece(const Vec6& target,
                            const Vec6& initial,
                            const Curve& first,
                            const Curve& second,
                            int edge_change1,
                            int edge_change2) {
    SearchResult best;
    best.score = rmse_at(target, initial, 0.0, 0.0, first, second);

    const vector<double> grid = build_grid(first, second);

    if (edge_change1 == 0 && edge_change2 == 0) {
        for (double x : grid) {
            for (double y : grid) {
                const double score =
                    rmse_at(target, initial, x, y, first, second);
                if (score < best.score) {
                    best = {x, y, score};
                }
            }
        }
    } else if ((edge_change1 < 0 && edge_change2 > 0) ||
               (edge_change1 > 0 && edge_change2 < 0)) {
        for (double x : grid) {
            const double y = abs(
                x * static_cast<double>(edge_change1) /
                static_cast<double>(edge_change2));
            const double score =
                rmse_at(target, initial, x, y, first, second);
            if (score < best.score) {
                best = {x, y, score};
            }
        }
    }

    // Preserve the intentional rule from the existing code: when exactly one
    // edge change is zero, ignore the pair and return the no-change scales.
    return best;
}

bool read_vec6(istream& input, Vec6& values) {
    for (double& value : values) {
        if (!(input >> value)) return false;
    }
    return true;
}

double safe_log(double value) {
    return value > 0.0 ? log(value) : LOG_FLOOR;
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0]
             << " <edge_change_1> <edge_change_2>\n";
        return 1;
    }

    int edge_change1;
    int edge_change2;
    try {
        edge_change1 = stoi(argv[1]);
        edge_change2 = stoi(argv[2]);
    } catch (const exception&) {
        cerr << "Error: invalid edge-change argument.\n";
        return 1;
    }

    Vec6 initial{};
    Vec6 target{};
    Curve first;
    Curve second;

    // Input layout, produced by rpll_rseas_then_rmse.sh:
    //   initial[6], target[6]
    //   first.ka first.kb, first.A[6], first.B[6]
    //   second.ka second.kb, second.A[6], second.B[6]
    if (!read_vec6(cin, initial) || !read_vec6(cin, target) ||
        !(cin >> first.ka >> first.kb) ||
        !read_vec6(cin, first.A) || !read_vec6(cin, first.B) ||
        !(cin >> second.ka >> second.kb) ||
        !read_vec6(cin, second.A) || !read_vec6(cin, second.B)) {
        cerr << "Error: incomplete bxk4f RMSE input.\n";
        return 1;
    }

    for (int i = 0; i < 6; ++i) {
        initial[i] = safe_log(initial[i]);
        target[i] = safe_log(target[i]);
    }

    first = normalize_curve(first);
    second = normalize_curve(second);

    const SearchResult best = minimize_piece(
        target, initial, first, second, edge_change1, edge_change2);

    cout << fixed << setprecision(0) << best.x << ' ' << best.y << ' ';
    cout << setprecision(17) << best.score << '\n';
    return 0;
}
