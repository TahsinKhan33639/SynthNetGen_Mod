#include <bits/stdc++.h>
using namespace std;

using Vec6 = array<double, 6>;

namespace {

constexpr double EPS = 1e-15;
constexpr double LOG_FLOOR = -15.0;

double safe_log(double value) {
    return value > 0.0 ? log(value) : LOG_FLOOR;
}

bool read_vec6(istream& input, Vec6& values) {
    for (double& value : values) {
        if (!(input >> value)) return false;
    }
    return true;
}

double find_best_multiplier(const Vec6& current,
                            const Vec6& target,
                            const Vec6& displacement) {
    double numerator = 0.0;
    double denominator = 0.0;

    // Minimize sum_i (current_i + u*displacement_i - target_i)^2.
    for (int i = 0; i < 6; ++i) {
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
    if (argc != 1) {
        cerr << "Usage: " << argv[0] << " < scale_input\n";
        return 1;
    }

    Vec6 current{};
    Vec6 target{};
    Vec6 displacement{};

    if (!read_vec6(cin, current) ||
        !read_vec6(cin, target) ||
        !read_vec6(cin, displacement)) {
        cerr << "Error: incomplete bxk4one RMSE input.\n";
        return 1;
    }

    for (int i = 0; i < 6; ++i) {
        current[i] = safe_log(current[i]);
        target[i] = safe_log(target[i]);

        // rpll writes log(current)-log(stage1).  Negating it gives the observed
        // movement from the current graph toward the stage-one graph.
        displacement[i] = -displacement[i];
    }

    double u = find_best_multiplier(current, target, displacement);
    if (abs(u) <= EPS) u = 0.0;

    cout << fixed << setprecision(12) << u << '\n';
    return 0;
}
