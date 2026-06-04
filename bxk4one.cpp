#include <bits/stdc++.h>
using namespace std;

// Finds u minimizing weighted squared error between:
// current + u * disp  and  target
//
// Error = sum_i w[i] * (current[i] + u*disp[i] - target[i])^2
//
// Returns:
//   u
//   weighted RSE = sqrt(weighted squared error / sum(weights))

struct Result {
    double u;
    double rse;
};

Result find_best_u(
    const vector<double>& current,
    const vector<double>& target,
    const vector<double>& disp,
    const vector<int>& weight
) {
    int n = 6;

    double numerator = 0.0;
    double denominator = 0.0;

    // Minimize:
    // sum w_i * (a_i + u d_i - t_i)^2
    //
    // Closed form:
    // u = sum w_i d_i (t_i - a_i) / sum w_i d_i^2

    for (int i = 0; i < n; i++) {
        numerator += weight[i] * disp[i] * (target[i] - current[i]);
        denominator += weight[i] * disp[i] * disp[i];
    }

    double u = 0.0;

    // avoid division by zero because apparently reality enjoys edge cases
    if (fabs(denominator) > 1e-15) {
        u = numerator / denominator;
    }

    double weighted_sq_error = 0.0;
    double weight_sum = 0.0;

    for (int i = 0; i < n; i++) {
        double diff = current[i] + u * disp[i] - target[i];
        weighted_sq_error += weight[i] * diff * diff;
        weight_sum += weight[i];
    }

    double rse = sqrt(weighted_sq_error / weight_sum);

    return {u, rse};
}

int main(int argc, char* argv[]) {
    int cm = stoi(argv[1]);

    ifstream file("weights4.txt");
    if (!file) {
        cerr << "Error: Could not open file.\n";
        return 1;
    }

    vector<vector<int>> all_weights(6, vector<int>(6));

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            file >> all_weights[i][j];
        }
    }

    if (cm < 0 || cm >= 6) {
        cerr << "Error: cm out of range.\n";
        return 1;
    }

    vector<double> current(6), target(6), disp(6);
	vector<int> weight(6);

    for (int i = 0; i < 6; i++) {cin >> current[i]; if (current[i] == 0) current[i] = -15; else current[i] = log(current[i]);}
    for (int i = 0; i < 6; i++) {cin >> target[i]; if (target[i] == 0) target[i] = -15; else target[i] = log(target[i]);}
    for (int i = 0; i < 6; i++) {cin >> disp[i]; disp[i] = -disp[i];}

    weight = all_weights[cm];

    Result ans = find_best_u(current, target, disp, weight);
	ans.u = clamp(ans.u, 0.0, 2.0);

    cout << fixed << setprecision(12);
    cout << ans.u << endl;
}
