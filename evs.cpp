#include <bits/stdc++.h>
using namespace std;

using ld = double;

const int MAX_ITERS = 10000;
const ld EPS = 1e-12;

pair<ld, vector<ld>> evlevc(vector<vector<int>> adj) {

	int n = adj.size();
    // Random initial vector
    vector<ld> x(n);

    mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
    uniform_real_distribution<ld> dist(0.0, 1.0);

    for (int i = 0; i < n; i++)
        x[i] = dist(rng);

    // Normalize initial vector
    ld norm = 0;
    for (ld val : x)
        norm += val * val;

    norm = sqrt(norm);

    for (ld &val : x)
        val /= norm;

    ld prev_lambda = 0;

    for (int iter = 1; iter <= MAX_ITERS; iter++) {
        vector<ld> y(n, 0.0);

        // y = A*x
        for (int i = 0; i < n; i++) {
            for (int nei : adj[i]) {
                y[i] += x[nei];
            }
        }

        // Rayleigh quotient λ = xᵀAx
        ld lambda = 0;
        for (int i = 0; i < n; i++) {
            lambda += x[i] * y[i];
        }

        // Normalize y
        norm = 0;
        for (ld val : y)
            norm += val * val;

        norm = sqrt(norm);

        if (norm < EPS) {
            cerr << "Zero vector encountered.\n";
            return {0, x};
        }

        for (ld &val : y)
            val /= norm;

		if (abs(lambda - prev_lambda) < EPS) {
		    return {lambda, x};
		}

        prev_lambda = lambda;
        x = move(y);
    }

    cout << fixed << setprecision(12);
    cout << "Largest eigenvalue (approx): " << prev_lambda << '\n';
    cout << "Reached max iterations.\n";

    return {prev_lambda, x};
}


int main(int argc, char* argv[]) {

    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " graph.el\n";
        return 1;
    }

    ifstream fin(argv[1]);
    if (!fin) {
        cerr << "Could not open file.\n";
        return 1;
    }

    vector<pair<int,int>> edges;
    unordered_map<int,int> compress;

    int u, v;
    int idx = 0;

    // Read edges + compress arbitrary labels
    while (fin >> u >> v) {
        if (u == v) continue; // ignore self-loops

        if (!compress.count(u))
            compress[u] = idx++;

        if (!compress.count(v))
            compress[v] = idx++;

        edges.push_back({compress[u], compress[v]});
    }

    int n = idx;

    vector<vector<int>> adj(n);

    for (auto &[a, b] : edges) {
        adj[a].push_back(b);
        adj[b].push_back(a); // undirected graph
    }

	auto [evl, evc] = evlevc(adj);
	cout << evl << endl;
}
