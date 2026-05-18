#include <bits/stdc++.h>
using namespace std;

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " input.el output.el\n";
        return 1;
    }

    string inputFile = argv[1];
    string outputFile = argv[2];

    ifstream fin(inputFile);
    if (!fin) {
        cerr << "Could not open input file.\n";
        return 1;
    }

    vector<pair<int,int>> edges;
    unordered_set<int> nodeSet;

    int u, v;
    while (fin >> u >> v) {
        if (u == v) continue;
        edges.push_back({u, v});
        nodeSet.insert(u);
        nodeSet.insert(v);
    }
    fin.close();

    vector<int> nodes(nodeSet.begin(), nodeSet.end());
    sort(nodes.begin(), nodes.end());

    unordered_map<int,int> id;
    for (int i = 0; i < (int)nodes.size(); i++) {
        id[nodes[i]] = i;
    }

    int n = nodes.size();
    if (n == 0) {
        cerr << "Empty graph.\n";
        return 1;
    }

    vector<int> deg(n, 0);

    for (auto &e : edges) {
        int a = id[e.first];
        int b = id[e.second];
        deg[a]++;
        deg[b]++;
    }

    for (int i = 0; i < n; i++) {
        if (deg[i] == 0 && n > 1) {
            cerr << "Node " << nodes[i]
                 << " is disconnected.\n";
            return 1;
        }
    }

    mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

    vector<vector<int>> adj(n);
    vector<int> currDeg(n, 0);

    // Build spanning tree
    vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), rng);

    vector<int> connected;
    connected.push_back(perm[0]);

    vector<bool> inTree(n, false);
    inTree[perm[0]] = true;

    for (int idx = 1; idx < n; idx++) {
        int u = perm[idx];

        vector<int> candidates;

        for (int v : connected) {
            if (currDeg[v] < deg[v]) {
                candidates.push_back(v);
            }
        }

        if (candidates.empty()) {
            cerr << "Failed to build spanning tree.\n";
            return 1;
        }

        int v = candidates[rng() % candidates.size()];

        adj[u].push_back(v);
        adj[v].push_back(u);

        currDeg[u]++;
        currDeg[v]++;

        connected.push_back(u);
        inTree[u] = true;
    }

    // Remaining degree
    vector<int> rem(n);

    long long totalRem = 0;

    for (int i = 0; i < n; i++) {
        rem[i] = deg[i] - currDeg[i];

        if (rem[i] < 0) {
            cerr << "Negative remaining degree at node "
                 << nodes[i] << "\n";
            return 1;
        }

        totalRem += rem[i];
    }


    set<pair<int,int>> existingEdges;

    for (int i = 0; i < n; i++) {
        for (int j : adj[i]) {
            if (i < j)
                existingEdges.insert({i, j});
        }
    }

    uniform_real_distribution<double> dist(0.0, 1.0);

	if (totalRem > 0) {
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {

                if (existingEdges.count({i, j}))
                    continue;

                if (rem[i] == 0 || rem[j] == 0)
                    continue;

                double p = (1.0 * rem[i] * rem[j]) / totalRem;
                p = min(1.0, p);

                if (dist(rng) < p) {
                    adj[i].push_back(j);
                    adj[j].push_back(i);

                    existingEdges.insert({i, j});

                    currDeg[i]++;
                    currDeg[j]++;
                }
            }
        }
    }

    // Output graph
    ofstream fout(outputFile);

    if (!fout) {
        cerr << "Could not open output file.\n";
        return 1;
    }

    for (auto &e : existingEdges) {
        fout << nodes[e.first] << " " << nodes[e.second] << "\n";
    }

    fout.close();

    return 0;
}
