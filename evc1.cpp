#include <bits/stdc++.h>
using namespace std;

vector<int> deg(20000);

bool cmp(int a, int b) {
    return deg[a] < deg[b];
}

using VP = pair<vector<int>, vector<int>>;
using RebuildResult = pair<VP, pair<VP, VP>>;

RebuildResult rebuildSets(const vector<int>& A, const vector<int>& B) {
    unordered_set<int> setA(A.begin(), A.end());
    unordered_set<int> setB(B.begin(), B.end());

    vector<int> common;
    for (int x : A) {
        if (setB.count(x)) common.push_back(x);
    }

    vector<int> onlyA;
    for (int x : A) {
        if (!setB.count(x)) onlyA.push_back(x);
    }

    vector<int> onlyB;
    for (int x : B) {
        if (!setA.count(x)) onlyB.push_back(x);
    }

    int cntA = (int)onlyA.size();

    vector<int> combined = onlyA;
    combined.insert(combined.end(), onlyB.begin(), onlyB.end());
    sort(combined.begin(), combined.end(), cmp);

    vector<int> newOnlyA;
    vector<int> newOnlyB;
    int total = (int)combined.size();

    newOnlyA.insert(newOnlyA.end(), combined.end() - cntA, combined.end());
    newOnlyB.insert(newOnlyB.end(), combined.begin(), combined.end() - cntA);

    vector<int> newA = newOnlyA;
    vector<int> newB = newOnlyB;
    newA.insert(newA.end(), common.begin(), common.end());
    newB.insert(newB.end(), common.begin(), common.end());

    return {{newA, newB}, {{onlyA, onlyB}, {newOnlyA, newOnlyB}}};
}

void sample_edge_prob(vector<vector<int>>& adj, int num_samples) {
    int n = (int)adj.size();
    if (n < 4) {
        cerr << "Graph must have at least 4 nodes.\n";
        return;
    }

    random_device rd;
    mt19937 gen(rd());
    vector<bool> visited(adj.size(), false);

    for (int sample = 0; sample < num_samples; sample++) {
        vector<int> nodes;
        uniform_int_distribution<int> node_dist(0, n - 1);
        uniform_int_distribution<int> prob_dist(0, 100);

        int a = node_dist(gen);
        uniform_int_distribution<int> neighbor_dist(0, (int)adj[a].size() - 1);
        int b = adj[a][neighbor_dist(gen)];

        adj[a].erase(remove(adj[a].begin(), adj[a].end(), b), adj[a].end());
        adj[b].erase(remove(adj[b].begin(), adj[b].end(), a), adj[b].end());

        if (adj[a].size() < adj[b].size()) {
            swap(a, b);
        }

        auto rebuilt = rebuildSets(adj[a], adj[b]);
        auto newSets = rebuilt.first;
        auto changes = rebuilt.second;

        vector<int> newA = newSets.first;
        vector<int> newB = newSets.second;
        vector<int> onlyA = changes.first.first;
        vector<int> onlyB = changes.first.second;
        vector<int> newOnlyA = changes.second.first;
        vector<int> newOnlyB = changes.second.second;

        for (int x : onlyA) {
            adj[x].erase(remove(adj[x].begin(), adj[x].end(), a), adj[x].end());
        }
        for (int x : onlyB) {
            adj[x].erase(remove(adj[x].begin(), adj[x].end(), b), adj[x].end());
        }
        for (int x : newOnlyA) {
            adj[x].push_back(a);
        }
        for (int x : newOnlyB) {
            adj[x].push_back(b);
        }

        adj[a] = newA;
        adj[b] = newB;
        adj[a].push_back(b);
        adj[b].push_back(a);
    }
}

void readGraph(const string& filename, vector<vector<int>>& adj) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: cannot open file '" << filename << "'\n";
        exit(1);
    }

    map<string, int> nodeMap;
    vector<pair<int, int>> edges;
    string line, u, v;
    int nodeCount = 0;

    while (getline(file, line)) {
        istringstream iss(line);
        if (!(iss >> u >> v)) continue;

        if (!nodeMap.count(u)) nodeMap[u] = nodeCount++;
        if (!nodeMap.count(v)) nodeMap[v] = nodeCount++;

        int a = nodeMap[u];
        int b = nodeMap[v];
        edges.push_back({a, b});
    }

    adj.assign(nodeCount, {});
    for (auto [a, b] : edges) {
        adj[a].push_back(b);
        adj[b].push_back(a);
    }
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <input_file1> <num_samples>\n";
        return 1;
    }

    vector<vector<int>> adj;
    readGraph(argv[1], adj);

    for (int i = 0; i < (int)adj.size(); i++) {
        deg[i] = (int)adj[i].size();
    }

    int num_samples = stoi(argv[2]);

    for (int i = 0; i < 1; i++) {
        sample_edge_prob(adj, num_samples);
    }

    for (int i = 0; i < (int)adj.size(); i++) {
        for (int v : adj[i]) {
            if (i < v) cout << i << " " << v << "\n";
        }
    }
    return 0;
}
