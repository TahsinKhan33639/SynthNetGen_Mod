#include <bits/stdc++.h>
using namespace std;

int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " input.el output.el\n";
        return 1;
    }

    string inputFile = argv[1];
    string outputFile = argv[2];

    ifstream in(inputFile);
    ofstream out(outputFile);

    if (!in.is_open()) {
        cerr << "cannot open input file " << inputFile << "\n";
        return 1;
    }
    if (!out.is_open()) {
        cerr << "cannot open output file " << outputFile << "\n";
        return 1;
    }

    unordered_map<string, int> id_map;
    int id_counter = 0;

    set<pair<int,int>> seen_edges;
    vector<pair<int,int>> edges;

    string u, v;
    while (in >> u >> v) {
        if (u == v) continue;

        if (!id_map.count(u)) id_map[u] = id_counter++;
        if (!id_map.count(v)) id_map[v] = id_counter++;

        int a = id_map[u];
        int b = id_map[v];

        if (a > b) swap(a, b);

        if (!seen_edges.count({a, b})) {
            seen_edges.insert({a, b});
            edges.push_back({a, b});
        }
    }

    int n = id_counter;

    vector<vector<int>> adj(n);
    for (auto &[a, b] : edges) {
        adj[a].push_back(b);
        adj[b].push_back(a);
    }

    vector<int> comp(n, -1);
    vector<int> comp_size;
    int comp_id = 0;

    for (int i = 0; i < n; i++) {
        if (comp[i] != -1) continue;

        queue<int> q;
        q.push(i);
        comp[i] = comp_id;
        int sz = 0;

        while (!q.empty()) {
            int cur = q.front();
            q.pop();
            sz++;

            for (int nei : adj[cur]) {
                if (comp[nei] == -1) {
                    comp[nei] = comp_id;
                    q.push(nei);
                }
            }
        }

        comp_size.push_back(sz);
        comp_id++;
    }

    int largest_comp = max_element(comp_size.begin(), comp_size.end()) - comp_size.begin();

    unordered_map<int,int> new_id;
    int new_counter = 0;

    for (int i = 0; i < n; i++) {
        if (comp[i] == largest_comp) {
            new_id[i] = new_counter++;
        }
    }

    int output_edges = 0;
    for (auto &[a, b] : edges) {
        if (comp[a] == largest_comp && comp[b] == largest_comp) {
            out << new_id[a] << " " << new_id[b] << "\n";
            output_edges++;
        }
    }

    return 0;
}
