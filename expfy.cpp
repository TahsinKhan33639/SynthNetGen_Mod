#include <chrono>
#include <bits/stdc++.h>
using namespace std;
using ld = double;
const int MAX_ITERS = 10000;
const ld EPS = 1e-12;

// Bound work by attempted graphlet samples rather than elapsed wall time.
// This makes each expfy process do a comparable amount of work even when
// several processes are competing for CPU time in parallel.
const long long MIN_SAMPLE_ATTEMPTS = 800000;
const long long MAX_SAMPLE_ATTEMPTS = 1000000;
const long long ATTEMPTS_PER_REQUESTED_CHANGE = 500;

// A single graphlet-sampling attempt must also be bounded internally.
// Otherwise pick_neighbor could keep retrying forever before the outer sample
// counter gets a chance to advance.
const int MAX_NEIGHBOR_PICK_ATTEMPTS = 100;

struct SampleStats {
    int done1 = 0;
    int done2 = 0;
    long long attempts = 0;
};

unordered_map<int, int> adj6bit_to_id = {
    {0b000000, 0},
    {0b000001, 1},
    {0b000010, 1},
    {0b000011, 3},
    {0b000100, 1},
    {0b000101, 3},
    {0b000110, 3},
    {0b000111, 4},
    {0b001000, 1},
    {0b001001, 3},
    {0b001010, 3},
    {0b001011, 5},
    {0b001100, 2},
    {0b001101, 6},
    {0b001110, 6},
    {0b001111, 7},
    {0b010000, 1},
    {0b010001, 3},
    {0b010010, 2},
    {0b010011, 6},
    {0b010100, 3},
    {0b010101, 5},
    {0b010110, 6},
    {0b010111, 7},
    {0b011000, 3},
    {0b011001, 4},
    {0b011010, 6},
    {0b011011, 7},
    {0b011100, 6},
    {0b011101, 7},
    {0b011110, 8},
    {0b011111, 9},
    {0b100000, 1},
    {0b100001, 2},
    {0b100010, 3},
    {0b100011, 6},
    {0b100100, 3},
    {0b100101, 6},
    {0b100110, 5},
    {0b100111, 7},
    {0b101000, 3},
    {0b101001, 6},
    {0b101010, 4},
    {0b101011, 7},
    {0b101100, 6},
    {0b101101, 8},
    {0b101110, 7},
    {0b101111, 9},
    {0b110000, 3},
    {0b110001, 6},
    {0b110010, 6},
    {0b110011, 8},
    {0b110100, 4},
    {0b110101, 7},
    {0b110110, 7},
    {0b110111, 9},
    {0b111000, 5},
    {0b111001, 7},
    {0b111010, 7},
    {0b111011, 9},
    {0b111100, 7},
    {0b111101, 9},
    {0b111110, 9},
    {0b111111, 10},
};

int get_graph_id(const vector<vector<int>>& adj) {
    int bits = 0;
    bits |= adj[2][3] ? (1 << 0) : 0;
    bits |= adj[1][3] ? (1 << 1) : 0;
    bits |= adj[1][2] ? (1 << 2) : 0;
    bits |= adj[0][3] ? (1 << 3) : 0;
    bits |= adj[0][2] ? (1 << 4) : 0;
    bits |= adj[0][1] ? (1 << 5) : 0;
    if (adj6bit_to_id.count(bits))
        return adj6bit_to_id[bits];
    else
        return -1;
}

SampleStats sample_edge_prob(vector<vector<int>>& adj, vector<vector<char>>& bvv, int num_samples, int num_samples2, int m1, int m2,
					 vector<ld> x, ld evl, ld evbound, vector<int> degb, double degbound, vector<int> gdnd, int hdgbound) {
	vector<double> cntr(6, 0);
    int n = adj.size();
	vector<bool> gnif(n, false);
	for (auto i : gdnd) gnif[i] = true;
    if (n < 4) {
        cerr << "Graph must have at least 4 nodes.\n";
        return {};
    }
    random_device rd;
    mt19937 rng(rd());

	vector<ld> y(n, 0), z(n, 0);
	for(int i = 0; i<n; i++){
		for(auto j : adj[i]){
			y[i] += x[j];
		}
		y[i] = 1/evl*y[i];
	}
	int edgec = 0, edgeb = 0;
	for(auto de : degb) edgeb += de;
	edgec = edgeb;
	vector<int> deg(n, 0), degz(n, 0);
	for(int i = 0; i<n; i++){
		deg[i] = adj[i].size();
	}

	vector<int> hdg(n, 0), hdgz(n, 0), hdgb(n, 0);
	for(int i=0; i<n; i++){
		hdg[deg[i]]++;
		hdgb[degb[i]]++;
	}

	auto j = [&](int u, int v) {
	    adj[v].push_back(u);
	    adj[u].push_back(v);
	    bvv[v][u] = bvv[u][v] = 1;
		y[u] = y[u] + 1/evl*x[v];
		y[v] = y[v] + 1/evl*x[u];
		hdg[deg[u]]--; hdg[deg[v]]--;
		deg[u]++; deg[v]++;
		hdg[deg[u]]++; hdg[deg[v]]++;
		edgec++;
	};

	auto d = [&](int u, int v) {
	    adj[v].erase(remove(adj[v].begin(), adj[v].end(), u), adj[v].end());
	    adj[u].erase(remove(adj[u].begin(), adj[u].end(), v), adj[u].end());
	    bvv[v][u] = bvv[u][v] = 0;
		y[u] = y[u] - 1/evl*x[v];
		y[v] = y[v] - 1/evl*x[u];
		hdg[deg[u]]--; hdg[deg[v]]--;
		deg[u]--; deg[v]--;
		hdg[deg[u]]++; hdg[deg[v]]++;
		edgec--;
	};

	int donereps = 0;
	int donereps2 = 0;
	int active_slot = 1;
	int m = m1;

	int lklk1=0, lklk2=0, lklk3=0, lklk4=0;

    auto op = [&](const vector<tuple<int,int,int>>& ops) {
        if (ops.empty()) return;

        auto valid = [&](int node, ld new_val) {
            ld old_val = y[node];
            if (new_val <= 0 || old_val <= 0 || x[node] <= 0)
                return false;

            ld target = log(x[node]);

            ld new_dist = abs(log(new_val) - target);
            if (new_dist <= evbound)
                return true;

            ld old_dist = abs(log(old_val) - target);
            if (new_dist < old_dist)
                return true;
            lklk1++;
            return false;
        };

        auto valideg = [&](int node, int new_val) {
            int old_val = deg[node];
            if (new_val <= 0 || old_val <= 0 || x[node] <= 0)
                return false;

            int target = degb[node];

            int new_dist = abs(new_val - target);
            double new_distl = abs(log(new_val) - log(target));
            if (new_distl <= degbound)
                return true;

			if (new_dist <= max(degbound*n*3/100, (double)5))
				return true;

            int old_dist = abs(old_val - target);
            if (new_dist < old_dist)
                return true;
            lklk2++;
            return false;
        };

        auto validhdg = [&](int dgX, int new_val) {
            if (dgX < 0 || dgX >= n)
                return false;

            int old_val = hdg[dgX];
            // Zero nodes in a degree bin is valid. Only negative counts are invalid.
            if (new_val < 0 || old_val < 0)
                return false;

            int target = hdgb[dgX];

            int new_dist = abs(new_val - target);
            if (new_dist <= hdgbound)
                return true;

            int old_dist = abs(old_val - target);
            if (new_dist < old_dist)
                return true;

            lklk3++;
            return false;
        };

        set<int> nodesinop;
        set<pair<int,int>> touched_edges;

        // Validate the complete operation list against the current simple graph
        // before changing temporary degrees or indexing histogram arrays.
        for (auto [type, n1, n2] : ops) {
            if (type != 1 && type != 2)
                return;
            if (n1 < 0 || n1 >= n || n2 < 0 || n2 >= n || n1 == n2)
                return;

            int a = min(n1, n2);
            int b = max(n1, n2);
            if (!touched_edges.insert({a, b}).second)
                return;

            bool edge_exists = bvv[n1][n2] != 0;
            if ((type == 1 && edge_exists) ||
                (type == 2 && !edge_exists)) {
                return;
            }

            nodesinop.insert(n1);
            nodesinop.insert(n2);
        }

        for (int node : nodesinop) {
            z[node] = y[node];
            degz[node] = deg[node];
        }

        for (auto [type, n1, n2] : ops) {
            if (type == 1) {
                z[n1] += 1 / evl * x[n2];
                z[n2] += 1 / evl * x[n1];
                degz[n1]++;
                degz[n2]++;
            } else {
                z[n1] -= 1 / evl * x[n2];
                z[n2] -= 1 / evl * x[n1];
                degz[n1]--;
                degz[n2]--;
            }
        }

        // A simple n-node graph has degrees in [0, n-1]. Reject an operation
        // before any degree is used as a histogram index if it violates that.
        for (int node : nodesinop) {
            if (deg[node] < 0 || deg[node] >= n ||
                degz[node] < 0 || degz[node] >= n) {
                return;
            }
        }

        bool open = true;
        for (int node : nodesinop) {
            open &= valid(node, z[node]);
            open &= valideg(node, degz[node]);
        }

        set<int> related_deg;
        for (int node : nodesinop) {
            related_deg.insert(deg[node]);
            related_deg.insert(degz[node]);
        }

        for (int degree : related_deg) {
            hdgz[degree] = hdg[degree];
        }
        for (int node : nodesinop) {
            hdgz[deg[node]]--;
            hdgz[degz[node]]++;
        }
        for (int degree : related_deg) {
            open &= validhdg(degree, hdgz[degree]);
        }

        if (!open) return;
        lklk4++;
        if (active_slot == 1 && donereps < num_samples) donereps++;
        if (active_slot == 2 && donereps2 < num_samples2) donereps2++;
        for (auto [type, n1, n2] : ops) {
            if (type == 1) j(n1, n2);
            if (type == 2) d(n1, n2);
        }
    };

	num_samples = max(0, num_samples);
	num_samples2 = max(0, num_samples2);

	const long long requested_changes =
	    static_cast<long long>(num_samples) + num_samples2;
	long long max_sample_attempts = 0;
	if (requested_changes > 0) {
		const long long scaled_attempt_limit =
		    requested_changes > MAX_SAMPLE_ATTEMPTS / ATTEMPTS_PER_REQUESTED_CHANGE
		        ? MAX_SAMPLE_ATTEMPTS
		        : requested_changes * ATTEMPTS_PER_REQUESTED_CHANGE;
		max_sample_attempts = clamp(
		    scaled_attempt_limit, MIN_SAMPLE_ATTEMPTS, MAX_SAMPLE_ATTEMPTS);
	}

	long long smpl = 0;
	for (; smpl < max_sample_attempts &&
	       ((donereps < num_samples) || (donereps2 < num_samples2));
	     ++smpl) {

		const bool need1 = donereps < num_samples;
		const bool need2 = donereps2 < num_samples2;

		if (!need1) {
			active_slot = 2;
		} else if (!need2) {
			active_slot = 1;
		} else {
			const double progress1 = static_cast<double>(donereps) / num_samples;
			const double progress2 = static_cast<double>(donereps2) / num_samples2;
			if (progress1 > progress2 + 0.02) active_slot = 2;
			if (progress1 < progress2 - 0.02) active_slot = 1;
		}
		m = active_slot == 1 ? m1 : m2;

        vector<int> nodes;

		uniform_int_distribution<int> dist0(0, n - 1);
		uniform_real_distribution<double> prob(0.0, 1.0);

		int node1 = dist0(rng);

		while (!gnif[node1]) {
			if (prob(rng) < 0.2) break;
			node1 = dist0(rng);
		}
        nodes.push_back(node1);

		vector<char> blocked(adj.size(), false);
		for (int x : nodes) blocked[x] = true;
		blocked[node1] = true;
	    vector<int> options;
	    options.reserve(32);

        auto pick_neighbor = [&](const vector<int>& cds) -> int {
            for (int attempt = 0; attempt < MAX_NEIGHBOR_PICK_ATTEMPTS;
                 ++attempt) {
                int u = cds[rng() % cds.size()];
                if (adj[u].empty()) continue;

                int v = adj[u][rng() % adj[u].size()];
                if (find(cds.begin(), cds.end(), v) != cds.end()) continue;

                // Preserve the old 70% rejection bias for nodes not listed in
                // data_middle1.txt, but count the rejection as an attempt.
                if (!gnif[v] && prob(rng) < 0.7) continue;

                return v;
            }
            return -1;
        };

        int node2 = pick_neighbor({nodes[0]});
        if (node2 == -1) continue;
        nodes.push_back(node2);

        int node3 = pick_neighbor(nodes);
        if (node3 == -1) continue;
        nodes.push_back(node3);

        int node4 = pick_neighbor(nodes);
        if (node4 == -1) continue;
        nodes.push_back(node4);

        shuffle (nodes.begin(), nodes.end(), rng);

        vector<vector<int>> sub(4, vector<int>(4, 0));
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < a; b++) {
                int u = nodes[a], v = nodes[b];
                sub[a][b] = sub[b][a] = bvv[u][v];
            }
        }
        int gid = get_graph_id(sub);
		cntr[gid-5]++;

		if (m == 60) {
	        int node5 = pick_neighbor(nodes);
	        if (node5 == -1) continue;
    	    nodes.push_back(node5);

	        vector<vector<int>> sub5(5, vector<int>(5, 0));
	        for (int a = 0; a < 5; a++) {
	            for (int b = 0; b < a; b++) {
	                int u = nodes[a], v = nodes[b];
	                sub5[a][b] = sub5[b][a] = bvv[u][v];
	            }
	        }

			if (sub5[0][1]+sub5[0][2]+sub5[0][3]+sub5[0][4]+sub5[1][2]+sub5[1][3]+sub5[1][4]+sub5[2][3]+sub5[2][4]+sub5[3][4] == 9) {
				int cpt1 = 1, cpt2 = 0;
				for (int cptl1 = 0; cptl1 < 5; cptl1++){
					for (int cptl2 = 0; cptl2 < cptl1; cptl2++){
						if (sub5[cptl1][cptl2] == 0){
							cpt1 = cptl1; cpt2 = cptl2;
						}
					}
				}
				op({{1,nodes[cpt1],nodes[cpt2]}});
			}
		}

		if (gid == 10 && ((43 <= m && m <= 47) || m == 59)){
			int u = -1, v = -1, w = -1, z = -1;
			u = nodes[0];
			v = nodes[1];
			w = nodes[2];
			z = nodes[3];
			if (m==43) { op({{2,u,v}}); }
			if (m==44) { op({{2,u,v}, {2,w,z}}); }
			if (m==45) { op({{2,u,z}, {2,v,z}}); }
			if (m==46) { op({{2,u,v}, {2,w,z}, {2,u,z}}); }
			if (m==47) { op({{2,u,v}, {2,u,w}, {2,v,w}}); }
			if (m==59) {
				uniform_int_distribution<int> dist(0, adj[u].size() - 1);
				int uuu = adj[u][dist(rng)];
				op({{1,uuu,v}, {1,uuu,w}, {1,uuu,z}});
			}
		}
		if (gid == 9 && (m == 24 || m == 31 || m == 55 || (48 <= m && m <= 49))){
			int u = -1, v = -1, w = -1, z = -1;
			vector<int> blanks;
			for (int i = 0; i < 4; i++){
				if (sub[i][0] + sub[i][1] + sub[i][2] + sub[i][3] == 2) blanks.push_back(i);
			}
			u = nodes[blanks[0]];
			v = nodes[blanks[1]];
			for (int i = 0; i < 4; i++){
				if (sub[i][0] + sub[i][1] + sub[i][2] + sub[i][3] == 3) blanks.push_back(i);
			}
			w = nodes[blanks[2]];
			z = nodes[blanks[3]];
			if (m==24) { op({{1,u,v}, {2,u,z}}); }
			if (m==31) { op({{1,u,v}}); }
			if (m==48) { op({{2,z,w}}); }
			if (m==49) { op({{2,u,w}}); }
			if (m==55) {
			    uniform_int_distribution<int> dist(0, adj[u].size() - 1);
				int uuu = adj[u][dist(rng)];
				if (deg[uuu] != 1) op({{1,u,v}, {2,u,uuu}});
			}
		}
        if (gid == 8 && ( (21 <= m && m <= 23) || m == 27 || m == 40 || m == 32 || m == 50 || m == 57 || m == 58)){
			int u = -1, v = -1, w = -1, z = -1;
            vector<int> dig;
            w = nodes[0];
            for (int i = 1; i < 4; i++){
                if (sub[0][i] == 1) dig.push_back(i);
                else z = nodes[i];
            }
            v = nodes[dig[0]];
            u = nodes[dig[1]];
			if (m==21) { op({{1,u,v}, {1,w,z}, {2,u,z}, {2,w,v}}); }
			if (m==22) { op({{1,v,u}, {1,z,w}, {2,z,v}, {2,w,v}}); }
			if (m==23) { op({{1,u,v}, {2,v,w}}); }
			if (m==32) { op({{1,u,v}}); }
			if (m==27) { op({{1,u,v}, {1,z,w}}); }
			if (m==40) { op({{2,u,z}}); }
			if (m==50) { op({{1,w,z}, {2,v,w}, {2,u,w}}); }
            if (m==57) {
                uniform_int_distribution<int> dist(0, adj[w].size() - 1);
                int ww1 = adj[w][dist(rng)];
                uniform_int_distribution<int> dist2(0, adj[z].size() - 1);
                int zz1 = adj[z][dist2(rng)];
                if (ww1 != zz1 && deg[ww1] != 1 && deg[zz1] != 1) op({{1,z,w}, {1,u,v}, {2,w,ww1}, {2,z,zz1}});
            }
			if (m==58) {
				op({{2,u,z}, {2,u,w}});
			}
		}
        if (gid == 7 && ( (10 <= m && m <= 17) || m == 28 || m == 33 || m == 41 || m == 42 || m == 51 || m == 56)){
			int u = -1, v = -1, w = -1, z = -1;
			vector<int> subsum(4, 0), par;
 			for (int i = 0; i < 4; i++){
				subsum[i] = sub[i][0] + sub[i][1] + sub[i][2] + sub[i][3];
                if (subsum[i] == 1) w = nodes[i];
                else if (subsum[i] == 3) z = nodes[i];
				else par.push_back(i);
            }
			u = nodes[par[0]];
			v = nodes[par[1]];
			if (m==10) { op({{1,w,v}, {2,u,v}}); }
			if (m==11) { op({{1,w,u}, {1,w,v}, {2,z,v}, {2,z,u}}); }
			if (m==12) { op({{1,w,u}, {1,w,v}, {2,v,u}, {2,z,u}}); }
			if (m==13) { op({{1,u,w}, {1,v,w}, {2,z,w}, {2,v,z}}); }
			if (m==14) { op({{1,u,w}, {2,v,z}}); }
			if (m==15) { op({{1,u,w}, {2,w,z}}); }
			if (m==28) { op({{1,u,w}, {1,v,w}}); }
			if (m==16) { op({{1,u,w}, {1,v,w}, {2,w,z}, {2,u,v}}); }
			if (m==17) { op({{1,w,u}, {2,z,u}}); }
			if (m==33) { op({{1,u,w}}); }
			if (m==41) { op({{2,u,z}}); }
			if (m==42) { op({{2,u,v}}); }
			if (m==51 && deg[w] != 1) { op({{2,z,w}}); }
        	if (m==56) {
			    uniform_int_distribution<int> dist(0, adj[w].size() - 1);
				int ww1 = adj[w][dist(rng)];
				int ww2 = adj[w][dist(rng)];
				if (ww1 != ww2 && deg[ww1] != 1 && deg[ww2] != 1)op({{1,v,w}, {1,u,w}, {2,w,ww1}, {2,w,ww2}});
			}
		}


		if (gid == 6 && ( (1 <= m && m <= 9) || 29 == m || m == 34 || m == 35 || m == 37 || m == 38 || m == 52)){
			int u = -1, v = -1, w = -1, z = -1;
			vector<int> singles;
			for (int i = 0; i < 4; i++){
				if (sub[i][0] + sub[i][1] + sub[i][2] + sub[i][3] == 1) singles.push_back(i);
			}
			z = nodes[singles[1]];
			u = nodes[singles[0]];
			for (int i = 0; i < 4; i++){
				if (sub[i][singles[1]] == 1) v = nodes[i];
			}
			for (int i = 0; i < 4; i++){
				if (sub[i][singles[0]] == 1) w = nodes[i];
			}
			if (m==1)  { op({{1,w,z}, {2,w,v}}); }
			if (m==2)  { op({{1,u,z}, {2,w,u}}); }
			if (m==3)  { op({{1,u,z}, {1,w,z}, {2,u,w}, {2,v,z}}); }
			if (m==4)  { op({{1,u,v}, {1,w,z}, {2,u,w}, {2,v,z}}); }
			if (m==5)  { op({{1,u,v}, {1,w,z}, {2,u,w}, {2,v,w}}); }
			if (m==6)  { op({{1,u,z}, {2,w,v}}); }
			if (m==7)  { op({{1,w,z}, {1,z,u}, {2,w,v}, {2,w,u}}); }
			if (m==8)  { op({{1,w,z}, {2,v,z}}); }
			if (m==9)  { op({{1,u,v}, {1,u,z}, {2,v,z}, {2,w,v}}); }
			if (m==29) { op({{1,u,v}, {1,u,z}, {1,w,z}}); }
			if (m==34) { op({{1,u,v}, {1,w,z}}); }
			if (m==35) { op({{1,u,z}, {1,w,z}}); }
			if (m==37) { op({{1,u,z}}); }
			if (m==38) { op({{1,u,v}}); }
			if (m==52 && deg[z] != 1) { op({{2,z,v}}); }
		}

		if (gid == 5 && ( (18 <= m && m <= 20) || m == 36 || m == 30 || m == 39 || m == 53 || m == 54)) {
			int u = -1, v = -1, w = -1, z = -1, center = -1;
            vector<int> subsum(4, 0), lone;
			for (int i = 0; i < 4; i++){
                subsum[i] = sub[i][0] + sub[i][1] + sub[i][2] + sub[i][3];
                if (subsum[i] == 3) {center = i;}
				else if (subsum[i] == 1) lone.push_back(i);
            }
			u = nodes[lone[0]];
			v = nodes[lone[1]];
			w = nodes[lone[2]];
			z = nodes[center];
			if (m==18) { op({{1,w,u}, {1,u,v}, {2,z,w}, {2,z,v}}); }
			if (m==19) { op({{1,u,v}, {1,w,v}, {2,w,z}, {2,z,v}}); }
			if (m==20) { op({{1,v,w}, {2,z,w}}); }
			if (m==30) { op({{1,u,v}, {1,w,v}, {1,u,w}}); }
			if (m==36) { op({{1,u,v}, {1,w,v}}); }
			if (m==39) { op({{1,u,v}}); }
			if (m==53 && deg[w] != 1) { op({{2,z,w}}); }
			if (m==54 && deg[w] != 1) { op({{1,u,v}, {2,z,w}}); }
		}
    }

	if (donereps < num_samples || donereps2 < num_samples2) {
		cerr << "Warning: sample-attempt limit reached after " << smpl
		     << " attempts before all requested changes were completed.\n";
	}

	double sumsss = 0;
	for (auto i : cntr) sumsss += i;
	//for (int i = 0; i<6; i++) {cntr[i] /= sumsss; cerr << cntr[i] << " ";} cerr << sumsss << endl;
	//cerr << lklk1 << " " << lklk2 << " " << lklk3 << " " << lklk4 << " " << endl;

	return {donereps, donereps2, smpl};
}

void readGraph(const string& filename,
               vector<vector<int>>& adj,
               vector<vector<char>>& adjMat) {
    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "Error: cannot open file '" << filename << "'\n";
        exit(1);
    }

    vector<pair<int, int>> edgess;
    string line;
    int u, v;
    int mx = -1;

    while (getline(infile, line)) {
        istringstream iss(line);
        if (!(iss >> u >> v)) continue;

        edgess.emplace_back(u, v);
        mx = max(mx, max(u, v));
    }

    int n = mx + 1;

    adj.assign(n, {});
    adjMat.assign(n, vector<char>(n, 0));

    for (auto [u, v] : edgess) {
        adj[u].push_back(v);
        adj[v].push_back(u);
        adjMat[u][v] = adjMat[v][u] = 1;
    }
}
// Return every connected component as a list of node IDs.
vector<vector<int>> get_connected_components(
    const vector<vector<int>>& adj) {
    const int n = static_cast<int>(adj.size());
    vector<char> visited(n, 0);
    vector<vector<int>> components;

    for (int start = 0; start < n; ++start) {
        if (visited[start]) continue;

        vector<int> component;
        queue<int> q;
        q.push(start);
        visited[start] = 1;

        while (!q.empty()) {
            const int u = q.front();
            q.pop();
            component.push_back(u);

            for (int v : adj[u]) {
                if (!visited[v]) {
                    visited[v] = 1;
                    q.push(v);
                }
            }
        }

        components.push_back(move(component));
    }

    return components;
}

// This runs only after all requested transformations are finished. If the
// graph has k components, it adds exactly k-1 bridge edges. On each step it
// chooses a random node X in the current largest component and connects X to a
// random node in one randomly chosen remaining component.
void reconnect_graph(vector<vector<int>>& adj,
                     vector<vector<char>>& bvv) {
    vector<vector<int>> components = get_connected_components(adj);
    if (components.size() <= 1) return;

    const size_t initial_component_count = components.size();
    cerr << "Graph is disconnected: " << initial_component_count
         << " connected components.\n";

    random_device rd;
    mt19937 rng(rd());

    while (components.size() > 1) {
        // Put the current largest component at index 0.
        const auto largest_it = max_element(
            components.begin(), components.end(),
            [](const vector<int>& a, const vector<int>& b) {
                return a.size() < b.size();
            });
        iter_swap(components.begin(), largest_it);

        uniform_int_distribution<size_t> other_component_dist(
            1, components.size() - 1);
        const size_t other_index = other_component_dist(rng);

        uniform_int_distribution<size_t> x_dist(
            0, components[0].size() - 1);
        uniform_int_distribution<size_t> y_dist(
            0, components[other_index].size() - 1);

        const int x = components[0][x_dist(rng)];
        const int y = components[other_index][y_dist(rng)];

        // Nodes in different components cannot already share an edge.
        adj[x].push_back(y);
        adj[y].push_back(x);
        bvv[x][y] = bvv[y][x] = 1;

        // The bridge merges these two components. The merged component remains
        // at index 0 and is reconsidered as the largest on the next iteration.
        components[0].insert(
            components[0].end(),
            components[other_index].begin(),
            components[other_index].end());

        if (other_index + 1 != components.size()) {
            swap(components[other_index], components.back());
        }
        components.pop_back();
    }

    cerr << "Added " << (initial_component_count - 1)
         << " bridge edges to reconnect the graph.\n";
}

ld rmsewc(vector<vector<int>> adj, vector<ld> x, ld evl) {
	int n = adj.size();
	vector<ld> y(n);
	ld rmse = 0;
	for (int i = 0; i < n; i++) {
	    for (int nei : adj[i]) {
	        y[i] += x[nei];
	    }
		y[i] /= evl;
		ld d = log(y[i]) - log(x[i]);
		rmse += d*d;
	}
	return sqrt(rmse/n);
}



pair<ld, vector<ld>> evlevc(vector<vector<int>> adj) {

	int n = adj.size();
    vector<ld> x(n);

    mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
    uniform_real_distribution<ld> dist(0.0, 1.0);

    for (int i = 0; i < n; i++)
        x[i] = dist(rng);

    ld norm = 0;
    for (ld val : x)
        norm += val * val;

    norm = sqrt(norm);

    for (ld &val : x)
        val /= norm;

    ld prev_lambda = 0;

    for (int iter = 1; iter <= MAX_ITERS; iter++) {
        vector<ld> y(n, 0.0);

        for (int i = 0; i < n; i++) {
            for (int nei : adj[i]) {
                y[i] += x[nei];
            }
        }

        ld lambda = 0;
        for (int i = 0; i < n; i++) {
            lambda += x[i] * y[i];
        }

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

    // stdout is reserved exclusively for the generated edge list.
    // Diagnostics must go to stderr or they become invalid lines in the graph.
    cerr << fixed << setprecision(12);
    cerr << "Largest eigenvalue (approx): " << prev_lambda << '\n';
    cerr << "Reached max iterations.\n";

    return {prev_lambda, x};
}

vector<ld> matchsortld(vector<ld> a, vector<ld> b) {
    int n = b.size();

    sort(a.begin(), a.end());

    vector<int> idx(n);
    iota(idx.begin(), idx.end(), 0);

    sort(idx.begin(), idx.end(), [&](int i, int j) {
        return b[i] < b[j];
    });

    vector<ld> res(n);
    for (int k = 0; k < n; k++) {
        res[idx[k]] = a[k];
    }

    return res;
}


vector<int> matchsortint(vector<int> a, vector<int> b) {
    int n = b.size();

    sort(a.begin(), a.end());

    vector<int> idx(n);
    iota(idx.begin(), idx.end(), 0);

    sort(idx.begin(), idx.end(), [&](int i, int j) {
        return b[i] < b[j];
    });

    vector<int> res(n);
    for (int k = 0; k < n; k++) {
        res[idx[k]] = a[k];
    }

    return res;
}

vector<int> readGoodNodes(int k) {
    ifstream infile("data_middle1.txt");
    if (!infile.is_open()) {
        cerr << "Error: cannot open file data_middle.txt \n";
        exit(1);
    }

    vector<int> goodnodes;

    for (int i = 0; ; i++) {
        int a[15];

        bool ok = true;
        for (int j = 0; j < 15; j++) {
            if (!(infile >> a[j])) {
                ok = false;
                break;
            }
        }

        if (!ok) break;

        if (a[k] > 0) {
            goodnodes.push_back(i);
        }
    }

    return goodnodes;
}

int main(int argc, char* argv[]) {
    if (argc < 10) {
        cerr << "Usage: " << argv[0]
             << " <target_file> <synth_file> <num_samples> <num_samples2>"
             << " <mode> <mode2> <evbound> <degbound> <hdgbound>"
             << " [completed_samples_file]\n";
        return 1;
    }
	vector<vector<char>> bvv, bvvt;
    vector<vector<int>> adj, adjt;
    readGraph(argv[1], adjt, bvvt);
    readGraph(argv[2], adj, bvv);
	int num_samples = stoi(argv[3]);
	int num_samples2 = stoi(argv[4]);
	int m = stoi(argv[5]);
	int m2 = stoi(argv[6]);
	ld evbound = stod(argv[7]);
	int degbound = stoi(argv[8]);
	double degboundf = (double)degbound/140;
	int hdgbound = stoi(argv[9]);
	const string completed_samples_file =
	    argc >= 11 ? argv[10] : "expfy_samples_done.tmp";

	auto [evlt, evct] = evlevc(adjt);
	auto [evl, evc] = evlevc(adj);

	int n = adj.size();
	vector<int> deg(n), degt(n);
	for(int i=0; i<n; i++){
		deg[i] = adj[i].size();
		degt[i] = adjt[i].size();
	}

	vector<int> gdnd = readGoodNodes(14);
	//cerr << gdnd.size() << endl;
	vector<ld> evcb = matchsortld(evct, evc);
	vector<int> degb = matchsortint(degt, deg);

    const SampleStats stats = sample_edge_prob(
        adj, bvv, num_samples, num_samples2, m, m2, evcb, evlt,
        evbound, degb, degboundf, gdnd, hdgbound);

    ofstream completed_out(completed_samples_file);
    if (!completed_out) {
        cerr << "Error: cannot write completed sample counts to '"
             << completed_samples_file << "'.\n";
        return 1;
    }
    completed_out << stats.done1 << " " << stats.done2 << "\n";
    completed_out.close();

    // Connectivity repair is deliberately done after transformation sampling,
    // so these bridge edges do not count as completed graphlet transformations.
    reconnect_graph(adj, bvv);

    for (size_t u = 0; u < adj.size(); u++) {
        for (int v : adj[u]) {
            if (u < static_cast<size_t>(v)) cout << u << " " << v << "\n";
        }
    }

    return 0;
}
