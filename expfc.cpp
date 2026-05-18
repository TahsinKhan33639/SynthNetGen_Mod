#include <chrono>
#include <bits/stdc++.h>
using namespace std;
using ld = double;
const int MAX_ITERS = 10000;
const ld EPS = 1e-12;

using Clock = std::chrono::steady_clock;

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

// adj -> gid
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

// Main loop
void sample_edge_prob(vector<vector<int>>& adj, vector<vector<char>>& bvv, int num_samples, int num_samples2, int m1, int m2,
					 vector<ld> x, ld evl, ld evbound, vector<int> degb, int degbound) {

    int n = adj.size();
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
	int leaft, leafc, leafz;
	for(int i=0; i<n; i++){
		if (degb[i] == 1) leaft++;
		if (deg[i] == 1) leafc++;
	}

	auto j = [&](int u, int v) {
	    adj[v].push_back(u);
	    adj[u].push_back(v);
	    bvv[v][u] = bvv[u][v] = 1;
		y[u] = y[u] + 1/evl*x[v];
		y[v] = y[v] + 1/evl*x[u];
		deg[u]++; deg[v]++;
		if(deg[u] == 2) leafc--;
		if(deg[v] == 2) leafc--;
		edgec++;
	};

	auto d = [&](int u, int v) {
	    adj[v].erase(remove(adj[v].begin(), adj[v].end(), u), adj[v].end());
	    adj[u].erase(remove(adj[u].begin(), adj[u].end(), v), adj[u].end());
	    bvv[v][u] = bvv[u][v] = 0;
		y[u] = y[u] - 1/evl*x[v];
		y[v] = y[v] - 1/evl*x[u];
		deg[u]--; deg[v]--;
		if(deg[u] == 1) leafc++;
		if(deg[v] == 1) leafc++;
		edgec--;
	};

	int donereps = 0;
	int donereps2 = 0;
	int m = m1;

    auto op = [&](const vector<tuple<int,int,int>>& ops) {
        bool open = true;
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

            return false;
        };

        auto valideg = [&](int node, int new_val) {
            int old_val = deg[node];
            if (new_val <= 0 || old_val <= 0 || x[node] <= 0)
                return false;

            int target = degb[node];

            int new_dist = abs(new_val - target);
            if (new_dist <= degbound)
                return true;

            int old_dist = abs(old_val - target);
            if (new_dist < old_dist)
                return true;

            return false;
        };

        set<int> nodesinop;
        for (auto [type, n1, n2] : ops) {
            nodesinop.insert(n1);
            nodesinop.insert(n2);
        }

        for (auto ndX : nodesinop){
            z[ndX] = y[ndX];
            degz[ndX] = deg[ndX];
        }

        for (auto [type, n1, n2] : ops) {
            if (type == 1){
                z[n1] += 1/evl*x[n2];
                z[n2] += 1/evl*x[n1];
                degz[n1]++; degz[n2]++;
            }
            if (type == 2){
                z[n1] -= 1/evl*x[n2];
                z[n2] -= 1/evl*x[n1];
                degz[n1]--; degz[n2]--;
            }
        }

        for (auto ndX : nodesinop){
            open &= valid(ndX, z[ndX]);
            open &= valideg(ndX, degz[ndX]);
        }

        leafz = leafc;
        for(auto ndX : nodesinop){
            if (deg[ndX] == 1 && degz[ndX] != 1) leafz--;
            if (deg[ndX] != 1 && degz[ndX] == 1) leafz++;
        }
        bool subopen = false;
        if (abs(leafz - leaft) <= abs(leafc - leaft)) subopen = true;
        if (leafz >= leaft*(1-0.02*degbound) && leafz <= leaft*(1+0.02*degbound)) subopen = true;
        open &= subopen;

        if (!open) return;
        if (m == m1) donereps++;
        if (m == m2) donereps2++;
        for (auto [type, n1, n2] : ops) {
            if (type == 1) j(n1, n2);
            if (type == 2) d(n1, n2);
        }
    };


	auto start = Clock::now();

	num_samples = max(1, num_samples);
	num_samples2 = max(1, num_samples2);

	for (int smpl = 0; (donereps < num_samples) || (donereps2 < num_samples2); smpl++) {

		if ((double)donereps / num_samples > (double)donereps2 / num_samples2 + 0.05) m = m2;
		if ((double)donereps / num_samples < (double)donereps2 / num_samples2 - 0.05) m = m1;

		if (donereps2 == num_samples2) m = m1;
		if (donereps == num_samples) m = m2;

		//if (smpl > 1000*num_samples) {cerr << "Too rare initial\n"; break;}
		auto now = Clock::now();
		if (std::chrono::duration_cast<std::chrono::seconds>(now - start).count() >= 20) {
        	break;
		}
        vector<int> nodes;

        uniform_int_distribution<int> dist0(0, n - 1);
        int node1 = dist0(rng);
		nodes.push_back(node1);

		vector<char> blocked(adj.size(), false);
		for (int x : nodes) blocked[x] = true;
		blocked[node1] = true;
	    vector<int> options;
	    options.reserve(32);
		auto pick_neighbor = [&](const vector<int>& candidates) -> int {
		    for (int u : candidates) {
		        for (int v : adj[u]) {
		            if (!blocked[v]) {
		                blocked[v] = true;
		                options.push_back(v);
		            }
		        }
		    }
		    if (options.empty()) return -1;

		    uniform_int_distribution<int> dist(0, options.size() - 1);
		    int newx = options[dist(rng)];
			options.erase(remove(options.begin(), options.end(), newx),options.end());
			return newx;
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

        // Build adj
        vector<vector<int>> sub(4, vector<int>(4, 0));
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < a; b++) {
                int u = nodes[a], v = nodes[b];
                sub[a][b] = sub[b][a] = bvv[u][v];
            }
        }
		// if (sub[3][1] != sub[1][3] || sub[1][2] != sub[2][1] || sub[2][3] != sub[3][2] || sub[0][1] != sub[1][0] || sub[0][2] != sub[2][0] || sub[0][3] != sub[3][0]) cerr << "Something wrong happaned but idk what\n";
        int gid = get_graph_id(sub);

		if (gid == 10 && (43 <= m && m <= 47)){
			int u = -1, v = -1, w = -1, z = -1;
			u = nodes[0];
			v = nodes[1];
			w = nodes[2];
			z = nodes[3];
			// 43: disconnect u-v
            // 44: disconnect z-w, u-v
            // 45: disconnect w-v, w-u
            // 46: to 6
            // 47: to 5
			if (m==43) { op({{2,u,v}}); }
			if (m==44) { op({{2,u,v}, {2,w,z}}); }
			if (m==45) { op({{2,u,z}, {2,v,z}}); }
			if (m==46) { op({{2,u,v}, {2,w,z}, {2,u,z}}); }
			if (m==47) { op({{2,u,v}, {2,u,w}, {2,v,w}}); }
		}
		if (gid == 9 && (m == 24 || m == 31 || m == 55 || (48 <= m && m <= 49))){
			int u = -1, v = -1, w = -1, z = -1;
			bool found = false;
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
			// 24: Disconnected u z
			// 31: to 10
			// 48: to 8
			// 49: to 7
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
        if (gid == 8 && ( (21 <= m && m <= 23) || m == 27 || m == 40 || m == 32 || m == 50)){
			int u = -1, v = -1, w = -1, z = -1;
            vector<int> dig;
            w = nodes[0];
            for (int i = 1; i < 4; i++){
                if (sub[0][i] == 1) dig.push_back(i);
                else z = nodes[i];
            }
            v = nodes[dig[0]];
            u = nodes[dig[1]];
			// diagonal
			// 21: u z
			// center leaf
			// 22: u v
			// 23: u w
			// 27: to 10
			// 32: add uv
			// 40: to 6 (remove uz)
			if (m==21) { op({{1,u,v}, {1,w,z}, {2,u,z}, {2,w,v}}); }
			if (m==22) { op({{1,v,u}, {1,z,w}, {2,z,v}, {2,w,v}}); }
			if (m==23) { op({{1,u,v}, {2,v,w}}); }
			if (m==32) { op({{1,u,v}}); }
			if (m==27) { op({{1,u,v}, {1,z,w}}); }
			if (m==40) { op({{2,u,z}}); }
			if (m==50) { op({{1,w,z}, {2,v,w}, {2,u,w}}); }
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
			// center, leaf
			// 10: z u
			// 11: w z
			// 12: w u
			// 13: u z
			// 14: u v
			// 15: u w
			// to 8 Diagonal
			// 16: u v
			// 17: u z
			// 28: to 10
			// 33: to 9 (add wv)
			// 41: to 6 (remove zu)
			// 42: to 5 (remove uv)
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
			// u--w--v--z
			// 1: u--w--z--v
			// 2: u--z--v--w
			// 3: u--z--w--v
			// 4: u--v--w--z
			// 5: u--v--z--w
			// 6: w--u--z--v
			// 7: w--z--u--v
			// 8: w center
			// 9: u center
			// 29: to 10
			// 34: to 9 (inverted u-z)
			// 35: to 9 (inverted u-v)
			// 37: to 8 (add uz)
			// 38: to 7 (add uv)
			// 52: disconnect z
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
			//   z with u,v,w
			// 18: center u
			// 19: z-u-v-w
			// 20: u-z-v-w
			// 30: to 10
			// 36: to 9 (inverted u-w)
			// 39: to 7 (add wv)
			// 53: disconnect w
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
	// cerr << donereps << " " << donereps2 << endl;
}

// Read graph
void readGraph(const string& filename, vector<vector<int>>& adj, vector<vector<char>>& adjMat) {
    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "cannot open file '" << filename << "'\n";
        exit(1);
    }

    unordered_map<string, int> nodeToId;
    vector<pair<int, int>> edges;
    string line, ustr, vstr;
    int nextId = 0;

    while (getline(infile, line)) {
        istringstream iss(line);
        if (!(iss >> ustr >> vstr)) continue;

        if (!nodeToId.count(ustr)) nodeToId[ustr] = nextId++;
        if (!nodeToId.count(vstr)) nodeToId[vstr] = nextId++;

        int u = nodeToId[ustr];
        int v = nodeToId[vstr];
        edges.emplace_back(u, v);
    }

    adj.assign(nextId, {});
    adjMat.assign(nextId, vector<char>(nextId, 0));

    for (auto [u, v] : edges) {
        adj[u].push_back(v);
        adj[v].push_back(u);
        adjMat[u][v] = adjMat[v][u] = 1;
    }
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
            cerr << "Zero vector\n";
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

vector<ld> matchsortld(vector<ld> a, vector<ld> b) {
    int n = a.size();

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
    int n = a.size();

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


int main(int argc, char* argv[]) {
    if (argc < 8) {
        cerr << "Usage: " << argv[0] << " <target_file> <synth_file> <num_samples> <num_samples2> <mode> <mode2> <evbound> <degbound>\n";
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

	auto [evlt, evct] = evlevc(adjt);
	auto [evl, evc] = evlevc(adj);

	int n = adj.size();
	vector<int> deg(n), degt(n);
	for(int i=0; i<n; i++){
		deg[i] = adj[i].size();
		degt[i] = adjt[i].size();
	}

	vector<ld> evcb = matchsortld(evct, evc);
	vector<int> degb = matchsortint(degt, deg);

    sample_edge_prob(adj, bvv, num_samples, num_samples2, m, m2, evcb, evlt, evbound, degb, degbound);

    for(int u = 0; u < adj.size(); u++){
        for(int v : adj[u]){
            if(u < v) cout << u << " " << v << "\n";
        }
    }

    return 0;
}
