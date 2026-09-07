#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace std;

#ifndef EXPFY_MIN_SAMPLE_ATTEMPTS
#define EXPFY_MIN_SAMPLE_ATTEMPTS 8000000LL
#endif
#ifndef EXPFY_MAX_SAMPLE_ATTEMPTS
#define EXPFY_MAX_SAMPLE_ATTEMPTS 10000000LL
#endif
#ifndef EXPFY_ATTEMPTS_PER_REQUESTED_CHANGE
#define EXPFY_ATTEMPTS_PER_REQUESTED_CHANGE 500LL
#endif
#ifndef EXPFY_MAX_NEIGHBOR_PICK_ATTEMPTS
#define EXPFY_MAX_NEIGHBOR_PICK_ATTEMPTS 100
#endif

namespace {

constexpr long long MIN_SAMPLE_ATTEMPTS = EXPFY_MIN_SAMPLE_ATTEMPTS;
constexpr long long MAX_SAMPLE_ATTEMPTS = EXPFY_MAX_SAMPLE_ATTEMPTS;
constexpr long long ATTEMPTS_PER_REQUESTED_CHANGE =
    EXPFY_ATTEMPTS_PER_REQUESTED_CHANGE;
constexpr int MAX_NEIGHBOR_PICK_ATTEMPTS =
    EXPFY_MAX_NEIGHBOR_PICK_ATTEMPTS;
constexpr int ORCA4_COLUMNS = 15;
constexpr int GOOD_NODE_COLUMN = 14;
constexpr int MAX_OPS = 4;
constexpr int MAX_TOUCHED_NODES = 2 * MAX_OPS;
constexpr int MAX_RELATED_DEGREES = 2 * MAX_TOUCHED_NODES;

struct SampleStats {
    int done1 = 0;
    int done2 = 0;
    long long attempts = 0;
};

// xoshiro256**. It is substantially cheaper than repeatedly constructing and
// invoking standard-library distribution objects inside the sampling loop.
class FastRng {
public:
    explicit FastRng(uint64_t seed) {
        for (uint64_t& value : state_) {
            value = splitmix64(seed);
        }
        if ((state_[0] | state_[1] | state_[2] | state_[3]) == 0) {
            state_[0] = 1;
        }
    }

    uint64_t next_u64() {
        const uint64_t result = rotl(state_[1] * 5, 7) * 9;
        const uint64_t t = state_[1] << 17;

        state_[2] ^= state_[0];
        state_[3] ^= state_[1];
        state_[1] ^= state_[2];
        state_[0] ^= state_[3];
        state_[2] ^= t;
        state_[3] = rotl(state_[3], 45);
        return result;
    }

    uint64_t bounded(uint64_t bound) {
        if (bound <= 1) return 0;

#if defined(__SIZEOF_INT128__)
        // Lemire's unbiased multiply-high reduction.
        uint64_t x = next_u64();
        __uint128_t product = static_cast<__uint128_t>(x) * bound;
        uint64_t low = static_cast<uint64_t>(product);
        if (low < bound) {
            const uint64_t threshold = static_cast<uint64_t>(-bound) % bound;
            while (low < threshold) {
                x = next_u64();
                product = static_cast<__uint128_t>(x) * bound;
                low = static_cast<uint64_t>(product);
            }
        }
        return static_cast<uint64_t>(product >> 64);
#else
        const uint64_t threshold = static_cast<uint64_t>(-bound) % bound;
        uint64_t value;
        do {
            value = next_u64();
        } while (value < threshold);
        return value % bound;
#endif
    }

    bool chance(unsigned numerator, unsigned denominator) {
        return bounded(denominator) < numerator;
    }

private:
    array<uint64_t, 4> state_{};

    static uint64_t rotl(uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }

    static uint64_t splitmix64(uint64_t& x) {
        uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }
};

uint64_t make_seed() {
    if (const char* text = getenv("EXPFY_SEED")) {
        try {
            size_t used = 0;
            const string value(text);
            const uint64_t seed = stoull(value, &used, 0);
            if (used == value.size()) return seed;
        } catch (const exception&) {
            // Fall through to a nondeterministic seed.
        }
        cerr << "Warning: ignoring invalid EXPFY_SEED='" << text << "'.\n";
    }

    random_device rd;
    uint64_t seed =
        static_cast<uint64_t>(chrono::high_resolution_clock::now()
                                  .time_since_epoch()
                                  .count());
    seed ^= static_cast<uint64_t>(rd()) << 32;
    seed ^= static_cast<uint64_t>(rd());
    return seed;
}

// One packed bit row per vertex. This uses about n^2/8 bytes rather than the
// old vector<vector<char>> matrix's n^2 bytes plus n row allocations.
class BitAdjacency {
public:
    BitAdjacency() = default;
    explicit BitAdjacency(int n) { reset(n); }

    void reset(int n) {
        n_ = n;
        words_per_row_ = (static_cast<size_t>(n) + 63) / 64;
        bits_.assign(static_cast<size_t>(n) * words_per_row_, 0);
    }

    [[nodiscard]] bool get(int u, int v) const {
        const size_t index =
            static_cast<size_t>(u) * words_per_row_ +
            static_cast<size_t>(v >> 6);
        return (bits_[index] >> (v & 63)) & 1ULL;
    }

    void set(int u, int v, bool present) {
        const size_t index =
            static_cast<size_t>(u) * words_per_row_ +
            static_cast<size_t>(v >> 6);
        const uint64_t mask = 1ULL << (v & 63);
        if (present) {
            bits_[index] |= mask;
        } else {
            bits_[index] &= ~mask;
        }
    }

    void set_undirected(int u, int v, bool present) {
        set(u, v, present);
        set(v, u, present);
    }

private:
    int n_ = 0;
    size_t words_per_row_ = 0;
    vector<uint64_t> bits_;
};

struct EdgeListData {
    int n = 0;
    vector<pair<int, int>> edges;
    vector<int> degree;
};

EdgeListData read_edge_list(const string& filename) {
    ifstream input(filename);
    if (!input) {
        throw runtime_error("cannot open graph file '" + filename + "'");
    }

    vector<pair<int, int>> edges;
    string line;
    int maximum_vertex = -1;

    while (getline(input, line)) {
        istringstream parser(line);
        int u = -1;
        int v = -1;
        if (!(parser >> u >> v)) continue;
        if (u < 0 || v < 0) {
            throw runtime_error("negative vertex ID in '" + filename + "'");
        }
        if (u == v) continue;
        if (u > v) swap(u, v);
        edges.emplace_back(u, v);
        maximum_vertex = max(maximum_vertex, v);
    }

    if (maximum_vertex < 0) {
        throw runtime_error("graph file '" + filename + "' has no edges");
    }

    sort(edges.begin(), edges.end());
    edges.erase(unique(edges.begin(), edges.end()), edges.end());

    EdgeListData result;
    result.n = maximum_vertex + 1;
    result.edges = move(edges);
    result.degree.assign(result.n, 0);
    for (const auto& [u, v] : result.edges) {
        result.degree[u]++;
        result.degree[v]++;
    }
    return result;
}

struct Graph {
    int n = 0;
    long long edge_count = 0;
    vector<vector<int>> adj;
    vector<int> degree;
    BitAdjacency matrix;

    explicit Graph(const EdgeListData& data)
        : n(data.n),
          edge_count(static_cast<long long>(data.edges.size())),
          adj(data.n),
          degree(data.degree),
          matrix(data.n) {
        for (int u = 0; u < n; ++u) {
            adj[u].reserve(static_cast<size_t>(degree[u]) + 4);
        }
        for (const auto& [u, v] : data.edges) {
            adj[u].push_back(v);
            adj[v].push_back(u);
            matrix.set_undirected(u, v, true);
        }
    }

    static bool erase_neighbor_swap(vector<int>& neighbors, int value) {
        const auto it = find(neighbors.begin(), neighbors.end(), value);
        if (it == neighbors.end()) return false;
        *it = neighbors.back();
        neighbors.pop_back();
        return true;
    }

    void add_edge_unchecked(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
        matrix.set_undirected(u, v, true);
        degree[u]++;
        degree[v]++;
        edge_count++;
    }

    void remove_edge_unchecked(int u, int v) {
        const bool removed_u = erase_neighbor_swap(adj[u], v);
        const bool removed_v = erase_neighbor_swap(adj[v], u);
        if (!removed_u || !removed_v) {
            throw runtime_error("adjacency list/matrix inconsistency while deleting edge");
        }
        matrix.set_undirected(u, v, false);
        degree[u]--;
        degree[v]--;
        edge_count--;
    }
};

constexpr array<uint8_t, 64> GRAPH_ID_BY_BITS = {
    0, 1, 1, 3, 1, 3, 3, 4,
    1, 3, 3, 5, 2, 6, 6, 7,
    1, 3, 2, 6, 3, 5, 6, 7,
    3, 4, 6, 7, 6, 7, 8, 9,
    1, 2, 3, 6, 3, 6, 5, 7,
    3, 6, 4, 7, 6, 8, 7, 9,
    3, 6, 6, 8, 4, 7, 7, 9,
    5, 7, 7, 9, 7, 9, 9, 10,
};

constexpr array<int8_t, 61> make_required_gid_table() {
    array<int8_t, 61> table{};
    for (int i = 0; i <= 60; ++i) table[i] = -1;

    for (int i = 1; i <= 9; ++i) table[i] = 6;
    table[29] = table[34] = table[35] = table[37] = table[38] = table[52] = 6;

    for (int i = 10; i <= 17; ++i) table[i] = 7;
    table[28] = table[33] = table[41] = table[42] = table[51] = table[56] = 7;

    table[18] = table[19] = table[20] = 5;
    table[30] = table[36] = table[39] = table[53] = table[54] = 5;

    table[21] = table[22] = table[23] = 8;
    table[27] = table[32] = table[40] = table[50] = table[57] = table[58] = 8;

    table[24] = table[31] = table[48] = table[49] = table[55] = 9;

    for (int i = 43; i <= 47; ++i) table[i] = 10;
    table[59] = 10;

    // 25 and 26 are unused. Mode 60 has its own five-node test.
    return table;
}

constexpr auto REQUIRED_GID = make_required_gid_table();

using Sub4 = array<array<uint8_t, 4>, 4>;

int graph_id_from_sub4(const Sub4& sub) {
    unsigned bits = 0;
    bits |= static_cast<unsigned>(sub[2][3]) << 0;
    bits |= static_cast<unsigned>(sub[1][3]) << 1;
    bits |= static_cast<unsigned>(sub[1][2]) << 2;
    bits |= static_cast<unsigned>(sub[0][3]) << 3;
    bits |= static_cast<unsigned>(sub[0][2]) << 4;
    bits |= static_cast<unsigned>(sub[0][1]) << 5;
    return GRAPH_ID_BY_BITS[bits];
}

vector<int> match_target_degrees_by_current_rank(
    vector<int> target_degrees,
    const vector<int>& current_degrees) {
    if (target_degrees.size() != current_degrees.size()) {
        throw runtime_error("target and synthetic graphs have different node counts");
    }

    sort(target_degrees.begin(), target_degrees.end());
    vector<int> indices(current_degrees.size());
    iota(indices.begin(), indices.end(), 0);
    sort(indices.begin(), indices.end(), [&](int lhs, int rhs) {
        return current_degrees[lhs] < current_degrees[rhs];
    });

    vector<int> matched(current_degrees.size());
    for (size_t rank = 0; rank < indices.size(); ++rank) {
        matched[indices[rank]] = target_degrees[rank];
    }
    return matched;
}

vector<uint8_t> read_good_node_mask(int n) {
    ifstream input("data_middle1.txt");
    if (!input) {
        throw runtime_error("cannot open data_middle1.txt");
    }

    vector<uint8_t> mask(n, 0);
    for (int row = 0; row < n; ++row) {
        for (int column = 0; column < ORCA4_COLUMNS; ++column) {
            long long value = 0;
            if (!(input >> value)) {
                throw runtime_error(
                    "data_middle1.txt ended before " + to_string(n) +
                    " complete 15-column rows were read");
            }
            if (column == GOOD_NODE_COLUMN && value > 0) mask[row] = 1;
        }
    }

    long long extra = 0;
    if (input >> extra) {
        throw runtime_error(
            "data_middle1.txt contains more rows than the synthetic graph");
    }
    return mask;
}

enum class EdgeAction : uint8_t { Add = 1, Remove = 2 };

struct EdgeOp {
    EdgeAction action;
    int u;
    int v;
};

class Sampler {
public:
    Sampler(Graph& graph,
            vector<int> target_degree,
            double degree_log_bound,
            int histogram_bound,
            vector<uint8_t> good_mask,
            uint64_t seed)
        : graph_(graph),
          target_degree_(move(target_degree)),
          degree_log_bound_(degree_log_bound),
          histogram_bound_(histogram_bound),
          good_mask_(move(good_mask)),
          rng_(seed),
          current_histogram_(graph.n, 0),
          target_histogram_(graph.n, 0),
          target_log_degree_(graph.n, 0.0),
          log_degree_(graph.n, -numeric_limits<double>::infinity()) {
        if (graph_.n < 4) {
            throw runtime_error("graph must have at least 4 nodes");
        }
        if (target_degree_.size() != static_cast<size_t>(graph_.n) ||
            good_mask_.size() != static_cast<size_t>(graph_.n)) {
            throw runtime_error("internal sampler vector-size mismatch");
        }
        if (!isfinite(degree_log_bound_) || degree_log_bound_ < 0.0) {
            throw runtime_error("degbound must encode a finite nonnegative log bound");
        }
        if (histogram_bound_ < 0) {
            throw runtime_error("hdgbound must be nonnegative");
        }

        for (int degree = 1; degree < graph_.n; ++degree) {
            log_degree_[degree] = log(static_cast<double>(degree));
        }

        for (int node = 0; node < graph_.n; ++node) {
            const int current = graph_.degree[node];
            const int target = target_degree_[node];
            if (current < 0 || current >= graph_.n ||
                target <= 0 || target >= graph_.n) {
                throw runtime_error(
                    "target/current degrees must lie in [1,n-1]");
            }
            current_histogram_[current]++;
            target_histogram_[target]++;
            target_log_degree_[node] = log_degree_[target];

            if (good_mask_[node]) {
                good_nodes_.push_back(node);
            } else {
                ordinary_nodes_.push_back(node);
            }
        }

        // Preserve the old fallback:
        // abs(new_degree-target_degree) <= max(degbound*n*3/100, 5).
        degree_absolute_fallback_ =
            max(degree_log_bound_ * graph_.n * 3.0 / 100.0, 5.0);
    }

    SampleStats run(int requested1, int requested2, int mode1, int mode2) {
        requested1 = max(0, requested1);
        requested2 = max(0, requested2);
        validate_mode(mode1);
        validate_mode(mode2);

        requested1_ = requested1;
        requested2_ = requested2;
        mode1_ = mode1;
        mode2_ = mode2;

        const long long requested_total =
            static_cast<long long>(requested1_) + requested2_;
        long long attempt_limit = 0;
        if (requested_total > 0) {
            const long long scaled =
                requested_total >
                        MAX_SAMPLE_ATTEMPTS /
                            ATTEMPTS_PER_REQUESTED_CHANGE
                    ? MAX_SAMPLE_ATTEMPTS
                    : requested_total * ATTEMPTS_PER_REQUESTED_CHANGE;
            attempt_limit = clamp(
                scaled, MIN_SAMPLE_ATTEMPTS, MAX_SAMPLE_ATTEMPTS);
        }

        long long attempts = 0;
        for (; attempts < attempt_limit &&
               (done1_ < requested1_ || done2_ < requested2_);
             ++attempts) {
            choose_active_slot();
            const int mode = active_slot_ == 1 ? mode1_ : mode2_;
            sample_one(mode);
        }

        if (done1_ < requested1_ || done2_ < requested2_) {
            cerr << "Warning: sample-attempt limit reached after " << attempts
                 << " attempts before all requested changes were completed.\n";
        }

        return {done1_, done2_, attempts};
    }

private:
    Graph& graph_;
    vector<int> target_degree_;
    double degree_log_bound_ = 0.0;
    double degree_absolute_fallback_ = 5.0;
    int histogram_bound_ = 0;
    vector<uint8_t> good_mask_;
    vector<int> good_nodes_;
    vector<int> ordinary_nodes_;
    FastRng rng_;
    vector<int> current_histogram_;
    vector<int> target_histogram_;
    vector<double> target_log_degree_;
    vector<double> log_degree_;

    int requested1_ = 0;
    int requested2_ = 0;
    int done1_ = 0;
    int done2_ = 0;
    int mode1_ = 1;
    int mode2_ = 1;
    int active_slot_ = 1;

    static void validate_mode(int mode) {
        if (mode < 1 || mode > 60 || mode == 25 || mode == 26) {
            throw runtime_error("transformation mode must be 1..60 except 25 and 26");
        }
    }

    void choose_active_slot() {
        const bool need1 = done1_ < requested1_;
        const bool need2 = done2_ < requested2_;

        if (!need1) {
            active_slot_ = 2;
        } else if (!need2) {
            active_slot_ = 1;
        } else {
            const double progress1 =
                static_cast<double>(done1_) / requested1_;
            const double progress2 =
                static_cast<double>(done2_) / requested2_;
            if (progress1 > progress2 + 0.02) active_slot_ = 2;
            if (progress1 < progress2 - 0.02) active_slot_ = 1;
        }
    }

    int pick_first_node() {
        // The old rejection loop gave each good node weight 1 and each other
        // node weight 0.2. Multiplying by five yields the exact integer weights
        // 5:1, avoiding a geometric retry loop without changing the distribution.
        const uint64_t good_weight =
            static_cast<uint64_t>(good_nodes_.size()) * 5ULL;
        const uint64_t total_weight =
            good_weight + static_cast<uint64_t>(ordinary_nodes_.size());
        if (total_weight == 0) return -1;

        const uint64_t choice = rng_.bounded(total_weight);
        if (choice < good_weight) {
            return good_nodes_[choice / 5ULL];
        }
        return ordinary_nodes_[choice - good_weight];
    }

    static bool contains_node(
        const array<int, 5>& nodes, int count, int value) {
        for (int i = 0; i < count; ++i) {
            if (nodes[i] == value) return true;
        }
        return false;
    }

    int pick_neighbor(const array<int, 5>& nodes, int count) {
        for (int attempt = 0;
             attempt < MAX_NEIGHBOR_PICK_ATTEMPTS;
             ++attempt) {
            const int u = nodes[rng_.bounded(static_cast<uint64_t>(count))];
            if (graph_.adj[u].empty()) continue;

            const int v = graph_.adj[u][rng_.bounded(graph_.adj[u].size())];
            if (contains_node(nodes, count, v)) continue;

            // Preserve the old 70% rejection bias for non-good nodes.
            if (!good_mask_[v] && rng_.chance(7, 10)) continue;
            return v;
        }
        return -1;
    }

    void shuffle_first_four(array<int, 5>& nodes) {
        for (int i = 3; i > 0; --i) {
            const int j = static_cast<int>(
                rng_.bounded(static_cast<uint64_t>(i + 1)));
            swap(nodes[i], nodes[j]);
        }
    }

    int random_neighbor(int node) {
        if (graph_.adj[node].empty()) return -1;
        return graph_.adj[node][rng_.bounded(graph_.adj[node].size())];
    }

    Sub4 build_sub4(const array<int, 5>& nodes) const {
        Sub4 sub{};
        for (int a = 0; a < 4; ++a) {
            for (int b = 0; b < a; ++b) {
                const uint8_t edge =
                    graph_.matrix.get(nodes[a], nodes[b]) ? 1 : 0;
                sub[a][b] = edge;
                sub[b][a] = edge;
            }
        }
        return sub;
    }

    static int sub_degree(const Sub4& sub, int node) {
        return sub[node][0] + sub[node][1] +
               sub[node][2] + sub[node][3];
    }

    bool degree_allowed(int node, int old_degree, int new_degree) const {
        if (new_degree <= 0 || new_degree >= graph_.n ||
            old_degree <= 0 || old_degree >= graph_.n) {
            return false;
        }

        const int target = target_degree_[node];
        const double new_log_distance =
            abs(log_degree_[new_degree] - target_log_degree_[node]);
        if (new_log_distance <= degree_log_bound_) return true;

        const int new_distance = abs(new_degree - target);
        if (new_distance <= degree_absolute_fallback_) return true;

        const int old_distance = abs(old_degree - target);
        return new_distance < old_distance;
    }

    bool histogram_allowed(int degree, int new_count) const {
        if (degree < 0 || degree >= graph_.n || new_count < 0) {
            return false;
        }

        const int old_count = current_histogram_[degree];
        const int target_count = target_histogram_[degree];
        const int new_distance = abs(new_count - target_count);
        if (new_distance <= histogram_bound_) return true;

        const int old_distance = abs(old_count - target_count);
        return new_distance < old_distance;
    }

    bool try_operation(initializer_list<EdgeOp> operations) {
        if (operations.size() == 0 || operations.size() > MAX_OPS) {
            return false;
        }

        array<pair<int, int>, MAX_OPS> touched_edges{};
        int edge_count = 0;
        array<int, MAX_TOUCHED_NODES> nodes{};
        int node_count = 0;

        auto add_unique_node = [&](int node) {
            for (int i = 0; i < node_count; ++i) {
                if (nodes[i] == node) return;
            }
            nodes[node_count++] = node;
        };

        for (const EdgeOp& operation : operations) {
            const int u = operation.u;
            const int v = operation.v;
            if (u < 0 || u >= graph_.n ||
                v < 0 || v >= graph_.n || u == v) {
                return false;
            }

            const pair<int, int> edge = minmax(u, v);
            for (int i = 0; i < edge_count; ++i) {
                if (touched_edges[i] == edge) return false;
            }
            touched_edges[edge_count++] = edge;

            const bool exists = graph_.matrix.get(u, v);
            if ((operation.action == EdgeAction::Add && exists) ||
                (operation.action == EdgeAction::Remove && !exists)) {
                return false;
            }

            add_unique_node(u);
            add_unique_node(v);
        }

        array<int, MAX_TOUCHED_NODES> old_degree{};
        array<int, MAX_TOUCHED_NODES> new_degree{};
        for (int i = 0; i < node_count; ++i) {
            old_degree[i] = graph_.degree[nodes[i]];
            new_degree[i] = old_degree[i];
        }

        auto node_position = [&](int node) {
            for (int i = 0; i < node_count; ++i) {
                if (nodes[i] == node) return i;
            }
            return -1;
        };

        for (const EdgeOp& operation : operations) {
            const int delta =
                operation.action == EdgeAction::Add ? 1 : -1;
            new_degree[node_position(operation.u)] += delta;
            new_degree[node_position(operation.v)] += delta;
        }

        for (int i = 0; i < node_count; ++i) {
            if (!degree_allowed(nodes[i], old_degree[i], new_degree[i])) {
                return false;
            }
        }

        array<int, MAX_RELATED_DEGREES> related_degree{};
        array<int, MAX_RELATED_DEGREES> proposed_histogram{};
        int related_count = 0;

        auto related_position = [&](int degree) {
            for (int i = 0; i < related_count; ++i) {
                if (related_degree[i] == degree) return i;
            }
            related_degree[related_count] = degree;
            proposed_histogram[related_count] =
                current_histogram_[degree];
            return related_count++;
        };

        for (int i = 0; i < node_count; ++i) {
            related_position(old_degree[i]);
            related_position(new_degree[i]);
        }
        for (int i = 0; i < node_count; ++i) {
            proposed_histogram[related_position(old_degree[i])]--;
            proposed_histogram[related_position(new_degree[i])]++;
        }
        for (int i = 0; i < related_count; ++i) {
            if (!histogram_allowed(
                    related_degree[i], proposed_histogram[i])) {
                return false;
            }
        }

        // Validation is complete. Mutating the graph can no longer leave a
        // half-applied operation.
        for (const EdgeOp& operation : operations) {
            if (operation.action == EdgeAction::Add) {
                graph_.add_edge_unchecked(operation.u, operation.v);
            } else {
                graph_.remove_edge_unchecked(operation.u, operation.v);
            }
        }

        for (int i = 0; i < node_count; ++i) {
            current_histogram_[old_degree[i]]--;
            current_histogram_[new_degree[i]]++;
            if (graph_.degree[nodes[i]] != new_degree[i]) {
                throw runtime_error("internal degree update mismatch");
            }
        }

        if (active_slot_ == 1 && done1_ < requested1_) {
            done1_++;
        } else if (active_slot_ == 2 && done2_ < requested2_) {
            done2_++;
        }
        return true;
    }

    void sample_one(int mode) {
        array<int, 5> nodes{};
        int count = 0;

        nodes[count++] = pick_first_node();
        if (nodes[0] < 0) return;

        for (; count < 4; ++count) {
            const int next = pick_neighbor(nodes, count);
            if (next < 0) return;
            nodes[count] = next;
        }
        shuffle_first_four(nodes);

        if (mode == 60) {
            const int fifth = pick_neighbor(nodes, 4);
            if (fifth < 0) return;
            nodes[4] = fifth;

            int induced_edges = 0;
            int missing_a = -1;
            int missing_b = -1;
            for (int a = 0; a < 5; ++a) {
                for (int b = 0; b < a; ++b) {
                    if (graph_.matrix.get(nodes[a], nodes[b])) {
                        induced_edges++;
                    } else {
                        missing_a = a;
                        missing_b = b;
                    }
                }
            }
            if (induced_edges == 9) {
                try_operation({
                    {EdgeAction::Add, nodes[missing_a], nodes[missing_b]},
                });
            }
            return;
        }

        const Sub4 sub = build_sub4(nodes);
        const int gid = graph_id_from_sub4(sub);
        if (REQUIRED_GID[mode] != gid) return;

        switch (gid) {
            case 10:
                apply_from_clique(mode, nodes);
                break;
            case 9:
                apply_from_diamond(mode, nodes, sub);
                break;
            case 8:
                apply_from_cycle(mode, nodes, sub);
                break;
            case 7:
                apply_from_paw(mode, nodes, sub);
                break;
            case 6:
                apply_from_path(mode, nodes, sub);
                break;
            case 5:
                apply_from_claw(mode, nodes, sub);
                break;
            default:
                break;
        }
    }

    void apply_from_clique(int mode, const array<int, 5>& nodes) {
        const int u = nodes[0];
        const int v = nodes[1];
        const int w = nodes[2];
        const int z = nodes[3];

        switch (mode) {
            case 43:
                try_operation({{EdgeAction::Remove, u, v}});
                break;
            case 44:
                try_operation({
                    {EdgeAction::Remove, u, v},
                    {EdgeAction::Remove, w, z},
                });
                break;
            case 45:
                try_operation({
                    {EdgeAction::Remove, u, z},
                    {EdgeAction::Remove, v, z},
                });
                break;
            case 46:
                try_operation({
                    {EdgeAction::Remove, u, v},
                    {EdgeAction::Remove, w, z},
                    {EdgeAction::Remove, u, z},
                });
                break;
            case 47:
                try_operation({
                    {EdgeAction::Remove, u, v},
                    {EdgeAction::Remove, u, w},
                    {EdgeAction::Remove, v, w},
                });
                break;
            case 59: {
                const int outside = random_neighbor(u);
                if (outside >= 0) {
                    try_operation({
                        {EdgeAction::Add, outside, v},
                        {EdgeAction::Add, outside, w},
                        {EdgeAction::Add, outside, z},
                    });
                }
                break;
            }
            default:
                break;
        }
    }

    void apply_from_diamond(
        int mode, const array<int, 5>& nodes, const Sub4& sub) {
        array<int, 2> degree2{};
        array<int, 2> degree3{};
        int count2 = 0;
        int count3 = 0;
        for (int i = 0; i < 4; ++i) {
            const int degree = sub_degree(sub, i);
            if (degree == 2 && count2 < 2) degree2[count2++] = i;
            if (degree == 3 && count3 < 2) degree3[count3++] = i;
        }
        if (count2 != 2 || count3 != 2) return;

        const int u = nodes[degree2[0]];
        const int v = nodes[degree2[1]];
        const int w = nodes[degree3[0]];
        const int z = nodes[degree3[1]];

        switch (mode) {
            case 24:
                try_operation({
                    {EdgeAction::Add, u, v},
                    {EdgeAction::Remove, u, z},
                });
                break;
            case 31:
                try_operation({{EdgeAction::Add, u, v}});
                break;
            case 48:
                try_operation({{EdgeAction::Remove, z, w}});
                break;
            case 49:
                try_operation({{EdgeAction::Remove, u, w}});
                break;
            case 55: {
                const int outside = random_neighbor(u);
                if (outside >= 0 && graph_.degree[outside] != 1) {
                    try_operation({
                        {EdgeAction::Add, u, v},
                        {EdgeAction::Remove, u, outside},
                    });
                }
                break;
            }
            default:
                break;
        }
    }

    void apply_from_cycle(
        int mode, const array<int, 5>& nodes, const Sub4& sub) {
        const int w = nodes[0];
        array<int, 2> adjacent{};
        int adjacent_count = 0;
        int z = -1;
        for (int i = 1; i < 4; ++i) {
            if (sub[0][i]) {
                if (adjacent_count < 2) adjacent[adjacent_count++] = i;
            } else {
                z = nodes[i];
            }
        }
        if (adjacent_count != 2 || z < 0) return;

        const int v = nodes[adjacent[0]];
        const int u = nodes[adjacent[1]];

        switch (mode) {
            case 21:
                try_operation({
                    {EdgeAction::Add, u, v},
                    {EdgeAction::Add, w, z},
                    {EdgeAction::Remove, u, z},
                    {EdgeAction::Remove, w, v},
                });
                break;
            case 22:
                try_operation({
                    {EdgeAction::Add, v, u},
                    {EdgeAction::Add, z, w},
                    {EdgeAction::Remove, z, v},
                    {EdgeAction::Remove, w, v},
                });
                break;
            case 23:
                try_operation({
                    {EdgeAction::Add, u, v},
                    {EdgeAction::Remove, v, w},
                });
                break;
            case 27:
                try_operation({
                    {EdgeAction::Add, u, v},
                    {EdgeAction::Add, z, w},
                });
                break;
            case 32:
                try_operation({{EdgeAction::Add, u, v}});
                break;
            case 40:
                try_operation({{EdgeAction::Remove, u, z}});
                break;
            case 50:
                try_operation({
                    {EdgeAction::Add, w, z},
                    {EdgeAction::Remove, v, w},
                    {EdgeAction::Remove, u, w},
                });
                break;
            case 57: {
                const int outside_w = random_neighbor(w);
                const int outside_z = random_neighbor(z);
                if (outside_w >= 0 && outside_z >= 0 &&
                    outside_w != outside_z &&
                    graph_.degree[outside_w] != 1 &&
                    graph_.degree[outside_z] != 1) {
                    try_operation({
                        {EdgeAction::Add, z, w},
                        {EdgeAction::Add, u, v},
                        {EdgeAction::Remove, w, outside_w},
                        {EdgeAction::Remove, z, outside_z},
                    });
                }
                break;
            }
            case 58:
                try_operation({
                    {EdgeAction::Remove, u, z},
                    {EdgeAction::Remove, u, w},
                });
                break;
            default:
                break;
        }
    }

    void apply_from_paw(
        int mode, const array<int, 5>& nodes, const Sub4& sub) {
        int w = -1;
        int z = -1;
        array<int, 2> middle{};
        int middle_count = 0;
        for (int i = 0; i < 4; ++i) {
            const int degree = sub_degree(sub, i);
            if (degree == 1) {
                w = nodes[i];
            } else if (degree == 3) {
                z = nodes[i];
            } else if (degree == 2 && middle_count < 2) {
                middle[middle_count++] = i;
            }
        }
        if (w < 0 || z < 0 || middle_count != 2) return;

        const int u = nodes[middle[0]];
        const int v = nodes[middle[1]];

        switch (mode) {
            case 10:
                try_operation({
                    {EdgeAction::Add, w, v},
                    {EdgeAction::Remove, u, v},
                });
                break;
            case 11:
                try_operation({
                    {EdgeAction::Add, w, u},
                    {EdgeAction::Add, w, v},
                    {EdgeAction::Remove, z, v},
                    {EdgeAction::Remove, z, u},
                });
                break;
            case 12:
                try_operation({
                    {EdgeAction::Add, w, u},
                    {EdgeAction::Add, w, v},
                    {EdgeAction::Remove, v, u},
                    {EdgeAction::Remove, z, u},
                });
                break;
            case 13:
                try_operation({
                    {EdgeAction::Add, u, w},
                    {EdgeAction::Add, v, w},
                    {EdgeAction::Remove, z, w},
                    {EdgeAction::Remove, v, z},
                });
                break;
            case 14:
                try_operation({
                    {EdgeAction::Add, u, w},
                    {EdgeAction::Remove, v, z},
                });
                break;
            case 15:
                try_operation({
                    {EdgeAction::Add, u, w},
                    {EdgeAction::Remove, w, z},
                });
                break;
            case 16:
                try_operation({
                    {EdgeAction::Add, u, w},
                    {EdgeAction::Add, v, w},
                    {EdgeAction::Remove, w, z},
                    {EdgeAction::Remove, u, v},
                });
                break;
            case 17:
                try_operation({
                    {EdgeAction::Add, w, u},
                    {EdgeAction::Remove, z, u},
                });
                break;
            case 28:
                try_operation({
                    {EdgeAction::Add, u, w},
                    {EdgeAction::Add, v, w},
                });
                break;
            case 33:
                try_operation({{EdgeAction::Add, u, w}});
                break;
            case 41:
                try_operation({{EdgeAction::Remove, u, z}});
                break;
            case 42:
                try_operation({{EdgeAction::Remove, u, v}});
                break;
            case 51:
                if (graph_.degree[w] != 1) {
                    try_operation({{EdgeAction::Remove, z, w}});
                }
                break;
            case 56: {
                const int outside1 = random_neighbor(w);
                const int outside2 = random_neighbor(w);
                if (outside1 >= 0 && outside2 >= 0 &&
                    outside1 != outside2 &&
                    graph_.degree[outside1] != 1 &&
                    graph_.degree[outside2] != 1) {
                    try_operation({
                        {EdgeAction::Add, v, w},
                        {EdgeAction::Add, u, w},
                        {EdgeAction::Remove, w, outside1},
                        {EdgeAction::Remove, w, outside2},
                    });
                }
                break;
            }
            default:
                break;
        }
    }

    void apply_from_path(
        int mode, const array<int, 5>& nodes, const Sub4& sub) {
        array<int, 2> leaves{};
        int leaf_count = 0;
        for (int i = 0; i < 4; ++i) {
            if (sub_degree(sub, i) == 1 && leaf_count < 2) {
                leaves[leaf_count++] = i;
            }
        }
        if (leaf_count != 2) return;

        const int u = nodes[leaves[0]];
        const int z = nodes[leaves[1]];
        int v = -1;
        int w = -1;
        for (int i = 0; i < 4; ++i) {
            if (sub[i][leaves[1]]) v = nodes[i];
            if (sub[i][leaves[0]]) w = nodes[i];
        }
        if (v < 0 || w < 0) return;

        switch (mode) {
            case 1:
                try_operation({
                    {EdgeAction::Add, w, z},
                    {EdgeAction::Remove, w, v},
                });
                break;
            case 2:
                try_operation({
                    {EdgeAction::Add, u, z},
                    {EdgeAction::Remove, w, u},
                });
                break;
            case 3:
                try_operation({
                    {EdgeAction::Add, u, z},
                    {EdgeAction::Add, w, z},
                    {EdgeAction::Remove, u, w},
                    {EdgeAction::Remove, v, z},
                });
                break;
            case 4:
                try_operation({
                    {EdgeAction::Add, u, v},
                    {EdgeAction::Add, w, z},
                    {EdgeAction::Remove, u, w},
                    {EdgeAction::Remove, v, z},
                });
                break;
            case 5:
                try_operation({
                    {EdgeAction::Add, u, v},
                    {EdgeAction::Add, w, z},
                    {EdgeAction::Remove, u, w},
                    {EdgeAction::Remove, v, w},
                });
                break;
            case 6:
                try_operation({
                    {EdgeAction::Add, u, z},
                    {EdgeAction::Remove, w, v},
                });
                break;
            case 7:
                try_operation({
                    {EdgeAction::Add, w, z},
                    {EdgeAction::Add, z, u},
                    {EdgeAction::Remove, w, v},
                    {EdgeAction::Remove, w, u},
                });
                break;
            case 8:
                try_operation({
                    {EdgeAction::Add, w, z},
                    {EdgeAction::Remove, v, z},
                });
                break;
            case 9:
                try_operation({
                    {EdgeAction::Add, u, v},
                    {EdgeAction::Add, u, z},
                    {EdgeAction::Remove, v, z},
                    {EdgeAction::Remove, w, v},
                });
                break;
            case 29:
                try_operation({
                    {EdgeAction::Add, u, v},
                    {EdgeAction::Add, u, z},
                    {EdgeAction::Add, w, z},
                });
                break;
            case 34:
                try_operation({
                    {EdgeAction::Add, u, v},
                    {EdgeAction::Add, w, z},
                });
                break;
            case 35:
                try_operation({
                    {EdgeAction::Add, u, z},
                    {EdgeAction::Add, w, z},
                });
                break;
            case 37:
                try_operation({{EdgeAction::Add, u, z}});
                break;
            case 38:
                try_operation({{EdgeAction::Add, u, v}});
                break;
            case 52:
                if (graph_.degree[z] != 1) {
                    try_operation({{EdgeAction::Remove, z, v}});
                }
                break;
            default:
                break;
        }
    }

    void apply_from_claw(
        int mode, const array<int, 5>& nodes, const Sub4& sub) {
        array<int, 3> leaves{};
        int leaf_count = 0;
        int center = -1;
        for (int i = 0; i < 4; ++i) {
            const int degree = sub_degree(sub, i);
            if (degree == 3) {
                center = i;
            } else if (degree == 1 && leaf_count < 3) {
                leaves[leaf_count++] = i;
            }
        }
        if (center < 0 || leaf_count != 3) return;

        const int u = nodes[leaves[0]];
        const int v = nodes[leaves[1]];
        const int w = nodes[leaves[2]];
        const int z = nodes[center];

        switch (mode) {
            case 18:
                try_operation({
                    {EdgeAction::Add, w, u},
                    {EdgeAction::Add, u, v},
                    {EdgeAction::Remove, z, w},
                    {EdgeAction::Remove, z, v},
                });
                break;
            case 19:
                try_operation({
                    {EdgeAction::Add, u, v},
                    {EdgeAction::Add, w, v},
                    {EdgeAction::Remove, w, z},
                    {EdgeAction::Remove, z, v},
                });
                break;
            case 20:
                try_operation({
                    {EdgeAction::Add, v, w},
                    {EdgeAction::Remove, z, w},
                });
                break;
            case 30:
                try_operation({
                    {EdgeAction::Add, u, v},
                    {EdgeAction::Add, w, v},
                    {EdgeAction::Add, u, w},
                });
                break;
            case 36:
                try_operation({
                    {EdgeAction::Add, u, v},
                    {EdgeAction::Add, w, v},
                });
                break;
            case 39:
                try_operation({{EdgeAction::Add, u, v}});
                break;
            case 53:
                if (graph_.degree[w] != 1) {
                    try_operation({{EdgeAction::Remove, z, w}});
                }
                break;
            case 54:
                if (graph_.degree[w] != 1) {
                    try_operation({
                        {EdgeAction::Add, u, v},
                        {EdgeAction::Remove, z, w},
                    });
                }
                break;
            default:
                break;
        }
    }
};

vector<vector<int>> connected_components(const Graph& graph) {
    vector<uint8_t> visited(graph.n, 0);
    vector<vector<int>> components;
    vector<int> stack;
    stack.reserve(graph.n);

    for (int start = 0; start < graph.n; ++start) {
        if (visited[start]) continue;

        vector<int> component;
        stack.clear();
        stack.push_back(start);
        visited[start] = 1;

        while (!stack.empty()) {
            const int u = stack.back();
            stack.pop_back();
            component.push_back(u);
            for (int v : graph.adj[u]) {
                if (!visited[v]) {
                    visited[v] = 1;
                    stack.push_back(v);
                }
            }
        }
        components.push_back(move(component));
    }
    return components;
}

void reconnect_graph(Graph& graph, FastRng& rng) {
    vector<vector<int>> components = connected_components(graph);
    if (components.size() <= 1) return;

    const size_t initial_count = components.size();
    cerr << "Graph is disconnected: " << initial_count
         << " connected components.\n";

    while (components.size() > 1) {
        const auto largest = max_element(
            components.begin(), components.end(),
            [](const vector<int>& lhs, const vector<int>& rhs) {
                return lhs.size() < rhs.size();
            });
        iter_swap(components.begin(), largest);

        const size_t other_index = 1 +
            rng.bounded(static_cast<uint64_t>(components.size() - 1));
        const int x = components[0][rng.bounded(components[0].size())];
        const int y = components[other_index][
            rng.bounded(components[other_index].size())];

        graph.add_edge_unchecked(x, y);
        components[0].insert(
            components[0].end(),
            components[other_index].begin(),
            components[other_index].end());

        if (other_index + 1 != components.size()) {
            swap(components[other_index], components.back());
        }
        components.pop_back();
    }

    cerr << "Added " << (initial_count - 1)
         << " bridge edges to reconnect the graph.\n";
}

int parse_int(const char* text, const char* name) {
    try {
        size_t used = 0;
        const string value(text);
        const long long parsed = stoll(value, &used);
        if (used != value.size() ||
            parsed < numeric_limits<int>::min() ||
            parsed > numeric_limits<int>::max()) {
            throw invalid_argument("out of range");
        }
        return static_cast<int>(parsed);
    } catch (const exception&) {
        throw runtime_error(string("invalid integer for ") + name + ": '" + text + "'");
    }
}

}  // namespace

int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc < 10) {
        cerr << "Usage: " << argv[0]
             << " <target_file> <synth_file> <num_samples> <num_samples2>"
             << " <mode> <mode2> <ignored_evb> <degbound_code> <hdgbound>"
             << " [completed_samples_file]\n";
        return 1;
    }

    try {
        const string target_file = argv[1];
        const string synth_file = argv[2];
        const int num_samples = parse_int(argv[3], "num_samples");
        const int num_samples2 = parse_int(argv[4], "num_samples2");
        const int mode1 = parse_int(argv[5], "mode");
        const int mode2 = parse_int(argv[6], "mode2");

        // argv[7] intentionally remains in the interface so existing rpll.sh
        // calls continue to work. The eigenvector system and evb constraint are
        // completely removed; this text is neither parsed nor used.
        (void)argv[7];

        const int degbound_code = parse_int(argv[8], "degbound_code");
        const int hdgbound = parse_int(argv[9], "hdgbound");
        if (degbound_code < 0 || hdgbound < 0) {
            throw runtime_error("degbound_code and hdgbound must be nonnegative");
        }
        const double degree_log_bound =
            static_cast<double>(degbound_code) / 140.0;
        const string completed_samples_file =
            argc >= 11 ? argv[10] : "expfy_samples_done.tmp";

        const EdgeListData target_data = read_edge_list(target_file);
        const EdgeListData synth_data = read_edge_list(synth_file);
        if (target_data.n != synth_data.n) {
            throw runtime_error(
                "target has " + to_string(target_data.n) +
                " nodes but synthetic graph has " +
                to_string(synth_data.n));
        }

        Graph graph(synth_data);
        const vector<int> target_degree =
            match_target_degrees_by_current_rank(
                target_data.degree, graph.degree);
        const vector<uint8_t> good_mask =
            read_good_node_mask(graph.n);

        const uint64_t seed = make_seed();
        Sampler sampler(
            graph,
            target_degree,
            degree_log_bound,
            hdgbound,
            good_mask,
            seed);
        const SampleStats stats =
            sampler.run(num_samples, num_samples2, mode1, mode2);

        ofstream completed_output(completed_samples_file);
        if (!completed_output) {
            throw runtime_error(
                "cannot write completed sample counts to '" +
                completed_samples_file + "'");
        }
        completed_output << stats.done1 << ' ' << stats.done2 << '\n';

        // Preserve the current expfy behavior: connectivity is repaired only
        // after sampling, and the bridge edges are not counted as transformations.
        FastRng reconnect_rng(seed ^ 0xd1b54a32d192ed03ULL);
        reconnect_graph(graph, reconnect_rng);

        for (int u = 0; u < graph.n; ++u) {
            for (int v : graph.adj[u]) {
                if (u < v) cout << u << ' ' << v << '\n';
            }
        }
    } catch (const exception& error) {
        cerr << "Error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
