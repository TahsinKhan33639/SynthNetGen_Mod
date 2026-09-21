// expfyb: modes 1--83 retained; learned inverses 84--95 added.
// Mode 73: compare existing G6 paths against hypothetical G5->G6 paths.
// Native continuous model; automatic training per invocation; C++17.
// Cutoff: maximize calibration existing-path gate with source precision >= target.
// The source labels do not indicate whether the reverse edit improves graphlets.
// T73--T95: exact graphlet verification and auditing are absent.
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
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

#include <filesystem>
#include <iomanip>
#include <map>
// Mode 73 collects, trains, and predicts natively within this executable.

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
// Modes 61--64 retain their conservative rooted-statistics work limit.
// The clique rewrites 65--67 do not use this work or these caches.
#ifndef EXPFY_G67_WORK_LIMIT
#define EXPFY_G67_WORK_LIMIT 200000LL
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
constexpr int MAX_TRANSFORMATION = 95;
// Mode 69 reads at most 32 + 31 neighbor degrees per proposal.
constexpr int PATH_LEAF_ROOT_DEGREE_LIMIT = 32;
// Mode 70 bounds both its branch checks and its common-neighbor scans.
constexpr int BOX_SWITCH_MAX_BRANCH_DEGREE = 32;
constexpr int BOX_SWITCH_NEIGHBOR_WORK_LIMIT = 128;
// Mode 71 uses four bounded common-neighbor checks and one bounded box sum.
constexpr int PAW_BOX_MAX_LOCAL_DEGREE = 32;
constexpr int PAW_BOX_NEIGHBOR_WORK_LIMIT = 128;
constexpr int TRIANGLE_BRANCH_MAX_DEGREE = 32;
constexpr long long G67_WORK_LIMIT = EXPFY_G67_WORK_LIMIT;
static_assert(G67_WORK_LIMIT > 0, "EXPFY_G67_WORK_LIMIT must be positive");
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

constexpr array<int8_t, MAX_TRANSFORMATION + 1> make_required_gid_table() {
    array<int8_t, MAX_TRANSFORMATION + 1> table{};
    for (int i = 0; i <= MAX_TRANSFORMATION; ++i) table[i] = -1;

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

    table[61] = 6;  // Leaf slide selected from an induced P4.
    table[62] = 7;  // Paw-tail relocation.
    table[63] = 7;  // Open a pendant triangle.
    table[64] = 6;  // P4 leaf slide plus a disjoint pendant triangle opening.
    table[65] = 10; // Exchange a clique edge across a shared triangle.
    table[66] = 10; // Switch edges of a clique and an adjoining box.
    table[67] = 10; // Move a clique edge into an adjoining box diagonal.
    table[68] = 6;  // Triangle-guided relocation of a true P4 endpoint leaf.
    table[69] = 6;  // Degree-histogram-preserving P4 reduction.
    table[70] = 6;  // Degree-preserving leaf/branch switch creating boxes.
    table[71] = 7;  // Open an isolated triangle into boxes using a leaf swap.
    table[72] = 7;  // Open a triangle and exchange a branch, keeping hub degrees.
    table[73] = 6;  // Convert an existing P4 that scores like a star-derived P4.
    table[74] = 5;  // Learned inverse of T9: star -> path.
    table[75] = 6;  // Learned inverse of T4: degree-preserving path switch.
    table[76] = 7;  // Learned inverse of T22: paw -> cycle.
    table[77] = 7;  // Learned inverse of T23: paw -> cycle.
    table[78] = 6;  // Learned inverse of T1: path -> path.
    table[79] = 5;  // Learned inverse of T8: star -> path.
    table[80] = 8;  // Learned inverse of T16: cycle -> paw.
    table[81] = 7;  // Learned inverse of T13: paw -> paw.
    table[82] = 7;  // Learned inverse of T10: paw -> paw.
    table[83] = 7;  // Learned inverse of T14: paw -> paw.
    table[84] = 6;  // Learned inverse of T2: path -> path.
    table[85] = 6;  // Learned inverse of T3: path -> path.
    table[86] = 6;  // Learned inverse of T5: path -> path.
    table[87] = 6;  // Learned inverse of T6: path -> path.
    table[88] = 5;  // Learned inverse of T7: star -> path.
    table[89] = 7;  // Learned inverse of T11: paw -> paw.
    table[90] = 7;  // Learned inverse of T12: paw -> paw.
    table[91] = 7;  // Learned inverse of T15: paw -> paw.
    table[92] = 8;  // Learned inverse of T17: cycle -> paw.
    table[93] = 5;  // Learned inverse of T18: star -> star.
    table[94] = 6;  // Learned inverse of T19: path -> star.
    table[95] = 8;  // Learned inverse of T21: cycle -> cycle.
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

vector<uint8_t> read_good_node_mask(
    int n, vector<uint8_t>& clique_sources, vector<long long>& triangle_counts) {
    ifstream input("data_middle1.txt");
    if (!input) {
        throw runtime_error("cannot open data_middle1.txt");
    }

    vector<uint8_t> mask(n, 0);
    clique_sources.assign(n, 0);
    triangle_counts.assign(n, 0);
    for (int row = 0; row < n; ++row) {
        long long degree = 0, triangles = 0;
        for (int column = 0; column < ORCA4_COLUMNS; ++column) {
            long long value = 0;
            if (!(input >> value)) {
                throw runtime_error(
                    "data_middle1.txt ended before " + to_string(n) +
                    " complete 15-column rows were read");
            }
            if (column == GOOD_NODE_COLUMN && value > 0) mask[row] = 1;
            if (column == 0) degree = value;
            if (column == 3) triangles = value;
        }
        triangle_counts[row] = triangles;
        // ORCA orbit 3 counts triangles at a vertex; orbit 0 is its degree.
        // High local closure favors deleting clique edges that create fewer
        // induced boxes. This input-derived heuristic is fixed for the run;
        // proposals use O(1) lookups and never scan common neighborhoods.
        const long double pairs = static_cast<long double>(degree) * (static_cast<long double>(degree) - 1) / 2;
        clique_sources[row] = degree >= 3 && triangles >= 0.75L * pairs;
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

// ---------------------------------------------------------------------------
// Mode 73: native smooth source classifier. No Python or generated model header.
// The labels describe provenance, NOT whether an edit improves graphlet counts.
// ---------------------------------------------------------------------------
using DegreeTuple73 = array<int, 4>;

// Star C-A, C-B, C-D -> path A-B-C-D by removing C-A and adding A-B.
// Input: (degree(A), degree(B), degree(C), degree(D)), with C the center.
// Output is hypothetical: no live graph is modified during collection.
DegreeTuple73 path_from_star73(const DegreeTuple73& star) {
    if (star[0] < 1 || star[1] < 1 || star[2] < 3 || star[3] < 1 ||
        star[1] == numeric_limits<int>::max())
        throw runtime_error("invalid degree tuple for an induced star");
    return {star[0], star[1]+1, star[2]-1, star[3]};
}

struct Example73 {
    DegreeTuple73 degree{}; // Ordered path A-B-C-D degrees; never sorted.
    bool from_star = false;  // false=defdeg; true=moddeg, the desired source.
    int partition = 0;  // 0=fit, 1=cutoff calibration, 2=held-out diagnostic
};

long long env_integer73(const char* name, long long fallback,
                        long long minimum, long long maximum) {
    const char* text = getenv(name);
    if (!text) return fallback;
    try {
        size_t used = 0;
        const string value(text);
        const long long parsed = stoll(value, &used);
        if (used != value.size() || parsed < minimum || parsed > maximum)
            throw invalid_argument("out of range");
        return parsed;
    } catch (const exception&) {
        throw runtime_error(string("invalid ") + name + "='" + text + "'");
    }
}

double env_real73(const char* name, double fallback,
                   double minimum, double maximum) {
    const char* text = getenv(name);
    if (!text) return fallback;
    try {
        size_t used = 0;
        const string value(text);
        const double parsed = stod(value, &used);
        if (used != value.size() || !isfinite(parsed) ||
            parsed < minimum || parsed > maximum)
            throw invalid_argument("out of range");
        return parsed;
    } catch (const exception&) {
        throw runtime_error(string("invalid ") + name + "='" + text + "'");
    }
}

// T73--T95 have no graphlet-effect counting, rejection, or audit path.
// Obsolete environment settings cannot re-enable that removed functionality.
void warn_removed_effect_options(int mode) {
    const string prefix = "EXPFY" + to_string(mode) + "_";
    for (const char* suffix : {"VERIFY", "AUDIT_EDITS", "DELTA_WORK",
                               "INCREASE_GID", "DECREASE_GID"}) {
        const string local = prefix + suffix;
        const string shared = string("EXPFY_ML_") + suffix;
        const string legacy = string("EXPFY73_") + suffix;
        if (getenv(local.c_str()) || getenv(shared.c_str()) || getenv(legacy.c_str())) {
            cerr << '[' << mode << "] NOTE: exact graphlet verification and auditing "
                    "were removed; old effect-check options are ignored.\n";
            return;
        }
    }
}

struct Config73 {
    long long data_attempts = 50000;
    double min_source_precision = 0.90;  // Flagged star-derived paths / all flagged source rows.
    double C = 0.5;
    int max_iterations = 500;
    bool both_orientations = false;
    string dump_directory;

    static Config73 from_environment() {
        Config73 c;
        warn_removed_effect_options(73);
        c.data_attempts = env_integer73("EXPFY73_DATA_ATTEMPTS", 50000, 100, 2000000);
        c.min_source_precision = env_real73("EXPFY73_MIN_PRECISION", 0.90, 0.0, 1.0);
        if (getenv("EXPFY73_MIN_PATH_RATE") || getenv("EXPFY73_PATH_RATE")) {
            cerr << "[73] WARNING: EXPFY73_MIN_PATH_RATE and EXPFY73_PATH_RATE are ignored. "
                    "This version maximizes path_gate_rate subject to the precision target; "
                    "there is no path-rate minimum or maximum.\n";
        }
        c.C = env_real73("EXPFY73_C", 0.5, 1e-8, 1e8);
        c.max_iterations = static_cast<int>(
            env_integer73("EXPFY73_MAX_ITERS", 500, 1, 10000));
        c.both_orientations = env_integer73("EXPFY73_BOTH", 0, 0, 1) != 0;
        if (const char* value = getenv("EXPFY73_DUMP_DIR")) c.dump_directory = value;
        return c;
    }
};

uint64_t mix73(uint64_t x) {
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

int grouped_partition73(DegreeTuple73 t, uint64_t seed) {
    DegreeTuple73 reversed = {t[3], t[2], t[1], t[0]};
    if (reversed < t) t = reversed;
    uint64_t h = mix73(seed ^ 0x6a09e667f3bcc909ULL);
    for (int value : t) h = mix73(h ^ static_cast<uint64_t>(value));
    const int bucket = static_cast<int>(h % 10);
    return bucket < 6 ? 0 : bucket < 8 ? 1 : 2;
}

class SmoothSourceModel73 {
public:
    static constexpr int VARIABLES = 10;
    static constexpr int USED_BASES = 6;  // 5 uniform knots + degree 3 - 2
    static constexpr int FEATURES = VARIABLES * USED_BASES;
    static constexpr int PARAMETERS = FEATURES + 1;
    using Feature = array<double, FEATURES>;
    using Parameter = array<double, PARAMETERS>;

    array<double, VARIABLES> lower{};
    array<double, VARIABLES> upper{};
    array<double, VARIABLES> inverse_step{};
    Feature weight{};
    double intercept = 0.0;
    double cutoff = numeric_limits<double>::infinity();
    bool fitted = false;
    bool converged = false;
    int iterations = 0;
    double final_gradient = numeric_limits<double>::infinity();
    double final_objective = numeric_limits<double>::infinity();

    static array<double, VARIABLES> variables_from_logs(const array<double, 4>& l) {
        return {l[0], l[1], l[2], l[3],
                l[0]-l[1], l[0]-l[2], l[0]-l[3],
                l[1]-l[2], l[1]-l[3], l[2]-l[3]};
    }

    static array<double, VARIABLES> variables(const DegreeTuple73& t) {
        array<double, 4> l{};
        for (int i = 0; i < 4; ++i) {
            if (t[i] <= 0) throw runtime_error("mode 73 requires positive degrees");
            l[i] = log(static_cast<double>(t[i]));
        }
        return variables_from_logs(l);
    }

    // Uniform cubic B-splines over 5 knots, constant extrapolation, last basis
    // dropped. Four basis functions can be nonzero in a single interval.
    array<double, USED_BASES> basis(int j, double value) const {
        array<double, USED_BASES> b{};
        const double q = clamp((value-lower[j]) * inverse_step[j], 0.0, 4.0);
        const int interval = min(3, static_cast<int>(q));
        const double t = q - interval;
        const double u = 1.0-t;
        const double t2 = t*t, t3 = t2*t;
        const array<double, 4> local = {
            u*u*u / 6.0,
            (3.0*t3-6.0*t2+4.0) / 6.0,
            (-3.0*t3+3.0*t2+3.0*t+1.0) / 6.0,
            t3 / 6.0
        };
        for (int k = 0; k < 4; ++k)
            if (interval+k < USED_BASES) b[interval+k] = local[k];
        return b;
    }

    Feature features(const DegreeTuple73& t) const {
        const auto v = variables(t);
        Feature out{};
        for (int j = 0; j < VARIABLES; ++j) {
            const auto b = basis(j, v[j]);
            for (int k = 0; k < USED_BASES; ++k)
                out[j*USED_BASES+k] = b[k];
        }
        return out;
    }

    static double sigmoid(double s) {
        if (s >= 0) return 1.0/(1.0+exp(-s));
        const double e = exp(s);
        return e/(1.0+e);
    }

    // Larger scores mean more similar to a star-derived path. These are
    // source-classification scores, not probabilities of a beneficial edit.
    double score_variables(const array<double, VARIABLES>& v) const {
        double s = intercept;
        for (int j = 0; j < VARIABLES; ++j) {
            const double q = clamp((v[j]-lower[j])*inverse_step[j], 0.0, 4.0);
            const int interval = min(3, static_cast<int>(q));
            const double t = q-interval, u = 1.0-t;
            const double t2 = t*t, t3 = t2*t;
            const array<double, 4> b = {
                u*u*u/6.0, (3.0*t3-6.0*t2+4.0)/6.0,
                (-3.0*t3+3.0*t2+3.0*t+1.0)/6.0, t3/6.0
            };
            for (int k = 0; k < 4; ++k)
                if (interval+k < USED_BASES)
                    s += weight[j*USED_BASES+interval+k]*b[k];
        }
        return s;
    }

    double score(const DegreeTuple73& t) const {
        return score_variables(variables(t));
    }

    double score_logs(const array<double, 4>& l) const {
        return score_variables(variables_from_logs(l));
    }

    bool train(const vector<Example73>& examples, double C, int max_iterations) {
        fitted = false;
        converged = false;
        iterations = 0;
        // Aggregate exact duplicates while preserving every observation's
        // weight, including conflicting labels on the same degree tuple.
        map<DegreeTuple73, array<long long, 2>> counts;
        long long total = 0, positives = 0;
        for (const auto& e : examples) {
            if (e.partition != 0) continue;
            counts[e.degree][e.from_star ? 1 : 0]++;
            ++total;
            positives += e.from_star;
        }
        if (positives < 20 || total-positives < 20 || counts.size() < 10) return false;

        lower.fill(numeric_limits<double>::infinity());
        upper.fill(-numeric_limits<double>::infinity());
        for (const auto& entry : counts) {
            const auto v = variables(entry.first);
            for (int j = 0; j < VARIABLES; ++j) {
                lower[j] = min(lower[j], v[j]);
                upper[j] = max(upper[j], v[j]);
            }
        }
        for (int j = 0; j < VARIABLES; ++j) {
            // Degenerate coordinates contribute no fitted information. Give
            // their basis a finite interval, then zero their constant columns.
            if (upper[j]-lower[j] < 1e-12) upper[j] = lower[j]+1.0;
            inverse_step[j] = 4.0/(upper[j]-lower[j]);
        }

        struct Row { Feature x; double mass; double target; };
        vector<Row> rows;
        rows.reserve(counts.size());
        Feature mean{}, scale{};
        for (const auto& entry : counts) {
            const double n = static_cast<double>(entry.second[0]+entry.second[1]);
            Row row{features(entry.first), n/total, entry.second[1]/n};
            for (int j = 0; j < FEATURES; ++j) mean[j] += row.mass*row.x[j];
            rows.push_back(row);
        }
        for (const auto& row : rows)
            for (int j = 0; j < FEATURES; ++j) {
                const double d = row.x[j]-mean[j];
                scale[j] += row.mass*d*d;
            }
        for (int j = 0; j < FEATURES; ++j)
            scale[j] = scale[j] < 1e-24 ? 1.0 : sqrt(scale[j]);
        for (auto& row : rows)
            for (int j = 0; j < FEATURES; ++j)
                row.x[j] = (row.x[j]-mean[j])/scale[j];

        const double lambda = 1.0/(C*static_cast<double>(total));
        auto objective = [&](const Parameter& p, Parameter& gradient) {
            gradient.fill(0.0);
            double value = 0.0;
            for (const Row& row : rows) {
                double s = p[FEATURES];
                for (int j = 0; j < FEATURES; ++j) s += p[j]*row.x[j];
                const double loss = s >= 0 ? (1.0-row.target)*s+log1p(exp(-s))
                                           : -row.target*s+log1p(exp(s));
                value += row.mass*loss;
                const double residual = row.mass*(sigmoid(s)-row.target);
                for (int j = 0; j < FEATURES; ++j) gradient[j] += residual*row.x[j];
                gradient[FEATURES] += residual;
            }
            for (int j = 0; j < FEATURES; ++j) {
                value += 0.5*lambda*p[j]*p[j];
                gradient[j] += lambda*p[j];
            }
            return value;
        };
        auto dot = [](const Parameter& a, const Parameter& b) {
            double s = 0;
            for (int j = 0; j < PARAMETERS; ++j) s += a[j]*b[j];
            return s;
        };
        auto infinity_norm = [](const Parameter& p) {
            double result = 0;
            for (double value : p) result = max(result, abs(value));
            return result;
        };

        Parameter p{}, gradient{};
        p[FEATURES] = log(static_cast<double>(positives)/(total-positives));
        double value = objective(p, gradient);
        vector<Parameter> history_s, history_y;
        vector<double> history_rho;
        constexpr size_t MEMORY = 12;

        for (int iteration = 0; iteration < max_iterations; ++iteration) {
            final_gradient = infinity_norm(gradient);
            if (final_gradient < 1e-7) { converged = true; break; }
            Parameter q = gradient;
            vector<double> alpha(history_s.size());
            for (int k = static_cast<int>(history_s.size())-1; k >= 0; --k) {
                alpha[k] = history_rho[k]*dot(history_s[k], q);
                for (int j = 0; j < PARAMETERS; ++j) q[j] -= alpha[k]*history_y[k][j];
            }
            double initial_scale = 1.0;
            if (!history_s.empty()) {
                const double yy = dot(history_y.back(), history_y.back());
                if (yy > 0) initial_scale = dot(history_s.back(), history_y.back())/yy;
            }
            Parameter direction{};
            for (int j = 0; j < PARAMETERS; ++j) direction[j] = q[j]*initial_scale;
            for (size_t k = 0; k < history_s.size(); ++k) {
                const double beta = history_rho[k]*dot(history_y[k], direction);
                for (int j = 0; j < PARAMETERS; ++j)
                    direction[j] += history_s[k][j]*(alpha[k]-beta);
            }
            for (double& d : direction) d = -d;
            double slope = dot(gradient, direction);
            if (!isfinite(slope) || slope >= -1e-20) {
                history_s.clear(); history_y.clear(); history_rho.clear();
                for (int j = 0; j < PARAMETERS; ++j) direction[j] = -gradient[j];
                slope = -dot(gradient, gradient);
            }
            Parameter next{}, next_gradient{};
            double step = 1.0, next_value = value;
            bool accepted = false;
            for (int backtrack = 0; backtrack < 50; ++backtrack) {
                for (int j = 0; j < PARAMETERS; ++j) next[j] = p[j]+step*direction[j];
                next_value = objective(next, next_gradient);
                if (isfinite(next_value) && next_value <= value+1e-4*step*slope) {
                    accepted = true;
                    break;
                }
                step *= 0.5;
            }
            if (!accepted) break;
            Parameter s{}, y{};
            for (int j = 0; j < PARAMETERS; ++j) {
                s[j] = next[j]-p[j];
                y[j] = next_gradient[j]-gradient[j];
            }
            const double sy = dot(s,y);
            if (sy > 1e-14*sqrt(dot(s,s)*dot(y,y))) {
                if (history_s.size() == MEMORY) {
                    history_s.erase(history_s.begin());
                    history_y.erase(history_y.begin());
                    history_rho.erase(history_rho.begin());
                }
                history_s.push_back(s); history_y.push_back(y); history_rho.push_back(1.0/sy);
            }
            p = next; gradient = next_gradient; value = next_value;
            iterations = iteration+1;
        }
        final_gradient = infinity_norm(gradient);
        final_objective = value;
        converged = final_gradient < 1e-7;
        if (!isfinite(value) || !isfinite(final_gradient)) return false;
        intercept = p[FEATURES];
        for (int j = 0; j < FEATURES; ++j) {
            weight[j] = p[j]/scale[j];
            intercept -= weight[j]*mean[j];
        }
        fitted = true;
        return true;
    }

    void save(ostream& out) const {
        out << "EXPFY_NATIVE_LINE73_V1\n" << setprecision(17);
        out << intercept << ' ' << cutoff << '\n';
        for (int j = 0; j < VARIABLES; ++j)
            out << lower[j] << ' ' << upper[j] << '\n';
        for (double w : weight) out << w << '\n';
    }
};

// A row is one sampled example, not a unique tuple. Repeated rows keep their
// observation weight. Larger scores mean more star-derived-path-like.
struct ScoredSource73 {
    double score = 0.0;
    bool from_star = false;  // false=existing G6 (defdeg); true=star-derived G6 (moddeg).
};

struct ThresholdPoint73 {
    double cutoff = numeric_limits<double>::infinity();
    long long path_calls = 0;
    long long star_calls = 0;

    long long calls() const { return path_calls + star_calls; }
    double precision() const {
        return calls() ? static_cast<double>(star_calls)/calls()
                       : numeric_limits<double>::quiet_NaN();
    }
};

struct ThresholdChoice73 {
    // feasible means at least one calibration path passes at the requested
    // source precision. A best possible gate of zero leaves mode 73 disabled.
    bool feasible = false;
    long long path_total = 0;
    long long star_total = 0;
    ThresholdPoint73 chosen;
    bool has_best_with_precision = false;
    ThresholdPoint73 best_with_precision;
    // Diagnostic only: best precision among thresholds accepting any path.
    // This is never used as an automatic fallback.
    bool has_best_nonzero_path = false;
    ThresholdPoint73 best_nonzero_path;
    vector<ThresholdPoint73> curve;
};

ThresholdChoice73 choose_max_gate_cutoff73(vector<ScoredSource73> rows,
                                           double min_precision) {
    if (!isfinite(min_precision) || min_precision < 0.0 || min_precision > 1.0)
        throw runtime_error("mode73 source precision must lie in [0,1]");
    ThresholdChoice73 result;
    for (const auto& row : rows) {
        if (!isfinite(row.score))
            throw runtime_error("nonfinite mode73 calibration score; refusing to select a cutoff");
        if (row.from_star) ++result.star_total;
        else ++result.path_total;
    }
    // Both sources are necessary for a meaningful calibration experiment.
    if (result.path_total == 0 || result.star_total == 0) return result;
    sort(rows.begin(),rows.end(),[](const auto& a,const auto& b) {
        return a.score > b.score;
    });
    result.curve.reserve(rows.size());
    long long paths = 0, stars = 0;
    for (size_t i = 0; i < rows.size();) {
        size_t j = i;
        // A deterministic threshold cannot split a tied score group.
        while (j < rows.size() && rows[j].score == rows[i].score) {
            if (rows[j].from_star) ++stars; else ++paths;
            ++j;
        }
        double cutoff = rows[i].score;
        // Place the cutoff in the gap where representable; otherwise keep
        // the included score. Always reproduce this prefix exactly with >=.
        if (j < rows.size()) {
            const double lower = rows[j].score;
            const double midpoint = lower/2.0 + cutoff/2.0;
            if (midpoint > lower && midpoint <= cutoff) cutoff = midpoint;
        }
        const ThresholdPoint73 point{cutoff,paths,stars};
        result.curve.push_back(point);
        const double precision = point.precision();

        // Inspect ALL thresholds: precision need not be monotone. Do not
        // stop at the first valid threshold or at the first precision failure.
        // Primary objective: most paths. Ties: higher precision, then higher
        // cutoff (the latter is retained automatically by descending order).
        if (precision >= min_precision &&
            (!result.has_best_with_precision ||
             paths > result.best_with_precision.path_calls ||
             (paths == result.best_with_precision.path_calls &&
              precision > result.best_with_precision.precision()))) {
            result.has_best_with_precision = true;
            result.best_with_precision = point;
        }
        if (paths > 0 &&
            (!result.has_best_nonzero_path ||
             precision > result.best_nonzero_path.precision() ||
             (precision == result.best_nonzero_path.precision() &&
              paths > result.best_nonzero_path.path_calls))) {
            result.has_best_nonzero_path = true;
            result.best_nonzero_path = point;
        }
        i = j;
    }
    if (result.has_best_with_precision && result.best_with_precision.path_calls > 0) {
        result.feasible = true;
        result.chosen = result.best_with_precision;
    }
    return result;
}

void print_threshold_point73(ostream& out, const char* name,
                             const ThresholdPoint73& point,
                             long long path_total) {
    const double rate = path_total ? static_cast<double>(point.path_calls)/path_total
                                      : numeric_limits<double>::quiet_NaN();
    out << "[73] " << name << ": cutoff=" << point.cutoff
        << " defdeg_flagged=" << point.path_calls << '/' << path_total
        << " path_gate_rate=" << rate
        << " source_precision=" << point.precision()
        << " (moddeg_flagged=" << point.star_calls
        << ", all_flagged=" << point.calls() << ")\n";
}

void print_max_gate_choice73(ostream& out, const ThresholdChoice73& choice,
                              double min_precision) {
    out << "[73] cutoff_policy=max_path_gate_at_precision\n"
        << "[73] objective: maximize path_gate_rate subject to source_precision >= "
        << min_precision << "; no path-rate minimum or maximum\n"
        << "[73] calibration_paths=" << choice.path_total
        << " calibration_star_derived_paths=" << choice.star_total << '\n';
    if (choice.feasible) {
        print_threshold_point73(out,"CHOSEN calibration",choice.chosen,choice.path_total);
        out << "[73] Maximum qualifying path gate among all whole-score thresholds "
                "for this fitted model on calibration data. "
                "Ties prefer higher precision, then higher cutoff.\n"
            << "[73] Held-out/live precision and gate rates may differ.\n";
        if (choice.chosen.calls() < 30)
            out << "[73] WARNING: fewer than 30 flagged calibration rows; "
                    "source precision has little support (duplicates may reduce it further).\n";
    } else {
        out << "[73] DISABLED: no threshold accepts a calibration path while meeting "
                "the source-precision target. Precision was not relaxed.\n";
        if (choice.path_total == 0 || choice.star_total == 0)
            out << "[73] Both calibration sources must be nonempty.\n";
        else if (choice.has_best_with_precision) {
            out << "[73] Best qualifying calibration path_gate_rate=0. "
                    "Only star-derived-path rows can be flagged at the requested precision.\n";
            print_threshold_point73(out,"zero-path qualifying threshold (NOT deployed)",
                                    choice.best_with_precision,choice.path_total);
        }
        if (choice.has_best_nonzero_path)
            print_threshold_point73(out,"best precision with any path accepted (NOT used)",
                                    choice.best_nonzero_path,choice.path_total);
    }
}

struct Counters73 {
    long long training_attempts = 0;
    long long complete_training_samples = 0;
    long long sampling_attempts = 0;
    long long path_candidates = 0;
    long long orientations_scored = 0;
    long long model_pass = 0;
    long long model_reject = 0;
    long long structure_reject = 0;
    long long degree_reject = 0;
    long long histogram_reject = 0;
    long long committed = 0;
};

enum class OperationFailure73 { None, Structure, Degree, Histogram };

// ---------------------------------------------------------------------------
// Learned inverses 74--95. A,B,C,D are FIXED transformation roles, not sorted
// degrees. Role matching, hypothetical degrees, dumps, and commits all use the
// same masks and the same edge-operation specification.
// ---------------------------------------------------------------------------
using InverseTuple = array<int, 4>;
struct InverseFitExample {
    InverseTuple degree{};
    bool positive = false;
    int partition = 0;
};

constexpr unsigned inverse_edge_bit(int a, int b) {
    return a > b ? inverse_edge_bit(b,a) :
           a == 0 ? (b == 1 ? 32U : b == 2 ? 16U : 8U) :
           a == 1 ? (b == 2 ? 4U : 2U) : 1U;
}
unsigned inverse_mask(const Sub4& sub) {
    unsigned result = 0;
    for (int i=0;i<4;++i) for (int j=i+1;j<4;++j)
        if (sub[i][j]) result |= inverse_edge_bit(i,j);
    return result;
}
unsigned inverse_permuted_mask(unsigned mask, const array<int,4>& p) {
    unsigned result=0;
    for (int i=0;i<4;++i) for (int j=i+1;j<4;++j)
        if (mask & inverse_edge_bit(p[i],p[j])) result |= inverse_edge_bit(i,j);
    return result;
}

struct InverseMatch {
    bool valid=false;
    array<int,4> role_index{};
};
struct InverseSpec {
    int mode=0, forward=0, source_gid=0, target_gid=0;
    unsigned source_mask=0, target_mask=0;
    array<pair<int,int>,2> forward_remove{}, forward_add{};
    int remove_count=0, add_count=0;
    InverseTuple forward_delta{};
    array<InverseMatch,64> source_matches{}, target_matches{};
    vector<array<int,4>> source_automorphisms;

    InverseSpec(int m, int f, initializer_list<pair<int,int>> source,
                initializer_list<pair<int,int>> remove,
                initializer_list<pair<int,int>> add):mode(m),forward(f) {
        if (remove.size()>2 || add.size()>2 || remove.size()+add.size()>4)
            throw runtime_error("inverse spec exceeds atomic operation capacity");
        for (const auto& e:source) source_mask |= inverse_edge_bit(e.first,e.second);
        target_mask=source_mask;
        for (const auto& e:remove) {
            const unsigned bit=inverse_edge_bit(e.first,e.second);
            if (!(target_mask&bit)) throw runtime_error("forward spec removes missing edge");
            target_mask &= ~bit;
            forward_remove[remove_count++]=e;
            --forward_delta[e.first]; --forward_delta[e.second];
        }
        for (const auto& e:add) {
            const unsigned bit=inverse_edge_bit(e.first,e.second);
            if (target_mask&bit) throw runtime_error("forward spec adds existing edge");
            target_mask |= bit;
            forward_add[add_count++]=e;
            ++forward_delta[e.first]; ++forward_delta[e.second];
        }
        source_gid=GRAPH_ID_BY_BITS[source_mask];
        target_gid=GRAPH_ID_BY_BITS[target_mask];
        if (source_gid<5 || target_gid<5)
            throw runtime_error("learned inverse templates must both be connected");
        for (unsigned mask=0;mask<64;++mask) {
            array<int,4> p={0,1,2,3};
            do {
                const unsigned mapped=inverse_permuted_mask(mask,p);
                if (mapped==source_mask && !source_matches[mask].valid)
                    source_matches[mask]={true,p};
                if (mapped==target_mask && !target_matches[mask].valid)
                    target_matches[mask]={true,p};
            } while (next_permutation(p.begin(),p.end()));
        }
        array<int,4> p={0,1,2,3};
        do {
            if (inverse_permuted_mask(source_mask,p)==source_mask)
                source_automorphisms.push_back(p);
        } while (next_permutation(p.begin(),p.end()));
    }
};

bool is_new_inverse(int mode) { return mode>=74 && mode<=95; }
const InverseSpec& inverse_spec(int mode) {
    // Fixed role conventions copied from the original forward transforms:
    // T1--T9:   A=u, B=w, C=v, D=z  (source path A-B-C-D).
    // T10--T17: A=u, B=v, C=z, D=w  (source paw: triangle A-B-C, tail C-D).
    // T18--T20: A=u, B=v, C=w, D=z  (source claw centered at D).
    // T21--T23: A=w, B=v, C=z, D=u  (source cycle A-B-C-D-A).
    static const array<InverseSpec,22> specs={
        InverseSpec(74,9, {{0,1},{1,2},{2,3}},
                         {{1,2},{2,3}}, {{0,2},{0,3}}),
        InverseSpec(75,4, {{0,1},{1,2},{2,3}},
                         {{0,1},{2,3}}, {{0,2},{1,3}}),
        InverseSpec(76,22,{{0,1},{1,2},{2,3},{3,0}},
                         {{1,2},{0,1}}, {{1,3},{0,2}}),
        InverseSpec(77,23,{{0,1},{1,2},{2,3},{3,0}},
                         {{0,1}}, {{1,3}}),
        InverseSpec(78,1, {{0,1},{1,2},{2,3}},
                         {{1,2}}, {{1,3}}),
        InverseSpec(79,8, {{0,1},{1,2},{2,3}},
                         {{2,3}}, {{1,3}}),
        InverseSpec(80,16,{{0,1},{0,2},{1,2},{2,3}},
                         {{2,3},{0,1}}, {{0,3},{1,3}}),
        InverseSpec(81,13,{{0,1},{0,2},{1,2},{2,3}},
                         {{2,3},{1,2}}, {{0,3},{1,3}}),
        InverseSpec(82,10,{{0,1},{0,2},{1,2},{2,3}},
                         {{0,1}}, {{1,3}}),
        InverseSpec(83,14,{{0,1},{0,2},{1,2},{2,3}},
                         {{1,2}}, {{0,3}}),

        // Remaining inverses from T1--T23, added in forward-transformation order.
        InverseSpec(84,2, {{0,1},{1,2},{2,3}},
                         {{0,1}}, {{0,3}}),
        InverseSpec(85,3, {{0,1},{1,2},{2,3}},
                         {{0,1},{2,3}}, {{0,3},{1,3}}),
        InverseSpec(86,5, {{0,1},{1,2},{2,3}},
                         {{0,1},{1,2}}, {{0,2},{1,3}}),
        InverseSpec(87,6, {{0,1},{1,2},{2,3}},
                         {{1,2}}, {{0,3}}),
        InverseSpec(88,7, {{0,1},{1,2},{2,3}},
                         {{1,2},{0,1}}, {{1,3},{0,3}}),
        InverseSpec(89,11,{{0,1},{0,2},{1,2},{2,3}},
                         {{0,2},{1,2}}, {{0,3},{1,3}}),
        InverseSpec(90,12,{{0,1},{0,2},{1,2},{2,3}},
                         {{0,1},{0,2}}, {{0,3},{1,3}}),
        InverseSpec(91,15,{{0,1},{0,2},{1,2},{2,3}},
                         {{2,3}}, {{0,3}}),
        InverseSpec(92,17,{{0,1},{0,2},{1,2},{2,3}},
                         {{0,2}}, {{0,3}}),
        InverseSpec(93,18,{{0,3},{1,3},{2,3}},
                         {{2,3},{1,3}}, {{0,2},{0,1}}),
        InverseSpec(94,19,{{0,3},{1,3},{2,3}},
                         {{2,3},{1,3}}, {{0,1},{1,2}}),
        InverseSpec(95,21,{{0,1},{1,2},{2,3},{3,0}},
                         {{2,3},{0,1}}, {{1,3},{0,2}})
    };
    if (!is_new_inverse(mode)) throw runtime_error("unknown learned inverse mode");
    return specs[mode-74];
}
InverseTuple inverse_shift(const InverseTuple& t, const InverseSpec& spec, int direction) {
    InverseTuple out{};
    for (int i=0;i<4;++i) {
        const long long value=static_cast<long long>(t[i])+direction*spec.forward_delta[i];
        if (value<=0 || value>numeric_limits<int>::max())
            throw runtime_error("nonpositive or overflowing hypothetical degree");
        out[i]=static_cast<int>(value);
    }
    return out;
}
array<int,4> inverse_role_nodes(const array<int,5>& nodes,const InverseMatch& match) {
    if (!match.valid) throw runtime_error("invalid inverse role match");
    return {nodes[match.role_index[0]],nodes[match.role_index[1]],
            nodes[match.role_index[2]],nodes[match.role_index[3]]};
}
InverseTuple inverse_degrees(const Graph& graph,const array<int,4>& nodes) {
    return {graph.degree[nodes[0]],graph.degree[nodes[1]],
            graph.degree[nodes[2]],graph.degree[nodes[3]]};
}
struct InverseOperations {
    array<EdgeOp,4> data{};
    int count=0;
};
InverseOperations inverse_operations(const InverseSpec& spec,const array<int,4>& nodes) {
    InverseOperations ops;
    for (int i=0;i<spec.add_count;++i) {
        const auto& e=spec.forward_add[i];
        ops.data[ops.count++]={EdgeAction::Remove,nodes[e.first],nodes[e.second]};
    }
    for (int i=0;i<spec.remove_count;++i) {
        const auto& e=spec.forward_remove[i];
        ops.data[ops.count++]={EdgeAction::Add,nodes[e.first],nodes[e.second]};
    }
    return ops;
}

// Specific EXPFY74_* ... EXPFY95_* overrides EXPFY_ML_*, which overrides compatible
// legacy EXPFY73_* options. BOTH is deliberately NOT inherited for new modes:
// the new modes use one randomized legal role assignment, not a max over roles.
string inverse_option(int mode,const char* suffix) {
    const string a="EXPFY"+to_string(mode)+"_"+suffix;
    if (getenv(a.c_str())) return a;
    const string b=string("EXPFY_ML_")+suffix;
    if (getenv(b.c_str())) return b;
    return string("EXPFY73_")+suffix;
}
struct InverseConfig {
    long long data_attempts=50000;
    double precision=0.9,C=0.5;
    int max_iterations=500;
    string dump_directory;
    static InverseConfig load(const InverseSpec& s) {
        InverseConfig c;
        warn_removed_effect_options(s.mode);
        auto integer=[&](const char* name,long long fallback,long long lo,long long hi) {
            return env_integer73(inverse_option(s.mode,name).c_str(),fallback,lo,hi);
        };
        auto real=[&](const char* name,double fallback,double lo,double hi) {
            return env_real73(inverse_option(s.mode,name).c_str(),fallback,lo,hi);
        };
        c.data_attempts=integer("DATA_ATTEMPTS",50000,100,2000000);
        c.precision=real("MIN_PRECISION",0.9,0.0,1.0);
        c.C=real("C",0.5,1e-8,1e8);
        c.max_iterations=static_cast<int>(integer("MAX_ITERS",500,1,10000));
        const string both="EXPFY"+to_string(s.mode)+"_BOTH";
        if (env_integer73(both.c_str(),0,0,1)!=0)
            throw runtime_error(both+" is not supported: new modes use one random legal role assignment");
        const string name=inverse_option(s.mode,"DUMP_DIR");
        if (const char* value=getenv(name.c_str())) {
            if (*value) c.dump_directory=(std::filesystem::path(value)/("mode"+to_string(s.mode))).string();
        }
        return c;
    }
};

struct InverseExample {
    InverseFitExample fit;
    InverseTuple raw{}, source_space{};
    array<int,4> nodes{};
    size_t group=0;
};
// Partition in the SAME source-space coordinates for b and c. The mapping
// between their descriptors is the fixed invertible shift forward_delta.
// Group automorphisms only for splitting; NEVER reorder model input degrees.
size_t inverse_assign_partitions(vector<InverseExample>& rows,const InverseSpec& spec,uint64_t seed) {
    const size_t n=rows.size();
    vector<size_t> parent(n),size(n,1); iota(parent.begin(),parent.end(),size_t{0});
    auto root=[&](size_t x) {
        while(parent[x]!=x) {parent[x]=parent[parent[x]];x=parent[x];} return x;
    };
    auto join=[&](size_t a,size_t b) {
        a=root(a);b=root(b);if(a==b)return;
        if(size[a]<size[b]) swap(a,b);
        parent[b]=a;size[a]+=size[b];
    };
    map<InverseTuple,size_t> degree_owner;
    map<array<int,4>,size_t> node_owner;
    vector<InverseTuple> keys(n);
    for(size_t i=0;i<n;++i) {
        auto key=rows[i].source_space;
        for(const auto& p:spec.source_automorphisms) {
            InverseTuple v={rows[i].source_space[p[0]],rows[i].source_space[p[1]],
                            rows[i].source_space[p[2]],rows[i].source_space[p[3]]};
            key=min(key,v);
        }
        keys[i]=key;
        auto a=degree_owner.emplace(key,i); if(!a.second)join(i,a.first->second);
        auto nodes=rows[i].nodes;sort(nodes.begin(),nodes.end());
        auto b=node_owner.emplace(nodes,i);if(!b.second)join(i,b.first->second);
    }
    map<size_t,InverseTuple> group_keys;
    for(size_t i=0;i<n;++i) {
        auto a=group_keys.emplace(root(i),keys[i]);if(!a.second)a.first->second=min(a.first->second,keys[i]);
    }
    for(size_t i=0;i<n;++i) {
        rows[i].group=root(i);
        uint64_t h=mix73(seed ^ 0x6a09e667f3bcc909ULL ^ static_cast<uint64_t>(spec.mode));
        for(int d:group_keys.at(root(i)))h=mix73(h^static_cast<uint64_t>(d));
        const int bucket=static_cast<int>(h%10);
        rows[i].fit.partition=bucket<6?0:bucket<8?1:2;
    }
    return group_keys.size();
}
class InverseSmoothModel {
public:
    static constexpr int VARIABLES = 10;
    static constexpr int USED_BASES = 6;  // 5 uniform knots + degree 3 - 2
    static constexpr int FEATURES = VARIABLES * USED_BASES;
    static constexpr int PARAMETERS = FEATURES + 1;
    using Feature = array<double, FEATURES>;
    using Parameter = array<double, PARAMETERS>;

    array<double, VARIABLES> lower{};
    array<double, VARIABLES> upper{};
    array<double, VARIABLES> inverse_step{};
    Feature weight{};
    double intercept = 0.0;
    double cutoff = numeric_limits<double>::infinity();
    bool fitted = false;
    bool converged = false;
    int iterations = 0;
    double final_gradient = numeric_limits<double>::infinity();
    double final_objective = numeric_limits<double>::infinity();

    static array<double, VARIABLES> variables_from_logs(const array<double, 4>& l) {
        return {l[0], l[1], l[2], l[3],
                l[0]-l[1], l[0]-l[2], l[0]-l[3],
                l[1]-l[2], l[1]-l[3], l[2]-l[3]};
    }

    static array<double, VARIABLES> variables(const InverseTuple& t) {
        array<double, 4> l{};
        for (int i = 0; i < 4; ++i) {
            if (t[i] <= 0) throw runtime_error("learned inverse requires positive degrees");
            l[i] = log(static_cast<double>(t[i]));
        }
        return variables_from_logs(l);
    }

    // Uniform cubic B-splines over 5 knots, constant extrapolation, last basis
    // dropped. Four basis functions can be nonzero in a single interval.
    array<double, USED_BASES> basis(int j, double value) const {
        array<double, USED_BASES> b{};
        const double q = clamp((value-lower[j]) * inverse_step[j], 0.0, 4.0);
        const int interval = min(3, static_cast<int>(q));
        const double t = q - interval;
        const double u = 1.0-t;
        const double t2 = t*t, t3 = t2*t;
        const array<double, 4> local = {
            u*u*u / 6.0,
            (3.0*t3-6.0*t2+4.0) / 6.0,
            (-3.0*t3+3.0*t2+3.0*t+1.0) / 6.0,
            t3 / 6.0
        };
        for (int k = 0; k < 4; ++k)
            if (interval+k < USED_BASES) b[interval+k] = local[k];
        return b;
    }

    Feature features(const InverseTuple& t) const {
        const auto v = variables(t);
        Feature out{};
        for (int j = 0; j < VARIABLES; ++j) {
            const auto b = basis(j, v[j]);
            for (int k = 0; k < USED_BASES; ++k)
                out[j*USED_BASES+k] = b[k];
        }
        return out;
    }

    static double sigmoid(double s) {
        if (s >= 0) return 1.0/(1.0+exp(-s));
        const double e = exp(s);
        return e/(1.0+e);
    }

    // Larger scores mean more similar to a positive-source example. These are
    // source-classification scores, not probabilities of a beneficial edit.
    double score_variables(const array<double, VARIABLES>& v) const {
        double s = intercept;
        for (int j = 0; j < VARIABLES; ++j) {
            const double q = clamp((v[j]-lower[j])*inverse_step[j], 0.0, 4.0);
            const int interval = min(3, static_cast<int>(q));
            const double t = q-interval, u = 1.0-t;
            const double t2 = t*t, t3 = t2*t;
            const array<double, 4> b = {
                u*u*u/6.0, (3.0*t3-6.0*t2+4.0)/6.0,
                (-3.0*t3+3.0*t2+3.0*t+1.0)/6.0, t3/6.0
            };
            for (int k = 0; k < 4; ++k)
                if (interval+k < USED_BASES)
                    s += weight[j*USED_BASES+interval+k]*b[k];
        }
        return s;
    }

    double score(const InverseTuple& t) const {
        return score_variables(variables(t));
    }

    double score_logs(const array<double, 4>& l) const {
        return score_variables(variables_from_logs(l));
    }

    bool train(const vector<InverseFitExample>& examples, double C, int max_iterations) {
        fitted = false;
        converged = false;
        iterations = 0;
        // Aggregate exact duplicates while preserving every observation's
        // weight, including conflicting labels on the same degree tuple.
        map<InverseTuple, array<long long, 2>> counts;
        long long total = 0, positives = 0;
        for (const auto& e : examples) {
            if (e.partition != 0) continue;
            counts[e.degree][e.positive ? 1 : 0]++;
            ++total;
            positives += e.positive;
        }
        if (positives < 20 || total-positives < 20 || counts.size() < 10) return false;

        lower.fill(numeric_limits<double>::infinity());
        upper.fill(-numeric_limits<double>::infinity());
        for (const auto& entry : counts) {
            const auto v = variables(entry.first);
            for (int j = 0; j < VARIABLES; ++j) {
                lower[j] = min(lower[j], v[j]);
                upper[j] = max(upper[j], v[j]);
            }
        }
        for (int j = 0; j < VARIABLES; ++j) {
            // Degenerate coordinates contribute no fitted information. Give
            // their basis a finite interval, then zero their constant columns.
            if (upper[j]-lower[j] < 1e-12) upper[j] = lower[j]+1.0;
            inverse_step[j] = 4.0/(upper[j]-lower[j]);
        }

        struct Row { Feature x; double mass; double target; };
        vector<Row> rows;
        rows.reserve(counts.size());
        Feature mean{}, scale{};
        for (const auto& entry : counts) {
            const double n = static_cast<double>(entry.second[0]+entry.second[1]);
            Row row{features(entry.first), n/total, entry.second[1]/n};
            for (int j = 0; j < FEATURES; ++j) mean[j] += row.mass*row.x[j];
            rows.push_back(row);
        }
        for (const auto& row : rows)
            for (int j = 0; j < FEATURES; ++j) {
                const double d = row.x[j]-mean[j];
                scale[j] += row.mass*d*d;
            }
        for (int j = 0; j < FEATURES; ++j)
            scale[j] = scale[j] < 1e-24 ? 1.0 : sqrt(scale[j]);
        for (auto& row : rows)
            for (int j = 0; j < FEATURES; ++j)
                row.x[j] = (row.x[j]-mean[j])/scale[j];

        const double lambda = 1.0/(C*static_cast<double>(total));
        auto objective = [&](const Parameter& p, Parameter& gradient) {
            gradient.fill(0.0);
            double value = 0.0;
            for (const Row& row : rows) {
                double s = p[FEATURES];
                for (int j = 0; j < FEATURES; ++j) s += p[j]*row.x[j];
                const double loss = s >= 0 ? (1.0-row.target)*s+log1p(exp(-s))
                                           : -row.target*s+log1p(exp(s));
                value += row.mass*loss;
                const double residual = row.mass*(sigmoid(s)-row.target);
                for (int j = 0; j < FEATURES; ++j) gradient[j] += residual*row.x[j];
                gradient[FEATURES] += residual;
            }
            for (int j = 0; j < FEATURES; ++j) {
                value += 0.5*lambda*p[j]*p[j];
                gradient[j] += lambda*p[j];
            }
            return value;
        };
        auto dot = [](const Parameter& a, const Parameter& b) {
            double s = 0;
            for (int j = 0; j < PARAMETERS; ++j) s += a[j]*b[j];
            return s;
        };
        auto infinity_norm = [](const Parameter& p) {
            double result = 0;
            for (double value : p) result = max(result, abs(value));
            return result;
        };

        Parameter p{}, gradient{};
        p[FEATURES] = log(static_cast<double>(positives)/(total-positives));
        double value = objective(p, gradient);
        vector<Parameter> history_s, history_y;
        vector<double> history_rho;
        constexpr size_t MEMORY = 12;

        for (int iteration = 0; iteration < max_iterations; ++iteration) {
            final_gradient = infinity_norm(gradient);
            if (final_gradient < 1e-7) { converged = true; break; }
            Parameter q = gradient;
            vector<double> alpha(history_s.size());
            for (int k = static_cast<int>(history_s.size())-1; k >= 0; --k) {
                alpha[k] = history_rho[k]*dot(history_s[k], q);
                for (int j = 0; j < PARAMETERS; ++j) q[j] -= alpha[k]*history_y[k][j];
            }
            double initial_scale = 1.0;
            if (!history_s.empty()) {
                const double yy = dot(history_y.back(), history_y.back());
                if (yy > 0) initial_scale = dot(history_s.back(), history_y.back())/yy;
            }
            Parameter direction{};
            for (int j = 0; j < PARAMETERS; ++j) direction[j] = q[j]*initial_scale;
            for (size_t k = 0; k < history_s.size(); ++k) {
                const double beta = history_rho[k]*dot(history_y[k], direction);
                for (int j = 0; j < PARAMETERS; ++j)
                    direction[j] += history_s[k][j]*(alpha[k]-beta);
            }
            for (double& d : direction) d = -d;
            double slope = dot(gradient, direction);
            if (!isfinite(slope) || slope >= -1e-20) {
                history_s.clear(); history_y.clear(); history_rho.clear();
                for (int j = 0; j < PARAMETERS; ++j) direction[j] = -gradient[j];
                slope = -dot(gradient, gradient);
            }
            Parameter next{}, next_gradient{};
            double step = 1.0, next_value = value;
            bool accepted = false;
            for (int backtrack = 0; backtrack < 50; ++backtrack) {
                for (int j = 0; j < PARAMETERS; ++j) next[j] = p[j]+step*direction[j];
                next_value = objective(next, next_gradient);
                if (isfinite(next_value) && next_value <= value+1e-4*step*slope) {
                    accepted = true;
                    break;
                }
                step *= 0.5;
            }
            if (!accepted) break;
            Parameter s{}, y{};
            for (int j = 0; j < PARAMETERS; ++j) {
                s[j] = next[j]-p[j];
                y[j] = next_gradient[j]-gradient[j];
            }
            const double sy = dot(s,y);
            if (sy > 1e-14*sqrt(dot(s,s)*dot(y,y))) {
                if (history_s.size() == MEMORY) {
                    history_s.erase(history_s.begin());
                    history_y.erase(history_y.begin());
                    history_rho.erase(history_rho.begin());
                }
                history_s.push_back(s); history_y.push_back(y); history_rho.push_back(1.0/sy);
            }
            p = next; gradient = next_gradient; value = next_value;
            iterations = iteration+1;
        }
        final_gradient = infinity_norm(gradient);
        final_objective = value;
        converged = final_gradient < 1e-7;
        if (!isfinite(value) || !isfinite(final_gradient)) return false;
        intercept = p[FEATURES];
        for (int j = 0; j < FEATURES; ++j) {
            weight[j] = p[j]/scale[j];
            intercept -= weight[j]*mean[j];
        }
        fitted = true;
        return true;
    }

    void save(ostream& out) const {
        out << "EXPFY_NATIVE_INVERSE_V1\n" << setprecision(17);
        out << intercept << ' ' << cutoff << '\n';
        for (int j = 0; j < VARIABLES; ++j)
            out << lower[j] << ' ' << upper[j] << '\n';
        for (double w : weight) out << w << '\n';
    }
};


struct InverseState {
    bool prepared=false,enabled=false;
    InverseConfig config;
    InverseSmoothModel model;
    long long attempts=0,candidates=0,model_pass=0,model_reject=0;
    long long degree_reject=0,histogram_reject=0,structure_reject=0;
    long long committed=0;
    ofstream accepted;
};

void inverse_print_point(ostream& out,int mode,const char* name,
                         const ThresholdPoint73& p,long long candidates) {
    out<<'['<<mode<<"] "<<name<<": cutoff="<<p.cutoff
       <<" candidate_flagged="<<p.path_calls<<'/'<<candidates
       <<" candidate_gate_rate="<<(candidates?static_cast<double>(p.path_calls)/candidates:0.0)
       <<" source_precision="<<p.precision()<<" positive_flagged="<<p.star_calls
       <<" all_flagged="<<p.calls()<<'\n';
}
void inverse_dump_roles(ostream& out,const InverseSpec& s,bool source_space) {
    out<<"mode="<<s.mode<<" forward=T"<<s.forward<<" variant="<<(source_space?'c':'b')<<'\n';
    out<<"Printed degree order is ALWAYS A B C D; degrees are never sorted.\n";
    if (s.forward >= 1 && s.forward <= 9)
        out<<"Original code roles: A=u, B=w, C=v, D=z.\n";
    else if (s.forward >= 10 && s.forward <= 17)
        out<<"Original code roles: A=u, B=v, C=z, D=w.\n";
    else if (s.forward >= 18 && s.forward <= 20)
        out<<"Original code roles: A=u, B=v, C=w, D=z.\n";
    else if (s.forward >= 21 && s.forward <= 23)
        out<<"Original code roles: A=w, B=v, C=z, D=u.\n";
    auto edges=[&](const char* label,unsigned mask) {
        out<<label;
        for(int i=0;i<4;++i)for(int j=i+1;j<4;++j)
            if(mask&inverse_edge_bit(i,j))out<<' '<<char('A'+i)<<char('A'+j);
        out<<'\n';
    };
    edges("forward_source_edges:",s.source_mask);
    edges("forward_output_edges:",s.target_mask);
    out<<"forward_degree_delta:";for(int d:s.forward_delta)out<<' '<<d;out<<'\n';
    out<<"inverse_degree_delta:";for(int d:s.forward_delta)out<<' '<<-d;out<<'\n';
    if(source_space)
        out<<"positive=existing forward-source graphlet (raw degrees)\n"
              "negative=existing forward-output graphlet after hypothetical inverse degree changes\n";
    else
        out<<"positive=existing forward-source graphlet after hypothetical forward degree changes\n"
              "negative=existing forward-output graphlet (raw degrees)\n";
    out<<"These are source/provenance labels, NOT edit-success labels.\n";
}

class Sampler {
public:
    Sampler(Graph& graph,
            vector<int> target_degree,
            double degree_log_bound,
            int histogram_bound,
            vector<uint8_t> good_mask,
            uint64_t seed,
            vector<uint8_t> clique_sources,
            vector<long long> triangle_counts)
        : graph_(graph),
          target_degree_(move(target_degree)),
          degree_log_bound_(degree_log_bound),
          histogram_bound_(histogram_bound),
          good_mask_(move(good_mask)),
          clique_sources_(move(clique_sources)),
          triangle_counts_(move(triangle_counts)),
          rng_(seed),
          current_histogram_(graph.n, 0),
          target_histogram_(graph.n, 0),
          target_log_degree_(graph.n, 0.0),
          log_degree_(graph.n, -numeric_limits<double>::infinity()) {
        seed73_ = seed;
        if (graph_.n < 4) {
            throw runtime_error("graph must have at least 4 nodes");
        }
        if (target_degree_.size() != static_cast<size_t>(graph_.n) ||
            good_mask_.size() != static_cast<size_t>(graph_.n) ||
            clique_sources_.size() != static_cast<size_t>(graph_.n) ||
            triangle_counts_.size() != static_cast<size_t>(graph_.n)) {
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
        if ((mode1_ == 73 && requested1_ > 0) ||
            (mode2_ == 73 && requested2_ > 0)) prepare_model73();
        if (is_new_inverse(mode1_) && requested1_ > 0) prepare_inverse(mode1_);
        if (is_new_inverse(mode2_) && requested2_ > 0) prepare_inverse(mode2_);
        // Build a cheap degree index only when mode 69 is requested. This is
        // a sampling snapshot: each sampled node's current degree is checked
        // again, including after edits made by another mode in a mixed run.
        if ((mode1_ == 69 && requested1_ > 0) ||
            (mode2_ == 69 && requested2_ > 0)) {
            for (auto& bucket : path_leaf_targets_) bucket.clear();
            for (int node = 0; node < graph_.n; ++node) {
                const int degree = graph_.degree[node];
                if (degree < PATH_LEAF_ROOT_DEGREE_LIMIT) {
                    path_leaf_targets_[degree].push_back(node);
                }
            }
        }

        if ((mode1_ == 70 && requested1_ > 0) ||
            (mode2_ == 70 && requested2_ > 0)) {
            // Mode 70 needs a true endpoint leaf. Index these once rather
            // than repeatedly sampling four-sets with no suitable endpoint.
            // Degrees stay fixed under this mode; mixed runs recheck samples.
            box_switch_leaves_.clear();
            for (int node = 0; node < graph_.n; ++node) {
                if (graph_.degree[node] == 1) box_switch_leaves_.push_back(node);
            }
        }

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
               (slot_pending(0) || slot_pending(1));
             ++attempts) {
            choose_active_slot();
            const int mode = active_slot_ == 1 ? mode1_ : mode2_;
            sample_one(mode);
        }

        if (done1_ < requested1_ || done2_ < requested2_) {
            cerr << (attempts >= attempt_limit
                      ? "Warning: sample-attempt limit reached after "
                      : "Warning: a requested mode was disabled after ") << attempts
                 << " attempts before all requested changes were completed.\n";
        }

        report73();
        report_inverse_modes();
        return {done1_, done2_, attempts};
    }

private:
    // Independent learner/cutoff per requested inverse; mode 73 remains intact.
    map<int,InverseState> inverse_states_;
    static constexpr bool INVERSE_SOURCE_SPACE = false;

    bool apply_inverse_operations(const InverseOperations& o,bool commit,
                                  OperationFailure73* failure=nullptr) {
        switch(o.count) {
            case 1:return try_operation({o.data[0]},commit,failure);
            case 2:return try_operation({o.data[0],o.data[1]},commit,failure);
            case 3:return try_operation({o.data[0],o.data[1],o.data[2]},commit,failure);
            case 4:return try_operation({o.data[0],o.data[1],o.data[2],o.data[3]},commit,failure);
            default:throw runtime_error("invalid inverse operation count");
        }
    }

    void prepare_inverse(int mode) {
        auto& state=inverse_states_[mode];if(state.prepared)return;state.prepared=true;
        const auto& spec=inverse_spec(mode);state.config=InverseConfig::load(spec);
        const auto& config=state.config;
        const auto started=chrono::steady_clock::now();
        const auto old_precision=cerr.precision();cerr<<setprecision(10);
        cerr<<'['<<mode<<"] inverse_of=T"<<spec.forward<<" candidate_G"<<spec.target_gid
            <<" -> G"<<spec.source_gid<<" comparison="<<(INVERSE_SOURCE_SPACE?"c_source_space":"b_output_space")
            <<" order=A,B,C,D (never sorted)\n";
        cerr<<'['<<mode<<"] collecting frozen data: attempts="<<config.data_attempts
            <<" min_source_precision="<<config.precision<<" graphlet_effect_checks=removed"
            <<" role_policy=one_random_compatible_assignment\n";
        if(spec.source_gid==spec.target_gid)
            cerr<<'['<<mode<<"] same-type forward transform: each matching sample contributes "
                           "one row to EACH source, not an if/else.\n";
        if(mode==75)
            cerr<<"[75] T4 preserves all four degrees. Roles are essential; b and c use "
                  "identical numerical training data here and should agree for the same seed.\n";
        const FastRng saved=rng_;rng_=FastRng(seed73_^0x243f6a8885a308d3ULL);
        vector<InverseExample> examples;
        examples.reserve(static_cast<size_t>(min(2*config.data_attempts,400000LL)));
        auto record=[&](const array<int,4>& nodes,bool positive) {
            InverseExample e;e.nodes=nodes;e.raw=inverse_degrees(graph_,nodes);
            e.fit.positive=positive;
            e.source_space=positive?e.raw:inverse_shift(e.raw,spec,-1);
            e.fit.degree=INVERSE_SOURCE_SPACE?e.source_space:inverse_shift(e.source_space,spec,+1);
            examples.push_back(e);
        };
        for(long long i=0;i<config.data_attempts;++i) {
            array<int,5> nodes{};nodes[0]=pick_first_node();if(nodes[0]<0)continue;
            int n=1;for(;n<4;++n){int v=pick_neighbor(nodes,n);if(v<0)break;nodes[n]=v;}
            if(n<4)continue;
            shuffle_first_four(nodes);
            const unsigned mask=inverse_mask(build_sub4(nodes));
            if(spec.source_matches[mask].valid)
                record(inverse_role_nodes(nodes,spec.source_matches[mask]),true);
            if(spec.target_matches[mask].valid)
                record(inverse_role_nodes(nodes,spec.target_matches[mask]),false);
        }
        rng_=saved;
        const size_t groups=inverse_assign_partitions(examples,spec,seed73_);
        vector<InverseFitExample> fit;fit.reserve(examples.size());
        array<array<long long,2>,3> counts{};
        for(const auto& e:examples){fit.push_back(e.fit);++counts[e.fit.partition][e.fit.positive?1:0];}
        for(int i=0;i<3;++i) {
            static const char* names[]={"fit","calibration","heldout"};
            cerr<<'['<<mode<<"] "<<names[i]<<": candidates="<<counts[i][0]<<" positives="<<counts[i][1]<<'\n';
        }
        cerr<<'['<<mode<<"] groups="<<groups<<"; grouping uses common source-space degrees "
                         "and sampled four-sets. No refit follows cutoff selection.\n";
        const bool dump=!config.dump_directory.empty();
        std::filesystem::path dir(config.dump_directory);
        if(dump) {
            std::filesystem::create_directories(dir);
            ofstream real(dir/"existing.txt"),modified(dir/"modified.txt"),samples(dir/"samples.csv"),roles(dir/"roles.txt");
            if(!real||!modified||!samples||!roles)throw runtime_error("cannot create inverse data dump");
            inverse_dump_roles(roles,spec,INVERSE_SOURCE_SPACE);
            samples<<"sample_id,A,B,C,D,raw_dA,raw_dB,raw_dC,raw_dD,feature_dA,feature_dB,feature_dC,feature_dD,positive,partition,group\n";
            for(size_t i=0;i<examples.size();++i) {
                const auto& e=examples[i];
                const bool existing=INVERSE_SOURCE_SPACE?e.fit.positive:!e.fit.positive;
                ostream& file=existing?static_cast<ostream&>(real):static_cast<ostream&>(modified);
                file<<e.fit.degree[0]<<' '<<e.fit.degree[1]<<' '<<e.fit.degree[2]<<' '<<e.fit.degree[3]<<'\n';
                samples<<i;for(int v:e.nodes)samples<<','<<v;for(int d:e.raw)samples<<','<<d;
                for(int d:e.fit.degree)samples<<','<<d;
                samples<<','<<e.fit.positive<<','<<e.fit.partition<<','<<e.group<<'\n';
            }
            real.close();modified.close();samples.close();roles.close();
            if(!real||!modified||!samples||!roles)throw runtime_error("cannot finish inverse data dump");
        }
        auto disable=[&](const string& why) {
            cerr<<'['<<mode<<"] DISABLED: "<<why<<'\n';
            if(dump){ofstream r(dir/"threshold_report.txt");r<<"DISABLED: "<<why<<'\n';}
            cerr.precision(old_precision);
        };
        if(counts[1][0]<20 || counts[1][1]<20) {
            disable("need >=20 calibration observations per source; inspect counts or increase DATA_ATTEMPTS");return;
        }
        if(!state.model.train(fit,config.C,config.max_iterations)) {
            disable("need >=20 fit observations per source, >=10 distinct descriptors, and finite fit");return;
        }
        const auto trained=chrono::steady_clock::now();
        cerr<<'['<<mode<<"] fit: converged="<<state.model.converged<<" iterations="<<state.model.iterations
            <<" gradient_inf="<<state.model.final_gradient<<'\n';
        if(!state.model.converged)cerr<<'['<<mode<<"] WARNING: optimizer tolerance not reached.\n";
        vector<ScoredSource73> calibration,test;
        ofstream scored;if(dump){scored.open(dir/"scored_samples.csv");scored<<"sample_id,positive,partition,gate_score\n"<<setprecision(17);}
        for(size_t i=0;i<examples.size();++i) {
            const auto& e=examples[i];if(e.fit.partition==0)continue;
            const double score=state.model.score(e.fit.degree);
            if(!isfinite(score))throw runtime_error("nonfinite learned inverse score");
            (e.fit.partition==1?calibration:test).push_back({score,e.fit.positive});
            if(dump)scored<<i<<','<<e.fit.positive<<','<<e.fit.partition<<','<<score<<'\n';
        }
        if(dump){scored.close();if(!scored)throw runtime_error("cannot finish score dump");}
        const auto choice=choose_max_gate_cutoff73(calibration,config.precision);
        state.enabled=choice.feasible;
        state.model.cutoff=state.enabled?choice.chosen.cutoff:numeric_limits<double>::infinity();
        auto report=[&](ostream& out) {
            out<<'['<<mode<<"] cutoff_policy=max_candidate_gate_at_precision target="<<config.precision<<'\n';
            if(state.enabled)inverse_print_point(out,mode,"CHOSEN calibration",choice.chosen,choice.path_total);
            else {
                out<<'['<<mode<<"] DISABLED: no threshold accepts a candidate at the requested source precision. No relaxation.\n";
                if(choice.has_best_nonzero_path)
                    inverse_print_point(out,mode,"best_precision_with_candidates (NOT used)",choice.best_nonzero_path,choice.path_total);
            }
            if(state.enabled && choice.chosen.calls()<30)
                out<<'['<<mode<<"] WARNING: <30 flagged calibration rows; duplicates can reduce support further.\n";
        };
        report(cerr);
        ThresholdPoint73 held;held.cutoff=state.model.cutoff;
        long long negatives=0,positives=0;
        for(const auto& e:test) {
            if(e.from_star)++positives;else ++negatives;
            if(state.enabled && e.score>=state.model.cutoff){if(e.from_star)++held.star_calls;else ++held.path_calls;}
        }
        inverse_print_point(cerr,mode,"HELDOUT (not used for selection)",held,negatives);
        if(state.enabled && (held.calls()==0 || held.precision()<config.precision))
            cerr<<'['<<mode<<"] WARNING: held-out source precision unavailable/below target; cutoff unchanged.\n";
        sort(test.begin(),test.end(),[](const auto& a,const auto& b){return a.score<b.score;});
        long double pairs=0,prev_negative=0;
        for(size_t i=0;i<test.size();) {
            size_t j=i;long long p=0,n=0;
            while(j<test.size()&&test[j].score==test[i].score){if(test[j].from_star)++p;else ++n;++j;}
            pairs+=p*(prev_negative+0.5L*n);prev_negative+=n;i=j;
        }
        const double auc=positives&&negatives?static_cast<double>(pairs/(static_cast<long double>(positives)*negatives)):numeric_limits<double>::quiet_NaN();
        cerr<<'['<<mode<<"] heldout_source_AUC="<<auc<<"; source precision is NOT edit-success precision.\n";
        if(dump) {
            ofstream r(dir/"threshold_report.txt"),curve(dir/"thresholds.csv"),model(dir/"model.txt");
            if(!r||!curve||!model)throw runtime_error("cannot create inverse reports");
            r<<setprecision(17);report(r);inverse_print_point(r,mode,"HELDOUT",held,negatives);
            curve<<"cutoff,candidates_flagged,positives_flagged,candidate_gate_rate,source_precision,meets_precision,chosen\n"<<setprecision(17);
            for(const auto& p:choice.curve)curve<<p.cutoff<<','<<p.path_calls<<','<<p.star_calls<<','
                <<static_cast<double>(p.path_calls)/choice.path_total<<','<<p.precision()<<','
                <<(p.precision()>=config.precision)<<','<<(state.enabled&&p.cutoff==choice.chosen.cutoff)<<'\n';
            inverse_dump_roles(model,spec,INVERSE_SOURCE_SPACE);state.model.save(model);
            r.close();curve.close();model.close();if(!r||!curve||!model)throw runtime_error("cannot finish inverse reports");
            state.accepted.open(dir/"accepted.csv");
            if(!state.accepted)throw runtime_error("cannot create inverse accepted log");
            state.accepted<<"edit,A,B,C,D,pre_dA,pre_dB,pre_dC,pre_dD,feature_dA,feature_dB,feature_dC,feature_dD,post_dA,post_dB,post_dC,post_dD,score,cutoff,exact_known,delta_G5,delta_G6,delta_G7,delta_G8,delta_G9,delta_G10\n"<<setprecision(17);
        }
        cerr<<'['<<mode<<"] seconds: collect_and_fit="<<chrono::duration<double>(trained-started).count()
            <<" calibration="<<chrono::duration<double>(chrono::steady_clock::now()-trained).count()<<'\n';
        cerr.precision(old_precision);
    }

    void apply_learned_inverse(int mode,const array<int,5>& nodes,const Sub4& sub) {
        auto& state=inverse_states_.at(mode);if(!state.enabled)return;
        const auto& spec=inverse_spec(mode);const auto& match=spec.target_matches[inverse_mask(sub)];
        if(!match.valid){++state.structure_reject;return;}
        ++state.candidates;
        const auto roles=inverse_role_nodes(nodes,match);
        const auto pre=inverse_degrees(graph_,roles),post=inverse_shift(pre,spec,-1);
        const auto feature=INVERSE_SOURCE_SPACE?post:pre;
        array<double,4> logs{};
        for(int j=0;j<4;++j) {
            if(feature[j]<=0 || feature[j]>=graph_.n){++state.structure_reject;return;}
            logs[j]=log_degree_[feature[j]];
        }
        const double score=state.model.score_logs(logs);
        if(!isfinite(score)||score<state.model.cutoff){++state.model_reject;return;}
        ++state.model_pass;
        const auto ops=inverse_operations(spec,roles);OperationFailure73 failure;
        if(!apply_inverse_operations(ops,true,&failure)) {
            if(failure==OperationFailure73::Degree)++state.degree_reject;
            else if(failure==OperationFailure73::Histogram)++state.histogram_reject;
            else ++state.structure_reject;
            return;
        }
        if(inverse_degrees(graph_,roles)!=post)throw runtime_error("inverse post-degree mapping mismatch");
        ++state.committed;
        if(state.accepted.is_open()) {
            state.accepted<<state.committed;for(int v:roles)state.accepted<<','<<v;
            for(int d:pre)state.accepted<<','<<d;
            for(int d:feature)state.accepted<<','<<d;
            for(int d:post)state.accepted<<','<<d;
            // Keep legacy CSV columns: exact_known=0; deltas are unmeasured.
            state.accepted<<','<<score<<','<<state.model.cutoff<<",0,,,,,,\n";
            if(!state.accepted)throw runtime_error("cannot write inverse accepted edit");
        }
    }

    void report_inverse_modes() {
        for(auto& entry:inverse_states_) {
            const int mode=entry.first;auto& s=entry.second;
            cerr<<'['<<mode<<"] summary: sampling_attempts="<<s.attempts<<" candidates="<<s.candidates
                <<" model_pass="<<s.model_pass<<" model_reject="<<s.model_reject
                <<" degree_reject="<<s.degree_reject<<" histogram_reject="<<s.histogram_reject
                <<" structure_reject="<<s.structure_reject
                <<" committed="<<s.committed<<'\n';
            if(s.accepted.is_open()){s.accepted.flush();if(!s.accepted)throw runtime_error("cannot flush inverse log");}
        }
    }

    Graph& graph_;
    vector<int> target_degree_;
    double degree_log_bound_ = 0.0;
    double degree_absolute_fallback_ = 5.0;
    int histogram_bound_ = 0;
    vector<uint8_t> good_mask_;
    vector<uint8_t> clique_sources_;
    // ORCA orbit 3 is a proposal heuristic, like good_mask_. Leaf relocations
    // preserve these counts exactly when the input table is current. Other
    // modes in a mixed run, or a frozen table, can make the estimates stale.
    vector<long long> triangle_counts_;
    vector<int> good_nodes_;
    vector<int> ordinary_nodes_;
    // Snapshot buckets avoid repeatedly drawing a node of the wrong degree.
    array<vector<int>, PATH_LEAF_ROOT_DEGREE_LIMIT> path_leaf_targets_;
    vector<int> box_switch_leaves_;
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

    // Caches are allocated lazily only when modes 61--64 need them. Every
    // successful edit (including legacy modes in a mixed pair) invalidates
    // cached decisions by advancing graph_revision_. Failed samples do not.
    uint64_t graph_revision_ = 1;
    struct RootStatistics {
        uint64_t revision = 0;
        long long triangles = 0;
        long long open_two_paths = 0;
        bool available = false;
    };
    vector<RootStatistics> root_cache_;
    vector<array<int, 3>> pendant_triangles_;
    uint64_t pendant_triangles_revision_ = 0;

    uint64_t seed73_ = 0;
    bool prepared73_ = false;
    bool enabled73_ = false;
    Config73 config73_;
    SmoothSourceModel73 model73_;
    Counters73 counters73_;
    ofstream accepted_log73_;
    double tuple_score73(const DegreeTuple73& t) const {
        // Inputs are already paths. Do not apply a path-to-star degree shift.
        for (int d : t)
            if (d <= 0 || d >= graph_.n) return -numeric_limits<double>::infinity();
        return model73_.score_logs({log_degree_[t[0]], log_degree_[t[1]],
                                    log_degree_[t[2]], log_degree_[t[3]]});
    }

    double path_gate_score73(const DegreeTuple73& t) const {
        const double first = tuple_score73(t);
        if (!config73_.both_orientations) return first;
        return max(first, tuple_score73({t[3],t[2],t[1],t[0]}));
    }

    void prepare_model73() {
        if (prepared73_) return;
        prepared73_ = true;
        config73_ = Config73::from_environment();
        const auto started = chrono::steady_clock::now();
        const streamsize previous_precision = cerr.precision();
        cerr << setprecision(8);
        cerr << "[73] collecting frozen data: attempts=" << config73_.data_attempts
             << " optimize_path_gate_rate=1"
             << " min_source_precision=" << config73_.min_source_precision
             << " both=" << config73_.both_orientations
             << " graphlet_effect_checks=removed\n";
        cerr << "[73] representation=line_space order=A-B-C-D (not sorted)\n"
             << "[73] labels: existing_P4=0, hypothetical_star_to_P4=1. "
             << "These are provenance labels, NOT edit-improvement labels.\n";

        // Use the same biased four-set sampler as expfy, with an independent
        // RNG stream, and leave both the graph and operation RNG unchanged.
        const FastRng saved_rng = rng_;
        rng_ = FastRng(seed73_ ^ 0x243f6a8885a308d3ULL);
        vector<Example73> examples;
        examples.reserve(static_cast<size_t>(min(config73_.data_attempts, 200000LL)));
        array<array<long long,2>,3> sizes{};
        auto record = [&](const DegreeTuple73& t, bool from_star) {
            const int partition = grouped_partition73(t, seed73_);
            examples.push_back({t,from_star,partition});
            sizes[partition][from_star ? 1 : 0]++;
        };
        for (long long attempt = 0; attempt < config73_.data_attempts; ++attempt) {
            ++counters73_.training_attempts;
            array<int,5> nodes{};
            nodes[0] = pick_first_node();
            if (nodes[0] < 0) continue;
            int count = 1;
            for (; count < 4; ++count) {
                const int next = pick_neighbor(nodes,count);
                if (next < 0) break;
                nodes[count] = next;
            }
            if (count != 4) continue;
            shuffle_first_four(nodes);
            ++counters73_.complete_training_samples;
            const Sub4 sub = build_sub4(nodes);
            const int gid = graph_id_from_sub4(sub);
            if (gid == 5) {
                array<int,3> leaves{};
                int center = -1, found = 0;
                for (int i = 0; i < 4; ++i) {
                    if (sub_degree(sub,i) == 3) center = i;
                    else if (sub_degree(sub,i) == 1) leaves[found++] = i;
                }
                if (center < 0 || found != 3) continue;
                // Star C-A,C-B,C-D -> path A-B-C-D after removing C-A,
                // adding A-B. Exactly the degree convention of expfydc.
                record(path_from_star73({graph_.degree[nodes[leaves[0]]],
                        graph_.degree[nodes[leaves[1]]],
                        graph_.degree[nodes[center]],
                        graph_.degree[nodes[leaves[2]]]}), true);
            } else if (gid == 6) {
                array<int,2> leaves{};
                int found = 0;
                for (int i = 0; i < 4; ++i)
                    if (sub_degree(sub,i) == 1) leaves[found++] = i;
                if (found != 2) continue;
                int b = -1, c = -1;
                for (int i = 0; i < 4; ++i) {
                    if (sub[i][leaves[0]]) b = i;
                    if (sub[i][leaves[1]]) c = i;
                }
                if (b < 0 || c < 0) continue;
                record({graph_.degree[nodes[leaves[0]]], graph_.degree[nodes[b]],
                        graph_.degree[nodes[c]], graph_.degree[nodes[leaves[1]]]}, false);
            }
        }
        rng_ = saved_rng;
        const auto collected = chrono::steady_clock::now();
        for (int p = 0; p < 3; ++p) {
            static const char* names[] = {"fit", "calibration", "heldout"};
            cerr << "[73] " << names[p] << ": existing_paths=" << sizes[p][0]
                 << " star_derived=" << sizes[p][1] << '\n';
        }
        cerr << "[73] all copies of a degree tuple, including its reverse, "
                "share one partition. No all-data refit follows calibration.\n";

        if (!config73_.dump_directory.empty()) {
            namespace fs = std::filesystem;
            fs::create_directories(config73_.dump_directory);
            const fs::path dir(config73_.dump_directory);
            ofstream existing(dir/"defdeg"), derived(dir/"moddeg"), all(dir/"samples.csv");
            if (!existing || !derived || !all) throw runtime_error("cannot write mode73 data dump");
            all << "w,x,y,z,from_star,partition\n";
            for (const auto& e : examples) {
                ostream& out = e.from_star ? static_cast<ostream&>(derived) : static_cast<ostream&>(existing);
                out << e.degree[0] << ' ' << e.degree[1] << ' ' << e.degree[2] << ' ' << e.degree[3] << '\n';
                all << e.degree[0] << ',' << e.degree[1] << ',' << e.degree[2] << ','
                    << e.degree[3] << ',' << e.from_star << ',' << e.partition << '\n';
            }
            existing.close(); derived.close(); all.close();
            if (!existing || !derived || !all) throw runtime_error("cannot finish mode73 data dump");
            accepted_log73_.open(dir/"accepted.csv");
            if (!accepted_log73_) throw runtime_error("cannot write mode73 edit log");
            accepted_log73_ << "edit,a,b,c,d,degree_a,degree_b,degree_c,degree_d,source_score,cutoff,exact_known,delta_G5,delta_G6,delta_G7,delta_G8,delta_G9,delta_G10\n";
            accepted_log73_ << setprecision(17);
        }
        if (sizes[1][0] < 20 || sizes[1][1] < 20 || !model73_.train(examples, config73_.C, config73_.max_iterations)) {
            cerr << "[73] DISABLED: need >=20 fitting examples per source, >=10 distinct fitting "
                    "tuples and >=20 calibration rows per source. Increase DATA_ATTEMPTS "
                    "or inspect whether both graphlets exist.\n";
            cerr.precision(previous_precision);
            return;
        }
        const auto trained = chrono::steady_clock::now();
        cerr << "[73] fit: iterations=" << model73_.iterations
             << " converged=" << model73_.converged
             << " gradient_inf=" << model73_.final_gradient
             << " objective=" << model73_.final_objective << '\n';
        if (!model73_.converged)
            cerr << "[73] WARNING: optimizer tolerance not reached; using last finite iterate.\n";

        // Score both sources using the SAME rule that defines a path gate.
        // With BOTH=1 this is max(score(t),score(reverse(t))) for BOTH classes.
        // Each row still counts once, not once per successful orientation.
        vector<ScoredSource73> calibration_scores;
        vector<ScoredSource73> source_test;
        ofstream score_dump;
        if (!config73_.dump_directory.empty()) {
            score_dump.open(std::filesystem::path(config73_.dump_directory)/"scored_samples.csv");
            if (!score_dump) throw runtime_error("cannot write mode73 score dump");
            score_dump << "w,x,y,z,from_star,partition,gate_score\n" << setprecision(17);
        }
        for (const auto& e : examples) {
            if (e.partition == 0) continue;
            const double score = path_gate_score73(e.degree);
            if (!isfinite(score)) throw runtime_error("nonfinite mode73 held-out score");
            if (e.partition == 1) calibration_scores.push_back({score,e.from_star});
            else source_test.push_back({score,e.from_star});
            if (score_dump.is_open()) {
                for (int value : e.degree) score_dump << value << ',';
                score_dump << e.from_star << ',' << e.partition << ',' << score << '\n';
            }
        }
        if (score_dump.is_open()) {
            score_dump.close();
            if (!score_dump) throw runtime_error("cannot finish mode73 score dump");
        }
        const auto choice = choose_max_gate_cutoff73(
            calibration_scores, config73_.min_source_precision);
        model73_.cutoff = choice.feasible ? choice.chosen.cutoff : numeric_limits<double>::infinity();
        enabled73_ = choice.feasible;
        print_max_gate_choice73(cerr,choice,config73_.min_source_precision);

        // Test labels are for reporting only. Never search or change the
        // chosen cutoff on this partition, and never refit the model afterward.
        ThresholdPoint73 test_point;
        test_point.cutoff = model73_.cutoff;
        long long test_existing = 0, test_derived = 0;
        for (const auto& row : source_test) {
            if (row.from_star) ++test_derived; else ++test_existing;
            if (enabled73_ && row.score >= model73_.cutoff) {
                if (row.from_star) ++test_point.star_calls;
                else ++test_point.path_calls;
            }
        }
        print_threshold_point73(cerr,"HELDOUT (not used for selection)",test_point,test_existing);
        const bool heldout_precision_ok = test_point.calls() &&
            test_point.precision() >= config73_.min_source_precision;
        cerr << "[73] heldout_precision_target_met=" << heldout_precision_ok
             << " heldout_moddeg=" << test_derived << '\n';
        if (enabled73_ && !heldout_precision_ok)
            cerr << "[73] WARNING: the precision target did not hold on held-out rows "
                    "(or there were no held-out flags). The cutoff is unchanged; "
                    "calibration precision is not a deployment guarantee.\n";
        if (test_existing < 20 || test_derived < 20)
            cerr << "[73] WARNING: fewer than 20 held-out rows in one source; "
                    "held-out metrics are based on very little support.\n";
        sort(source_test.begin(), source_test.end(), [](const auto& a, const auto& b) {
            return a.score < b.score;
        });
        double pairs = 0, previous_negatives = 0;
        long long positives = 0, negatives = 0;
        for (size_t i = 0; i < source_test.size();) {
            size_t j = i;
            long long p = 0, n = 0;
            while (j < source_test.size() && source_test[j].score == source_test[i].score) {
                if (source_test[j].from_star) ++p; else ++n;
                ++j;
            }
            pairs += p*(previous_negatives+0.5*n);
            previous_negatives += n; positives += p; negatives += n;
            i = j;
        }
        const double auc = positives && negatives ? pairs/(static_cast<double>(positives)*negatives)
                                                   : numeric_limits<double>::quiet_NaN();
        cerr << "[73] heldout_source_AUC_gate_score=" << auc
             << " gate=" << (config73_.both_orientations ? "max_of_two_orientations" : "single_orientation")
             << '\n';
        cerr << "[73] Source precision means moddeg among flagged defdeg+moddeg rows. "
                "It does NOT measure successful graphlet-improving edits.\n";
        if (isfinite(auc) && auc < 0.55)
            cerr << "[73] WARNING: weak held-out source discrimination; validate actual edit effects.\n";
        if (!config73_.dump_directory.empty()) {
            const auto dir = std::filesystem::path(config73_.dump_directory);
            ofstream report(dir/"threshold_report.txt"), curve(dir/"thresholds.csv");
            if (!report || !curve) throw runtime_error("cannot write mode73 threshold diagnostics");
            report << setprecision(17);
            print_max_gate_choice73(report,choice,config73_.min_source_precision);
            print_threshold_point73(report,"HELDOUT (not used for selection)",test_point,test_existing);
            curve << "cutoff,defdeg_flagged,moddeg_flagged,path_gate_rate,source_precision,accepts_path,meets_min_precision,chosen\n";
            curve << setprecision(17);
            for (const auto& p : choice.curve) {
                curve << p.cutoff << ',' << p.path_calls << ',' << p.star_calls << ','
                      << static_cast<double>(p.path_calls)/choice.path_total << ','
                      << p.precision() << ','
                      << (p.path_calls > 0) << ','
                      << (p.precision() >= config73_.min_source_precision) << ','
                      << (choice.feasible && p.cutoff == choice.chosen.cutoff) << '\n';
            }
            report.close(); curve.close();
            if (!report || !curve) throw runtime_error("cannot finish mode73 threshold diagnostics");
        }
        if (!config73_.dump_directory.empty()) {
            ofstream model_file(std::filesystem::path(config73_.dump_directory)/"model.txt");
            if (!model_file) throw runtime_error("cannot write mode73 model dump");
            model73_.save(model_file);
            model_file.close();
            if (!model_file) throw runtime_error("cannot finish mode73 model dump");
        }
        const auto finished = chrono::steady_clock::now();
        cerr << "[73] seconds: collect=" << chrono::duration<double>(collected-started).count()
             << " fit_and_dump=" << chrono::duration<double>(trained-collected).count()
             << " calibration_diagnostics=" << chrono::duration<double>(finished-trained).count() << '\n';
        cerr.precision(previous_precision);
    }

    bool try_model_relocation73(int a, int b, int c, int d) {
        const DegreeTuple73 degrees = {graph_.degree[a], graph_.degree[b],
                                        graph_.degree[c], graph_.degree[d]};
        ++counters73_.orientations_scored;
        const double score = tuple_score73(degrees);
        if (!isfinite(score) || score < model73_.cutoff) {
            ++counters73_.model_reject;
            return false;
        }
        ++counters73_.model_pass;
        const initializer_list<EdgeOp> operation = {
            {EdgeAction::Add, a, c}, {EdgeAction::Remove, a, b}
        };
        OperationFailure73 failure;
        if (!try_operation(operation, true, &failure)) {
            if (failure == OperationFailure73::Degree) ++counters73_.degree_reject;
            else if (failure == OperationFailure73::Histogram) ++counters73_.histogram_reject;
            else ++counters73_.structure_reject;
            return false;
        }
        ++counters73_.committed;
        if (accepted_log73_.is_open()) {
            accepted_log73_ << counters73_.committed << ',' << a << ',' << b << ',' << c << ',' << d;
            for (int value : degrees) accepted_log73_ << ',' << value;
            // Keep legacy CSV columns: exact_known=0; deltas are unmeasured.
            accepted_log73_ << ',' << score << ',' << model73_.cutoff << ",0,,,,,,\n";
            if (!accepted_log73_) throw runtime_error("cannot write mode73 accepted edit log");
        }
        return true;
    }

    void report73() {
        if (!prepared73_) return;
        const auto& c = counters73_;
        cerr << "[73] summary: sampling_attempts=" << c.sampling_attempts
             << " path_candidates=" << c.path_candidates
             << " orientations_scored=" << c.orientations_scored
             << " model_pass=" << c.model_pass << " model_reject=" << c.model_reject
             << " structure_reject=" << c.structure_reject
             << " degree_reject=" << c.degree_reject << " histogram_reject=" << c.histogram_reject
             << " committed=" << c.committed << '\n';
        if (accepted_log73_.is_open()) {
            accepted_log73_.flush();
            if (!accepted_log73_) throw runtime_error("cannot flush mode73 accepted edit log");
        }
    }

    bool slot_pending(int slot) const {
        const int mode = slot == 0 ? mode1_ : mode2_;
        if (mode == 73 && prepared73_ && !enabled73_) return false;
        if (is_new_inverse(mode)) {
            const auto it = inverse_states_.find(mode);
            if (it != inverse_states_.end() && it->second.prepared && !it->second.enabled) return false;
        }
        return slot == 0 ? done1_ < requested1_ : done2_ < requested2_;
    }

    static void validate_mode(int mode) {
        if (mode < 1 || mode > MAX_TRANSFORMATION || mode == 25 || mode == 26) {
            throw runtime_error("transformation mode must be 1..95 except 25 and 26");
        }
    }

    void choose_active_slot() {
        const bool need1 = slot_pending(0);
        const bool need2 = slot_pending(1);

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
        return new_distance <= old_distance;
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
        return new_distance <= old_distance;
    }

    bool try_operation(initializer_list<EdgeOp> operations,
                       bool commit = true, OperationFailure73* failure = nullptr) {
        if (failure) *failure = OperationFailure73::Structure;
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

        if (failure) *failure = OperationFailure73::Degree;
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
        if (failure) *failure = OperationFailure73::Histogram;
        for (int i = 0; i < related_count; ++i) {
            if (!histogram_allowed(
                    related_degree[i], proposed_histogram[i])) {
                return false;
            }
        }

        if (failure) *failure = OperationFailure73::None;
        if (!commit) return true;
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

        ++graph_revision_;
        if (active_slot_ == 1 && done1_ < requested1_) {
            done1_++;
        } else if (active_slot_ == 2 && done2_ < requested2_) {
            done2_++;
        }
        return true;
    }

    static long long choose_two(long long count) {
        return count < 2 ? 0 : count * (count - 1) / 2;
    }

    // Triangles containing root, and induced length-two paths root-b-c.
    // There is no approximation or sampled count in this helper.
    bool rooted_statistics(int root, RootStatistics& result) {
        if (root_cache_.empty()) root_cache_.resize(graph_.n);
        RootStatistics& cached = root_cache_[root];
        if (cached.revision == graph_revision_) {
            result = cached;
            return cached.available;
        }
        cached = {};
        cached.revision = graph_revision_;
        const auto& neighbors = graph_.adj[root];
        const long long degree = static_cast<long long>(neighbors.size());
        if (degree + choose_two(degree) > G67_WORK_LIMIT) {
            result = cached;
            return false;
        }
        long long triangles = 0;
        long long walks = 0;
        for (size_t i = 0; i < neighbors.size(); ++i) {
            walks += graph_.degree[neighbors[i]] - 1LL;
            for (size_t j = i + 1; j < neighbors.size(); ++j) {
                triangles += graph_.matrix.get(neighbors[i], neighbors[j]);
            }
        }
        cached.triangles = triangles;
        cached.open_two_paths = walks - 2 * triangles;
        cached.available = true;
        result = cached;
        return true;
    }

    // Reattach a true degree-one vertex. All changed four-node graphlets
    // contain that leaf, so only G5 (claw), G6 (P4), and G7 (paw) can change.
    // The deltas below account for EVERY affected four-set, not just the
    // graphlet from which the proposal was sampled.
    bool leaf_relocation_allowed(int leaf, int from, int to, int mode) {
        if ((mode != 61 && mode != 62) || leaf < 0 || from < 0 || to < 0 ||
            leaf >= graph_.n || from >= graph_.n || to >= graph_.n ||
            leaf == from || leaf == to || from == to ||
            graph_.degree[leaf] != 1 || !graph_.matrix.get(leaf, from) ||
            graph_.matrix.get(leaf, to)) {
            return false;
        }
        RootStatistics old_root, new_root;
        if (!rooted_statistics(from, old_root) ||
            !rooted_statistics(to, new_root)) return false;

        // Work in H = G minus the leaf. Removing the leaf does not change
        // triangles. It removes one open root-from-leaf path at 'to' exactly
        // when from and to are adjacent.
        const long long delta7 = new_root.triangles - old_root.triangles;
        const long long delta6 = new_root.open_two_paths -
            static_cast<long long>(graph_.matrix.get(from, to)) -
            old_root.open_two_paths;
        const long long delta5 =
            choose_two(graph_.degree[to]) - new_root.triangles -
            (choose_two(graph_.degree[from] - 1LL) - old_root.triangles);
        if (delta6 > 0 || delta7 > 0 || delta5 + delta6 + delta7 < 0) {
            return false;
        }
        if ((mode == 61 && delta6 >= 0) || (mode == 62 && delta7 >= 0)) {
            return false;
        }
        // The total number of connected four-sets does not fall, while both
        // G6 and G7 counts do not rise. Their percentages therefore do not
        // rise either. G6 drops strictly in mode 61, G7 in mode 62.
        return true;
    }

    bool try_leaf_relocation(int leaf, int from, int to, int mode) {
        if (!leaf_relocation_allowed(leaf, from, to, mode)) return false;
        return try_operation({
            {EdgeAction::Add, leaf, to},
            {EdgeAction::Remove, leaf, from},
        });
    }

    bool try_triangle_guided_leaf_relocation(int leaf, int from, int to) {
        if (leaf < 0 || from < 0 || to < 0 || leaf >= graph_.n ||
            from >= graph_.n || to >= graph_.n || leaf == to || from == to ||
            graph_.degree[leaf] != 1 || !graph_.matrix.get(leaf, from) ||
            graph_.matrix.get(leaf, to)) return false;

        // Only four-sets containing the leaf change. Its paw count equals
        // the number of triangles at its root; its star count equals the
        // number of nonadjacent pairs among the root's other neighbors.
        // Cached orbit counts make this an O(1) proposal rule. This is not
        // an exact acceptance certificate when another mode changes triangles
        // or the caller deliberately freezes the orbit table.
        const long long paw_delta = triangle_counts_[to] - triangle_counts_[from];
        if (paw_delta <= 0) return false;
        const long long star_delta = choose_two(graph_.degree[to]) -
            choose_two(graph_.degree[from] - 1LL) - paw_delta;
        // Allow some star growth while favoring triangle extensions (paws).
        // Actual log-frequency improvement is measured by the outer optimizer.
        if (static_cast<long double>(star_delta) > 2.0L * paw_delta) return false;

        // A leaf can be reattached anywhere in its connected component without
        // disconnecting the remaining graph. Raw G8/G9/G10 stay unchanged.
        return try_operation({
            {EdgeAction::Remove, leaf, from},
            {EdgeAction::Add, leaf, to},
        });
    }

    bool has_only_common_neighbor(int a, int b, int known_common) {
        // The caller supplies a known common neighbor from its sampled paw.
        // Inspect the smaller adjacency list, or decline if it is too large.
        if (graph_.degree[a] > graph_.degree[b]) swap(a, b);
        if (graph_.degree[a] > PAW_BOX_MAX_LOCAL_DEGREE) return false;
        for (int x : graph_.adj[a]) {
            if (x != known_common && graph_.matrix.get(b, x)) return false;
        }
        return true;
    }

    bool try_triangle_branch_exchange(int u, int v, int z, int w) {
        // Sampled paw: triangle u-v-z and tail z-w. Keep the larger hub u
        // and hub w at their current degrees while moving a nonleaf branch
        // toward u. This favors path growth when paths are underrepresented.
        if (graph_.degree[u] <= graph_.degree[w]) return false;
        const int leaf = random_neighbor(u);
        if (leaf < 0 || graph_.degree[leaf] != 1) return false;
        const int branch = random_neighbor(w);
        if (branch < 0 || branch == u || branch == v || branch == z ||
            branch == leaf || graph_.matrix.get(u, branch) ||
            graph_.degree[branch] < graph_.degree[v] ||
            graph_.degree[branch] > TRIANGLE_BRANCH_MAX_DEGREE) return false;

        // Protect diamonds and cliques: no source triangle edge may belong
        // to another triangle. Each check reads at most 32 neighbor entries.
        if (!has_only_common_neighbor(u, v, z) ||
            !has_only_common_neighbor(u, z, v) ||
            !has_only_common_neighbor(v, z, u)) return false;

        // Neither removing w-branch nor adding u-branch may change a
        // triangle. v is allowed as a neighbor because u-v is being removed.
        // The new w-leaf edge cannot form a triangle since u-w is absent.
        // All checks use current adjacency, even with a frozen orbit table.
        for (int x : graph_.adj[branch]) {
            if (x == w) continue;
            if (graph_.matrix.get(w, x) ||
                (x != v && graph_.matrix.get(u, x))) return false;
        }

        // The selected triangle opens into u-z-w-leaf-u. Only v loses a
        // degree and the old leaf gains one; raw diamonds/cliques stay fixed
        // and raw paws fall. Global path/box effects are measured outside
        // expfy, without enumerating surrounding four-sets here.
        // u-z-v replaces the deleted triangle edge. The surviving u-z-w
        // path keeps both sides connected when branch changes attachment.
        return try_operation({
            {EdgeAction::Remove, u, v},
            {EdgeAction::Remove, w, branch},
            {EdgeAction::Add, w, leaf},
            {EdgeAction::Add, u, branch},
        });
    }

    bool try_paw_box_switch(int u, int v, int z, int w) {
        // Sampled paw: triangle u-v-z and tail z-w. Find a true leaf at w,
        // then exchange u-z and w-leaf for u-w and z-leaf. The triangle
        // becomes the box u-v-z-w-u, and every vertex keeps its degree.
        const int leaf = random_neighbor(w);
        if (leaf < 0 || graph_.degree[leaf] != 1 ||
            graph_.degree[u] > PAW_BOX_MAX_LOCAL_DEGREE) return false;

        // No triangle edge may belong to another triangle. Thus removing
        // u-z cannot destroy a diamond or clique. u and w currently share
        // only z, so deleting u-z also ensures that adding u-w creates no
        // triangle. These checks use current adjacency, not cached orbits.
        if (!has_only_common_neighbor(u, z, v) ||
            !has_only_common_neighbor(u, v, z) ||
            !has_only_common_neighbor(v, z, u) ||
            !has_only_common_neighbor(u, w, z)) return false;

        int work = 0;
        for (int t : graph_.adj[u]) {
            if (t == z) continue;
            if (graph_.degree[t] > PAW_BOX_NEIGHBOR_WORK_LIMIT - work) return false;
            work += graph_.degree[t];
        }

        // With C(a,b) denoting the current common-neighbor count, the exact
        // raw box delta is C(v,w) + sum_{t in N(u)-{z,v}}[C(t,w)-C(t,z)+1].
        // The correction excludes u itself from old common-neighbor counts.
        // Decline moves whose newly formed box is outweighed by losses
        // elsewhere. At most 128 neighbor entries are read in this sum.
        int boxes = 0;
        for (int t : graph_.adj[u]) {
            if (t == z) continue;
            int common_w = 0, common_z = 0;
            for (int x : graph_.adj[t]) {
                common_w += graph_.matrix.get(w, x);
                common_z += graph_.matrix.get(z, x);
            }
            boxes += t == v ? common_w : common_w - common_z + 1;
        }
        if (boxes <= 0) return false;

        // Exactly degree(u)+degree(v)+degree(z)-6 paws become stars; raw
        // diamonds and cliques stay fixed. Paths and the frequency denominator
        // can increase, so the outer optimizer still measures log error.
        // u-v-z survives the deletion and the detached leaf is reattached,
        // preserving connectivity without a bridge search.
        return try_operation({
            {EdgeAction::Remove, u, z},
            {EdgeAction::Remove, w, leaf},
            {EdgeAction::Add, u, w},
            {EdgeAction::Add, z, leaf},
        });
    }

    bool try_box_leaf_switch(int leaf, int a, int middle, int far) {
        // Start with the induced path leaf-a-middle-far. Exchange a-leaf
        // with b-branch for a-branch and b-leaf, preserving every degree.
        if (graph_.degree[leaf] != 1) return false;
        const int branch = random_neighbor(far);
        if (branch < 0 || graph_.degree[branch] < 2 ||
            graph_.degree[branch] > BOX_SWITCH_MAX_BRANCH_DEGREE ||
            graph_.matrix.get(a, branch)) return false;
        const int b = random_neighbor(branch);
        if (b < 0 || b == far || b == a || b == leaf) return false;

        // Both the removed b-branch edge and the new a-branch edge must
        // have no triangles. Check current adjacency, independent of the
        // freshness of the orbit table. No triangle, diamond, or K4 changes.
        int work = 0;
        for (int x : graph_.adj[branch]) {
            if (x == b) continue;
            if (graph_.matrix.get(a, x) || graph_.matrix.get(b, x)) return false;
            work += graph_.degree[x];
            if (work > BOX_SWITCH_NEIGHBOR_WORK_LIMIT) return false;
        }

        // C(r,x) denotes the current number of common neighbors. The exact
        // raw box increment is sum_x [C(a,x)-C(b,x)+1], x in N(branch)-{b}.
        // The +1 removes branch itself from the old common-neighbor count.
        // At most 128 neighbor entries are visited, never surrounding four-sets.
        int boxes = 0;
        bool connected = graph_.matrix.get(b, a) || graph_.matrix.get(b, middle);
        for (int x : graph_.adj[branch]) {
            if (x == b) continue;
            int common_a = 0, common_b = 0;
            for (int neighbor : graph_.adj[x]) {
                common_a += graph_.matrix.get(a, neighbor);
                common_b += graph_.matrix.get(b, neighbor);
            }
            boxes += common_a - common_b + 1;
            // A second common neighbor gives b-y-x-branch as an alternate
            // path. The two earlier tests instead use the sampled P4.
            if (common_b >= 2) connected = true;
        }

        // With degrees and triangles fixed, raw stars, paws, diamonds, and
        // cliques stay fixed. The raw P4 delta is walk_delta - 4*boxes.
        // Requiring walk_delta <= 3*boxes also keeps the total connected
        // four-set count from increasing, so the box percentage rises and
        // the P4 percentage falls. The outer optimizer measures log error.
        const long long walk_delta = (graph_.degree[a] - graph_.degree[b]) *
            (graph_.degree[branch] - 1LL);
        if (boxes <= 0 || walk_delta > 3LL * boxes || !connected) return false;
        return try_operation({
            {EdgeAction::Remove, a, leaf},
            {EdgeAction::Remove, b, branch},
            {EdgeAction::Add, a, branch},
            {EdgeAction::Add, b, leaf},
        });
    }

    bool try_balanced_path_leaf_relocation(int leaf, int from) {
        if (leaf < 0 || from < 0 || leaf >= graph_.n || from >= graph_.n ||
            graph_.degree[leaf] != 1 || !graph_.matrix.get(leaf, from)) return false;
        const int degree = graph_.degree[from];
        if (degree < 2 || degree > PATH_LEAF_ROOT_DEGREE_LIMIT) return false;
        const auto& candidates = path_leaf_targets_[degree - 1];
        if (candidates.empty()) return false;
        const int to = candidates[rng_.bounded(candidates.size())];
        if (to == leaf || to == from || graph_.degree[to] != degree - 1 ||
            graph_.matrix.get(leaf, to) ||
            triangle_counts_[to] != triangle_counts_[from]) return false;

        // The roots exchange degrees d and d-1, preserving the histogram.
        // Matching triangle counts also preserves raw stars and paws. As in
        // mode 68, triangle counts are exact for a fresh table followed only
        // by leaf moves; a frozen table or other modes make this a heuristic.
        // Boxes, diamonds, and cliques cannot contain the moved true leaf.
        long long from_walks = 0, to_walks = 0;
        for (int neighbor : graph_.adj[from]) {
            from_walks += graph_.degree[neighbor] - 1LL;
        }
        for (int neighbor : graph_.adj[to]) {
            to_walks += graph_.degree[neighbor] - 1LL;
        }
        // P(root) = sum_{neighbor}(degree(neighbor)-1) - 2*T(root).
        // Equal T cancels; adjacent roots need one correction for the leaf.
        // Only bounded one-hop degree sums are read, never neighbor pairs.
        const long long path_delta = to_walks - from_walks -
            static_cast<long long>(graph_.matrix.get(from, to));
        if (path_delta >= 0) return false;
        return try_operation({
            {EdgeAction::Remove, leaf, from},
            {EdgeAction::Add, leaf, to},
        });
    }

    bool try_pendant_triangle_opening(int u, int v, int hub) {
        if (u < 0 || v < 0 || hub < 0 || u >= graph_.n ||
            v >= graph_.n || hub >= graph_.n ||
            u == v || u == hub || v == hub ||
            graph_.degree[u] != 2 || graph_.degree[v] != 2 ||
            graph_.degree[hub] < 3 ||
            !graph_.matrix.get(u, v) ||
            !graph_.matrix.get(u, hub) || !graph_.matrix.get(v, hub)) {
            return false;
        }
        // Exactly degree(hub)-2 paws become claws. All other connected
        // four-node counts, and the total connected-four-set count, stay fixed.
        // u-hub-v keeps the removed edge's endpoints connected.
        return try_operation({{EdgeAction::Remove, u, v}});
    }

    void refresh_pendant_triangles() {
        if (pendant_triangles_revision_ == graph_revision_) return;
        pendant_triangles_.clear();
        for (int u = 0; u < graph_.n; ++u) {
            if (graph_.degree[u] != 2) continue;
            for (int index = 0; index < 2; ++index) {
                const int v = graph_.adj[u][index];
                const int hub = graph_.adj[u][1 - index];
                if (u < v && graph_.degree[v] == 2 &&
                    graph_.degree[hub] >= 3 && graph_.matrix.get(v, hub)) {
                    pendant_triangles_.push_back({u, v, hub});
                }
            }
        }
        pendant_triangles_revision_ = graph_revision_;
    }

    bool try_combined_leaf_opening(int leaf, int from, int to) {
        if (leaf < 0 || leaf >= graph_.n || graph_.degree[leaf] != 1) return false;
        refresh_pendant_triangles();
        if (pendant_triangles_.empty()) return false;
        if (!leaf_relocation_allowed(leaf, from, to, 61)) return false;

        const size_t count = pendant_triangles_.size();
        const size_t start = rng_.bounded(count);
        // Bound per-proposal retries. Failure simply consumes a sampling attempt.
        const size_t limit = min(count, size_t{16});
        for (size_t attempt = 0; attempt < limit; ++attempt) {
            const auto triangle = pendant_triangles_[(start + attempt) % count];
            const int u = triangle[0], v = triangle[1];
            bool disjoint = true;
            for (int node : triangle) {
                if (node == leaf || node == from || node == to) disjoint = false;
            }
            if (!disjoint) continue;

            // The leaf slide is already certified to lower G6 without raising
            // G7 or shrinking the connected-four-set denominator. None of the
            // triangle's vertices is touched by that slide. Removing u-v then
            // converts exactly degree(hub)-2 more paws to claws, with no other
            // four-graphlet changes. Both G6/G7 percentages drop strictly.
            // All three edits share the usual atomic degree/histogram checks.
            if (try_operation({
                    {EdgeAction::Add, leaf, to},
                    {EdgeAction::Remove, leaf, from},
                    {EdgeAction::Remove, u, v},
                })) return true;
        }
        return false;
    }

    void sample_one(int mode) {
        if (is_new_inverse(mode)) {
            auto& state = inverse_states_.at(mode);
            if (!state.enabled) return;
            ++state.attempts;
        }
        if (mode == 73) {
            if (!enabled73_) return;
            ++counters73_.sampling_attempts;
        }
        array<int, 5> nodes{};
        int count = 0;

        if (mode == 70) {
            if (box_switch_leaves_.empty()) return;
            nodes[count++] = box_switch_leaves_[rng_.bounded(box_switch_leaves_.size())];
            if (graph_.degree[nodes[0]] != 1) return;
        } else {
            nodes[count++] = pick_first_node();
        }
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
        if (is_new_inverse(mode)) {
            apply_learned_inverse(mode, nodes, sub);
            return;
        }

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

    void apply_from_clique(int mode, const array<int, 5>& sampled_nodes) {
        array<int, 5> nodes = sampled_nodes;
        if (mode >= 65) {
            // The clique was shuffled already. Put its eligible endpoints
            // first instead of rejecting a usable clique for its orientation.
            // Both loops inspect exactly four vertices, independent of degree.
            int count = 0;
            for (int i = 0; i < 4; ++i) {
                if (clique_sources_[sampled_nodes[i]]) nodes[count++] = sampled_nodes[i];
            }
            if (count < 2) return;
            for (int i = 0; i < 4; ++i) {
                if (!clique_sources_[sampled_nodes[i]]) nodes[count++] = sampled_nodes[i];
            }
        }
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
            case 65: {
                // Exchange uv for ux, where x shares the triangle edge wz.
                // The five-node witness is isomorphic before/after: its K4
                // changes from uvwz to uxwz. Prefer a lower-degree recipient
                // to reduce cliques outside the witness, without counting them.
                const int x = random_neighbor(w);
                if (x < 0 || contains_node(nodes, 4, x) ||
                    graph_.degree[x] >= graph_.degree[v] ||
                    !graph_.matrix.get(x, z) || graph_.matrix.get(x, u) ||
                    graph_.matrix.get(x, v)) break;
                try_operation({
                    {EdgeAction::Remove, u, v},
                    {EdgeAction::Add, u, x},
                });
                break;
            }
            case 66: {
                // K4 uvwz and induced C4 w-x-y-z-w share edge wz.
                // uv + xy -> ux + vy preserves every degree. On this six-node
                // witness it removes one K4 and one C4, adding three diamonds.
                const int x = random_neighbor(w);
                if (x < 0 || contains_node(nodes, 4, x) ||
                    graph_.matrix.get(x, u) || graph_.matrix.get(x, z)) break;
                const int y = random_neighbor(x);
                if (y < 0 || contains_node(nodes, 4, y) ||
                    !graph_.matrix.get(y, z) || graph_.matrix.get(y, v) ||
                    graph_.matrix.get(y, w)) break;
                // x-u-w-v-y replaces xy; u-w-v replaces uv. No bridge search.
                try_operation({
                    {EdgeAction::Remove, u, v},
                    {EdgeAction::Remove, x, y},
                    {EdgeAction::Add, u, x},
                    {EdgeAction::Add, v, y},
                });
                break;
            }
            case 67: {
                // Transfer uv into diagonal wy of the adjoining C4 w-x-y-z-w.
                // The source K4 and the box both become diamonds. x-u and x-v
                // may exist; the sampled six-node witness still loses K4s and
                // gains no boxes. Effects outside it are measured, not scanned.
                const int x = random_neighbor(w);
                if (x < 0 || contains_node(nodes, 4, x) ||
                    graph_.matrix.get(x, z)) break;
                const int y = random_neighbor(x);
                if (y < 0 || contains_node(nodes, 4, y) ||
                    !graph_.matrix.get(y, z) || graph_.matrix.get(y, u) ||
                    graph_.matrix.get(y, v) || graph_.matrix.get(y, w)) break;
                // u-w-v survives the deletion, preserving connectivity.
                try_operation({
                    {EdgeAction::Remove, u, v},
                    {EdgeAction::Add, w, y},
                });
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
            case 72:
                if (try_triangle_branch_exchange(u, v, z, w)) break;
                try_triangle_branch_exchange(v, u, z, w);
                break;
            case 71:
                if (try_paw_box_switch(u, v, z, w)) break;
                try_paw_box_switch(v, u, z, w);
                break;
            case 62: {
                if (graph_.degree[w] != 1) break;
                // Try one nearby root and one degree-biased global root. Both
                // use the exact global four-graphlet-delta check above.
                const int nearby = random_neighbor(z);
                if (nearby >= 0 && try_leaf_relocation(w, z, nearby, 62)) break;
                const int seed = static_cast<int>(rng_.bounded(graph_.n));
                const int destination = random_neighbor(seed);
                if (destination >= 0) try_leaf_relocation(w, z, destination, 62);
                break;
            }
            case 63:
                try_pendant_triangle_opening(u, v, z);
                break;
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
            case 73: {
                ++counters73_.path_candidates;
                // Turn u-w-v-z into a star centered at v. Provenance scores
                // are only a heuristic; no graphlet-effect check is performed.
                if (try_model_relocation73(u,w,v,z)) break;
                if (config73_.both_orientations) try_model_relocation73(z,v,w,u);
                break;
            }
            case 70:
                if (try_box_leaf_switch(u, w, v, z)) break;
                try_box_leaf_switch(z, v, w, u);
                break;
            case 69:
                // Try the true endpoint leaves of the sampled induced P4.
                if (try_balanced_path_leaf_relocation(u, w)) break;
                try_balanced_path_leaf_relocation(z, v);
                break;
            case 68: {
                // Start with a true endpoint leaf of the sampled P4, then
                // draw one global destination using the existing orbit bias.
                // A global draw avoids restricting moves to adjacent roots.
                const int destination = pick_first_node();
                if (try_triangle_guided_leaf_relocation(u, w, destination)) break;
                try_triangle_guided_leaf_relocation(z, v, destination);
                break;
            }
            case 64:
                if (try_combined_leaf_opening(z, v, w)) break;
                try_combined_leaf_opening(u, w, v);
                break;
            case 61:
                // On u-w-v-z, lift a true endpoint leaf one step inward.
                // Test both orientations but commit at most one change.
                if (try_leaf_relocation(z, v, w, 61)) break;
                try_leaf_relocation(u, w, v, 61);
                break;
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
	    static ofstream moddeg("moddeg");
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
            case 20:{
	                bool ididityay = try_operation({
	                    {EdgeAction::Add, v, w},
	                    {EdgeAction::Remove, z, w},
	                });
					if (ididityay) moddeg << graph_.degree[u] << ' ' << graph_.degree[z] << ' '
	                       << graph_.degree[v] << ' ' << graph_.degree[w] << '\n';
	            }
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
        vector<uint8_t> clique_sources;
        vector<long long> triangle_counts;
        const vector<uint8_t> good_mask =
            read_good_node_mask(graph.n, clique_sources, triangle_counts);

        const uint64_t seed = make_seed();
        Sampler sampler(
            graph,
            target_degree,
            degree_log_bound,
            hdgbound,
            good_mask,
            seed,
            move(clique_sources),
            move(triangle_counts));
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
