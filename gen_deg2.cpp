// Construct a fresh connected simple graph using only its requested degrees.
// Build: g++ -O3 -std=c++17 gen_deg_new.cpp -o gen_deg_new
#include <algorithm>
#include <chrono>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <random>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {
using Node = uint32_t;
using Edge = std::pair<Node, Node>;
using Clock = std::chrono::steady_clock;
// Guard against tiny malformed histogram files requesting enormous arrays.
constexpr uint64_t MAX_NODES = 10000000;

uint64_t integer(const std::string& text, const char* name) {
    if (text.empty()) throw std::runtime_error(std::string(name) + " is empty");
    uint64_t value = 0;
    for (char c : text) {
        if (c < '0' || c > '9')
            throw std::runtime_error(std::string(name) + " must be a nonnegative decimal integer");
        const unsigned digit = static_cast<unsigned>(c - '0');
        if (value > (std::numeric_limits<uint64_t>::max() - digit) / 10)
            throw std::runtime_error(std::string(name) + " is out of range");
        value = value * 10 + digit;
    }
    return value;
}

class FastRng {
    uint64_t state_;
public:
    explicit FastRng(uint64_t seed) : state_(seed) {}
    uint64_t next() {
        uint64_t z = (state_ += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }
    uint64_t bounded(uint64_t bound) {
        __uint128_t product = static_cast<__uint128_t>(next()) * bound;
        uint64_t low = static_cast<uint64_t>(product);
        if (low < bound) {
            const uint64_t threshold = -bound % bound;
            while (low < threshold) {
                product = static_cast<__uint128_t>(next()) * bound;
                low = static_cast<uint64_t>(product);
            }
        }
        return static_cast<uint64_t>(product >> 64);
    }
};

uint64_t random_seed() {
    std::random_device rd;
    return static_cast<uint64_t>(Clock::now().time_since_epoch().count()) ^
           (static_cast<uint64_t>(rd()) << 32) ^ static_cast<uint64_t>(rd());
}

struct EdgeHash {
    size_t operator()(uint64_t z) const {
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return static_cast<size_t>(z ^ (z >> 31));
    }
};
using EdgeSet = std::unordered_set<uint64_t, EdgeHash>;

Edge ordered(Node u, Node v) { return {std::min(u, v), std::max(u, v)}; }
uint64_t key(const Edge& edge) {
    return (static_cast<uint64_t>(edge.first) << 32) | edge.second;
}

// Parse exactly two unsigned integers per nonblank line, without allocating
// token strings for each of the millions of input edges.
bool pair_on_line(const std::string& line, uint64_t& a, uint64_t& b, size_t number) {
    size_t pos = 0;
    auto skip_space = [&] {
        while (pos < line.size() && std::isspace(static_cast<unsigned char>(line[pos]))) ++pos;
    };
    skip_space();
    if (pos == line.size()) return false;
    uint64_t values[2]{};
    for (int column = 0; column < 2; ++column) {
        if (pos == line.size() || line[pos] < '0' || line[pos] > '9')
            throw std::runtime_error("expected two nonnegative integers on line " + std::to_string(number));
        while (pos < line.size() && line[pos] >= '0' && line[pos] <= '9') {
            const unsigned digit = static_cast<unsigned>(line[pos++] - '0');
            if (values[column] > (std::numeric_limits<uint64_t>::max() - digit) / 10)
                throw std::runtime_error("integer out of range on line " + std::to_string(number));
            values[column] = values[column] * 10 + digit;
        }
        if (pos < line.size() && !std::isspace(static_cast<unsigned char>(line[pos])))
            throw std::runtime_error("invalid integer on line " + std::to_string(number));
        skip_space();
    }
    if (pos != line.size())
        throw std::runtime_error("extra data after two integers on line " + std::to_string(number));
    a = values[0]; b = values[1];
    return true;
}

enum class InputMode { Edges, Degrees, Histogram };

std::vector<Node> read_degrees(const std::string& path, InputMode mode) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open input '" + path + "'");
    std::vector<Node> degrees;
    if (mode == InputMode::Degrees) {
        std::string token;
        while (input >> token) {
            const uint64_t degree = integer(token, "degree");
            if (degree >= MAX_NODES || degrees.size() >= MAX_NODES)
                throw std::runtime_error("degree input exceeds the 10000000-node limit");
            degrees.push_back(static_cast<Node>(degree));
        }
    } else {
        std::unordered_map<uint64_t, Node> ids;
        EdgeSet seen;
        seen.max_load_factor(0.85F);
        auto node_id = [&](uint64_t label) {
            const auto found = ids.find(label);
            if (found != ids.end()) return found->second;
            if (degrees.size() >= MAX_NODES)
                throw std::runtime_error("input exceeds the 10000000-node limit");
            const Node id = static_cast<Node>(degrees.size());
            ids.emplace(label, id);
            degrees.push_back(0);
            return id;
        };
        std::string line;
        size_t number = 0;
        while (std::getline(input, line)) {
            ++number;
            uint64_t a, b;
            if (!pair_on_line(line, a, b, number)) continue;
            if (mode == InputMode::Histogram) {
                if (a >= MAX_NODES || b > MAX_NODES - degrees.size())
                    throw std::runtime_error("histogram exceeds the 10000000-node limit");
                degrees.insert(degrees.end(), static_cast<size_t>(b), static_cast<Node>(a));
            } else {
                if (a == b) continue;
                const Node u = node_id(a), v = node_id(b);
                if (seen.insert(key(ordered(u, v))).second) {
                    ++degrees[u]; ++degrees[v];
                }
            }
        }
        // Input topology and labels are discarded here. Canonicalizing the
        // multiset makes equal degree distributions + equal seed reproducible
        // even when supplied through different source graphs or histograms.
        std::sort(degrees.begin(), degrees.end());
    }
    if (input.bad()) throw std::runtime_error("error reading input '" + path + "'");
    if (degrees.empty()) throw std::runtime_error("input specifies no nodes");
    return degrees;
}

uint64_t validate_degrees(const std::vector<Node>& degrees, bool header,
                          bool allow_disconnected) {
    uint64_t sum = 0;
    bool isolated = false;
    for (Node degree : degrees) {
        if (degree >= degrees.size())
            throw std::runtime_error("every degree must be smaller than the node count");
        sum += degree;
        isolated = isolated || degree == 0;
    }
    if (sum % 2 != 0) throw std::runtime_error("degree sum must be even");
    if (isolated && !header)
        throw std::runtime_error("isolated nodes require --header; an edge list alone cannot encode them");
    if (!allow_disconnected && degrees.size() > 1) {
        if (isolated)
            throw std::runtime_error("a connected graph cannot have isolated nodes; use --allow-disconnected");
        if (sum / 2 < degrees.size() - 1)
            throw std::runtime_error("a connected graph needs at least n-1 edges; use --allow-disconnected");
    }
    return sum / 2;
}

struct HeapNode {
    uint64_t tie;
    Node vertex;
    Node degree;
    bool operator<(const HeapNode& other) const {
        if (degree != other.degree) return degree < other.degree;
        if (tie != other.tie) return tie < other.tie;
        return vertex < other.vertex;
    }
};

// With existing != nullptr, complete residual budgets while excluding the
// edges of the initial tree. That constrained greedy attempt can fail even
// for a graphical target; callers retry or use the unconstrained fallback.
bool havel_hakimi(const std::vector<Node>& degrees, FastRng& rng,
                  std::vector<Edge>& edges, EdgeSet* existing = nullptr) {
    std::vector<HeapNode> initial;
    initial.reserve(degrees.size());
    for (size_t u = 0; u < degrees.size(); ++u)
        if (degrees[u] != 0)
            initial.push_back({rng.next(), static_cast<Node>(u), degrees[u]});
    std::priority_queue<HeapNode> heap(std::less<HeapNode>(), std::move(initial));
    std::vector<HeapNode> pending;
    pending.reserve(*std::max_element(degrees.begin(), degrees.end()));
    std::vector<HeapNode> blocked;
    while (!heap.empty()) {
        const HeapNode u = heap.top();
        heap.pop();
        if (u.degree > heap.size()) return false;
        pending.clear();
        blocked.clear();
        // Keep these vertices outside the heap until all neighbors have been
        // chosen, so each is distinct. Havel--Hakimi guarantees a realization
        // for every graphical sequence; ties may be broken randomly.
        Node chosen = 0;
        while (chosen < u.degree && !heap.empty()) {
            HeapNode v = heap.top();
            heap.pop();
            const Edge edge = ordered(u.vertex, v.vertex);
            if (existing && existing->count(key(edge))) {
                blocked.push_back(v);
                continue;
            }
            edges.push_back(edge);
            if (existing) existing->insert(key(edge));
            ++chosen;
            --v.degree;
            if (v.degree != 0) pending.push_back(v);
        }
        if (chosen != u.degree) return false;
        for (const HeapNode& v : pending) heap.push(v);
        for (const HeapNode& v : blocked) heap.push(v);
    }
    return true;
}

struct ConstructionStats {
    uint64_t tree_attempts = 0;
    bool used_tree = false;
};

std::vector<Edge> realize(const std::vector<Node>& degrees, uint64_t m,
                          uint64_t tree_attempts, FastRng& rng,
                          ConstructionStats& stats) {
    std::vector<Edge> edges;
    if (m > edges.max_size()) throw std::runtime_error("too many edges for this system");
    edges.reserve(static_cast<size_t>(m));
    for (uint64_t attempt = 0; attempt < tree_attempts && degrees.size() > 1; ++attempt) {
        ++stats.tree_attempts;
        edges.clear();
        std::vector<Node> remaining = degrees;
        std::vector<Node> order(degrees.size());
        std::iota(order.begin(), order.end(), Node{0});
        std::vector<uint64_t> ties(degrees.size());
        for (uint64_t& tie : ties) tie = rng.next();
        std::sort(order.begin(), order.end(), [&](Node a, Node b) {
            if (degrees[a] != degrees[b]) return degrees[a] > degrees[b];
            if (ties[a] != ties[b]) return ties[a] > ties[b];
            return a < b;
        });

        std::vector<Node> available;
        available.reserve(degrees.size());
        available.push_back(order[0]);
        for (size_t position = 1; position < order.size(); ++position) {
            // Descending degrees ensure the pool cannot empty before all
            // nodes have joined, given positive degrees and m >= n-1.
            if (available.empty()) throw std::runtime_error("internal tree capacity error");
            const size_t parent_index = static_cast<size_t>(rng.bounded(available.size()));
            const Node parent = available[parent_index], child = order[position];
            edges.push_back(ordered(parent, child));
            if (--remaining[parent] == 0) {
                available[parent_index] = available.back();
                available.pop_back();
            }
            if (--remaining[child] != 0) available.push_back(child);
        }

        EdgeSet existing;
        existing.max_load_factor(0.85F);
        existing.reserve(static_cast<size_t>(m));
        for (const Edge& edge : edges) existing.insert(key(edge));
        // Only initial-tree edges can block two still-active HH vertices.
        // Each is examined at most once this way: there is no all-pairs scan.
        if (havel_hakimi(remaining, rng, edges, &existing)) {
            if (edges.size() != m) throw std::runtime_error("internal tree completion error");
            stats.used_tree = true;
            return edges;
        }
    }
    edges.clear();
    if (!havel_hakimi(degrees, rng, edges))
        throw std::runtime_error("degree sequence is not graphical: no simple graph realizes it");
    if (edges.size() != m) throw std::runtime_error("internal degree realization error");
    return edges;
}

struct SwapStats { uint64_t attempts = 0, accepted = 0; };

SwapStats randomize(std::vector<Edge>& edges, uint64_t passes, FastRng& rng) {
    SwapStats stats;
    if (edges.size() < 2 || passes == 0) return stats;
    if (passes > std::numeric_limits<uint64_t>::max() / edges.size())
        throw std::runtime_error("swap attempt budget is out of range");
    const uint64_t budget = passes * edges.size();
    EdgeSet existing;
    existing.max_load_factor(0.85F);
    existing.reserve(edges.size());
    for (const Edge& edge : edges) existing.insert(key(edge));
    for (uint64_t attempt = 0; attempt < budget; ++attempt) {
        ++stats.attempts;
        const size_t i = static_cast<size_t>(rng.bounded(edges.size()));
        const size_t j = static_cast<size_t>(rng.bounded(edges.size()));
        if (i == j) continue;
        Edge a = edges[i], b = edges[j];
        if (rng.next() & 1) std::swap(a.first, a.second);
        if (rng.next() & 1) std::swap(b.first, b.second);
        if (a.first == b.first || a.first == b.second ||
            a.second == b.first || a.second == b.second) continue;
        const Edge next_a = ordered(a.first, b.first);
        const Edge next_b = ordered(a.second, b.second);
        if (existing.count(key(next_a)) || existing.count(key(next_b))) continue;
        // Reuse the two hash nodes instead of allocating on every accepted
        // swap. Each switch preserves every degree and keeps the graph simple.
        auto first = existing.extract(key(edges[i]));
        auto second = existing.extract(key(edges[j]));
        first.value() = key(next_a);
        second.value() = key(next_b);
        existing.insert(std::move(first));
        existing.insert(std::move(second));
        edges[i] = next_a;
        edges[j] = next_b;
        ++stats.accepted;
    }
    return stats;
}

class DisjointSet {
    std::vector<Node> parent_, size_;
public:
    size_t count;
    explicit DisjointSet(size_t n) : parent_(n), size_(n, 1), count(n) {
        std::iota(parent_.begin(), parent_.end(), Node{0});
    }
    Node find(Node v) {
        while (parent_[v] != v) {
            parent_[v] = parent_[parent_[v]];
            v = parent_[v];
        }
        return v;
    }
    bool unite(Node a, Node b) {
        a = find(a); b = find(b);
        if (a == b) return false;
        if (size_[a] < size_[b]) std::swap(a, b);
        parent_[b] = a;
        size_[a] += size_[b];
        --count;
        return true;
    }
};

// A random switch may disconnect even a successfully completed initial tree.
// Merge components without changing degrees: remove a cycle edge (a,b) in
// one component and a forest edge (c,d) in another, then add (a,c),(b,d).
// The first component stays connected; both sides of the second reattach.
// New edges cross components and hence cannot be duplicates or self-loops.
size_t connect_components(std::vector<Edge>& edges, size_t n) {
    if (n == 1) return 0;
    DisjointSet components(n);
    std::vector<size_t> forest, spare;
    forest.reserve(n - 1);
    spare.reserve(n - 1);
    for (size_t i = 0; i < edges.size(); ++i) {
        if (components.unite(edges[i].first, edges[i].second)) forest.push_back(i);
        else if (spare.size() < n - 1) spare.push_back(i);
        if (components.count == 1) return 0;
    }

    std::vector<Node> roots, root_position(n);
    std::vector<size_t> representative(n);
    for (size_t u = 0; u < n; ++u) {
        const Node v = static_cast<Node>(u);
        if (components.find(v) == v) {
            root_position[v] = static_cast<Node>(roots.size());
            roots.push_back(v);
        }
    }
    for (size_t i : forest) representative[components.find(edges[i].first)] = i;
    const size_t merges = roots.size() - 1;
    // There are m-n+c cycle edges for c components. Because m>=n-1,
    // at least c-1 donors exist; retaining at most n-1 is sufficient.
    if (spare.size() < merges) throw std::runtime_error("internal connectivity budget error");
    auto remove_root = [&](Node root) {
        const size_t position = root_position[root];
        roots[position] = roots.back();
        root_position[roots[position]] = static_cast<Node>(position);
        roots.pop_back();
    };
    while (roots.size() > 1) {
        const size_t donor = spare.back();
        spare.pop_back();
        const Edge a = edges[donor];
        const Node root_a = components.find(a.first);
        const Node root_b = roots.front() == root_a ? roots.back() : roots.front();
        const size_t recipient = representative[root_b];
        const Edge b = edges[recipient];
        edges[donor] = ordered(a.first, b.first);
        edges[recipient] = ordered(a.second, b.second);
        remove_root(root_a);
        remove_root(root_b);
        components.unite(root_a, root_b);
        const Node merged = components.find(root_a);
        representative[merged] = donor;
        root_position[merged] = static_cast<Node>(roots.size());
        roots.push_back(merged);
    }
    return merges;
}

void check_distinct_paths(const char* input, const char* output) {
    struct stat a{}, b{};
    if (::stat(input, &a) == 0 && ::stat(output, &b) == 0 &&
        a.st_dev == b.st_dev && a.st_ino == b.st_ino)
        throw std::runtime_error("input and output must be different files");
}

void write_graph(const char* path, size_t n, const std::vector<Edge>& edges, bool header) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error(std::string("cannot open output '") + path + "'");
    if (header) output << n << ' ' << edges.size() << '\n';
    for (const Edge& edge : edges) output << edge.first << ' ' << edge.second << '\n';
    output.close();
    if (!output) throw std::runtime_error(std::string("error writing output '") + path + "'");
}

void usage(std::ostream& out, const char* name) {
    out << "Usage: " << name << " INPUT.el OUTPUT.el [OPTIONS]\n"
        << "   or: " << name << " --degrees DEGREES.txt OUTPUT.el [OPTIONS]\n"
        << "   or: " << name << " --histogram HISTOGRAM.txt OUTPUT.el [OPTIONS]\n"
        << "  --degrees: whitespace-separated degrees, one per node in ID order.\n"
        << "  --histogram: lines containing degree and number of nodes.\n"
        << "  --seed UINT64          Reproducible random seed.\n"
        << "  --swaps-per-edge UINT  Random switch attempts per edge (default 5).\n"
        << "  --tree-attempts UINT   Tree-first attempts before HH fallback (default 4).\n"
        << "  --allow-disconnected  Disable the default connected-output requirement.\n"
        << "  --header               Write an n m header; required for isolated nodes.\n"
        << "Fresh node IDs are 0..n-1. Output is connected by default.\n"
        << "Input graph mode uses only its simple-graph degree distribution.\n"
        << "With default budgets: O((n+m) log n) time, O(n+m) memory; at most 10000000 nodes.\n";
}

double seconds(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration<double>(b - a).count();
}
} // namespace

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    if (argc == 2 && std::string(argv[1]) == "--help") {
        usage(std::cout, argv[0]); return 0;
    }
    try {
        InputMode mode = InputMode::Edges;
        int position = 1;
        if (argc > 1 && std::string(argv[1]) == "--degrees") {
            mode = InputMode::Degrees; ++position;
        } else if (argc > 1 && std::string(argv[1]) == "--histogram") {
            mode = InputMode::Histogram; ++position;
        }
        if (argc < position + 2) { usage(std::cerr, argv[0]); return 1; }
        const char* input = argv[position++];
        const char* output = argv[position++];
        bool header = false, seeded = false, allow_disconnected = false;
        uint64_t seed = 0, swaps = 5, tree_attempts = 4;
        while (position < argc) {
            const std::string option = argv[position++];
            if (option == "--header") header = true;
            else if (option == "--allow-disconnected") allow_disconnected = true;
            else if (option == "--seed" || option == "--swaps-per-edge" || option == "--tree-attempts") {
                if (position == argc) throw std::runtime_error("missing value for " + option);
                const uint64_t value = integer(argv[position++], option.c_str());
                if (option == "--seed") { seed = value; seeded = true; }
                else if (option == "--tree-attempts") tree_attempts = value;
                else swaps = value;
            } else throw std::runtime_error("unknown option '" + option + "'");
        }
        if (!seeded) seed = random_seed();
        check_distinct_paths(input, output);
        const auto start = Clock::now();
        std::vector<Node> degrees = read_degrees(input, mode);
        const uint64_t m = validate_degrees(degrees, header, allow_disconnected);
        FastRng rng(seed);
        const auto input_done = Clock::now();
        ConstructionStats construction;
        std::vector<Edge> edges = realize(degrees, m, allow_disconnected ? 0 : tree_attempts,
                                          rng, construction);
        const auto constructed = Clock::now();
        const SwapStats stats = randomize(edges, swaps, rng);
        const auto randomized = Clock::now();
        const size_t component_merges = allow_disconnected ? 0 : connect_components(edges, degrees.size());
        const auto connected = Clock::now();
        write_graph(output, degrees.size(), edges, header);
        const auto written = Clock::now();
        std::cerr << "gen_deg_new: nodes=" << degrees.size() << " edges=" << edges.size()
                  << " seed=" << seed << " swap_attempts=" << stats.attempts
                  << " swaps_accepted=" << stats.accepted
                  << " construction=" << (construction.used_tree ? "tree_hh" : "hh")
                  << " tree_attempts=" << construction.tree_attempts
                  << " connectivity=" << (allow_disconnected ? "not_enforced" : "connected")
                  << " component_merges=" << component_merges << std::fixed << std::setprecision(3)
                  << " input_seconds=" << seconds(start, input_done)
                  << " construct_seconds=" << seconds(input_done, constructed)
                  << " randomize_seconds=" << seconds(constructed, randomized)
                  << " connect_seconds=" << seconds(randomized, connected)
                  << " write_seconds=" << seconds(connected, written) << '\n';
    } catch (const std::exception& error) {
        std::cerr << "gen_deg_new: " << error.what() << '\n'; return 1;
    }
    return 0;
}
