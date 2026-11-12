// Sample binary presence of ORCA's connected node orbits on 2--4 vertices.
// Each 1 is witnessed by an induced graphlet; sampling can miss rare orbits.
// Build: g++ -O2 -std=c++17 sampled_orbits.cpp -o sampled_orbits

#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <utility>
#include <vector>

namespace {

using Node = uint32_t;
using Edge = std::pair<Node, Node>;
using Clock = std::chrono::steady_clock;

constexpr uint64_t DEFAULT_SAMPLES = 10000000;
// Rows preserve numeric IDs. Refuse enormous ID gaps before allocating CSR
// offsets or the result table (10 million rows already produce 300 MB text).
constexpr uint64_t MAX_NODES = 10000000;
constexpr uint64_t MAX_NODE_ID =
    static_cast<uint64_t>(std::numeric_limits<int>::max()) - 1;

uint64_t decimal_integer(const std::string& text, const char* name,
                         uint64_t maximum) {
    if (text.empty()) throw std::runtime_error(std::string(name) + " is empty");
    uint64_t value = 0;
    for (char c : text) {
        if (c < '0' || c > '9') {
            throw std::runtime_error(std::string(name) +
                                     " must be an unsigned decimal integer");
        }
        const unsigned digit = static_cast<unsigned>(c - '0');
        if (value > maximum / 10 ||
            (value == maximum / 10 && digit > maximum % 10)) {
            throw std::runtime_error(std::string(name) + " is out of range");
        }
        value = value * 10 + digit;
    }
    return value;
}

class FastRng {
public:
    explicit FastRng(uint64_t seed) : state_(seed) {}

    uint64_t next() {
        // SplitMix64; deterministic even for seed 0.
        uint64_t z = (state_ += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }

    uint64_t bounded(uint64_t bound) {
        // Multiply-high reduction with rejection eliminates modulo bias.
        // Callers always supply a positive bound.
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

private:
    uint64_t state_;
};

uint64_t random_seed() {
    std::random_device device;
    uint64_t seed = static_cast<uint64_t>(Clock::now().time_since_epoch().count());
    seed ^= static_cast<uint64_t>(device()) << 32;
    seed ^= static_cast<uint64_t>(device());
    return seed;
}

struct Graph {
    std::vector<Edge> edges;
    std::vector<size_t> offsets;
    std::vector<Node> neighbors;

    size_t node_count() const { return offsets.size() - 1; }
    size_t degree(Node v) const { return offsets[v + 1] - offsets[v]; }

    bool adjacent(Node u, Node v) const {
        // Search the smaller list; every CSR row is sorted.
        if (degree(u) > degree(v)) std::swap(u, v);
        return std::binary_search(neighbors.begin() + offsets[u],
                                  neighbors.begin() + offsets[u + 1], v);
    }
};

bool whitespace(char c) {
    return std::isspace(static_cast<unsigned char>(c)) != 0;
}

bool read_endpoints(const std::string& line, Node& u, Node& v, size_t number) {
    size_t pos = 0;
    while (pos < line.size() && whitespace(line[pos])) ++pos;
    if (pos == line.size()) return false;

    uint64_t endpoints[2]{};
    for (int column = 0; column < 2; ++column) {
        if (pos == line.size() || line[pos] < '0' || line[pos] > '9') {
            throw std::runtime_error("expected two nonnegative integer node IDs on line " +
                                     std::to_string(number));
        }
        uint64_t value = 0;
        while (pos < line.size() && line[pos] >= '0' && line[pos] <= '9') {
            const unsigned digit = static_cast<unsigned>(line[pos++] - '0');
            if (value > MAX_NODE_ID / 10 ||
                (value == MAX_NODE_ID / 10 && digit > MAX_NODE_ID % 10)) {
                throw std::runtime_error("node ID out of range on line " +
                                         std::to_string(number));
            }
            value = value * 10 + digit;
        }
        endpoints[column] = value;
        if (pos < line.size() && !whitespace(line[pos])) {
            throw std::runtime_error("invalid node ID on line " + std::to_string(number));
        }
        while (pos < line.size() && whitespace(line[pos])) ++pos;
    }
    if (pos != line.size()) {
        throw std::runtime_error("extra data after two node IDs on line " +
                                 std::to_string(number));
    }
    u = static_cast<Node>(endpoints[0]);
    v = static_cast<Node>(endpoints[1]);
    return true;
}

Graph read_graph(const std::string& filename) {
    std::ifstream input(filename);
    if (!input) throw std::runtime_error("cannot open graph '" + filename + "'");

    Graph graph;
    std::string line;
    size_t line_number = 0;
    size_t nodes = 0;
    while (std::getline(input, line)) {
        ++line_number;
        Node u, v;
        if (!read_endpoints(line, u, v, line_number) || u == v) continue;
        if (u > v) std::swap(u, v);
        if (static_cast<uint64_t>(v) >= MAX_NODES) {
            throw std::runtime_error("node IDs require more than " +
                std::to_string(MAX_NODES) + " rows (MAX_NODES) on line " +
                std::to_string(line_number));
        }
        nodes = std::max(nodes, static_cast<size_t>(v) + 1);
        graph.edges.emplace_back(u, v);
    }
    if (input.bad()) throw std::runtime_error("error reading graph '" + filename + "'");
    if (graph.edges.empty()) throw std::runtime_error("graph has no non-loop edges");

    std::sort(graph.edges.begin(), graph.edges.end());
    graph.edges.erase(std::unique(graph.edges.begin(), graph.edges.end()), graph.edges.end());
    if (graph.edges.size() > std::numeric_limits<size_t>::max() / 2 ||
        2 * graph.edges.size() > graph.neighbors.max_size()) {
        throw std::runtime_error("too many graph edges");
    }
    graph.offsets.assign(nodes + 1, 0);
    for (const Edge& edge : graph.edges) {
        ++graph.offsets[edge.first + 1];
        ++graph.offsets[edge.second + 1];
    }
    for (size_t v = 1; v <= nodes; ++v) graph.offsets[v] += graph.offsets[v - 1];
    graph.neighbors.resize(2 * graph.edges.size());
    std::vector<size_t> cursor(graph.offsets.begin(), graph.offsets.end() - 1);
    // Canonically sorted edges append each row's smaller neighbors first,
    // followed by its larger neighbors, both in increasing order.
    for (const Edge& edge : graph.edges) {
        graph.neighbors[cursor[edge.first]++] = edge.second;
        graph.neighbors[cursor[edge.second]++] = edge.first;
    }
    return graph;
}

constexpr unsigned EDGE_BITS[4][4] = {
    {0, 1, 2, 4}, {1, 0, 8, 16}, {2, 8, 0, 32}, {4, 16, 32, 0}
};

using OrbitRow = std::array<uint16_t, 4>;
using OrbitTable = std::array<OrbitRow, 64>;

void mark(OrbitRow& row, int vertex, int orbit) {
    row[vertex] |= static_cast<uint16_t>(1u << orbit);
}

OrbitTable make_orbit_table() {
    OrbitTable table{};
    for (unsigned mask = 0; mask < 64; ++mask) {
        OrbitRow& row = table[mask];
        int degree[4]{};
        int edges = 0;
        for (int i = 0; i < 4; ++i) {
            for (int j = i + 1; j < 4; ++j) {
                if (!(mask & EDGE_BITS[i][j])) continue;
                ++edges;
                ++degree[i];
                ++degree[j];
                mark(row, i, 0);
                mark(row, j, 0);
            }
        }
        // All connected induced triples, including those that were not a
        // prefix of the sampled expansion. ORCA 1/2 are P3 end/center; 3 is K3.
        for (int omitted = 0; omitted < 4; ++omitted) {
            int triple_degree[4]{};
            int triple_edges = 0;
            for (int i = 0; i < 4; ++i) if (i != omitted) {
                for (int j = i + 1; j < 4; ++j) if (j != omitted) {
                    if (!(mask & EDGE_BITS[i][j])) continue;
                    ++triple_edges;
                    ++triple_degree[i];
                    ++triple_degree[j];
                }
            }
            if (triple_edges < 2) continue;
            for (int i = 0; i < 4; ++i) if (i != omitted) {
                mark(row, i, triple_edges == 3 ? 3 : (triple_degree[i] == 1 ? 1 : 2));
            }
        }
        unsigned reached = 1;
        for (int pass = 0; pass < 3; ++pass) {
            for (int i = 0; i < 4; ++i) if (reached & (1u << i)) {
                for (int j = 0; j < 4; ++j) {
                    if (mask & EDGE_BITS[i][j]) reached |= 1u << j;
                }
            }
        }
        if (reached != 15) continue;

        const int max_degree = *std::max_element(degree, degree + 4);
        for (int i = 0; i < 4; ++i) {
            int orbit = -1;
            if (edges == 3) { // P4 or K1,3.
                orbit = max_degree == 3 ? (degree[i] == 3 ? 7 : 6)
                                        : (degree[i] == 1 ? 4 : 5);
            } else if (edges == 4) { // C4 or paw.
                orbit = max_degree == 2 ? 8
                    : (degree[i] == 1 ? 9 : (degree[i] == 2 ? 10 : 11));
            } else if (edges == 5) { // Diamond.
                orbit = degree[i] == 2 ? 12 : 13;
            } else if (edges == 6) { // K4.
                orbit = 14;
            }
            if (orbit < 0) throw std::runtime_error("invalid connected graphlet table entry");
            mark(row, i, orbit);
        }
    }
    return table;
}

bool selected(Node candidate, const std::array<Node, 4>& nodes, unsigned count) {
    for (unsigned i = 0; i < count; ++i) if (nodes[i] == candidate) return true;
    return false;
}

struct Extension {
    Node node;
    unsigned source;
};

Extension grow(const Graph& graph, const std::array<Node, 4>& nodes,
               unsigned count, uint64_t degree_sum, uint64_t external, FastRng& rng) {
    // Each incidence leaving the selected vertices is equally likely. Reject
    // edges back into the tuple, as in BLANT's edge-based expansion. With at
    // most three selected vertices and a nonempty boundary, success probability
    // is at least 1/7. A bounded fallback also handles dense tiny components.
    for (unsigned attempt = 0; attempt < 32; ++attempt) {
        uint64_t index = rng.bounded(degree_sum);
        unsigned source = 0;
        while (index >= graph.degree(nodes[source])) {
            index -= graph.degree(nodes[source++]);
        }
        const Node candidate = graph.neighbors[graph.offsets[nodes[source]] + index];
        if (!selected(candidate, nodes, count)) return {candidate, source};
    }
    // Conditional on reaching this fallback, choose a uniformly random actual
    // boundary incidence. This preserves the same sampling distribution.
    uint64_t index = rng.bounded(external);
    for (unsigned source = 0; source < count; ++source) {
        for (size_t pos = graph.offsets[nodes[source]];
             pos < graph.offsets[nodes[source] + 1]; ++pos) {
            const Node candidate = graph.neighbors[pos];
            if (selected(candidate, nodes, count)) continue;
            if (index == 0) return {candidate, source};
            --index;
        }
    }
    throw std::runtime_error("inconsistent graph boundary");
}

struct SampleResult {
    std::vector<uint16_t> presence;
    std::array<uint64_t, 5> sizes{};
};

SampleResult sample(const Graph& graph, uint64_t samples, uint64_t seed) {
    SampleResult result;
    result.presence.assign(graph.node_count(), 0);
    const OrbitTable table = make_orbit_table();
    FastRng rng(seed);
    for (uint64_t iteration = 0; iteration < samples; ++iteration) {
        const Edge& edge = graph.edges[rng.bounded(graph.edges.size())];
        std::array<Node, 4> nodes{{edge.first, edge.second, 0, 0}};
        unsigned count = 2;
        unsigned mask = EDGE_BITS[0][1];
        uint64_t degree_sum = graph.degree(nodes[0]) + graph.degree(nodes[1]);
        uint64_t internal_incidences = 2;
        while (count < 4) {
            const uint64_t external = degree_sum - internal_incidences;
            if (external == 0) break; // The entire component has size 2 or 3.
            const Extension extension = grow(graph, nodes, count, degree_sum, external, rng);
            nodes[count] = extension.node;
            for (unsigned i = 0; i < count; ++i) {
                if (i == extension.source || graph.adjacent(nodes[i], extension.node)) {
                    mask |= EDGE_BITS[i][count];
                    internal_incidences += 2;
                }
            }
            degree_sum += graph.degree(extension.node);
            ++count;
        }
        ++result.sizes[count];
        const OrbitRow& row = table[mask];
        for (unsigned i = 0; i < count; ++i) result.presence[nodes[i]] |= row[i];
    }
    return result;
}

void refuse_input_as_output(const char* input, const char* output) {
    struct stat input_stat{}, output_stat{};
    if (::stat(input, &input_stat) == 0 && ::stat(output, &output_stat) == 0 &&
        input_stat.st_dev == output_stat.st_dev && input_stat.st_ino == output_stat.st_ino) {
        throw std::runtime_error("graph and output must be different files");
    }
}

void write_presence(const std::string& filename, const std::vector<uint16_t>& presence) {
    std::ofstream output(filename, std::ios::binary);
    if (!output) throw std::runtime_error("cannot open output '" + filename + "'");
    char row[30];
    for (const uint16_t mask : presence) {
        for (unsigned orbit = 0; orbit < 15; ++orbit) {
            row[2 * orbit] = (mask & (1u << orbit)) ? '1' : '0';
            row[2 * orbit + 1] = orbit == 14 ? '\n' : ' ';
        }
        output.write(row, sizeof(row));
    }
    output.close();
    if (!output) throw std::runtime_error("error writing output '" + filename + "'");
}

double seconds(Clock::time_point before, Clock::time_point after) {
    return std::chrono::duration<double>(after - before).count();
}

} // namespace

int main(int argc, char* argv[]) {
    if (argc < 3 || argc > 5) {
        std::cerr << "Usage: " << argv[0] << " GRAPH OUTPUT [SAMPLES [SEED]]\n"
                  << "  Headerless undirected integer edge list; default SAMPLES=10000000.\n"
                  << "  Writes 15 binary ORCA4 node-orbit columns for each node ID.\n";
        return 1;
    }
    try {
        const uint64_t samples = argc >= 4
            ? decimal_integer(argv[3], "SAMPLES", std::numeric_limits<uint64_t>::max())
            : DEFAULT_SAMPLES;
        if (samples == 0) throw std::runtime_error("SAMPLES must be positive");
        const uint64_t seed = argc == 5
            ? decimal_integer(argv[4], "SEED", std::numeric_limits<uint64_t>::max())
            : random_seed();
        refuse_input_as_output(argv[1], argv[2]);
        const auto start = Clock::now();
        const Graph graph = read_graph(argv[1]);
        const auto read_done = Clock::now();
        const SampleResult result = sample(graph, samples, seed);
        const auto sample_done = Clock::now();
        write_presence(argv[2], result.presence);
        const auto write_done = Clock::now();
        std::cerr << "sampled_orbits: nodes=" << graph.node_count()
                  << " edges=" << graph.edges.size() << " samples=" << samples
                  << " seed=" << seed << " sampled2=" << result.sizes[2]
                  << " sampled3=" << result.sizes[3] << " sampled4=" << result.sizes[4]
                  << std::fixed << std::setprecision(3)
                  << " read_seconds=" << seconds(start, read_done)
                  << " sample_seconds=" << seconds(read_done, sample_done)
                  << " write_seconds=" << seconds(sample_done, write_done) << '\n';
    } catch (const std::exception& error) {
        std::cerr << "sampled_orbits: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
