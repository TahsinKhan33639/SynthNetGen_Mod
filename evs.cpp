#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using Adjacency = std::vector<std::vector<int>>;
constexpr int MAX_ITERATIONS = 10000;
constexpr double RESIDUAL_TOLERANCE = 1e-10;

struct EigenvalueResult {
    double value = 0.0;
    double residual = 0.0;
    int iterations = 0;
    bool converged = true;
};

Adjacency read_graph(const std::string& filename) {
    std::ifstream input(filename);
    if (!input) throw std::runtime_error("cannot open graph file '" + filename + "'");

    // Compress labels as before; missing labels represent isolated vertices and
    // do not affect the largest eigenvalue. Count each undirected edge once,
    // matching expfy's simple-graph interpretation.
    std::unordered_map<long long, int> labels;
    std::vector<std::pair<int, int>> edges;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        std::istringstream parser(line);
        parser >> std::ws;
        if (parser.eof()) continue;
        long long first, second;
        if (!(parser >> first >> second)) {
            throw std::runtime_error("expected two integer vertex labels on line " +
                                     std::to_string(line_number));
        }
        parser >> std::ws;
        if (!parser.eof()) {
            throw std::runtime_error("extra data on edge-list line " +
                                     std::to_string(line_number));
        }
        if (first == second) continue;
        const auto index = [&](long long label) {
            const auto found = labels.find(label);
            if (found != labels.end()) return found->second;
            if (labels.size() >= static_cast<std::size_t>(std::numeric_limits<int>::max())) {
                throw std::runtime_error("too many vertices");
            }
            const int id = static_cast<int>(labels.size());
            labels.emplace(label, id);
            return id;
        };
        int u = index(first), v = index(second);
        if (u > v) std::swap(u, v);
        edges.emplace_back(u, v);
    }
    if (input.bad()) throw std::runtime_error("error reading graph file '" + filename + "'");
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    Adjacency adjacency(labels.size());
    for (const auto& edge : edges) {
        adjacency[edge.first].push_back(edge.second);
        adjacency[edge.second].push_back(edge.first);
    }
    return adjacency;
}

// A graph of maximum degree two is a union of paths, cycles and isolated
// vertices. Their exact eigenvalues avoid the O(n^2) iteration counts caused
// by the small spectral gap of long paths.
double degree_two_eigenvalue(const Adjacency& adjacency) {
    std::vector<bool> seen(adjacency.size(), false);
    std::vector<int> pending;
    double largest = 0.0;
    const double pi = std::acos(-1.0);
    for (std::size_t root = 0; root < adjacency.size(); ++root) {
        if (seen[root]) continue;
        pending.clear();
        pending.push_back(static_cast<int>(root));
        seen[root] = true;
        std::size_t count = 0;
        bool cycle = true;
        while (!pending.empty()) {
            const int vertex = pending.back();
            pending.pop_back();
            ++count;
            cycle = cycle && adjacency[vertex].size() == 2;
            for (int neighbor : adjacency[vertex]) {
                if (!seen[neighbor]) {
                    seen[neighbor] = true;
                    pending.push_back(neighbor);
                }
            }
        }
        const double value = cycle ? 2.0 :
            (count == 1 ? 0.0 : 2.0 * std::cos(pi / (static_cast<double>(count) + 1.0)));
        largest = std::max(largest, value);
    }
    return largest;
}

EigenvalueResult largest_eigenvalue(const Adjacency& adjacency,
                                   int max_iterations = MAX_ITERATIONS) {
    if (adjacency.empty()) return {};
    std::size_t maximum_degree = 0;
    double degree_sum = 0.0;
    for (const auto& neighbors : adjacency) {
        maximum_degree = std::max(maximum_degree, neighbors.size());
        degree_sum += static_cast<double>(neighbors.size());
    }
    if (maximum_degree <= 2) return {degree_two_eigenvalue(adjacency), 0.0, 0, true};

    // Both quantities are lower bounds on the spectral radius. This positive
    // shift makes the Perron eigenvalue strictly larger in magnitude than its
    // negative partner, including on bipartite graphs. Unlike a max-degree
    // shift, it also works efficiently for stars with a very high-degree hub.
    const double shift = std::max(std::sqrt(static_cast<double>(maximum_degree)),
                                  degree_sum / static_cast<double>(adjacency.size()));
    std::vector<double> x(adjacency.size(), 1.0 / std::sqrt(static_cast<double>(adjacency.size())));
    std::vector<double> ax(adjacency.size(), 0.0);
    EigenvalueResult result;
    result.converged = false;

    for (int iteration = 0; iteration <= max_iterations; ++iteration) {
        double value = 0.0;
        for (std::size_t i = 0; i < adjacency.size(); ++i) {
            double sum = 0.0;
            for (int neighbor : adjacency[i]) sum += x[neighbor];
            ax[i] = sum;
            value += x[i] * sum;
        }
        double residual_squared = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i) {
            const double error = ax[i] - value * x[i];
            residual_squared += error * error;
        }
        result = {value, std::sqrt(residual_squared), iteration, false};
        if (!std::isfinite(result.value) || !std::isfinite(result.residual)) {
            throw std::runtime_error("nonfinite eigenvalue iteration");
        }
        if (result.residual <= RESIDUAL_TOLERANCE * std::max(1.0, std::abs(value))) {
            result.converged = true;
            return result;
        }
        if (iteration == max_iterations) break;

        // Reuse ax for (A + shift*I)x and swap buffers after normalization.
        double norm_squared = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i) {
            ax[i] += shift * x[i];
            norm_squared += ax[i] * ax[i];
        }
        const double norm = std::sqrt(norm_squared);
        if (!std::isfinite(norm) || norm == 0.0) {
            throw std::runtime_error("invalid norm in eigenvalue iteration");
        }
        for (double& value_next : ax) value_next /= norm;
        x.swap(ax);
    }
    return result;
}

int run_graph(const std::string& filename, int max_iterations = MAX_ITERATIONS) {
    const Adjacency adjacency = read_graph(filename);
    const EigenvalueResult result = largest_eigenvalue(adjacency, max_iterations);
    if (!result.converged) {
        std::cerr << std::setprecision(std::numeric_limits<double>::max_digits10)
                  << "evs: eigenpair residual did not converge after " << result.iterations
                  << " iterations (estimate=" << result.value
                  << ", residual=" << result.residual << ").\n";
        return 2;
    }
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
              << result.value << '\n';
    if (!std::cout) throw std::runtime_error("cannot write eigenvalue output");
    return 0;
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " graph.el\n";
        return 1;
    }
    try {
        return run_graph(argv[1]);
    } catch (const std::exception& error) {
        std::cerr << "evs: " << error.what() << '\n';
        return 1;
    }
}
