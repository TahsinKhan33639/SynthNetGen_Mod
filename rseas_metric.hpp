#ifndef RSEAS_METRIC_HPP
#define RSEAS_METRIC_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace rseas_metric {

constexpr std::size_t DIM = 6;
constexpr double LOG_FLOOR = -15.0;

using Vec6 = std::array<double, DIM>;

struct Config {
    double error_range = 0.0;
    Vec6 weights{1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
};

inline double safe_log(double value) {
    return value > 0.0 ? std::log(value) : LOG_FLOOR;
}

// signed_errors[i] may have either sign. RSEAS uses its absolute value.
inline double score(const Vec6& signed_errors, const Config& config) {
    Vec6 distances{};
    double maximum_distance = 0.0;
    bool outside_error_range = false;

    for (std::size_t i = 0; i < DIM; ++i) {
        distances[i] = std::abs(signed_errors[i]);
        maximum_distance = std::max(maximum_distance, distances[i]);
        if (distances[i] > config.error_range) {
            outside_error_range = true;
        }
    }

    if (!outside_error_range) {
        return maximum_distance - config.error_range;
    }

    double result = 0.0;
    for (std::size_t i = 0; i < DIM; ++i) {
        const double excess = std::max(0.0, distances[i] - config.error_range);
        result += config.weights[i] * excess * excess;
    }
    return result;
}

}  // namespace rseas_metric

#endif  // RSEAS_METRIC_HPP
