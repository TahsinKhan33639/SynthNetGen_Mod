#ifndef RSEAS_METRIC_K45_HPP
#define RSEAS_METRIC_K45_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace rseas_metric {

constexpr double LOG_FLOOR = -15.0;

using Vec = std::vector<double>;

inline std::size_t dimension_for_k(int k) {
    if (k == 4) return 6;
    if (k == 5) return 21;
    throw std::invalid_argument("graphlet size k must be 4 or 5");
}

struct Config {
    double error_range = 0.0;
    Vec weights;
};

inline double safe_log(double value) {
    return value > 0.0 ? std::log(value) : LOG_FLOOR;
}

// signed_errors[i] may have either sign. RSEAS uses its absolute value.
inline double score(const Vec& signed_errors, const Config& config) {
    if (signed_errors.empty()) {
        throw std::invalid_argument("RSEAS requires at least one coordinate");
    }
    if (signed_errors.size() != config.weights.size()) {
        throw std::invalid_argument(
            "RSEAS error and weight vectors must have the same length");
    }
    if (!std::isfinite(config.error_range) || config.error_range < 0.0) {
        throw std::invalid_argument(
            "RSEAS error range must be finite and nonnegative");
    }

    double maximum_distance = 0.0;
    bool outside_error_range = false;

    for (std::size_t i = 0; i < signed_errors.size(); ++i) {
        if (!std::isfinite(signed_errors[i]) ||
            !std::isfinite(config.weights[i]) || config.weights[i] <= 0.0) {
            throw std::invalid_argument(
                "RSEAS errors must be finite and weights must be positive");
        }

        const double distance = std::abs(signed_errors[i]);
        maximum_distance = std::max(maximum_distance, distance);
        if (distance > config.error_range) {
            outside_error_range = true;
        }
    }

    if (!outside_error_range) {
        return maximum_distance - config.error_range;
    }

    double result = 0.0;
    for (std::size_t i = 0; i < signed_errors.size(); ++i) {
        const double distance = std::abs(signed_errors[i]);
        const double excess =
            std::max(0.0, distance - config.error_range);
        result += config.weights[i] * excess * excess;
    }
    return result;
}

}  // namespace rseas_metric

#endif  // RSEAS_METRIC_K45_HPP
