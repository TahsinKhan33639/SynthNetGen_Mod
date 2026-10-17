#ifndef BXK4F_DISABLED_HPP
#define BXK4F_DISABLED_HPP

#include <array>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace bxk4f_disabled {

struct Config {
    std::array<bool, 96> transformations{};
    std::array<bool, 4> evc{};
};

template <std::size_t Size>
void parse_line(const std::string& line, std::array<bool, Size>& disabled,
                bool allow_zero, const std::string& label) {
    std::istringstream input(line);
    std::string token;
    bool saw_zero = false;
    std::size_t tokens = 0;
    while (input >> token) {
        unsigned int value = 0;
        for (char digit : token) {
            if (digit < '0' || digit > '9') {
                throw std::runtime_error(label + " must contain unsigned decimal integers");
            }
            value = value * 10 + static_cast<unsigned int>(digit - '0');
            if (value >= Size) {
                throw std::runtime_error(label + " integers must be in [1," +
                                         std::to_string(Size - 1) + "]");
            }
        }
        ++tokens;
        if (value == 0) {
            if (!allow_zero) {
                throw std::runtime_error(label + " integers must be in [1," +
                                         std::to_string(Size - 1) + "]");
            }
            saw_zero = true;
        } else {
            disabled[value] = true;
        }
    }
    if (saw_zero && tokens != 1) {
        throw std::runtime_error(label + " must use a single 0 for no disabled EVC transformations");
    }
}

inline Config load(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open disabled configuration: " + path);
    std::string transformations, evc, extra;
    if (!std::getline(input, transformations) || !std::getline(input, evc)) {
        throw std::runtime_error("disabled configuration must contain exactly two lines: " + path);
    }
    if (std::getline(input, extra)) {
        throw std::runtime_error("disabled configuration must contain exactly two lines: " + path);
    }
    if (input.bad()) throw std::runtime_error("failed to read disabled configuration: " + path);
    Config config;
    parse_line(transformations, config.transformations, false, "disabled configuration line 1");
    // EVC exclusions are applied by the shell, but validate the complete file
    // here so direct bxk4f invocations use the same configuration contract.
    parse_line(evc, config.evc, true, "disabled configuration line 2");
    return config;
}

}  // namespace bxk4f_disabled

#endif
