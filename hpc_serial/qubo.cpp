#include "qubo.h"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

static bool read_line(std::istream& in, std::string& line) {
    if (!std::getline(in, line)) return false;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    return true;
}

Qubo load_qubo(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open QUBO input: " + path);
    std::string line, word;
    if (!read_line(in, line) || line != "QUBO1") throw std::runtime_error("expected QUBO1 header");
    Qubo q;
    if (!read_line(in, line)) throw std::runtime_error("missing variable count");
    {
        std::istringstream row(line);
        if (!(row >> word >> q.n) || word != "n" || q.n < 1 || (row >> word))
            throw std::runtime_error("invalid variable count");
    }
    if (!read_line(in, line)) throw std::runtime_error("missing offset");
    {
        std::istringstream row(line);
        if (!(row >> word >> q.offset) || word != "offset" || !std::isfinite(q.offset) || (row >> word))
            throw std::runtime_error("invalid offset");
    }
    if (!read_line(in, line) || line.rfind("meta ", 0) != 0)
        throw std::runtime_error("missing metadata record");
    q.diagonal.assign(q.n, 0);
    q.neighbors.resize(q.n);
    std::unordered_set<std::uint64_t> seen;
    while (read_line(in, line)) {
        std::istringstream row(line);
        int i, j;
        double coefficient;
        if (!(row >> word >> i >> j >> coefficient) || word != "q" || (row >> word)
            || i < 0 || i > j || j >= q.n || !std::isfinite(coefficient))
            throw std::runtime_error("invalid QUBO coefficient record");
        std::uint64_t key = std::uint64_t(i) * std::uint64_t(q.n) + std::uint64_t(j);
        if (!seen.insert(key).second) throw std::runtime_error("duplicate QUBO coefficient");
        if (i == j) q.diagonal[i] = coefficient;
        else {
            q.neighbors[i].push_back({j, coefficient});
            q.neighbors[j].push_back({i, coefficient});
        }
    }
    return q;
}

double energy(const Qubo& q, const std::vector<unsigned char>& bits) {
    double result = q.offset;
    for (int i = 0; i < q.n; ++i) {
        if (!bits[i]) continue;
        result += q.diagonal[i];
        for (const auto& [j, weight] : q.neighbors[i])
            if (j > i && bits[j]) result += weight;
    }
    return result;
}

std::vector<double> local_fields(const Qubo& q, const std::vector<unsigned char>& bits) {
    // A field is the coefficient of one bit, given all the other current bits.
    std::vector<double> result = q.diagonal;
    for (int i = 0; i < q.n; ++i)
        for (const auto& [j, weight] : q.neighbors[i])
            if (bits[j]) result[i] += weight;
    return result;
}

double flip_delta(const std::vector<unsigned char>& bits, const std::vector<double>& fields, int index) {
    // Turning 0 to 1 adds its field; turning 1 to 0 subtracts it.
    return (1 - 2 * int(bits[index])) * fields[index];
}

void update_fields_after_flip(const Qubo& q, std::vector<double>& fields, int index, int change) {
    // Only neighbors' fields change after an accepted flip; the QUBO stays read-only.
    for (const auto& [neighbor, weight] : q.neighbors[index]) fields[neighbor] += weight * change;
}
