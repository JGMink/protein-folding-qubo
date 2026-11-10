#include "qubo.h"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

// Read one line into `line`, dropping a trailing '\r' so Windows-style files parse too.
// Returns false at end of file. (`static` = file-private, same as in C.)
static bool read_line(std::istream& in, std::string& line) {
    if (!std::getline(in, line)) return false;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    return true;
}

static FoldLayout parse_fold_layout(const std::string& meta, int n);  // defined below

// Parse a QUBO1 file: header, "n <count>", "offset <value>", a "meta {...}" line (its "sequence" and
// "lattice_edges" become the fold layout), then "q i j coefficient" lines with i <= j. i == j fills the diagonal;
// i < j is stored in both bits' neighbor lists. Rejects duplicates, out-of-range indices,
// and non-finite numbers. `throw` jumps out like a goto-error; run_experiment catches it.
// (`istringstream row(line); row >> a >> b` is roughly sscanf on one line.)
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
    q.fold = parse_fold_layout(line.substr(5), q.n);  // everything after "meta "
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

// Offset + diagonal terms of set bits + each pair of set bits once (the `j > i` test avoids counting a pair twice).
// (`const auto& [j, weight]` unpacks each (int, double) pair into two names.)
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

// Compute every field from scratch; afterwards update_fields_after_flip maintains them incrementally.
std::vector<double> local_fields(const Qubo& q, const std::vector<unsigned char>& bits) {
    // A field is the coefficient of one bit, given all the other current bits.
    std::vector<double> result = q.diagonal;
    for (int i = 0; i < q.n; ++i)
        for (const auto& [j, weight] : q.neighbors[i])
            if (bits[j]) result[i] += weight;
    return result;
}

// Energy change for flipping bit `index`, in O(1). Sign trick: a 0 bit gives (1 - 0) = +1, a 1 bit gives -1.
double flip_delta(const std::vector<unsigned char>& bits, const std::vector<double>& fields, int index) {
    // Turning 0 to 1 adds its field; turning 1 to 0 subtracts it.
    return (1 - 2 * int(bits[index])) * fields[index];
}

// After bit `index` flips (change = +1 or -1), each neighbor's field shifts by weight * change; the caller toggles the bit.
void update_fields_after_flip(const Qubo& q, std::vector<double>& fields, int index, int change) {
    // Only neighbors' fields change after an accepted flip; the QUBO stays read-only.
    for (const auto& [neighbor, weight] : q.neighbors[index]) fields[neighbor] += weight * change;
}

// Position just after `"key":` in a compact JSON line, or npos. Just enough JSON handling for two meta fields.
static std::size_t value_start(const std::string& meta, const std::string& key) {
    const std::string quoted = "\"" + key + "\"";
    std::size_t at = meta.find(quoted);
    if (at == std::string::npos) return std::string::npos;
    at += quoted.size();
    while (at < meta.size() && meta[at] == ' ') ++at;
    if (at >= meta.size() || meta[at] != ':') return std::string::npos;
    ++at;
    while (at < meta.size() && meta[at] == ' ') ++at;
    return at < meta.size() ? at : std::string::npos;
}

// Build the FoldLayout from the text of the meta record: chain length = letters in "sequence", sites = n / chain
// length, "lattice_edges" fills next_to. Strict: anything unexpected throws, so a wrong layout can't mislabel a fold.
static FoldLayout parse_fold_layout(const std::string& meta, int n) {
    FoldLayout layout;
    // Chain length: count the characters between the quotes of "sequence":"ABC...".
    std::size_t at = value_start(meta, "sequence");
    if (at == std::string::npos || meta[at] != '"')
        throw std::runtime_error("the meta record needs a \"sequence\" string (the lattice is required)");
    const std::size_t close = meta.find('"', at + 1);
    if (close == std::string::npos || close == at + 1)
        throw std::runtime_error("invalid \"sequence\" in the meta record");
    layout.residues = static_cast<int>(close - at - 1);
    if (n % layout.residues != 0) throw std::runtime_error("variable count is not residues * sites");
    layout.sites = n / layout.residues;

    // "lattice_edges":[[u,v],...] nests two brackets deep: collect every integer, rejecting any other character.
    at = value_start(meta, "lattice_edges");
    if (at == std::string::npos || meta[at] != '[')
        throw std::runtime_error("the meta record needs \"lattice_edges\" (the lattice is required)");
    std::vector<long long> numbers;
    int depth = 0;
    for (std::size_t i = at; i < meta.size(); ++i) {
        const char c = meta[i];
        if (c == '[') ++depth;
        else if (c == ']') { if (--depth == 0) break; }
        else if (c >= '0' && c <= '9') {
            long long value = 0;
            for (; i < meta.size() && meta[i] >= '0' && meta[i] <= '9'; ++i) {
                value = value * 10 + (meta[i] - '0');
                if (value >= layout.sites) throw std::runtime_error("lattice edge index out of range");
            }
            --i;  // the outer loop's ++i steps past the last digit
            numbers.push_back(value);
        } else if (c != ',' && c != ' ') throw std::runtime_error("invalid \"lattice_edges\" in the meta record");
    }
    if (depth != 0 || numbers.empty() || numbers.size() % 2 != 0)
        throw std::runtime_error("invalid \"lattice_edges\" in the meta record");

    // Each edge makes the two sites neighbors of each other (stored both ways, like Qubo::neighbors).
    layout.next_to.resize(layout.sites);
    for (std::size_t k = 0; k < numbers.size(); k += 2) {
        const int u = static_cast<int>(numbers[k]), v = static_cast<int>(numbers[k + 1]);
        if (u == v) throw std::runtime_error("lattice edge joins a site to itself");
        layout.next_to[u].push_back(v);
        layout.next_to[v].push_back(u);
    }
    return layout;
}

// Three rules, as in qubo_bench.model.validate: one 1-bit per residue block, no shared site, consecutive residues
// adjacent. (Contact energy isn't needed here: the penalty terms are zero for a valid fold.)
bool is_valid_fold(const FoldLayout& layout, const std::vector<unsigned char>& bits) {
    const int residues = layout.residues, sites = layout.sites;
    if (residues < 1 || sites < 1 || bits.size() != static_cast<std::size_t>(residues) * sites ||
        layout.next_to.size() != static_cast<std::size_t>(sites))
        return false;
    std::vector<int> site_of(residues);          // site chosen by each residue
    std::vector<unsigned char> taken(sites, 0);  // which sites are already occupied
    for (int i = 0; i < residues; ++i) {
        int found = -1;
        for (int u = 0; u < sites; ++u) {
            if (!bits[i * sites + u]) continue;
            if (found != -1) return false;       // a second site for the same residue
            found = u;
        }
        if (found == -1 || taken[found]) return false;  // no site, or two residues on one site
        taken[found] = 1;
        site_of[i] = found;
    }
    for (int i = 0; i + 1 < residues; ++i) {
        bool adjacent = false;
        for (int s : layout.next_to[site_of[i]])  // range-for: visit each neighbor site in turn
            if (s == site_of[i + 1]) { adjacent = true; break; }
        if (!adjacent) return false;
    }
    return true;
}
